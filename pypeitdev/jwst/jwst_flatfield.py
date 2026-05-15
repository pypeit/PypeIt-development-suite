"""
THIS IS GENERATED WITH CLAUDE.

Standalone reproduction of the calwebb_spec2 NIRSpec flat-field interpolation.

This module reimplements the core flat-field interpolation algorithm from the
JWST pipeline's ``flat_field`` step (``jwst.flatfield.flat_field``) so that
PypeIt can generate the interpolated flat images directly from CRDS reference
files and a wavelength image, **without** running ``calwebb_spec2``.

The three flat-field components for NIRSpec are:

- **D-flat** (detector flat): Pixel-to-pixel sensitivity variations on the
  detector.  Stored as a 3-D cube (wavelength × y × x) in the CRDS ``dflat``
  reference file, plus a 1-D "fast variation" table.
- **S-flat** (spectrograph flat): Throughput variations in the spectrograph
  optics.  Stored as a 3-D cube (or 2-D image) in the CRDS ``sflat`` reference
  file, plus a 1-D "fast variation" table.
- **F-flat** (fore-optics flat): Filter-dependent throughput of the
  fore-optics.  For MSA data, stored per-quadrant in a
  ``NirspecQuadFlatModel``; for fixed-slit data, stored in a
  ``NirspecFlatModel``.  Contains a 1-D "fast variation" table and, for MSA,
  a per-shutter 3-D image cube.

At each pixel, the total flat is:

.. math::
    \\text{flat}(y, x) = D(y, x, \\lambda) \\times S(y, x, \\lambda) \\times F(\\lambda)

where the interpolation in :math:`\\lambda` is performed using the per-pixel
wavelength grid from the ``assign_wcs`` step.

.. include:: ../include/links.rst
"""
import math

import numpy as np
from astropy.io import fits

from IPython import embed

from pypeit import log


# Dispersion direction constants
HORIZONTAL = 1
VERTICAL = 2

MICRONS_100 = 1.0e-4  # 100 microns, in meters

# DQ flag bit values matching JWST convention
DQ_DO_NOT_USE = 1         # Bit 0
DQ_NO_FLAT_FIELD = 1024   # Bit 10  (2**10)
DQ_UNRELIABLE_FLAT = 2048 # Bit 11  (2**11)
BADFLAT = DQ_NO_FLAT_FIELD | DQ_DO_NOT_USE | DQ_UNRELIABLE_FLAT


def create_interpolated_flat(waveimg, slit_xstart, slit_xsize, slit_ystart,
                             slit_ysize, fflat_file=None, sflat_file=None,
                             dflat_file=None, exposure_type='NRS_MSASPEC',
                             slit_name=None, quadrant=None, msa_x=None,
                             msa_y=None, dispaxis=HORIZONTAL,
                             subarray_xstart=1, subarray_ystart=1):
    """
    Generate the interpolated flat-field image for a single NIRSpec slit.

    This is the top-level function that combines D-flat, S-flat, and F-flat
    components, reproducing the behavior of ``calwebb_spec2``'s ``flat_field``
    step for a single slit.

    Parameters
    ----------
    waveimg : `numpy.ndarray`_
        2-D wavelength image for this slit, in microns. Shape ``(ny, nx)``
        matching the extracted slit dimensions (JWST convention: row = y,
        col = x).
    slit_xstart : :obj:`int`
        1-indexed X start position of the slit on the full detector.
    slit_xsize : :obj:`int`
        Number of columns in the slit cutout.
    slit_ystart : :obj:`int`
        1-indexed Y start position of the slit on the full detector.
    slit_ysize : :obj:`int`
        Number of rows in the slit cutout.
    fflat_file : :obj:`str`, optional
        Path to the CRDS ``fflat`` reference file. If ``None``, the F-flat
        component is set to 1.0.
    sflat_file : :obj:`str`, optional
        Path to the CRDS ``sflat`` reference file. If ``None``, the S-flat
        component is set to 1.0.
    dflat_file : :obj:`str`, optional
        Path to the CRDS ``dflat`` reference file. If ``None``, the D-flat
        component is set to 1.0.
    exposure_type : :obj:`str`, optional
        JWST exposure type, e.g. ``'NRS_MSASPEC'``, ``'NRS_FIXEDSLIT'``.
        Default is ``'NRS_MSASPEC'``.
    slit_name : :obj:`str`, optional
        Name of the slit (e.g. ``'S200A1'``).  Needed for fixed-slit
        flat table lookup.
    quadrant : :obj:`int`, optional
        MSA quadrant number (0-indexed: 0–3).  Only needed for MSA
        ``fflat`` files.
    msa_x : :obj:`int`, optional
        MSA shutter X index (1-indexed).  Only needed for MSA ``fflat``.
    msa_y : :obj:`int`, optional
        MSA shutter Y index (1-indexed).  Only needed for MSA ``fflat``.
    dispaxis : :obj:`int`, optional
        Dispersion direction: 1 = horizontal, 2 = vertical.
    subarray_xstart : :obj:`int`, optional
        1-indexed X start of the subarray on the full detector. Default 1.
    subarray_ystart : :obj:`int`, optional
        1-indexed Y start of the subarray on the full detector. Default 1.

    Returns
    -------
    flat_2d : `numpy.ndarray`_
        The combined interpolated flat-field correction, same shape as
        ``waveimg``.
    flat_dq : `numpy.ndarray`_
        DQ array (uint32) for the flat.
    flat_err : `numpy.ndarray`_
        Error array for the flat.
    wl_out : `numpy.ndarray`_
        Copy of the wavelength image (for storage in the output product).
    """

    # Convert slit position to 0-indexed full-frame coordinates
    xstart = slit_xstart - 1 + subarray_xstart - 1
    ystart = slit_ystart - 1 + subarray_ystart - 1
    xstop = xstart + slit_xsize
    ystop = ystart + slit_ysize

    # Prepare the wavelength image
    wl = waveimg.copy()
    nan_mask = np.isnan(wl)
    if nan_mask.any():
        wl[nan_mask] = 0.0
    max_wl = np.nanmax(wl)
    if 0.0 < max_wl < MICRONS_100:
        log.warning('Wavelengths appear to be in meters; expected microns.')

    # --- F-flat (fore-optics) component ---
    f_flat, f_flat_dq, f_flat_err = _fore_optics_flat(
        wl, fflat_file, exposure_type, dispaxis, slit_name,
        quadrant=quadrant, msa_x=msa_x, msa_y=msa_y)

    # --- S-flat (spectrograph) component ---
    s_flat, s_flat_dq, s_flat_err = _spectrograph_flat(
        wl, sflat_file, xstart, xstop, ystart, ystop,
        exposure_type, dispaxis, slit_name)

    # --- D-flat (detector) component ---
    d_flat, d_flat_dq, d_flat_err = _detector_flat(
        wl, dflat_file, xstart, xstop, ystart, ystop,
        exposure_type, dispaxis, slit_name)

    # --- Combine the three components ---
    flat_2d = f_flat * s_flat * d_flat
    flat_dq = _combine_dq(f_flat_dq, s_flat_dq, d_flat_dq,
                          default_shape=flat_2d.shape)

    # Combine errors in quadrature (relative errors)
    sum_var = np.zeros_like(flat_2d)
    if f_flat_err is not None:
        safe_f = np.where(f_flat == 0, 1.0, f_flat)
        sum_var += (f_flat_err / safe_f) ** 2
    if s_flat_err is not None:
        safe_s = np.where(s_flat == 0, 1.0, s_flat)
        sum_var += (s_flat_err / safe_s) ** 2
    if d_flat_err is not None:
        safe_d = np.where(d_flat == 0, 1.0, d_flat)
        sum_var += (d_flat_err / safe_d) ** 2
    flat_err = flat_2d * np.sqrt(sum_var)

    # Set DO_NOT_USE pixels to flat = 1.0
    bad = np.bitwise_and(flat_dq, DQ_DO_NOT_USE).astype(bool)
    flat_2d[bad] = 1.0

    # Mask bad flat values
    mask = flat_2d <= 0.0
    if mask.any():
        log.info(f'{mask.sum()} flat-field values <= 0; masking.')
        flat_2d[mask] = np.nan
        flat_dq[mask] = np.bitwise_or(flat_dq[mask], BADFLAT)

    return flat_2d.astype(np.float32), flat_dq, flat_err.astype(np.float32), wl


# ========================================================================
#  Internal helpers — each component follows the same pattern as the
#  calwebb code but reads the reference file via plain astropy.io.fits
# ========================================================================

def _read_flat_reffile(filepath):
    """
    Open a NIRSpec flat-field CRDS reference file.

    Parameters
    ----------
    filepath : :obj:`str`
        Path to the FITS reference file.

    Returns
    -------
    hdul : `~astropy.io.fits.HDUList`
        Opened FITS file.
    """
    return fits.open(filepath)


def _read_image_wl_from_hdu(hdul, ext_name='WAVELENGTH'):
    """
    Read the 1-D wavelength grid for image planes from the reference file.

    The reference file stores wavelengths in a binary table extension named
    ``WAVELENGTH``.  This corresponds to ``read_image_wl()`` in the JWST
    pipeline.

    Parameters
    ----------
    hdul : `~astropy.io.fits.HDUList`
        Opened FITS reference file.
    ext_name : :obj:`str`, optional
        Extension name.  Default ``'WAVELENGTH'``.

    Returns
    -------
    wavelength : `numpy.ndarray`_
        1-D array of wavelengths, one per plane of the SCI cube.
    """
    try:
        tbl = hdul[ext_name].data
    except KeyError:
        log.warning(f'No {ext_name} extension found in reference file.')
        return np.array([1.0])

    wl = tbl['wavelength'].flatten().copy()
    # Remove NaN and non-positive values (assumed to be trailing padding)
    good = np.isfinite(wl) & (wl > 0.0)
    return wl[good]


def _read_flat_table_from_hdu(hdul, ext_name='FAST_VARIATION',
                              exposure_type='NRS_MSASPEC', slit_name=None):
    """
    Read the "fast variation" flat-field table from the reference file.

    This corresponds to ``read_flat_table()`` in the JWST pipeline.  The table
    contains a wavelength-dependent 1-D flat-field correction that varies
    rapidly with wavelength.

    Parameters
    ----------
    hdul : `~astropy.io.fits.HDUList`
        Opened FITS reference file.
    ext_name : :obj:`str`, optional
        Extension name.  Default ``'FAST_VARIATION'``.
    exposure_type : :obj:`str`, optional
        Exposure type for slit-name row selection.
    slit_name : :obj:`str`, optional
        Slit name for fixed-slit row selection.

    Returns
    -------
    tab_wl : `numpy.ndarray`_
        1-D wavelength array.
    tab_flat : `numpy.ndarray`_
        1-D flat-field values.
    tab_flat_err : `numpy.ndarray`_
        1-D flat-field error values.
    """
    FIXED_SLIT_TYPES = ['NRS_LAMP', 'NRS_AUTOWAVE', 'NRS_BRIGHTOBJ',
                        'NRS_FIXEDSLIT']

    try:
        data = hdul[ext_name].data
    except KeyError:
        # No fast-variation table — return unity flat
        log.warning(f'No {ext_name} extension in reference file; '
                    f'returning unity for fast-variation component.')
        return (np.array([0.5, 100.0]),
                np.array([1.0, 1.0]),
                np.array([0.0, 0.0]))

    colnames = [c.lower() for c in data.dtype.names]
    has_slit_col = 'slit_name' in colnames
    has_nelem = 'nelem' in colnames

    # Determine which row to use
    row = None
    if exposure_type in FIXED_SLIT_TYPES and has_slit_col and slit_name is not None:
        slit_name_lc = slit_name.lower()
        for i in range(len(data)):
            col_val = data['slit_name'][i].strip().lower()
            if col_val == 'any' or col_val == slit_name_lc:
                row = i
                break
        if row is None:
            log.warning(f'Slit name {slit_name} not found in flat table; '
                        f'using first row.')
            row = 0

    # Extract arrays
    if row is not None:
        tab_wl = data['wavelength'][row].copy().flatten()
        tab_flat = data['data'][row].copy().flatten()
        if 'error' in colnames:
            tab_flat_err = data['error'][row].copy().flatten()
        else:
            tab_flat_err = np.zeros_like(tab_flat)
        nelem = data['nelem'][row] if has_nelem else len(tab_wl)
    elif len(data['wavelength'].shape) > 1:
        tab_wl = data['wavelength'][0].copy().flatten()
        tab_flat = data['data'][0].copy().flatten()
        if 'error' in colnames:
            tab_flat_err = data['error'][0].copy().flatten()
        else:
            tab_flat_err = np.zeros_like(tab_flat)
        nelem = data['nelem'][0] if has_nelem else len(tab_wl)
    else:
        tab_wl = data['wavelength'].copy().flatten()
        tab_flat = data['data'].copy().flatten()
        if 'error' in colnames:
            tab_flat_err = data['error'].copy().flatten()
        else:
            tab_flat_err = np.zeros_like(tab_flat)
        nelem = len(tab_wl)

    # Truncate to nelem
    tab_wl = tab_wl[:nelem]
    tab_flat = tab_flat[:nelem]
    tab_flat_err = tab_flat_err[:nelem]

    # Filter out NaN, zero, and negative values
    good = np.isfinite(tab_wl) & np.isfinite(tab_flat) & (tab_wl > 0) & (tab_flat != 0)
    tab_wl = tab_wl[good]
    tab_flat = tab_flat[good]
    tab_flat_err = tab_flat_err[good]

    return tab_wl, tab_flat, tab_flat_err


def _interpolate_flat_cube(image_flat, image_dq, image_err, image_wl, wl):
    """
    Interpolate a 3-D flat-field cube to a 2-D image at the given wavelengths.

    Linear interpolation along the wavelength axis (first axis) of the cube.
    This corresponds to ``interpolate_flat()`` in the JWST pipeline.

    Parameters
    ----------
    image_flat : `numpy.ndarray`_
        3-D flat array, shape ``(nz, ny, nx)``.
    image_dq : `numpy.ndarray`_
        2-D or 3-D DQ array.
    image_err : `numpy.ndarray`_
        2-D or 3-D error array.
    image_wl : `numpy.ndarray`_
        1-D wavelength grid for the cube planes, length ``nz``.
    wl : `numpy.ndarray`_
        2-D target wavelength array, shape ``(ny, nx)``.

    Returns
    -------
    flat_2d : `numpy.ndarray`_
        Interpolated 2-D flat, same shape as ``wl``.
    flat_dq : `numpy.ndarray`_
        DQ array.
    flat_err : `numpy.ndarray`_
        Error array.
    """
    if image_flat.ndim < 3:
        return image_flat.copy(), image_dq.copy(), image_err.copy()

    nz, ysize, xsize = image_flat.shape
    if nz == 1:
        return (image_flat[0].copy(),
                image_dq[0].copy() if image_dq.ndim == 3 else image_dq.copy(),
                image_err[0].copy() if image_err.ndim == 3 else image_err.copy())

    grid = np.indices((ysize, xsize), dtype=np.intp)
    iy = grid[0]
    ix = grid[1]

    # Find the interpolation interval for each pixel
    k = np.full(wl.shape, -1, dtype=np.intp)
    k = np.where(wl <= image_wl[0], 0, k)
    k = np.where(wl >= image_wl[-1], nz - 2, k)

    for k_test in range(nz - 1):
        in_interval = (wl >= image_wl[k_test]) & (wl < image_wl[k_test + 1])
        unassigned = (k == -1)
        k = np.where(unassigned & in_interval, k_test, k)
        if np.all(k >= 0):
            break

    # Linear interpolation
    denom = image_wl[k + 1] - image_wl[k]
    zero_denom = (denom == 0.0)
    denom = np.where(zero_denom, 1.0, denom)
    p = np.where(zero_denom, 0.0, (wl - image_wl[k]) / denom)
    q = 1.0 - p

    flat_2d = q * image_flat[k, iy, ix] + p * image_flat[k + 1, iy, ix]

    if image_err.ndim == 2:
        flat_err = image_err.copy()
    else:
        flat_err = q * image_err[k, iy, ix] + p * image_err[k + 1, iy, ix]

    if image_dq.ndim == 2:
        flat_dq = image_dq.copy()
    else:
        flat_dq = np.where(
            p == 0.0,
            image_dq[k, iy, ix],
            np.bitwise_or(image_dq[k, iy, ix], image_dq[k + 1, iy, ix]))
        flat_bad = np.bitwise_and(flat_dq, DQ_DO_NOT_USE).astype(bool)
        flat_2d[flat_bad] = 1.0

    # Flag pixels outside the reference wavelength range
    flat_dq[wl < image_wl[0]] |= DQ_DO_NOT_USE
    flat_dq[wl > image_wl[-1]] |= DQ_DO_NOT_USE
    flat_bad = np.bitwise_and(flat_dq, DQ_DO_NOT_USE).astype(bool)
    flat_2d[flat_bad] = 1.0

    return flat_2d.astype(image_flat.dtype), flat_dq, flat_err


def _clean_wl(wl, dispaxis):
    """
    Replace zeros and NaNs in the wavelength array with nearest-neighbor
    averages along the cross-dispersion direction.

    Parameters
    ----------
    wl : `numpy.ndarray`_
        2-D wavelength array.
    dispaxis : :obj:`int`
        1 = horizontal, 2 = vertical.

    Returns
    -------
    wl_c : `numpy.ndarray`_
        Cleaned copy of ``wl``.
    """
    wl_c = wl.copy().astype(np.float64)
    wl_c[wl_c <= 0.0] = np.nan
    shape = wl_c.shape

    if dispaxis == HORIZONTAL:
        for i in range(shape[1]):
            col = wl_c[:, i]
            if np.any(np.isfinite(col)):
                replacement = np.nanmean(col)
                col[np.isnan(col)] = replacement
        # Fill leading/trailing all-NaN columns by edge extension
        finite_cols = np.any(np.isfinite(wl_c), axis=0)
        if finite_cols.any():
            first = np.argmax(finite_cols)
            last = shape[1] - 1 - np.argmax(finite_cols[::-1])
            if first > 0:
                wl_c[:, :first] = wl_c[:, first:first + 1]
            if last < shape[1] - 1:
                wl_c[:, last + 1:] = wl_c[:, last:last + 1]
    elif dispaxis == VERTICAL:
        for j in range(shape[0]):
            row = wl_c[j, :]
            if np.any(np.isfinite(row)):
                replacement = np.nanmean(row)
                row[np.isnan(row)] = replacement
        finite_rows = np.any(np.isfinite(wl_c), axis=1)
        if finite_rows.any():
            first = np.argmax(finite_rows)
            last = shape[0] - 1 - np.argmax(finite_rows[::-1])
            if first > 0:
                wl_c[:first, :] = wl_c[first:first + 1, :]
            if last < shape[0] - 1:
                wl_c[last + 1:, :] = wl_c[last:last + 1, :]

    return wl_c


def _combine_fast_slow(wl, flat_2d, flat_dq, flat_err, tab_wl, tab_flat,
                       tab_flat_err, dispaxis):
    """
    Multiply the 2-D "slow" flat by the 1-D "fast variation" table.

    Uses 3-point Gaussian quadrature to average the fast-variation table over
    each pixel's wavelength bin.  Reproduces ``combine_fast_slow()`` from the
    JWST pipeline.

    Parameters
    ----------
    wl : `numpy.ndarray`_
        2-D wavelength array (microns).
    flat_2d : `numpy.ndarray`_ or :obj:`float`
        The slowly-varying 2-D flat component.
    flat_dq : `numpy.ndarray`_ or ``None``
        DQ array for the slow component.
    flat_err : `numpy.ndarray`_ or ``None``
        Error array for the slow component.
    tab_wl : `numpy.ndarray`_
        1-D wavelength grid from the fast-variation table.
    tab_flat : `numpy.ndarray`_
        1-D flat values from the fast-variation table.
    tab_flat_err : `numpy.ndarray`_
        1-D flat errors.
    dispaxis : :obj:`int`
        1 = horizontal, 2 = vertical.

    Returns
    -------
    combined_flat : `numpy.ndarray`_
        ``flat_2d * interpolated_tab_flat``.
    combined_dq : `numpy.ndarray`_
        Updated DQ array.
    combined_err : `numpy.ndarray`_
        Updated error array.
    """
    wl_c = _clean_wl(wl, dispaxis)
    dwl = np.zeros_like(wl_c)

    if flat_dq is None:
        combined_dq = np.zeros(wl.shape, dtype=np.uint32)
    else:
        combined_dq = flat_dq.copy()

    if flat_err is None:
        combined_err = np.zeros(wl.shape, dtype=np.float64)
    else:
        combined_err = flat_err.copy().astype(np.float64)

    if dispaxis == HORIZONTAL:
        dwl[:, :-1] = wl_c[:, 1:] - wl_c[:, :-1]
        dwl[:, -1] = dwl[:, -2]
    elif dispaxis == VERTICAL:
        dwl[:-1, :] = wl_c[1:, :] - wl_c[:-1, :]
        dwl[-1, :] = dwl[-2, :]

    # 3-point Gaussian quadrature (averaged over the pixel's wavelength bin)
    d = math.sqrt(0.6) / 2.0
    dx = np.array([-d, 0.0, d])
    wgt = np.array([5.0, 8.0, 5.0]) / 18.0

    values = np.zeros_like(wl_c)
    for offset, weight in zip(dx, wgt):
        wavelengths = wl_c + dwl * offset
        values += weight * np.interp(wavelengths, tab_wl, tab_flat,
                                     left=np.nan, right=np.nan)

    # Simple interpolation for errors
    error_value = np.interp(wl_c, tab_wl, tab_flat_err,
                            left=np.nan, right=np.nan)

    # Handle bad wavelength values
    bad = wl <= 0
    values[bad] = 1.0
    error_value[bad] = 0.0

    # Handle missing (out-of-range) values
    missing = np.isnan(values)
    values[missing] = 1.0
    error_value[missing] = 0.0
    combined_dq[missing] |= DQ_NO_FLAT_FIELD
    combined_dq[missing] |= DQ_DO_NOT_USE

    # Add new 1-D errors and input 2-D errors in quadrature
    v1 = np.square(error_value * flat_2d)
    v2 = np.square(combined_err * values)
    combined_err = np.sqrt(np.nansum([v1, v2], axis=0))

    return (flat_2d * values).astype(np.float32), combined_dq, combined_err.astype(np.float32)


def _fore_optics_flat(wl, fflat_file, exposure_type, dispaxis, slit_name,
                      quadrant=None, msa_x=None, msa_y=None):
    """
    Compute the fore-optics (F-flat) component of the NIRSpec flat field.

    Parameters
    ----------
    wl : `numpy.ndarray`_
        2-D wavelength image (microns).
    fflat_file : :obj:`str` or ``None``
        Path to the ``fflat`` CRDS reference file.
    exposure_type : :obj:`str`
        Exposure type.
    dispaxis : :obj:`int`
        Dispersion direction (1=horizontal, 2=vertical).
    slit_name : :obj:`str` or ``None``
        Slit name for fixed-slit lookups.
    quadrant : :obj:`int` or ``None``
        0-indexed MSA quadrant (for MSA fflat only).
    msa_x : :obj:`int` or ``None``
        1-indexed MSA shutter X index (for MSA fflat only).
    msa_y : :obj:`int` or ``None``
        1-indexed MSA shutter Y index (for MSA fflat only).

    Returns
    -------
    f_flat : `numpy.ndarray`_
        F-flat correction array.
    f_flat_dq : `numpy.ndarray`_ or ``None``
        DQ array.
    f_flat_err : `numpy.ndarray`_ or ``None``
        Error array.
    """
    if fflat_file is None:
        return (np.ones(wl.shape, dtype=np.float32), None, None)

    hdul = fits.open(fflat_file)

    # For MSA data with quadrant-based fflat (NirspecQuadFlatModel),
    # there are multiple SCI/DQ/ERR/WAVELENGTH/FAST_VARIATION extensions,
    # one per quadrant. The naming convention uses EXTVER for quadrant indexing.
    # For fixed-slit, there is a single set of extensions.
    if exposure_type == 'NRS_MSASPEC' and quadrant is not None:
        # MSA: read the per-quadrant fast_variation table
        # The FAST_VARIATION extensions are numbered by quadrant (EXTVER)
        ext_ver = quadrant + 1  # 1-indexed
        tab_wl, tab_flat, tab_flat_err = _read_flat_table_from_hdu_quadrant(
            hdul, quadrant=quadrant, exposure_type=exposure_type,
            slit_name=slit_name)

        # MSA: combine the per-shutter 1-D image flat with the table
        if msa_x is not None and msa_y is not None:
            # Read the per-quadrant image cube and wavelength
            sci_ext = _get_quadrant_ext(hdul, 'SCI', quadrant)
            err_ext = _get_quadrant_ext(hdul, 'ERR', quadrant)
            image_wl = _read_image_wl_quadrant(hdul, quadrant)

            if sci_ext is not None and len(sci_ext.shape) >= 3:
                msa_y_idx = msa_y - 1  # 0-indexed
                msa_x_idx = msa_x - 1
                one_d_flat = sci_ext[:, msa_y_idx, msa_x_idx]
                one_d_err = err_ext[:, msa_y_idx, msa_x_idx]
                # Interpolate the 1-D MSA flat at each table wavelength and
                # multiply into tab_flat
                tab_flat = tab_flat * np.interp(tab_wl, image_wl, one_d_flat,
                                               1.0, 1.0)
                f_flat_err_1d = np.interp(wl, image_wl, one_d_err, 0.0, 0.0)
            else:
                f_flat_err_1d = None
        else:
            f_flat_err_1d = None

        # The F-flat "image" for MSA is effectively 1-D (applied via table),
        # so flat_2d = 1.0 (scalar)
        flat_2d = 1.0
        f_flat_dq = None
        f_flat_err = f_flat_err_1d

        # Combine 2-D (=1.0) and 1-D components
        f_flat, f_flat_dq, f_flat_err = _combine_fast_slow(
            wl, flat_2d, f_flat_dq, f_flat_err,
            tab_wl, tab_flat, tab_flat_err, dispaxis)

    else:
        # Fixed-slit or other non-quadrant mode
        tab_wl, tab_flat, tab_flat_err = _read_flat_table_from_hdu(
            hdul, exposure_type=exposure_type, slit_name=slit_name)
        flat_2d = 1.0
        f_flat_dq = None
        f_flat_err = None

        f_flat, f_flat_dq, f_flat_err = _combine_fast_slow(
            wl, flat_2d, f_flat_dq, f_flat_err,
            tab_wl, tab_flat, tab_flat_err, dispaxis)

    hdul.close()

    # Flag NaN and zero values
    _flag_bad_flat(f_flat, f_flat_dq)

    # For F-flat, bad pixels get NaN (not 1.0) because F-flat provides
    # flux calibration scaling
    if f_flat_dq is not None:
        bad = np.bitwise_and(f_flat_dq, DQ_DO_NOT_USE).astype(bool)
        f_flat[bad] = np.nan

    return f_flat, f_flat_dq, f_flat_err


def _spectrograph_flat(wl, sflat_file, xstart, xstop, ystart, ystop,
                       exposure_type, dispaxis, slit_name):
    """
    Compute the spectrograph (S-flat) component of the NIRSpec flat field.

    Parameters
    ----------
    wl : `numpy.ndarray`_
        2-D wavelength image (microns).
    sflat_file : :obj:`str` or ``None``
        Path to the ``sflat`` CRDS reference file.
    xstart, xstop, ystart, ystop : :obj:`int`
        0-indexed pixel coordinates for the slit cutout on the full detector.
    exposure_type : :obj:`str`
        Exposure type.
    dispaxis : :obj:`int`
        Dispersion direction.
    slit_name : :obj:`str` or ``None``
        Slit name.

    Returns
    -------
    s_flat : `numpy.ndarray`_
    s_flat_dq : `numpy.ndarray`_ or ``None``
    s_flat_err : `numpy.ndarray`_ or ``None``
    """
    if sflat_file is None:
        return (np.ones(wl.shape, dtype=np.float32), None, None)

    hdul = fits.open(sflat_file)

    # Read the fast-variation table
    tab_wl, tab_flat, tab_flat_err = _read_flat_table_from_hdu(
        hdul, exposure_type=exposure_type, slit_name=slit_name)

    if tab_wl.max() < MICRONS_100:
        log.warning('Wavelengths in S-flat table appear to be in meters.')

    # Read the image (2-D or 3-D)
    sci = hdul['SCI'].data
    dq = hdul['DQ'].data
    err = hdul['ERR'].data

    if sci.ndim == 3:
        # 3-D cube: extract slit region and interpolate along wavelength
        image_flat = sci[:, ystart:ystop, xstart:xstop]
        image_dq = dq[:, ystart:ystop, xstart:xstop] if dq.ndim == 3 \
            else dq[ystart:ystop, xstart:xstop]
        image_err = err[:, ystart:ystop, xstart:xstop] if err.ndim == 3 \
            else err[ystart:ystop, xstart:xstop]
        image_wl = _read_image_wl_from_hdu(hdul)

        if image_wl.max() < MICRONS_100:
            log.warning('Wavelengths in S-flat image appear to be in meters.')

        flat_2d, s_flat_dq, s_flat_err = _interpolate_flat_cube(
            image_flat, image_dq, image_err, image_wl, wl)
    else:
        flat_2d = sci[ystart:ystop, xstart:xstop].copy()
        s_flat_dq = dq[ystart:ystop, xstart:xstop].copy()
        s_flat_err = err[ystart:ystop, xstart:xstop].copy()

    hdul.close()

    # Flag NaN and zero
    _flag_bad_flat(flat_2d, s_flat_dq)

    # Reset bad pixels to 1.0
    if s_flat_dq is not None:
        bad = np.bitwise_and(s_flat_dq, DQ_DO_NOT_USE).astype(bool)
        flat_2d[bad] = 1.0

    # Combine slow (2-D) and fast (1-D table) components
    s_flat, s_flat_dq, s_flat_err = _combine_fast_slow(
        wl, flat_2d, s_flat_dq, s_flat_err,
        tab_wl, tab_flat, tab_flat_err, dispaxis)

    return s_flat, s_flat_dq, s_flat_err


def _detector_flat(wl, dflat_file, xstart, xstop, ystart, ystop,
                   exposure_type, dispaxis, slit_name):
    """
    Compute the detector (D-flat) component of the NIRSpec flat field.

    Parameters
    ----------
    wl : `numpy.ndarray`_
        2-D wavelength image (microns).
    dflat_file : :obj:`str` or ``None``
        Path to the ``dflat`` CRDS reference file.
    xstart, xstop, ystart, ystop : :obj:`int`
        0-indexed pixel coordinates for the slit cutout on the full detector.
    exposure_type : :obj:`str`
        Exposure type.
    dispaxis : :obj:`int`
        Dispersion direction.
    slit_name : :obj:`str` or ``None``
        Slit name.

    Returns
    -------
    d_flat : `numpy.ndarray`_
    d_flat_dq : `numpy.ndarray`_ or ``None``
    d_flat_err : `numpy.ndarray`_ or ``None``
    """
    if dflat_file is None:
        return (np.ones(wl.shape, dtype=np.float32), None, None)

    hdul = fits.open(dflat_file)

    # Read the fast-variation table
    tab_wl, tab_flat, tab_flat_err = _read_flat_table_from_hdu(
        hdul, exposure_type=exposure_type, slit_name=slit_name)

    if tab_wl.max() < MICRONS_100:
        log.warning('Wavelengths in D-flat table appear to be in meters.')

    # Read the 3-D image cube
    sci = hdul['SCI'].data
    dq = hdul['DQ'].data
    err = hdul['ERR'].data

    image_flat = sci[:, ystart:ystop, xstart:xstop]
    image_dq = dq[..., ystart:ystop, xstart:xstop]
    image_err = err[..., ystart:ystop, xstart:xstop]
    image_wl = _read_image_wl_from_hdu(hdul)

    if image_wl.max() < MICRONS_100:
        log.warning('Wavelengths in D-flat image appear to be in meters.')

    flat_2d, d_flat_dq, d_flat_err = _interpolate_flat_cube(
        image_flat, image_dq, image_err, image_wl, wl)

    hdul.close()

    # Flag NaN and zero
    _flag_bad_flat(flat_2d, d_flat_dq)

    # Reset bad pixels to 1.0
    if d_flat_dq is not None:
        bad = np.bitwise_and(d_flat_dq, DQ_DO_NOT_USE).astype(bool)
        flat_2d[bad] = 1.0

    # Combine slow (2-D) and fast (1-D table)
    d_flat, d_flat_dq, d_flat_err = _combine_fast_slow(
        wl, flat_2d, d_flat_dq, d_flat_err,
        tab_wl, tab_flat, tab_flat_err, dispaxis)

    return d_flat, d_flat_dq, d_flat_err


def _flag_bad_flat(flat_2d, flat_dq):
    """
    In-place flag NaN and zero-valued flat pixels in the DQ array.

    Parameters
    ----------
    flat_2d : `numpy.ndarray`_
        The flat-field image.
    flat_dq : `numpy.ndarray`_ or ``None``
        The DQ array (modified in-place).
    """
    if flat_dq is None:
        return
    bad_flag = DQ_DO_NOT_USE | DQ_NO_FLAT_FIELD
    nan_mask = np.isnan(flat_2d)
    flat_dq[nan_mask] |= bad_flag
    zero_mask = (flat_2d == 0.0)
    flat_dq[zero_mask] |= bad_flag


def _combine_dq(f_flat_dq, s_flat_dq, d_flat_dq, default_shape):
    """
    Combine the three component DQ arrays via bitwise OR.

    Parameters
    ----------
    f_flat_dq, s_flat_dq, d_flat_dq : `numpy.ndarray`_ or ``None``
        Individual DQ arrays.
    default_shape : :obj:`tuple`
        Shape for a zero-filled fallback DQ array.

    Returns
    -------
    flat_dq : `numpy.ndarray`_
        Combined DQ array (uint32).
    """
    dq_list = [dq for dq in [f_flat_dq, s_flat_dq, d_flat_dq]
               if dq is not None]

    flat_dq = np.zeros(default_shape, dtype=np.uint32)
    if len(dq_list) == 0:
        flat_dq |= BADFLAT
    else:
        for dq in dq_list:
            flat_dq = np.bitwise_or(flat_dq, dq.astype(np.uint32))

    # Flag DO_NOT_USE where NO_FLAT_FIELD is set
    no_flat = np.bitwise_and(flat_dq, DQ_NO_FLAT_FIELD).astype(bool)
    flat_dq[no_flat] |= DQ_DO_NOT_USE

    # Flag full BADFLAT where DO_NOT_USE is set
    bad = np.bitwise_and(flat_dq, DQ_DO_NOT_USE).astype(bool)
    flat_dq[bad] |= BADFLAT

    return flat_dq


# ========================================================================
#  Helpers for NirspecQuadFlatModel (MSA fflat) — read by EXTVER
# ========================================================================

def _get_quadrant_ext(hdul, extname, quadrant):
    """
    Get an extension from a NirspecQuadFlatModel FITS file for a specific
    quadrant.

    The convention is multiple extensions with the same EXTNAME but
    different EXTVER (1-indexed quadrant + 1).

    Parameters
    ----------
    hdul : `~astropy.io.fits.HDUList`
        Opened FITS file.
    extname : :obj:`str`
        Extension name (e.g. ``'SCI'``, ``'ERR'``).
    quadrant : :obj:`int`
        0-indexed quadrant number.

    Returns
    -------
    data : `numpy.ndarray`_ or ``None``
        The extension data, or ``None`` if not found.
    """
    extver = quadrant + 1
    try:
        return hdul[extname, extver].data
    except KeyError:
        log.warning(f'Extension {extname} ver={extver} not found.')
        return None


def _read_image_wl_quadrant(hdul, quadrant):
    """
    Read the wavelength grid for a specific quadrant.

    Parameters
    ----------
    hdul : `~astropy.io.fits.HDUList`
        Opened FITS file.
    quadrant : :obj:`int`
        0-indexed quadrant number.

    Returns
    -------
    wavelength : `numpy.ndarray`_
        1-D wavelength array.
    """
    extver = quadrant + 1
    try:
        tbl = hdul['WAVELENGTH', extver].data
    except KeyError:
        return np.array([1.0])

    wl = tbl['wavelength'].flatten().copy()
    good = np.isfinite(wl) & (wl > 0.0)
    return wl[good]


def _read_flat_table_from_hdu_quadrant(hdul, quadrant, exposure_type,
                                       slit_name=None):
    """
    Read the fast-variation flat table for a specific MSA quadrant.

    Parameters
    ----------
    hdul : `~astropy.io.fits.HDUList`
        Opened FITS file.
    quadrant : :obj:`int`
        0-indexed quadrant number.
    exposure_type : :obj:`str`
        Exposure type.
    slit_name : :obj:`str` or ``None``
        Slit name.

    Returns
    -------
    tab_wl, tab_flat, tab_flat_err : `numpy.ndarray`_
        1-D arrays.
    """
    extver = quadrant + 1
    ext_name = ('FAST_VARIATION', extver)
    try:
        data = hdul[ext_name].data
    except KeyError:
        log.warning(f'No FAST_VARIATION extension for quadrant {quadrant}; '
                    f'returning unity.')
        return (np.array([0.5, 100.0]),
                np.array([1.0, 1.0]),
                np.array([0.0, 0.0]))

    colnames = [c.lower() for c in data.dtype.names]
    tab_wl = data['wavelength'][0].copy().flatten()
    tab_flat = data['data'][0].copy().flatten()
    if 'error' in colnames:
        tab_flat_err = data['error'][0].copy().flatten()
    else:
        tab_flat_err = np.zeros_like(tab_flat)

    if 'nelem' in colnames:
        nelem = data['nelem'][0]
        tab_wl = tab_wl[:nelem]
        tab_flat = tab_flat[:nelem]
        tab_flat_err = tab_flat_err[:nelem]

    # Filter
    good = (np.isfinite(tab_wl) & np.isfinite(tab_flat)
            & (tab_wl > 0) & (tab_flat != 0))
    return tab_wl[good], tab_flat[good], tab_flat_err[good]


