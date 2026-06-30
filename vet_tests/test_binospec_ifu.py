"""
Tests for the Binospec IFU cube building and fiber extraction tools.
"""
import re
from pathlib import Path

import numpy as np
from astropy.io import fits

from pypeit import spec2dobj
from pypeit.onespec import OneSpec
from pypeit.specobjs import SpecObjs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts.binospec_ifu_cube import BinospecIFUCube
from pypeit.core.datacube import (
    load_fibers,
    project_to_sky,
    resample_and_combine,
    write_onespec,
)


def test_binospec_ifu_cube_from_spec1d(redux_out):
    """Test building a datacube from a Binospec IFU spec1d file."""
    sci_dir = Path(redux_out) / 'mmt_binospec_ifu' / 'G270' / 'Science'
    spec1d_files = sorted(sci_dir.glob('spec1d_*.fits'))
    assert len(spec1d_files) > 0, f"No spec1d files found in {sci_dir}"

    spec1d_file = spec1d_files[0]
    output_file = sci_dir / 'test_cube_spec1d.fits'

    # Build cube from spec1d
    pargs = BinospecIFUCube.parse_args([str(spec1d_file), '-o', str(output_file)])
    BinospecIFUCube.main(pargs)

    # Validate the output datacube
    assert output_file.exists(), f"Datacube not created: {output_file}"

    with fits.open(output_file) as hdu:
        # Check expected HDU structure
        assert 'FLUX' in hdu, "FLUX HDU missing from datacube"
        assert 'VAR' in hdu, "VAR HDU missing from datacube"

        flux = hdu['FLUX'].data
        var = hdu['VAR'].data

        # Check 3D shape (n_wave, ny, nx)
        assert flux.ndim == 3, f"FLUX is not 3D: {flux.ndim}D"
        assert var.shape == flux.shape, "VAR shape does not match FLUX shape"

        # Check that the cube contains real data
        assert np.any(flux != 0), "FLUX cube is all zeros"
        assert np.any(var > 0), "VAR cube has no positive values"

        # Check WCS keywords
        hdr = hdu['FLUX'].header
        assert hdr['CTYPE1'] == 'RA---TAN'
        assert hdr['CTYPE2'] == 'DEC--TAN'
        assert hdr['CTYPE3'] == 'WAVE'

    # Clean up
    output_file.unlink()


def test_binospec_ifu_extract(redux_out):
    """Test the Binospec IFU fiber-extraction pipeline on a real spec1d file.

    ``binospec_ifu_extract`` is an interactive GUI, so its ``main()`` cannot
    run headless.  This drives the same data-loading and combination helpers
    that ``main()`` uses up to the point the GUI is launched: load the science
    fibers, project them to sky, combine a subset, and write a OneSpec file.
    """
    sci_dir = Path(redux_out) / 'mmt_binospec_ifu' / 'G270' / 'Science'
    spec1d_files = sorted(sci_dir.glob('spec1d_*.fits'))
    assert len(spec1d_files) > 0, f"No spec1d files found in {sci_dir}"

    spec1d_file = spec1d_files[0]

    # Load the spec1d and its raw header (mirrors binospec_ifu_extract.main).
    sobjs = SpecObjs.from_fitsfile(str(spec1d_file))
    assert sobjs.nobj > 0, f"No objects in {spec1d_file}"
    with fits.open(spec1d_file) as hdu:
        raw_hdr = hdu[0].header.copy()

    spectrograph = load_spectrograph(raw_hdr['PYP_SPEC'])
    assert spectrograph.name == 'mmt_binospec_ifu'

    targetx, targety = spectrograph.load_sky_layout()
    fibers = load_fibers(sobjs, spectrograph, targetx, targety, prefix='OPT')
    assert len(fibers) > 0, "No science fibers loaded from spec1d"

    # Project fiber instrument-frame positions to sky coordinates.
    x = np.array([f['x'] for f in fibers])
    y = np.array([f['y'] for f in fibers])
    ra, dec = project_to_sky(x, y, raw_hdr, spectrograph)
    assert ra.shape == (len(fibers),)
    assert dec.shape == (len(fibers),)
    assert np.all(np.isfinite(ra)) and np.all(np.isfinite(dec))

    # Combine a handful of fibers into a single 1D spectrum.
    subset = fibers[:min(5, len(fibers))]
    out_wave, out_flux, out_ivar = resample_and_combine(
        [f['wave'] for f in subset],
        [f['flux'] for f in subset],
        [f['ivar'] for f in subset])
    assert out_wave.size > 0, "Combined spectrum is empty"
    assert out_flux.shape == out_wave.shape
    assert out_ivar.shape == out_wave.shape
    assert np.any(np.isfinite(out_flux)), "Combined flux is all non-finite"

    # Write and round-trip the extracted OneSpec.
    output_file = sci_dir / 'test_extract1d.fits'
    write_onespec(out_wave, out_flux, out_ivar, raw_hdr,
                  spectrograph.name, str(output_file))
    assert output_file.exists(), f"OneSpec not created: {output_file}"

    one = OneSpec.from_file(str(output_file))
    np.testing.assert_allclose(one.wave, out_wave)
    np.testing.assert_allclose(one.flux, out_flux)
    np.testing.assert_allclose(one.ivar, out_ivar)
    assert one.PYP_SPEC == 'mmt_binospec_ifu'

    # Clean up
    output_file.unlink()


def test_binospec_ifu_sky_fiber_skymodel(redux_out):
    """Validate the sky-fibers-only joint sky model.

    The ``G270_Sky`` setup reduces the same data as ``G270`` but with
    ``[reduce][skysub][joint_fit_use_sci] = False``, so the joint sky model is
    built from the dedicated sky fibers alone rather than from all fibers.  This
    confirms that reduction completes and produces a sane sky model.
    """
    redux_path = Path(redux_out) / 'mmt_binospec_ifu' / 'G270_Sky'
    sci_dir = redux_path / 'Science'
    spec2d_files = sorted(sci_dir.glob('spec2d_*.fits'))
    assert len(spec2d_files) > 0, f"No spec2d files found in {sci_dir}"

    spec2d_file = spec2d_files[0]

    spec = load_spectrograph('mmt_binospec_ifu')
    detname = spec.get_det_name(1)
    spec2DObj = spec2dobj.Spec2DObj.from_file(str(spec2d_file), detname,
                                              chk_version=True)

    skymodel = spec2DObj.skymodel
    assert skymodel is not None, "SKYMODEL missing from spec2d"
    assert skymodel.shape == spec2DObj.sciimg.shape, \
        "SKYMODEL shape does not match science image"
    assert np.all(np.isfinite(skymodel)), "SKYMODEL has non-finite values"
    assert np.any(skymodel > 0), "SKYMODEL has no positive sky counts"


def _joint_sky_fiber_counts(redux_dir):
    """Parse the joint sky-fit fiber counts from a reduction's run log.

    ``FiberFindObjects._fit_fiber_sky_bspline`` logs one line per pass per
    detector recording how many sky and science fibers entered the joint
    bspline fit, e.g. ``... 1470817 pts from 40 sky + 320 sci fibers ...``.
    Returns the list of ``(n_sky, n_sci)`` tuples found in the ``run_pypeit``
    log under ``redux_dir`` (the ``*.log`` file written next to the ``.pypeit``
    file, excluding the dev-suite ``*.test.log`` wrapper log).
    """
    redux_dir = Path(redux_dir)
    logs = [p for p in redux_dir.glob('*.log')
            if not p.name.endswith('.test.log')]
    assert len(logs) == 1, \
        f"Expected exactly one run log in {redux_dir}, found {logs}"
    pattern = re.compile(r'(\d+)\s+sky\s+\+\s+(\d+)\s+sci\s+fibers')
    counts = []
    for line in logs[0].read_text().splitlines():
        if '_fit_fiber_sky_bspline' not in line:
            continue
        m = pattern.search(line)
        if m:
            counts.append((int(m.group(1)), int(m.group(2))))
    return counts


def test_binospec_ifu_joint_fit_fiber_selection(redux_out):
    """Confirm ``joint_fit_use_sci`` selects the intended fibers.

    ``G270`` reduces with the default ``joint_fit_use_sci = True`` so the joint
    sky bspline is fit over all sky *and* science fibers; ``G270_Sky`` reduces
    the *same* data with ``joint_fit_use_sci = False`` so only the dedicated
    sky fibers contribute.  This pins down that the parameter actually changes
    which fibers enter the fit, reading the per-pass provenance the sky fit
    logs for each detector.
    """
    base = Path(redux_out) / 'mmt_binospec_ifu'

    all_counts = _joint_sky_fiber_counts(base / 'G270')
    sky_counts = _joint_sky_fiber_counts(base / 'G270_Sky')

    assert all_counts, "No joint sky-fit log lines found for G270"
    assert sky_counts, "No joint sky-fit log lines found for G270_Sky"

    # Default reduction: both sky and science fibers participate.
    assert all(n_sky > 0 for n_sky, _ in all_counts), \
        f"G270 joint fit used no sky fibers: {all_counts}"
    assert all(n_sci > 0 for _, n_sci in all_counts), \
        f"G270 joint fit used no science fibers despite " \
        f"joint_fit_use_sci=True: {all_counts}"

    # Sky-only reduction: dedicated sky fibers only, zero science fibers.
    assert all(n_sky > 0 for n_sky, _ in sky_counts), \
        f"G270_Sky joint fit used no sky fibers: {sky_counts}"
    assert all(n_sci == 0 for _, n_sci in sky_counts), \
        f"G270_Sky joint fit used science fibers despite " \
        f"joint_fit_use_sci=False: {sky_counts}"
