# selected scripts from:
# https://github.com/gbrammer/grizli/blob/master/grizli/utils.py
# https://github.com/gbrammer/grizli/blob/master/grizli/jwst_utils.py

import numpy as np

def unset_dq_bits(value, okbits=32+64+512, verbose=False):
    """
    Unset bit flags from a DQ array
    For WFC3/IR, the following DQ bits can usually be unset:
    32, 64: these pixels usually seem OK
       512: blobs not relevant for grism exposures
    Parameters
    ----------
    value : int, `~numpy.ndarray`
        Input DQ value
    okbits : int
        Bits to unset
    verbose : bool
        Print some information
    Returns
    -------
    new_value : int, `~numpy.ndarray`
    """
    bin_bits = np.binary_repr(okbits)
    n = len(bin_bits)
    for i in range(n):
        if bin_bits[-(i+1)] == '1':
            if verbose:
                print(2**i)

            value -= (value & 2**i)

    return value

def exposure_oneoverf_correction(file, axis=None, thresholds=[5, 4, 3], erode_mask=None, dilate_iterations=3,
                                 deg_pix=64, make_plot=True, init_model=0, in_place=False, skip_miri=True, verbose=True,
                                 **kwargs):
    """
    1/f correction for individual exposure

    1. Create a "background" mask with `sep`
    2. Identify sources above threshold limit in the background-subtracted
       image
    3. Iterate a row/column correction on threshold-masked images.  A
       chebyshev polynomial is fit to the correction array to try to isolate
       just the high-frequency oscillations.

    Parameters
    ----------
    file : str
        JWST raw image filename

    axis : int
        Axis over which to calculated the correction. If `None`, then defaults
        to ``axis=1`` (rows) for NIRCam and ``axis=1`` (columns) for NIRISS.

    thresholds : list
        List of source identification thresholds

    erode_mask : bool
        Erode the source mask to try to remove individual pixels that satisfy
        the S/N threshold.  If `None`, then set to False if the exposure is a
        NIRISS dispersed image to avoid clipping compact high-order spectra
        from the mask and True otherwise (for NIRISS imaging and NIRCam
        generally).

    dilate_iterations : int
        Number of `binary_dilation` iterations of the source mask

    deg_pix : int
        Scale in pixels for each degree of the smooth chebyshev polynomial

    make_plot : bool
        Make a diagnostic plot

    init_model : scalar, array-like
        Initial correction model, e.g., for doing both axes

    in_place : bool
        If True, remove the model from the 'SCI' extension of ``file``

    skip_miri : bool
        Don't run on MIRI exposures

    verbose : bool
        Print status messages

    Returns
    -------
    fig : `~matplotlib.figure.Figure`, None
        Diagnostic figure if `make_plot=True`

    model : array-like
        The row- or column-average correction array
    """
    import numpy as np
    import scipy.ndimage as nd
    import matplotlib.pyplot as plt

    import astropy.io.fits as pyfits
    import sep

    im = pyfits.open(file)
    if (im[0].header['INSTRUME'] in 'MIRI') & (skip_miri):
        im.close()

        msg = 'exposure_oneoverf_correction: Skip for MIRI'
        #utils.log_comment(utils.LOGFILE, msg, verbose=verbose)

        return None, 0

    if axis is None:
        if im[0].header['INSTRUME'] in ('NIRISS', 'NIRSPEC'):
            axis = 0
        else:
            axis = 1

    elif axis < 0:
        # Opposite axis
        if im[0].header['INSTRUME'] in ('NIRISS', 'NIRSPEC'):
            axis = 1
        else:
            axis = 0

    msg = f'exposure_oneoverf_correction: {file} axis={axis} deg_pix={deg_pix}'
    #utils.log_comment(utils.LOGFILE, msg, verbose=verbose)

    if im[0].header['INSTRUME'] in ('NIRSPEC'):
        erode_mask = False

    if erode_mask is None:
        if im[0].header['FILTER'].startswith('GR150'):
            erode_mask = False
        elif im[0].header['PUPIL'].startswith('GRISM'):
            erode_mask = False
        else:
            erode_mask = True

    dq = unset_dq_bits(im['DQ'].data, 4)
    dqmask = dq == 0
    mask = dqmask

    err = im['ERR'].data
    dqmask &= (err > 0) & np.isfinite(err)

    sci = im['SCI'].data.astype(np.float32) - init_model
    if deg_pix == 0:
        bw = sci.shape[0] // 64
    else:
        bw = deg_pix

    bkg = sep.Background(sci, mask=~dqmask, bw=bw, bh=bw)
    back = bkg.back()

    sn_mask = (sci - back) / err > thresholds[0]
    if erode_mask:
        sn_mask = nd.binary_erosion(sn_mask)

    sn_mask = nd.binary_dilation(sn_mask, iterations=dilate_iterations)

    mask = dqmask & ~sn_mask

    cheb = 0

    if make_plot:
        fig, ax = plt.subplots(1, 1, figsize=(6, 3))
    else:
        fig = None

    for _iter, thresh in enumerate(thresholds):

        if deg_pix == 0:
            sci = im['SCI'].data * 1. - init_model
        else:
            sci = im['SCI'].data * 1. - back - init_model

        sci[~mask] = np.nan
        med = np.nanmedian(sci, axis=axis)

        if axis == 0:
            model = np.zeros_like(sci) + (med - cheb)
        else:
            model = (np.zeros_like(sci) + (med - cheb)).T

        sn_mask = ((sci - model) / err > thresh) & dqmask
        if erode_mask:
            sn_mask = nd.binary_erosion(sn_mask)

        sn_mask = nd.binary_dilation(sn_mask, iterations=dilate_iterations)
        mask = dqmask & ~sn_mask
        mask &= (sci - model) / err > -thresh

        if make_plot:
            ax.plot(med, alpha=0.5)

    nx = med.size
    xarr = np.linspace(-1, 1, nx)
    ok = np.isfinite(med)

    if deg_pix == 0:
        # Don't remove anything from the median profile
        cheb = 0.
        deg = -1

    elif deg_pix >= nx:
        # Remove constant component
        cheb = np.nanmedian(med)
        deg = 0

    else:
        # Remove smooth component
        deg = nx // deg_pix

        for _iter in range(3):
            coeffs = np.polynomial.chebyshev.chebfit(xarr[ok], med[ok], deg)
            cheb = np.polynomial.chebyshev.chebval(xarr, coeffs)
            ok = np.isfinite(med) & (np.abs(med - cheb) < 0.05)

        if make_plot:
            ax.plot(np.arange(nx)[ok], cheb[ok], color='r')

    if axis == 0:
        model = np.zeros_like(sci) + (med - cheb)
    else:
        model = (np.zeros_like(sci) + (med - cheb)).T

    if make_plot:
        ax.set_title(f'{file} axis={axis}')
        ax.grid()

        fig.tight_layout(pad=0)

    data = im['SCI'].data
    im.close()

    if in_place:
        msg = f'exposure_oneoverf_correction: {file} apply to file'
        #utils.log_comment(utils.LOGFILE, msg, verbose=verbose)

        with pyfits.open(file, mode='update') as im:
            im[0].header['ONEFEXP'] = True, 'Exposure 1/f correction applied'
            im[0].header['ONEFAXIS'] = axis, 'Axis for 1/f correction'
            im[0].header['ONEFDEG'] = deg, 'Degree of smooth component'
            im[0].header['ONEFNPIX'] = deg_pix, 'Pixels per smooth degree'

            model[~np.isfinite(model)] = 0

            im['SCI'].data -= model
            im.flush()

        if make_plot:
            fig.savefig(file.split('.fits')[0] + f'_onef_axis{axis}.png')
            plt.close('all')

    return fig, model, data