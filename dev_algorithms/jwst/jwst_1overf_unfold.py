# modified from https://github.com/JWST-EREBUS/unfold_jwst/blob/master/unfold_jwst/background.py

import os, time, sys
import copy
import warnings
import scipy
import numpy as np
from scipy import ndimage
from scipy.stats import norm

import matplotlib.pyplot as plt

import multiprocessing

import astropy
from astropy.time import Time

from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.stats import biweight_location, biweight_midvariance
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.utils.exceptions import AstropyWarning

from photutils import detect_sources
from photutils import StdBackgroundRMS, MADStdBackgroundRMS, BiweightScaleBackgroundRMS
from photutils import Background2D, MeanBackground, MedianBackground, SExtractorBackground
from photutils import MMMBackground, BiweightLocationBackground, ModeEstimatorBackground

import jwst

warnings.simplefilter('ignore', category=AstropyWarning)

def fnoise_sub(data, bpm=None, error=None, namp=4, minimum_pixels=10, rej_nsigma=3, maxiters=5,
               brightstar_nsigma=2, npixels=4, fwhm=3, back_type='sextractor', back_rms_type='biweight',
               back_size=(51, 51), back_filter_size=(3, 3), back_rej_nsigma=3, back_maxiters=5,
               sub_bkg=True, mask_brightstar=True, evenOdd=True, skip_col=False, skip_row=False, inst='NIRCam'):
    """
    Subtract amplifier noise (i.e., the so called 1/f noise caused by readout) row-by-row
    This function was modified by Feige Wang and Jinyi Yang based on rowamp_sub.py in tshirt package (https://github.com/eas342/tshirt)
    We add sigma clipping and more proper object mask and bad pixel mask. We estimated and subtracted a
    smooth background before subtracting 1/f noise to avoid the over/under-subtraction for pixels with different local
    background (i.e., close to bright star halo). The smooth background was added back to the data after 1/f noise
    subtraction. It supports stripe subtraction along both row and column.

    Parameters
    ----------
    data: numpy array
        The input image

    bpm: numpy array (bool) or None
        Bad pixel mask of your data. True means bad

    error: numpy array or None
        Error array of your input image

    namp: int
        Number of amplifier. 4 for NIRCam stripe mode and 1 is for 1 output amplifier

    objmask: numpy array (bool) or None
        Mask for the bright sources. True pixels will not be used for background estimation.

    minimum_pixels: int or float
        Minimum number of pixes must be in a row. It will interpolate the background of
        other amplifiers to estimate the background in this row if the number of good pixels
        is smaller than minimum_pixels.

    rej_nsigma: int or float
        The number of standard deviations to use for clipping limit.

    maxiters: int or float
        The maximum number of sigma-clipping iterations to perform or
        `None` to clip until convergence is achieved (i.e., iterate
        until the last iteration clips nothing). If convergence is
        achieved prior to ``maxiters`` iterations, the clipping
        iterations will stop. The default is 5.

    brightstar_nsigma: int or float
        The nsigma for identifying bright stars.

    npixels: int
        The number of connected pixels used for identifying bright stars.

    fwhm: int
        FWHM in unit of pixel used for build a Gaussian kernel

    back_type: str
        Background type

    back_rms_type: str
        Background rms type

    back_size: int or array_like (int)
        The box size along each axis.  If ``box_size`` is a scalar then
        a square box of size ``box_size`` will be used.  If ``box_size``
        has two elements, they should be in ``(ny, nx)`` order.  For
        best results, the box shape should be chosen such that the
        ``data`` are covered by an integer number of boxes in both
        dimensions.

    back_filter_size: int or array_like (int), optional
        The window size of the 2D median filter to apply to the
        low-resolution background map.  If ``filter_size`` is a scalar
        then a square box of size ``filter_size`` will be used.  If
        ``filter_size`` has two elements, they should be in ``(ny, nx)``
        order.  A filter size of ``1`` (or ``(1, 1)``) means no
        filtering.

    back_rej_nsigma: int or float
        The number of standard deviations to use when estimating background for both the lower
        and upper clipping limit.

    back_maxiters: int or float
        The maximum number of sigma-clipping iterations used for estimating background to perform or
        `None` to clip until convergence is achieved (i.e., iterate
        until the last iteration clips nothing). If convergence is
        achieved prior to ``maxiters`` iterations, the clipping
        iterations will stop. The default is 5.

    sub_bkg: bool
        Subtract a smooth background before subtracting the stripe.

    mask_brightstar: bool
        Mask brightstar before estimating background and 1/f noise?

    evenOdd: bool
        Remove the even and odd offsets before doing row-by-row medians?

    skip_col: bool
        skip stripe along column subtraction?

    inst: str
        Instrument name. If inst!='NIRCam', only namp=1 is supported.

    Returns
    -------
    outimg: numpy array with the same shape of data
        Corrected data

    modelimg: numpy array with the same shape of data
        model image
    """

    ## Make a good pixel mask for pixels that will be used for background estimation.
    gpm = np.isfinite(data) & (data != 0.)
    if bpm is not None:
        gpm &= np.invert(bpm)

    ## Estimate a smooth background and make a bright star mask
    if sub_bkg or mask_brightstar:
        objmask, background_array, _ = get_objmask_bkg(data, bpm=np.invert(gpm), error=error,
                                                       brightstar_nsigma=brightstar_nsigma,
                                                       npixels=npixels, fwhm=fwhm, back_type=back_type,
                                                       back_rms_type=back_rms_type,
                                                       back_size=back_size, back_filter_size=back_filter_size,
                                                       back_rej_nsigma=back_rej_nsigma,
                                                       back_maxiters=back_maxiters)
        if mask_brightstar:
            gpm &= np.invert(objmask)
        if sub_bkg:
            data_masked = np.ma.masked_where(np.invert(gpm), data - background_array)
        else:
            data_masked = np.ma.masked_where(np.invert(gpm), data)
    else:
        data_masked = np.ma.masked_where(np.invert(gpm), data)

    data_masked.fill_value = np.nan

    if skip_row:
        modelimg = np.zeros_like(data)
        ## do the subtraction
        outimg = data - modelimg
    else:
        ## Get the model image
        modelimg = get_rowamp_model(data_masked, namp, minimum_pixels=minimum_pixels,
                                    rej_nsigma=rej_nsigma, maxiters=maxiters, evenOdd=evenOdd, inst=inst)

        modelimg[np.isnan(modelimg)] = 0.
        ## do the subtraction
        outimg = data - modelimg

    ## Also subtract stripe along column?
    # I set namp = 1 for the column stripe the readout is along row but not column.
    if not skip_col:
        modelimg_col = get_rowamp_model(data_masked.T - modelimg.T, 1, minimum_pixels=minimum_pixels,
                                        rej_nsigma=rej_nsigma, maxiters=maxiters, evenOdd=evenOdd, inst=inst)
        outimg -= modelimg_col.T
        modelimg += modelimg_col.T

    # I do not subtract anything for zero pixels.
    outimg[data == 0.] = 0.

    return outimg, modelimg

def get_objmask_bkg(data, bpm=None, error=None, brightstar_nsigma=2, npixels=4, fwhm=3, back_type='sextractor',
                    back_rms_type='biweight', back_size=(51, 51), back_filter_size=(3, 3), back_rej_nsigma=3,
                    back_maxiters=5, qa=None, show=False):
    """
    Derive object mask, background and background rms.

    Parameters
    ----------
    data
    bpm
    error
    brightstar_nsigma
    npixels
    fwhm : float or int
    back_type
    back_rms_type
    back_size
    back_filter_size
    back_rej_nsigma
    back_maxiters

    Returns
    -------

    """

    ## Build a Gaussian kernel
    sigma = fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()

    this_data = copy.deepcopy(data)
    if error is not None:
        this_error = copy.deepcopy(error)
    else:
        this_error = None

    if bpm is None:
        this_bpm_bkg = np.zeros_like(data, dtype='bool')
    else:
        this_bpm_bkg = copy.deepcopy(bpm)

    gpm = np.invert(copy.deepcopy(this_bpm_bkg))

    ## Estimate background and objmask with two iterations
    objmask, background_array = np.zeros_like(data, dtype='bool'), np.zeros_like(data)
    iiter = 0
    while iiter < 2:
        if iiter == 0:
            this_back_size = (201, 201)
        else:
            this_back_size = back_size
        # Background estimation
        background_array, background_rms = BKG2D(this_data, this_back_size, bpm=this_bpm_bkg,
                                                 filter_size=back_filter_size,
                                                 sigclip=back_rej_nsigma, back_type=back_type,
                                                 back_rms_type=back_rms_type,
                                                 back_maxiters=back_maxiters, verbose=False)
        # Mask bright star
        if this_error is None:
            threshold = background_array + brightstar_nsigma * background_rms
        else:
            threshold = background_array + brightstar_nsigma * this_error
        convolved_data = convolve(data * gpm, kernel)
        convolved_data[data == 0.] = 0.

        # Do the detection using Image Segmentation
        # The return is a Segmentation image
        # msgs.info('Making bright star mask with iter={:}'.format(niter+1))
        segm = detect_sources(convolved_data, threshold, npixels=npixels)

        # grow mask for bright stars
        bright_blob = ndimage.binary_erosion(segm.data > 0., iterations=11)
        mask_bright_blob = ndimage.binary_dilation(bright_blob, iterations=50)
        objmask = np.logical_or(mask_bright_blob, (segm.data > 0.))
        this_bpm_bkg &= objmask
        # replace bright star values with background values
        this_data[objmask] = background_array[objmask]

        if error is None:
            this_error = background_rms
        else:
            this_error = copy.deepcopy(error)
        iiter += 1

    if qa is not None or show:
        bsub = data - background_array
        gpm_stat = np.logical_and(np.invert(objmask), np.invert(this_bpm_bkg))
        data_stat = bsub[gpm_stat].flatten()
        error_stat = error[gpm_stat].flatten()
        # data_median, data_std = np.median(data_stat), np.std(data_stat)

        # Plot Chi distribution
        sig_range = 6.
        n_bins = 50
        binsize = 2.0 * sig_range / n_bins
        bins_histo = -sig_range + np.arange(n_bins) * binsize + binsize / 2.0
        xvals = np.arange(-10.0, 10, 0.02)
        gauss = norm(loc=0.0, scale=1.0)

        plt.hist(data_stat * utils.inverse(error_stat), bins_histo, density=True, histtype='step', align='mid',
                 color='k', linewidth=1)
        plt.plot(xvals, gauss.pdf(xvals), 'r-', lw=1)
        plt.xlim(-5.5, 5.5)
        plt.ylim(0, 0.55)
        plt.xlabel(r'$\chi$', fontsize=13, labelpad=0)
        plt.ylabel(r'Distribution', fontsize=13)
        if qa is not None:
            plt.savefig(qa, dpi=300)
        if show:
            plt.show()
        plt.close()

    return objmask, background_array, this_error


def get_rowamp_model(data_masked, namp, minimum_pixels=10, rej_nsigma=3, maxiters=5, evenOdd=True, inst='NIRCam'):
    """
    Get a 1/f noise model
    """

    slowread_model = np.zeros_like(data_masked.data)
    fastread_model = np.zeros_like(data_masked.data)

    if namp == 4:
        ## mask to keep track of which rows have enough pixels to use
        amp_badrow_mask = np.zeros_like(data_masked.data, dtype=bool)

        ## list where the amplifiers are
        if inst == 'NIRCam':
            ampStarts = [0, 512, 1024, 1536]
            ampEnds = [512, 1024, 1536, 2048]
            ampWidth = 512
        else:
            #msgs.error('Only NIRCam is supported at this moment.')
            print('Only NIRCam is supported at this moment.')

        ## loop through amps
        for amp in np.arange(4):
            this_amp_tmp = data_masked[:, ampStarts[amp]:ampEnds[amp]]
            if evenOdd == True:
                thisAmp, even_odd_model = do_even_odd(this_amp_tmp, rej_nsigma=rej_nsigma, maxiters=maxiters)
                slowread_model[:, ampStarts[amp]:ampEnds[amp]] = even_odd_model
            else:
                ## even if not doing an even/odd correction, still do an overall median
                _, this_median, _ = sigma_clipped_stats(this_amp_tmp, sigma=rej_nsigma, maxiters=maxiters,
                                                        cenfunc='median', stdfunc='std')
                slowread_model[:, ampStarts[amp]:ampEnds[amp]] = np.nanmedian(this_median)
                thisAmp = this_amp_tmp - slowread_model[:, ampStarts[amp]:ampEnds[amp]]

            bad_rows = np.sum(np.isfinite(thisAmp.filled()), axis=1) <= minimum_pixels
            _, medVals, _ = sigma_clipped_stats(thisAmp, axis=1, sigma=rej_nsigma, maxiters=maxiters, cenfunc='median',
                                                stdfunc='mad_std')
            medVals[bad_rows] = 0.
            ## tile this to make a model across the fast-read direction
            fastread_model[:, ampStarts[amp]:ampEnds[amp]] = np.tile(medVals, [ampWidth, 1]).T
            ## check if the mask leaves too few pixels in some of the rows in this amplifier
            amp_badrow_mask[bad_rows, ampStarts[amp]:ampEnds[amp]] = True

        ## Let's replace the bad rows with that row in other amplifiers
        replace_rows = (np.sum(amp_badrow_mask, axis=1) > 0) & (
                np.sum(amp_badrow_mask, axis=1) < amp_badrow_mask.shape[1])
        if np.sum(replace_rows) > 0:
            _, rowModel, _ = sigma_clipped_stats(fastread_model, axis=1, sigma=rej_nsigma, maxiters=maxiters,
                                                 cenfunc='median', stdfunc='std')
            for amp in np.arange(4):
                this_amp_replace = replace_rows & (np.sum(amp_badrow_mask[:, ampStarts[amp]:ampEnds[amp]], axis=1) > 0)
                if np.sum(this_amp_replace) > 0:
                    this_shape = ampEnds[amp] - ampStarts[amp]
                    fastread_model[this_amp_replace, ampStarts[amp]:ampEnds[amp]] = np.tile(rowModel[this_amp_replace],
                                                                                            [this_shape, 1]).T

    elif namp == 1:
        if evenOdd == True:
            thisAmp, slowread_model = do_even_odd(data_masked, rej_nsigma=rej_nsigma, maxiters=maxiters)
        else:
            _, slowread_model, _ = sigma_clipped_stats(data_masked, sigma=rej_nsigma, maxiters=maxiters,
                                                       cenfunc='median', stdfunc='std')
            thisAmp = data_masked - slowread_model

        _, medVals, _ = sigma_clipped_stats(thisAmp, axis=1, sigma=rej_nsigma, maxiters=maxiters, cenfunc='median',
                                            stdfunc='std')
        ## tile this to make a constant model
        tiled_med = np.tile(medVals, [data_masked.data.shape[1], 1]).T
        ## put the results in the model image
        fastread_model[:, :] = tiled_med
    else:
        #msgs.error('{:} amplifiers is not implemented yet.'.format(namp))
        print('{:} amplifiers is not implemented yet.'.format(namp))

    ## put the results in the model image
    modelimg = slowread_model + fastread_model

    return modelimg

def BKG2D(data, back_size, bpm=None, filter_size=(3, 3), sigclip=5, back_type='sextractor', back_rms_type='biweight',
          back_maxiters=5, verbose=True):
    """
    Estimate background for a given data

    Parameters
    ----------
    data:
    back_size:
    bpm:
    filter_size:
    sigclip:
    back_type:
    back_rms_type:
    back_maxiters:
    verbose:

    Returns
    -------
    bkg_map
    rms_map
    """

    ## Sky background estimation
    if 'median' in back_type.lower():
        bkg_estimator = MedianBackground()
    elif back_type.lower() == 'mean':
        bkg_estimator = MeanBackground()
    elif back_type.lower() == 'sextractor':
        bkg_estimator = SExtractorBackground()
    elif back_type.lower() == 'mmm':
        bkg_estimator = MMMBackground()
    elif back_type.lower() == 'biweight':
        bkg_estimator = BiweightLocationBackground()
    elif back_type.lower() == 'mode':
        bkg_estimator = ModeEstimatorBackground()
    else:
        msgs.warn('{:} Background is not found, using MedianBackground Instead.'.format(back_type))
        back_type = 'median'
        bkg_estimator = MedianBackground()

    if back_rms_type.lower() == 'std':
        bkgrms_estimator = StdBackgroundRMS()
    elif back_rms_type.lower() == 'mad':
        bkgrms_estimator = MADStdBackgroundRMS()
    elif back_rms_type.lower() == 'biweight':
        bkgrms_estimator = BiweightScaleBackgroundRMS()
    else:
        msgs.warn('{:} Background RMS type is not found, using STD Instead.'.format(back_rms_type))
        bkgrms_estimator = StdBackgroundRMS()

    if verbose:
        msgs.info('Estimating {:} BACKGROUND with Photutils Background2D.'.format(back_type))

    tmp = copy.deepcopy(data)
    Sigma_Clip = SigmaClip(sigma=sigclip, maxiters=back_maxiters)
    bkg = Background2D(tmp, back_size, mask=bpm, filter_size=filter_size, sigma_clip=Sigma_Clip,
                       bkg_estimator=bkg_estimator, bkgrms_estimator=bkgrms_estimator)
    bkg_map, rms_map = bkg.background, bkg.background_rms
    bkg_map[tmp == 0.] = 0.

    if back_type == 'GlobalMedian':
        bkg_map = np.ones_like(bkg_map) * np.nanmedian(bkg_map[np.invert(bpm)])

    return bkg_map, rms_map

def do_even_odd(this_amp, rej_nsigma=3, maxiters=10):
    """
    Do an even-odd correction for a given amplifier
    If only one amplifier is used, it can be the whole image
    This function was modified by Feige Wang based on rowamp_sub.py in tshirt package (https://github.com/eas342/tshirt)
    """

    even_odd_model = np.zeros_like(this_amp.data)
    _, even_offset, _ = sigma_clipped_stats(this_amp[:, 0::2], sigma=rej_nsigma, maxiters=maxiters, cenfunc='median',
                                            stdfunc='std')
    _, odd_offset, _ = sigma_clipped_stats(this_amp[:, 1::2], sigma=rej_nsigma, maxiters=maxiters, cenfunc='median',
                                           stdfunc='std')
    even_odd_model[:, 0::2] = even_offset
    even_odd_model[:, 1::2] = odd_offset

    return this_amp - even_odd_model, even_odd_model
