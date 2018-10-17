


""" Module for basic utilties with holy grail
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import numba as nb

from scipy.ndimage.filters import gaussian_filter
from scipy.signal import resample
import scipy
from scipy.optimize import curve_fit
from astropy.io import fits
from matplotlib import pyplot as plt
from pypeit import msgs
from pypeit.core import arc




from pypeit.core import arc


def shift_and_stretch_old(spec, shift, stretch):

    x1 = np.arange(spec.shape[0])
    spec_resamp = scipy.signal.resample(spec, int(spec.size + stretch))
    x2 = np.arange(spec_resamp.shape[0]) + shift
    spec_out = np.interp(x1, x2, spec_resamp,left = 0.0,right=0.0)
    return spec_out


def shift_and_stretch(spec, shift, stretch):

    # Positive value of shift means features shift to larger pixel values

    nspec = spec.shape[0]
    # pad the spectrum on both sizes
    x1 = np.arange(nspec)/float(nspec)
    nspec_stretch = int(nspec*stretch)
    x2 = np.arange(nspec_stretch)/float(nspec_stretch)
    spec_str = (scipy.interpolate.interp1d(x1, spec, kind = 'quadratic', bounds_error = False, fill_value = 0.0))(x2)
    # Now create a shifted version
    ind_shift = np.arange(nspec_stretch) - shift
    spec_str_shf = (scipy.interpolate.interp1d(np.arange(nspec_stretch), spec_str, kind = 'quadratic', bounds_error = False, fill_value = 0.0))(ind_shift)
    # Now interpolate onto the original grid
    spec_out = (scipy.interpolate.interp1d(np.arange(nspec_stretch), spec_str_shf, kind = 'quadratic', bounds_error = False, fill_value = 0.0))(np.arange(nspec))

    return spec_out


# JFH I think this should be done with scipy.optimize to find the maximum value of the cc correlation as a function of
# shift and stretch, rather than with curve_fit
def xcorr_shift_stretch(theta, y1, y2):

    #shift = theta[0]
    #stretch = theta[1]
    y2_corr = shift_and_stretch(y2, theta[0], theta[1])
    # Zero lag correlation
    corr_zero = np.sum(y1*y2_corr)
    corr_denom = np.sqrt(np.sum(y1*y1)*np.sum(y2*y2))
    corr_norm = corr_zero/corr_denom
    #corr = scipy.signal.correlate(y1, y2_corr, mode='same')
    return -corr_norm


# This function only applies a stretch
def xcorr_stretch(theta, y1, y2):

    #shift = theta[0]
    #stretch = theta[1]
    y2_corr = shift_and_stretch(y2, 0.0, theta[0])
    # Zero lag correlation
    corr_zero = np.sum(y1*y2_corr)
    corr_denom = np.sqrt(np.sum(y1*y1)*np.sum(y2*y2))
    corr_norm = corr_zero/corr_denom
    return -corr_norm


def xcross_shift(inspec1,inspec2, smooth = 5.0):

    nspec = inspec1.shape[0]
    y1 = scipy.ndimage.filters.gaussian_filter(inspec1, smooth)
    y2 = scipy.ndimage.filters.gaussian_filter(inspec2, smooth)
    lags = np.arange(-nspec + 1, nspec)
    corr = scipy.signal.correlate(y1, y2, mode='full')
    corr_denom = np.sqrt(np.sum(y1*y1)*np.sum(y2*y2))
    corr_norm = corr/corr_denom
    pix_max = corr_norm.argmax()
    corr_max_fit, lag_max, width, cent_err = arc.fit_arcspec(lags, corr_norm, np.array([pix_max]), 7)

    return lag_max[0], corr_max_fit[0]

def xcross_shift_nosm(y1,y2, debug = False):

    nspec = y1.shape[0]
    lags = np.arange(-nspec + 1, nspec)
    corr = scipy.signal.correlate(y1, y2, mode='full')
    corr_denom = np.sqrt(np.sum(y1*y1)*np.sum(y2*y2))
    corr_norm = corr/corr_denom
    output = arc.detect_lines(corr_norm, nfitpix=7, sigdetect=5.0, FWHM=20.0, cont_samp=30, nfind = 1)
    pix_max = output[1]
    corr_max = np.interp(pix_max, np.arange(lags.shape[0]),corr_norm)
    lag_max  = np.interp(pix_max, np.arange(lags.shape[0]),lags)
    #pix_max = corr_norm.argmax()
    #corr_max_fit, lag_max, width, cent_err = arc.fit_arcspec(lags, corr_norm, np.array([pix_max]), 7)
    if debug:
        # Interpolate for bad lines since the fitting code often returns nan
        plt.figure(figsize=(14, 6))
        plt.plot(lags, corr_norm, color='black', drawstyle = 'steps-mid', lw=3, label = 'x-corr', linewidth = 1.0)
        plt.plot(lag_max[0], corr_max[0],'g+', markersize =6.0, label = 'peak')
        plt.title('Best shift = {:5.3f}'.format(lag_max[0]) + ',  corr_max = {:5.3f}'.format(corr_max[0]))
        plt.legend()
        plt.show()

    return lag_max[0], corr_max[0]


def fit_shift_stretch(inspec1, inspec2, smooth = 3.0, shift_mnmx = (0.8,1.2), stretch_mnmx = (0.8,1.2), debug = True):

    nspec = inspec1.size
    y1 = scipy.ndimage.filters.gaussian_filter(inspec1, smooth)
    y2 = scipy.ndimage.filters.gaussian_filter(inspec2, smooth)

    # Do the cross-correlation first and determine the
    shift_cc, cc_val = xcross_shift_nosm(y1, y2, debug = debug)

    bounds = [(shift_mnmx[0]*shift_cc,shift_mnmx[1]*shift_cc), stretch_mnmx]
    guess = np.array([shift_cc,1.0])
    #result = scipy.optimize.minimize(xcorr_shift_stretch, guess, args=(y1,y2), bounds=bounds)
    result = scipy.optimize.differential_evolution(xcorr_shift_stretch, args=(y1,y2), tol = 1e-4,
                                                   bounds=bounds, disp=False, polish=True) # ToDO Should we polish?

    #corr_val = xcorr_stretch((shift,stretch), y1,y2)
    #corr = scipy.signal.correlate(y1, y2, mode='same')

    if not result.success:
        msgs.warn('Fit for shift and stretch did not converge!')


    if debug:
        x1 = np.arange(nspec)
        inspec2_trans = shift_and_stretch(inspec2, result.x[0], result.x[1])
        plt.plot(x1,inspec1, 'k-', drawstyle='steps')
        plt.plot(x1,inspec2_trans, 'r-', drawstyle='steps')
        plt.show()


    return result.success, result.x[0], result.x[1], -result.fun, shift_cc, cc_val


def fit_shift_stretch_iter(inspec1, inspec2, smooth = 5.0, shift_mnmx = (-1.0,1.0), stretch_mnmx = (0.9,1.1), debug = True):



    shift_tot = 0.0
    stretch_tot = 1.0
    niter = 10

    stretch_vec = np.linspace(stretch_mnmx[0],stretch_mnmx[1],endpoint=True, num = 100)
    nspec = inspec1.size
    y1 = scipy.ndimage.filters.gaussian_filter(inspec1, smooth)
    y2 = scipy.ndimage.filters.gaussian_filter(inspec2, smooth)
    bounds = [stretch_mnmx]
    guess = np.array([1.0])

    for iter in range(niter):
        this_shift, corr_val = xcross_shift_nosm(y1,y2)
        shift_tot += shift
        y2_trans = shift_and_stretch(y2,shift,1.0)
        corr_vec = np.zeros_like(stretch_vec)
        for ii in np.arange(stretch_vec.shape[0]):
            corr_vec[ii] = xcorr_stretch([stretch_vec[ii]],y1,y2_trans)
        pix_max = corr_vec.argmax()
        corr_max_fit, stretch_max, width, cent_err = arc.fit_arcspec(stretch_vec, -corr_vec, np.array([pix_max]), 7)

        #result = scipy.optimize.minimize(xcorr_shift, guess, args=(y1,y2), bounds=bounds)
        #result = scipy.optimize.differential_evolution(xcorr_shift_stretch, args=(y1,y2), bounds=bounds, popsize=25, recombination=0.7,
        #                                           disp=False, polish=True) # ToDO Should we polish?

    #corr_val = xcorr_stretch((shift,stretch), y1,y2)
    #corr = scipy.signal.correlate(y1, y2, mode='same')

    if not result.success:
        msgs.warn('Fit for shift and stretch did not converge!')


    if debug:
        x1 = np.arange(nspec)
        inspec2_trans = shift_and_stretch(inspec2, result.x[0], result.x[1])
        plt.plot(x1,inspec1, 'k-', drawstyle='steps')
        plt.plot(x1,inspec2_trans, 'r-', drawstyle='steps')
        plt.show()


    return result.success, result.x[0], result.x[1], -result.fun





hdu = fits.open('./spec_array.fits')
spec = hdu[0].data.astype('float64')
inspec1=spec[:,0]
nslits = spec.shape[1]
nspec = spec.shape[0]
#inspec2 = spec[:,4]

# TESTING
#shift = 200.0
#stretch = 0.9
#smooth = 5.0

# Apply the shift and stretch to a copy of inspec1
#inspec2 = shift_and_stretch(inspec1, shift, stretch)
# evaluate the correlation coefficient at this value
success = np.zeros(nslits,dtype=bool)
shift = np.zeros(nslits)
shift_cc = np.zeros(nslits)
stretch = np.zeros(nslits)
xcorr = np.zeros(nslits)
xcorr_cc = np.zeros(nslits)
smooth = 5.0

y1 = scipy.ndimage.filters.gaussian_filter(inspec1, smooth)

for islit in range(nslits):
    output = fit_shift_stretch(inspec1,spec[:,islit], debug = True)
    success[islit], shift[islit], stretch[islit], xcorr[islit], shift_cc[islit], xcorr_cc[islit] = output

