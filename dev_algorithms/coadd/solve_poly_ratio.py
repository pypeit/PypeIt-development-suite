

import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from pypeit import utils
from pypeit.core import coadd
from scipy import interpolate
from coadd_new import interp_spec
from scipy import stats
from pypeit import msgs
import IPython

#def solve_poly_fn(theta, xvector, polyfunc, nback = None):
#
#    if nback is None:
#        acoeff = theta
#    else:
#        acoeff = theta[0:-nback]
#        bcoeff = theta[-nback:]
#
#    ymult = utils.func_val(acoeff, xvector, polyfunc, minx=wave_min, maxx=wave_max)
#    if nback is not None:
#        ymult = utils.func_val(acoeff, xvector, polyfunc, minx=wave_min, maxx=wave_max)


def error_renormalize(data, model, ivar, mask=None, clip = 6.0, max_corr = 5.0):

    if mask is None:
        mask = (ivar > 0)

    chi2 = ivar*(data-model)**2
    igood = (chi2 < clip**2) & mask
    if (np.sum(igood) > 0):
        chi2_good = chi2[igood]
        gauss_prob = 1.0 - 2.0 * stats.norm.cdf(-1.0)
        chi2_sigrej = np.percentile(chi2_good[igood], 100.0*gauss_prob)
        sigma_corr = np.sqrt(chi2_sigrej)
        if sigma_corr < 1.0:
            msgs.warn("Error renormalization found correction factor sigma_corr = {:f} < 1." + msgs.newline() +
                      " Errors are overestimated so not applying correction".format(sigma_corr))
            sigma_corr = 1.0
        if sigma_corr > max_corr:
            msgs.warn("Error renormalization found sigma_corr/sigma = {:f} > {:f}." + msgs.newline() +
                      "Errors are severely underestimated." + msgs.newline() +
                      "Setting correction to sigma_corr = {:f}".format(sigma_corr, max_corr, max_corr))
            sigma_corr = max_corr

        return sigma_corr

def poly_ratio_fitfunc_chi2(theta, flux_ref, thismask, arg_dict):
    """
    Function to be optimized for poly_ratio rescaling

    Args:
        theta:
        flux_ref:
        ivar_ref:
        thismask:
        arg_dict:

    Returns:

    """

    # Unpack the data to be rescaled, the mask for the reference spectrum, and the wavelengths
    flux = arg_dict['flux']
    ivar = arg_dict['ivar']
    mask = arg_dict['mask']
    ivar_ref = arg_dict['ivar_ref']
    wave = arg_dict['wave']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    func = arg_dict['func']
    # Evaluate the polynomial for rescaling
    ymult = utils.func_val(theta, wave, func, minx=wave_min, maxx=wave_max)
    flux_scale = ymult*flux
    mask_both = mask & thismask
    # This is the formally correct ivar
    #totvar = utils.calc_ivar(ivar_ref) + ymult**2 * utils.calc_ivar(ivar)
    #ivartot = mask_both*utils.calc_ivar(totvar)

    # The errors are rescaled at every function evaluation, but we only allow the errors to get smaller by up to a
    # factor of 1e4, and we only allow them to get larger slowly (as the square root).  This should very strongly
    # constrain the flux-corrrection vectors from going too small (or negative), or too large.
    ## Schlegel's version here, this will not work with robust_optimize since the chi^2 is not a real chi^2
    vmult = np.fmax(ymult,1e-4)*(ymult <= 1.0) + np.sqrt(ymult)*(ymult > 1.0)
    ivartot = mask_both/(1.0/(ivar + np.invert(mask_both)) + np.square(vmult)/(ivar_ref + np.invert(mask_both)))
    chi_vec = mask_both * (flux_ref - flux_scale) * np.sqrt(ivartot)
    chi2 = np.sum(np.square(chi_vec))
    return chi2

def poly_ratio_fitfunc(flux_ref, thismask, arg_dict, **kwargs_opt):

    # flux_ref, ivar_ref act like the 'data', the rescaled flux will be the 'model'

    # Function that we are optimizing
    #result = scipy.optimize.differential_evolution(poly_ratio_fitfunc_chi2, args=(flux_ref, ivar_ref, thismask, arg_dict,), **kwargs_opt)
    guess = arg_dict['guess']
    result = scipy.optimize.minimize(poly_ratio_fitfunc_chi2, guess, args=(flux_ref, thismask, arg_dict),  **kwargs_opt)
    flux = arg_dict['flux']
    ivar = arg_dict['ivar']
    mask = arg_dict['mask']
    ivar_ref = arg_dict['ivar_ref']
    wave = arg_dict['wave']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    func = arg_dict['func']
    # Evaluate the polynomial for rescaling
    ymult = utils.func_val(result.x, wave, func, minx=wave_min, maxx=wave_max)
    flux_scale = ymult*flux
    mask_both = mask & thismask
    totvar = utils.calc_ivar(ivar_ref) + ymult**2*utils.calc_ivar(ivar)
    ivartot = mask_both*utils.calc_ivar(totvar)
    # Now rescale the errors

    return result, flux_scale, ivartot




from astropy.table import Table
from astropy.io import fits
from matplotlib import pyplot as plt
import sys

fitfunc = poly_ratio_fitfunc
stackfile  = '/Users/joe/REDUX/lris_redux/Nov_2004/Final/SDSSJ073522.43+295710.1_N.fits'
hdu = fits.open(stackfile)
flux_ref = hdu[0].data
sig_ref = hdu[1].data
ivar_ref = utils.calc_ivar(sig_ref)
mask_ref = (ivar_ref > 0.0)
wave_ref = hdu[2].data

infiles = ['/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0063.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0064.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0065.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0066.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0219.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0220.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0221.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0222.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0063.fits.gz']
nfiles = len(infiles)
objid = 0
for idx, file in enumerate(infiles):
    obj = Table.read(file, hdu=5)
    flux = obj[objid]['FLUX_OPT']
    ivar = obj[objid]['IVAR_OPT']
    wave = obj[objid]['WAVE_OPT']
    if idx == 0:
        nspec = flux.size
        flux_arr = np.zeros((nfiles, nspec))
        wave_arr = np.zeros((nfiles, nspec))
        ivar_arr = np.zeros((nfiles, nspec))
    flux_arr[idx,:] = flux
    ivar_arr[idx,:] = ivar
    wave_arr[idx,:] = wave

mask_arr = (ivar_arr > 0.0)
#plt.plot(wave_ref, flux_ref, drawstyle='steps-mid')
#for ii in range(nfiles):
#    plt.plot(wave_arr[:,ii], flux_arr[:,ii],drawstyle='steps-mid')

flux_inter, ivar_inter, mask_inter = interp_spec(wave_ref, wave_arr, flux_arr, ivar_arr, mask_arr)

idx = 5
wave = wave_ref
flux = flux_inter[idx, :]
ivar = ivar_inter[idx, :]
mask = ivar > 0


#if norder is None:
#    # If user did not specify an order to use, determine it automatically based on S/N
#    rms_snr, _ = coadd.sn_weights(wave, flux, ivar, mask=mask)
#    rms_snr_ref, _ = coadd.sn_weights(wave, flux_ref, ivar_ref, mask=mask_ref)
#sys.exit(-1)

#def solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, norder, mask = None, mask_ref = None, func='legendre',
#                     maxiter=3, sticky=True, use_mad=True,
#                     lower=3.0, upper=3.0, min_good =0.05):
func ='legendre'
norder = 3
maxiter = 3
sticky = True
lower = 3.0 #stats.norm.cdf(-3.0)
upper = 3.0 #stats.norm.cdf(3.0)
min_good = 0.05
use_mad=True
debug=True

if mask is None:
    mask = (ivar > 0.0)
if mask_ref is None:
    mask_ref = (ivar_ref > 0.0)

nspec = wave.size
# Determine an initial guess
if ((np.sum(mask) > min_good*nspec) & (np.sum(mask_ref) > min_good*nspec)):
    flux_25 = np.percentile(flux[mask],25)
    flux_ref_25 = np.percentile(flux_ref[mask_ref],25)
    ratio = flux_ref_25/flux_25
else:
    ratio = 1.0
guess = np.append(ratio, np.zeros(norder-1))
arg_dict = dict(flux = flux, ivar = ivar, mask = mask,
                ivar_ref = ivar_ref,
                wave = wave, wave_min = wave.min(),
                wave_max = wave.max(), func = func, norder = norder, guess = guess)

result, ymodel, ivartot, outmask = utils.robust_optimize(flux_ref, poly_ratio_fitfunc, arg_dict,
                                                         inmask=mask_ref, use_mad=False,
                                                         maxiter=maxiter, lower=lower, upper=upper, sticky=sticky)

if debug:
    plt.plot(wave,flux_ref, color='black', drawstyle='steps-mid', zorder=3, label='Reference spectrum')
    plt.plot(wave,flux, color='dodgerblue', drawstyle='steps-mid', zorder = 10, alpha = 0.5, label='Original spectrum')
    plt.plot(wave,ymodel, color='red', drawstyle='steps-mid', alpha=0.7, zorder=1, linewidth=2, label='Rescaled spectrum')
    plt.legend()
    plt.show()
