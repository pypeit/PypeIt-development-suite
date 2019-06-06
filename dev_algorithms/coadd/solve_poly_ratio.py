

import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from pypeit import utils
from pypeit.core import coadd
from scipy import interpolate

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



def interp_spec(wave_ref, waves, fluxes, ivars, masks):
    '''
    Interpolate all spectra to the page of the iref spectrum
    Args:
        waves:
        fluxes:
        sigs:
        masks:
        iref:
    Returns:
    '''

    masks_float = np.zeros_like(fluxes)
    masks_float[masks] = 1.0

    #wave_iref,flux_iref,ivar_iref = waves[iref,:], fluxes[iref,:],ivars[iref,:]
    #sig_iref = np.sqrt(utils.calc_ivar(ivar_iref))

    # Interpolate spectra to have the same wave grid with the iexp spectrum.
    # And scale spectra to the same flux level with the iexp spectrum.
    fluxes_inter = np.zeros_like(fluxes)
    ivars_inter = np.zeros_like(ivars)
    masks_inter = np.zeros_like(masks)

    nexp = np.shape(fluxes)[0]
    for ii in range(nexp):

        mask_ii = masks[ii,:]

        if np.sum(wave_ref == waves[ii,:]) == np.size(wave_ref):
            # do not interpolate if the wavelength is exactly same with wave_ref
            fluxes_inter[ii,:] = fluxes[ii,:].copy()
            ivars_inter[ii,:] = ivars[ii, :].copy()
            masks_inter[ii, :] = mask_ii.copy()
        else:
            flux_inter_ii = interpolate.interp1d(waves[ii,:][mask_ii],fluxes[ii,:][mask_ii],kind='cubic',\
                                     bounds_error=False,fill_value=0.)(wave_ref)
            ivar_inter_ii = interpolate.interp1d(waves[ii,:][mask_ii],ivars[ii, :][mask_ii], kind='cubic',\
                                     bounds_error=False,fill_value=0.)(wave_ref)
            mask_inter_ii = interpolate.interp1d(waves[ii,:][mask_ii],masks_float[ii,:][mask_ii],kind='cubic',\
                                         bounds_error=False,fill_value=0.)(wave_ref)
            mask_inter_ii = (mask_inter_ii>0.5) & (ivar_inter_ii>0.) & (flux_inter_ii!=0.)

            fluxes_inter[ii,:] = flux_inter_ii #* ratio_ii
            ivars_inter[ii,:] = ivar_inter_ii #* ratio_ii
            masks_inter[ii,:] = mask_inter_ii

    return fluxes_inter,ivars_inter,masks_inter



def poly_ratio_fitfunc_chi2(theta, flux_ref, ivar_ref, thismask, arg_dict):
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
    mask_ref = arg_dict['mask_ref']
    wave = arg_dict['wave']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    func = arg_dict['func']
    # Evaluate the polynomial for rescaling
    ymult = utils.func_val(theta, wave, func, minx=wave_min, maxx=wave_max)
    flux_scale = ymult*flux
    mask_both = mask & mask_ref & thismask
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

def poly_ratio_fitfunc(flux_ref, ivar_ref, thismask, arg_dict, **kwargs_opt):

    # flux_ref, ivar_ref act like the 'data', the rescaled flux will be the 'model'

    # Function that we are optimizing
    result = scipy.optimize.differential_evolution(poly_ratio_fitfunc_chi2, args=(flux_ref, ivar_ref, thismask, arg_dict,), **kwargs_opt)
    flux = arg_dict['flux']
    ivar = arg_dict['ivar']
    mask = arg_dict['mask']
    mask_ref = arg_dict['mask_ref']
    wave = arg_dict['wave']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    func = arg_dict['func']
    # Evaluate the polynomial for rescaling
    ymult = utils.func_val(result.x, wave, func, minx=wave_min, maxx=wave_max)
    flux_scale = ymult*flux
    mask_both = mask & mask_ref & thismask
    totvar = utils.calc_ivar(ivar_ref) + ymult**2*utils.calc_ivar(ivar)
    ivartot = mask_both*utils.calc_ivar(totvar)

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
maks_ref = (ivar_ref > 0.0)
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
sys.exit(-1)

idx = 5
wave = wave_ref
flux = flux_inter[:,idx]
ivar = ivar_inter[:,idx]
#mask =

#if norder is None:
#    # If user did not specify an order to use, determine it automatically based on S/N
#    rms_snr, _ = coadd.sn_weights(wave, flux, ivar, mask=mask)
#    rms_snr_ref, _ = coadd.sn_weights(wave, flux_ref, ivar_ref, mask=mask_ref)


def solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, norder, mask = None, mask_ref = None, func='legendre',
                     maxiter=3, sticky=True, use_mad=True,
                     lower=3.0, upper=3.0):

    if mask is None:
        mask = (ivar > 0.0)
    if mask_ref is None:
        mask_ref = (ivar_ref > 0.0)

    arg_dict = dict(flux = flux, ivar = ivar, mask = mask, mask_ref = mask_ref,
                    wave = wave, wave_min = wave.min(),
                    wave_max = wave.max(), func = func, norder = norder)
    # Determine an initial guess


    result, ymodel, outmask = utils.robust_optimize(flux_ref, poly_ratio_fitfunc, arg_dict, invvar=invvar, inmask=inmask,
                                                    maxiter=maxiter, lower=lower, upper=upper, sticky=sticky,
                                                    use_mad=use_mad)


