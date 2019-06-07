


import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from pypeit import utils
from pypeit.core import coadd
from scipy import interpolate
from astropy import stats


def median_ratio_flux(flux,ivar,flux_ref,ivar_ref,mask=None,mask_ref=None,
                      cenfunc='median',snr_cut=2.0, maxiters=5,sigma=3):
    '''
    Calculate the ratio between reference spectrum and your spectrum.
    Need to be in the same wave grid !!!
    Args:
        flux:
        sig:
        flux_iref:
        sig_iref:
        mask:
        mask_iref:
        snr_cut:
        maxiters:
        sigma:
    Returns:
        ratio_median: median ratio
    '''

    ## Mask for reference spectrum and your spectrum
    if mask_ref is None:
        mask_ref = (ivar_ref > 0.0) & (flux_ref*np.sqrt(ivar_ref) > snr_cut)
    if mask is None:
        mask = (ivar > 0.0) & (flux*np.sqrt(ivar) > snr_cut)

    ## Calculate the ratio
    ratio = flux_ref/(flux + (flux == 0.0))
    mask_all = mask & mask_ref & (flux > 0.0) & (flux_ref > 0.0)
    ratio_mean,ratio_median,ratio_std = stats.sigma_clipped_stats(ratio,np.invert(mask_all),cenfunc=cenfunc,
                                                                  maxiters=maxiters,sigma=sigma)
    #ratio_array = np.ones_like(flux) * ratio_median

    return ratio_median


def robust_median(flux, ivar, mask = None, maxiters=5,sigrej=3):

    if mask is None:
        mask = (ivar > 0.0)

    spec_mean1, spec_median1, sped_std1 = stats.sigma_clipped_stats(flux, np.invert(mask),
                                                                 cenfunc='median',maxiters=maxiters, sigma=sigrej)
    med_mask = mask & (flux > 0.5*spec_median1)

    spec_mean, spec_median, sped_std = stats.sigma_clipped_stats(flux, np.invert(med_mask),
                                                                 cenfunc='median',maxiters=maxiters, sigma=sigrej)

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
