


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
        flux_ref:
        sig_ref:
        mask:
        mask_ref:
        snr_cut:
        maxiters:
        sigma:
    Returns:
        ratio_median: median ratio
    '''

    ## Mask for reference spectrum and your spectrum
    if mask_ref is None:
        mask_ref = (ivar_ref > 0.0) & (flux_ref*np.sqrt(ivar_ref) > snr_cut)
    if np.sum(mask_ref)<1:
        msgs.warn('Not a single pixel has SNR>{:}, estimate median ratio based on data with \
                   20-80 percentile'.format(snr_cut))
        p20 = np.percentile(flux_ref, 20)
        p80 = np.percentile(flux_ref, 80)
        mask_ref = (ivar_ref > 0.0) & (flux_ref>p20) & (flux_ref<p80)
    if mask is None:
        mask = (ivar > 0.0) & (flux*np.sqrt(ivar) > snr_cut)
    if np.sum(mask_ref)<1:
        msgs.warn('Not a single pixel has SNR>{:}, estimate median ratio based on data with \
                   20-80 percentile'.format(snr_cut))
        p20 = np.percentile(flux, 20)
        p80 = np.percentile(flux, 80)
        mask_ref = (ivar > 0.0) & (flux > p20) & (flux < p80)

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

def interp_oned(wave_new, wave_old, flux_old, ivar_old, mask_old):
    '''
    Args:
       wave_new: (one-D array) New wavelength
       wave_old: (one-D array) Old wavelength
       flux_old: (one-D array) Old flux
       ivar_old: (one-D array) Old ivar
       mask_old: (one-D array) Old float mask
    Returns :
       flux_new, ivar_new, mask_new (bool)
    '''

    # make the mask array to be float, used for interpolation
    masks_float = np.zeros_like(flux_old)
    masks_float[mask_old] = 1.0

    flux_new = interpolate.interp1d(wave_old[mask_old], flux_old[mask_old], kind='cubic', \
                                         bounds_error=False, fill_value=0.)(wave_new)
    ivar_new = interpolate.interp1d(wave_old[mask_old], ivar_old[mask_old], kind='cubic', \
                                         bounds_error=False, fill_value=0.)(wave_new)
    mask_new_tmp = interpolate.interp1d(wave_old[mask_old], masks_float[mask_old], kind='cubic', \
                                         bounds_error=False, fill_value=0.)(wave_new)
    mask_new = (mask_new_tmp > 0.5) & (ivar_new > 0.) & (flux_new != 0.)
    return flux_new,ivar_new,mask_new

def interp_spec(wave_ref, waves, fluxes, ivars, masks):
    '''
    Interpolate all spectra to the page of wave_ref
    Args:
        waves:
        fluxes:
        sigs:
        masks:
        iref:
    Returns:
    '''

    if (fluxes.ndim==2) and (wave_ref.ndim==1):

        nexp = np.shape(fluxes)[0]
        # Interpolate spectra to have the same wave grid with the iexp spectrum.
        # And scale spectra to the same flux level with the iexp spectrum.
        fluxes_inter = np.zeros((nexp, len(wave_ref)))
        ivars_inter = np.zeros((nexp, len(wave_ref)))
        masks_inter = np.zeros((nexp, len(wave_ref)), dtype=bool)

        for ii in range(nexp):
            mask_ii = masks[ii, :]
            if np.sum(wave_ref == waves[ii, :]) == np.size(wave_ref):
                # do not interpolate if the wavelength is exactly same with wave_ref
                fluxes_inter[ii, :] = fluxes[ii, :].copy()
                ivars_inter[ii, :] = ivars[ii, :].copy()
                masks_inter[ii, :] = mask_ii.copy()
            else:
                flux_inter_ii, ivar_inter_ii, mask_inter_ii = interp_oned(wave_ref, waves[ii, :], \
                                                                          fluxes[ii, :],ivars[ii, :], \
                                                                          masks[ii, :])
                fluxes_inter[ii, :] = flux_inter_ii  # * ratio_ii
                ivars_inter[ii, :] = ivar_inter_ii  # * ratio_ii
                masks_inter[ii, :] = mask_inter_ii

    elif (fluxes.ndim==1) and (wave_ref.ndim==1):
        fluxes_inter, ivars_inter, masks_inter = interp_oned(wave_ref,fluxes,ivars,masks)

    elif (fluxes.ndim==1) and (wave_ref.ndim==2):
        nexp = np.shape(wave_ref)[0]
        fluxes_inter = np.zeros_like(wave_ref)
        ivars_inter = np.zeros_like(wave_ref)
        masks_inter = np.zeros_like(wave_ref)

        for ii in range(nexp):
            if np.sum(wave_ref[ii, :] == waves) == np.size(waves):
                # do not interpolate if the wavelength is exactly same with wave_ref
                fluxes_inter[ii, :] = fluxes.copy()
                ivars_inter[ii, :] = ivars.copy()
                masks_inter[ii, :] = masks.copy()
            else:
                flux_inter_ii, ivar_inter_ii, mask_inter_ii = interp_oned(wave_ref[ii, :], waves, \
                                                                          fluxes, ivars, masks)
                fluxes_inter[ii, :] = flux_inter_ii  # * ratio_ii
                ivars_inter[ii, :] = ivar_inter_ii  # * ratio_ii
                masks_inter[ii, :] = mask_inter_ii

    return fluxes_inter,ivars_inter,masks_inter
