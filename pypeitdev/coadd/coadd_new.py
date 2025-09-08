


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
