

import numpy as np
from scipy.interpolate import interp1d
from astropy import stats
from astropy.stats import sigma_clip,sigma_clipped_stats
from pypeit.core import load
from pypeit import utils
import matplotlib.pyplot as plt

def load_1dspec_to_array(fnames,gdobj=None,order=None,ex_value='OPT',flux_value=True):
    '''
    Load the spectra from the 1d fits file into arrays.
    If Echelle, you need to specify which order you want to load.
    It can NOT load all orders for Echelle data.
    Args:
        fnames:
        gdobj:
        extensions:
        order: set to None if longslit data
        ex_value:
        flux_value:
    Returns:
        waves:
        fluxes:
        ivars:
        masks:
    '''

    nimg = len(fnames)
    sobjs0,header0 = load.load_specobjs(fnames[0], order=order)
    nspec = sobjs0[0].optimal['COUNTS'].size

    waves = np.zeros((nimg,nspec))
    fluxes = np.zeros_like(waves)
    ivars = np.zeros_like(waves)
    masks = np.zeros_like(waves,dtype=bool)

    for ispec in range(nimg):
        specobjs, headers = load.load_specobjs(fnames[ispec], order=order)

        for indx, spobj in enumerate(specobjs):
            if gdobj[ispec] in spobj.idx:
                ext = indx

        ## unpack wave/flux/mask
        if ex_value == 'OPT':
            wave = specobjs[ext].optimal['WAVE']
            mask = specobjs[ext].optimal['MASK']
            if flux_value:
                flux = specobjs[ext].optimal['FLAM']
                ivar = specobjs[ext].optimal['FLAM_IVAR']
            else:
                flux = specobjs[ext].optimal['COUNTS']
                ivar = specobjs[ext].optimal['COUNTS_IVAR']
        else:
            wave = specobjs[ext].boxcar['WAVE']
            if flux_value:
                flux = specobjs[ext].boxcar['FLAM']
                ivar = specobjs[ext].boxcar['FLAM_IVAR']
            else:
                flux = specobjs[ext].boxcar['COUNTS']
                ivar = specobjs[ext].boxcar['COUNTS_IVAR']

        waves[ispec,:] = wave
        fluxes[ispec,:] = flux
        ivars[ispec,:] = ivar
        masks[ispec,:] = mask

    return waves,fluxes,ivars,masks

def interp_spec(waves,fluxes,ivars,masks,iref=0):
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

    wave_iref,flux_iref,ivar_iref = waves[iref,:], fluxes[iref,:],ivars[iref,:]
    #sig_iref = np.sqrt(utils.calc_ivar(ivar_iref))

    # Interpolate spectra to have the same wave grid with the ispec spectrum.
    # And scale spectra to the same flux level with the ispec spectrum.
    waves_inter = np.zeros_like(waves)
    fluxes_inter = np.zeros_like(fluxes)
    ivars_inter = np.zeros_like(ivars)
    masks_inter = np.zeros_like(masks)

    nspec = np.shape(fluxes)[0]
    for ii in range(nspec):

        waves_inter[ii, :] = wave_iref.copy()
        mask_ii = masks[ii,:]

        if ii == iref:
            fluxes_inter[ii,:] = fluxes[ii,:].copy()
            ivars_inter[ii,:] = ivars[ii, :].copy()
            masks_inter[ii, :] = mask_ii.copy()
        else:
            flux_inter_ii = interp1d(waves[ii,:][mask_ii],fluxes[ii,:][mask_ii],kind='cubic',\
                                     bounds_error=False,fill_value=0.)(wave_iref)
            ivar_inter_ii = interp1d(waves[ii,:][mask_ii],ivars[ii, :][mask_ii], kind='cubic',\
                                     bounds_error=False,fill_value=0.)(wave_iref)
            mask_inter_ii = interp1d(waves[ii,:][mask_ii],masks_float[ii,:][mask_ii],kind='cubic',\
                                         bounds_error=False,fill_value=0.)(wave_iref)
            mask_inter_ii = (mask_inter_ii>0.5) & (ivar_inter_ii>0.) & (flux_inter_ii!=0.)

            #sig_inter_ii = np.sqrt(utils.calc_ivar(ivar_inter_ii))
            #ratio_ii = cal_ratio(wave_iref, flux_iref, flux_inter_ii, sig_iref=sig_iref, sig=sig_inter_ii, \
            #                     cenfunc=cenfunc, snr_cut=snr_cut, maxiters=maxiters, sigma=sigma)

            fluxes_inter[ii,:] = flux_inter_ii #* ratio_ii
            ivars_inter[ii,:] = ivar_inter_ii #* ratio_ii
            masks_inter[ii,:] = mask_inter_ii

            #plt.plot(waves[ii,:],masks_float[ii,:],'k-')
            #plt.plot(waves_inter[ii,:],masks_inter[ii,:],'b-')
            #plt.plot(waves_inter[ii,:],fluxes_inter[ii,:]/np.max(fluxes_inter[ii,:]),'r-')

            #plt.plot(waves[ii,:],fluxes[ii,:],'k-')
            #plt.plot(waves_inter[ii,:][masks_inter[ii,:]], fluxes_inter[ii,:][masks_inter[ii,:]],'b-')
            #plt.show()

    return waves_inter,fluxes_inter,ivars_inter,masks_inter

def cal_ratio(wave,flux_iref,flux,sig_iref=None,sig=None,cenfunc='median',snr_cut=3.0, maxiters=5,sigma=3):
    '''
    Calculate the ratio between reference spectrum and your spectrum.
    Args:
        wave:
        flux_iref:
        flux:
        sig_iref:
        sig:
        snr_cut:
        maxiters:
        sigma:
    Returns:
        Ratio_array: ratio array with the same size with your spectrum
    '''

    ## If no sigma array then set to sigma = flux*snr_cut/2.0 of the spectrum
    if sig_iref is None:
        sig_iref = np.ones_like(flux_iref)*sigma_clipped_stats(flux_iref,cenfunc=cenfunc,maxiters=maxiters,\
                                                               sigma=sigma)[1]*snr_cut/2.0
    if sig is None:
        sig = np.ones_like(flux)*sigma_clipped_stats(flux,cenfunc=cenfunc,maxiters=maxiters,sigma=sigma)[1]*snr_cut/2.0

    ## Mask for reference spectrum and your spectrum
    mask_iref = (sig_iref>0.) & (flux_iref/sig_iref>snr_cut)
    mask = (sig>0.) & (flux/sig>snr_cut)

    ## Calculate the ratio
    ratio = flux_iref / flux
    mask_all = mask & mask_iref & (np.isfinite(ratio))
    ratio_mean,ratio_median,ratio_std = sigma_clipped_stats(ratio,np.invert(mask_all),cenfunc=cenfunc,\
                                                            maxiters=maxiters,sigma=sigma)
    ratio_array = np.ones_like(flux) * ratio_median

    return ratio_array

def scale_spec(waves,fluxes,ivars,flux_iref,ivar_iref,masks=None,mask_iref=None,\
               cenfunc='median',snr_cut=3.0, maxiters=5,sigma=3):


datapath = '/Users/feige/Dropbox/OBS_DATA/GMOS/GS_2018B_FT202/R400_Flux/'
fnames = [datapath+'spec1d_flux_S20180903S0136-J0252-0503_GMOS-S_1864May27T160716.387.fits',\
          datapath+'spec1d_flux_S20180903S0137-J0252-0503_GMOS-S_1864May27T160719.968.fits',\
          datapath+'spec1d_flux_S20180903S0138-J0252-0503_GMOS-S_1864May27T160723.353.fits',\
          datapath+'spec1d_flux_S20180903S0141-J0252-0503_GMOS-S_1864May27T160727.033.fits',\
          datapath+'spec1d_flux_S20180903S0142-J0252-0503_GMOS-S_1864May27T160730.419.fits',\
          datapath+'spec1d_flux_S20181015S0140-J0252-0503_GMOS-S_1864May27T185252.770.fits']
gdobj = ['SPAT1073-SLIT0001-DET03','SPAT1167-SLIT0001-DET03','SPAT1071-SLIT0001-DET03','SPAT1072-SLIT0001-DET03',\
         'SPAT1166-SLIT0001-DET03','SPAT1073-SLIT0001-DET03']
ex_value = 'OPT'
flux_value = True
waves,fluxes,ivars,masks = load_1dspec_to_array(fnames,gdobj=gdobj,order=None,ex_value=ex_value,flux_value=flux_value)
waves_inter,fluxes_inter,ivars_inter,masks_inter = interp_spec(waves,fluxes,ivars,masks,iref=0)
from IPython import embed
embed()

