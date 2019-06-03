

import numpy as np
from scipy.interpolate import interp1d
from astropy import stats
from astropy.stats import sigma_clip,sigma_clipped_stats
from pypeit.core import load
from pypeit import utils
from pypeit import msgs
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

    #ToDo: make it also works for single fits frame.
    nexp = len(fnames)
    sobjs0,header0 = load.load_specobjs(fnames[0], order=order)
    nspec = sobjs0[0].optimal['COUNTS'].size

    waves = np.zeros((nexp,nspec))
    fluxes = np.zeros_like(waves)
    ivars = np.zeros_like(waves)
    masks = np.zeros_like(waves,dtype=bool)

    for iexp in range(nexp):
        specobjs, headers = load.load_specobjs(fnames[iexp], order=order)

        for indx, spobj in enumerate(specobjs):
            if gdobj[iexp] in spobj.idx:
                ext = indx
        # Error message if this fails

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

        waves[iexp,:] = wave
        fluxes[iexp,:] = flux
        ivars[iexp,:] = ivar
        masks[iexp,:] = mask

    return waves,fluxes,ivars,masks

def interp_spec(wave_ref, waves,fluxes,ivars,masks):
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

    # Interpolate spectra to have the same wave grid with the iexp spectrum.
    # And scale spectra to the same flux level with the iexp spectrum.
    waves_inter = np.zeros_like(waves)
    fluxes_inter = np.zeros_like(fluxes)
    ivars_inter = np.zeros_like(ivars)
    masks_inter = np.zeros_like(masks)

    nexp = np.shape(fluxes)[0]
    for ii in range(nexp):

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

            fluxes_inter[ii,:] = flux_inter_ii #* ratio_ii
            ivars_inter[ii,:] = ivar_inter_ii #* ratio_ii
            masks_inter[ii,:] = mask_inter_ii

    return waves_inter,fluxes_inter,ivars_inter,masks_inter

def sn_weights(waves, fluxes, ivars, masks, dv_smooth=10000.0, const_weights=False, debug=False, verbose=False):
    """ Calculate the S/N of each input spectrum and create an array of (S/N)^2 weights to be used
    for coadding.

    Parameters
    ----------
    fluxes: float ndarray, shape = (nexp, nspec)
        Stack of (nexp, nspec) spectra where nexp = number of exposures, and nspec is the length of the spectrum.
    sigs: float ndarray, shape = (nexp, nspec)
        1-sigm noise vectors for the spectra
    masks: bool ndarray, shape = (nexp, nspec)
        Mask for stack of spectra. True=Good, False=Bad.
    waves: flota ndarray, shape = (nspec,) or (nexp, nspec)
        Reference wavelength grid for all the spectra. If wave is a 1d array the routine will assume
        that all spectra are on the same wavelength grid. If wave is a 2-d array, it will use the individual

    Optional Parameters:
    --------------------
    dv_smooth: float, 10000.0
         Velocity smoothing used for determining smoothly varying S/N ratio weights.

    Returns
    -------
    rms_sn : array
        Root mean square S/N value for each input spectra
    weights : ndarray
        Weights to be applied to the spectra. These are signal-to-noise squared weights.
    """

    sigs = np.sqrt(utils.calc_ivar(ivars))

    if fluxes.ndim == 1:
        nstack = 1
        nspec = fluxes.shape[0]
        flux_stack = fluxes.reshape((nstack, nspec))
        sig_stack = sigs.reshape((nstack,nspec))
        mask_stack = masks.reshape((nstack, nspec))
    elif fluxes.ndim == 2:
        nstack = fluxes.shape[0]
        nspec = fluxes.shape[1]
        flux_stack = fluxes
        sig_stack = sigs
        mask_stack = masks
    else:
        msgs.error('Unrecognized dimensionality for flux')

    # if the wave
    if waves.ndim == 1:
        wave_stack = np.outer(np.ones(nstack), waves)
    elif waves.ndim == 2:
        wave_stack = waves
    else:
        msgs.error('wavelength array has an invalid size')

    ivar_stack = utils.calc_ivar(sig_stack**2)
    # Calculate S/N
    sn_val = flux_stack*np.sqrt(ivar_stack)
    sn_val_ma = np.ma.array(sn_val, mask = np.invert(mask_stack))
    sn_sigclip = stats.sigma_clip(sn_val_ma, sigma=3, maxiters=5)
    sn2 = (sn_sigclip.mean(axis=1).compressed())**2 #S/N^2 value for each spectrum
    rms_sn = np.sqrt(sn2) # Root Mean S/N**2 value for all spectra
    rms_sn_stack = np.sqrt(np.mean(sn2))

    if rms_sn_stack <= 3.0 or const_weights:
        if verbose:
            msgs.info("Using constant weights for coadding, RMS S/N = {:g}".format(rms_sn_stack))
        weights = np.outer(sn2, np.ones(nspec))
        return rms_sn, weights
    else:
        if verbose:
            msgs.info("Using wavelength dependent weights for coadding")
        weights = np.ones_like(flux_stack) #((fluxes.shape[0], fluxes.shape[1]))
        spec_vec = np.arange(nspec)
        for iexp in range(nstack):
            imask = mask_stack[iexp,:]
            wave_now = wave_stack[iexp, imask]
            spec_now = spec_vec[imask]
            dwave = (wave_now - np.roll(wave_now,1))[1:]
            dv = (dwave/wave_now[1:])*c_kms
            dv_pix = np.median(dv)
            med_width = int(np.round(dv_smooth/dv_pix))
            sn_med1 = scipy.ndimage.filters.median_filter(sn_val[iexp,imask]**2, size=med_width, mode='reflect')
            sn_med2 = np.interp(spec_vec, spec_now, sn_med1)
            #sn_med2 = np.interp(wave_stack[iexp,:], wave_now,sn_med1)
            sig_res = np.fmax(med_width/10.0, 3.0)
            gauss_kernel = convolution.Gaussian1DKernel(sig_res)
            sn_conv = convolution.convolve(sn_med2, gauss_kernel)
            weights[iexp,:] = sn_conv

        # Finish
        return rms_sn, weights


def median_ratio_flux(flux,sig,flux_iref,sig_iref,mask=None,mask_iref=None,
                      cenfunc='median',snr_cut=3.0, maxiters=5,sigma=3):
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
    if mask_iref is None:
        mask_iref = (sig_iref>0.) & (flux_iref/sig_iref>snr_cut)
    if mask is None:
        mask = (sig>0.) & (flux/sig>snr_cut)

    ## Calculate the ratio
    ratio = flux_iref / flux
    mask_all = mask & mask_iref & (flux>0.) & (flux_iref>0.)
    ratio_mean,ratio_median,ratio_std = sigma_clipped_stats(ratio,np.invert(mask_all),cenfunc=cenfunc,\
                                                            maxiters=maxiters,sigma=sigma)
    #ratio_array = np.ones_like(flux) * ratio_median

    return ratio_median

def scale_spec(waves,fluxes,ivars,masks=None,flux_iref=None,ivar_iref=None,mask_iref=None,
               iref=None,cenfunc='median',snr_cut=3.0, maxiters=5,sigma=3,
               scale_method='median',hand_scale=None,SN_MAX_MEDSCALE=20.,SN_MIN_MEDSCALE=0.5,
               dv_smooth=10000.0,const_weights=False,debug=False,verbose=False):
    '''
    Scale the spectra into the same page with the reference spectrum.

    Args:
        waves:
        fluxes:
        ivars:
        masks:
        flux_iref:
        ivar_iref:
        mask_iref:
        iref:
        cenfunc:
        snr_cut:
        maxiters:
        sigma:
    Returns:
    '''

    # ToDo: Actually do we need wavelength to be here? maybe not.
    # Just need to make sure all the spectra in the same wave grid.
    if waves.ndim == 2:
        if abs(np.mean(waves)-np.mean(waves[0,:]))>1e-6:
            raise ValueError("Your spectra are not in the same wave grid. Please run 'interp_spec' first!")
        else:
            wave = waves[0,:]
    else:
        wave = waves.copy()

    if masks is None:
        masks = (ivars>0.) & np.isfinite(fluxes)

    # estimates the SNR of each spectrum and the stacked mean SNR
    rms_sn, weights = sn_weights(waves, fluxes, ivars, masks, dv_smooth=dv_smooth, \
                                 const_weights=const_weights, debug=debug, verbose=verbose)
    rms_sn_stack = np.sqrt(np.mean(rms_sn**2))

    nexp = waves.shape[0]

    # if flux_iref is None, then use the iref spectrum as the reference.
    if flux_iref is None:
        if iref is None:
            # find the highest SNR spectrum as the reference
            iref = np.argmax(rms_sn)
        flux_iref = fluxes[iref, :]
        ivar_iref = ivars[iref, :]
        mask_iref = masks[iref, :]
        msgs.info('Using the {:} spectrum to scale your spectra'.format(iref))
    else:
        iref = None
        # ToDo: Need to think about how to deal with ivar_iref=None.
        if ivar_iref is None:
            sig_iref = flux_iref.copy() / snr_cut / 2.0
            ivar_iref = utils.calc_ivar(sig_iref**2)

    if mask_iref is None:
        mask_iref = (ivar_iref>0.) & np.isfinite(flux_iref)

    sig_iref = np.sqrt(utils.calc_ivar(ivar_iref))

    # Estimate the scale factors
    fluxes_scale= np.zeros_like(fluxes)
    ivars_scale= np.zeros_like(ivars)
    scales = []
    for iexp in range(nexp):
        if scale_method == 'hand':
            omethod = 'hand'
            # Input?
            if hand_scale is None:
                msgs.error("Need to provide hand_scale parameter, one value per spectrum")
            fluxes_scale[iexp,:] = fluxes[iexp,:] * hand_scale[iexp]
            ivars_scale[iexp,:] = ivars[iexp,:] * 1.0/hand_scale[iexp]**2
            scales.append(hand_scale[qq])

        elif ((rms_sn_stack <= SN_MAX_MEDSCALE) and (rms_sn_stack > SN_MIN_MEDSCALE)) or (scale_method=='median'):
            omethod = 'median_flux'
            if iexp == iref:
                fluxes_scale[iexp, :] = fluxes[iexp, :].copy()
                ivars_scale[iexp, :] = ivars[iexp, :].copy()
                scales.append(1.)
                continue
            # Median ratio (reference to spectrum)
            #med_scale = median_ratio_flux(spectra, smask, qq, iref)
            flux = fluxes[iexp,:]
            sig = np.sqrt(utils.calc_ivar(ivars[iexp,:]))
            mask = masks[iexp,:]
            med_scale = median_ratio_flux(flux,sig,flux_iref,sig_iref,mask=mask,mask_iref=mask_iref,
                      cenfunc=cenfunc,snr_cut=snr_cut, maxiters=maxiters,sigma=sigma)
            # Apply
            med_scale= np.minimum(med_scale, 10.0)
            fluxes_scale[iexp,:] = fluxes[iexp,:] * med_scale
            ivars_scale[iexp,:] = ivars[iexp,:] * 1.0/med_scale**2
            scales.append(med_scale)

        elif rms_sn_stack <= SN_MIN_MEDSCALE:
            omethod = 'none_SN'
        elif (rms_sn_stack > SN_MAX_MEDSCALE) or scale_method=='poly':
            msgs.work("Should be using poly here, not median")
            omethod = 'median_flux'
            if iexp == iref:
                fluxes_scale[iexp, :] = fluxes[iexp, :].copy()
                ivars_scale[iexp, :] = ivars[iexp, :].copy()
                scales.append(1.)
                continue
            # Median ratio (reference to spectrum)
            #med_scale = median_ratio_flux(spectra, smask, qq, iref)
            flux = fluxes[iexp,:]
            sig = np.sqrt(utils.calc_ivar(ivars[iexp,:]))
            mask = masks[iexp,:]
            med_scale = median_ratio_flux(flux,sig,flux_iref,sig_iref,mask=mask,mask_iref=mask_iref,
                      cenfunc=cenfunc,snr_cut=snr_cut, maxiters=maxiters,sigma=sigma)
            # Apply
            med_scale= np.minimum(med_scale, 10.0)
            fluxes_scale[iexp,:] = fluxes[iexp,:] * med_scale
            ivars_scale[iexp,:] = ivars[iexp,:] * 1.0/med_scale**2
            scales.append(med_scale)
        else:
            msgs.error("Scale method not recognized! Check documentation for available options")
    # Finish
    return fluxes_scale,ivars_scale,scales, omethod


def long_clean(waves,fluxes,ivars,masks=None,cenfunc='median', snr_cut=3.0, maxiters=5, sigma=2,
               scale_method='median',hand_scale=None, SN_MAX_MEDSCALE=20., SN_MIN_MEDSCALE=0.5,
               dv_smooth=10000.0,const_weights=False, debug=False, verbose=False):

    if masks is None:
        masks = (ivars>0.) & np.isfinite(fluxes)

    masks_out = masks.copy()

    nexp = np.shape(fluxes)[0]
    nspec  = np.shape(fluxes)[1]

    for iexp in range(nexp):
        wave_iexp,flux_iexp,ivar_iexp = waves[iexp,:],fluxes[iexp,:],ivars[iexp,:]

        # Interpolate spectra into the native wave grid of the iexp spectrum
        waves_inter, fluxes_inter, ivars_inter, masks_inter = interp_spec(waves[iexp], waves, fluxes, ivars, masks)

        # avsigclip the interpolated spectra to obtain a high SNR average -- This now comes from coadd2d
        flux_iref,flux_median,flux_std = stats.sigma_clipped_stats(fluxes_inter, mask=np.invert(masks_inter), mask_value=0.,
                                                             sigma=sigma,maxiters=maxiters,cenfunc=cenfunc,axis=0)

        # This stuff dissappears
        nused = np.sum(masks_inter,axis=0)
        mask_iref = nused == nexp
        sig2 = 1.0/(ivars_inter+(ivars_inter<=0))
        newsig2 = np.sum(sig2*masks_inter,axis=0)/(nused**2+(nused==0))
        ivar_iref = mask_iref/(newsig2+(newsig2<=0.0))

        # scale the spectra to the reference spectrum
        fluxes_scale, ivars_scale, scales, omethod = scale_spec(waves_inter,fluxes_inter,ivars_inter,masks=masks_inter,\
                                                        flux_iref=flux_iref,ivar_iref=ivar_iref,mask_iref=mask_iref,\
                                                        iref=None,cenfunc=cenfunc,snr_cut=snr_cut,maxiters=maxiters,\
                                                        sigma=sigma,scale_method=scale_method,hand_scale=hand_scale,\
                                                        SN_MAX_MEDSCALE=SN_MAX_MEDSCALE,SN_MIN_MEDSCALE=SN_MIN_MEDSCALE,\
                                                        dv_smooth=dv_smooth,const_weights=const_weights,\
                                                        debug=debug,verbose=verbose)

        from IPython import embed
        embed()

        # do rejections and get mask for the iexp spectrum (import from long_combspec.pro)

        # QA plots for the residual of the new masked iexp spectrum and the high SNR average.



def long_combspec(waves,fluxes,ivars,masks=None):

    ## Some numbers and variables
    ndim_wave = np.ndim(waves)
    ndim_spec = np.ndim(fluxes)
    spec_shape = np.shape(fluxes)
    nexp = spec_shape[0]
    nspec  = spec_shape[1]

    ## Some judgements here
    if (ndim_wave == 1) and (ndim_spec==1):
        msgs.warn('Only one spectrum, no coadding will be performed.')
        return
    if (ndim_wave == 1) and (ndim_spec>1):
        waves = np.reshape(np.tile(waves,nexp),np.shape(fluxes))

    ## Doing rejections and get the new masks
    masks_out = long_clean(waves, fluxes, ivars, masks=masks)

    ## resample all spectra (including mask) into the same grid
    msgs.work('Write a function to resample. No interpolation please !')

    ## Now all spectra with the new masks have been resampled to the sample grid. Coadd them!



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

long_clean(waves,fluxes,ivars,masks=masks,cenfunc='median', snr_cut=3.0, maxiters=5, sigma=2,
               scale_method='median',hand_scale=None, SN_MAX_MEDSCALE=20., SN_MIN_MEDSCALE=0.5,
               dv_smooth=10000.0,const_weights=False, debug=False, verbose=False)



aaaa
for i in range(fluxes_inter.shape[0]):
    plt.plot(waves_inter[i,:],fluxes_inter[i,:])
plt.show()

for i in range(fluxes_inter.shape[0]):
    plt.plot(waves_inter[i,:],fluxes_scale[i,:])
plt.show()


waves_inter,fluxes_inter,ivars_inter,masks_inter = interp_spec(waves,fluxes,ivars,masks,iref=0)

#rms_sn, weights = sn_weights(waves_inter, fluxes_inter, ivars_inter, masks_inter, dv_smooth=10000.0, \
#                             const_weights=False, debug=False, verbose=False)

fluxes_scale,ivars_scale,scales, omethod =  scale_spec(waves_inter,fluxes_inter,ivars_inter,masks=masks_inter,
                                                       flux_iref=None,ivar_iref=None,mask_iref=None,iref=None,
                                                       cenfunc='median',snr_cut=3.0, maxiters=5,sigma=3,
                                                       scale_method='median')