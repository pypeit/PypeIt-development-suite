'''
This is NOT a code.

Based on the discussions with Joe, Fred on Apr 22, 2019.

Basic longslit strategy:
  1. Doing rejections
     interpolate all other spectra to the grid of the iref's grid
     doing rejections and generate a mask for the iref spectrum.
     also solve the scaling factor (scale the iref spectrum to the median spectrum)
     Loop over all the spectra then you get masks for all individual spectrum with the native wave grid
  2. Coadd spectra with a fixed grid.
     Generate a wavelength grid, constant pixel for longslit
     The wavelength might be shifted by up to half pixel.
  3. coadd spectra with a non-fixed grid wavelength (i.e. compute the median/mean wavelength of
     all pixels within a dpix/dv relative to the iref's spectra)

Basic Echelle strategy:
  1. Loop over orders doing rejections and solve the scaling factors for each spectrum
  2. coadd each order and solve the overlap sticking factor for each order
  3. generate a giant wavelength grid and put all individual orders/spectra into the same grid
     all the individual orders/spectra are scaled to the median based on the scaling factors and the overlap sticking factors
  4. Coadding all individual orders/spectra at the same time with rejections.


A core function for generating mask, rescaling, dealing with overlap regions between different orders and coadding
'''

import numpy as np
from scipy.interpolate import interp1d
from astropy import stats
from astropy.stats import sigma_clip,sigma_clipped_stats

def poly_ratio():
    '''
    This is a function to calculate the polynomial ratio
    :return:
    '''

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
    ## ToDo: Remove the following few lines
    import numpy as np
    from scipy.interpolate import interp1d
    from astropy import stats
    from astropy.stats import sigma_clip, sigma_clipped_stats

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

def long_clean(waves,fluxes,sigs,sigma=3,maxiters=5,snr_cut=3.0,iref=None,cenfunc='median'):
    '''
    This will be a function to generate the mask and scaling factor/array for the spectra
    :param waves:
    :param fluxes:
    :param sigs:
    :param iref:
    :param sigma:
    :param maxiters:
    :param cenfunc:
    :param stdfunc:
    :return:
    '''

    # if iref is None then find the highest SNR one as the reference
    if iref is None:
        snrs = fluxes / sigs
        snr = np.nanmedian(snrs, axis=1)
        iref = np.argmax(snr)

    # empty array for mask and rescaling
    mask_iref_new = np.ones_like(fluxes,dtype=bool) # good is True, bad is false
    scale_array = np.ones_like(fluxes,dtype=float)

    wave_iref,flux_iref,sig_iref = waves[iref,:], fluxes[iref,:],sigs[iref,:]
    #mask_iref = sigs[iref, :] > 0
    #med_iref = fluxes[iref,:][mask_iref]/sigs[iref,:][mask_iref]
    #sn_iref = stats.sigma_clip(med_iref, sigma=3, maxiters=5)#np.median(med_iref)

    nspec = np.shape(fluxes)[0]
    for ispec in range(nspec):
        wave_ispec = waves[ispec,:]
        flux_ispec = fluxes[ispec,:]
        sig_ispec = sigs[ispec,:]

        flux_iref_inter = interp1d(wave_iref,flux_iref,kind='cubic',bounds_error=False,fill_value=0.)(wave_ispec)
        sig_iref_inter = interp1d(wave_iref,sig_iref,kind='cubic',bounds_error=False,fill_value=0.)(wave_ispec)

        ## ToDo: change the following function to poly_ratio
        ratio_ispec = cal_ratio(wave_ispec,flux_iref_inter,flux_ispec,sig_iref=sig_iref_inter,sig=sig_ispec, \
                                cenfunc=cenfunc, snr_cut=snr_cut, maxiters=maxiters,sigma=sigma)

        # Interpolate spectra to have the same wave grid with the ispec spectrum.
        # And scale spectra to the same flux level with the iref spectrum.
        flux_inter = np.zeros_like(fluxes)
        sig_inter = np.zeros_like(sigs)
        for ii in range(nspec):
            if ii == ispec:
                flux_inter[ii,:] = fluxes[ii,:].copy()*ratio_ispec
                sig_inter[ii,:] = sigs[ii, :].copy()*ratio_ispec
            else:
                flux_inter_ii = interp1d(waves[ii,:],fluxes[ii,:],kind='cubic',bounds_error=False,fill_value=0.)(wave_ispec)
                sig_inter_ii = interp1d(waves[ii,:],sigs[ii, :], kind='cubic',bounds_error=False,fill_value=0.)(wave_ispec)
                ratio_ii = cal_ratio(wave_ispec, flux_iref_inter, flux_inter_ii, sig_iref=sig_iref_inter, sig=sig_inter_ii, \
                                     cenfunc=cenfunc, snr_cut=snr_cut, maxiters=maxiters, sigma=sigma)
                flux_inter[ii,:] = flux_inter_ii * ratio_ii
                sig_inter[ii,:] = sig_inter_ii * ratio_ii

        ## Now all spectra should have been scaled to the flux level same with iref and in the ispec wave grid
        #for ii in range(nspec):
        #    plt.plot(wave_ispec,flux_inter[ii,:])
        #plt.show()
        #plt.close('all')
        mean,median,std = sigma_clipped_stats(flux_inter, mask=None, mask_value=0., sigma=sigma,maxiters=maxiters,
                            cenfunc=cenfunc,axis=0)
        mask_ispec = abs(flux_inter[ispec,:] - median)<sigma*std
    ## Import the rejection part from the long_combspec.pro to do the rejections
    ## and get the mask array for the iref spectrum .


    ## The following three lines won't be used
    #med = np.median(fluxes,axis=1)
    #med_all = np.reshape(np.repeat(med,npix),spec_shape)
    #flux_med = np.median(fluxes*med_all,axis=0)

    ## Make a plot of the residual distribution of the iref spectrum.

    return mask_iref_new


def long_comb(waves,fluxes,sigs,iref=0,wave_method='pixel'):
    '''
    This function will coadd the longslit spectra or one specific order of Echelle spectra
    :param waves:
    :param fluxes:
    :param sigs:
    :return:
    '''

    ## Some numbers and variables
    ndim_wave = np.ndim(waves)
    ndim_spec = np.ndim(fluxes)
    spec_shape = np.shape(fluxes)
    nspec = spec_shape[0]
    npix  = spec_shape[1]
    masks = np.ones_like(fluxes,dtype='bool') # True means good, False means masked
    rescales = np.zeros_like(fluxes) # Array for storing masks

    ## Some judgements here
    if (ndim_wave == 1) and (ndim_spec==1):
        msgs.warn('Only one spectrum, no coadding will be performed.')
        return
    if (ndim_wave == 1) and (ndim_spec>1):
        waves = np.reshape(np.tile(waves,3),np.shape(fluxes))

    ## Doing rejections
    for ispec in range(nspec):
        # calling the long_clean function to do the rejection and to get the mask
        # and rescaling factor/function for the ispec
        masks[ispec,:],rescales[ispec,:] = long_clean(waves,fluxes,sigs,iref=ispec)

    # Compute weights (either single values or smoothed SN2 weights)

    ## Coadd spectra with two different methods
    # Coadd spectra with a fixed grid.
    # The wavelength might be shifted by up to half pixel.
    # Generate a wavelength grid, constant pixel for longslit, constant velocity for Echelle

    # I should write a function to coadd

    # coadd spectra with a non-fixed grid wavelength (i.e. compute the median/mean wavelength of
    # all pixels within a dpix/dv relative to the iref's spectra)

    ## Save coadded spectra. What's the data model?

    return wave_fix,flux_fix,sig_fix, wave, flux, sig


def ech_comb(specfiles,objid):
    '''
    This function will be used to coadd Echelle spectra
    :param specfiles: I would prefer using the list of fits files rather than arrays here since different orders
            might have different sizes.
    :return:
    '''

    ## coadd spectra order by order using long_comb
    wave_fix_all, flux_fix_all, sig_fix_all = np.zeros((norder,npix)),np.zeros((norder,npix)),np.zeros((norder,npix))
    wave, flux, sig = np.zeros((norder,npix)),np.zeros((norder,npix)),np.zeros((norder,npix))
    for iord in range(norder):
        waves,fluxes,sigs = # read the iord spectra from all files of a specific objid

        # run long_comb return the coadded spectrum of this order with the fixed velocity grid method
        wave_fix, flux_fix, sig_fix, wave, flux, sig = long_comb(waves, fluxes, sigs, iref=0,wave_method='velocity')

    ## Scale different orders to the same page. If no overlap just times by 1.0

    ## Merge orders directly or do one more iteration with the rejection and then coadd
    if giantcoadd:
        # generate a giant wavelength grid

        # rectification the spectra to the above wave grid using the nearest grid point interpolation

        # Put all the individual orders of all spectra into a giant array.
        # Need to take the differences between different orders derived from the overlap regions into account

        # call the long_clean to do the rejection

        # compute weights (either single value or a smoothed SN2 weights)

        # coadd the giant spectra array

    else:
        ## merge orders
        # should use the sensfunc*(smoothed SN2) as the weights