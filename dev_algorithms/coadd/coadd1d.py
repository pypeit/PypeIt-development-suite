
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy import stats
from astropy.io import fits
from astropy import convolution

from pkg_resources import resource_filename
from pypeit import utils
from pypeit import msgs
from pypeit.core import load
from pypeit.core.wavecal import wvutils

## Plotting parameters
plt.rcdefaults()
plt.rcParams['font.family'] = 'times new roman'
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["ytick.direction"] = 'in'
plt.rcParams["xtick.direction"] = 'in'
plt.rcParams["xtick.labelsize"] = 15
plt.rcParams["ytick.labelsize"] = 15
plt.rcParams["axes.labelsize"] = 17


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

        # Initialize ext
        ext = None
        for indx, spobj in enumerate(specobjs):
            if gdobj[iexp] in spobj.idx:
                ext = indx
        if ext is None:
            msgs.error('Can not find extension {:} in {:}.'.format(gdobj[iexp],fnames[iexp]))

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
        elif ex_value == 'BOX':
            wave = specobjs[ext].boxcar['WAVE']
            if flux_value:
                flux = specobjs[ext].boxcar['FLAM']
                ivar = specobjs[ext].boxcar['FLAM_IVAR']
            else:
                flux = specobjs[ext].boxcar['COUNTS']
                ivar = specobjs[ext].boxcar['COUNTS_IVAR']
        else:
            msgs.error('{:} is not recognized. Please change to either BOX or OPT.'.format(ex_value))

        waves[iexp,:] = wave
        fluxes[iexp,:] = flux
        ivars[iexp,:] = ivar
        masks[iexp,:] = mask

    return waves,fluxes,ivars,masks

def new_wave_grid(waves,wave_method='iref',iref=0,wave_grid_min=None,wave_grid_max=None,
                  A_pix=None,v_pix=None,samp_fact=1.0,**kwargs):
    """ Create a new wavelength grid for the spectra to be rebinned and coadded on

    Parameters
    ----------
    waves : masked ndarray
        Set of N original wavelength arrays
        nexp, nspec
    wave_method : str, optional
        Desired method for creating new wavelength grid.
        'iref' -- Use the first wavelength array (default)
        'velocity' -- Constant velocity
        'pixel' -- Constant pixel grid
        'concatenate' -- Meld the input wavelength arrays
    iref : int, optional
      Reference spectrum
    wave_grid_min: float, optional
      min wavelength value for the final grid
    wave_grid_max: float, optional
      max wavelength value for the final grid
    A_pix : float
      Pixel size in same units as input wavelength array (e.g. Angstroms)
      If not input, the median pixel size is calculated and used
    v_pix : float
      Pixel size in km/s for velocity method
      If not input, the median km/s per pixel is calculated and used
    samp_fact: float
      sampling factor to make the wavelength grid finer or coarser.  samp_fact > 1.0 oversamples (finer),
      samp_fact < 1.0 undersamples (coarser)

    Returns
    -------
    wave_grid : ndarray
        New wavelength grid, not masked
    """
    if not isinstance(waves, np.ma.MaskedArray):
        waves = np.ma.array(waves,mask=waves<10.0)

    if wave_method == 'velocity':  # Constant km/s
        spl = 299792.458
        if v_pix is None:
            # Find the median velocity of a pixel in the input
            dv = spl * np.abs(waves - np.roll(waves,1)) / waves   # km/s
            v_pix = np.median(dv)

        # to make the wavelength grid finer or coarser
        v_pix = v_pix/samp_fact

        # Generate wavelength array
        if wave_grid_min is None:
            wave_grid_min = np.min(waves)
        if wave_grid_max is None:
            wave_grid_max = np.max(waves)
        x = np.log10(v_pix/spl + 1)
        npix = int(np.log10(wave_grid_max/wave_grid_min) / x) + 1
        wave_grid = wave_grid_min * 10**(x*np.arange(npix))

    elif wave_method == 'pixel': # Constant Angstrom
        if A_pix is None:
            dA =  np.abs(waves - np.roll(waves,1))
            A_pix = np.median(dA)

        # Generate wavelength array
        if wave_grid_min is None:
            wave_grid_min = np.min(waves)
        if wave_grid_max is None:
            wave_grid_max = np.max(waves)
        wave_grid = wvutils.wavegrid(wave_grid_min, wave_grid_max + A_pix, \
                                     A_pix,samp_fact=samp_fact)

    elif wave_method == 'loggrid':
        dloglam_n = np.log10(waves) - np.roll(np.log10(waves), 1)
        dloglam = np.median(dloglam_n.compressed())
        wave_grid_max = np.max(waves)
        wave_grid_min = np.min(waves)
        loglam_grid = wvutils.wavegrid(np.log10(wave_grid_min), np.log10(wave_grid_max)+dloglam, \
                                       dloglam,samp_fact=samp_fact)
        wave_grid = 10**loglam_grid

    elif wave_method == 'concatenate':  # Concatenate
        # Setup
        loglam = np.log10(waves) # This deals with padding (0's) just fine, i.e. they get masked..
        nexp = waves.shape[0]
        newloglam = loglam[iref, :].compressed()  # Deals with mask
        # Loop
        for j in range(nexp):
            if j == iref:
                continue
            #
            iloglam = loglam[j,:].compressed()
            dloglam_0 = (newloglam[1]-newloglam[0])
            dloglam_n =  (newloglam[-1] - newloglam[-2]) # Assumes sorted
            if (newloglam[0] - iloglam[0]) > dloglam_0:
                kmin = np.argmin(np.abs(iloglam - newloglam[0] - dloglam_0))
                newloglam = np.concatenate([iloglam[:kmin], newloglam])
            #
            if (iloglam[-1] - newloglam[-1]) > dloglam_n:
                kmin = np.argmin(np.abs(iloglam - newloglam[-1] - dloglam_n))
                newloglam = np.concatenate([newloglam, iloglam[kmin:]])
        # Finish
        wave_grid = 10**newloglam

    elif wave_method == 'iref':
        wave_grid = waves[iref, :].compressed()

    else:
        msgs.error("Bad method for scaling: {:s}".format(wave_method))

    return wave_grid


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

def interp_spec(wave_new, waves, fluxes, ivars, masks):
    '''
    Interpolate all spectra to the page of wave_new
    Args:
        waves:
        fluxes:
        sigs:
        masks:
        iref:
    Returns:
    '''

    if (fluxes.ndim==2) and (wave_new.ndim==1):
        nexp = np.shape(fluxes)[0]
        # Interpolate spectra to have the same wave grid with the iexp spectrum.
        # And scale spectra to the same flux level with the iexp spectrum.
        fluxes_inter = np.zeros((nexp, wave_new.size))
        ivars_inter  = np.zeros((nexp, wave_new.size))
        masks_inter  = np.zeros((nexp, wave_new.size), dtype=bool)
        for ii in range(nexp):
            mask_ii = masks[ii, :]
            if np.sum(wave_new == waves[ii, :]) == np.size(wave_new):
                # do not interpolate if the wavelength is exactly same with wave_new
                fluxes_inter[ii, :] = fluxes[ii, :].copy()
                ivars_inter[ii, :] = ivars[ii, :].copy()
                masks_inter[ii, :] = mask_ii.copy()
            else:
                flux_inter_ii, ivar_inter_ii, mask_inter_ii = \
                    interp_oned(wave_new, waves[ii, :],fluxes[ii, :],ivars[ii, :], masks[ii, :])
                fluxes_inter[ii, :] = flux_inter_ii  # * ratio_ii
                ivars_inter[ii, :] = ivar_inter_ii  # * ratio_ii
                masks_inter[ii, :] = mask_inter_ii

    elif (fluxes.ndim==1) and (wave_new.ndim==1):
        fluxes_inter, ivars_inter, masks_inter = interp_oned(wave_new,fluxes,ivars,masks)

    elif (fluxes.ndim==1) and (wave_new.ndim==2):
        nexp = np.shape(wave_new)[0]
        fluxes_inter = np.zeros_like(wave_new)
        ivars_inter = np.zeros_like(wave_new)
        masks_inter = np.zeros_like(wave_new)

        for ii in range(nexp):
            if np.sum(wave_new[ii, :] == waves) == np.size(waves):
                # do not interpolate if the wavelength is exactly same with wave_new
                fluxes_inter[ii, :] = fluxes.copy()
                ivars_inter[ii, :] = ivars.copy()
                masks_inter[ii, :] = masks.copy()
            else:
                flux_inter_ii, ivar_inter_ii, mask_inter_ii = \
                    interp_oned(wave_new[ii, :], waves, fluxes, ivars, masks)
                fluxes_inter[ii, :] = flux_inter_ii  # * ratio_ii
                ivars_inter[ii, :] = ivar_inter_ii  # * ratio_ii
                masks_inter[ii, :] = mask_inter_ii

    return fluxes_inter, ivars_inter, masks_inter

def sn_weights(waves, fluxes, ivars, masks, dv_smooth=10000.0, const_weights=False, verbose=False):
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

def robust_median_ratio(flux,ivar,flux_ref,ivar_ref, ref_percentile=20.0, min_good=0.05, mask=None, mask_ref=None,
                        cenfunc='median', maxiters=5, max_factor = 10.0, sigrej = 3.0):
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
    if mask is None:
        mask = ivar > 0.0
    if mask_ref is None:
        mask_ref = ivar_ref > 0.0

    nspec = flux.size
    snr_ref = flux_ref * np.sqrt(ivar_ref)
    snr_ref_best = np.percentile(snr_ref[mask_ref], ref_percentile)
    calc_mask = (snr_ref > snr_ref_best) & mask_ref & mask
    if (np.sum(calc_mask) > min_good*nspec):
        # Take the best part of the higher SNR reference spectrum
        flux_ref_mean, flux_ref_median, flux_ref_std = \
            stats.sigma_clipped_stats(flux_ref,np.invert(calc_mask),cenfunc=cenfunc,maxiters=maxiters,sigma=sigrej)
        flux_dat_mean, flux_dat_median, flux_dat_std = \
            stats.sigma_clipped_stats(flux,np.invert(calc_mask),cenfunc=cenfunc,maxiters=maxiters,sigma=sigrej)
        if (flux_ref_median < 0.0) or (flux_dat_mean < 0.0):
            msgs.warn('Negative median flux found. Not rescaling')
            ratio = 1.0
        else:
            ratio = np.fmax(np.fmin(flux_ref_median/flux_dat_median, max_factor), 1.0/max_factor)
    else:
        msgs.warn('Found only {%d} good pixels for computing median flux ratio.' + msgs.newline() +
                  'No median rescaling applied'.format(np.sum(calc_mask)))
        ratio = 1.0

    return ratio

def scale_spec(wave, flux, ivar, flux_ref, ivar_ref, mask=None, mask_ref=None,
               cenfunc='median',ref_percentile=20.0, maxiters=5,sigrej=3,
               scale_method=None,hand_scale=None, sn_max_medscale=2.0,sn_min_medscale=0.5,
               dv_smooth=10000.0,const_weights=False,verbose=False):
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

    if mask is None:
        mask = ivar > 0.0
    if mask_ref is None
        mask_ref = ivar_ref > 0.0

    # estimates the SNR of each spectrum and the stacked mean SNR
    rms_sn, weights = sn_weights(wave, flux, ivar, mask, dv_smooth=dv_smooth, verbose) const_weights=const_weights, verbose=verbose)
    rms_sn_stack = np.sqrt(np.mean(rms_sn**2))

    nexp = fluxes.shape[0]

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
            med_scale = median_ratio_flux(flux,sig,flux_iref,sig_iref,mask=mask,mask_ref=mask_iref,
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


def compute_stack(waves,fluxes,ivars,masks,wave_grid,weights):
    '''
    Compute the stacked spectrum based on spectra and wave_grid with weights being taken into account.
    Args:
        waves:
        fluxes:
        ivars:
        masks:
        wave_grid:
        weights:
    Returns:
        weighted stacked wavelength, flux and ivar
    '''

    waves_flat = waves[masks].ravel()
    fluxes_flat = fluxes[masks].ravel()
    ivars_flat = ivars[masks].ravel()
    vars_flat = utils.calc_ivar(ivars_flat)
    weights_flat = weights[masks].ravel()

    # Counts how many pixels in each wavelength bin
    num_stack, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False)

    # Calculate the stacked wavelength
    wave_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=waves_flat*weights_flat)
    weights_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=weights_flat)
    wave_stack = wave_stack_total / (weights_total+(weights_total==0.))

    # Calculate the stacked flux
    flux_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=fluxes_flat*weights_flat)
    flux_stack = flux_stack_total / (weights_total+(weights_total==0.))

    # Calculate the stacked ivar
    var_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=vars_flat*weights_flat**2)
    var_stack = var_stack_total / (weights_total+(weights_total==0.))**2
    ivar_stack = utils.calc_ivar(var_stack)
    #sig_stack = np.sqrt(var_stack)

    # New mask for the stack
    mask_stack = (wave_stack>0) & (ivar_stack>0.) & (flux_stack!=0.)

    return wave_stack, flux_stack, ivar_stack, mask_stack

def resid_gauss_plot(chi,one_sigma):

    max = 6.0
    min = -6.0

    n_bins = 50
    n_tot = len(chi)

    binsize = (max - min) / n_bins
    bins_histo = min + np.arange(n_bins)*binsize+binsize/2.0

    def gauss1(x, mean, sigma, area):
        ygauss = np.exp(-np.power(x - mean, 2.) / (2 * np.power(sigma, 2.)))
        norm = area / (sigma * np.sqrt(2 * np.pi))

        return norm * ygauss

    xvals = np.arange(-10.0,10,0.02)
    ygauss = gauss1(xvals,0.0,1.0,1.0)
    ygauss_new = gauss1(xvals,0.0,one_sigma,1.0)
    plt.hist(chi,bins=bins_histo,normed=True,histtype='step', align='mid',color='k',linewidth=3,label='Chi distribution')
    plt.plot(xvals,ygauss,'c-',lw=3,label='sigma=1')
    plt.plot(xvals,ygauss_new,'m--',lw=2,label='new sigma={:}'.format(round(one_sigma,2)))
    plt.xlabel('Residual distribution')
    plt.xlim([-6.05,6.05])
    plt.legend(fontsize=13,loc=2)
    plt.show()

    return

def coaddspec_qa(waves,fluxes,ivars,masks,wave_stack,flux_stack,ivar_stack,mask_stack,
                 qafile=None,verbose=False):
    """
    QA plot for 1D coadd of spectra
    """

    plt.figure(figsize=(10,6))
    ax1 = plt.axes([0.07, 0.13, 0.9, 0.4])
    ax2 = plt.axes([0.07, 0.55,0.9, 0.4])
    plt.setp(ax2.get_xticklabels(), visible=False)

    # Plot individual exposures
    ymin = np.percentile(flux_stack[mask_stack], 5)
    ymax = 2.0 * np.percentile(flux_stack[mask_stack], 95)
    ylim = [ymin,ymax]
    xlim = [np.min(wave_stack[mask_stack]),np.max(wave_stack[mask_stack])]

    nexp = np.shape(fluxes)[0]
    cmap = plt.get_cmap('RdYlBu_r')
    for iexp in range(nexp):
        color = cmap(float(iexp) / nexp)
        ax1.plot(waves[iexp,:], fluxes[iexp,:], color=color,alpha=0.5)
        ax1.plot(waves[iexp,:][np.invert(masks[iexp,:])], fluxes[iexp,:][np.invert(masks[iexp,:])],\
                 's',mfc='None',mec='k')
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.set_xlabel('Wavelength (Angstrom)')
    ax1.set_ylabel('Flux')

    # Plot coadded spectrum
    if (np.max(wave_stack[mask_stack])>9000.0):
        skytrans_file = resource_filename('pypeit', '/data/skisim/atm_transmission_secz1.5_1.6mm.dat')
        skycat = np.genfromtxt(skytrans_file,dtype='float')
        scale = 0.8*ylim[1]
        ax2.plot(skycat[:,0]*1e4,skycat[:,1]*scale,'m-',alpha=0.5)

    ax2.plot(wave_stack[mask_stack], np.sqrt(utils.calc_ivar(ivar_stack))[mask_stack], ls='steps-',color='0.7')
    ax2.plot(wave_stack[mask_stack], flux_stack[mask_stack], ls='steps-',color='b')
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    ax2.set_ylabel('Flux')

    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.2)
    if qafile is not None:
        if len(qafile.split('.'))==1:
            msgs.info("No fomat given for the qafile, save to PDF format.")
            qafile = qafile+'.pdf'
        plt.savefig(qafile,dpi=300)
        msgs.info("Wrote coadd QA: {:s}".format(qafile))
    if verbose:
        plt.show()
    plt.close()

    return


def long_reject(waves, fluxes, ivars, masks, fluxes_stack, ivars_stack, do_offset=True,
                sigrej_final=3.,do_var_corr=True, qafile=None, SN_MAX = 20.0, check=False):

    nexp = np.shape(fluxes)[0]

    outmasks = np.copy(masks)

    # Loop on images to update noise model for rejection
    for iexp in range(nexp):

        # Grab the spectrum
        iflux = fluxes[iexp,:]
        ivar = ivars[iexp,:]
        imask = outmasks[iexp,:]

        # Grab the stack with the same grid with the iexp
        newflux_now = fluxes_stack[iexp,:]
        newvar = utils.calc_ivar(ivars_stack[iexp,:])

        # var_tot
        var_tot = newvar + utils.calc_ivar(ivar)
        ivar_real = utils.calc_ivar(var_tot)
        # smooth out possible outliers in noise
        var_med = scipy.ndimage.filters.median_filter(var_tot, size=5, mode='reflect')
        var_smooth = scipy.ndimage.filters.median_filter(var_tot, size=99, mode='reflect')
        # conservatively always take the largest variance
        var_final = np.maximum(var_med, var_smooth)
        ivar_final = utils.calc_ivar(var_final)
        # Cap S/N ratio at SN_MAX to prevent overly aggressive rejection
        ivar_cap = np.minimum(ivar_final, (SN_MAX / (newflux_now + (newflux_now <= 0.0))) ** 2)
        # ; adjust rejection to reflect the statistics of the distribtuion
        # ; of errors. This fixes cases where for not totally understood
        # ; reasons the noise model is not quite right and
        # ; many pixels are rejected.

        # ; Is the model offset relative to the data? If so take it out
        if do_offset:
            diff1 = iflux - newflux_now
            # idum = np.where(arrmask[*, j] EQ 0, nnotmask)
            nnotmask = np.sum(imask)
            nmed_diff = np.maximum(nnotmask // 20, 10)
            # ; take out the smoothly varying piece
            # ; JXP -- This isnt going to work well if the data has a bunch of
            # ; null values in it
            w = np.ones(5, 'd')
            diff_med = scipy.ndimage.filters.median_filter(diff1 * (imask), size=nmed_diff, mode='reflect')
            diff_sm = np.convolve(diff_med, w / w.sum(), mode='same')
            chi2 = (diff1 - diff_sm) ** 2 * ivar_real
            goodchi = (imask) & (ivar_real > 0.0) & (chi2 <= 36.0)  # AND masklam, ngd)
            if np.sum(goodchi) == 0:
                goodchi = np.array([True] * iflux.size)

            offset_mean, offset, offset_std = stats.sigma_clipped_stats(diff1[goodchi],sigma=3., maxiters=5,cenfunc='median')
            # djs_iterstat, (arrflux[goodchi, j]-newflux_now[goodchi]) $
            #   , invvar = ivar_real[goodchi], mean = offset_mean $
            #   , median = offset $
        else:
            offset = 0.

        chi2 = (iflux - newflux_now - offset) ** 2 * ivar_real
        goodchi = imask & (ivar_real > 0.0) & (chi2 <= 36.0)  # AND masklam, ngd)
        ngd = np.sum(goodchi)
        if ngd == 0:
            goodchi = np.array([True] * iflux.size)
        # ; evalute statistics of chi2 for good pixels and excluding
        # ; extreme 6-sigma outliers
        chi2_good = chi2[goodchi]
        chi2_srt = chi2_good.copy()
        chi2_srt.sort()
        # ; evaluate at 1-sigma and then scale
        gauss_prob = 1.0 - 2.0 * (1. - scipy.stats.norm.cdf(1.))  # gaussint(-double(1.0d))
        sigind = int(np.round(gauss_prob * ngd))
        chi2_sigrej = chi2_srt[sigind]
        one_sigma = np.minimum(np.maximum(np.sqrt(chi2_sigrej), 1.0), 5.0)
        sigrej_eff = sigrej_final * one_sigma
        chi2_cap = (iflux - newflux_now - offset) ** 2 * ivar_cap
        # Grow??
        # Is this correct? This is not growing mask
        # chi_mask = (chi2_cap > sigrej_eff**2) & (~rmask[qq,:])
        chi_mask = (chi2_cap > sigrej_eff ** 2) | np.invert(imask)
        nrej = np.sum(chi_mask)
        # Apply
        if nrej > 0:
            msgs.info("Rejecting {:d} pixels in exposure {:d}".format(nrej, iexp))
            # print(rspec.data['wave'][qq,chi_mask])
            outmasks[iexp, chi_mask] = False
            # rspec.select = qq
            # rspec.add_to_mask(chi_mask)
            # outmask[*, j] = (arrmask[*, j] EQ 1) OR (chi2_cap GT sigrej_eff^2)

        if check:
            msgs.info('Measured effective rejection from distribution of chi^2')
            msgs.info('Instead of rejecting sigrej={:}. Use threshold sigrej_eff={:}'\
                      .format(sigrej_final,np.round(sigrej_eff,2)))

            gdtmp = (outmasks[iexp, :] >0) & (ivar_real>0.)
            chi = (iflux[gdtmp] - newflux_now[gdtmp] - offset) * np.sqrt(ivar_cap[gdtmp])
            # Plot Chi distribution
            resid_gauss_plot(chi, one_sigma)
            # Compare individual exposoures with stack spectrum.
            plt.plot(waves[iexp,:],newflux_now,'k-',lw=2,label='Coadd model')
            plt.plot(waves[iexp,:],iflux,'b-',alpha=0.7,label='{:}th exposure'.format(iexp+1))
            plt.plot(waves[iexp,:],np.sqrt(utils.calc_ivar(ivar)),'c:')
            plt.plot(waves[iexp,:][np.invert(gdtmp)],iflux[np.invert(gdtmp)],'s',mfc='None',
                        mec='r',label='Rejected pixels')
            ymin = np.percentile(newflux_now,5)
            ymax = 2.0*np.percentile(newflux_now,95)
            plt.ylim([ymin,ymax])
            plt.xlim([waves[iexp,:].min(),waves[iexp,:].max()])
            plt.xlabel('Wavelength (Angstrom)')
            plt.ylabel('Flux')
            plt.legend(fontsize=13)
            plt.show()

    return outmasks

def long_comb(waves, fluxes, ivars, masks,wave_method='pixel', wave_grid_min=None, wave_grid_max=None, \
              A_pix=None, v_pix=None, samp_fact = 1.0, cenfunc='median', snr_cut=2.0, maxiters=5, sigma=3, \
              scale_method='median', hand_scale=None, SN_MAX_MEDSCALE=20., SN_MIN_MEDSCALE=0.5, \
              dv_smooth=10000.0, const_weights=False, maxiter_reject = 5, SN_MAX_REJECT=20., \
              fill_val=None,qafile=None,outfile=None,verbose=False):

    # Define a common fixed wavegrid
    wave_grid = new_wave_grid(waves,wave_method=wave_method,wave_grid_min=wave_grid_min,wave_grid_max=wave_grid_max,
                      A_pix=A_pix,v_pix=v_pix,samp_fact=samp_fact)

    # Evaluate the sn_weights. This is done once at the beginning
    rms_sn, weights = sn_weights(waves,fluxes,ivars,masks, dv_smooth=dv_smooth, \
                                 const_weights=const_weights, verbose=verbose)

    # get_scale_factors == this will use the common wavegrid and interpolate spectra on it. We can be sinful about
    # covariance for the purposes of determining scale factors. This should return:
    # the scale factors evaluated on each native wavegrid, the prelim coadd evaluated on each native wavegrid

    # Interpolate spectra into the native wave grid of the iexp spectrum
    fluxes_inter, ivars_inter, masks_inter = interp_spec(wave_grid, waves, fluxes, ivars, masks)

    # avsigclip the interpolated spectra to obtain a high SNR average -- ToDo: This now comes from coadd2d
    flux_iref, flux_median, flux_std = stats.sigma_clipped_stats(fluxes_inter, mask=np.invert(masks_inter), mask_value=0.,
                                                                 sigma=sigma, maxiters=maxiters, cenfunc=cenfunc, axis=0)
    # ToDo: This stuff disappears
    nexp = np.shape(fluxes)[0]
    nused = np.sum(masks_inter, axis=0)
    mask_iref = nused == nexp
    sig2 = 1.0 / (ivars_inter + (ivars_inter <= 0))
    newsig2 = np.sum(sig2 * masks_inter, axis=0) / (nused ** 2 + (nused == 0))
    ivar_iref = mask_iref / (newsig2 + (newsig2 <= 0.0))

    fluxes_scale, ivars_scale, scales, omethod = scale_spec(wave_grid, fluxes_inter, ivars_inter, masks=masks_inter, \
                                                            flux_iref=flux_iref, ivar_iref=ivar_iref, mask_iref=mask_iref, \
                                                            iref=None, cenfunc=cenfunc, snr_cut=snr_cut, maxiters=maxiters, \
                                                            sigma=sigma, scale_method=scale_method, hand_scale=hand_scale, \
                                                            SN_MAX_MEDSCALE=SN_MAX_MEDSCALE,
                                                            SN_MIN_MEDSCALE=SN_MIN_MEDSCALE, \
                                                            dv_smooth=dv_smooth, const_weights=const_weights, \
                                                            verbose=verbose)

    scale_array = np.transpose(np.ones_like(fluxes.T)*scales)
    fluxes_native_scale = fluxes*scale_array
    ivars_native_scale = ivars * 1.0/scale_array**2

    # Doing rejections and coadding based on the scaled spectra
    iIter = 0
    thismask = np.copy(masks)
    while iIter < maxiter_reject:
        wave_stack, flux_stack, ivar_stack, mask_stack = compute_stack(waves, fluxes_native_scale, ivars_native_scale, \
                                                           thismask, wave_grid, weights)
        fluxes_native_stack, ivars_native_stack, masks_native_stack = interp_spec(waves, wave_stack, flux_stack, \
                                                                                  ivar_stack,mask_stack)
        if iIter == maxiter_reject -1:
            thismask = long_reject(waves, fluxes_native_scale, ivars_native_scale, thismask, fluxes_native_stack, \
                                   ivars_native_stack, SN_MAX=SN_MAX_REJECT, check=verbose)
        else:
            thismask = long_reject(waves, fluxes_native_scale, ivars_native_scale, thismask, fluxes_native_stack, \
                                   ivars_native_stack, SN_MAX=SN_MAX_REJECT, check=False)

        iIter = iIter +1

    # Plot the final coadded spectrum
    coaddspec_qa(waves, fluxes_native_scale, ivars_native_scale, thismask, \
                 wave_stack, flux_stack, ivar_stack, mask_stack,qafile=qafile, verbose=verbose)

    # Write to disk?
    if outfile is not None:
        if len(outfile.split('.'))==1:
            msgs.info("No fomat given for the outfile, save to fits format.")
            outfile = outfile+'.fits'
        write_to_fits(wave_stack, flux_stack, ivar_stack, mask_stack, outfile, clobber=True, fill_val=fill_val)

    return wave_stack, flux_stack, ivar_stack, mask_stack, scale_array


def write_to_fits(wave, flux, ivar, mask, outfil, clobber=True, fill_val=None):
    """ Write to a multi-extension FITS file.
    unless select=True.
    Otherwise writes 1D arrays
    Parameters
    ----------
    outfil : str
      Name of the FITS file
    select : int, optional
      Write only the select spectrum.
      This will always trigger if there is only 1 spectrum
      in the data array
    clobber : bool (True)
      Clobber existing file?
    add_wave : bool (False)
      Force writing of wavelengths as array, instead of using FITS
      header keywords to specify a wcs.
    fill_val : float, optional
      Fill value for masked pixels
    """

    # Flux
    if fill_val is None:
        flux = flux[mask]
    else:
        flux[mask] = fill_val
    prihdu = fits.PrimaryHDU(flux)
    hdu = fits.HDUList([prihdu])
    prihdu.name = 'FLUX'

    # Error  (packing LowRedux style)
    sig = np.sqrt(utils.calc_ivar(ivar))
    if fill_val is None:
        sig = sig[mask]
    else:
        sig[mask] = fill_val
    sighdu = fits.ImageHDU(sig)
    sighdu.name = 'ERROR'
    hdu.append(sighdu)

    # Wavelength
    if fill_val is None:
        wave = wave[mask]
    else:
        wave[mask] = fill_val
    wvhdu = fits.ImageHDU(wave)
    wvhdu.name = 'WAVELENGTH'
    hdu.append(wvhdu)

    hdu.writeto(outfil, overwrite=clobber)
    msgs.info('Wrote spectrum to {:s}'.format(outfil))


import os
datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/GMOS/R400_Flux/')
fnames = [datapath+'spec1d_flux_S20180903S0136-J0252-0503_GMOS-S_1864May27T160716.387.fits',\
          datapath+'spec1d_flux_S20180903S0137-J0252-0503_GMOS-S_1864May27T160719.968.fits',\
          datapath+'spec1d_flux_S20180903S0138-J0252-0503_GMOS-S_1864May27T160723.353.fits',\
          datapath+'spec1d_flux_S20180903S0141-J0252-0503_GMOS-S_1864May27T160727.033.fits',\
          datapath+'spec1d_flux_S20180903S0142-J0252-0503_GMOS-S_1864May27T160730.419.fits',\
          datapath+'spec1d_flux_S20181015S0140-J0252-0503_GMOS-S_1864May27T185252.770.fits']
gdobj = ['SPAT1073-SLIT0001-DET03','SPAT1167-SLIT0001-DET03','SPAT1071-SLIT0001-DET03','SPAT1072-SLIT0001-DET03',\
         'SPAT1166-SLIT0001-DET03','SPAT1073-SLIT0001-DET03']

# parameters for load_1dspec_to_array
ex_value = 'OPT'
flux_value = True

# Reading data
waves,fluxes,ivars,masks = load_1dspec_to_array(fnames,gdobj=gdobj,order=None,ex_value=ex_value,flux_value=flux_value)

# Coadding
wave_stack, flux_stack, ivar_stack, mask_stack, scale_array = \
    long_comb(waves, fluxes, ivars, masks,wave_method='pixel', scale_method='median', \
              maxiter_reject = 5, qafile='J0252_gmos', outfile='J0252_gmos.fits', verbose=True)

from IPython import embed
embed()
# Now do some QA. We need the option to call the same QA routines while we are iterating just in case for debuggin
# Loop over exposures and show the comparison of the the chi_distribution to a Gaussian, and the updated errors

# Loop over exposures show the individual exposure, compared to the stack, compared to the single object noise, and maybe
# also show the modified noise vector (not sure about this). Label the masked pixels.

#newmask = reject_update_mask(waves, fluxes, ivars, masks, sn_weights, scale_factors, waves_stack, flux_stack, ivar_stack)

# Once you have the outmask, we update the



'''
long_clean(waves,fluxes,ivars,masks=masks,cenfunc='median', snr_cut=2.0, maxiters=5, sigma=2,
               scale_method='median',hand_scale=None, SN_MAX_MEDSCALE=20., SN_MIN_MEDSCALE=0.5,
               dv_smooth=10000.0,const_weights=False, debug=False, verbose=False)

for i in range(fluxes_inter.shape[0]):
    plt.plot(waves_inter[i,:],fluxes_inter[i,:])
plt.show()

for i in range(fluxes_inter.shape[0]):
    plt.plot(waves_inter[i,:],fluxes_scale[i,:])
plt.show()


waves_inter,fluxes_inter,ivars_inter,masks_inter = interp_spec(waves[0,:],waves,fluxes,ivars,masks)

#rms_sn, weights = sn_weights(waves_inter, fluxes_inter, ivars_inter, masks_inter, dv_smooth=10000.0, \
#                             const_weights=False, debug=False, verbose=False)

fluxes_scale,ivars_scale,scales, omethod =  scale_spec(waves_inter,fluxes_inter,ivars_inter,masks=masks_inter,
                                                       flux_iref=None,ivar_iref=None,mask_iref=None,iref=None,
                                                       cenfunc='median',snr_cut=2.0, maxiters=5,sigma=3,
                                                       scale_method='median')
                                                       
                                                       
                                                       
    ivars_new, ivars_tot = update_errors(waves, fluxes, ivars, flux_stack_native, ivar_stack_native)
    thismask, qdone = pydl.djs_reject(fluxes, flux_stack_native, outmask=thismask, inmask=inmask, invvar=ivars_tot,
                                      lower=lower, upper=upper, maxdev=maxdev, maxrej=maxrej,
                                      groupdim=groupdim, groupsize=groupsize, groupbadpix=groupbadpix, grow=grow,
                                      use_mad=use_mad, sticky=sticky)
    iIter += 1

if (iIter == maxiter) & (maxiter != 0):
    msgs.warn('Maximum number of iterations maxiter={:}'.format(maxiter) + ' reached in robust_polyfit_djs')
outmask = np.copy(thismask)
if np.sum(outmask) == 0:
    msgs.warn('All points were rejected!!! The fits will be zero everywhere.')
waves_stack, flux_stack, ivar_stack = compute_stack(waves, fluxes, ivars, newmasks, wave_grid, sn_weights)


def long_clean(waves,fluxes,ivars,masks=None,cenfunc='median', snr_cut=2.0, maxiters=5, sigma=2,
               scale_method='median',hand_scale=None, SN_MAX_MEDSCALE=20., SN_MIN_MEDSCALE=0.5,
               dv_smooth=10000.0,const_weights=False, debug=False, verbose=False):

    if masks is None:
        masks = (ivars>0.) & np.isfinite(fluxes)

    masks_out = masks.copy()

    nexp = np.shape(fluxes)[0]
    for iexp in range(nexp):
        msgs.info('Cleaning the {:}th exposures.'.format(iexp))
        wave_iexp,flux_iexp,ivar_iexp = waves[iexp,:],fluxes[iexp,:],ivars[iexp,:]

        # Interpolate spectra into the native wave grid of the iexp spectrum
        fluxes_inter, ivars_inter, masks_inter = interp_spec(waves[iexp], waves, fluxes, ivars, masks)

        # avsigclip the interpolated spectra to obtain a high SNR average -- ToDo: This now comes from coadd2d
        flux_iref,flux_median,flux_std = stats.sigma_clipped_stats(fluxes_inter, mask=np.invert(masks_inter), mask_value=0.,
                                                             sigma=sigma,maxiters=maxiters,cenfunc=cenfunc,axis=0)
        # ToDo: This stuff disappears
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

        # do rejections and get mask for the iexp spectrum (import from long_combspec.pro)
        # doing rejections on fluxes_scale, ivars_scale, masks_inter
        msgs.info('Cleaning the {:}th exposure.'.format(iexp))

        # QA plots for the residual of the new masked iexp spectrum and the high SNR average.
        if verbose:
            msgs.info('QA plot for the {:}th exposure.'.format(iexp))
            #call the QA plots here.

    return masks_out


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
'''