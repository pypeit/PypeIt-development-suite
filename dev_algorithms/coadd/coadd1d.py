
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



def gauss1(x, mean, sigma, area):
    ygauss = np.exp(-np.power(x - mean, 2.) / (2 * np.power(sigma, 2.)))
    norm = area / (sigma * np.sqrt(2 * np.pi))
    return norm * ygauss

def renormalize_errors_qa(chi,sigma_corr, sig_range = 6.0):

    n_bins = 50
    binsize = 2.0*sig_range/n_bins
    bins_histo = -sig_range + np.arange(n_bins)*binsize+binsize/2.0

    xvals = np.arange(-10.0,10,0.02)
    ygauss = gauss1(xvals,0.0,1.0,1.0)
    ygauss_new = gauss1(xvals,0.0,sigma_corr,1.0)
    plt.figure(figsize=(10,6))
    plt.hist(chi,bins=bins_histo,normed=True,histtype='step', align='mid',color='k',linewidth=3,label='Chi distribution')
    plt.plot(xvals,ygauss,'c-',lw=3,label='sigma=1')
    plt.plot(xvals,ygauss_new,'m--',lw=2,label='new sigma={:}'.format(round(sigma_corr,2)))
    plt.xlabel('Residual distribution')
    plt.xlim([-6.05,6.05])
    plt.legend(fontsize=13,loc=2)
    plt.show()
    plt.close()

    return

def renormalize_errors(chi, mask, clip = 6.0, max_corr = 5.0, debug=False):

    chi2 = chi**2
    igood = (chi2 < clip**2) & mask
    if (np.sum(igood) > 0):
        gauss_prob = 1.0 - 2.0 * scipy.stats.norm.cdf(-1.0)
        chi2_sigrej = np.percentile(chi2[igood], 100.0*gauss_prob)
        sigma_corr = np.sqrt(chi2_sigrej)
        if sigma_corr < 1.0:
            msgs.warn("Error renormalization found correction factor sigma_corr = {:f}".format(sigma_corr) +
                      " < 1." + msgs.newline() +
                      " Errors are overestimated so not applying correction".format(sigma_corr))
            sigma_corr = 1.0
        if sigma_corr > max_corr:
            msgs.warn("Error renormalization found sigma_corr/sigma = {:f} > {:f}." + msgs.newline() +
                      "Errors are severely underestimated." + msgs.newline() +
                      "Setting correction to sigma_corr = {:f}".format(sigma_corr, max_corr, max_corr))
            sigma_corr = max_corr

    else:
        msgs.warn('No good pixels in error_renormalize. Threre are probably issues with your data')
        sigma_corr = 1.0

    if debug:
        renormalize_errors_qa(chi, sigma_corr)

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
    ivartot1 = mask_both*utils.calc_ivar(totvar)
    # Now rescale the errors
    chi = (flux_scale - flux_ref)*np.sqrt(ivartot1)
    try:
        debug = arg_dict['debug']
    except KeyError:
        debug = False

    sigma_corr = renormalize_errors(chi, mask=thismask, debug=debug)
    ivartot = ivartot1/sigma_corr

    return result, flux_scale, ivartot


def solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, norder, mask = None, mask_ref = None,
                     scale_min = 0.05, scale_max = 100.0, func='legendre',
                     maxiter=3, sticky=True, lower=3.0, upper=3.0, min_good=0.05, debug=False):


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
    wave_min = wave.min()
    wave_max = wave.max()
    arg_dict = dict(flux = flux, ivar = ivar, mask = mask,
                    ivar_ref = ivar_ref, wave = wave, wave_min = wave_min,
                    wave_max = wave_max, func = func, norder = norder, guess = guess, debug=debug)

    result, ymodel, ivartot, outmask = utils.robust_optimize(flux_ref, poly_ratio_fitfunc, arg_dict,inmask=mask_ref,
                                                             maxiter=maxiter, lower=lower, upper=upper, sticky=sticky)
    ymult1 = utils.func_val(result.x, wave, func, minx=wave_min, maxx=wave_max)
    ymult = np.fmin(np.fmax(ymult1, scale_min), scale_max)
    flux_rescale = ymult*flux
    ivar_rescale = ivar/ymult

    if debug:
        plt.plot(wave,flux_ref, color='black', drawstyle='steps-mid', zorder=3, label='Reference spectrum')
        plt.plot(wave,flux, color='dodgerblue', drawstyle='steps-mid', zorder = 10, alpha = 0.5, label='Original spectrum')
        plt.plot(wave,ymodel, color='red', drawstyle='steps-mid', alpha=0.7, zorder=1, linewidth=2, label='Rescaled spectrum')
        plt.legend()
        plt.show()

    return ymult, flux_rescale, ivar_rescale, outmask






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
        fluxes_inter, ivars_inter, masks_inter = interp_oned(wave_new,waves,fluxes,ivars,masks)

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

def scale_spec(wave, flux, ivar, flux_ref, ivar_ref, mask=None, mask_ref=None, min_good=0.05,
               cenfunc='median',ref_percentile=20.0, maxiters=5, sigrej=3, max_median_factor=10,
               npoly = None, scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5, debug=True):
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
    if mask_ref is None:
        mask_ref = ivar_ref > 0.0

    # estimates the SNR of each spectrum and the stacked mean SNR
    rms_sn, weights = sn_weights(wave, flux, ivar, mask)
    rms_sn_stack = np.sqrt(np.mean(rms_sn**2))

    if scale_method is None:
        if rms_sn_stack > sn_max_medscale:
            scale_method = 'poly'
        elif ((rms_sn_stack <= sn_max_medscale) and (rms_sn_stack > sn_min_medscale)):
            scale_method = 'median'
        else:
            scale_method = 'none'

    # Estimate the scale factor
    if scale_method == 'hand':
        # Input?
        if hand_scale is None:
            msgs.error("Need to provide hand_scale parameter, single value")
        flux_scale = flux * hand_scale
        ivar_scale = ivar * 1.0/hand_scale**2
        scale = np.full(flux.size,hand_scale)
    elif scale_method == 'median':
        # Median ratio (reference to spectrum)
        med_scale = robust_median_ratio(flux,ivar,flux_ref,ivar_ref,ref_percentile=ref_percentile,min_good=min_good,\
                                        mask=mask, mask_ref=mask_ref,cenfunc=cenfunc, maxiters=maxiters,\
                                        max_factor=max_median_factor,sigrej=sigrej)
        # Apply
        flux_scale = flux * med_scale
        ivar_scale = ivar * 1.0/med_scale**2
        scale = np.full(flux.size,med_scale)
    elif scale_method == 'none':
        flux_scale = flux.copy()
        ivar_scale = ivar.copy()
        scale = np.ones_like(flux)
    elif scale_method == 'poly':
        # Decide on the order of the polynomial rescaling
        if npoly is None:
            if rms_sn_stack > 25.0:
                npoly = 5 # Is this stable?
            elif rms_sn_stack > 8.0:
                npoly = 3
            elif rms_sn_stack >= 5.0:
                npoly = 2
            else:
                npoly = 1
        scale, flux_scale, ivar_scale, outmask = solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, npoly,
                                                                      mask=mask, mask_ref=mask_ref, debug=debug)
    else:
        msgs.error("Scale method not recognized! Check documentation for available options")
    # Finish
    return flux_scale, ivar_scale, scale, scale_method


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

def resid_gauss_plot(chi,one_sigma,debug=False):

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
    plt.figure(figsize=(10,6))
    plt.hist(chi,bins=bins_histo,normed=True,histtype='step', align='mid',color='k',linewidth=3,label='Chi distribution')
    plt.plot(xvals,ygauss,'c-',lw=3,label='sigma=1')
    plt.plot(xvals,ygauss_new,'m--',lw=2,label='new sigma={:}'.format(round(one_sigma,2)))
    plt.xlabel('Residual distribution')
    plt.xlim([-6.05,6.05])
    plt.legend(fontsize=13,loc=2)

    if debug:
        plt.show()
    plt.close()

    return

def coadd_qa(wave, flux, ivar, mask=None, wave_coadd=None, flux_coadd=None, ivar_coadd=None, mask_coadd=None,
             qafile=None, debug=False):

    if mask is None:
        mask = ivar>0.

    plt.figure(figsize=(10,6))
    plt.plot(wave[np.invert(mask)], flux[np.invert(mask)],'s',zorder=10,mfc='None', mec='r', label='Rejected pixels')
    plt.plot(wave, flux, color='dodgerblue', linestyle='steps-mid',zorder=2, alpha=0.7,label='Single exposure')
    plt.plot(wave, np.sqrt(utils.calc_ivar(ivar)),zorder=3, color='0.7', linestyle='steps-mid')
    ymin = np.percentile(flux, 5)
    ymax = 2.0 * np.percentile(flux, 95)

    # plot coadd
    if (wave_coadd is not None) and (flux_coadd is not None) and (ivar_coadd is not None):
        if mask_coadd is None:
            mask_coadd = ivar_coadd>0.
        plt.plot(wave_coadd, flux_coadd, color='k', linestyle='steps-mid', lw=2, zorder=1, label='Coadd')
        ymin = np.percentile(flux_coadd, 5)
        ymax = 2.0 * np.percentile(flux_coadd, 95)

    # Plot transmission
    if (np.max(wave[mask])>9000.0):
        skytrans_file = resource_filename('pypeit', '/data/skisim/atm_transmission_secz1.5_1.6mm.dat')
        skycat = np.genfromtxt(skytrans_file,dtype='float')
        scale = 0.8*ymax
        plt.plot(skycat[:,0]*1e4,skycat[:,1]*scale,'m-',alpha=0.5,zorder=11)

    plt.ylim([ymin, ymax])
    plt.xlim([wave.min(), wave.max()])
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Flux')
    plt.legend(fontsize=13)

    if qafile is not None:
        if len(qafile.split('.'))==1:
            msgs.info("No fomat given for the qafile, save to PDF format.")
            qafile = qafile+'.pdf'
        plt.savefig(qafile,dpi=300)
        msgs.info("Wrote QA: {:s}".format(qafile))
    if debug:
        plt.show()

    return

def coaddspec_qa_old(waves,fluxes,ivars,masks,wave_stack,flux_stack,ivar_stack,mask_stack,
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

    return


def long_reject(waves, fluxes, ivars, masks, fluxes_stack, ivars_stack, do_offset=True,
                sigrej_final=3.,do_var_corr=True, qafile=None, SN_MAX = 20.0, debug=False):

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

        msgs.info('Measured effective rejection from distribution of chi^2')
        msgs.info('Instead of rejecting sigrej={:}. Use threshold sigrej_eff={:}'\
                  .format(sigrej_final,np.round(sigrej_eff,2)))

        gdtmp = (outmasks[iexp, :] >0) & (ivar_real>0.)
        chi = (iflux[gdtmp] - newflux_now[gdtmp] - offset) * np.sqrt(ivar_cap[gdtmp])
        # Plot Chi distribution
        if debug:
            renormalize_errors_qa(chi, one_sigma)
            # Compare individual exposoures with stack spectrum.
            coadd_qa(waves[iexp,:],iflux,ivar,mask=gdtmp,wave_coadd=waves[iexp,:],flux_coadd=newflux_now,\
                      ivar_coadd=ivars_stack[iexp,:], mask_coadd=None,qafile=qafile,debug=debug)

    return outmasks

def long_comb(waves, fluxes, ivars, masks,wave_method='pixel', wave_grid_min=None, wave_grid_max=None, \
              A_pix=None, v_pix=None, samp_fact = 1.0, cenfunc='median', ref_percentile=20.0, maxiters=5, sigrej=3, \
              scale_method='median', hand_scale=None, sn_max_medscale=20., sn_min_medscale=0.5, \
              dv_smooth=10000.0, const_weights=False, maxiter_reject = 5, sn_max_reject=20., \
              fill_val=None,qafile=None,outfile=None,verbose=False,debug=False):

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
    flux_ref, flux_median, flux_std = stats.sigma_clipped_stats(fluxes_inter, mask=np.invert(masks_inter), mask_value=0.,
                                                                 sigma=sigrej, maxiters=maxiters, cenfunc=cenfunc, axis=0)
    # ToDo: This stuff disappears
    nexp = np.shape(fluxes)[0]
    nused = np.sum(masks_inter, axis=0)
    mask_ref = nused == nexp
    sig2 = 1.0 / (ivars_inter + (ivars_inter <= 0))
    newsig2 = np.sum(sig2 * masks_inter, axis=0) / (nused ** 2 + (nused == 0))
    ivar_ref = mask_ref / (newsig2 + (newsig2 <= 0.0))

    nexp = np.shape(fluxes)[0]
    fluxes_scale = np.copy(fluxes)
    ivars_scale = np.copy(ivars)
    scales = np.ones_like(fluxes)
    for iexp in range(nexp):
        flux_iref, ivar_iref, mask_iref = interp_spec(waves[iexp,:], wave_grid, flux_ref, ivar_ref, mask_ref)
        flux_scale,ivar_scale,scale,omethod = scale_spec(waves[iexp,:],fluxes[iexp,:],ivars[iexp,:],\
                                                         flux_iref,ivar_iref,mask=masks[iexp,:],mask_ref=mask_iref,\
                                                         cenfunc=cenfunc,ref_percentile=ref_percentile,maxiters=maxiters, \
                                                         sigrej=sigrej,scale_method=scale_method,hand_scale=hand_scale, \
                                                         sn_max_medscale=sn_max_medscale,sn_min_medscale=sn_min_medscale,\
                                                         debug=debug)
        fluxes_scale[iexp,:] = flux_scale
        ivars_scale[iexp,:] = ivar_scale
        scales[iexp,:] = scale

    # Doing rejections and coadding based on the scaled spectra
    iIter = 0
    thismask = np.copy(masks)
    while iIter < maxiter_reject:
        wave_stack, flux_stack, ivar_stack, mask_stack = compute_stack(waves, fluxes_scale, ivars_scale,thismask,\
                                                                       wave_grid, weights)
        fluxes_native_stack, ivars_native_stack, masks_native_stack = interp_spec(waves, wave_stack, flux_stack, \
                                                                                  ivar_stack,mask_stack)
        if iIter == maxiter_reject -1:
            thismask = long_reject(waves, fluxes_scale, ivars_scale, thismask, fluxes_native_stack, \
                                   ivars_native_stack, SN_MAX=sn_max_reject, debug=debug)
        else:
            thismask = long_reject(waves, fluxes_scale, ivars_scale, thismask, fluxes_native_stack, \
                                   ivars_native_stack, SN_MAX=sn_max_reject, debug=False)

        iIter = iIter +1

    # Plot the final coadded spectrum
    coadd_qa(waves, fluxes_scale, ivars_scale, thismask, wave_stack, flux_stack, \
             ivar_stack, mask_stack,qafile=qafile, debug=debug)

    # Write to disk?
    if outfile is not None:
        if len(outfile.split('.'))==1:
            msgs.info("No fomat given for the outfile, save to fits format.")
            outfile = outfile+'.fits'
        write_to_fits(wave_stack, flux_stack, ivar_stack, mask_stack, outfile, clobber=True, fill_val=fill_val)

    return wave_stack, flux_stack, ivar_stack, mask_stack, scales

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


# Now do some QA. We need the option to call the same QA routines while we are iterating just in case for debuggin
# Loop over exposures and show the comparison of the the chi_distribution to a Gaussian, and the updated errors

# Loop over exposures show the individual exposure, compared to the stack, compared to the single object noise, and maybe
# also show the modified noise vector (not sure about this). Label the masked pixels.

#newmask = reject_update_mask(waves, fluxes, ivars, masks, sn_weights, scale_factors, waves_stack, flux_stack, ivar_stack)

# Once you have the outmask, we update the

