


import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from sklearn import mixture
from astropy.io import fits
from pypeit.core import pydl
import astropy.units as u
from astropy.io import fits
from astropy import stats
from pypeit.core import flux
from pypeit.core import load
from astropy import table
from pypeit.core import save
from pypeit.core import coadd2d
from pypeit.core import coadd1d
from pypeit.spectrographs import util
from pypeit import utils
from pypeit import msgs
import pickle
PYPEIT_FLUX_SCALE = 1e-17
from astropy.io import fits
import copy
from IPython import embed
import qso_pca
from pypeit.spectrographs.spectrograph import Spectrograph
from pypeit.spectrographs.util import load_spectrograph


##############################
#  Telluric model functions  #
##############################

def get_bounds_tell(tell_dict, resln_guess, resln_frac_bounds, pix_shift_bounds):

    # Set the bounds for the optimization
    bounds_tell = [(tell_dict['pressure_grid'].min(), tell_dict['pressure_grid'].max()),
                   (tell_dict['temp_grid'].min(), tell_dict['temp_grid'].max()),
                   (tell_dict['h2o_grid'].min(), tell_dict['h2o_grid'].max()),
                   (tell_dict['airmass_grid'].min(), tell_dict['airmass_grid'].max()),
                   (resln_guess * resln_frac_bounds[0], resln_guess * resln_frac_bounds[1]),
                   pix_shift_bounds]

    return bounds_tell

def read_telluric_grid(filename, wave_min=None, wave_max=None, pad = 0):

    hdul = fits.open(filename)
    wave_grid_full = 10.0*hdul[1].data
    model_grid_full = hdul[0].data
    nspec_full = wave_grid_full.size

    if wave_min is not None:
        ind_lower = np.argmin(np.abs(wave_grid_full - wave_min)) - pad
    else:
        ind_lower = 0
    if wave_max is not None:
        ind_upper = np.argmin(np.abs(wave_grid_full - wave_max)) + pad
    else:
        ind_upper=nspec_full
    wave_grid = wave_grid_full[ind_lower:ind_upper]
    model_grid = model_grid_full[:,:,:,:, ind_lower:ind_upper]

    pg = hdul[0].header['PRES0']+hdul[0].header['DPRES']*np.arange(0,hdul[0].header['NPRES'])
    tg = hdul[0].header['TEMP0']+hdul[0].header['DTEMP']*np.arange(0,hdul[0].header['NTEMP'])
    hg = hdul[0].header['HUM0']+hdul[0].header['DHUM']*np.arange(0,hdul[0].header['NHUM'])
    if hdul[0].header['NAM'] > 1:
        ag = hdul[0].header['AM0']+hdul[0].header['DAM']*np.arange(0,hdul[0].header['NAM'])
    else:
        ag = hdul[0].header['AM0']+1*np.arange(0,1)

    loglam = np.log10(wave_grid)
    dloglam = np.median(loglam[1:] - loglam[:-1])
    # Guess resolution from wavelength sampling of telluric grid if it is not provided
    resln_guess = 1.0/(3.0 * dloglam * np.log(10.0))  # assume roughly Nyquist sampling
    pix_per_R = 1.0/resln_guess / (dloglam * np.log(10.0)) / (2.0 * np.sqrt(2.0 * np.log(2)))
    tell_pad_pix = int(np.ceil(10.0 * pix_per_R))

    tell_dict = dict(wave_grid=wave_grid, dloglam=dloglam,
                     resln_guess=resln_guess, pix_per_R=pix_per_R, tell_pad_pix=tell_pad_pix,
                     pressure_grid=pg, temp_grid=tg, h2o_grid=hg, airmass_grid=ag, tell_grid=model_grid)
    return tell_dict


def interp_telluric_grid(theta,tell_dict):

    pg = tell_dict['pressure_grid']
    tg = tell_dict['temp_grid']
    hg = tell_dict['h2o_grid']
    ag = tell_dict['airmass_grid']
    model_grid = tell_dict['tell_grid']
    press,temp,hum,airmass = theta
    if len(pg) > 1:
        p_ind = int(np.round((press-pg[0])/(pg[1]-pg[0])))
    else:
        p_ind = 0
    if len(tg) > 1:
        t_ind = int(np.round((temp-tg[0])/(tg[1]-tg[0])))
    else:
        t_ind = 0
    if len(hg) > 1:
        h_ind = int(np.round((hum-hg[0])/(hg[1]-hg[0])))
    else:
        h_ind = 0
    if len(ag) > 1:
        a_ind = int(np.round((airmass-ag[0])/(ag[1]-ag[0])))
    else:
        a_ind = 0

    return model_grid[p_ind,t_ind,h_ind,a_ind]

def conv_telluric(tell_model,dloglam, res):

    pix_per_sigma = 1.0/res/(dloglam*np.log(10.0))/(2.0 * np.sqrt(2.0 * np.log(2))) # number of dloglam pixels per 1 sigma dispersion
    sig2pix = 1.0/pix_per_sigma # number of sigma per 1 pix
    #conv_model = scipy.ndimage.filters.gaussian_filter1d(tell_model, pix)
    # x = loglam/sigma on the wavelength grid from -4 to 4, symmetric, centered about zero.
    x = np.hstack([-1*np.flip(np.arange(sig2pix,4,sig2pix)),np.arange(0,4,sig2pix)])
    # g = Gaussian evaluated at x, sig2pix multiplied in to properly normalize the convolution
    g = (1.0/(np.sqrt(2*np.pi)))*np.exp(-0.5*(x)**2)*sig2pix
    conv_model = scipy.signal.convolve(tell_model,g,mode='same')
    return conv_model

def shift_telluric(tell_model, loglam, dloglam, shift):
    loglam_shift = loglam + shift*dloglam
    tell_model_shift = np.interp(loglam_shift, loglam, tell_model)
    return tell_model_shift


def eval_telluric(theta_tell, tell_dict):

    ntheta = len(theta_tell)
    tellmodel_hires = interp_telluric_grid(theta_tell[:4], tell_dict)
    tellmodel_conv = conv_telluric(tellmodel_hires, tell_dict['dloglam'], theta_tell[4])
    tell_pad = tell_dict['tell_pad']

    if ntheta == 6:
        tellmodel_out = shift_telluric(tellmodel_conv, tell_dict['loglam'], tell_dict['dloglam'], theta_tell[5])
        return tellmodel_out[tell_pad[0]:-tell_pad[1]]
    else:
        return tellmodel_conv[tell_pad[0]:-tell_pad[1]]


def trim_tell_dict(tell_dict, ind_lower, ind_upper):

    wave_grid = tell_dict['wave_grid']
    ind_lower_pad = np.fmax(ind_lower - tell_dict['tell_pad_pix'], 0)
    ind_upper_pad = np.fmin(ind_upper + tell_dict['tell_pad_pix'], wave_grid.size - 1)
    tell_pad_tuple = (ind_lower - ind_lower_pad, ind_upper_pad - ind_upper)
    tell_wave_grid = wave_grid[ind_lower_pad:ind_upper_pad + 1]
    tell_model_grid = tell_dict['tell_grid'][:, :, :, :, ind_lower_pad:ind_upper_pad + 1]
    tell_dict_fit = dict(pressure_grid=tell_dict['pressure_grid'], temp_grid=tell_dict['temp_grid'],
                         h2o_grid=tell_dict['h2o_grid'], airmass_grid=tell_dict['airmass_grid'],
                         tell_grid=tell_model_grid, tell_pad=tell_pad_tuple,
                         resln_guess=tell_dict['resln_guess'],pix_per_R=tell_dict['pix_per_R'],
                         dloglam=tell_dict['dloglam'], loglam=np.log10(tell_wave_grid))
    return tell_dict_fit

def sort_telluric(wave, wave_mask, tell_dict):

    norders = wave.shape[1]
    tell_med = np.zeros(norders)
    # Do a quick loop over all the orders to sort them in order of strongest to weakest telluric absorption
    for iord in range(norders):
        ind_lower = np.argmin(np.abs(tell_dict['wave_grid'] - np.min(wave[wave_mask[:,iord],iord])))
        ind_upper = np.argmin(np.abs(tell_dict['wave_grid'] - np.max(wave[wave_mask[:,iord],iord])))
        tm_grid = tell_dict['tell_grid'][:, :, :, :,ind_lower:ind_upper]
        tell_model_mid = tm_grid[tm_grid.shape[0]//2, tm_grid.shape[1]//2,tm_grid.shape[2]//2,tm_grid.shape[3]//2,:]
        tell_med[iord] = np.mean(tell_model_mid)

    # Perform fits in order of telluric strength
    srt_order_tell = tell_med.argsort()

    return srt_order_tell

############################
#  Object model functions  #
############################

# Sensitivity function evaluation function. Model for counts is flam_true/sensfunc
def sensfunc_obj_model(theta, obj_dict):

    wave_star = obj_dict['wave']
    wave_min = obj_dict['wave_min']
    wave_max = obj_dict['wave_max']
    flam_true = obj_dict['flam_true']
    func = obj_dict['func']

    sensfunc = np.exp(utils.func_val(theta, wave_star, func, minx=wave_min, maxx=wave_max))
    counts_model = flam_true/(sensfunc + (sensfunc == 0.0))

    return counts_model, (sensfunc > 0.0)

# QSO evaluation function. Model for QSO is a PCA spectrum
def qso_obj_model(theta, obj_dict):

    pca_model = qso_pca.pca_eval(theta, obj_dict['pca_dict'])
    # TODO Is the prior evaluation slowing things down??
    # TODO Disablingthe prior for now as I think it slows things down for no big gain
    #ln_pca_pri = qso_pca.pca_lnprior(theta_PCA, arg_dict['pca_dict'])
    #ln_pca_pri = 0.0
    #flux_model, tell_model, spec_model, modelmask
    return pca_model, (pca_model > 0.0)


############################
#  Fitting routines        #
############################

def tellfit_chi2(theta, flux, thismask, arg_dict):

    obj_model_func = arg_dict['obj_model_func']
    flux_ivar = arg_dict['ivar']

    theta_obj = theta[:-6]
    theta_tell = theta[-6:]
    tell_model = eval_telluric(theta_tell, arg_dict['tell_dict'])
    obj_model, modelmask = obj_model_func(theta_obj, arg_dict['obj_dict'])

    if not np.any(modelmask):
        return np.inf
    else:
        totalmask = thismask & modelmask
        chi_vec = totalmask * (flux - tell_model*obj_model) * np.sqrt(flux_ivar)
        robust_scale = 2.0
        huber_vec = scipy.special.huber(robust_scale, chi_vec)
        loss_function = np.sum(np.square(huber_vec * totalmask))
        return loss_function

def tellfit(flux, thismask, arg_dict, **kwargs_opt):

    # Unpack arguments
    obj_model_func = arg_dict['obj_model_func'] # Evaluation function
    flux_ivar = arg_dict['ivar'] # Inverse variance of flux or counts
    bounds = arg_dict['bounds']  # bounds for differential evolution optimizaton
    seed = arg_dict['seed']      # Seed for differential evolution optimizaton
    result = scipy.optimize.differential_evolution(tellfit_chi2, bounds, args=(flux, thismask, arg_dict,), seed=seed,
                                                   **kwargs_opt)

    theta_obj  = result.x[:-6]
    theta_tell = result.x[-6:]
    tell_model = eval_telluric(theta_tell, arg_dict['tell_dict'])
    obj_model, modelmask = obj_model_func(theta_obj, arg_dict['obj_dict'])
    totalmask = thismask & modelmask
    chi_vec = totalmask*(flux - tell_model*obj_model)*np.sqrt(flux_ivar)

    try:
        debug = arg_dict['debug']
    except KeyError:
        debug = False

    # Name of function for title in case QA requested
    obj_model_func_name = getattr(obj_model_func, '__name__', repr(obj_model_func))
    sigma_corr, maskchi = coadd1d.renormalize_errors(chi_vec, mask=totalmask, title = obj_model_func_name,
                                                     debug=debug)
    ivartot = flux_ivar/sigma_corr**2

    return result, tell_model*obj_model, ivartot



def unpack_orders(sobjs, ret_flam=False):

    # TODO This should be a general reader:
    #  For echelle:  read in all the orders into a (nspec, nporders) array
    #  FOr longslit: read in the stanard into a (nspec, 1) array
    # read in the


    # Read in the spec1d file
    norders = len(sobjs) # ToDO: This is incorrect if you have more than one object in the sobjs
    if ret_flam:
        nspec = sobjs[0].optimal['FLAM'].size
    else:
        nspec = sobjs[0].optimal['COUNTS'].size
    # Allocate arrays and unpack spectrum
    wave = np.zeros((nspec, norders))
    #wave_mask = np.zeros((nspec, norders),dtype=bool)
    flam = np.zeros((nspec, norders))
    flam_ivar = np.zeros((nspec, norders))
    flam_mask = np.zeros((nspec, norders),dtype=bool)
    for iord in range(norders):
        wave[:,iord] = sobjs[iord].optimal['WAVE']
        #wave_mask[:,iord] = sobjs[iord].optimal['WAVE'] > 0.0
        flam_mask[:,iord] = sobjs[iord].optimal['MASK']
        if ret_flam:
            flam[:,iord] = sobjs[iord].optimal['FLAM']
            flam_ivar[:,iord] = sobjs[iord].optimal['FLAM_IVAR']
        else:
            flam[:,iord] = sobjs[iord].optimal['COUNTS']
            flam_ivar[:,iord] = sobjs[iord].optimal['COUNTS_IVAR']

    return wave, flam, flam_ivar, flam_mask


def sensfunc_guess(wave, counts_ps, inmask, flam_true, tell_dict_now, resln_guess, airmass_guess, polyorder, func,
                   lower=3.0, upper=3.0, debug=False):

    # Model parameter guess for starting the optimizations
    tell_guess = (np.median(tell_dict_now['pressure_grid']), np.median(tell_dict_now['temp_grid']),
                  np.median(tell_dict_now['h2o_grid']), airmass_guess, resln_guess, 0.0)
    tell_model1 = eval_telluric(tell_guess, tell_dict_now)
    sensguess_arg = tell_model1 * flam_true/(counts_ps + (counts_ps < 0.0))
    sensguess = np.log(sensguess_arg)
    fitmask = inmask & np.isfinite(sensguess) & (sensguess_arg > 0.0)
    # Perform an initial fit to the sensitivity function to set the starting point for optimization
    mask, coeff = utils.robust_polyfit_djs(wave, sensguess, polyorder, function=func,
                                       minx=wave.min(), maxx=wave.max(), inmask=fitmask, lower=lower, upper=upper,
                                       use_mad=True)
    sensfit_guess = np.exp(utils.func_val(coeff, wave, func, minx=wave.min(), maxx=wave.max()))

    if debug:
        plt.plot(wave, sensguess_arg)
        plt.plot(wave, sensfit_guess)
        plt.ylim(-0.1 * sensfit_guess.min(), 1.3 * sensfit_guess.max())
        plt.show()

    return coeff

def qso_guess(flux, ivar, mask, pca_dict, tell_dict, airmass, resln_guess, tell_norm_thresh):

    # Estimate the normalization for setting the bounds
    tell_guess = (np.median(tell_dict['pressure_grid']), np.median(tell_dict['temp_grid']),
                  np.median(tell_dict['h2o_grid']), airmass, resln_guess, 0.0)
    tell_model = eval_telluric(tell_guess, tell_dict)
    # Just use the mean PCA model for estimating the normalization
    pca_mean = np.exp(pca_dict['components'][0,:])
    tell_mask = tell_model > tell_norm_thresh
    # Create a reference model and bogus noise
    flux_ref = pca_mean*tell_model
    ivar_ref = utils.inverse((pca_mean/100.0)**2)
    flam_norm_inv = coadd1d.robust_median_ratio(flux, ivar, flux_ref, ivar_ref, mask=mask, mask_ref = tell_mask)
    #data_norm = np.sum(tell_mask*mask_norm*flux_med) # Data has absorption in it so we don't multilply by tell_model
    # Model has no absorption in it, so we use the average telluric absorption in the grid
    #pca_norm = np.sum(tell_mask*mask_norm*tell_model*pca_mean)
    #flam_norm = data_norm/pca_norm

    return 1.0/flam_norm_inv


def ind_lower_upper(wave_grid, mask):

    # This presumes that the data has been interpolated onto the telluric model grid
    wave_grid_ma = np.ma.array(np.copy(wave_grid))
    wave_grid_ma.mask = np.invert(mask)
    ind_lower = np.ma.argmin(wave_grid_ma)
    ind_upper = np.ma.argmax(wave_grid_ma)
    return ind_lower, ind_upper


def trim_spectrum(wave_grid, flux, ivar, mask):

    # Slice out the parts of the data that are not masked
    ind_lower, ind_upper = ind_lower_upper(wave_grid, mask)
    wave_fit = wave_grid[ind_lower:ind_upper + 1]
    flux_fit = flux[ind_lower:ind_upper + 1]
    ivar_fit = ivar[ind_lower:ind_upper + 1]
    mask_fit = mask[ind_lower:ind_upper + 1]

    return wave_fit, flux_fit, ivar_fit, mask_fit, ind_lower, ind_upper

def general_spec_reader(specfile, ret_flam=False):

    # Place holder routine that provides a generic spectrum reader

    bonus = {}
    try:
        # Read in the standard spec1d file produced by Pypeit
        sobjs, head = load.load_specobjs(specfile)
        wave, counts, counts_ivar, counts_mask = unpack_orders(sobjs, ret_flam=ret_flam)
        bonus['ECH_ORDER'] = (sobjs.ech_order).astype(int)
        bonus['ECH_ORDERINDX'] = (sobjs.ech_orderindx).astype(int)
        bonus['ECH_SNR'] = (sobjs.ech_snr).astype(float)
        bonus['NORDERS'] = wave.shape[1]
    except:
        # Read in the coadd 1d spectra file
        hdu = fits.open(specfile)
        head = hdu[0].header
        data = hdu[1].data
        wave_in, flux_in, flux_ivar_in, mask_in = data['OPT_WAVE'], data['OPT_FLAM'], data['OPT_FLAM_IVAR'], data[
            'OPT_MASK']
        wave = wave_in
        counts = flux_in
        counts_ivar = flux_ivar_in
        counts_mask = mask_in
        #wave = np.reshape(wave_in,(wave_in.size,1))
        #counts = np.reshape(flux_in,(wave_in.size,1))
        #counts_ivar = np.reshape(flux_ivar_in,(wave_in.size,1))
        #counts_mask = np.reshape(mask_in,(wave_in.size,1))

    try:
        spectrograph = load_spectrograph(head['INSTRUME'])
    except:
        # This is a hack until a generic spectrograph is implemented.
        spectrograph = load_spectrograph('shane_kast_blue')

    meta_spec = dict(core={}, bonus=bonus)
    core_keys = spectrograph.header_cards_for_spec()
    for key in core_keys:
        try:
            meta_spec['core'][key.upper()] = head[key.upper()]
        except KeyError:
            pass

    return wave, counts, counts_ivar, counts_mask, meta_spec

def interpolate_inmask(wave_grid, mask, wave_inmask, inmask):

    if inmask is not None:
        if wave_inmask is None:
            msgs.error('If you are specifying a mask you need to pass in the corresponding wavelength grid')
        # TODO we shoudld consider refactoring the interpolator to take a list of images and masks to remove the
        # the fake zero images in the call below
        _, _, inmask_int = coadd1d.interp_spec(wave_grid, wave_inmask, np.ones_like(wave_inmask), np.ones_like(wave_inmask), inmask)
        if mask.ndim == 2:
            norders = mask.shape[1]
            inmask_out = np.tile(inmask_int, (norders,1)).T
        elif mask.ndim == 1:
            inmask_out = inmask_int
        else:
            msgs.error('Unrecognized shape for data mask')
        return (mask & inmask_out)
    else:
        return mask

## TODO Add a function called telluric_fit_loop which loops over the orders/spectra to be fit. Pack the obj_dict
## and arg_dict into lists
#def telluric_fit_loop(arg_dict_list, obj_dict_list):

# TODO Make all of this a class

def sensfunc_telluric(spec1dfile, telgridfile, outfile, star_type=None, star_mag=None, star_ra=None, star_dec=None,
                      inmask=None, wave_inmask=None,
                      resln_guess=None, resln_frac_bounds=(0.5,1.5), pix_shift_bounds = (-2.0,2.0),
                      delta_coeff_bounds=(-20.0, 20.0), minmax_coeff_bounds=(-5.0, 5.0),
                      polyorder=7, func='legendre', maxiter=3, sticky=True, lower=3.0, upper=3.0, ret_flam=False,
                      seed=None, tol=1e-4, popsize=40, recombination=0.7, polish=True, disp=True, debug=False):

    """
    Loop over orders to jointly fit a sensitivity function and telluric correction for a standard star spectrum for each
    order individually


    Args:
        spec1dfile:
        telgridfile:
        star_type:
        star_mag:
        ra:
        dec:
        resln_guess:
        resln_frac_bounds:
        delta_coeff_bounds:
        polyorder:
        func:
        maxiter:
        sticky:
        use_mad:
        lower:
        upper:
        debug:
        seed:
        tol:
        popsize:
        recombination:
        disp:
        polish:

    Returns:

    """

    # Read in the telluric grid
    tell_dict = read_telluric_grid(telgridfile)
    wave_grid = tell_dict['wave_grid']
    ngrid = wave_grid.size

    # Guess resolution from wavelength sampling of telluric grid if it is not provided
    if resln_guess is None:
        resln_guess = tell_dict['resln_guess']

    # optimizer requires a seed. This guarantees that the fit will be deterministic and hence reproducible
    if seed is None:
        seed_data = np.fmin(int(np.abs(np.sum(wave_grid))), 2 ** 32 - 1)
        seed = np.random.RandomState(seed=seed_data)

    # Read in the data and interpolate data onto the regular telluric wave_grid
    wave, counts, counts_ivar, counts_mask, meta_spec = general_spec_reader(spec1dfile, ret_flam=ret_flam)
    nspec, norders = wave.shape
    counts_int, counts_ivar_int, counts_mask_int = coadd1d.interp_spec(wave_grid, wave, counts, counts_ivar, counts_mask)
    counts_ps = counts_int/meta_spec['core']['EXPTIME']
    counts_ps_ivar = counts_ivar_int*meta_spec['core']['EXPTIME']**2

    # If an input mask was provided apply it by interpolating onto wave_grid
    counts_ps_mask = interpolate_inmask(wave_grid, counts_mask_int, wave_inmask, inmask)

    # Read in standard star dictionary and interpolate onto regular telluric wave_grid
    star_ra = meta_spec['core']['RA'] if star_ra is None else star_ra
    star_dec = meta_spec['core']['DEC'] if star_dec is None else star_dec
    std_dict = flux.get_standard_spectrum(star_type=star_type, star_mag=star_mag, ra=star_ra, dec=star_dec)
    flam_true, flam_ivar, flam_mask = coadd1d.interp_spec(wave_grid, std_dict['wave'].value, std_dict['flux'].value,
                                                          np.ones_like(std_dict['flux'].value),
                                                          np.ones_like(std_dict['flux'].value, dtype=bool))
    if np.size(polyorder) > 1:
        if np.size(polyorder) != norders:
            msgs.error('polyorder must have either have norder elements or be a scalar')
        polyorder_vec = np.array(polyorder)
    else:
        polyorder_vec = np.full(norders, polyorder)

    # Allocate the meta parameter tables
    meta_table = table.Table(meta={'name': 'Parameter Values'})
    meta_table['WAVE_GRID'] = [wave_grid]
    meta_table['TOL'] = tol
    meta_table['POPSIZE'] = popsize
    meta_table['RECOMBINATION'] = recombination
    meta_table['SPEC1DFILE'] = os.path.basename(spec1dfile)
    meta_table['TELGRIDFILE'] = os.path.basename(telgridfile)
    for key,value in meta_spec['core'].items():
        meta_table[key.upper()] = value
    # Add sensitivity function specific quantities
    meta_table['STAR_TYPE'] = star_type if star_type is not None else ''
    meta_table['STAR_MAG'] = star_mag if star_mag is not None else 0.0
    meta_table['STAR_RA'] = star_ra
    meta_table['STAR_DEC'] = star_dec
    meta_table['STAR_CAL_FILE'] = std_dict['cal_file']
    meta_table['FUNCTION'] = func
    meta_table['STD_NAME'] = std_dict['name']

    out_table = table.Table(meta={'name': 'Sensfunc and Telluric Correction'})
    # TODO Check the Pypeline so that this will also work with multislit reductions
    out_table['IND_LOWER'] = np.zeros(norders, dtype=int)
    out_table['IND_UPPER'] = np.zeros(norders, dtype=int)
    out_table['WAVE_MIN'] = np.zeros(norders)
    out_table['WAVE_MAX'] = np.zeros(norders)
    out_table['TELLURIC'] = np.zeros((norders, ngrid))
    out_table['TELL_PRESS'] = np.zeros(norders)
    out_table['TELL_TEMP'] = np.zeros(norders)
    out_table['TELL_H2O'] = np.zeros(norders)
    out_table['TELL_AIRMASS'] = np.zeros(norders)
    out_table['TELL_RESLN'] = np.zeros(norders)
    out_table['TELL_SHIFT'] = np.zeros(norders)
    out_table['CHI2'] = np.zeros(norders)
    out_table['SUCCESS'] = np.zeros(norders,dtype=bool)
    out_table['NITER'] = np.zeros(norders, dtype=int)
    out_table['SENS_ORDER'] = polyorder_vec
    out_table['SENS_COEFF'] = np.zeros((norders, polyorder_vec.max() + 1))
    out_table['SENSFUNC'] = np.zeros((norders, ngrid))
    try:
        out_table['ECH_ORDER'] = meta_spec['bonus']['ECH_ORDER']
        out_table['ECH_ORDERINDX'] =  meta_spec['bonus']['ECH_ORDERINDX']
        out_table['ECH_SNR'] = meta_spec['bonus']['ECH_SNR']
    except:
        pass


    # Sort order by the strength of their telluric absorption
    srt_order_tell = sort_telluric(wave, counts_mask, tell_dict)
    wave_all_min = np.inf
    wave_all_max = -np.inf
    for iord in srt_order_tell:
        msgs.info("Fitting sensitivity function for order: {:d}/{:d}".format(iord, norders))

        # Slice out the parts of the data that are not masked, this deals with wavelenghts that are zero
        wave_fit, counts_ps_fit, counts_ps_ivar_fit, counts_ps_mask_fit, ind_lower, ind_upper = \
            trim_spectrum(wave_grid, counts_ps[:,iord], counts_ps_ivar[:,iord], counts_ps_mask[:,iord])
        flam_true_fit = flam_true[ind_lower:ind_upper + 1]
        wave_min = wave_grid[ind_lower]
        wave_max = wave_grid[ind_upper]

        # This presumes that the data has been interpolated onto the telluric model grid
        tell_dict_fit = trim_tell_dict(tell_dict, ind_lower, ind_upper)

        # Guess the coefficients by doing a preliminary fit to the sensitivity function with the average telluric behavior
        guess_obj = sensfunc_guess(wave_fit, counts_ps_fit, counts_ps_mask_fit, flam_true_fit, tell_dict_fit, resln_guess,
                                   meta_spec['core']['AIRMASS'], polyorder_vec[iord], func, lower=lower, upper=upper,
                                   debug=debug)
        # Polynomial coefficient bounds
        bounds_obj = [(np.fmin(np.abs(this_coeff)*delta_coeff_bounds[0], minmax_coeff_bounds[0]),
                         np.fmax(np.abs(this_coeff)*delta_coeff_bounds[1], minmax_coeff_bounds[1])) for this_coeff in guess_obj]

        # Set the bounds for the optimization
        bounds_tell = get_bounds_tell(tell_dict_fit, resln_guess, resln_frac_bounds, pix_shift_bounds)
        bounds = bounds_obj + bounds_tell
        # Create the arg_dict
        obj_dict = dict(wave=wave_fit, wave_min=wave_min, wave_max=wave_max, flam_true=flam_true_fit, func=func,
                        polyorder=polyorder_vec[iord])
        arg_dict = dict(ivar=counts_ps_ivar_fit, tell_dict=tell_dict_fit,
                        obj_model_func=sensfunc_obj_model, obj_dict=obj_dict,
                        bounds=bounds, seed=seed, debug=debug)
        result, ymodel, ivartot, outmask = utils.robust_optimize(counts_ps_fit, tellfit, arg_dict,
                                                                 inmask=counts_ps_mask_fit,
                                                                 maxiter=maxiter, lower=lower, upper=upper,
                                                                 sticky=sticky,
                                                                 tol=tol, popsize=popsize, recombination=recombination,
                                                                 polish=polish, disp=disp)
        # TODO move to sensfunc_telluric and individual routines
        theta_obj = result.x[:-6]
        theta_tell = result.x[-6:]
        telluric_fit = eval_telluric(theta_tell, arg_dict['tell_dict'])
        sensfit = np.exp(utils.func_val(theta_obj, wave_fit, obj_dict['func'], minx=wave_min, maxx=wave_max))
        counts_model_fit = telluric_fit * flam_true_fit/ (sensfit + (sensfit == 0.0))
        if debug:
            ## TODO Add rejected points to QA plot
            plt.plot(wave_fit, counts_ps_fit * sensfit, drawstyle='steps-mid')
            plt.plot(wave_fit, counts_ps_fit * sensfit / (telluric_fit + (telluric_fit == 0.0)), drawstyle='steps-mid')
            plt.plot(wave_fit, flam_true_fit, drawstyle='steps-mid')
            plt.ylim(-0.1 * flam_true_fit.max(), 1.5 * flam_true_fit.max())
            plt.show()

            plt.plot(wave_fit, counts_ps_fit, drawstyle='steps-mid', color='k', label='star spectrum', alpha=0.7)
            plt.plot(wave_fit, counts_model_fit, drawstyle='steps-mid', color='red', linewidth=1.0, label='model', zorder=3,
                     alpha=0.7)
            plt.ylim(-0.1 * counts_ps_fit.max(), 1.5 * counts_ps_fit.max())
            plt.legend()
            plt.show()

        out_table['CHI2'][iord] = result.fun
        out_table['SUCCESS'][iord] = result.success
        out_table['NITER'][iord] = result.nit
        out_table['IND_LOWER'][iord] = ind_lower
        out_table['IND_UPPER'][iord] = ind_upper
        out_table['WAVE_MIN'][iord] = wave_grid[ind_lower]
        out_table['WAVE_MAX'][iord] = wave_grid[ind_upper]
        out_table['SENS_COEFF'][iord][0:polyorder_vec[iord]+1] = theta_obj
        out_table['TELLURIC'][iord][ind_lower:ind_upper+1] = telluric_fit
        out_table['SENSFUNC'][iord][ind_lower:ind_upper+1] = sensfit
        out_table['TELL_PRESS'][iord] = theta_tell[0]
        out_table['TELL_TEMP'][iord] = theta_tell[1]
        out_table['TELL_H2O'][iord] = theta_tell[2]
        out_table['TELL_AIRMASS'][iord] = theta_tell[3]
        out_table['TELL_RESLN'][iord] = theta_tell[4]
        out_table['TELL_SHIFT'][iord] = theta_tell[5]
        wave_all_min = np.fmin(out_table['WAVE_MIN'][iord], wave_all_min)
        wave_all_max = np.fmax(out_table['WAVE_MAX'][iord], wave_all_max)

    meta_table['WAVE_MIN'] = wave_all_min
    meta_table['WAVE_MAX'] = wave_all_max
    # Write to outfile
    msgs.info('Writing sensitivity function and telluric to file: {:}'.format(outfile))
    hdu_param = fits.table_to_hdu(meta_table)
    hdu_table = fits.table_to_hdu(out_table)
    hdulist = fits.HDUList()
    hdulist.append(hdu_param)
    hdulist.append(hdu_table)
    hdulist.writeto(outfile, overwrite=True)


def telluric_qso(spec1dfile, telgridfile, pcafile, npca, z_qso, inmask=None, wave_inmask=None,
                 sn_cap = 30.0,
                 delta_zqso = 0.1, bounds_norm = (0.1,3.0), resln_guess=None,
                 resln_frac_bounds=(0.5,1.5), pix_shift_bounds = (-2.0,2.0),
                 tell_norm_thresh=0.9, maxiter=3, sticky=True, lower=3.0, upper=3.0,
                 seed=None, tol=1e-3, popsize=25, recombination=0.7, disp=True, polish=True,
                 debug=True, show=True):


    # Read in the telluric grid
    tell_dict = read_telluric_grid(telgridfile)
    wave_grid = tell_dict['wave_grid']
    ngrid = wave_grid.size

    # Guess resolution from wavelength sampling of telluric grid if it is not provided
    if resln_guess is None:
        resln_guess = tell_dict['resln_guess']

    # optimizer requires a seed. This guarantees that the fit will be deterministic and hence reproducible
    if seed is None:
        seed_data = np.fmin(int(np.abs(np.sum(wave_grid))), 2 ** 32 - 1)
        seed = np.random.RandomState(seed=seed_data)

    # Read in the data and interpolate data onto the regular telluric wave_grid
    wave_in, flux_in, flux_ivar_in, mask_in, meta_spec = general_spec_reader(spec1dfile, ret_flam=True)


    # Interpolate the qso onto the regular telluric wave_grid
    flux, flux_ivar, mask = coadd1d.interp_spec(wave_grid, wave_in, flux_in, flux_ivar_in, mask_in)
    # If an input mask was provided get into the wave_grid, and apply it
    mask_tot = interpolate_inmask(wave_grid, mask, wave_inmask, inmask)

    # Slice out the parts of the data that are not masked
    wave_fit, flux_fit, flux_ivar_fit1, mask_fit, ind_lower, ind_upper = trim_spectrum(wave_grid, flux, flux_ivar, mask_tot)

    # cap the inverse variance with a SN_cap to avoid excessive rejection for high S/N data. This amounts to adding
    # a tiny bit of error, 1.0/SN_CAP to the sigma of the spectrum.
    flux_ivar_fit = utils.clip_ivar(flux_fit, flux_ivar_fit1, sn_cap, mask=mask_fit)

    wave_min = wave_grid[ind_lower]
    wave_max = wave_grid[ind_upper]

    # This presumes that the data has been interpolated onto the telluric model grid
    tell_dict_fit = trim_tell_dict(tell_dict, ind_lower, ind_upper)

    #
    # fitfunc, arg_dict, bounds = instantiate_obj_model(model_type, wave_fit, model_params)
    pca_dict = qso_pca.init_pca(pcafile, wave_fit, z_qso, npca)

    flam_norm = qso_guess(flux_fit, flux_ivar_fit, mask_fit, pca_dict, tell_dict_fit, meta_spec['core']['AIRMASS'],
                          resln_guess, tell_norm_thresh)
    # Set the bounds for the PCA and truncate to the right dimension
    coeffs = pca_dict['coeffs'][:,1:npca]
    # Compute the min and max arrays of the coefficients which are not the norm, i.e. grab the coeffs that aren't the first one
    coeff_min = np.amin(coeffs, axis=0)  # only
    coeff_max = np.amax(coeffs, axis=0)
    bounds_z = [(z_qso - delta_zqso, z_qso + delta_zqso)]                # QSO redshift: can vary within delta_zqso
    bounds_flam = [(flam_norm*bounds_norm[0], flam_norm*bounds_norm[1])] # Norm: bounds determined from estimate above
    bounds_coeff = [(i, j) for i, j in zip(coeff_min, coeff_max)]        # Coefficients:  determined from PCA model
    # Set the bounds for the telluric
    bounds_tell = get_bounds_tell(tell_dict_fit, resln_guess, resln_frac_bounds, pix_shift_bounds)
    # Final bounds for the optimizaiton
    bounds =  bounds_z + bounds_flam + bounds_coeff + bounds_tell
    # Create the arg_dict
    obj_dict = dict(npca=npca, pca_dict=pca_dict)
    arg_dict = dict(ivar=flux_ivar_fit, tell_dict=tell_dict_fit,
                    obj_model_func= qso_obj_model, obj_dict=obj_dict,
                    bounds=bounds, seed=seed, debug=debug)
    result, ymodel, ivartot, outmask = utils.robust_optimize(flux_fit, tellfit, arg_dict,
                                                             inmask=mask_fit,
                                                             maxiter=maxiter, lower=lower, upper=upper,
                                                             sticky=sticky,
                                                             tol=tol, popsize=popsize, recombination=recombination,
                                                             polish=polish, disp=disp)
    theta_obj = result.x[:-6]
    theta_tell = result.x[-6:]
    # model_params, tell_params eval_obj_model(model_type, wave_fit, model_params, theta)
    pca_model, modelmask = qso_obj_model(theta_obj, obj_dict)
    telluric_fit = eval_telluric(theta_tell, arg_dict['tell_dict'])
    flux_model = telluric_fit*pca_model

    if debug:
        # Plot the data
        flux_med =scipy.ndimage.filters.median_filter(flux_fit, size=15)

        plt.plot(wave_fit, flux_fit, color='k', drawstyle='steps-mid', label='data', zorder=1, alpha=0.5)
        plt.plot(wave_fit, flux_model, color='r', drawstyle='steps-mid', label='model',zorder=2, alpha=0.5)
        plt.plot(wave_fit, pca_model, color='cornflowerblue', drawstyle='steps-mid', label='pca',zorder=3)
        plt.ylim((-0.2, 1.5*pca_model.max()))
        plt.legend()
        plt.show()
        embed()
        # Plot the telluric corrected and rescaled orders
        flux_fit_corr = flux_fit / (telluric_fit + (telluric_fit == 0.0))
        plt.plot(wave_fit, flux_fit_corr, color='k', drawstyle='steps-mid')
        plt.plot(wave_fit, pca_model, color='cornflowerblue', drawstyle='steps-mid', label='pca')
        plt.plot(wave_fit, pca_model.max()*0.9*telluric_fit, color='magenta', drawstyle='steps-mid', label='pca', alpha=0.4)
        plt.ylim((0.0, pca_model.max()*1.5))
        plt.legend()
        plt.show()

    # Interpolate the model to the grid of the input spectrum
    tell_model_out, _, _ = coadd1d.interp_spec(wave_in, wave_fit, telluric_fit, flux_ivar_fit, mask_fit)
    pca_model_out, _, _ = coadd1d.interp_spec(wave_in, wave_fit, pca_model, flux_ivar_fit, mask_fit)

    flux_out = flux_in / (tell_model_out + (tell_model_out == 0.0))
    flux_ivar_out = flux_ivar_in * tell_model_out**2
    mask_out = mask_in & (tell_model_out>0.)

    # Save the final corrected spectrum
    outfile = spec1dfile.replace('.fits', '_tell_corr.fits') # Final spectra
    head = fits.getheader(spec1dfile)
    save.save_coadd1d_to_fits(outfile, wave_in, flux_out, flux_ivar_out, mask_out,
                              header=head, ex_value='OPT', overwrite=True)

    if show:
        ymin, ymax = coadd1d.get_ylim(flux_out, flux_ivar_out, mask_out)

        fig = plt.figure(figsize=(12, 8))
        plt.plot(wave_in, flux_out, color='black', drawstyle='steps-mid',zorder=1,alpha=0.8, label='Telluric Correctted Spectrum')
        plt.plot(wave_in, np.sqrt(utils.calc_ivar(flux_ivar_out)),zorder=2, color='red', alpha=0.7,
                   drawstyle='steps-mid', linestyle=':')

        plt.ylim([ymin, ymax])
        plt.xlim([wave_in.min(), wave_in.max()])
        plt.xlabel('Wavelength ($\\rm\\AA$)')
        plt.ylabel('Flux')
        plt.show()

    # Add an analogous I/O to tables as with sensfunc_telluric
    # Allocate the output tables
    # Allocate the meta parameter tables
    meta_table = table.Table(meta={'name': 'Parameter Values'})
    meta_table['WAVE_GRID'] = [wave_grid]
    meta_table['TOL'] = tol
    meta_table['POPSIZE'] = popsize
    meta_table['RECOMBINATION'] = recombination
    meta_table['SPEC1DFILE'] = os.path.basename(spec1dfile)
    meta_table['TELGRIDFILE'] = os.path.basename(telgridfile)
    for key,value in meta_spec['core'].items():
        meta_table[key.upper()] = value

    out_table = table.Table(meta={'name': 'Sensfunc and Telluric Correction'})
    out_table['TELLURIC'] = [tell_model_out]
    out_table['TELL_PRESS'] = theta_tell[0]
    out_table['TELL_TEMP'] = theta_tell[1]
    out_table['TELL_H2O'] = theta_tell[2]
    out_table['TELL_AIRMASS'] = theta_tell[3]
    out_table['TELL_RESLN'] = theta_tell[4]
    out_table['TELL_SHIFT'] = theta_tell[5]
    out_table['PCA_THETA'] = [theta_obj]
    out_table['PCA_MODEL'] = [pca_model_out]
    out_table['CHI2'] = result.fun
    out_table['SUCCESS'] = result.success
    out_table['NITER'] = result.nit

    # Write the telluric and pca models to fits file,.
    tellfile = spec1dfile.replace('.fits', '_tell_model.fits') # Telluric and PCA model

    msgs.info('Writing telluric and PCA models to file: {:}'.format(tellfile))
    hdu_param = fits.table_to_hdu(meta_table)
    hdu_table = fits.table_to_hdu(out_table)
    hdulist = fits.HDUList()
    hdulist.append(hdu_param)
    hdulist.append(hdu_table)
    hdulist.writeto(tellfile, overwrite=True)


        # TODO Do we need software that can determine telluric absorption from a telluric star?
    #  There are two options here:
    #  1) run the sensfunc_telluric code as if it were a high-resolution staandard.
    #  2) Flux the telluric star with the sensfunc form a high-resolution standard, as if it were regular data, then
    #     run the sensfunc telluric code with a lower-order polynomial adjustment to tweak the model star spectrum
    #
    #
    # Option 2) seems like the better approach, but we need to make sure the code is compatible with this. I think it
    # would be exactly the same.

    #qso_pca_dict = dict(wave=wave, wave_mask=wave_mask, flam=flam, flam_mask=flam_mask, flam_ivar=flam_ivar,
    #                airmass=airmass,  ind=ind_tell,
    #                tell_params=tell_params, tell_model=tell_model, tell_model_orders=tell_model_orders,
    #                pca_params=pca_params, pca_model=pca_model, pca_model_orders=pca_model_orders, result = result)

#    return qso_pca_dict

