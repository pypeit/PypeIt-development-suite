


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
import IPython
import qso_pca


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

    tell_dict = dict(wave_grid=wave_grid, pressure_grid=pg, temp_grid=tg, h2o_grid=hg, airmass_grid=ag,
                     tell_grid=model_grid)
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

# QSO TELLFIT
def qso_tellfit_eval(theta, arg_dict):

    tell_dict = arg_dict['tell_dict']
    npca = arg_dict['npca']

    # There are npca-1 coefficients, norm and redshift for npca+1 parameters
    theta_PCA = theta[:npca + 1]
    # Telluric is 6-d for general case of airmass, temp, pressure, H2O, resln, shift
    theta_tell = theta[-6:]
    tell_model = eval_telluric(theta_tell, arg_dict['tell_dict'])
    pca_model = qso_pca.pca_eval(theta_PCA, arg_dict['pca_dict'])
    # TODO Is the prior evaluation slowing things down??
    # TODO Disablingthe prior for now as I think it slows things down for no big gain
    #ln_pca_pri = qso_pca.pca_lnprior(theta_PCA, arg_dict['pca_dict'])
    ln_pca_pri = 0.0
    return tell_model, pca_model, ln_pca_pri

def qso_tellfit_chi2(theta, flux, thismask, arg_dict):


    tell_model, pca_model, ln_pca_pri = qso_tellfit_eval(theta, arg_dict)
    chi_vec = thismask*(flux - tell_model*pca_model)*np.sqrt(arg_dict['flux_ivar'])

    robust_scale = 2.0
    huber_vec = scipy.special.huber(robust_scale, chi_vec)
    chi2_remap = np.sum(np.square(huber_vec * thismask))

    lnL = -chi2_remap/2.0
    lnptot = lnL + ln_pca_pri
    chi2_tot = -2.0*lnptot
    return chi2_tot


def qso_tellfit(flux, thismask, arg_dict, **kwargs_opt):

    # Function that we are optimizing
    chi2_func = arg_dict['chi2_func']
    flux_ivar = arg_dict['flux_ivar']
    bounds = arg_dict['bounds']
    seed = arg_dict['seed']
    result = scipy.optimize.differential_evolution(chi2_func, bounds, args=(flux, thismask, arg_dict,), seed=seed,
                                                   **kwargs_opt)

    tell_model, pca_model, ln_pca_pri = qso_tellfit_eval(result.x, arg_dict)
    chi_vec = thismask*(flux - tell_model*pca_model)*np.sqrt(flux_ivar)
    try:
        debug = arg_dict['debug']
    except KeyError:
        debug = False

    sigma_corr, maskchi = coadd1d.renormalize_errors(chi_vec, mask=thismask, title='qso_tellfit', debug=debug)
    ivartot = flux_ivar/sigma_corr ** 2

    return result, tell_model*pca_model, ivartot

## SENSFUNC_TELLFIT
def sensfunc_tellfit_eval(theta, arg_dict):

    wave_star = arg_dict['wave']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    flam_true = arg_dict['flam_true']
    tell_dict = arg_dict['tell_dict']
    tell_pad = tell_dict['tell_pad']
    order = arg_dict['order']
    func = arg_dict['func']

    theta_sens = theta[:order+1]
    theta_tell = theta[order+1:]
    sensmodel = np.exp(utils.func_val(theta_sens, wave_star, func, minx=wave_min, maxx=wave_max))
    tellmodel_conv = eval_telluric(theta_tell, tell_dict)
    counts_model = tellmodel_conv*flam_true/(sensmodel + (sensmodel == 0.0))

    return counts_model, sensmodel

def sensfunc_tellfit_chi2(theta, counts_ps, thismask, arg_dict):
    """

    Args:
        theta:
        counts_ps:
        thismask:
        arg_dict:

    Returns:

    """

    counts_model, sensmodel = sensfunc_tellfit_eval(theta, arg_dict)
    counts_ps_ivar = arg_dict['ivar']

    if np.sum(np.abs(sensmodel)) < 1e-6:
        return np.inf
    else:
        chi_vec = thismask*(sensmodel != 0.0)*(counts_model - counts_ps)*np.sqrt(counts_ps_ivar)
        # Robustly characterize the dispersion of this distribution
        #chi_mean, chi_median, chi_std = stats.sigma_clipped_stats(chi_vec, np.invert(thismask), cenfunc='median',
        #                                                          stdfunc=stats.mad_std, maxiters=5, sigma=2.0)
        robust_scale = 2.0
        huber_vec = scipy.special.huber(robust_scale, chi_vec)
        loss_function = np.sum(np.square(huber_vec * thismask))
        return loss_function


def sensfunc_tellfit(counts_ps, thismask, arg_dict, **kwargs_opt):

    # Function that we are optimizing
    chi2_func = arg_dict['chi2_func']
    counts_ps_ivar = arg_dict['ivar']
    #guess = arg_dict['guess']
    bounds = arg_dict['bounds']
    seed = arg_dict['seed']
    #result = scipy.optimize.minimize(chi2_func, guess, bounds = bounds, method='TNC', args=(counts_ps, thismask, arg_dict),  **kwargs_opt)
    result = scipy.optimize.differential_evolution(chi2_func, bounds, args=(counts_ps, thismask, arg_dict,), seed=seed,
                                                   **kwargs_opt)

    counts_model, sensmodel = sensfunc_tellfit_eval(result.x, arg_dict)
    chi_vec = thismask * (sensmodel != 0.0) * (counts_model - counts_ps) * np.sqrt(counts_ps_ivar)

    try:
        debug = arg_dict['debug']
    except KeyError:
        debug = False

    sigma_corr, maskchi = coadd1d.renormalize_errors(chi_vec, mask=thismask, title = 'sensfunc_tellfit', debug=debug)
    ivartot = counts_ps_ivar/sigma_corr**2

    return result, counts_model, ivartot


def get_bounds_tell(tell_dict, resln_guess, resln_frac_bounds, pix_shift_bounds):

    # Set the bounds for the optimization
    bounds_tell = [(tell_dict['pressure_grid'].min(), tell_dict['pressure_grid'].max()),
                   (tell_dict['temp_grid'].min(), tell_dict['temp_grid'].max()),
                   (tell_dict['h2o_grid'].min(), tell_dict['h2o_grid'].max()),
                   (tell_dict['airmass_grid'].min(), tell_dict['airmass_grid'].max()),
                   (resln_guess * resln_frac_bounds[0], resln_guess * resln_frac_bounds[1]),
                   pix_shift_bounds]

    return bounds_tell

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

def qso_norm_guess(pca_dict, tell_dict, airmass, resln_guess, tell_norm_thresh, flux, ivar, mask):

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

## Not used
def get_dloglam_data(wave):

    loglam = np.log10(wave)
    loglam_ma = np.ma.array(np.copy(loglam))
    loglam_ma.mask = np.invert(wave > 10.0)
    dloglam_arr = np.abs(loglam_ma - np.roll(loglam_ma, 1, axis=0))[1:, :]
    dloglam = np.ma.median(dloglam_arr)

    return dloglam


def trim_tell_dict(tell_dict, ind_lower, ind_upper, dloglam, tell_pad):

    wave_grid = tell_dict['wave_grid']
    ind_lower_pad = np.fmax(ind_lower - tell_pad, 0)
    ind_upper_pad = np.fmin(ind_upper + tell_pad, wave_grid.size - 1)
    tell_pad_tuple = (ind_lower - ind_lower_pad, ind_upper_pad - ind_upper)
    tell_wave_grid = wave_grid[ind_lower_pad:ind_upper_pad + 1]
    tell_model_grid = tell_dict['tell_grid'][:, :, :, :, ind_lower_pad:ind_upper_pad + 1]
    tell_dict_fit = dict(pressure_grid=tell_dict['pressure_grid'], temp_grid=tell_dict['temp_grid'],
                         h2o_grid=tell_dict['h2o_grid'], airmass_grid=tell_dict['airmass_grid'],
                         tell_grid=tell_model_grid, tell_pad=tell_pad_tuple, dloglam=dloglam,
                         loglam=np.log10(tell_wave_grid))

    return tell_dict_fit

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



def sensfunc_telluric(spec1dfile, telgridfile, outfile, star_type=None, star_mag=None, ra=None, dec=None,
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

    # Read in the standard star spectrum and interpolate it onto the regular telluric wave grid.
    ## ToDo: currently just using try except to deal with different format. We should fix the data model
    try:
        # Read in the standard spec1d file produced by Pypeit
        sobjs, head = load.load_specobjs(spec1dfile)
        wave, counts, counts_ivar, counts_mask = unpack_orders(sobjs, ret_flam=ret_flam)
    except:
        # Read in the coadd 1d spectra file
        hdu = fits.open(spec1dfile)
        head = hdu[0].header
        data = hdu[1].data
        wave_in, flux_in, flux_ivar_in, mask_in = data['OPT_WAVE'], data['OPT_FLAM'], data['OPT_FLAM_IVAR'], data[
            'OPT_MASK']
        wave = np.reshape(wave_in,(wave_in.size,1))
        counts = np.reshape(flux_in,(wave_in.size,1))
        counts_ivar = np.reshape(flux_ivar_in,(wave_in.size,1))
        counts_mask = np.reshape(mask_in,(wave_in.size,1))

    exptime = head['EXPTIME']
    airmass = head['AIRMASS']

    loglam = np.log10(wave_grid)
    dloglam = np.median(loglam[1:] - loglam[:-1])
    # Guess resolution from wavelength sampling of telluric grid if it is not provided
    if resln_guess is None:
        resln_guess = 1.0/(3.0 * dloglam * np.log(10.0))  # assume roughly Nyquist sampling

    pix_per_R = 1.0 / resln_guess / (dloglam * np.log(10.0)) / (2.0 * np.sqrt(2.0 * np.log(2)))
    tell_pad_pix = int(np.ceil(10.0 * pix_per_R))

    # Interpolate the data and standard star spectrum onto the regular telluric wave_grid
    counts_int, counts_ivar_int, counts_mask_int = coadd1d.interp_spec(wave_grid, wave, counts, counts_ivar, counts_mask)
    counts_ps = counts_int/exptime
    counts_ps_ivar = counts_ivar_int* exptime ** 2
    counts_ps_mask = counts_mask_int

    # Read in standard star dictionary
    std_dict = flux.get_standard_spectrum(star_type=star_type, star_mag=star_mag, ra=ra, dec=dec)
    flam_true, flam_ivar, flam_mask = coadd1d.interp_spec(wave_grid, std_dict['wave'].value, std_dict['flux'].value,
                                                          0.0*std_dict['flux'].value,
                                                          np.ones_like(std_dict['flux'].value, dtype=bool))

    # This guarantees that the fit will be deterministic and hence reproducible
    if seed is None:
        seed_data = np.fmin(int(np.abs(np.sum(std_dict['flux'].value))), 2 ** 32 - 1)
        seed = np.random.RandomState(seed=seed_data)


    nspec, norders = wave.shape
    if np.size(polyorder) > 1:
        if np.size(polyorder) != norders:
            msgs.error('polyorder must have either have norder elements or be a scalar')
        polyorder_vec = np.array(polyorder)
    else:
        polyorder_vec = np.full(norders, polyorder)

    # Allocate the output tables
    meta_table = table.Table(meta={'name': 'Parameter Values'})
    meta_table['WAVE_GRID'] = [wave_grid]
    meta_table['TOL'] = tol
    meta_table['POPSIZE'] = popsize
    meta_table['RECOMBINATION'] = recombination
    meta_table['STAR_TYPE'] = star_type if star_type is not None else ''
    meta_table['STAR_MAG'] = star_mag if star_mag is not None else 0.0
    meta_table['STAR_RA'] = ra if ra is not None else 0.0
    meta_table['STAR_DEC'] = dec if dec is not None else 0.0
    meta_table['FUNCTION'] = func
    meta_table['PYPELINE'] = head['PYPELINE']
    meta_table['EXPTIME'] = head['EXPTIME']
    meta_table['AIRMASS'] = head['AIRMASS']
    meta_table['NORDERS'] = norders
    meta_table['CAL_FILE'] = std_dict['cal_file']
    meta_table['STD_NAME'] = std_dict['name']
    meta_table['TELGRIDFILE'] = os.path.basename(telgridfile)
    meta_table['SPEC1DFILE'] = os.path.basename(spec1dfile)


    out_table = table.Table(meta={'name': 'Sensfunc and Telluric Correction'})
    # TODO Check the Pypeline so that this will also work with multislit reductions
    out_table['ECH_ORDER'] = (sobjs.ech_order).astype(int)
    out_table['ECH_ORDERINDX'] = (sobjs.ech_orderindx).astype(int)
    out_table['ECH_SNR'] = (sobjs.ech_snr).astype(float)
    out_table['IND_LOWER'] = np.zeros(norders, dtype=int)
    out_table['IND_UPPER'] = np.zeros(norders, dtype=int)
    out_table['WAVE_MIN'] = np.zeros(norders)
    out_table['WAVE_MAX'] = np.zeros(norders)
    out_table['SENS_ORDER'] = polyorder_vec
    out_table['SENS_COEFF'] = np.zeros((norders, polyorder_vec.max() + 1))
    out_table['TELLURIC'] = np.zeros((norders, ngrid))
    out_table['SENSFUNC'] = np.zeros((norders, ngrid))
    out_table['TELL_PRESS'] = np.zeros(norders)
    out_table['TELL_TEMP'] = np.zeros(norders)
    out_table['TELL_H2O'] = np.zeros(norders)
    out_table['TELL_AIRMASS'] = np.zeros(norders)
    out_table['TELL_RESLN'] = np.zeros(norders)
    out_table['TELL_SHIFT'] = np.zeros(norders)
    out_table['CHI2'] = np.zeros(norders)
    out_table['SUCCESS'] = np.zeros(norders,dtype=bool)
    out_table['NITER'] = np.zeros(norders, dtype=int)

    # Sort order by the strength of their telluric absorption
    srt_order_tell = sort_telluric(wave, counts_mask, tell_dict)
    wave_all_min = np.inf
    wave_all_max = -np.inf
    for iord in srt_order_tell:
        msgs.info("Fitting sensitivity function for order: {:d}/{:d}".format(iord, norders))

        # Slice out the parts of the data that are not masked, this deals with wavelenghth that are zero
        wave_fit, counts_ps_fit, counts_ps_ivar_fit, counts_ps_mask_fit, ind_lower, ind_upper = \
            trim_spectrum(wave_grid, counts_ps[:,iord], counts_ps_ivar[:,iord], counts_ps_mask[:,iord])
        flam_true_fit = flam_true[ind_lower:ind_upper + 1]
        wave_min = wave_grid[ind_lower]
        wave_max = wave_grid[ind_upper]

        # This presumes that the data has been interpolated onto the telluric model grid
        tell_dict_fit = trim_tell_dict(tell_dict, ind_lower, ind_upper, dloglam, tell_pad_pix)

        # Guess the coefficients by doing a preliminary fit to the sensitivity function with the average telluric behavior
        guess_coeff = sensfunc_guess(wave_fit, counts_ps_fit, counts_ps_mask_fit, flam_true_fit, tell_dict_fit, resln_guess, airmass,
                               polyorder_vec[iord], func, lower=lower, upper=upper, debug=debug)
        # Polynomial coefficient bounds
        bounds_coeff = [(np.fmin(np.abs(this_coeff)*delta_coeff_bounds[0], minmax_coeff_bounds[0]),
                         np.fmax(np.abs(this_coeff)*delta_coeff_bounds[1], minmax_coeff_bounds[1])) for this_coeff in guess_coeff]

        # Set the bounds for the optimization
        bounds_tell = get_bounds_tell(tell_dict_fit, resln_guess, resln_frac_bounds, pix_shift_bounds)
        bounds = bounds_coeff + bounds_tell
        # Create the arg_dict
        arg_dict = dict(wave=wave_fit, counts_ps=counts_ps_fit, ivar=counts_ps_ivar_fit,
                        wave_min=wave_min, wave_max=wave_max, flam_true=flam_true_fit,
                        tell_dict=tell_dict_fit,  bounds=bounds,order=polyorder_vec[iord], func=func,
                        chi2_func = sensfunc_tellfit_chi2, seed=seed, debug=debug)
        result, ymodel, ivartot, outmask = utils.robust_optimize(counts_ps_fit, sensfunc_tellfit, arg_dict,
                                                                 inmask=counts_ps_mask_fit,
                                                                 maxiter=maxiter, lower=lower, upper=upper,
                                                                 sticky=sticky,
                                                                 tol=tol, popsize=popsize, recombination=recombination,
                                                                 polish=polish, disp=disp)
        # TODO move to sensfunc_telluric and individual routines
        sens_coeff = result.x[:polyorder_vec[iord] + 1]
        tell_params = result.x[polyorder_vec[iord] + 1:]
        tell_pad_tuple = arg_dict['tell_dict']['tell_pad']
        telluric_fit = eval_telluric(tell_params, arg_dict['tell_dict'])
        sensfit = np.exp(utils.func_val(sens_coeff, wave_fit, arg_dict['func'], minx=wave_min, maxx=wave_max))
        counts_model_fit = telluric_fit * flam_true_fit/ (sensfit + (sensfit == 0.0))
        IPython.embed()
        if debug:
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
        out_table['SENS_COEFF'][iord][0:polyorder_vec[iord]+1] = sens_coeff
        out_table['TELLURIC'][iord][ind_lower:ind_upper+1] = telluric_fit
        out_table['SENSFUNC'][iord][ind_lower:ind_upper+1] = sensfit
        out_table['TELL_PRESS'][iord] = tell_params[0]
        out_table['TELL_TEMP'][iord] = tell_params[1]
        out_table['TELL_H2O'][iord] = tell_params[2]
        out_table['TELL_AIRMASS'][iord] = tell_params[3]
        out_table['TELL_RESLN'][iord] = tell_params[4]
        out_table['TELL_SHIFT'][iord] = tell_params[5]
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

    # Read in the data
    # TODO we need to sort out our data model!!
    hdu = fits.open(spec1dfile)
    head = hdu[0].header
    data = hdu[1].data
    wave_in, flux_in, flux_ivar_in, mask_in = data['OPT_WAVE'], data['OPT_FLAM'], data['OPT_FLAM_IVAR'], data['OPT_MASK']
    mask_in = mask_in.astype(bool) ## TODO hack fix this
    airmass = head['AIRMASS']

    # This guarantees that the fit will be deterministic and hence reproducible
    if seed is None:
        seed_data = np.fmin(int(np.abs(np.sum(flux_in))), 2 ** 32 - 1)
        seed = np.random.RandomState(seed=seed_data)

    # Read in the telluric grid
    tell_dict = read_telluric_grid(telgridfile)
    wave_grid = tell_dict['wave_grid']
    ngrid = wave_grid.size

    loglam = np.log10(wave_grid)
    dloglam = np.median(loglam[1:] - loglam[:-1])
    # Guess resolution from wavelength sampling of telluric grid if it is not provided
    if resln_guess is None:
        resln_guess = 1.0/(3.0 * dloglam * np.log(10.0))  # assume roughly Nyquist sampling

    pix_per_R = 1.0 / resln_guess / (dloglam * np.log(10.0)) / (2.0 * np.sqrt(2.0 * np.log(2)))
    tell_pad_pix = int(np.ceil(10.0 * pix_per_R))

    # Interpolate the qso onto the regular telluric wave_grid
    flux, flux_ivar, mask = coadd1d.interp_spec(wave_grid, wave_in, flux_in, flux_ivar_in, mask_in)
    # If an input mask was provided get into the wave_grid, and apply it
    if inmask is not None:
        if wave_inmask is None:
            msgs.error('If you are specifying a mask you need to pass in the corresponding wavelength grid')
        # TODO we shoudld consider refactoring the interpolator to take a list of images and masks to remove the
        # the fake zero images in the call below
        _, _, inmask_int = coadd1d.interp_spec(wave_grid, wave_inmask,
                                               np.ones_like(wave_inmask), np.ones_like(wave_inmask), inmask)
        mask_tot = mask & inmask_int
    else:
        mask_tot = mask

    # Slice out the parts of the data that are not masked
    wave_fit, flux_fit, flux_ivar_fit1, mask_fit, ind_lower, ind_upper = trim_spectrum(wave_grid, flux, flux_ivar, mask_tot)

    # cap the inverse variance with a SN_cap to avoid excessive rejection for high S/N data. This amounts to adding
    # a tiny bit of error, 1.0/SN_CAP to the sigma of the spectrum.
    flux_ivar_fit = utils.cap_ivar(flux_fit, flux_ivar_fit1, sn_cap, mask=mask_fit)

    wave_min = wave_grid[ind_lower]
    wave_max = wave_grid[ind_upper]

    # This presumes that the data has been interpolated onto the telluric model grid
    tell_dict_fit = trim_tell_dict(tell_dict, ind_lower, ind_upper, dloglam, tell_pad_pix)

    #
    # fitfunc, arg_dict, bounds = instantiate_obj_model(model_type, wave_fit, model_params)
    pca_dict = qso_pca.init_pca(pcafile, wave_fit, z_qso, npca)

    flam_norm = qso_norm_guess(pca_dict, tell_dict_fit, airmass, resln_guess, tell_norm_thresh, flux_fit, flux_ivar_fit, mask_fit)
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
    arg_dict = dict(npca=npca, flux_ivar=flux_ivar_fit, tell_dict=tell_dict_fit, pca_dict=pca_dict, debug=debug,
                    bounds=bounds, seed=seed, chi2_func = qso_tellfit_chi2)
    result, ymodel, ivartot, outmask = utils.robust_optimize(flux_fit, qso_tellfit, arg_dict,
                                                             inmask=mask_fit,
                                                             maxiter=maxiter, lower=lower, upper=upper,
                                                             sticky=sticky,
                                                             tol=tol, popsize=popsize, recombination=recombination,
                                                             polish=polish, disp=disp)
    theta = result.x
    # model_params, tell_params eval_obj_model(model_type, wave_fit, model_params, theta)
    tell_model, pca_model, ln_pca_pri = qso_tellfit_eval(theta, arg_dict)
    flux_model = tell_model*pca_model
    pca_params = theta[:npca + 1]
    tell_params = theta[-6:]


    if debug:
        # Plot the data
        flux_med =scipy.ndimage.filters.median_filter(flux_fit, size=15)

        plt.plot(wave_fit, flux_fit, color='k', drawstyle='steps-mid', label='data', zorder=1, alpha=0.5)
        plt.plot(wave_fit, flux_model, color='r', drawstyle='steps-mid', label='model',zorder=2, alpha=0.5)
        plt.plot(wave_fit, pca_model, color='cornflowerblue', drawstyle='steps-mid', label='pca',zorder=3)
        plt.ylim((-0.2, 1.5*pca_model.max()))
        plt.legend()
        plt.show()
        IPython.embed()
        # Plot the telluric corrected and rescaled orders
        flux_fit_corr = flux_fit / (tell_model + (tell_model == 0.0))
        plt.plot(wave_fit, flux_fit_corr, color='k', drawstyle='steps-mid')
        plt.plot(wave_fit, pca_model, color='cornflowerblue', drawstyle='steps-mid', label='pca')
        plt.plot(wave_fit, pca_model.max()*0.9*tell_model, color='magenta', drawstyle='steps-mid', label='pca', alpha=0.4)
        plt.ylim((0.0, pca_model.max()*1.5))
        plt.legend()
        plt.show()

    # Interpolate the model to the grid of raw spectrum
    tell_model_out, _, _ = coadd1d.interp_spec(wave_in, wave_fit, tell_model, flux_ivar_fit, mask_fit)
    pca_model_out, _, _ = coadd1d.interp_spec(wave_in, wave_fit, pca_model, flux_ivar_fit, mask_fit)

    flux_out = flux_in / (tell_model_out + (tell_model_out == 0.0))
    flux_ivar_out = flux_ivar_in * tell_model_out**2
    mask_out = mask_in & (tell_model_out>0.)

    # Save the final corrected spectrum
    outfile = spec1dfile.replace('.fits', '_tell_corr.fits') # Final spectra
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
    meta_table = table.Table(meta={'name': 'Parameter Values'})
    meta_table['WAVE_GRID'] = [wave_in]
    meta_table['TOL'] = tol
    meta_table['POPSIZE'] = popsize
    meta_table['RECOMBINATION'] = recombination
    meta_table['PYPELINE'] = head['PYPELINE']
    meta_table['EXPTIME'] = head['EXPTIME']
    meta_table['AIRMASS'] = head['AIRMASS']
    meta_table['TELGRIDFILE'] = os.path.basename(telgridfile)
    meta_table['SPEC1DFILE'] = os.path.basename(spec1dfile)

    out_table = table.Table(meta={'name': 'Sensfunc and Telluric Correction'})
    out_table['TELLURIC'] = [tell_model_out]
    out_table['TELL_PRESS'] = tell_params[0]
    out_table['TELL_TEMP'] = tell_params[1]
    out_table['TELL_H2O'] = tell_params[2]
    out_table['TELL_AIRMASS'] = tell_params[3]
    out_table['TELL_RESLN'] = tell_params[4]
    out_table['TELL_SHIFT'] = tell_params[5]
    out_table['PCA_THETA'] = [pca_params]
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

