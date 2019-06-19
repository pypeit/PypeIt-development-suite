


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
    if ntheta == 6:
        tellmodel_out = shift_telluric(tellmodel_conv, tell_dict['loglam'], tell_dict['dloglam'], theta_tell[5])
        return tellmodel_out
    else:
        return tellmodel_conv


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
    counts_model = tellmodel_conv[tell_pad[0]:-tell_pad[1]]*flam_true/(sensmodel + (sensmodel == 0.0))

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


def tellfit_chi2(theta, flam, thismask, arg_dict):

    flam_ivar = arg_dict['ivar']
    flam_true = arg_dict['flam_true']
    tell_dict = arg_dict['tell_dict']
    tell_pad = tell_dict['tell_pad']
    tellmodel_conv = eval_telluric(theta, tell_dict)
    chi_vec = thismask*(tellmodel_conv[tell_pad[0]:-tell_pad[1]]*flam_true - flam)*np.sqrt(flam_ivar)
    chi2 = np.sum(np.square(chi_vec))

    return chi2



def tellfit(flam, thismask, arg_dict, **kwargs_opt):

    # Function that we are optimizing
    chi2_func= arg_dict['chi2_func']
    result = scipy.optimize.differential_evolution(chi2_func, args=(flam, thismask, arg_dict,), **kwargs_opt)
    tell_out = result.x
    tellfit_conv = eval_telluric(tell_out, arg_dict['tell_dict'])
    tell_pad = arg_dict['tell_dict']['tell_pad']
    flam_model = tellfit_conv[tell_pad[0]:-tell_pad[1]]*arg_dict['flam_true']

    return result, flam_model


def unpack_orders(sobjs, ret_flam=False):

    # TODO This should be a general reader:
    #  For echelle:  read in all the orders into a (nspec, nporders) array
    #  FOr longslit: read in the stanard into a (nspec, 1) array
    # read in the


    # Read in the spec1d file
    norders = len(sobjs)
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
    tell_pad = tell_dict_now['tell_pad']
    tell_guess = (np.median(tell_dict_now['pressure_grid']), np.median(tell_dict_now['temp_grid']),
                  np.median(tell_dict_now['h2o_grid']), airmass_guess, resln_guess, 0.0)
    tell_model1 = eval_telluric(tell_guess, tell_dict_now)
    sensguess_arg = tell_model1[tell_pad[0]:-tell_pad[1]] * flam_true/(counts_ps + (counts_ps < 0.0))
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



def fit_joint_telluric(counts_ps_in, counts_ps_ivar_in, counts_ps_mask_in, flam_true_in, tell_dict, sensfunc=True,
                       airmass=None, resln_guess=None, pix_shift_bounds = (-2.0,2.0), resln_frac_bounds=(0.5,1.5),
                       delta_coeff_bounds=(-20.0, 20.0), minmax_coeff_bounds=(-5.0, 5.0),
                       polyorder=7, func='legendre', maxiter=3, sticky=True, use_mad=False,
                       lower=3.0, upper=3.0, seed=None, debug=False, **kwargs_opt):
    """
    Jointly fit a sensitivity function and telluric correction for an input standart star spectrum.

    Args:
        wave:
        counts_ps:
        counts_ps_ivar:
        std_dict:
        tell_dict:
        inmask:
        airmass:
        resln_guess:
        resln_frac_range:
        delta_coeff:
        seed:
        polyorder:
        func:
        maxiter:
        sticky:
        use_mad:
        lower:
        upper:
        tol:
        popsize:
        recombination:
        disp:
        polish:
        debug:

    Returns:

    """

    airmass_guess = np.median(tell_dict['airmass_grid']) if airmass is None else airmass

    # This guarantees that the fit will be deterministic and hence reproducible
    if seed is None:
        seed_data = np.fmin(int(np.abs(np.sum(counts_ps_in[np.isfinite(counts_ps_in)]))), 2 ** 32 - 1)
        seed = np.random.RandomState(seed=seed_data)

    # Slice out the parts of the data that are not masked
    wave_grid = tell_dict['wave_grid']
    wave_grid_ma = np.ma.array(np.copy(wave_grid))
    wave_grid_ma.mask = np.invert(counts_ps_mask_in)
    ind_lower = np.ma.argmin(wave_grid_ma)
    ind_upper = np.ma.argmax(wave_grid_ma)
    wave = wave_grid[ind_lower:ind_upper+1]
    wave_min = wave_grid[ind_lower]
    wave_max = wave_grid[ind_upper]
    counts_ps = counts_ps_in[ind_lower:ind_upper+1]
    counts_ps_ivar = counts_ps_ivar_in[ind_lower:ind_upper+1]
    counts_ps_mask = counts_ps_mask_in[ind_lower:ind_upper+1]
    flam_true = flam_true_in[ind_lower:ind_upper+1]

    # Determine the padding and use a subset of the full tell_model_grid to make the convolutions faster
    loglam = np.log10(wave)
    dloglam = np.median(loglam[1:] - loglam[:-1])
    if resln_guess is None:
        resln_guess = 1.0/(3.0*dloglam*np.log(10.0)) # assume roughly Nyquist sampling
    pix = 1.0/resln_guess/(dloglam*np.log(10.0))/(2.0 * np.sqrt(2.0 * np.log(2))) # number of pixels per resolution element
    tell_pad = int(np.ceil(10.0 * pix))
    # This presumes that the star has been interpolated onto the telluric model grid
    ind_lower_pad = np.fmax(ind_lower - tell_pad, 0)
    ind_upper_pad = np.fmin(ind_upper + tell_pad, wave_grid.size-1)
    tell_pad_tuple = (ind_lower-ind_lower_pad, ind_upper_pad-ind_upper)
    tell_wave_grid = wave_grid[ind_lower_pad:ind_upper_pad+1]
    tell_model_grid = tell_dict['tell_grid'][:,:,:,:,ind_lower_pad:ind_upper_pad+1]
    tell_dict_now = dict(pressure_grid=tell_dict['pressure_grid'], temp_grid=tell_dict['temp_grid'],
                         h2o_grid=tell_dict['h2o_grid'], airmass_grid=tell_dict['airmass_grid'],
                         tell_grid=tell_model_grid, tell_pad=tell_pad_tuple, dloglam=dloglam,
                         loglam = np.log10(tell_wave_grid))

    if sensfunc:
        # Joint sensitivity function and telluric fits
        chi2_func = sensfunc_tellfit_chi2
        fitting_function=sensfunc_tellfit
        # Guess the coefficients by doing a fit to the sensitivity function with the average telluric behavior
        guess_coeff = sensfunc_guess(wave, counts_ps, counts_ps_mask, flam_true, tell_dict_now, resln_guess, airmass_guess,
                               polyorder, func, lower=lower, upper=upper, debug=debug)
        # Polynomial coefficient bounds
        bounds_coeff = [(np.fmin(np.abs(this_coeff)*delta_coeff_bounds[0], minmax_coeff_bounds[0]),
                         np.fmax(np.abs(this_coeff)*delta_coeff_bounds[1], minmax_coeff_bounds[1])) for this_coeff in guess_coeff]
    else:
        # Telluric only fits
        chi2_func = tellfit_chi2
        fitting_function=tellfit
        bounds_coeff = []
        guess_coeff = [] ## TODO this will break

    # Set the bounds for the optimization
    bounds_tell = [(tell_dict_now['pressure_grid'].min(), tell_dict_now['pressure_grid'].max()),
                   (tell_dict_now['temp_grid'].min(), tell_dict_now['temp_grid'].max()),
                   (tell_dict_now['h2o_grid'].min(), tell_dict_now['h2o_grid'].max()),
                   (tell_dict_now['airmass_grid'].min(), tell_dict_now['airmass_grid'].max()),
                   (resln_guess*resln_frac_bounds[0], resln_guess*resln_frac_bounds[1]),
                   pix_shift_bounds]
    bounds = bounds_coeff + bounds_tell
    bounds_tell_bar = np.mean(np.array(bounds_tell), axis=1)
    guess_tell = [bounds_tell_bar[0],
                  bounds_tell_bar[1],
                  bounds_tell_bar[2],
                  airmass_guess,
                  resln_guess,
                  0.0]
    guess = guess_coeff.tolist() + guess_tell

    arg_dict = dict(wave=wave, bounds=bounds, guess=guess, counts_ps=counts_ps, ivar=counts_ps_ivar,
                    wave_min=wave_min, wave_max=wave_max, flam_true=flam_true, tell_dict=tell_dict_now, order=polyorder,
                    sensfunc = sensfunc, func=func, chi2_func=chi2_func, seed=seed, debug=debug)
    result, ymodel, ivartot, outmask = utils.robust_optimize(counts_ps, fitting_function, arg_dict, inmask=counts_ps_mask,
                                                             maxiter=maxiter,lower=lower, upper=upper, sticky=sticky,
                                                             use_mad=use_mad, **kwargs_opt)
    #bounds = bounds, tol=tol, popsize=popsize,
    #recombination=recombination, disp=disp, polish=polish, seed=seed)

    sens_coeff = result.x[:polyorder + 1]
    tell_params = result.x[polyorder + 1:]
    tell_pad = arg_dict['tell_dict']['tell_pad']
    telluric_fit = eval_telluric(tell_params, arg_dict['tell_dict'])[tell_pad[0]:-tell_pad[1]]
    sensfit = np.exp(utils.func_val(sens_coeff, wave, arg_dict['func'], minx=wave_min, maxx=wave_max))
    counts_model = telluric_fit*arg_dict['flam_true']/(sensfit + (sensfit == 0.0))
    if debug:
        plt.plot(wave,counts_ps*sensfit, drawstyle='steps-mid')
        plt.plot(wave,counts_ps*sensfit/(telluric_fit + (telluric_fit == 0.0)), drawstyle='steps-mid')
        plt.plot(wave,flam_true, drawstyle='steps-mid')
        plt.ylim(-0.1*flam_true.max(),1.5*flam_true.max())
        plt.show()

        plt.plot(wave,counts_ps, drawstyle='steps-mid',color='k',label='star spectrum',alpha=0.7)
        plt.plot(wave,counts_model,drawstyle='steps-mid', color='red',linewidth=1.0,label='model',zorder=3,alpha=0.7)
        plt.ylim(-0.1*counts_ps.max(),1.5*counts_ps.max())
        plt.legend()
        plt.show()


    return result, tell_params, telluric_fit, sens_coeff, sensfit, ind_lower, ind_upper



def sensfunc_telluric(spec1dfile, telgridfile, outfile, star_type=None, star_mag=None, ra=None, dec=None,
                      resln_guess=None, resln_frac_bounds=(0.5,1.5), delta_coeff_bounds=(-20.0, 20.0),
                      polyorder=7, func='legendre', maxiter=3, sticky=True, use_mad=False, lower=3.0, upper=3.0,
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
    tell_model_dict = read_telluric_grid(telgridfile)
    wave_grid = tell_model_dict['wave_grid']
    ngrid = wave_grid.size

    # Read in the standard star spectrum and interpolate it onto the regular telluric wave grid.

    sobjs, head = load.load_specobjs(spec1dfile)
    wave, counts, counts_ivar, counts_mask = unpack_orders(sobjs)
    exptime = head['EXPTIME']
    airmass = head['AIRMASS']



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

    nspec, norders = wave.shape
    if np.size(polyorder) > 1:
        if np.size(polyorder) != norders:
            msgs.error('polyorder must have either have norder elements or be a scalar')
        polyorder_vec = np.array(polyorder)
    else:
        polyorder_vec = np.full(norders, polyorder)

    # This guarantees that the fit will be deterministic and hence reproducible
    if seed is None:
        seed_data = np.fmin(int(np.abs(np.sum(std_dict['flux'].value))), 2 ** 32 - 1)
        seed = np.random.RandomState(seed=seed_data)

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
    meta_table['SPEC1DFILE'] = spec1dfile
    meta_table['CAL_FILE'] = std_dict['cal_file']
    meta_table['STD_NAME'] = std_dict['name']


    out_table = table.Table(meta={'name': 'Sensfunc and Telluric Correction'})
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
    srt_order_tell = sort_telluric(wave, counts_mask, tell_model_dict)
    wave_all_min = np.inf
    wave_all_max = -np.inf
    for iord in srt_order_tell:
        msgs.info("Fitting sensitivity function for order: {:d}/{:d}".format(iord, norders))
        result, tell_params, tellfit, sens_coeff, sensfit, ind_lower, ind_upper = fit_joint_telluric(
            counts_ps[:, iord], counts_ps_ivar[:, iord], counts_ps_mask[:, iord], flam_true, tell_model_dict,
            airmass=airmass, resln_guess=resln_guess, resln_frac_bounds=resln_frac_bounds, delta_coeff_bounds=delta_coeff_bounds,
            polyorder=polyorder_vec[iord], func=func, maxiter=maxiter, sticky=sticky, use_mad=use_mad, lower=lower, upper=upper,
            seed = seed, debug=debug, tol=1e-4, popsize=40, recombination=0.7, polish=polish, disp=disp)
        out_table['CHI2'][iord] = result.fun
        out_table['SUCCESS'][iord] = result.success
        out_table['NITER'][iord] = result.nit
        out_table['IND_LOWER'][iord] = ind_lower
        out_table['IND_UPPER'][iord] = ind_upper
        out_table['WAVE_MIN'][iord] = wave_grid[ind_lower]
        out_table['WAVE_MAX'][iord] = wave_grid[ind_upper]
        out_table['SENS_COEFF'][iord][0:polyorder_vec[iord]+1] = sens_coeff
        out_table['TELLURIC'][iord][ind_lower:ind_upper+1] = tellfit
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

