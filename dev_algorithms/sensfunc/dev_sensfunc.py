


import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from sklearn import mixture
from astropy.io import fits
from pypeit.core import pydl
import astropy.units as u
from astropy.io import fits
from pypeit.core import flux
from pypeit.core import load
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

def init_pca(filename,wave_grid,redshift, npca):
    # Read in the pickle file from coarse_pca.create_coarse_pca
    # The relevant pieces are the wavelengths (wave_pca_c), the PCA components (pca_comp_c),
    # and the Gaussian mixture model prior (mix_fit)

    loglam = np.log10(wave_grid)
    dloglam = np.median(loglam[1:] - loglam[:-1])
    wave_pca_c, cont_all_c, pca_comp_c, coeffs_c, mean_pca, covar_pca, diff_pca, mix_fit, chi2, dof = pickle.load(open(filename,'rb'))
    num_comp = pca_comp_c.shape[0] # number of PCA components
    # Interpolate PCA components onto wave_grid
    pca_interp = scipy.interpolate.interp1d(wave_pca_c*(1+redshift),pca_comp_c, bounds_error=False, fill_value=0.0, axis=1)
    pca_comp_new = pca_interp(wave_grid)
    # Generate a mixture model for the coefficients prior, what should ngauss be?
    prior = mixture.GaussianMixture(n_components = npca-1).fit(coeffs_c[:, 1:npca])
    # Construct the PCA dict
    pca_dict = {'npca': npca, 'components': pca_comp_new, 'prior': prior, 'coeffs': coeffs_c,
                'z_fid': redshift, 'dloglam': dloglam}
    return pca_dict

def pca_eval(theta,pca_dict):
    C = pca_dict['components']
    z_fid = pca_dict['z_fid']
    dloglam = pca_dict['dloglam']
    npca = pca_dict['npca']  # Size of the PCA currently being used, original PCA in the dict could be larger
    z_qso = theta[0]
    norm = theta[1]
    A = theta[2:]
    dshift = int(np.round(np.log10((1.0 + z_qso)/(1.0 + z_fid))/dloglam))
    C_now = np.roll(C[:npca,:], dshift, axis=1)
    return norm*np.exp(np.dot(np.append(1.0,A),C_now))

def pca_lnprior(theta,pca_dict):
    gaussian_mixture_model = pca_dict['prior']
    A = theta[2:]
    return gaussian_mixture_model.score_samples(A.reshape(1,-1))

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

def update_bounds(bounds, delta_coeff, coeff):

    bounds_new = [(this_coeff * delta_coeff[0], this_coeff * delta_coeff[1]) for this_coeff in coeff]
    bounds_tell = bounds[len(coeff):]
    bounds_new.extend(bounds_tell)
    return bounds_new



def sensfunc_tellfit_chi2(theta, counts_ps, thismask, arg_dict):
    """

    Args:
        theta:
        counts_ps:
        thismask:
        arg_dict:

    Returns:

    """
    wave_star = arg_dict['wave']
    counts_ps_ivar = arg_dict['ivar']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    flam_true = arg_dict['flam_true']
    tell_dict = arg_dict['tell_dict']
    tell_pad = tell_dict['tell_pad']
    order = arg_dict['order']
    func = arg_dict['func']

    theta_sens = theta[:order+1]
    theta_tell = theta[order+1:]
    sensmodel = utils.func_val(theta_sens, wave_star, func, minx=wave_min, maxx=wave_max)
    tellmodel_conv = eval_telluric(theta_tell, tell_dict)
    if np.sum(np.abs(sensmodel)) < 1e-6:
        return np.inf
    else:
        chi_vec = thismask*(sensmodel != 0.0)*(tellmodel_conv[tell_pad[0]:-tell_pad[1]]*flam_true/(sensmodel + (sensmodel == 0.0)) -
                                               counts_ps)*np.sqrt(counts_ps_ivar)
        chi2 = np.sum(np.square(chi_vec))

    return chi2


def sensfunc_tellfit(counts_ps, thismask, arg_dict, **kwargs_opt):

    # Function that we are optimizing
    chi2_func = arg_dict['chi2_func']
    result = scipy.optimize.differential_evolution(chi2_func, args=(counts_ps, thismask, arg_dict,), **kwargs_opt)
    wave_star = arg_dict['wave']
    order = arg_dict['order']
    coeff_out = result.x[:order+1]
    tell_out = result.x[order+1:]
    tellfit_conv = eval_telluric(tell_out, arg_dict['tell_dict'])
    tell_pad = arg_dict['tell_dict']['tell_pad']
    sensfit = utils.func_val(coeff_out, wave_star, arg_dict['func'], minx=arg_dict['wave_min'], maxx=arg_dict['wave_max'])
    counts_model = tellfit_conv[tell_pad[0]:-tell_pad[1]]*arg_dict['flam_true']/(sensfit + (sensfit == 0.0))

    return result, counts_model


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


def tell_qso_evalmodel(theta, arg_dict):

    npca = arg_dict['npca']
    ind = arg_dict['ind']
    theta_PCA = theta[:npca + 1]
    theta_tell = theta[-5:]
    tell_model = eval_telluric(theta_tell, arg_dict['tell_dict'])
    tell_model_orders = populate_orders(tell_model, ind)
    pca_model = pca_eval(theta_PCA, arg_dict['pca_dict'])
    pca_model_orders = populate_orders(pca_model, ind)
    flam_model = tell_model_orders*pca_model_orders
    return flam_model

def tell_qso_chi2(theta, flam, thismask, arg_dict):

    #norders = arg_dict['norders']
    npca = arg_dict['npca']
    theta_PCA = theta[:npca + 1]
    flam_model = tell_qso_evalmodel(theta, arg_dict)
    chi_vec = thismask*(flam - flam_model)*np.sqrt(arg_dict['flam_ivar'])
    chi2 = np.sum(np.square(chi_vec))
    lnL = -chi2/2.0
    lnpri = pca_lnprior(theta_PCA, arg_dict['pca_dict'])
    lnptot = lnL + lnpri
    chi2_tot = -2.0*lnptot
    return chi2_tot

def tell_qso_fit(flam, thismask, arg_dict, **kwargs_opt):

    # Function that we are optimizing
    tell_qso_chi2 = arg_dict['tell_qso_chi2']
    result = scipy.optimize.differential_evolution(tell_qso_chi2, args=(flam, thismask, arg_dict,), **kwargs_opt)
    flam_model = tell_qso_evalmodel(result.x, arg_dict)
    return result, flam_model


def unpack_orders(spec1dfile, ret_flam=False):

    # Read in the spec1d file
    sobjs, head = load.load_specobjs(spec1dfile)
    norders = len(sobjs)
    nspec = sobjs[0].optimal['COUNTS'].size
    # Allocate arrays and unpack spectrum
    wave = np.zeros((nspec, norders))
    wave_mask = np.zeros((nspec, norders),dtype=bool)
    flam = np.zeros((nspec, norders))
    flam_ivar = np.zeros((nspec, norders))
    flam_mask = np.zeros((nspec, norders),dtype=bool)
    for iord in range(norders):
        wave[:,iord] = sobjs[iord].optimal['WAVE']
        wave_mask[:,iord] = sobjs[iord].optimal['WAVE'] > 0.0
        flam_mask[:,iord] = wave_mask[:, iord] & sobjs[iord].optimal['MASK']
        if ret_flam:
            flam[:,iord] = sobjs[iord].optimal['FLAM']
            flam_ivar[:,iord] = sobjs[iord].optimal['FLAM_IVAR']
        else:
            flam[:,iord] = sobjs[iord].optimal['COUNTS']
            flam_ivar[:,iord] = sobjs[iord].optimal['COUNTS_IVAR']

    return wave, wave_mask, flam, flam_ivar, flam_mask, head


def sensfunc_guess(wave, counts_ps, inmask, flam_true, tell_dict_now, resln_guess, airmass_guess, polyorder, func,
                   lower=3.0, upper=3.0, debug=False):

    # Model parameter guess for starting the optimizations
    tell_pad = tell_dict_now['tell_pad']
    tell_guess = (np.median(tell_dict_now['pressure_grid']), np.median(tell_dict_now['temp_grid']),
                  np.median(tell_dict_now['h2o_grid']), airmass_guess, resln_guess, 0.0)
    tell_model1 = eval_telluric(tell_guess, tell_dict_now)
    sensguess = tell_model1[tell_pad[0]:-tell_pad[1]] * flam_true/(counts_ps + (counts_ps < 0.0))
    fitmask = inmask & np.isfinite(sensguess)
    # Perform an initial fit to the sensitivity function to set the starting point for optimization
    mask, coeff = utils.robust_polyfit_djs(wave, sensguess, polyorder, function=func,
                                       minx=wave.min(), maxx=wave.max(), inmask=fitmask, lower=lower, upper=upper,
                                       use_mad=True)
    sensfit_guess = utils.func_val(coeff, wave, func, minx=wave.min(), maxx=wave.max())

    if debug:
        plt.plot(wave, sensguess)
        plt.plot(wave, sensfit_guess)
        plt.ylim(-0.1 * sensfit_guess.min(), 1.3 * sensfit_guess.max())
        plt.show()

    return coeff

def sensfunc_telluric_joint(wave, counts_ps, counts_ps_ivar, flam_true, tell_dict, inmask=None, sensfunc=True,
                            airmass=None, resln_guess=None, pix_shift_bounds = (-2.0,2.0), resln_frac_bounds=(0.5,1.5),
                            delta_coeff_bounds=(-20.0, 20.0), minmax_coeff_bounds=(-5.0, 5.0),
                            seed=None, polyorder=7, func='legendre', maxiter=3, sticky=True, use_mad=True,
                            lower=3.0, upper=3.0, tol=1e-4, popsize=30, recombination=0.7, disp=True, polish=True,
                            debug=False):
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
    if inmask is None:
        inmask = (counts_ps_ivar > 0.0)

    if use_mad:
        invvar = None
    else:
        invvar = counts_ps_ivar

    airmass_guess = np.median(tell_dict['airmass_grid']) if airmass is None else airmass

    # This guarantees that the fit will be deterministic and hence reproducible
    if seed is None:
        seed_data = np.fmin(int(np.abs(np.sum(counts_ps[np.isfinite(counts_ps)]))), 2 ** 32 - 1)
        seed = np.random.RandomState(seed=seed_data)


    # Determine the padding and use a subset of the full tell_model_grid to make the convolutions faster
    loglam = np.log10(wave)
    dloglam = np.median(loglam[1:] - loglam[:-1])
    if resln_guess is None:
        resln_guess = 1.0/(3.0*dloglam*np.log(10.0)) # assume roughly Nyquist sampling
    pix = 1.0/resln_guess/(dloglam*np.log(10.0))/(2.0 * np.sqrt(2.0 * np.log(2))) # number of pixels per resolution element
    tell_pad = int(np.ceil(10.0 * pix))
    ind_lower = np.argmin(np.abs(tell_dict['wave_grid']-np.min(wave)))
    ind_upper = np.argmin(np.abs(tell_dict['wave_grid']-np.max(wave)))
    # This presumes that the telluric model grid is evaluated on the exact same grid as the star.
    ngrid = ind_upper-ind_lower + 1
    if ngrid != wave.size:
        msgs.error('size of tell_wave_grid={:d} != wave = {:d}. '
                   'Something is wrong with your telluric grid'.format(ngrid, wave.size))
    ind_lower_pad = np.fmax(ind_lower - tell_pad, 0)
    ind_upper_pad = np.fmin(ind_upper + tell_pad, tell_dict['wave_grid'].size-1)
    tell_pad_tuple = (ind_lower-ind_lower_pad, ind_upper_pad-ind_upper)
    tell_wave_grid = tell_dict['wave_grid'][ind_lower_pad:ind_upper_pad+1]
    tell_model_grid = tell_dict['tell_grid'][:,:,:,:,ind_lower_pad:ind_upper_pad+1]
    tell_dict_now = dict(pressure_grid=tell_dict['pressure_grid'], temp_grid=tell_dict['temp_grid'],
                         h2o_grid=tell_dict['h2o_grid'], airmass_grid=tell_dict['airmass_grid'],
                         tell_grid=tell_model_grid, tell_pad=tell_pad_tuple, dloglam=dloglam,
                         loglam = np.log10(tell_wave_grid))

    wave_min = wave.min()
    wave_max = wave.max()
    if sensfunc:
        # Joint sensitivity function and telluric fits
        chi2_func = sensfunc_tellfit_chi2
        fitting_function=sensfunc_tellfit
        # Guess the coefficients by doing a fit to the sensitivity function with the average telluric behavior
        coeff = sensfunc_guess(wave, counts_ps, inmask, flam_true, tell_dict_now, resln_guess, airmass_guess,
                               polyorder, func, lower=lower, upper=upper, debug=debug)
        # Polynomial coefficient bounds
        bounds_coeff = [(np.fmin(np.abs(this_coeff)*delta_coeff_bounds[0], minmax_coeff_bounds[0]),
                         np.fmax(np.abs(this_coeff)*delta_coeff_bounds[1], minmax_coeff_bounds[1])) for this_coeff in coeff]
    else:
        # Telluric only fits
        chi2_func = tellfit_chi2
        fitting_function=tellfit
        bounds_coeff = []

    # Set the bounds for the optimization
    bounds_tell = [(tell_dict_now['pressure_grid'].min(), tell_dict_now['pressure_grid'].max()),
                   (tell_dict_now['temp_grid'].min(), tell_dict_now['temp_grid'].max()),
                   (tell_dict_now['h2o_grid'].min(), tell_dict_now['h2o_grid'].max()),
                   (tell_dict_now['airmass_grid'].min(), tell_dict_now['airmass_grid'].max()),
                   (resln_guess*resln_frac_bounds[0], resln_guess*resln_frac_bounds[1]),
                   pix_shift_bounds]
    bounds = bounds_coeff + bounds_tell

    arg_dict = dict(wave=wave, bounds=bounds, counts_ps=counts_ps, ivar=counts_ps_ivar,
                    wave_min=wave_min, wave_max=wave_max, flam_true=flam_true, tell_dict=tell_dict_now, order=polyorder,
                    sensfunc = sensfunc, func=func, chi2_func=chi2_func)
    result, ymodel, outmask = utils.robust_optimize(counts_ps, fitting_function, arg_dict, invvar=invvar, inmask=inmask,
                                                    maxiter=maxiter,lower=lower, upper=upper, sticky=sticky,
                                                    use_mad=use_mad, bounds = bounds, tol=tol, popsize=popsize,
                                                    recombination=recombination, disp=disp, polish=polish, seed=seed)

    sens_coeff = result.x[:polyorder + 1]
    tell_params = result.x[polyorder + 1:]
    tell_pad = arg_dict['tell_dict']['tell_pad']
    telluric_fit = eval_telluric(tell_params, arg_dict['tell_dict'])[tell_pad[0]:-tell_pad[1]]
    sensfit = utils.func_val(sens_coeff, wave, arg_dict['func'], minx=arg_dict['wave_min'], maxx=arg_dict['wave_max'])
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


    return tell_params, telluric_fit, sens_coeff, sensfit



def ech_sensfunc_telluric(spec1dfile, telgridfile, star_type=None, star_mag=None, ra=None, dec=None,
                          resln_guess=None, resln_frac_bounds=(0.5,1.5), delta_coeff_bounds=(-20.0, 20.0),
                          polyorder=7, func='legendre', maxiter=3, sticky=True, use_mad=True, lower=3.0, upper=3.0, seed=None,
                          tol=1e-4, popsize=30, recombination=0.7, disp=True, polish=True,
                          debug=False):
    """
    Loop over orders to jointly fit a sensitivity function and telluric correction for a standard star spectrum for each
    order individua


    Args:
        spec1dfile:
        telgridfile:
        star_type:
        star_mag:
        ra:
        dec:
        resln_guess:
        resln_frac_range:
        delta_coeff:
        polyorder:
        func:
        maxiter:
        sticky:
        use_mad:
        lower:
        upper:
        seed:
        tol:
        popsize:
        recombination:
        disp:
        polish:
        debug:

    Returns:

    """

    # Read in standard star dictionary
    std_dict = flux.get_standard_spectrum(star_type=star_type, star_mag=star_mag, ra=ra, dec=dec)
    # This guarantees that the fit will be deterministic and hence reproducible
    if seed is None:
        seed_data = np.fmin(int(np.abs(np.sum(std_dict['flux'].value))), 2 ** 32 - 1)
        seed = np.random.RandomState(seed=seed_data)

    # Read in the telluric grid
    tell_model_dict = read_telluric_grid(telgridfile)
    wave_grid = tell_model_dict['wave_grid']

    # Read in the standard star spectrum and interpolate it onto the wave grid
    wave, wave_mask, counts, counts_ivar, counts_mask, head = unpack_orders(spec1dfile)
    exptime = head['EXPTIME']
    airmass = head['AIRMASS']

    counts_int, counts_ivar_int, counts_mask_int = coadd1d.interp_spec(wave_grid, wave, counts, counts_ivar, counts_mask)

    nspec, norders = wave.shape
    if np.size(polyorder) > 1:
        if np.size(polyorder) != norders:
            msgs.error('polyorder must have either have norder elements or be a scalar')
        polyorder_vec = np.array(polyorder)
    else:
        polyorder_vec = np.full(norders, polyorder)

    # Sort order by the strength of their telluric absorption
    srt_order_tell = telluric_sort(wave, wave_mask, tell_model_dict)

    telluric_out = np.zeros((nspec, norders))
    sensfunc_out = np.zeros((nspec, norders))
    sens_dict = {}
    tell_dict = {}
    wave_all_min=np.inf
    wave_all_max=-np.inf
    for iord in srt_order_tell:
        IPython.embed()
        msgs.info("Fitting sensitiity function for order: {:d}/{:d}".format(iord, norders))
        # Interpolate the data onto the fixed telluric grid
        wave_mask_iord = wave_mask[:,iord]
        wave_iord = wave[wave_mask_iord, iord]
        wave_all_min = np.fmin(wave_iord.min(),wave_all_min)
        wave_all_max = np.fmax(wave_iord.max(),wave_all_max)
        counts_ps = counts[wave_mask_iord, iord]/exptime
        counts_ps_ivar = counts_ivar[wave_mask_iord,iord]*exptime ** 2
        counts_ps_mask = counts_mask[wave_mask_iord,iord]
        inmask = counts_ps_mask & (counts_ps_ivar > 0.0)
        # Interpolate standard star spectrum onto the data wavelength grid
        flam_true = scipy.interpolate.interp1d(std_dict['wave'], std_dict['flux'], bounds_error=False,
                                               fill_value='extrapolate')(wave_iord)
        tell_params, tellfit, sens_coeff, sensfit = sensfunc_telluric_joint(
            wave_iord, counts_ps, counts_ps_ivar, flam_true, tell_model_dict, inmask=inmask, airmass=airmass,
            resln_guess=resln_guess, resln_frac_bounds=resln_frac_bounds, delta_coeff_bounds=delta_coeff_bounds, seed=seed,
            polyorder=polyorder_vec[iord], func=func, maxiter=maxiter, sticky=sticky, use_mad=use_mad, lower=lower, upper=upper,
            tol=tol, popsize=popsize, recombination=recombination, disp=disp, polish=polish, debug=debug)
        telluric_out[wave_mask_iord,iord] = tellfit
        sensfunc_out[wave_mask_iord,iord] = sensfit
        sens_dict[str(iord)] = dict(polyorder=polyorder_vec[iord], wave_min=wave_iord.min(), wave_max=wave_iord.max(),
                                    sens_coeff=sens_coeff,
                                    wave=wave[:,iord], wave_mask=wave_mask[:,iord],
                                    sensfunc=sensfunc_out[:,iord])
        tell_dict[str(iord)] = dict(wave_min=wave_iord.min(), wave_max=wave_iord.max(),
                                    pressure=tell_params[0], temp=tell_params[1], h2o=tell_params[2],
                                    airmass=tell_params[3], resolution=tell_params[4], shift=0.0,
                                    wave=wave[:,iord], wave_mask=wave_mask[:,iord],
                                    telluric=telluric_out[:,iord])

    sens_dict['meta'] = dict(exptime=exptime, airmass=airmass, nslits=norders, std_file=spec1dfile,
                             cal_file=std_dict['cal_file'], std_name=std_dict['name'],
                             std_ra=std_dict['std_ra'], std_dec=std_dict['std_dec'],
                             wave_min=wave_all_min, wave_max=wave_all_max)
    tell_dict['meta'] = dict(airmass=airmass, nslits=norders, wave_min=wave_all_min, wave_max=wave_all_max)

    return sens_dict, tell_dict

def telluric_sort(wave, wave_mask, tell_dict):

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

def populate_orders(quantity, ind):

    norders = ind.shape[0]
    nspec_ord = ind[:,1] - ind[:, 0] + 1
    if len(set(nspec_ord.tolist()))!= 1:
        msgs.error('There is a problem with your upper and lower indices. They do not all imply spectrum of the same length')
    nspec = nspec_ord[0]
    quantity_orders = np.zeros((nspec, norders), dtype=quantity.dtype)
    for iord in range(norders):
        quantity_orders[:,iord] = quantity[ind[iord, 0]:ind[iord, 1]+1]

    return quantity_orders


def get_inmask_orders(wave, wavegrid_inmask, inmask):

    norders = wave.shape[1]
    inmask_orders = np.ones_like(wave, dtype=bool)
    if inmask is not None:
        for iord in range(norders):
            ind_lower = np.argmin(np.abs(wavegrid_inmask - wave[:,iord].min()))
            ind_upper = np.argmin(np.abs(wavegrid_inmask - wave[:,iord].max()))
            inmask_orders[:, iord] = inmask[ind_lower:ind_upper+1]

    return inmask_orders

def ech_telluric_qso(spec1dfile, telgridfile, pcafile, npca, z_qso, inmask=None, wavegrid_inmask=None,
                     delta_zqso = 0.1, bounds_norm = (0.7,1.3), resln_guess=None, resln_frac_range=(0.5,1.5),
                     tell_norm_thresh=0.9, maxiter=3, sticky=True, use_mad=True, lower=3.0, upper=3.0,
                     seed=None, tol=1e-3, popsize=30, recombination=0.65, disp=True, polish=True,
                     debug=True):

    # Read in the data
    wave, wave_mask, flam, flam_ivar, flam_mask, head = unpack_orders(spec1dfile, ret_flam=True)
    airmass = head['AIRMASS']
    nspec, norders = wave.shape

    # Populate the input mask along the orders
    inmask_orders = get_inmask_orders(wave, wavegrid_inmask, inmask)

    # This guarantees that the fit will be deterministic and hence reproducible
    if seed is None:
        seed_data = np.fmin(int(np.abs(np.sum(flam))), 2 ** 32 - 1)
        seed = np.random.RandomState(seed=seed_data)

    # Are we using median absolute deviation for the fits?
    if use_mad:
        invvar = None
    else:
        invvar = flam_ivar

    # Guess resolution from spectral sampling
    loglam = np.log10(wave[:, 0])
    dloglam = np.median(loglam[1:] - loglam[:-1])
    if resln_guess is None:
        resln_guess = 1.0/(3.0*dloglam*np.log(10))  # assume roughly Nyquist sampling

    # Determine the padding and use a subset of the full tell_model_grid to make the convolutions faster
    pix_per_sigma = 1.0 / resln_guess /(dloglam*np.log(10))/ (2.0 * np.sqrt(2.0 * np.log(2)))  # number of pixels per resolution element
    tell_pad = int(np.ceil(15.0*pix_per_sigma))

    # Read in the telluric grid
    tell_dict = read_telluric_grid(telgridfile,wave_min = wave.min(), wave_max=wave.max(), pad=tell_pad)
    # Determine the indices to populate telluric models onto our orders
    ind_tell = np.zeros((norders, 2), dtype=int)
    for iord in range(norders):
        ind_tell[iord, 0] = np.argmin(np.abs(tell_dict['wave_grid'] - np.min(wave[:,iord])))
        ind_tell[iord, 1] = np.argmin(np.abs(tell_dict['wave_grid'] - np.max(wave[:,iord])))
    tell_dict['dloglam'] = dloglam
    tell_dict['tell_pad'] = (0,0) # no padding since we already padded everything in determining the wavelength range above

    # Read in the PCA model information
    pca_dict = init_pca(pcafile,tell_dict['wave_grid'],z_qso, npca)
    pca_dict['ind'] = ind_tell # The PCA is on the same grid as the telluric models and thus has the same ind_lower/ind_upper

    # Estimate the normalization for setting the bounds
    tell_guess = (np.median(tell_dict['pressure_grid']), np.median(tell_dict['temp_grid']),
                  np.median(tell_dict['h2o_grid']), airmass, resln_guess)
    tell_model = eval_telluric(tell_guess, tell_dict)
    tell_model_orders = populate_orders(tell_model, ind_tell)
    # Just use the mean PCA model for estimating the normalization
    pca_mean_orders = populate_orders(np.exp(pca_dict['components'][0,:]), ind_tell)
    tell_mask = tell_model_orders > tell_norm_thresh
    inmask_tot = inmask_orders & flam_mask
    data_norm = np.sum(tell_mask*inmask_tot*flam) # Data has absorption in it so we don't multilply by tell_model
    # Model has no absorption in it, so we use the average telluric absorption in the grid
    pca_norm = np.sum(tell_mask*inmask_tot*tell_model_orders*pca_mean_orders)
    flam_norm = data_norm/pca_norm

    # Set the bounds for the PCA and truncate to the right dimension
    coeffs = pca_dict['coeffs'][:,1:npca]
    # Compute the min and max arrays of the coefficients which are not the norm, i.e. grab the coeffs that aren't the first one
    coeff_min = np.amin(coeffs, axis=0)  # only
    coeff_max = np.amax(coeffs, axis=0)
    bounds_z = [(z_qso - delta_zqso, z_qso + delta_zqso)]                # QSO redshift: can vary within delta_zqso
    bounds_flam = [(flam_norm*bounds_norm[0], flam_norm*bounds_norm[1])] # Norm: bounds determined from estimate above
    bounds_coeff = [(i, j) for i, j in zip(coeff_min, coeff_max)]        # Coefficients:  determined from PCA model
    # Set the bounds for the telluric
    bounds_tell = [(tell_dict['pressure_grid'].min(), tell_dict['pressure_grid'].max()),
                   (tell_dict['temp_grid'].min()    , tell_dict['temp_grid'].max()),
                   (tell_dict['h2o_grid'].min()     , tell_dict['h2o_grid'].max()),
                   (tell_dict['airmass_grid'].min() , tell_dict['airmass_grid'].max()),
                   (resln_guess*resln_frac_range[0] , resln_guess*resln_frac_range[1])]
    # Final bounds for the optimizaiton
    bounds =  bounds_z + bounds_flam + bounds_coeff + bounds_tell
    # Create the arg_dict
    arg_dict = dict(nspec=nspec, norders=norders, npca=npca, ind=ind_tell, flam_ivar=flam_ivar, tell_dict=tell_dict,
                    pca_dict=pca_dict, tell_qso_chi2 = tell_qso_chi2)
    result, ymodel, outmask = utils.robust_optimize(flam, tell_qso_fit, arg_dict, invvar=invvar, inmask=inmask_tot,
                                                    maxiter=maxiter,lower=lower, upper=upper, sticky=sticky,
                                                    use_mad=use_mad, bounds = bounds, tol=tol, popsize=popsize,
                                                    recombination=recombination, disp=disp, polish=polish, seed=seed)

    theta = result.x
    flam_model = tell_qso_evalmodel(theta, arg_dict)
    pca_params = theta[:npca + 1]
    tell_params = theta[-5:]
    tell_model = eval_telluric(tell_params, arg_dict['tell_dict'])
    tell_model_orders = populate_orders(tell_model, ind_tell)
    pca_model = pca_eval(pca_params, arg_dict['pca_dict'])
    pca_model_orders = populate_orders(pca_model, ind_tell)

    if debug:
        color_scheme = [('black', 'red'), ('cornflowerblue', 'magenta')]
        # Plot the data
        for iord in range(norders):
            colors = color_scheme[np.mod(iord, 2)]
            this_mask = wave_mask[:,iord]
            plt.plot(wave[this_mask,iord], scipy.ndimage.filters.median_filter(flam[this_mask,iord], size=15), color=colors[0], drawstyle='steps-mid')
        plt.ylim((-0.2, flam_model.max()))
        plt.legend()
        plt.show()


        # Plot what we actually fit
        for iord in range(norders):
            colors = color_scheme[np.mod(iord, 2)]
            this_mask = wave_mask[:,iord]
            plt.plot(wave[this_mask,iord], flam[this_mask,iord], color=colors[0], drawstyle='steps-mid')
            plt.plot(wave[this_mask,iord], flam_model[this_mask,iord], color=colors[1], drawstyle='steps-mid')
        plt.ylim((0.0, flam_model.max()))
        plt.legend()
        plt.show()

        # Plot the telluric corrected and rescaled orders
        for iord in range(norders):
            colors = color_scheme[np.mod(iord, 2)]
            this_mask = wave_mask[:,iord]
            spec_iord = flam[this_mask,iord]/(tell_model_orders[this_mask, iord] + (tell_model_orders[this_mask,iord] == 0.0))
            plt.plot(wave[this_mask,iord],scipy.ndimage.filters.median_filter(spec_iord, size=15),
                     color=colors[0], drawstyle='steps-mid')
        plt.plot(tell_dict['wave_grid'], pca_model, color='green', label='pca model')
        plt.ylim((0.0, pca_model.max()*1.5))
        plt.legend()
        plt.show()

    qso_pca_dict = dict(wave=wave, wave_mask=wave_mask, flam=flam, flam_mask=flam_mask, flam_ivar=flam_ivar,
                    airmass=airmass,  ind=ind_tell,
                    tell_params=tell_params, tell_model=tell_model, tell_model_orders=tell_model_orders,
                    pca_params=pca_params, pca_model=pca_model, pca_model_orders=pca_model_orders, result = result)

    return qso_pca_dict

def ech_telluric(wave, wave_mask, flam, flam_ivar, flam_mask, flam_true, airmass, telgridfile,
                 inmask=None, wavegrid_inmask=None,
                 resln_guess=None, resln_frac_bounds=(0.5,1.5), maxiter=3, sticky=True, use_mad=True, lower=3.0, upper=3.0,
                 seed=None, tol=1e-4, popsize=30, recombination=0.7, disp=True, polish=True,
                 debug=False):
    """
    Loop over orders to fit a telluric correction order by order given a known spectrum flam_true


    Args:
        spec1dfile:
        telgridfile:
        star_type:
        star_mag:
        ra:
        dec:
        resln_guess:
        resln_frac_range:
        delta_coeff:
        polyorder:
        func:
        maxiter:
        sticky:
        use_mad:
        lower:
        upper:
        seed:
        tol:
        popsize:
        recombination:
        disp:
        polish:
        debug:

    Returns:

    """

    # This guarantees that the fit will be deterministic and hence reproducible
    if seed is None:
        seed_data = np.fmin(int(np.abs(np.sum(flam))), 2 ** 32 - 1)
        seed = np.random.RandomState(seed=seed_data)


    nspec, norders = wave.shape

    # Populate the input mask along the orders
    inmask_orders = get_inmask_orders(wave, wavegrid_inmask, inmask)

    # Read in the telluric grid
    tell_model_dict = read_telluric_grid(telgridfile)
    # Sort orders by the strength of their telluric absorption
    srt_order_tell = telluric_sort(wave, wave_mask, tell_model_dict)
    telluric_out = np.zeros((nspec, norders))
    tell_dict = {}
    wave_all_min=np.inf
    wave_all_max=-np.inf
    for iord in srt_order_tell:
        msgs.info("Fitting telluric absorption for order: {:d}/{:d}".format(iord, norders))
        wave_mask_iord = wave_mask[:,iord]
        wave_iord = wave[wave_mask_iord, iord]
        wave_all_min = np.fmin(wave_iord.min(),wave_all_min)
        wave_all_max = np.fmax(wave_iord.max(),wave_all_max)
        flam_iord = flam[wave_mask_iord, iord]
        flam_ivar_iord = flam_ivar[wave_mask_iord,iord]
        flam_mask_iord = flam_mask[wave_mask_iord,iord]
        inmask_iord = inmask_orders[wave_mask_iord, iord]
        inmask_tot = flam_mask_iord & (flam_ivar_iord > 0.0) & inmask_iord
        flam_true_iord = flam_true[wave_mask_iord, iord]
        tell_params, tellfit, sens_coeff, sensfit = sensfunc_telluric_joint(
            wave_iord, flam_iord, flam_ivar_iord, flam_true_iord, tell_model_dict, inmask=inmask_tot, sensfunc=False,
            airmass=airmass, resln_guess=resln_guess, resln_frac_bounds=resln_frac_bounds, seed=seed, maxiter=maxiter,
            sticky=sticky, use_mad=use_mad, lower=lower, upper=upper,
            tol=tol, popsize=popsize, recombination=recombination, disp=disp, polish=polish, debug=debug)
        telluric_out[wave_mask_iord,iord] = tellfit
        tell_dict[str(iord)] = dict(wave_min=wave_iord.min(), wave_max=wave_iord.max(),
                                    pressure=tell_params[0], temp=tell_params[1], h2o=tell_params[2],
                                    airmass=tell_params[3], resolution=tell_params[4],
                                    wave=wave[:,iord], wave_mask=wave_mask[:,iord],
                                    telluric=telluric_out[:,iord])

    tell_dict['meta'] = dict(airmass=airmass, nslits=norders, wave_min=wave_all_min, wave_max=wave_all_max)

    return tell_dict


dev_path = os.getenv('PYPEIT_DEV')

# NIRES
#spec1dfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/NIRES/J0252/ut180903/Science_coadd/spec1d_HIP13917.fits')
#telgridfile = os.path.join(dev_path, 'dev_algorithms/sensfunc/TelFit_MK_NIR_9300_26100_AM1.000_R8000.fits')
#polyorder = [7,11,7,7,7]
#star_type='A0'
#star_mag = 8.63

# XSHOOTER
#star_mag  = None
#star_type = None
#spec1dfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/Pypeit_files/PISCO_nir_REDUCED/Science_coadd/spec1d_STD.fits')
#header = fits.getheader(spec1dfile)
#telgridfile =  os.path.join(dev_path, 'dev_algorithms/sensfunc/TelFit_Paranal_NIR_9800_25000_R25000.fits')
#polyorder=6
#sens_dict = ech_sensfunc_telluric(spec1dfile, telgridfile, polyorder=polyorder, ra=header['RA'], dec=header['DEC'],
#                          star_mag=star_mag, star_type=star_type)
#sensfuncfile = 'EG274_sensfunc.json'
#save.save_sens_dict(sens_dict,sensfuncfile)
#sens_dict = load.load_sens_dict(sensfuncfile)


# First fit the sensivity function using the standard star
star_mag  = None
star_type = None
spec1dfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Science/spec1d_XSHOO.2018-11-08T00:11:57.074-Feige110_XShooter_NIR_2018Nov08T001157.074.fits')
#spec1dfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439/vlt_xshooter_nir/Science/spec1d_Feige110.fits')
header = fits.getheader(spec1dfile)
telgridfile =  os.path.join(dev_path, 'dev_algorithms/sensfunc/TelFit_Paranal_NIR_9800_25000_R25000.fits')
polyorder=6
sens_dict, tell_dict = ech_sensfunc_telluric(spec1dfile, telgridfile, polyorder=polyorder, ra=header['RA'], dec=header['DEC'],
                                             star_mag=star_mag, star_type=star_type,tol=1e-4,popsize=100, debug=True)
# Write the sens_dict and tell_dict out to a file
sensfuncfile = 'Feige110_sensfunc.json'
save.save_sens_dict(sens_dict,sensfuncfile)
# Test loading
sens_dict1 = load.load_sens_dict(sensfuncfile)

telluricfile = 'Feige110_telluric.json'
save.save_sens_dict(tell_dict, telluricfile)
# Test loading
tell_dict1 = load.load_sens_dict(telluricfile)

# Now flux calibrate the data. At the moment this is a bit cumbersome becuase the interface to the fluxspec class needs to be improved. Probably
# the easiest way is to use the script for now. This could be a stack of 1d files, or a 1d file from a 2d co-add

## Now we co-add the fluxed data if it is not already co-added and get a high S/N order by order spectrum

spec1dfluxfile = '/Users/joe/Dropbox/PypeIt_Redux/XSHOOTER/Pypeit_files/PISCO_nir_REDUCED/Science_coadd_feb19/spec1d_flux_PSOJ205p09.fits'


# Now run the piece of code to rescale the individual orders to match each other in the overlap regions. Now we have
# high S/N fluxed, and matched spectra for each order.


# Now perform a global QSO PCA fit to all the orders combined with order by order telluric fits to each order.
spec1dfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439/vlt_xshooter_nir/Science_coadd/spec1d_J0439_flux_feige110.fits')
pcafile = os.path.join(dev_path, 'dev_algorithms/sensfunc/qso_pca_1200_3100.pckl')
npca = 8
z_qso = 6.51#7.54#6.51
vlt_xshooter_nir = util.load_spectrograph('vlt_xshooter_nir')
wavegrid = vlt_xshooter_nir.wavegrid(midpoint=True)
# inmask here includes both a Lya mask (unnecessary for lens), BAL mask, and cuts off at the end of the PCA wave grid
inmask = (wavegrid > (1.0 + z_qso)*1220.0) & ((wavegrid < 10770) | (wavegrid > 11148)) & (wavegrid < 3099*(1+z_qso))
qso_pca_dict = ech_telluric_qso(spec1dfile, telgridfile, pcafile, npca, z_qso, inmask=inmask, wavegrid_inmask=wavegrid, debug=True)

# Now that we have the PCA model perform an order by order fit to the telluric absorption as we did for the standard. This allows
# resolution and a small wavelength shift to be free parameters

# TODO we need to add an absorption threshold in median telluric absorption after which the code uses the average of the model
# coming from all the other orders
tell_qso_dict = ech_telluric(qso_pca_dict['wave'], qso_pca_dict['wave_mask'], qso_pca_dict['flam'], qso_pca_dict['flam_ivar'],
                             qso_pca_dict['flam_mask'], qso_pca_dict['pca_model_orders'], qso_pca_dict['airmass'],telgridfile,
                             tol=1e-4, popsize=30, recombination=0.7, disp=True, polish=True, debug=False)

# Apply this telluric to the data, now merge the data




#
#
# def sensfunc_old(theta, arg_dict):
#
#     wave_star = arg_dict['wave_star']
#     counts_ps = arg_dict['counts_ps']
#     counts_ps_ivar = arg_dict['counts_ps_ivar']
#     thismask = arg_dict['thismask']
#     wave_min = arg_dict['wave_min']
#     wave_max = arg_dict['wave_max']
#     flam_true = arg_dict['flam_true']
#     tell_dict = arg_dict['tell_dict']
#     order = arg_dict['order']
#     func = arg_dict['func']
#
#     theta_sens = theta[:order+1]
#     theta_tell = theta[order+1:]
#     sensmodel = utils.func_val(theta_sens, wave_star, func, minx=wave_min, maxx=wave_max)
#     tellmodel_conv = eval_telluric(theta_tell, wave_star, tell_dict)
#     if np.sum(sensmodel) < 1e-6:
#         return np.inf
#     else:
#         chi_vec = thismask*(sensmodel != 0.0)*(tellmodel_conv*flam_true/(sensmodel + (sensmodel == 0.0)) -
#                                                counts_ps)*np.sqrt(counts_ps_ivar)
#         chi2 = np.sum(np.square(chi_vec))
#     return chi2
#
#
# def sens_tellfit_old(optfunc, bounds, arg_dict, tol=1e-4, popsize=30, recombination=0.7, disp=True, polish=True, seed=None):
#
#     result = scipy.optimize.differential_evolution(optfunc, args=(arg_dict,), tol=tol,
#                                                    bounds=bounds, popsize=popsize,recombination=recombination,
#                                                    disp=disp, polish=polish, seed=seed)
#     wave_star = arg_dict['wave_star']
#     order = arg_dict['order']
#     coeff_out = result.x[:order+1]
#     tell_out = result.x[order+1:]
#     tellfit_conv = eval_telluric(tell_out, wave_star,arg_dict['tell_dict'])
#     sensfit = utils.func_val(coeff_out, wave_star, arg_dict['func'], minx=arg_dict['wave_min'], maxx=arg_dict['wave_max'])
#
#     return result, tellfit_conv, sensfit, coeff_out, tell_out

#std_dict = {}
#std_dict['std_ra'] = head['RA']
#std_dict['std_dec'] = head['DEC']
#std_dict['exptime'] = exptime
#std_dict['airmass'] = head['AIRMASS']
#std_dict['std_name'] = ['EG274']
#std_dict['cal_file'] = ['EG274']
#std_dict['wave'] = wave_std
#std_dict['flux'] = flux_std


#dev_path = os.getenv('PYPEIT_DEV')
#xshooter_file = os.path.join(dev_path,'dev_algorithms/sensfunc/xshooter_standards/fEG274.dat')
#output = np.loadtxt(xshooter_file)
#wave_std = output[:,0]
#flux_std = output[:,1]/PYPEIT_FLUX_SCALE

        # Create a fake std_dict for EG274
#std_dict = flux.get_standard_spectrum(star_type=star_type, star_mag=star_mag, ra=ra, dec=dec)

#
#
#
# def ech_telluric_qso_old(spec1dfile, telgridfile, pcafile, npca, z_qso, inmask=None, wavegrid_inmask=None, delta_zqso = 0.1,
#                             bounds_norm = (0.7,1.3), bounds_rescale = (0.95,1.05), resln_guess=None, resln_frac_range=(0.7,1.3),
#                             tell_norm_thresh=0.9, maxiter=3, sticky=True,
#                             use_mad=True, lower=3.0, upper=3.0, seed=None, tol=1e-2, popsize=30, recombination=0.65, disp=True, polish=True,
#                             debug=True):
#
#
#     # Read in the spec1d file
#     sobjs, head = load.load_specobjs(spec1dfile)
#     airmass = head['AIRMASS']
#     norders = len(sobjs)
#     nspec = sobjs[0].optimal['COUNTS'].size
#     # Allocate arrays and unpack spectrum
#     wave = np.zeros((nspec, norders))
#     wave_mask = np.zeros((nspec, norders),dtype=bool)
#     flux = np.zeros((nspec, norders))
#     flux_ivar = np.zeros((nspec, norders))
#     flux_mask = np.zeros((nspec, norders),dtype=bool)
#     inmask_orders = np.ones((nspec, norders), dtype=bool)
#     for iord in range(norders):
#         flux[:,iord] = sobjs[iord].optimal['FLAM']
#         flux_ivar[:,iord] = sobjs[iord].optimal['FLAM_IVAR']
#         wave[:,iord] = sobjs[iord].optimal['WAVE_GRID']
#         wave_mask[:,iord] = sobjs[iord].optimal['WAVE_GRID_MASK']
#         flux_mask[:,iord] = wave_mask[:, iord] & sobjs[iord].optimal['MASK']
#         # If the user input a mask, populate it onto the orders using wavegrid_inmask
#         if inmask is not None:
#             ind_lower = np.argmin(np.abs(wavegrid_inmask - wave[:,iord].min()))
#             ind_upper = np.argmin(np.abs(wavegrid_inmask - wave[:,iord].max()))
#             inmask_orders[:, iord] = inmask[ind_lower:ind_upper+1]
#
#     if use_mad:
#         invvar = None
#     else:
#         invvar = flux_ivar
#
#     # This guarantees that the fit will be deterministic and hence reproducible
#     if seed is None:
#         seed_data = np.fmin(int(np.abs(np.sum(flux))), 2 ** 32 - 1)
#         seed = np.random.RandomState(seed=seed_data)
#
#     loglam = np.log10(wave[:, 0])
#     dloglam = np.median(loglam[1:] - loglam[:-1])
#     # Guess resolution from spectral sampling
#     if resln_guess is None:
#         resln_guess = 1.0/(3.0*dloglam*np.log(10))  # assume roughly Nyquist sampling
#
#     # Determine the padding and use a subset of the full tell_model_grid to make the convolutions faster
#     pix_per_sigma = 1.0 / resln_guess /(dloglam*np.log(10))/ (2.0 * np.sqrt(2.0 * np.log(2)))  # number of pixels per resolution element
#     tell_pad = int(np.ceil(15.0*pix_per_sigma))
#
#     # Read in the telluric grid
#     tell_dict = read_telluric_grid(telgridfile,wave_min = wave.min(), wave_max=wave.max(), pad=tell_pad)
#     # Determine the indices for populate telluric models onto our orders
#     ind_tell = np.zeros((norders, 2), dtype=int)
#     for iord in range(norders):
#         ind_tell[iord, 0] = np.argmin(np.abs(tell_dict['wave_grid'] - np.min(wave[:,iord])))
#         ind_tell[iord, 1] = np.argmin(np.abs(tell_dict['wave_grid'] - np.max(wave[:,iord])))
#
#     tell_dict['dloglam'] = dloglam
#     tell_dict['tell_pad'] = (0,0) # no padding since we already padded everything in determining the wavelength range above
#
#     # Read in the PCA
#     pca_dict = init_pca(pcafile,tell_dict['wave_grid'],z_qso, npca)
#     pca_dict['ind'] = ind_tell # The PCA is on the same grid as the telluric models and thus has the same ind_lower/ind_upper
#
#     # Determine the normalization
#     tell_guess = (np.median(tell_dict['pressure_grid']), np.median(tell_dict['temp_grid']),
#                   np.median(tell_dict['h2o_grid']), airmass, resln_guess)
#     tell_model = eval_telluric(tell_guess, tell_dict)
#     tell_model_orders = populate_orders(tell_model, ind_tell)
#     # Just use the mean for estimating the normalization
#     pca_mean_orders = populate_orders(np.exp(pca_dict['components'][0,:]), ind_tell)
#     tell_mask = tell_model_orders > tell_norm_thresh
#     inmask_tot = inmask_orders & flux_mask
#     data_norm = np.sum(tell_mask*inmask_tot*flux) # Data has absorption in it so we don't multilply by tell_model
#     # Model has no absorption in it, so we model by average model absorption
#     pca_norm = np.sum(tell_mask*inmask_tot*tell_model_orders*pca_mean_orders)
#     flux_norm = data_norm/pca_norm
#
#     # Set the bounds for the PCA and truncate to the right dimension
#     coeffs = pca_dict['coeffs'][:,1:npca]
#     # Compute the min and max arrays of the coefficients which are not the norm, i.e. grab the coeffs that aren't the first one
#     coeff_min = np.amin(coeffs, axis=0)  # only
#     coeff_max = np.amax(coeffs, axis=0)
#     # Determine the norm bounds from the estimate above
#     bounds_scale = [bounds_rescale]*norders
#     bounds_z = [(z_qso - delta_zqso, z_qso + delta_zqso)]
#     bounds_flux = [(flux_norm*bounds_norm[0], flux_norm*bounds_norm[1])]
#     bounds_coeff = [(i, j) for i, j in zip(coeff_min, coeff_max)]
#     # Set the bounds for the optimization
#     bounds_tell = [(tell_dict['pressure_grid'].min(), tell_dict['pressure_grid'].max()),
#                    (tell_dict['temp_grid'].min()    , tell_dict['temp_grid'].max()),
#                    (tell_dict['h2o_grid'].min()     , tell_dict['h2o_grid'].max()),
#                    (tell_dict['airmass_grid'].min() , tell_dict['airmass_grid'].max()),
#                    (resln_guess*resln_frac_range[0] , resln_guess*resln_frac_range[1])]
#     bounds = bounds_scale + bounds_z + bounds_flux + bounds_coeff + bounds_tell
#
#     # Create the arg_dict
#     arg_dict = dict(nspec=nspec, norders=norders, npca=npca, ind=ind_tell, flux_ivar=flux_ivar, tell_dict=tell_dict, pca_dict=pca_dict,
#                     tell_qso_chi2 = tell_qso_chi2)
#
#     result, ymodel, outmask = utils.robust_optimize(flux, tell_qso_fit, arg_dict, invvar=invvar, inmask=inmask_tot,
#                                                     maxiter=maxiter,lower=lower, upper=upper, sticky=sticky,
#                                                     use_mad=use_mad, bounds = bounds, tol=tol, popsize=popsize,
#                                                     recombination=recombination, disp=disp, polish=polish, seed=seed)
#
#     theta = result.x
#     flux_model = tell_qso_evalmodel(theta, arg_dict)
#     theta_rescale = theta[0:norders]
#     theta_PCA = theta[norders:norders + npca + 1]
#     theta_tell = theta[-5:]
#     tell_model = eval_telluric(theta_tell, arg_dict['tell_dict'])
#     tell_model_orders = populate_orders(tell_model, ind_tell)
#     pca_model = pca_eval(theta_PCA, arg_dict['pca_dict'])
#     pca_model_orders = populate_orders(pca_model, ind_tell)
#     rescale_orders = np.outer(np.ones(nspec), theta_rescale)
#
#     color_scheme = [('black', 'red'), ('cornflowerblue', 'magenta')]
#     if debug:
#         # Plot the data
#         for iord in range(norders):
#             colors = color_scheme[np.mod(iord, 2)]
#             this_mask = wave_mask[:,iord]
#             plt.plot(wave[this_mask,iord], scipy.ndimage.filters.median_filter(flux[this_mask,iord], size=15), color=colors[0], drawstyle='steps-mid')
#         plt.ylim((-0.2, flux_model.max()))
#         plt.legend()
#         plt.show()
#
#
#         # Plot what we actually fit
#         for iord in range(norders):
#             colors = color_scheme[np.mod(iord, 2)]
#             this_mask = wave_mask[:,iord]
#             plt.plot(wave[this_mask,iord], flux[this_mask,iord], color=colors[0], drawstyle='steps-mid')
#             plt.plot(wave[this_mask,iord], flux_model[this_mask,iord], color=colors[1], drawstyle='steps-mid')
#         plt.ylim((0.0, flux_model.max()))
#         plt.legend()
#         plt.show()
#
#         # Plot the telluric corrected and rescaled orders
#         for iord in range(norders):
#             colors = color_scheme[np.mod(iord, 2)]
#             this_mask = wave_mask[:,iord]
#             spec_iord = flux[this_mask,iord]*rescale_orders[this_mask, iord]/\
#                         (tell_model_orders[this_mask, iord] + (tell_model_orders[this_mask,iord] == 0.0))
#             plt.plot(wave[this_mask,iord],scipy.ndimage.filters.median_filter(spec_iord, size=15),
#                      color=colors[0], drawstyle='steps-mid')
#         plt.plot(tell_dict['wave_grid'], pca_model, color='green', label='pca model')
#         plt.ylim((0.0, pca_model.max()*1.5))
#         plt.legend()
#         plt.show()
#
#     from IPython import embed
#     embed()
#     return result


