

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


def sensfunc_guess(wave, counts_ps, inmask, flam_true, tell_dict_now, resln_guess, airmass_guess, polyorder, func,
                   ind_lower, ind_upper,
                   lower=3.0, upper=3.0, debug=False):

    # Model parameter guess for starting the optimizations
    tell_guess = (np.median(tell_dict_now['pressure_grid']), np.median(tell_dict_now['temp_grid']),
                  np.median(tell_dict_now['h2o_grid']), airmass_guess, resln_guess, 0.0)
    tell_model1 = eval_telluric(tell_guess, tell_dict_now, ind_lower, ind_upper)
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

#

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


def qso_guess(flux, ivar, mask, pca_dict, tell_dict, airmass, resln_guess, tell_norm_thresh, ind_lower, ind_upper):

    # Estimate the normalization for setting the bounds
    tell_guess = (np.median(tell_dict['pressure_grid']), np.median(tell_dict['temp_grid']),
                  np.median(tell_dict['h2o_grid']), airmass, resln_guess, 0.0)
    tell_model = eval_telluric(tell_guess, tell_dict, ind_lower, ind_upper)
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


def update_bounds(bounds, delta_coeff, coeff):

    bounds_new = [(this_coeff * delta_coeff[0], this_coeff * delta_coeff[1]) for this_coeff in coeff]
    bounds_tell = bounds[len(coeff):]
    bounds_new.extend(bounds_tell)
    return bounds_new


## NOT YET USED
def tellfit_chi2(theta, flux, thismask, arg_dict):

    flux_ivar = arg_dict['ivar']
    flux_true = arg_dict['flux_true']
    tell_dict = arg_dict['tell_dict']
    tell_pad = tell_dict['tell_pad']
    tellmodel_conv = eval_telluric(theta, tell_dict)
    chi_vec = thismask*(tellmodel_conv*flux_true - flux)*np.sqrt(flux_ivar)
    chi2 = np.sum(np.square(chi_vec))

    return chi2

def tellfit(flam, thismask, arg_dict, **kwargs_opt):

    # Function that we are optimizing
    chi2_func= arg_dict['chi2_func']
    result = scipy.optimize.differential_evolution(chi2_func, args=(flam, thismask, arg_dict,), **kwargs_opt)
    tell_out = result.x
    tellfit_conv = eval_telluric(tell_out, arg_dict['tell_dict'])
    tell_pad = arg_dict['tell_dict']['tell_pad']
    flam_model = tellfit_conv*arg_dict['flam_true']

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
                          polyorder=7, func='legendre', maxiter=3, sticky=True, use_mad=True, lower=3.0, upper=3.0,
                          debug = False,
                          seed=None, tol=1e-4, popsize=30, recombination=0.7, disp=True, polish=True):
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

# Deprecated
#def trim_tell_dict(tell_dict, ind_lower, ind_upper):
#
#    wave_grid = tell_dict['wave_grid']
#    ind_lower_pad = np.fmax(ind_lower - tell_dict['tell_pad_pix'], 0)
#    ind_upper_pad = np.fmin(ind_upper + tell_dict['tell_pad_pix'], wave_grid.size - 1)
#    tell_pad_tuple = (ind_lower - ind_lower_pad, ind_upper_pad - ind_upper)
#    tell_wave_grid = wave_grid[ind_lower_pad:ind_upper_pad + 1]
#    tell_model_grid = tell_dict['tell_grid'][:, :, :, :, ind_lower_pad:ind_upper_pad + 1]
#    tell_dict_fit = dict(pressure_grid=tell_dict['pressure_grid'], temp_grid=tell_dict['temp_grid'],
#                         h2o_grid=tell_dict['h2o_grid'], airmass_grid=tell_dict['airmass_grid'],
#                         tell_grid=tell_model_grid, tell_pad_tuple=tell_pad_tuple,
#                         resln_guess=tell_dict['resln_guess'],pix_per_R=tell_dict['pix_per_R'],
#                         dloglam=tell_dict['dloglam'], loglam=np.log10(tell_wave_grid))
#    return tell_dict_fit



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


# deprecated
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
        chi_vec = thismask * (sensmodel != 0.0) * (counts_model - counts_ps) * np.sqrt(counts_ps_ivar)
        # Robustly characterize the dispersion of this distribution
        # chi_mean, chi_median, chi_std = stats.sigma_clipped_stats(chi_vec, np.invert(thismask), cenfunc='median',
        #                                                          stdfunc=stats.mad_std, maxiters=5, sigma=2.0)
        robust_scale = 2.0
        huber_vec = scipy.special.huber(robust_scale, chi_vec)
        loss_function = np.sum(np.square(huber_vec * thismask))
        return loss_function


# deprecated
def qso_tellfit_chi2(theta, flux, thismask, arg_dict):
    tell_model, pca_model, ln_pca_pri = qso_tellfit_eval(theta, arg_dict)
    chi_vec = thismask * (flux - tell_model * pca_model) * np.sqrt(arg_dict['flux_ivar'])

    robust_scale = 2.0
    huber_vec = scipy.special.huber(robust_scale, chi_vec)
    chi2_remap = np.sum(np.square(huber_vec * thismask))

    lnL = -chi2_remap / 2.0
    lnptot = lnL + ln_pca_pri
    chi2_tot = -2.0 * lnptot
    return chi2_tot


## Not used
def get_dloglam_data(wave):

    loglam = np.log10(wave)
    loglam_ma = np.ma.array(np.copy(loglam))
    loglam_ma.mask = np.invert(wave > 10.0)
    dloglam_arr = np.abs(loglam_ma - np.roll(loglam_ma, 1, axis=0))[1:, :]
    dloglam = np.ma.median(dloglam_arr)

    return dloglam


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



def sensfunc_tellfit(counts_ps, thismask, arg_dict, **kwargs_opt):

    # Function that we are optimizing
    chi2_func = arg_dict['chi2_func']
    counts_ps_ivar = arg_dict['ivar']
    bounds = arg_dict['bounds']
    seed = arg_dict['seed']
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



def qso_tellfit(flux, thismask, arg_dict, **kwargs_opt):

    # Function that we are optimizing
    chi2_func = arg_dict['chi2_func']
    flux_ivar = arg_dict['ivar']
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


    pix = 1.0/resln_guess/(dloglam*np.log(10.0))/(2.0 * np.sqrt(2.0 * np.log(2))) # number of pixels per resolution element
    tell_pad = int(np.ceil(10.0 * pix))

    # This presumes that the data has been interpolated onto the telluric model grid
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
    elif qso:
        pass
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

    # TODO move to sensfunc_telluric and individual routines
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
        #tell_dict_fit = trim_tell_dict(tell_dict, ind_lower, ind_upper)

        # Guess the coefficients by doing a preliminary fit to the sensitivity function with the average telluric behavior
        guess_obj = sensfunc_guess(wave_fit, counts_ps_fit, counts_ps_mask_fit, flam_true_fit, tell_dict, resln_guess,
                                   meta_spec['core']['AIRMASS'], polyorder_vec[iord], func, lower=lower, upper=upper,
                                   debug=debug)
        # Polynomial coefficient bounds
        bounds_obj = [(np.fmin(np.abs(this_coeff)*delta_coeff_bounds[0], minmax_coeff_bounds[0]),
                         np.fmax(np.abs(this_coeff)*delta_coeff_bounds[1], minmax_coeff_bounds[1])) for this_coeff in guess_obj]

        # Set the bounds for the optimization
        bounds_tell = get_bounds_tell(tell_dict, resln_guess, resln_frac_bounds, pix_shift_bounds)
        bounds = bounds_obj + bounds_tell
        # Create the arg_dict
        obj_dict = dict(wave=wave_fit, wave_min=wave_min, wave_max=wave_max, flam_true=flam_true_fit, func=func,
                        polyorder=polyorder_vec[iord])
        arg_dict = dict(ivar=counts_ps_ivar_fit, tell_dict=tell_dict, ind_lower=ind_lower, ind_upper=ind_upper,
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
        telluric_fit = eval_telluric(theta_tell, arg_dict['tell_dict'], arg_dict['ind_lower'], arg_dict['ind_upper'])
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
    #tell_dict_fit = trim_tell_dict(tell_dict, ind_lower, ind_upper)

    #
    # fitfunc, arg_dict, bounds = instantiate_obj_model(model_type, wave_fit, model_params)
    pca_dict = qso_pca.init_pca(pcafile, wave_fit, z_qso, npca)

    flam_norm = qso_guess(flux_fit, flux_ivar_fit, mask_fit, pca_dict, tell_dict, meta_spec['core']['AIRMASS'],
                          resln_guess, tell_norm_thresh, ind_lower, ind_upper)
    # Set the bounds for the PCA and truncate to the right dimension
    coeffs = pca_dict['coeffs'][:,1:npca]
    # Compute the min and max arrays of the coefficients which are not the norm, i.e. grab the coeffs that aren't the first one
    coeff_min = np.amin(coeffs, axis=0)  # only
    coeff_max = np.amax(coeffs, axis=0)
    bounds_z = [(z_qso - delta_zqso, z_qso + delta_zqso)]                # QSO redshift: can vary within delta_zqso
    bounds_flam = [(flam_norm*bounds_norm[0], flam_norm*bounds_norm[1])] # Norm: bounds determined from estimate above
    bounds_coeff = [(i, j) for i, j in zip(coeff_min, coeff_max)]        # Coefficients:  determined from PCA model
    # Set the bounds for the telluric
    bounds_tell = get_bounds_tell(tell_dict, resln_guess, resln_frac_bounds, pix_shift_bounds)
    # Final bounds for the optimizaiton
    bounds =  bounds_z + bounds_flam + bounds_coeff + bounds_tell
    # Create the arg_dict
    obj_dict = dict(npca=npca, pca_dict=pca_dict)
    arg_dict = dict(ivar=flux_ivar_fit, tell_dict=tell_dict,
                    ind_lower=ind_lower, ind_upper=ind_upper,
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
    telluric_fit = eval_telluric(theta_tell, arg_dict['tell_dict'], arg_dict['ind_lower'], arg_dict['ind_upper'])
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

#obj_params = dict(star_type=None, star_mag=None, star_ra=None, star_dec=None, delta_coeff_bounds=(-20.0, 20.0),
#                  minmax_coeff_bounds=(-5.0, 5.0), polyorder=7, func='legendre')


#flam_norm = qso_guess(flux_fit, flux_ivar_fit, mask_fit, pca_dict, tell_dict, meta_spec['core']['AIRMASS'],
#                      resln_guess, tell_norm_thresh, ind_lower, ind_upper)

# Guess the coefficients by doing a preliminary fit to the sensitivity function with the average telluric behavior
#guess_obj = sensfunc_guess(wave_fit, counts_ps_fit, counts_ps_mask_fit, flam_true_fit, tell_dict, resln_guess,
#                           meta_spec['core']['AIRMASS'], polyorder_vec[iord], func, lower=lower, upper=upper,
#                           debug=debug)

