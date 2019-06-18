
def update_bounds(bounds, delta_coeff, coeff):

    bounds_new = [(this_coeff * delta_coeff[0], this_coeff * delta_coeff[1]) for this_coeff in coeff]
    bounds_tell = bounds[len(coeff):]
    bounds_new.extend(bounds_tell)
    return bounds_new




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
