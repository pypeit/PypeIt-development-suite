


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
from telluric import sensfunc_telluric, telluric_qso
import IPython




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
spec1dfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Science/spec1d_XSHOO.2018-11-08T00:16:56.583-Feige110_XShooter_NIR_2018Nov08T001656.583.fits')
#spec1dfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439/vlt_xshooter_nir/Science/spec1d_Feige110.fits')
header = fits.getheader(spec1dfile)
telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
#telgridfile =  os.path.join(dev_path, 'dev_algorithms/sensfunc/TelFit_Paranal_NIR_9800_25000_R25000.fits')
polyorder=5 # changed from 6
outfile = 'Feige110_sens_tell.fits'
#
#sensfunc_telluric(spec1dfile, telgridfile, outfile, polyorder=polyorder, ra=header['RA'], dec=header['DEC'],
#                  star_mag=star_mag, star_type=star_type, debug=False)

### Test telluric star
spec1dfileflux = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0224-4711/Test_tell/spec1d_stack_TELL_B8IV_V5p8.fits')
star_mag  = 5.8
star_type = 'B8'

polyorder=5 # changed from 6
outfile = 'Tell_B8_V5p8_sens_tell.fits'

sensfunc_telluric(spec1dfileflux, telgridfile, outfile, polyorder=polyorder, ra=None, dec=None,
                  star_mag=star_mag, star_type=star_type, ret_flam=True, debug=True)
sys.exit(-1)

spec1dfileflux = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/NIR_Stack/spec1d_stack_Pisco_all.fits')
pcafile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux//qso_pca_1200_3100.pckl')
npca = 8
z_qso = 7.54#6.514 #7.54#6.51
# Create the input mask
hdu = fits.open(spec1dfileflux)
head = hdu[0].header
data = hdu[1].data
wave, flux, flux_ivar, mask = data['OPT_WAVE'], data['OPT_FLAM'], data['OPT_FLAM_IVAR'], data['OPT_MASK']
inmask = (wave > (1.0 + z_qso)*1220.0) & (wave < 3099*(1.0 + z_qso))

#vlt_xshooter_nir = util.load_spectrograph('vlt_xshooter_nir')
#wavegrid = vlt_xshooter_nir.wavegrid(midpoint=True)
# inmask here includes both a Lya mask (unnecessary for lens), BAL mask, and cuts off at the end of the PCA wave grid
#inmask = (wavegrid > (1.0 + z_qso)*1220.0) & ((wavegrid < 10770) | (wavegrid > 11148)) & (wavegrid < 3099*(1+z_qso))
## TODO add an inmask and inmask on its own wavelength grid
telluric_qso(spec1dfileflux, telgridfile, pcafile, npca, z_qso, inmask = inmask, wave_inmask=wave, debug=True)
sys.exit(-1)
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


