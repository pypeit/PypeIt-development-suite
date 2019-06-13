import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#from pypeit.core import coadd1d
#from coadd1d_old import *
from pypeit.core.coadd1d import *
from pypeit.core import load, save

def ech_combspec(fnames, objids, ex_value='OPT', flux_value=True, wave_method='loggrid', A_pix=None, v_pix=None,
                 samp_fact = 1.0, wave_grid_min=None, wave_grid_max=None, ref_percentile=20.0, maxiter_scale=5,
                 sigrej=3, scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5,
                 dv_smooth=10000.0, const_weights=False, maxiter_reject=5, sn_cap=20.0, lower=3.0, upper=3.0,
                 maxrej=None, max_factor=10.0, maxiters=5, min_good=0.05, phot_scale_dicts=None,
                 qafile=None, outfile = None, debug=False, show=True):

    # Loading Echelle data
    waves, fluxes, ivars, masks, header = load.load_1dspec_to_array(fnames, gdobj=gdobj, order=None, ex_value=ex_value,
                                                                    flux_value=flux_value)

    # data shape
    data_shape = np.shape(waves)
    npix = data_shape[0] # detector size in the wavelength direction
    norder = data_shape[1]
    nexp = data_shape[2]

    # create some arrays
    scales = np.zeros_like(waves)
    weights = np.zeros_like(waves)
    outmasks = np.zeros_like(waves,dtype=bool)

    # output name root for fits and QA plots
    if outfile is None:
        outfile = header['TARGET']
    outfile_order = 'spec1d_order_{:}'.format(outfile)
    outfile_stack = 'spec1d_stack_{:}'.format(outfile)
    outfile_giant_stack = 'spec1d_giant_stack_{:}'.format(outfile)

    if qafile is None:
        qafile = header['TARGET']
    qafile_stack = 'spec1d_stack_{:}'.format(outfile)
    qafile_giant_stack = 'spec1d_giant_stack_{:}'.format(outfile)

    # Generate a giant wave_grid
    wave_grid = new_wave_grid(waves, wave_method=wave_method, wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max,
                              A_pix=A_pix, v_pix=v_pix, samp_fact=samp_fact)

    # Arrays to store stacked individual order spectra.
    waves_stack_orders = np.zeros((np.size(wave_grid)-1, norder))
    fluxes_stack_orders = np.zeros_like(waves_stack_orders)
    ivars_stack_orders = np.zeros_like(waves_stack_orders)
    masks_stack_orders = np.zeros_like(waves_stack_orders,dtype=bool)

    # Loop over orders to get the initial stacks of each order
    for ii in range(norder):

        # get the slice of iord spectra of all exposures
        waves_iord, fluxes_iord, ivars_iord, masks_iord = waves[:,ii,:], fluxes[:,ii,:], ivars[:,ii,:], masks[:,ii,:]

        # Get the stacked spectrum for each order
        waves_stack_orders[:, ii], fluxes_stack_orders[:, ii], ivars_stack_orders[:, ii],  masks_stack_orders[:, ii], \
        outmask_iord, weights_iord, scales_iord, rms_sn_iord = \
                    long_combspec(wave_grid, waves_iord, fluxes_iord, ivars_iord, masks_iord, ref_percentile=ref_percentile,
                                  maxiter_scale=maxiter_scale, sigrej=sigrej, scale_method=scale_method, hand_scale=hand_scale,
                                  sn_max_medscale=sn_max_medscale, sn_min_medscale=sn_min_medscale, dv_smooth=dv_smooth,
                                  const_weights=const_weights, maxiter_reject=maxiter_reject, sn_cap=sn_cap, lower=lower,
                                  upper=upper, maxrej=maxrej, title='Stack_{:}th order'.format(ii), qafile=None,
                                  outfile=None, debug=debug, show=show)

        # store new masks, scales and weights, all of these arrays are in native wave grid
        scales[:,ii,:] = scales_iord
        weights[:,ii,:] = weights_iord
        outmasks[:,ii,:] = outmask_iord

    # scale different orders based on overlapped regions using median scale
    fluxes_stack_orders_scale, ivars_stack_orders_scale, order_ratios = \
        order_median_scale(waves_stack_orders, fluxes_stack_orders, ivars_stack_orders, masks_stack_orders,
                           min_good=min_good, maxiters=maxiters, max_factor=max_factor, sigrej=sigrej, debug=debug)

    ## Stack with the first method: combine the stacked individual order spectra directly
    # Get weights for individual order stacks
    rms_sn_stack, weights_stack = sn_weights(waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale,
                                             masks_stack_orders, dv_smooth=dv_smooth, const_weights=const_weights, verbose=True)

    # Reject and stack
    wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused = spec_reject_comb(
        wave_grid, waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale, masks_stack_orders,
        weights_stack, sn_cap=sn_cap, lower=lower, upper=upper, maxrej=maxrej, maxiter_reject=maxiter_reject,
        title='Order merged stacked spectra', qafile=qafile_stack, debug=debug, show=show)

    # Save stacked individual order spectra
    save.save_coadd1d_to_fits(waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale, masks_stack_orders,
                              header=header, outfile=outfile_order, ex_value = ex_value, overwrite=True)
    save.save_coadd1d_to_fits(wave_stack, flux_stack, ivar_stack, mask_stack, header=header, outfile=outfile_stack,
                              ex_value = ex_value, overwrite=True)

    # apply order_ratios to the scales array: order_ratio*scale
    scales_new = np.copy(scales)
    for ii in range(norder):
        scales_new[:,ii,:] *= order_ratios[ii]

    fluxes_scale = fluxes * scales_new
    ivars_scale = ivars * 1.0/scales_new**2

    # reshaping 3D arrays (npix, norder, nexp) to 2D arrays (npix, norder*nexp)
    waves_2d = np.reshape(waves,(npix, norder*nexp))
    fluxes_2d = np.reshape(fluxes_scale, np.shape(waves_2d))
    ivars_2d = np.reshape(ivars_scale, np.shape(waves_2d))
    masks_2d = np.reshape(masks, np.shape(waves_2d))
    weights_2d = np.reshape(weights, np.shape(waves_2d))

    wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack, outmask_giant_stack, nused_giant_stack = \
        spec_reject_comb(wave_grid, waves_2d, fluxes_2d, ivars_2d, masks_2d, weights_2d, sn_cap=sn_cap, lower=lower,
                         upper=upper, maxrej=maxrej, maxiter_reject=maxiter_reject, title='Finale stacked spectra',
                         qafile=qafile_stack, debug=debug, show=show)

    save.save_coadd1d_to_fits(wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack,
                              header=header, outfile=outfile_giant_stack, ex_value=ex_value, overwrite=True)

    return wave_stack, flux_stack, ivar_stack, mask_stack

datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/NIRES/NIRES_May19/Science/')
fnames = [datapath + 'spec1d_flux_s190519_0037-J1007+2115_NIRES_2019May19T055221.895.fits',
          datapath + 'spec1d_flux_s190519_0038-J1007+2115_NIRES_2019May19T055923.665.fits',
          datapath + 'spec1d_flux_s190519_0041-J1007+2115_NIRES_2019May19T062048.865.fits',
          datapath + 'spec1d_flux_s190519_0042-J1007+2115_NIRES_2019May19T062750.635.fits',
          datapath + 'spec1d_flux_s190519_0045-J1007+2115_NIRES_2019May19T064943.885.fits',
          datapath + 'spec1d_flux_s190519_0046-J1007+2115_NIRES_2019May19T065646.165.fits',
          datapath + 'spec1d_flux_s190519_0049-J1007+2115_NIRES_2019May19T071920.215.fits',
          datapath + 'spec1d_flux_s190519_0050-J1007+2115_NIRES_2019May19T072621.985.fits',
          datapath + 'spec1d_flux_s190519_0053-J1007+2115_NIRES_2019May19T074819.315.fits',
          datapath + 'spec1d_flux_s190519_0054-J1007+2115_NIRES_2019May19T075521.595.fits',
          datapath + 'spec1d_flux_s190519_0057-J1007+2115_NIRES_2019May19T081918.265.fits',
          datapath + 'spec1d_flux_s190519_0058-J1007+2115_NIRES_2019May19T082620.545.fits']

gdobj = ['OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
         'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
         'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001']

wave_stack, flux_stack, ivar_stack, mask_stack = ech_combspec(fnames, gdobj)
