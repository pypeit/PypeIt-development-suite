import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#from pypeit.core import coadd1d
#from coadd1d_old import *
from pypeit.core.coadd1d import *
from pypeit.core import load


def order_median_scale(waves_stack_orders, fluxes_stack_orders, ivars_stack_orders, masks_stack_orders,
                       min_good=0.05, maxiters=5, max_factor=10., sigrej=3, debug=False):

    norder = np.shape(waves_stack_orders)[1]
    order_ratios = np.ones(norder)

    ## re-scale bluer orders to match the reddest order.
    # scaling spectrum order by order. We use the reddest order as the reference since slit loss in redder is smaller
    for ii in range(norder - 1):
        iord = norder - ii - 1
        wave_blue, flux_blue, ivar_blue, mask_blue = waves_stack_orders[:, iord-1], fluxes_stack_orders[:, iord-1],\
                                                     ivars_stack_orders[:, iord-1], masks_stack_orders[:, iord-1]

        wave_red_tmp, flux_red_tmp = waves_stack_orders[:, iord], fluxes_stack_orders[:, iord]*order_ratios[iord]
        ivar_red_tmp, mask_red_tmp = ivars_stack_orders[:, iord]*1.0/order_ratios[iord]**2, masks_stack_orders[:, iord]
        wave_mask = wave_red_tmp>1.0
        wave_red, flux_red, ivar_red, mask_red = wave_red_tmp[wave_mask], flux_red_tmp[wave_mask], \
                                                 ivar_red_tmp[wave_mask], mask_red_tmp[wave_mask],

        # interpolate iord-1 (bluer) to iord-1 (redder)
        flux_blue_inter, ivar_blue_inter, mask_blue_inter = interp_spec(wave_red, wave_blue, flux_blue, ivar_blue, mask_blue)

        npix_overlap = np.sum(mask_blue_inter & mask_red)
        percentile_iord = np.fmax(100.0 * (npix_overlap / np.sum(mask_red)-0.05), 10)

        order_ratio_iord = robust_median_ratio(flux_blue_inter, ivar_blue_inter, flux_red, ivar_red, mask=mask_blue_inter,
                                               mask_ref=mask_red, ref_percentile=percentile_iord, min_good=min_good,
                                               maxiters=maxiters, max_factor=max_factor, sigrej=sigrej)

        order_ratios[iord - 1] = np.fmax(np.fmin(order_ratio_iord, max_factor), 1.0/max_factor)
        msgs.info('Scaled {}th order to {}th order by {:}'.format(iord-1, iord, order_ratios[iord-1]))

        if debug:
            plt.figure(figsize=(12, 8))
            plt.plot(wave_red[mask_red], flux_red[mask_red], 'k-', label='reference spectrum')
            plt.plot(wave_blue[mask_blue], flux_blue[mask_blue],color='dodgerblue', lw=3, label='raw spectrum')
            plt.plot(wave_blue[mask_blue], flux_blue[mask_blue]*order_ratios[iord-1], color='r',
                     alpha=0.5, label='re-scaled spectrum')

            med_width = (2.0 * np.ceil(0.03 / 2.0 * wave_blue[mask_blue].size) + 1).astype(int)
            flux_med, ivar_med = median_filt_spec(flux_blue, ivar_blue, mask_blue, med_width)
            ymax = 1.5 * flux_med.max()
            ymin = -0.15 * ymax

            plt.ylim([ymin, ymax])
            plt.xlim([np.min(wave_blue[mask_blue]), np.max(wave_red[mask_red])])
            plt.legend()
            plt.xlabel('wavelength')
            plt.ylabel('Flux')
            plt.show()

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


# parameters for load_1dspec_to_array
ex_value = 'OPT'
flux_value = True

# parameters for new_wave_grid
wave_method = 'loggrid'
A_pix = None
v_pix = None
samp_fact = 1.0
wave_grid_min=None
wave_grid_max=None

# parameters for long_combspec
ref_percentile=20.0
maxiter_scale=5
sigrej=3
scale_method=None
hand_scale=None
sn_max_medscale=2.0
sn_min_medscale=0.5
dv_smooth=10000.0
const_weights=False
maxiter_reject=5
sn_cap=20.0
lower=3.0
upper=3.0
maxrej=None
qafile=None
outfile = 'J1007_NIRES'
debug=False

# parameters for robust_median_ratio
max_factor = 10.0
maxiters = 5
min_good = 0.05

# Loading Echelle data
waves, fluxes, ivars, masks = load.load_1dspec_to_array(fnames, gdobj=gdobj, order=None, ex_value=ex_value, flux_value=flux_value)
scales = np.zeros_like(waves)
weights = np.zeros_like(waves)
outmasks = np.zeros_like(waves,dtype=bool)

# data shape
data_shape = np.shape(waves)
npix = data_shape[0] # detector size in the wavelength direction
norder = data_shape[1]
nexp = data_shape[2]

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
                              upper=upper, maxrej=maxrej, qafile=None, outfile=None, debug=debug)

    # store new masks, scales and weights, all of these arrays are in native wave grid
    scales[:,ii,:] = scales_iord
    weights[:,ii,:] = weights_iord
    outmasks[:,ii,:] = outmask_iord

if debug:
    plt.figure(figsize=(12,8))
    ymin = []
    ymax = []
    for ii in range(norder):
        wave_stack_iord = waves_stack_orders[:, ii]
        flux_stack_iord = fluxes_stack_orders[:, ii]
        ivar_stack_iord = ivars_stack_orders[:, ii]
        mask_stack_iord = masks_stack_orders[:, ii]
        plt.plot(wave_stack_iord[mask_stack_iord],flux_stack_iord[mask_stack_iord])
        #plt.plot(wave_stack_iord[mask_stack_iord],1.0/np.sqrt(ivar_stack_iord[mask_stack_iord]))

        med_width = (2.0 * np.ceil(0.1 / 2.0 * np.size(wave_stack_iord[mask_stack_iord])) + 1).astype(int)
        flux_med, ivar_med = median_filt_spec(flux_stack_iord, ivar_stack_iord, mask_stack_iord, med_width)
        ymax.append(1.5 * flux_med.max())
        ymin.append(-0.2 * flux_med.max())
    plt.xlim([np.min(waves_stack_orders[masks_stack_orders]),np.max(waves_stack_orders[masks_stack_orders])])
    plt.ylim([np.min(ymin),np.max(ymax)])
    plt.xlabel('Wavelength ($\\rm\\AA$)')
    plt.ylabel('Flux')
    plt.show()

order_ratios = order_median_scale(waves_stack_orders, fluxes_stack_orders, ivars_stack_orders, masks_stack_orders,
                       min_good=min_good, maxiters=maxiters, max_factor=max_factor, sigrej=sigrej, debug=debug)

from IPython import embed
embed()

# apply order_ratios to the scales array: order_ratio*scale
scales_new = np.copy(scales)
for ii in range(norder):
    scales_new[ii*npix:(ii+1)*npix,:] *= order_ratios[ii]

fluxes_scale = fluxes * scales_new
ivars_scale = ivars * 1.0/scales_new**2

# Todo: save stacked individual order spectra into one single fits
if outfile is not None:
    outfile_order = 'spec1d_stack_order{:04d}_{:}'.format(ii, outfile)

wave_grid = new_wave_grid(waves, wave_method=wave_method, wave_grid_min=None, wave_grid_max=None,
                          A_pix=None, v_pix=None, samp_fact=1.0)

#wave_ref, flux_ref, ivar_ref, mask_ref, nused = compute_stack(waves, fluxes_scale, ivars_scale, masks, wave_grid, weights)
wave_stack, flux_stack, ivar_stack, mask_stack, outmask, weights, scales = spec_reject(
    waves, fluxes_scale, ivars_scale, masks, weights, wave_grid, debug=True)
write_to_fits(wave_stack, flux_stack, ivar_stack, mask_stack, 'J1007_NIRES_Coadd.fits', clobber=True)


from IPython import embed
embed()

for ii in range(norder):
    plt.plot(waves_stack_orders[:, ii][masks_stack_orders[:, ii]], fluxes_stack_orders[:, ii][masks_stack_orders[:, ii]])
plt.ylim([ymin, ymax])
plt.show()

for ii in range(norder):
    plt.plot(waves_stack_orders[:, ii][masks_stack_orders[:, ii]], order_ratios[ii]*fluxes_stack_orders[:, ii][masks_stack_orders[:, ii]])
plt.ylim([ymin, ymax])
plt.show()