import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#from pypeit.core import coadd1d
import coadd1d_old as coadd1d
from pypeit.core import load


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

gdobjid = ['OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
         'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
         'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001']
# parameters for load_1dspec_to_array
ex_value = 'OPT'
flux_value = True
outfile = 'J1007'

order_vec = np.array([0,1,2,3,4])
norder = np.size(order_vec)
nexp = np.size(gdobjid)
npix = fits.getheader(fnames[0], 0)['NPIX']

# giant array for storing you data
waves = np.zeros((nexp,npix*norder))
fluxes = np.zeros((nexp,npix*norder))
ivars = np.zeros((nexp,npix*norder))
masks = np.zeros((nexp,npix*norder),dtype=bool)
scales = np.zeros((nexp,npix*norder))
weights = np.zeros((nexp,npix*norder))
outmasks = np.zeros((nexp,npix*norder),dtype=bool)

# arrays to store stacked individual order spectra. I set npix_stack = npix*2 to make sure
waves_stack_orders = np.zeros((norder,npix*2))
fluxes_stack_orders = np.zeros((norder,npix*2))
ivars_stack_orders = np.zeros((norder,npix*2))
masks_stack_orders = np.zeros((norder,npix*2),dtype=bool)

# Loop over orders to get the initial stack and scale factor or each order
for ii, iord in enumerate(order_vec):
    gdobjorder = []
    for jj in range(nexp):
        gdobjorder_jj = '{:}-ORDER{:04d}-DET01'.format(gdobjid[jj],iord)
        gdobjorder.append(gdobjorder_jj)

    # Reading data
    waves_iord, fluxes_iord, ivars_iord, masks_iord = load.load_1dspec_to_array(fnames, gdobj=gdobjorder, order=iord,
                                                                           ex_value=ex_value, flux_value=flux_value)
    # store the data to the giant array
    waves[:, ii*npix:(ii+1)*npix] = waves_iord.copy()
    fluxes[:, ii*npix:(ii+1)*npix] = fluxes_iord.copy()
    ivars[:, ii*npix:(ii+1)*npix] = ivars_iord.copy()
    masks[:, ii*npix:(ii+1)*npix] = masks_iord.copy()

    # Get the stacked spectrum for each order
    # Todo: save stacked individual order spectra into one single fits
    if outfile is not None:
        outfile_order = 'spec1d_stack_order{:04d}_{:}'.format(iord,outfile)
    wave_stack_iord, flux_stack_iord, ivar_stack_iord, mask_stack_iord, outmask_iord, weights_iord, scales_iord, rms_sn_iord = \
        coadd1d.combspec(waves_iord, fluxes_iord, ivars_iord, masks_iord, wave_grid_method='loggrid',wave_grid_min=None,
                         wave_grid_max=None, A_pix=None, v_pix=None, samp_fact = 1.0, ref_percentile=20.0,
                         maxiter_scale=5, sigrej=3, scale_method=None, hand_scale=None, sn_max_medscale=2.0,
                         sn_min_medscale=0.5, dv_smooth=10000.0, const_weights=False, maxiter_reject=5, sn_cap=20.0,
                         lower=3.0, upper=3.0, maxrej=None, qafile=None, outfile=outfile_order, debug=False)

    # store individual stacked spectrum
    npix_stack = np.size(wave_stack_iord)
    waves_stack_orders[ii,:npix_stack] = wave_stack_iord
    fluxes_stack_orders[ii,:npix_stack] = flux_stack_iord
    ivars_stack_orders[ii,:npix_stack] = ivar_stack_iord
    masks_stack_orders[ii,:npix_stack] = mask_stack_iord

    # store new masks, scales and weights, all of these arrays are in native wave grid
    scales[:, ii*npix:(ii+1)*npix] = scales_iord.copy()
    weights[:, ii*npix:(ii+1)*npix] = weights_iord.copy()
    outmasks[:, ii*npix:(ii+1)*npix] = outmask_iord.copy()

## re-scale bluer orders to match the reddest order.
# scaling spectrum order by order. We use the reddest order as the reference since slit loss in redder is smaller
for ii in range(norder - 1):
    iord = norder - ii - 1
    sn_iord_iref = fluxes[iord] * (1. / sigs[iord])
    sn_iord_scale = fluxes[iord - 1] * (1. / sigs[iord - 1])
    allok = (sigs[iord - 1, :] > 0) & (sigs[iord, :] > 0) & (sn_iord_iref > SN_MIN_MEDSCALE) & (
    sn_iord_scale > SN_MIN_MEDSCALE)

    if sum(allok) > np.maximum(num_min_pixels, len(wave) * overlapfrac):
        # Ratio
        med_flux = spectra.data['flux'][iord, allok] / spectra.data['flux'][iord - 1, allok]
        # Clip
        mn_scale, med_scale, std_scale = stats.sigma_clipped_stats(med_flux, sigma=nsig, iters=niter)
        med_scale = np.maximum(np.minimum(med_scale, 1.5),0.5)
        spectra.data['flux'][iord - 1, :] *= med_scale
        spectra.data['sig'][iord - 1, :] *= med_scale
        msgs.info('Scaled %s order by a factor of %s'%(iord,str(med_scale)))

        if debug:
            allok_iord = sigs[iord, :] > 0
            allok_iordm1 = sigs[iord-1, :] > 0
            plt.figure(figsize=(12, 6))
            plt.plot(wave[allok], spectra.data['flux'][iord, allok], '-',lw=10,color='0.7', label='Scale region')
            plt.plot(wave[allok_iord], spectra.data['flux'][iord,allok_iord], 'r-', label='reference spectrum')
            plt.plot(wave[allok_iordm1], fluxes_raw[iord - 1,allok_iordm1], 'k-', label='raw spectrum')
            plt.plot(spectra.data['wave'][iord - 1, allok_iordm1], spectra.data['flux'][iord - 1, allok_iordm1], 'b-',
                     label='scaled spectrum')
            mny, medy, stdy = stats.sigma_clipped_stats(fluxes[iord, allok], sigma=nsig, iters=niter)
            plt.ylim([0.1 * medy, 4.0 * medy])
            plt.xlim([np.min(wave[sigs[iord - 1, :] > 0]), np.max(wave[sigs[iord, :] > 0])])
            plt.legend()
            plt.xlabel('wavelength')
            plt.ylabel('Flux')
            plt.show()
    else:
        msgs.warn('Not enough overlap region for sticking different orders.')



from IPython import embed
embed()
