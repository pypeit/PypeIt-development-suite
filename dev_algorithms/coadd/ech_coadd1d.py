import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#from pypeit.core import coadd1d
#from coadd1d_old import *
from pypeit.core.coadd1d import *
from pypeit.core import load


def load_ext_to_array(hdulist, ext_id, ex_value='OPT', flux_value=True):
    '''
    It will be called by load_1dspec_to_array.
    Load one-d spectra from ext_id in the hdulist
    Args:
        hdulist: FITS HDU list
        ext_id: extension name, i.e., 'SPAT1073-SLIT0001-DET03', 'OBJID0001-ORDER0003', 'OBJID0001-ORDER0002-DET01'
        ex_value: 'OPT' or 'BOX'
        flux_value: if True load fluxed data, else load unfluxed data
    Returns:
         wave, flux, ivar, mask
    '''

    if (ex_value != 'OPT') and (ex_value != 'BOX'):
        msgs.error('{:} is not recognized. Please change to either BOX or OPT.'.format(ex_value))

    # get the order/slit information
    ntrace0 = np.size(hdulist)-1
    idx_names = []
    for ii in range(ntrace0):
        idx_names.append(hdulist[ii+1].name) # idx name

    # Initialize ext
    ext = None
    for indx in (idx_names):
        if ext_id in indx:
            ext = indx
    if ext is None:
        msgs.error('Can not find extension {:}.'.format(ext_id))
    else:
        hdu_iexp = hdulist[ext]

    wave = hdu_iexp.data['{:}_WAVE'.format(ex_value)]
    mask = hdu_iexp.data['{:}_MASK'.format(ex_value)]
    if flux_value:
        flux = hdu_iexp.data['{:}_FLAM'.format(ex_value)]
        ivar = hdu_iexp.data['{:}_FLAM_IVAR'.format(ex_value)]
    else:
        msgs.warn('Loading unfluxed spectra')
        flux = hdu_iexp.data['{:}_COUNTS'.format(ex_value)]
        ivar = hdu_iexp.data['{:}_COUNTS_IVAR'.format(ex_value)]

    return wave, flux, ivar, mask

def load_1dspec_to_array(fnames, gdobj=None, order=None, ex_value='OPT', flux_value=True):
    '''
    Load the spectra from the 1d fits file into arrays.
    If Echelle, you need to specify which order you want to load.
    It can NOT load all orders for Echelle data.
    Args:
        fnames:
        gdobj:
        extensions:
        order: set to None if longslit data
        ex_value:
        flux_value:
    Returns:
        waves:
        fluxes:
        ivars:
        masks:
    '''

    # read in the first fits file
    nexp = np.size(fnames)
    if nexp == 1:
        fname0 = fnames.copy()
    else:
        fname0 = fnames[0]
    hdulist = fits.open(fname0)
    npix = hdulist[0].header['NPIX']
    pypeline = hdulist[0].header['PYPELINE']

    # get the order/slit information
    ntrace0 = np.size(hdulist)-1
    idx_names = []
    idx_objids = []
    idx_orders = []
    for ii in range(ntrace0):
        idx_names.append(hdulist[ii+1].name) # idx name
        idx_objids.append(hdulist[ii+1].name.split('-')[0])
        idx_orders.append(int(hdulist[ii+1].name.split('-')[1][5:])) # slit ID or order ID

    if pypeline == "Echelle":
        order_vec = np.unique(idx_orders)
        norder = np.size(order_vec)
    else:
        order_vec = None
        norder = 1

    ## Loading data from a single fits file
    if nexp == 1:
        # initialize arrays
        if (order is None) and (pypeline == "Echelle"):
            waves = np.zeros((npix, norder))
            fluxes = np.zeros_like(waves)
            ivars = np.zeros_like(waves)
            masks = np.zeros_like(waves, dtype=bool)

            for ii, iord in enumerate(order_vec):
                ext_id = gdobj+'-ORDER{:04d}'.format(iord)
                wave_iord, flux_iord, ivar_iord, mask_iord = load_ext_to_array(hdulist, ext_id, ex_value=ex_value,
                                                                               flux_value=flux_value)
                waves[:,ii] = wave_iord
                fluxes[:,ii] = flux_iord
                ivars[:,ii] = ivar_iord
                masks[:,ii] = mask_iord
        else:
            if pypeline == "Echelle":
                ext_id = gdobj+'-ORDER{:04d}'.format(order)
            else:
                ext_id = gdobj
            waves, fluxes, ivars, masks = load_ext_to_array(hdulist, ext_id, ex_value=ex_value, flux_value=flux_value)

    ## Loading data from a list of fits files
    else:
        # initialize arrays
        if (order is None) and (pypeline == "Echelle"):
            # store all orders into one single array
            waves = np.zeros((npix * norder, nexp))
        else:
            # store a specific order or longslit
            waves = np.zeros((npix, nexp))
        fluxes = np.zeros_like(waves)
        ivars = np.zeros_like(waves)
        masks = np.zeros_like(waves,dtype=bool)

        for iexp in range(nexp):
            hdulist_iexp = fits.open(fnames[iexp])
            if (order is None) and (pypeline == "Echelle"):
                wave = np.zeros(npix*norder)
                flux = np.zeros_like(wave)
                ivar = np.zeros_like(wave)
                mask = np.zeros_like(wave,dtype=bool)
                for ii, iord in enumerate(order_vec):
                    ext_id = gdobj[iexp]+'-ORDER{:04d}'.format(iord)
                    wave_iord, flux_iord, ivar_iord, mask_iord = load_ext_to_array(hdulist_iexp, ext_id, ex_value=ex_value,
                                                                                   flux_value=flux_value)
                    wave[npix*ii:npix*(ii+1)] = wave_iord
                    flux[npix*ii:npix*(ii+1)] = flux_iord
                    ivar[npix*ii:npix*(ii+1)] = ivar_iord
                    mask[npix*ii:npix*(ii+1)] = mask_iord
            else:
                if pypeline == "Echelle":
                    ext_id = gdobj[iexp]+'-ORDER{:04d}'.format(order)
                else:
                    ext_id = gdobj[iexp]
                wave, flux, ivar, mask = load_ext_to_array(hdulist_iexp, ext_id, ex_value=ex_value, flux_value=flux_value)
            waves[:, iexp] = wave
            fluxes[:, iexp] = flux
            ivars[:, iexp] = ivar
            masks[:, iexp] = mask

    return waves,fluxes,ivars,masks


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

outfile = 'J1007'
norder = 5


from IPython import embed
embed()


# Load data
waves, fluxes, ivars, masks = load_1dspec_to_array(fnames, gdobj=gdobj, order=None, ex_value=ex_value, flux_value=flux_value)
scales = np.zeros_like(waves)
weights = np.zeros_like(waves)
outmasks = np.zeros_like(waves,dtype=bool)


npix = int(np.shape(waves)[0]/norder)
# arrays to store stacked individual order spectra. I set npix_stack = npix*2 to make sure
waves_stack_orders = np.zeros((npix*2, norder))
fluxes_stack_orders = np.zeros((npix*2, norder))
ivars_stack_orders = np.zeros((npix*2, norder))
masks_stack_orders = np.zeros((npix*2, norder),dtype=bool)

# Loop over orders to get the initial stack and scale factor or each order
for ii in range(norder):

    waves_iord, fluxes_iord, ivars_iord, masks_iord = waves[npix*ii:npix*(ii+1),:], fluxes[npix*ii:npix*(ii+1),:], \
                                                      ivars[npix*ii:npix*(ii+1),:], masks[npix*ii:npix*(ii+1),:]

    # Get the stacked spectrum for each order
    # Todo: save stacked individual order spectra into one single fits
    if outfile is not None:
        outfile_order = 'spec1d_stack_order{:04d}_{:}'.format(ii,outfile)
    wave_stack_iord, flux_stack_iord, ivar_stack_iord, mask_stack_iord, outmask_iord, weights_iord, scales_iord, rms_sn_iord = \
                combspec(waves_iord, fluxes_iord, ivars_iord, masks_iord, wave_grid_method='loggrid',wave_grid_min=None,
                         wave_grid_max=None, A_pix=None, v_pix=None, samp_fact = 1.0, ref_percentile=20.0,
                         maxiter_scale=5, sigrej=3, scale_method=None, hand_scale=None, sn_max_medscale=2.0,
                         sn_min_medscale=0.5, dv_smooth=10000.0, const_weights=False, maxiter_reject=5, sn_cap=20.0,
                         lower=3.0, upper=3.0, maxrej=None, qafile=None, outfile=outfile_order, debug=False)

    # store individual stacked spectrum
    npix_stack = np.size(wave_stack_iord)
    waves_stack_orders[:npix_stack, ii] = wave_stack_iord
    fluxes_stack_orders[:npix_stack, ii] = flux_stack_iord
    ivars_stack_orders[:npix_stack, ii] = ivar_stack_iord
    masks_stack_orders[:npix_stack, ii] = mask_stack_iord

    # store new masks, scales and weights, all of these arrays are in native wave grid
    scales[ii*npix:(ii+1)*npix, :] = scales_iord.copy()
    weights[ii*npix:(ii+1)*npix, :] = weights_iord.copy()
    outmasks[ii*npix:(ii+1)*npix, :] = outmask_iord.copy()


debug = True
max_factor = 10.0
wave_method = 'loggrid'


order_ratios = np.ones(norder)
## re-scale bluer orders to match the reddest order.
# scaling spectrum order by order. We use the reddest order as the reference since slit loss in redder is smaller
for ii in range(norder - 1):
    iord = norder - ii - 1

    wave_blue, flux_blue, ivar_blue, mask_blue = waves_stack_orders[:, iord-1], fluxes_stack_orders[:, iord-1],\
                                                 ivars_stack_orders[:, iord-1], masks_stack_orders[:, iord-1]

    wave_red, flux_red, ivar_red, mask_red = waves_stack_orders[:, iord], fluxes_stack_orders[:, iord]*order_ratios[iord],\
                                             ivars_stack_orders[:, iord]*1.0/order_ratios[iord]**2, masks_stack_orders[:, iord]

    # interpolate iord-1 (bluer) to iord-1 (redder)
    flux_tmp, ivar_tmp, mask_tmp = interp_spec(wave_red, wave_blue, flux_blue, ivar_blue, mask_blue)

    npix_overlap = np.sum(mask_tmp & mask_red)
    percentile_iord = np.fmax(100.0 * (npix_overlap / np.sum(mask_red)-0.05), 10)

    order_ratio_iord = robust_median_ratio(flux_tmp, ivar_tmp, flux_red, ivar_red, mask=mask_tmp,
                                               mask_ref=mask_red, ref_percentile=percentile_iord, min_good=0.05,
                                               maxiters=5, max_factor=10.0, sigrej=3.0)

    order_ratios[iord - 1] = np.fmax(np.fmin(order_ratio_iord, max_factor), 1.0/max_factor)
    msgs.info('Scaled {}th order to {}th order by {:}'.format(order_vec[iord-1],order_vec[iord],order_ratios[iord-1]))

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

# apply order_ratios to the scales array: order_ratio*scale
scales_new = np.copy(scales)
for ii in range(norder):
    scales_new[ii*npix:(ii+1)*npix,:] *= order_ratios[ii]

fluxes_scale = fluxes * scales_new
ivars_scale = ivars * 1.0/scales_new**2

wave_grid = new_wave_grid(waves, wave_method=wave_method, wave_grid_min=None, wave_grid_max=None,
                          A_pix=None, v_pix=None, samp_fact=1.0)

#wave_ref, flux_ref, ivar_ref, mask_ref, nused = compute_stack(waves, fluxes_scale, ivars_scale, masks, wave_grid, weights)
wave_stack, flux_stack, ivar_stack, mask_stack, outmask, weights, scales = spec_reject(
    waves, fluxes_scale, ivars_scale, masks, weights, wave_grid, debug=True)
write_to_fits(wave_stack, flux_stack, ivar_stack, mask_stack, 'J1007_NIRES_Coadd.fits', clobber=True)


def spec_reject(waves, fluxes, ivars, masks, weights, wave_grid, sn_cap=20.0, lower=3.0, upper=3.0, maxrej=None,
                maxiter_reject=5, qafile=None, debug=False):

    iIter = 0
    qdone = False
    thismask = np.copy(masks)
    while (not qdone) and (iIter < maxiter_reject):
        wave_stack, flux_stack, ivar_stack, mask_stack, nused = compute_stack(
            waves, fluxes, ivars, thismask, wave_grid, weights)
        flux_stack_nat, ivar_stack_nat, mask_stack_nat = interp_spec(
            waves, wave_stack, flux_stack, ivar_stack,mask_stack)
        rejivars, sigma_corrs, outchi, maskchi = update_errors(waves, fluxes, ivars, thismask,
                                                               flux_stack_nat, ivar_stack_nat, mask_stack_nat,
                                                               sn_cap=sn_cap)
        thismask, qdone = pydl.djs_reject(fluxes, flux_stack_nat, outmask=thismask,inmask=masks, invvar=rejivars,
                                          lower=lower,upper=upper, maxrej=maxrej, sticky=False)
        # print out how much was rejected
        for iexp in range(nexp):
            thisreject = thismask[:, iexp]
            nrej = np.sum(np.invert(thisreject))
            if nrej > 0:
                msgs.info("Rejecting {:d} pixels in exposure {:d}".format(nrej, iexp))

        iIter = iIter +1

    if (iIter == maxiter_reject) & (maxiter_reject != 0):
        msgs.warn('Maximum number of iterations maxiter={:}'.format(maxiter_reject) + ' reached in combspec')
    outmask = np.copy(thismask)

    # Compute the final stack using this outmask
    wave_stack, flux_stack, ivar_stack, mask_stack, nused = compute_stack(
        waves, fluxes, ivars, outmask, wave_grid, weights)
    # Used only for plotting below
    flux_stack_nat, ivar_stack_nat, mask_stack_nat = interp_spec(waves, wave_stack, flux_stack, ivar_stack, mask_stack)
    if debug:
        for iexp in range(nexp):
            # plot the residual distribution
            renormalize_errors_qa(outchi[:, iexp], maskchi[:, iexp], sigma_corrs[iexp])
            # plot the rejections for each exposures
            coadd_iexp_qa(waves[:, iexp], fluxes[:, iexp], ivars[:, iexp], flux_stack_nat[:, iexp],
                          ivar_stack_nat[:, iexp], mask=outmask[:, iexp], mask_stack=mask_stack_nat[:, iexp],
                          qafile=None, debug=debug)

    # Plot the final coadded spectrum
    if debug:
        weights_qa(waves, weights, outmask)
        coadd_qa(wave_stack,flux_stack,ivar_stack, nused, mask=mask_stack, qafile=qafile, debug=debug)

    return wave_stack, flux_stack, ivar_stack, mask_stack, outmask, weights, scales

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