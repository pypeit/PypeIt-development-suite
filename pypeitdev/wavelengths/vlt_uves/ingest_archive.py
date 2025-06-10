import os
import sys
import numpy as np
import itertools
from astropy.table import Table, vstack
from astropy.io import fits
from pkg_resources import resource_filename
from matplotlib import pyplot as plt
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.wavecal import templates
from pypeit.core.wavecal import wvutils
from pypeit.core import arc
from pypeit.core.fitting import robust_fit
from pypeit.core import coadd
from pypeit.core import fitting
from pypeit.core.wavecal import autoid, waveio, wv_fitting
from pypeit.core.wavecal.wvutils import get_xcorr_arc, xcorr_shift
from pypeit import utils
from pypeit import msgs
from pypeit import wavecalib
from astropy import table
from scipy import interpolate
from IPython import embed
from astropy import constants as const

c_kms = const.c.to('km/s').value


def load_archive(outfile, n_final=4, func='legendre'):
    """
    Load the pypeit format archived templates that were created with the `esorex_to_pypeit` function.

    Args:
        outfile (str):
            Name of the output file

    """
    # load the list of pypeit templates
    ptempl_table_file = os.path.join(os.getenv('PYPEIT_DEV'), 'dev_algorithms', 'wavelengths', 'vlt_uves',
                                     'uves_templ_pypeit.dat')
    ptbl = Table.read(ptempl_table_file, format='ascii')
    p_nrows = len(ptbl)

    # recomputing the final order min, max, order_vec, norders
    order_min = ptbl['IOrder'].min()
    order_max = ptbl['EOrder'].max()
    order_vec = np.arange(order_min, order_max + 1, 1)
    norders = order_vec.size
    nspec = 3000

    # xmin, xmax for wavelength vs pixel fits
    xmin, xmax = 0.0, 1.0
    params = Table([[order_min],[order_max],[norders],[n_final],[nspec],[func],[xmin],[xmax]]
                ,names=('order_min','order_max','norders','n_final','nspec','func','xmin','xmax'))

    # create new final table
    tot_nrows = p_nrows
    final_table = empty_design_table(tot_nrows, norders, n_final=n_final)
    final_table['wave'] = np.zeros((tot_nrows, norders, nspec))
    final_table['arcspec'] = np.zeros((tot_nrows, norders, nspec))

    # Now deal with the pypeit format
    # load the pypeit templates (they are WaveCalib files)
    for irow in np.arange(p_nrows):
        ifinal_row = irow
        templ_file = os.path.join(os.getenv('PYPEIT_DEV'), 'dev_algorithms', 'wavelengths', 'vlt_uves',
                                  ptbl[irow]['Name'])
        waveCalib = wavecalib.WaveCalib.from_file(templ_file, chk_version=False)
        # this is the order vector available for this wavecalib file
        this_order_vec_raw = np.arange(ptbl[irow]['IOrder'], ptbl[irow]['EOrder'] + 1, 1)
        # select the orders that we want to use
        if ptbl[irow]['Spatids'] == 'all':
            igood = np.ones(waveCalib.spat_ids.size, dtype=bool)
            if this_order_vec_raw.size != waveCalib.spat_ids.size:
                msgs.error('the number of order determined using IOrder and EOrder does not match the number of '
                           'orders in the WaveCalib file')
        elif isinstance(ptbl[irow]['Spatids'], (list, tuple, np.integer)):
            spat_ids = np.atleast_1d(ptbl[irow]['Spatids'])
            igood = np.isin(waveCalib.spat_ids, spat_ids)
        else:
            msgs.error('Unrecognized format for Spatids')
        this_order_vec = this_order_vec_raw[igood]
        this_arc = np.array([arc.resize_spec(ww, nspec) for ww in waveCalib.arc_spectra.T])
        this_wave = np.array([arc.resize_spec(wvfit.wave_soln, nspec) for wvfit in waveCalib.wv_fits])
        nsolns = this_order_vec_raw.size
        nsolns_good = np.sum(igood)
        indx = this_order_vec_raw - order_min
        indx_good = this_order_vec - order_min
        # Information for the file is stored for convenience, although this is redundant with the arrays below
        final_table['filename'][ifinal_row] = ptbl[irow]['Name']
        final_table['nsolns'][ifinal_row] = nsolns
        final_table['nsolns_good'][ifinal_row] = nsolns_good
        final_table['bluest_order'][ifinal_row] = this_order_vec_raw[-1]
        final_table['bluest_good_order'][ifinal_row] = this_order_vec[-1]
        final_table['reddest_order'][ifinal_row] = this_order_vec_raw[0]
        final_table['reddest_good_order'][ifinal_row] = this_order_vec[0]
        final_table['xdisp_file'][ifinal_row] = ptbl[irow]['XDISP']
        final_table['ech_angle_file'][ifinal_row] = ptbl[irow]['ECH']
        final_table['xd_angle_file'][ifinal_row] = ptbl[irow]['XDAng']
        final_table['det_file'][ifinal_row] = ptbl[irow]['Chip']
        # Arrays (nfile, norders)
        final_table['order'][ifinal_row, indx] = this_order_vec_raw
        final_table['populated'][ifinal_row, indx] = True
        final_table['populated_and_good'][ifinal_row, indx_good] = True
        final_table['ech_angle'][ifinal_row, indx] = ptbl[irow]['ECH']
        final_table['xd_angle'][ifinal_row, indx] = ptbl[irow]['XDAng']
        final_table['xdisp'][ifinal_row, indx] = ptbl[irow]['XDISP']
        final_table['det'][ifinal_row, indx] = ptbl[irow]['Chip']
        final_table['binspec'][ifinal_row, indx] = ptbl[irow]['Rbin']
        final_table['lambda_cen'][ifinal_row, indx] = np.median(this_wave, axis=1)
        final_table['wave'][ifinal_row, indx, :] = this_wave
        final_table['arcspec'][ifinal_row, indx, :] = this_arc
        # Fit the wavelengths
        xnspecmin1 = float(nspec - 1)
        xvec = np.arange(nspec) / xnspecmin1
        this_coeff_array = np.zeros((nsolns_good, n_final + 1))
        for ii, iwave in enumerate(this_wave[igood, :]):
            pypeitFit = fitting.robust_fit(xvec, iwave, n_final, function=function, maxiter=10,
                                           lower=1e10, upper=1e10, maxrej=0, sticky=True,
                                           minx=xmin, maxx=xmax, weights=None)
            this_coeff_array[ii, :] = pypeitFit.fitc
        final_table['coeff'][ifinal_row, indx_good, :] = this_coeff_array

    # Write out to multi-extension fits
    print(f'Writing UVES PypeIt wv_calib archive to file: {outfile}')
    hdu_param = fits.BinTableHDU(params.as_array())
    hdu_table = fits.BinTableHDU(final_table.as_array())

    hdulist = fits.HDUList()
    hdulist.append(hdu_param)
    hdulist.append(hdu_table)
    hdulist.writeto(outfile, overwrite=True)


def fit_wvcalib_vs_angles(arxiv_file, outfile, func='legendre',
                         ech_nmax = 3, ech_coeff_fit_order_min=1, ech_coeff_fit_order_max=2,
                         xd_reddest_fit_polyorder=2, sigrej=3.0, maxrej=1, debug=False):
    """
    Fit the coefficients of the wavelength solution vs. the ECH angle. Also fit the bluest_order as a function of XD angle
    for each XDISP

    Args:
        arxiv_file (str):
             File containing the XIDL archive of the HIRES calibration
        outfile (str):
            File to write the output to.
        func (str):
            Function to fit the coefficients.
        ech_nmax (int):
            Polynomial coefficients from nmax to n_final + 1 will be fit with the lower order coeff_fit_order_min
        ech_coeff_fit_order_min (int):
            Polynomial order to fit the set of coefficients from nmax to n_final+1. These coefficients have
            weaker trends with ech_angle so we fit them with lower order.
        ech_coeff_fit_order_max (int):
            Polynomial order to fit the set of coefficients from 0 to nmax. These coefficients show stronger
            trends with ech_angle so we fit them with higher order.
        sigrej (float):
             Rejection threshold for the coefficient fits.
        maxrej (int):
             Maximum number of rejections to allow for each iteration of the coefficient fits rejection.
        debug (bool):
             If True, show plots illustring the fits

    Returns:
        None
    """


    arxiv_params = Table.read(arxiv_file, hdu=1)[0]
    arxiv = Table.read(arxiv_file, hdu=2)

    ech_angle_fit_params, ech_angle_fit_coeffs = fit_coeffs_vs_ech_angle(
        arxiv_params, arxiv, func=func, nmax = ech_nmax, coeff_fit_order_min=ech_coeff_fit_order_min,
        coeff_fit_order_max=ech_coeff_fit_order_max, sigrej=sigrej, maxrej=maxrej, debug=debug)

    xd_angle_fit_params, xd_angle_fit_coeffs = fit_reddest_vs_xd_angle(
        arxiv, polyorder=xd_reddest_fit_polyorder, func=func, sigrej=sigrej, maxrej=maxrej, debug=debug)

    fit_params = table.hstack((ech_angle_fit_params, xd_angle_fit_params))

    hdulist = fits.HDUList()
    hdulist.append(fits.BinTableHDU(fit_params.as_array()))  # hdu = 1
    hdulist.append(fits.ImageHDU(np.array(ech_angle_fit_coeffs)))  # hdu = 2
    hdulist.append(fits.ImageHDU(np.array(xd_angle_fit_coeffs)))  # hdu = 3
    hdulist.writeto(outfile, overwrite=True)


def fit_coeffs_vs_ech_angle(arxiv_params, arxiv, func='legendre', nmax = 3, coeff_fit_order_min=1, coeff_fit_order_max=2,
                            sigrej=3.0, maxrej=1, debug=False):
    """
    Fit the coefficients of the wavelength solution vs. the ECH angle. Called by fit_coeffs_vs_angles

    Args:
        arxiv_file (str):
             File containing the XIDL archive of the HIRES calibration
        outfile (str):
            File to write the output to.
        func (str):
            Function to fit the coefficients.
        nmax (int):
            Polynomial coefficients from nmax to n_final + 1 will be fit with the lower order coeff_fit_order_min
        coeff_fit_order_min (int):
            Polynomial order to fit the set of coefficients from nmax to n_final+1. These coefficients have
            weaker trends with ech_angle so we fit them with lower order.
        coeff_fit_order_max (int):
            Polynomial order to fit the set of coefficients from 0 to nmax. These coefficients show stronger
            trends with ech_angle so we fit them with higher order.
        sigrej (float):
             Rejection threshold for the coefficient fits.
        maxrej (int):
             Maximum number of rejections to allow for each iteration of the coefficient fits rejection.
        debug (bool):
             If True, show plots illustring the fits

    Returns:
        ech_angle_fit_params (astropy.table.Table):
            Table containing the fit parameters.
        ech_angle_fit_coeffs (numpy.ndarray):
            Array containing the fit coefficients.
    """

    order_min, order_max = arxiv_params['order_min'], arxiv_params['order_max']
    order_vec = np.arange(order_min, order_max + 1, 1)
    norders = arxiv_params['norders'] # Total number of orders in the arxiv
    n_final = arxiv_params['n_final'] # order of wavelength solution fits
    ech_angles = arxiv["ech_angle"][arxiv["populated_and_good"]]
    # Determine the min,max params for the fits using all the echelle angles in the arxiv
    ech_min, ech_max = ech_angles.min(), ech_angles.max()
    ech_vec = ech_min + (ech_max - ech_min) * np.arange(100) / 99

    # Assign orders for each coefficient that we are fitting
    if nmax > n_final + 1:
        msgs.error(f'nmax={nmax} cannot be greater than n_final+1={n_final+1}. Reduce nmax')
    # This vector holds the polynomial order used to fit each coefficient
    coeff_fit_order_vec = np.full(n_final+1, coeff_fit_order_min)
    # DP: the fits look better if we remove the following line
    # coeff_fit_order_vec[0:nmax] = coeff_fit_order_max
    ## TODO This needs to be modified to be lower order for cases where there are very few fits in the arxiv. Right now
    # we fit the first 0:nmax coeffs always with a coeff_fit_order_max orderp polynomial.

    ech_angle_fit_params=Table([
        [ech_min],[ech_max],[norders], [order_min], [order_max], [n_final],[coeff_fit_order_vec], [func],[arxiv_params['func']],
        [arxiv_params['xmin']], [arxiv_params['xmax']]],
        names=('ech_xmin','ech_xmax','norders','order_min', 'order_max', 'ech_n_final','ech_coeff_fit_order', 'ech_func',
               'wave_func',
               'wave_xmin', 'wave_xmax'))

    ech_angle_fit_coeffs = np.zeros((norders, n_final + 1, coeff_fit_order_max + 1))

    for iord, this_order in enumerate(order_vec):
        populated = arxiv['populated_and_good'][:, iord] # & (arxiv['xdisp'][:, iord] == xdisp)
        nsolns_this_order = np.sum(populated)
        if nsolns_this_order > 0:
            ech_angle_this_order = arxiv['ech_angle'][:, iord][populated]
            coeff_this_order = arxiv['coeff'][:, iord, :][populated, :]
            #xd_angle_this_order = arxiv['xd_angle'][:, iord][populated]
            #lambda_cen_this_order = arxiv['lambda_cen'][:, iord][populated]
            for ic in range(coeff_this_order.shape[1]):
                pypeitFit = fitting.robust_fit(ech_angle_this_order, coeff_this_order[:, ic], coeff_fit_order_vec[ic], function=func,
                                               minx=ech_min, maxx=ech_max, maxiter=25,
                                               lower=sigrej, upper=sigrej, maxrej=maxrej, sticky=True, use_mad=True,
                                               weights=None)
                ech_angle_fit_coeffs[iord, ic, 0:coeff_fit_order_vec[ic]+1] = pypeitFit.fitc
                if debug:
                    this_fit = fitting.evaluate_fit(pypeitFit.fitc, func, ech_vec, minx=ech_min, maxx=ech_max)
                    plt.plot(ech_vec, this_fit, color='blue', label='fit')
                    fit_gpm = pypeitFit.bool_gpm
                    plt.plot(ech_angle_this_order[fit_gpm], coeff_this_order[fit_gpm, ic], marker='o', markersize=7.0, mfc='black',
                             mec='black', fillstyle='full', linestyle='None', zorder=5, label='used by fit')
                    plt.plot(ech_angle_this_order[np.logical_not(fit_gpm)], coeff_this_order[np.logical_not(fit_gpm), ic],
                             marker='s', markersize=9.0, mfc='red', mec='red', fillstyle='full', linestyle='None',
                             zorder=7, label='rejected')
                    plt.legend()
                    plt.title(
                        f'order={this_order}, cc_ii={ic}, nkept={np.sum(fit_gpm)}, nrej={np.sum(np.logical_not(fit_gpm))}')
                    plt.xlabel('ech_angle')
                    plt.ylabel('coeff')
                    plt.ylim(this_fit.min() - 0.05 * np.abs(this_fit.min()),
                             this_fit.max() + 0.05 * np.abs(this_fit.max()))
                    plt.show()


    return ech_angle_fit_params, ech_angle_fit_coeffs


def fit_reddest_vs_xd_angle(arxiv, func='legendre', polyorder = 2, sigrej=3.0, maxrej=1, debug=False):

    # Fit the reddest order on detector 3. We use this as the reference since the bluest order seems
    # to have a more noisy trend
    xd_angles = arxiv['xd_angle_file']
    # Determine the min,max params for the fits using all the XD angles in the arxiv
    xd_min, xd_max = xd_angles.min(), xd_angles.max()
    xd_vec = xd_min + (xd_max - xd_min) * np.arange(100) / 99

    xd_angle_fit_params=Table([[xd_min],[xd_max], [['UV', 'RED']], [polyorder], [func]]
                ,names=('xd_xmin','xd_xmax','xdisp_vec', 'xd_polyorder', 'xd_func'))

    # First dimension is UV or RED, second dimension is the set of polynomial coefficients
    xd_angle_fit_coeffs = np.zeros((2, polyorder + 1))

    for idisp, xdisp in enumerate(['UV', 'RED']):

        indx = (arxiv['det_file'] == 3) & (arxiv['xdisp_file'] == xdisp)
        xd_angles_this_disp = xd_angles[indx]
        reddest_order_this_disp = arxiv['reddest_order'][indx].astype(float)
        pypeitFit = fitting.robust_fit(xd_angles_this_disp, reddest_order_this_disp, polyorder,
                                       function=func, minx=xd_min, maxx=xd_max, maxiter=25,
                                       lower=sigrej, upper=sigrej, maxrej=maxrej, sticky=True, use_mad=True,
                                       weights=None)
        xd_angle_fit_coeffs[idisp, :] = pypeitFit.fitc
        if debug:
            this_fit = fitting.evaluate_fit(pypeitFit.fitc, func, xd_vec, minx=xd_min, maxx=xd_max)
            plt.plot(xd_vec, this_fit, color='green', label='fit')
            fit_gpm = pypeitFit.bool_gpm
            plt.plot(xd_angles_this_disp[fit_gpm], reddest_order_this_disp[fit_gpm], marker='o', markersize=7.0,
                     mfc='black', mec='black', fillstyle='full', linestyle='None', zorder=5, label='used by fit')
            plt.plot(xd_angles_this_disp[np.logical_not(fit_gpm)],reddest_order_this_disp[np.logical_not(fit_gpm)],
                     marker='s', markersize=9.0, mfc='red', mec='red', fillstyle='full', linestyle='None',
                     zorder=7, label='rejected')
            plt.legend()
            plt.title(f'XDISP={xdisp}, nkept={np.sum(fit_gpm)}, nrej={np.sum(np.logical_not(fit_gpm))}')
            plt.xlabel('xd_angle')
            plt.ylabel('reddest_order')
            plt.ylim(this_fit.min()-3, this_fit.max() +3)
            plt.show()

    return xd_angle_fit_params, xd_angle_fit_coeffs





def echelle_composite_arcspec(arxiv_file, outfile, show_individual_solns=False, do_total=False, show_orders=False, debug=False):

    color_tuple = ('green', 'cyan', 'magenta', 'blue', 'darkorange', 'yellow', 'dodgerblue', 'purple',
                   'lightgreen', 'cornflowerblue')
    colors = itertools.cycle(color_tuple)

    arxiv_params = Table.read(arxiv_file, hdu=1)[0]
    arxiv = Table.read(arxiv_file, hdu=2)
    norders = arxiv_params['norders']
    order_vec = np.arange(arxiv_params['order_min'], arxiv_params['order_max'] + 1, 1)


    # First loop over the orders to determine wavelength coverage and sampling of each order
    wave_grid_min = np.zeros(norders)
    wave_grid_max = np.zeros(norders)
    dwave_pix = np.zeros(norders)
    dloglam_pix = np.zeros(norders)
    nspec_per_order = np.zeros(norders, dtype=int)
    for iord, this_order in enumerate(order_vec):
        # TODO Should we make XDISP specific composites? I don't think so
        populated = arxiv['populated_and_good'][:, iord]
        nsolns_this_order = np.sum(populated)
        if nsolns_this_order > 0:
            this_wave = arxiv['wave'][:, iord, :][populated, :]
            this_arc = arxiv['arcspec'][:, iord, :][populated, :]
            this_gpm = (this_wave > 0.0) & (this_arc != 0.0)
            this_dwave = np.zeros_like(this_wave)
            for ii in range(nsolns_this_order):
                this_dwave[ii, :] = wvutils.get_delta_wave(this_wave[ii, :], this_gpm[ii, :])

            wave_grid_min[iord], wave_grid_max[iord] = this_wave.min(), this_wave.max()
            dwave_pix[iord] = np.median(this_dwave.min(axis=1, where=this_dwave!=0, initial=10))
            dloglam_pix[iord] = np.median((this_dwave/this_wave/np.log(10.0)).min(axis=1, where=this_dwave!=0, initial=10))
        else:
            msgs.error(f'No arc solutions contribute to order={iord}. There must be a bug')

    # Use the smallest value of dloglam across all orders for the spectral grid spacing
    dloglam_pix_final = dloglam_pix.min()
    dv_pix_final = np.log(10.0)*c_kms*dloglam_pix_final

    # Use the same wavelength grid for all orders
    wave_total_composite, wave_total_composite_mid, dsamp = wvutils.wavegrid(
        wave_grid_min.min(), wave_grid_max.max(), dloglam_pix_final, log10=True)

    nspec_composite = wave_total_composite.size
    ind_min =  np.zeros(norders, dtype=int)
    ind_max =  np.zeros(norders, dtype=int)
    nvec = np.arange(nspec_composite)

    # Determine the size of the output array by looping over all orders and finding the maximum grid we need to store
    for iord, this_order in enumerate(order_vec):
        indx = (wave_total_composite >= wave_grid_min[iord]) & (wave_total_composite <= wave_grid_max[iord])
        nspec_per_order[iord] = np.sum((wave_total_composite >= wave_grid_min[iord]) & (wave_total_composite <= wave_grid_max[iord]))
        ind_min[iord] = nvec[indx].min()
        ind_max[iord] = nvec[indx].max()

    # Allocate output arrays for composite arc
    nspec_max = nspec_per_order.max()
    wave_composite = np.zeros((nspec_max, norders))
    arc_composite = np.zeros((nspec_max, norders))
    gpm_composite = np.zeros((nspec_max, norders), dtype=bool)

    sn_smooth_npix = 1  # Should not matter since we use uniform weights
    for iord, this_order in enumerate(order_vec):
        populated = arxiv['populated_and_good'][:, iord]
        nsolns_this_order = np.sum(populated)
        if nsolns_this_order > 0:
            this_wave = arxiv['wave'][:, iord, :][populated, :]
            this_arc = arxiv['arcspec'][:, iord, :][populated, :]
            this_gpm = (this_wave > 0.0) & (this_arc != 0.0)
            this_wave_composite = wave_total_composite[ind_min[iord]:ind_max[iord]]
            this_nspec = this_wave_composite.size
            arc_interp_iord = np.zeros((this_nspec, nsolns_this_order))
            gpm_arc_iord = np.zeros((this_nspec, nsolns_this_order), dtype=bool)
            # Interpolate our arcs onto the new grid
            for ii in range(nsolns_this_order):
                this_arc_interp = interpolate.interp1d(
                    this_wave[ii, this_gpm[ii, :]], this_arc[ii, this_gpm[ii,:]], kind='cubic',
                    bounds_error=False, fill_value=-1e10)(this_wave_composite)
                arc_interp_iord[:, ii] = this_arc_interp
                gpm_arc_iord[:, ii] = arc_interp_iord[:, ii] > -1e9
                if show_individual_solns:
                    # plt.plot(iwave[in_gpm], this_arc[ii, in_gpm], color=next(colors), alpha=0.7)
                    plt.plot(this_wave_composite[gpm_arc_iord[ii, :]],this_arc_interp[gpm_arc_iord[ii, :]],
                             color=next(colors), alpha=0.7)
            if show_individual_solns:
                plt.title(f'Order={this_order}')
                plt.show()

            wave_grid_in = np.repeat(this_wave_composite[:, np.newaxis], nsolns_this_order, axis=1)
            ivar_arc_iord = utils.inverse(np.abs(arc_interp_iord) + 10.0)
            wave_grid_mid, wave_grid_stack, arcspec_stack, _, arcspec_gpm, = coadd.combspec(
                utils.array_to_explist(wave_grid_in), utils.array_to_explist(arc_interp_iord),
                utils.array_to_explist(ivar_arc_iord), utils.array_to_explist(gpm_arc_iord), sn_smooth_npix,
                wave_method='user_input', wave_grid_input=this_wave_composite,
                ref_percentile=70.0, maxiter_scale=5, sigrej_scale=3.0, scale_method='median',
                sn_min_polyscale=2.0, sn_min_medscale=0.5, const_weights=True, maxiter_reject=5, sn_clip=30.0,
                lower=5.0, upper=5.0,
                debug=debug, debug_scale=debug, show_scale=debug, show=show_orders, verbose=True)
            #ind_mid_min[iord] = np.argmin(np.abs(wave_grid_mid.min() - wave_total_composite_mid))
            #ind_mid_max[iord] = np.argmin(np.abs(wave_grid_mid.max() - wave_total_composite_mid))
            wave_composite[0:wave_grid_mid.size, iord] = wave_grid_mid
            arc_composite[0:wave_grid_mid.size, iord] = arcspec_stack
            gpm_composite[0:wave_grid_mid.size, iord] = arcspec_gpm

    # Now generate a final composite arc combining all the orders. Experimental. Not sure we need this.
    if do_total:
        show_total=False
        ivar_composite = utils.inverse(np.abs(arc_composite) + 10.0)
        # TODO this will crash since it is not taking lists.
        wave_grid_mid, wave_grid_stack, arcspec_stack, _, arcspec_gpm = coadd.combspec(
            wave_composite, arc_composite, ivar_composite, gpm_composite, sn_smooth_npix,
            wave_method='user_input', wave_grid_input=wave_total_composite, ref_percentile=70.0, maxiter_scale=5, sigrej_scale=3.0, scale_method='median',
            sn_min_polyscale=2.0, sn_min_medscale=0.5, const_weights=True, maxiter_reject=5, sn_clip=30.0,
            lower=5.0, upper=5.0,
            debug=debug, debug_scale=debug, show_scale=debug, show=show_total, verbose=True)

    params = Table([[os.path.basename(arxiv_file)], [arxiv_params['order_min']],[arxiv_params['order_max']],[norders],
                  [wave_composite[gpm_composite > 0.0].min()], [wave_composite[gpm_composite > 0.0].max()],
                  [dloglam_pix_final], [dv_pix_final]],
                  names=('arxiv_file','order_min', 'order_max', 'norders','wave_min','wave_max','dloglam','dv'))

    hdulist = fits.HDUList()
    hdulist.append(fits.BinTableHDU(params.as_array()))  # hdu = 1
    hdulist.append(fits.ImageHDU(np.array(wave_composite)))  # hdu = 2
    hdulist.append(fits.ImageHDU(np.array(arc_composite)))  # hdu = 3
    hdulist.append(fits.ImageHDU(np.array(gpm_composite.astype(float))))  # hdu = 3
    hdulist.writeto(outfile, overwrite=True)


if __init__ == '__main__':
    # xidl_arxiv_file = os.path.join(os.getenv('PYPEIT_DEV'), 'dev_algorithms', 'hires_wvcalib', 'hires_wvcalib_xidl.fits')
    # # Create the astropy table form of the xidl save file arxiv
    # if not os.path.isfile(xidl_arxiv_file):
    #     ingest_xidl_archive(xidl_arxiv_file)
    # # append the pypeit templates to the xidl archive
    # arxiv_file = os.path.join(os.getenv('PYPEIT_DEV'), 'dev_algorithms', 'hires_wvcalib', 'hires_wvcalib.fits')
    # if not os.path.isfile(arxiv_file):
    #     append_pypeit_archive(arxiv_file, xidl_arxiv_file)

    # Load the wavelength solutions


    # sys.exit(-1)
    # Perform fits to the coefficients vs ech angle
    # TODO see if pca works better here
    debug=False
    wvcalib_angle_fit_file = os.path.join(os.getenv('PYPEIT_DEV'), 'dev_algorithms', 'wavelengths', 'vlt_uves', 'wvcalib_angle_fits.fits')
    if not os.path.isfile(wvcalib_angle_fit_file):
        fit_wvcalib_vs_angles(arxiv_file, wvcalib_angle_fit_file, func='legendre',
                          ech_nmax = 3, ech_coeff_fit_order_min=1, ech_coeff_fit_order_max=2,
                          xd_reddest_fit_polyorder=2, sigrej=3.0, maxrej=1, debug=debug)

    # Compute a composite arc from the solution arxiv
    composite_arcfile = os.path.join(os.getenv('PYPEIT_DEV'), 'dev_algorithms', 'wavelengths', 'vlt_uves', 'UVES_composite_arc.fits')
    if not os.path.isfile(composite_arcfile):
        echelle_composite_arcspec(arxiv_file, composite_arcfile, show_orders=debug)

    sys.exit(-1)

    use_unknowns = True
    line_lists_all = waveio.load_line_lists(['ThAr'])
    line_lists = line_lists_all[np.where(line_lists_all['ion'] != 'UNKNWN')]
    unknwns = line_lists_all[np.where(line_lists_all['ion'] == 'UNKNWN')]
    tot_line_list = table.vstack([line_lists, unknwns]) if use_unknowns else line_lists
    spectrograph = load_spectrograph('keck_hires')
    par = spectrograph.default_pypeit_par()['calibrations']['wavelengths']
    n_final = 4
    # xmin, xmax for wavelength vs pixel fits
    fmin, fmax = 0.0, 1.0
    color_tuple = ('green', 'cyan', 'magenta', 'blue', 'darkorange', 'yellow', 'dodgerblue', 'purple',
                   'lightgreen', 'cornflowerblue')
    colors = itertools.cycle(color_tuple)
    #
    #
    # # Read template file
    # templ_table_file = os.path.join(
    #     resource_filename('pypeit', 'data'), 'arc_lines',
    #     'hires', 'hires_templ_xidl.dat')
    # tbl = Table.read(templ_table_file, format='ascii')
    # nrows = len(tbl)
    #
    # order_min = tbl['IOrder'].min()
    # order_max = 118
    #
    # order_vec = np.arange(order_min, order_max +1, 1)
    # norders = order_vec.size
    #
    # # Subset of orders in every file. Populated indicates whether a given order is populated
    # lambda_cen = np.zeros((norders, nrows))
    # ech_angle = np.zeros((norders, nrows))
    # populated = np.zeros((norders, nrows), dtype=bool)
    # XDISP_is_red = np.zeros((norders, nrows), dtype=bool)
    # binspec = np.zeros((norders, nrows), dtype=int)
    # det = np.zeros((norders, nrows), dtype=int)
    # xd_angle = np.zeros((norders, nrows))
    # coeff = np.zeros((norders, nrows, n_final+1))
    # bluest_order = np.zeros(nrows, dtype=int)
    # xd_angle_file = np.zeros(nrows)
    # ech_angle_file = np.zeros(nrows)
    # det_file = np.zeros(nrows)
    # XDISP_is_red_file = np.zeros(nrows, dtype=bool)
    #
    # for irow in np.arange(nrows):
    #     this_order_vec_raw, this_wave, this_arc = templates.xidl_hires(
    #         os.path.join(os.getenv('HIRES_CALIBS'), 'ARCS', tbl[irow]['Name']), specbin=tbl[irow]['Rbin'])
    #     if irow == 0:
    #         nspec = this_wave.shape[1]
    #         xnspecmin1 = float(nspec - 1)
    #         xvec = np.arange(nspec)/xnspecmin1
    #         wave = np.zeros((norders, nrows, nspec))
    #         arcspec = np.zeros((norders, nrows, nspec))
    #         table_xidl = Table([np.zeros((nrows, norders), dtype="<U30"),
    #                             np.zeros((nrows, norders)),
    #                             np.zeros((nrows, norders)),
    #                             np.zeros((nrows, norders), dtype=int),
    #                             np.zeros((nrows,), dtype=int),
    #                             np.zeros((nrows, norders), dtype="<U3"),
    #                             np.zeros((nrows, norders), dtype=int),
    #                             np.zeros((nrows, norders), dtype=int),
    #                             np.zeros((nrows, norders), dtype=bool),
    #                             np.zeros((nrows, norders, n_final+1)),
    #                             np.zeros((nrows, norders, nspec)),
    #                             np.zeros((nrows, norders, nspec)),],
    #                             names = ('filename', 'ech_angle', 'xd_angle', 'order', 'bluest_order',
    #                                      'xdisp', 'det', 'binspec', 'populated', 'coeff', 'wave', 'arcspec'))
    #     else:
    #         assert this_wave.shape[1] == nspec
    #     # Restrict to what is labeled as good in the Table
    #     igood = (this_order_vec_raw >= tbl[irow]['IOrder']) & (this_order_vec_raw <= tbl[irow]['EOrder'])
    #     nsolns = np.sum(igood)
    #     this_order_vec = this_order_vec_raw[igood]
    #     indx = this_order_vec - order_min
    #     populated[indx, irow] = True
    #     ech_angle[indx, irow] = tbl[irow]['ECH']
    #     xd_angle[indx, irow] = tbl[irow]['XDAng']
    #     XDISP_is_red[indx, irow] = tbl[irow]['XDISP'] == 'RED'
    #     binspec[indx, irow] =  tbl[irow]['Rbin']
    #     det[indx, irow] =  tbl[irow]['Chip']
    #
    #     wave[indx, irow, :] = this_wave[igood, :]
    #     arcspec[indx, irow, :] = this_arc[igood, :]
    #     lambda_cen[indx, irow] = np.median(this_wave[igood, :], axis=1)
    #     # Fit the wavelengths
    #     coeff_array = np.zeros((nsolns, n_final +1))
    #     for ii, iwave in enumerate(this_wave[igood, :]):
    #         pypeitFit = fitting.robust_fit(xvec, iwave, n_final, function=par['func'], maxiter=10,
    #                                        lower=1e10, upper=1e10, maxrej=0, sticky=True,
    #                                        minx=fmin, maxx=fmax, weights=None)
    #         coeff_array[ii, :] = pypeitFit.fitc
    #     coeff[indx, irow, :] = coeff_array
    #     # file specific
    #     bluest_order[irow] = this_order_vec[-1]
    #     ech_angle_file[irow] = tbl[irow]['ECH']
    #     xd_angle_file[irow] = tbl[irow]['XDAng']
    #     det_file[irow] = tbl[irow]['Chip']
    #     XDISP_is_red_file[irow] = tbl[irow]['XDISP'] == 'RED'

    #all_dlam = []
    #all_lam = []
    #all_orders = []


    pad_factor = 0.10
    # xvec_pad = np.arange(-int(np.round(pad_factor*nspec)), int(np.round((1.0 + pad_factor)*nspec)))/xnspecmin1
    #
    # # Plot the polynomial coefficients versus echelle angle order by order
    # debug_all=False
    # show_wv_grid=False
    # ncoeff_fit_order = 2
    # coeff_vs_order = np.zeros((norders, n_final + 1, ncoeff_fit_order+1))
    # ech_min, ech_max = ech_angle_file.min(), ech_angle_file.max()
    # ech_vec = ech_min + (ech_max-ech_min)*np.arange(100)/99
    # func = 'legendre'
    # debug_fits=False
    # for iord, this_order in enumerate(order_vec):
    #     if np.any(populated[iord, :]):
    #         nsolns = np.sum( populated[iord, :])
    #         this_ech = ech_angle[iord, populated[iord, :]]
    #         this_xd_angle = xd_angle[iord, populated[iord, :]]
    #         this_lambda_cen = lambda_cen[iord, populated[iord, :]]
    #         this_coeff = coeff[iord, populated[iord, :], :]
    #         for ic in range(n_final + 1):
    #             pypeitFit = fitting.robust_fit(this_ech, this_coeff[:, ic], ncoeff_fit_order, function=func,
    #                                            minx=ech_min, maxx=ech_max, maxiter=25,
    #                                            lower=3.0, upper=3.0, maxrej=2, sticky=True,use_mad=True, weights=None)
    #             coeff_vs_order[iord, ic, :] = pypeitFit.fitc
    #             if debug_fits:
    #                 this_fit = fitting.evaluate_fit(pypeitFit.fitc, func, ech_vec, minx=ech_min, maxx=ech_max)
    #                 plt.plot(ech_vec, this_fit, color='blue', label='fit')
    #                 fit_gpm = pypeitFit.bool_gpm
    #                 plt.plot(this_ech[fit_gpm], this_coeff[fit_gpm, ic], marker='o', markersize=7.0, mfc='black',
    #                          mec='black', fillstyle='full', linestyle='None', zorder=5, label='used by fit')
    #                 plt.plot(this_ech[np.logical_not(fit_gpm)], this_coeff[np.logical_not(fit_gpm), ic],
    #                          marker='s', markersize=9.0, mfc='red', mec='red',fillstyle='full', linestyle='None', zorder=7,label='rejected')
    #                 plt.legend()
    #                 plt.title(f'order={this_order}, cc_ii={ic}, nkept={np.sum(fit_gpm)}, nrej={np.sum(np.logical_not(fit_gpm))}')
    #                 plt.xlabel('ech_angle')
    #                 plt.ylabel('coeff')
    #                 plt.ylim(this_fit.min() - 0.05*np.abs(this_fit.min()), this_fit.max() + 0.05*np.abs(this_fit.max()))
    #                 plt.show()
    #         #for ii in range(n_final+1):
    #         #    plt.plot(this_xd_angle, this_coeff[:, ii], 'k.', label=f'order={iorder}, cc_ii={ii}')
    #         #    plt.legend()
    #         #    plt.xlabel('xd_angle')
    #         #    plt.ylabel('coeff')
    #         #    plt.show()


    # Plot lam vs dlam/lam for each order

    #for iord, this_order in enumerate(order_vec):
    for iord in np.arange(norders)[::-1]:
        this_order = order_vec[iord]
        if np.any(populated[iord, :]):
            nsolns = np.sum(populated[iord, :])
            this_ech = ech_angle[iord, populated[iord, :]]
            this_xd_angle = xd_angle[iord, populated[iord, :]]
            this_lambda_cen = lambda_cen[iord, populated[iord, :]]
            this_wave = wave[iord, populated[iord, :], :]
            this_arc = arcspec[iord, populated[iord, :], :]
            #this_coeff = coeff[iord, populated[iord, :], :]

            this_dwave = np.zeros_like(this_wave)
            for ii, iwave in enumerate(this_wave):
                this_dwave[ii, :] = wvutils.get_delta_wave(iwave, (iwave > 0.0))

            # Now try a fit. TODO any wavelength grid will work here
            med_dlam = np.median(this_dwave[this_wave > 1.0])
            fit = robust_fit(this_wave.flatten(), this_dwave.flatten(), 3, maxiter=25, maxdev = 0.10*med_dlam, groupbadpix=True)
            wave_grid_fit, wave_grid_fit_mid, dsamp = wvutils.get_wave_grid(this_wave.T,wave_method='log10')
            dwave_fit = fit.eval(wave_grid_fit)
            gpm = fit.bool_gpm.copy()
            gpm.resize(this_wave.shape)

            if show_wv_grid:
                for ii, iwave in enumerate(this_wave):
                    this_color=next(colors)
                    this_gpm = gpm[ii, :]
                    plt.plot(iwave[this_gpm], this_dwave[ii, this_gpm], marker='o', markersize=1.0, mfc=this_color,
                    fillstyle='full',  linestyle='None', zorder=1)
                    plt.plot(iwave[np.logical_not(this_gpm)], this_dwave[ii, np.logical_not(this_gpm)], marker='o',
                         markersize=2.0, mfc='red', fillstyle='full', zorder=3, linestyle='None')

                plt.plot(wave_grid_fit, dwave_fit, color='black', label='fit', zorder=10)
                plt.title(f'order={this_order}', fontsize=14)
                plt.legend()
                plt.show()

            lam_min, lam_max = this_wave[gpm].min(), this_wave[gpm].max()
            wave_grid = get_variable_dlam_wavegrid(lam_min, lam_max, wave_grid_fit, dwave_fit)
            nspec_tmpl = wave_grid.shape[0]
            # TESTING
            #dwave_chk = wvutils.get_delta_wave(wave_grid, (wave_grid > 0.0))
            #plt.plot(wave_grid, dwave_chk, color='red', label='our new grid')
            #plt.plot(wave_grid_fit, dwave_fit, color='black', label='fit')
            #plt.title('dwave compared to our fit')
            #plt.legend()
            #plt.show()
            tmpl_iord = np.zeros((nsolns, nspec_tmpl))
            gpm_tmpl = np.zeros((nsolns, nspec_tmpl), dtype=bool)
            # Interpolate our arcs onto the new grid
            for ii, iwave in enumerate(this_wave):
                in_gpm = this_arc[ii, :] != 0.0
                tmpl_iord[ii, :] = interpolate.interp1d(iwave[in_gpm], this_arc[ii, in_gpm], kind='cubic', bounds_error=False, fill_value=-1e10)(wave_grid)
                gpm_tmpl[ii, :] = tmpl_iord[ii, :] > -1e9
                if show_wv_grid:
                    #plt.plot(iwave[in_gpm], this_arc[ii, in_gpm], color=next(colors), alpha=0.7)
                    plt.plot(wave_grid[gpm_tmpl[ii, :]], tmpl_iord[ii, gpm_tmpl[ii, :]], color=next(colors), alpha=0.7)
                    plt.show()

            sn_smooth_npix = 1 # Should not matter since we use uniform weights
            wave_grid_in = np.repeat(wave_grid[:, np.newaxis], nsolns, axis=1)
            ivar_tmpl_iord = utils.inverse(np.abs(tmpl_iord) + 10.0)
            wave_grid_mid, wave_grid_stack, arcspec_tmpl, _, arcspec_tmpl_gpm = coadd.combspec(
                wave_grid_in, tmpl_iord.T, ivar_tmpl_iord.T, gpm_tmpl.T, sn_smooth_npix,
                wave_method='iref',  ref_percentile=70.0, maxiter_scale=5, sigrej_scale=3.0, scale_method='median',
                sn_min_polyscale=2.0, sn_min_medscale=0.5, const_weights=True, maxiter_reject=5, sn_clip=30.0, lower=5.0, upper=5.0,
                debug=debug_all, debug_scale=debug_all, show_scale=debug_all, show=True, verbose=True)

            all_patt_dict = {}
            detections = {}
            wv_calib = {}

            all_patt_dict_pad = {}
            detections_pad = {}
            wv_calib_pad = {}

            for slit in range(nsolns):
                print('Working on soln={:d}'.format(slit))
                # Trim the template to the relevant range. Hack this for now
                #itmpl = (wave_grid_mid >= 0.999*iwave.min()) & (wave_grid_mid <= 1.001*iwave.max())
                coeff_predict = np.zeros(n_final + 1)
                for ic in range(n_final+1):
                    coeff_predict[ic] = fitting.evaluate_fit(coeff_vs_order[iord, ic, :], func, this_ech[slit], minx=ech_min, maxx=ech_max)

                wave_predict = fitting.evaluate_fit(coeff_predict, par['func'], xvec, minx=fmin, maxx=fmax)
                wave_predict_pad = fitting.evaluate_fit(coeff_predict, par['func'], xvec_pad, minx=fmin, maxx=fmax)
                wave_true = this_wave[slit, :]
                # Substitute wave_true here as a test
                arcspec_templ_predict_pad =  interpolate.interp1d(wave_grid_stack[arcspec_tmpl_gpm], arcspec_tmpl[arcspec_tmpl_gpm],
                                                              kind='cubic', bounds_error=False, fill_value=0.0)(wave_predict_pad)
                arcspec_templ_predict =  interpolate.interp1d(wave_grid_stack[arcspec_tmpl_gpm], arcspec_tmpl[arcspec_tmpl_gpm],
                                                              kind='cubic', bounds_error=False, fill_value=0.0)(wave_predict)
                #arcspec_tmpl_trim = arcspec_tmpl[itmpl]
                #wave_grid_mid_trim = wave_grid_mid[itmpl]
                #arc_in_pad = np.zeros_like(arcspec_tmpl_trim)
                #in_gpm = this_arc[slit, :] != 0.0
                #npix = np.sum(in_gpm)
                #arc_in_pad[:npix] = this_arc[slit, in_gpm]
                #xcorr_poly(this_wave[slit, in_gpm], this_arc[slit, in_gpm], wave_grid_mid, arcspec_tmpl, smooth=1.0, percent_ceil=50.0, use_raw_arc=False,
                #           sigdetect=10.0, fwhm=4.0, debug=True, seed=42)

                # WITH PADDING
                #arc_pad = np.zeros_like(xvec_pad)
                #arc_pad[0:nspec] = this_arc[slit, :]
                #detections_pad[str(slit)], spec_cont_sub_pad, all_patt_dict_pad[str(slit)] = autoid.reidentify(
                #    arc_pad, arcspec_templ_predict_pad, wave_predict_pad,  tot_line_list, par['nreid_min'],
                #    cc_thresh=par['cc_thresh'], match_toler=par['match_toler'], cc_local_thresh=par['cc_local_thresh'],
                #    nlocal_cc=par['nlocal_cc'], nonlinear_counts=1e10,
                #    sigdetect=par['sigdetect'], fwhm=par['fwhm'], debug_peaks=True, debug_xcorr=True, debug_reid=True)

                # WITHOUT PADDING
                detections[str(slit)], spec_cont_sub, all_patt_dict[str(slit)] = autoid.reidentify(
                    this_arc[slit, :], arcspec_templ_predict, wave_predict,  tot_line_list, par['nreid_min'],
                    cc_thresh=par['cc_thresh'], match_toler=par['match_toler'], cc_local_thresh=par['cc_local_thresh'],
                    nlocal_cc=par['nlocal_cc'], nonlinear_counts=1e10,
                    sigdetect=par['sigdetect'], fwhm=par['fwhm'], debug_peaks=True, debug_xcorr=True, debug_reid=True)

                # Check if an acceptable reidentification solution was found
                if not all_patt_dict[str(slit)]['acceptable']:
                    wv_calib[str(slit)] = None
                    continue

                final_fit = wv_fitting.fit_slit(spec_cont_sub, all_patt_dict[str(slit)], detections[str(slit)],
                                                tot_line_list, match_toler=par['match_toler'],func=par['func'], n_first=par['n_first'],
                                                sigrej_first=par['sigrej_first'], n_final=n_final,sigrej_final=par['sigrej_final'])

                #autoid.arc_fit_qa(final_fit, title='Silt: {}'.format(str(slit)))


    #        dlam = []
    #        lam = []
    #        for iwave in this_wave:
    #            dlam += list(wvutils.get_delta_wave(iwave, np.ones_like(iwave,dtype=bool)))
    #            lam += iwave.tolist()
    #            #xlam += ((np.array(iwave) - np.array(iwave).min())/(np.array(iwave).max() - np.array(iwave).min())).tolist()
    #        plt.plot(lam, 3.0e5*np.array(dlam)/np.array(lam), '.', label=f'order={iorder}')
    #        plt.legend()
    #        plt.show()

     #       all_dlam += dlam
     #       all_lam += lam
     #       all_orders += [iorder]*len(lam)


    #plt.plot(all_lam, 3.0e5*np.array(all_dlam)/np.array(all_lam), '.')
    #plt.legend()
    #plt.show()

    # Plot the central wavelength vs echelle angle order by order
    for iord, iorder in enumerate(order_vec):
        if np.any(populated[iord, :]):
            this_ech = ech_angle[iord, populated[iord, :]]
            this_xd_angle = xd_angle[iord, populated[iord, :]]
            this_lambda_cen = lambda_cen[iord, populated[iord, :]]
            plt.plot(this_ech, this_lambda_cen, 'k.', label=f'order={iorder}')
            plt.legend()
            plt.show()


    for xdisp in ['UV', 'RED']:
        for idet in [1,2,3]:
            iord = (XDISP_is_red_file == (xdisp == 'RED')) & (det_file == idet)
            plt.plot(xd_angle_file[iord], bluest_order[iord], 'k.', label=f'XDISP={xdisp}, det={idet}')
            plt.legend()
            plt.show()


