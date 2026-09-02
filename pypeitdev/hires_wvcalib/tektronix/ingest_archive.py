import os
import itertools
from IPython import embed

import numpy as np
from scipy import interpolate
from astropy import constants as const
from astropy.table import Table
from astropy.io import fits
from astropy import table
from matplotlib import pyplot as plt

from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.wavecal import wvutils
from pypeit.core import arc
from pypeit.core.fitting import robust_fit
from pypeit.core import coadd
from pypeit.core import fitting
from pypeit.core.wave import airtovac
from pypeit.core.wavecal import autoid, waveio, wv_fitting
from pypeit.core.wavecal.wvutils import  get_xcorr_arc, xcorr_shift
from pypeit import utils
from pypeit import wavecalib
from pypeit import PypeItError

import pypeit
pypeit_path = os.path.dirname(os.path.realpath(pypeit.__file__))

c_kms = const.c.to('km/s').value


def load_archive(outfile, n_final=2, func='legendre'):
    """
    Load the pypeit format archived templates that were created with the `esorex_to_pypeit` function.

    Args:
        outfile (str):
            Name of the output file
    """
    # load the list of MAKEE templates
    ptempl_table_file = os.path.join(os.getenv('PYPEIT_DEV'), 'pypeitdev', 'hires_wvcalib', 'tektronix',
                                     'hires_orig_templ_makee.dat')
    ptbl = Table.read(ptempl_table_file, format='ascii.fixed_width')

    # Number of rows in the pypeit template table for this xdisp
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

    # Load the HIRES (new) wavelength solution... the MAKEE solutions seem wrong (and doesn't seem to be just an air->vacuum issue, but that may contribute)
    hires_file = os.path.join(pypeit_path, 'data', 'arc_lines', 'reid_arxiv',
                             'keck_hires_composite_arc.fits')
    print("Loading template file: ", hires_file)
    hires_new = fits.open(hires_file, chk_version=False)
    hires_orders = np.arange(hires_new[1].data['order_min'][0], hires_new[1].data['order_max'][0]+1)
    hires_waves = hires_new[2].data
    hires_specs = hires_new[3].data
    mskarr = np.ma.array(hires_waves, mask=(hires_waves == 0))
    mid_waves = np.ma.mean(mskarr, axis=0).data
    # Make a fitting function to get the order number as a function of wavelength
    order_fit = np.polyfit(mid_waves, hires_orders, 6)

    # Now deal with the MAKEE format
    # load the MAKEE templates
    for irow in np.arange(p_nrows):
        ifinal_row = irow
        templ_file = os.path.join(os.getenv('PYPEIT_DEV'), 'pypeitdev', 'hires_wvcalib', 'tektronix',
                                  'templates_makee', ptbl[irow]['Name'])
        print("Loading template file: ", templ_file)
        hdu = fits.open(templ_file, chk_version=False)
        # this is the order vector available for this template
        this_order_vec_raw = np.arange(ptbl[irow]['IOrder'], ptbl[irow]['EOrder'] + 1, 1)[::-1]
        # For each order, get the arc spectrum and wavelength solution
        nord = hdu[0].header['NAXIS2']
        if this_order_vec_raw.size != nord:
            embed()
            raise PypeItError(f'Number of orders in the header ({nord}) does not match the number of orders in the table ({this_order_vec_raw.size}) for file {templ_file}')
        all_wave, all_arc, all_shift = [], [], np.zeros(nord)
        # Repeat twice (get a better cross-correlation shift by iterating)
        for rep in range(2):
            for ii in range(nord):
                coeff = np.array(hdu[0].header['WV_0_{0:02d}'.format(ii+1)].split() + \
                         hdu[0].header['WV_4_{0:02d}'.format(ii+1)].split()).astype(float)[::-1]
                pixarr = np.arange(hdu[0].header['NAXIS1'])
                wavarr = np.polyval(coeff, pixarr)

                # Generate a spectrum from the hires template for this order and interpolate it onto the MAKEE wavelength solution
                this_ord = this_order_vec_raw[ii]
                ww_ord = np.where(hires_orders==this_ord)[0][0]
                this_hspec = hires_specs[:, ww_ord]
                this_hwave = hires_waves[:, ww_ord]
                wgud = (this_hwave > 0)
                this_hires_spec = np.interp(wavarr, this_hwave[wgud], this_hspec[wgud])
                # Plot to see if the MAKEE and HIRES spectra are aligned.
                # They are not, so we need to cross-correlate and shift the MAKEE spectrum to match the HIRES spectrum
                if False:
                    plt.plot(wavarr, hdu[0].data[ii], color='blue', label='MAKEE')
                    plt.plot(wavarr, this_hires_spec, color='red', label='HIRES')
                    plt.legend()
                    plt.show()
                if rep == 0:
                    # Calculate the shift
                    shift, corr = xcorr_shift(hdu[0].data[ii], this_hires_spec, 100, debug=False)
                    all_shift[ii] = 0.0#shift
                else:
                    # On the second iteration, apply the shift to the MAKEE wavelength solution and store it in the final table
                    wavarr = np.polyval(coeff, pixarr-fit_shift[ii])

                    all_wave.append(arc.resize_spec(wavarr, nspec))
                    # all_wave.append(arc.resize_spec(airtovac(wavarr*units.AA).value, nspec))
                    all_arc.append(arc.resize_spec(hdu[0].data[ii], nspec))
            if rep == 0:
                robfit = robust_fit(np.arange(nord), all_shift, 1, function='legendre', maxiter=10,
                                      lower=5, upper=5, maxrej=1, sticky=True,
                                        minx=0, maxx=nord-1, weights=None)
                fit_shift = robfit.eval(np.arange(nord))

        # Convert to numpy arrays
        this_wave = np.array(all_wave)
        this_arc = np.array(all_arc)

        # Assume all are good
        igood = np.ones(nord, dtype=bool)
        embed()

        this_order_vec = this_order_vec_raw[igood]
        nsolns = this_order_vec_raw.size
        nsolns_good = np.sum(igood)
        indx = this_order_vec_raw - order_min
        indx_good = this_order_vec - order_min
        # Information for the file is stored for convenience, although this is redundant with the arrays below
        final_table['filename'][ifinal_row] = ptbl[irow]['Name']
        final_table['nsolns'][ifinal_row] = nsolns
        final_table['nsolns_good'][ifinal_row] = nsolns_good
        final_table['bluest_order'][ifinal_row] = this_order_vec_raw[0]
        final_table['bluest_good_order'][ifinal_row] = this_order_vec[0]
        final_table['reddest_order'][ifinal_row] = this_order_vec_raw[-1]
        final_table['reddest_good_order'][ifinal_row] = this_order_vec[-1]
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
            pypeitFit = fitting.robust_fit(xvec, iwave, n_final, function=func, maxiter=10,
                                           lower=1e10, upper=1e10, maxrej=0, sticky=True,
                                           minx=xmin, maxx=xmax, weights=None)
            this_coeff_array[ii, :] = pypeitFit.fitc
        final_table['coeff'][ifinal_row, indx_good, :] = this_coeff_array

    # Write out to multi-extension fits
    print(f'Writing HIRES (original) PypeIt wv_calib archive to file: {outfile}')
    hdu_param = fits.BinTableHDU(params.as_array())
    hdu_table = fits.BinTableHDU(final_table.as_array())

    hdulist = fits.HDUList()
    hdulist.append(hdu_param)
    hdulist.append(hdu_table)
    hdulist.writeto(outfile, overwrite=True)


def empty_design_table(nrows, norders, n_final=4):
    """
    Construct an empty arxiv table.

    Args:
        nrows (:obj:`int`):
            Number of rows in the table.
        norders (:obj:`int`):
            Number of orders in the table.
        n_final (:obj:`int`, optional):
            Number of final coefficients in the wavelength solution.

    Returns:
        `astropy.table.Table`_: Instance of the empty arxiv table.
    """
    return Table([np.zeros(nrows, dtype="<U30"),
                        np.zeros(nrows, dtype=int),
                        np.zeros(nrows, dtype=int),
                        np.zeros(nrows, dtype=int),
                        np.zeros(nrows, dtype=int),
                        np.zeros(nrows, dtype=int),
                        np.zeros(nrows, dtype=int),
                        np.zeros(nrows, dtype="<U5"),
                        np.zeros(nrows, dtype=float),
                        np.zeros(nrows, dtype=float),
                        np.zeros(nrows, dtype=int),
                        np.zeros((nrows, norders), dtype=int),
                        np.zeros((nrows, norders), dtype=bool),
                        np.zeros((nrows, norders), dtype=bool),
                        np.full((nrows, norders), -1e10),
                        np.full((nrows, norders), -1e10),
                        np.zeros((nrows, norders), dtype="<U5"),
                        np.zeros((nrows, norders), dtype=int),
                        np.zeros((nrows, norders), dtype=int),
                        np.zeros((nrows, norders)),
                        np.zeros((nrows, norders, n_final + 1)),],
                        names=('filename', 'nsolns', 'nsolns_good', 'bluest_order', 'bluest_good_order',
                               'reddest_order', 'reddest_good_order', 'xdisp_file',
                               'ech_angle_file', 'xd_angle_file', 'det_file',
                               'order', 'populated', 'populated_and_good', 'ech_angle', 'xd_angle',
                               'xdisp', 'det', 'binspec','lambda_cen', 'coeff'))


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
    print(f"order_min = {order_min},\t order_max = {order_max}")
    print(f"order_vec = {order_vec}")
    print(f"norders = {norders}")
    print(f"n_final = {n_final}")
    print(f"ech_angle = {ech_angles}")

    # Determine the min,max params for the fits using all the echelle angles in the arxiv
    ech_min, ech_max = ech_angles.min(), ech_angles.max()
    ech_vec = ech_min + (ech_max - ech_min) * np.arange(100) / 99
    print(f"ech_min = {ech_min},\t ech_max = {ech_max}")
    print(f"ech_vec = {ech_vec}")

    # Assign orders for each coefficient that we are fitting
    if nmax > n_final + 1:
        raise PypeItError(f'nmax={nmax} cannot be greater than n_final+1={n_final+1}. Reduce nmax')
    # This vector holds the polynomial order used to fit each coefficient
    coeff_fit_order_vec = np.full(n_final+1, coeff_fit_order_min)
    print(f"coeff_fit_order_vec = {coeff_fit_order_vec}")
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
        print(f"iord = {iord}, this_order = {this_order}, nsolns_this_order = {nsolns_this_order}")
        if nsolns_this_order > 0:
            ech_angle_this_order = arxiv['ech_angle'][:, iord][populated]
            coeff_this_order = arxiv['coeff'][:, iord, :][populated, :]
            print(f"ech_angle_this_order = {ech_angle_this_order}, coeff_this_order = {coeff_this_order}")
            #xd_angle_this_order = arxiv['xd_angle'][:, iord][populated]
            #lambda_cen_this_order = arxiv['lambda_cen'][:, iord][populated]
            for ic in range(coeff_this_order.shape[1]):
                print(f"working for iord = {iord}, this_order = {this_order}, nsolns_this_order = {nsolns_this_order}, ic = {ic} =============================================")
                print(f"xarray = ech_angle_this_order = {ech_angle_this_order}")
                print(f"yarray = coeff_this_order = {coeff_this_order[:, ic]}")
                print(f"order = coeff_fit_order_vec[ic] = {coeff_fit_order_vec[ic]}")
                print(f"function = {func}")
                print(f"minx = {ech_min}")
                print(f"maxx = {ech_max}")
                print(f"lower = upper = {sigrej}")
                print(f"maxrej = {maxrej}")
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

    xd_angle_fit_params=Table([[xd_min],[xd_max], [['UV', 'RED', 'RED97']], [polyorder], [func]]
                ,names=('xd_xmin','xd_xmax','xdisp_vec', 'xd_polyorder', 'xd_func'))

    # First dimension is UV or RED, second dimension is the set of polynomial coefficients
    xd_angle_fit_coeffs = np.zeros((3, polyorder + 1))

    #for idisp, xdisp in enumerate(['UV', 'RED']):
    for idisp, xdisp in enumerate(['UV','RED','RED97']):


        print(f"xdisp = {xdisp}")
#        indx = (arxiv['det_file'] == 3) & (arxiv['xdisp_file'] == xdisp)
        indx = (arxiv['det_file'] == 1) & (arxiv['xdisp_file'] == xdisp)
        if len(np.where(indx)[0])==0:
            continue
        xd_angles_this_disp = xd_angles[indx]
        reddest_order_this_disp = arxiv['reddest_order'][indx].astype(float)
        print(f"{indx}")
        pypeitFit = fitting.robust_fit(xd_angles_this_disp, reddest_order_this_disp, polyorder,
                                       function=func, minx=xd_min, maxx=xd_max, maxiter=25,
                                       lower=sigrej, upper=sigrej, maxrej=maxrej, sticky=True, use_mad=True,
                                       weights=None)
        xd_angle_fit_coeffs[idisp, :] = pypeitFit.fitc
        debug = True
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
    print(f"order_vec = {order_vec}")


    # First loop over the orders to determine wavelength coverage and sampling of each order
    wave_grid_min = np.zeros(norders)
    wave_grid_max = np.zeros(norders)
    dwave_pix = np.zeros(norders)
    dloglam_pix = np.zeros(norders)
    nspec_per_order = np.zeros(norders, dtype=int)
    for iord, this_order in enumerate(order_vec):
        print(f"iord = {iord},  this_order = {this_order}")
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
            raise PypeItError(f'No arc solutions contribute to order={iord}. There must be a bug')
    print("end the loop of enumerate(order_vec)")

    # Use the smallest value of dloglam across all orders for the spectral grid spacing
    dloglam_pix_final = dloglam_pix.min()
    dv_pix_final = np.log(10.0)*c_kms*dloglam_pix_final
    print(f"dloglam_pix_final = {dloglam_pix_final}")
    print(f"dv_pix_final = {dv_pix_final}")
    print(f"wave_grid_min.min = {wave_grid_min.min()}")
    print(f"wave_grid_max.max = {wave_grid_max.max()}")


    # Use the same wavelength grid for all orders
    wave_total_composite, wave_total_composite_mid, dsamp = wvutils.wavegrid(
        wave_grid_min.min(), wave_grid_max.max(), dloglam_pix_final, log10=True)

    nspec_composite = wave_total_composite.size
    ind_min =  np.zeros(norders, dtype=int)
    ind_max =  np.zeros(norders, dtype=int)
    nvec = np.arange(nspec_composite)
    print(f"nspec_composite = {nspec_composite}")

    # Determine the size of the output array by looping over all orders and finding the maximum grid we need to store
    for iord, this_order in enumerate(order_vec):
        indx = (wave_total_composite >= wave_grid_min[iord]) & (wave_total_composite <= wave_grid_max[iord])
        nspec_per_order[iord] = np.sum((wave_total_composite >= wave_grid_min[iord]) & (wave_total_composite <= wave_grid_max[iord]))
        ind_min[iord] = nvec[indx].min()
        ind_max[iord] = nvec[indx].max()
        print(f"iord = {iord},  this_order = {this_order}, ind_min[iord] = {ind_min[iord]}, ind_max[iord] = {ind_max[iord]}")

    # Allocate output arrays for composite arc
    nspec_max = nspec_per_order.max()
    wave_composite = np.zeros((nspec_max, norders))
    arc_composite = np.zeros((nspec_max, norders))
    gpm_composite = np.zeros((nspec_max, norders), dtype=bool)

    sn_smooth_npix = 1  # Should not matter since we use uniform weights
    for iord, this_order in enumerate(order_vec):
        print(f"iord = {iord},  this_order = {this_order}")
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
            # wave_grid_mid, wave_grid_stack, arcspec_stack, _, arcspec_gpm, = coadd.combspec(
            #     utils.array_to_explist(wave_grid_in), utils.array_to_explist(arc_interp_iord),
            #     utils.array_to_explist(ivar_arc_iord), utils.array_to_explist(gpm_arc_iord), sn_smooth_npix,
            #     wave_method='user_input', wave_grid_input=this_wave_composite,
            #     ref_percentile=70.0, maxiter_scale=5, sigrej_scale=3.0, scale_method='median',
            #     sn_min_polyscale=2.0, sn_min_medscale=0.5, const_weights=True, maxiter_reject=5, sn_clip=30.0,
            #     lower=5.0, upper=5.0,
            #     debug=debug, debug_scale=debug, show_scale=debug, show=show_orders, verbose=True)
            wave_grid_mid, wave_grid_stack, arcspec_stack, _, arcspec_gpm, = coadd.combspec(
                utils.array_to_explist(wave_grid_in), utils.array_to_explist(arc_interp_iord),
                utils.array_to_explist(ivar_arc_iord), utils.array_to_explist(gpm_arc_iord), sn_smooth_npix,
                wave_method='user_input', wave_grid_input=this_wave_composite,
                ref_percentile=70.0, maxiter_scale=5, sigrej_scale=3.0, scale_method='median',
                sn_min_polyscale=2.0, sn_min_medscale=0.5, weight_method='constant', maxiter_reject=5, sn_clip=30.0,
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


    params=Table([[os.path.basename(arxiv_file)], [arxiv_params['order_min']],[arxiv_params['order_max']],[norders],
                  [wave_composite[gpm_composite > 0.0].min()], [wave_composite[gpm_composite > 0.0].max()],
                  [dloglam_pix_final], [dv_pix_final]],
                  names=('arxiv_file','order_min', 'order_max', 'norders','wave_min','wave_max','dloglam','dv'))


    print(f'Writing HIRES xidl+pypeit composite_arc archive to file: {outfile}')
    hdulist = fits.HDUList()
    hdulist.append(fits.BinTableHDU(params.as_array()))  # hdu = 1
    hdulist.append(fits.ImageHDU(np.array(wave_composite)))  # hdu = 2
    hdulist.append(fits.ImageHDU(np.array(arc_composite)))  # hdu = 3
    hdulist.append(fits.ImageHDU(np.array(gpm_composite.astype(float))))  # hdu = 3
    hdulist.writeto(outfile, overwrite=True)


if __name__ == '__main__':
    # Allow for command line arguments to specify an overwrite option
    import argparse
    parser = argparse.ArgumentParser(description='Ingest HIRES (original detector) makee wavelength calibrations and perform fits to the wavelength solution coefficients vs. XDISP and ECH angle')
    parser.add_argument('-o', '--overwrite', action='store_true', help='Whether to overwrite existing files')
    args = parser.parse_args()

    overwrite = args.overwrite

    # Store info in a file
    arxiv_file = os.path.join(os.getenv('PYPEIT_DEV'), 'pypeitdev', 'hires_wvcalib', 'tektronix', 'hires_wvcalib.fits')

    # Load the wavelength solutions
    load_archive(arxiv_file)

    # Perform fits to the coefficients vs ech angle
    # TODO see if pca works better here
    debug = False
    wvcalib_angle_fit_file = os.path.join(os.getenv('PYPEIT_DEV'), 'pypeitdev', 'hires_wvcalib', 'tektronix', f'keck_hires_orig_angle_fits.fits')
    if os.path.isfile(wvcalib_angle_fit_file) and not overwrite:
        print(f'File {wvcalib_angle_fit_file} already exists. Use --overwrite or -o to overwrite it.')
    else:
        fit_wvcalib_vs_angles(arxiv_file, wvcalib_angle_fit_file, func='legendre',
                          ech_nmax = 3, ech_coeff_fit_order_min=1, ech_coeff_fit_order_max=2,
                          xd_reddest_fit_polyorder=2, sigrej=3.0, maxrej=1, debug=debug)

    # Compute a composite arc from the solution arxiv
    print("A composite arc is not required for the wavelength calibration. We can just use the pre-existing one.")
    # composite_arcfile = os.path.join(os.getenv('PYPEIT_DEV'), 'pypeitdev', 'hires_wvcalib', 'tektronix', f'keck_hires_orig_composite_arc.fits')
    # if os.path.isfile(composite_arcfile) and not overwrite:
    #     print(f'File {composite_arcfile} already exists. Use --overwrite or -o to overwrite it.')
    # else:
    #     echelle_composite_arcspec(arxiv_file, composite_arcfile, show_orders=debug)

    print("All done!")
