""" Codes to help port ESI wavelengths from XIDL """
import os

import numpy as np

from astropy.table import Table

def ingest_xidl_archive(outfile, n_final=4, func='legendre'):

    # Read archive file
    tbl = Table.read('esi_archive_orders.fits')

    # Orders
    order_vec = np.arange(15, 6 - 1, -1)
    norders = order_vec.size

    # xmin, xmax for wavelength vs pixel fits
    fmin, fmax = 0.0, 1.0

    wv_dicts = []

    for ss, order in enumerate(order_vec):
        wv_dict = {}
        wv_dict['order'] = order

        #
        this_order_vec_raw, this_wave, this_arc = templates.xidl_hires(
            os.path.join(os.getenv('HIRES_CALIBS'), 'ARCS', tbl[irow]['Name']), specbin=tbl[irow]['Rbin'])
        if irow == 0:
            nspec = this_wave.shape[1]
            xnspecmin1 = float(nspec - 1)
            xvec = np.arange(nspec)/xnspecmin1
            params['nspec'] = nspec
            table_xidl['wave'] = np.zeros((nrows, norders, nspec))
            table_xidl['arcspec'] = np.zeros((nrows, norders, nspec))
        else:
            assert this_wave.shape[1] == nspec

        # Good is what is labeled as good in the Table, we also store everything in the file
        igood = (this_order_vec_raw >= tbl[irow]['IOrder']) & (this_order_vec_raw <= tbl[irow]['EOrder'])
        nsolns_good = np.sum(igood)
        this_order_vec = this_order_vec_raw[igood]
        indx = this_order_vec_raw - order_min
        indx_good = this_order_vec - order_min
        nsolns = this_order_vec_raw.size
        # Information for the file is stored for convenience, although this is redundant with the arrays below
        table_xidl['filename'][irow] = tbl[irow]['Name']
        table_xidl['nsolns'][irow] = nsolns
        table_xidl['nsolns_good'][irow] = nsolns_good
        table_xidl['bluest_order'][irow] = this_order_vec[-1]
        table_xidl['bluest_good_order'][irow] = this_order_vec_raw[-1]
        table_xidl['reddest_order'][irow] = this_order_vec[0]
        table_xidl['reddest_good_order'][irow] = this_order_vec_raw[0]
        table_xidl['xdisp_file'][irow] = tbl[irow]['XDISP']
        table_xidl['ech_angle_file'][irow] = tbl[irow]['ECH']
        table_xidl['xd_angle_file'][irow] = tbl[irow]['XDAng']
        table_xidl['det_file'][irow] = tbl[irow]['Chip']
        # Arrays (nfile, norders)
        table_xidl['order'][irow, indx] = this_order_vec_raw
        table_xidl['populated'][irow, indx] = True
        table_xidl['populated_and_good'][irow, indx_good] = True
        table_xidl['ech_angle'][irow, indx] = tbl[irow]['ECH']
        table_xidl['xd_angle'][irow, indx] = tbl[irow]['XDAng']
        table_xidl['xdisp'][irow, indx] = tbl[irow]['XDISP']
        table_xidl['det'][irow, indx] = tbl[irow]['Chip']
        table_xidl['binspec'][irow, indx] = tbl[irow]['Rbin']
        table_xidl['lambda_cen'][irow, indx] = np.median(this_wave, axis=1)
        table_xidl['wave'][irow, indx, :] = this_wave
        table_xidl['arcspec'][irow, indx, :] = this_arc
        # Fit the wavelengths
        this_coeff_array = np.zeros((nsolns, n_final + 1))
        for ii, iwave in enumerate(this_wave[igood, :]):
            pypeitFit = fitting.robust_fit(xvec, iwave, n_final, function=func, maxiter=10,
                                           lower=1e10, upper=1e10, maxrej=0, sticky=True,
                                           minx=fmin, maxx=fmax, weights=None)
            this_coeff_array[ii, :] = pypeitFit.fitc
        table_xidl['coeff'][irow, indx, :] = this_coeff_array

    # Write out to multi-extension fits
    print(f'Writing HIRES xidl wv_calib archive to file: {outfile}')
    hdu_param = fits.BinTableHDU(params.as_array())
    hdu_table = fits.BinTableHDU(table_xidl.as_array())

    hdulist = fits.HDUList()
    hdulist.append(hdu_param)
    hdulist.append(hdu_table)
    hdulist.writeto(outfile, overwrite=True)

    return
