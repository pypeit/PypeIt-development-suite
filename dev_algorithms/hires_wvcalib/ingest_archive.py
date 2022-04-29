

import os
import numpy as np
import itertools
from astropy.table import Table
from pkg_resources import resource_filename
from matplotlib import pyplot as plt
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.wavecal import templates
from pypeit.core.wavecal import wvutils
from pypeit.core.fitting import robust_fit
from pypeit.core import coadd
from pypeit.core.wavecal import autoid, waveio, wv_fitting
from pypeit import utils
from astropy import table
from scipy import interpolate
from IPython import embed

def get_variable_dlam_wavegrid(lam_min, lam_max, wave_grid_fit, dwave_fit):

    lam_out = [lam_min]
    while lam_out[-1] < lam_max:
        lam_next = lam_out[-1] + np.interp(lam_out[-1], wave_grid_fit, dwave_fit)
        lam_out.append(lam_next)

    return np.array(lam_out, dtype=float)

# Read template file
templ_table_file = os.path.join(
    resource_filename('pypeit', 'data'), 'arc_lines',
    'hires', 'hires_templ.dat')
tbl = Table.read(templ_table_file, format='ascii')
nrows = len(tbl)

order_min = tbl['IOrder'].min()
order_max = 118

order_vec = np.arange(order_min, order_max +1, 1)
norders = order_vec.size

# Subset of orders in every file. Populated indicates whether a given order is populated
lambda_cen = np.zeros((norders, nrows))
ech_angle = np.zeros((norders, nrows))
populated = np.zeros((norders, nrows), dtype=bool)
XDISP_is_red = np.zeros((norders, nrows), dtype=bool)
binspec = np.zeros((norders, nrows), dtype=int)
det = np.zeros((norders, nrows), dtype=int)
xd_angle = np.zeros((norders, nrows))

bluest_order = np.zeros(nrows, dtype=int)
xd_angle_file = np.zeros(nrows)
ech_angle_file = np.zeros(nrows)
det_file = np.zeros(nrows)
XDISP_is_red_file = np.zeros(nrows, dtype=bool)

for irow in np.arange(nrows):
    this_order_vec_raw, this_wave, this_arc = templates.xidl_hires(
        os.path.join(os.getenv('HIRES_CALIBS'), 'ARCS', tbl[irow]['Name']), specbin=tbl[irow]['Rbin'])
    if irow == 0:
        nspec = this_wave.shape[1]
        wave = np.zeros((norders, nrows, nspec))
        arcspec = np.zeros((norders, nrows, nspec))
    else:
        assert this_wave.shape[1] == nspec
    # Restrict to what is labeled as good in the Table
    igood = (this_order_vec_raw >= tbl[irow]['IOrder']) & (this_order_vec_raw <= tbl[irow]['EOrder'])
    this_order_vec = this_order_vec_raw[igood]
    indx = this_order_vec - order_min
    populated[indx, irow] = True
    ech_angle[indx, irow] = tbl[irow]['ECH']
    xd_angle[indx, irow] = tbl[irow]['XDAng']
    XDISP_is_red[indx, irow] = tbl[irow]['XDISP'] == 'RED'
    binspec[indx, irow] =  tbl[irow]['Rbin']
    det[indx, irow] =  tbl[irow]['Chip']

    wave[indx, irow, :] = this_wave[igood, :]
    arcspec[indx, irow, :] = this_arc[igood, :]
    lambda_cen[indx, irow] = np.median(this_wave[igood, :], axis=1)
    # file specific
    bluest_order[irow] = this_order_vec[-1]
    ech_angle_file[irow] = tbl[irow]['ECH']
    xd_angle_file[irow] = tbl[irow]['XDAng']
    det_file[irow] = tbl[irow]['Chip']
    XDISP_is_red_file[irow] = tbl[irow]['XDISP'] == 'RED'

#all_dlam = []
#all_lam = []
#all_orders = []

color_tuple = ('green', 'cyan', 'magenta', 'blue', 'darkorange', 'yellow', 'dodgerblue', 'purple',
               'lightgreen', 'cornflowerblue')
colors = itertools.cycle(color_tuple)

use_unknowns = True
line_lists_all = waveio.load_line_lists(['ThAr'])
line_lists = line_lists_all[np.where(line_lists_all['ion'] != 'UNKNWN')]
unknwns = line_lists_all[np.where(line_lists_all['ion'] == 'UNKNWN')]
tot_line_list = table.vstack([line_lists, unknwns]) if use_unknowns else line_lists
spectrograph = load_spectrograph('keck_hires')
par = spectrograph.default_pypeit_par()['calibrations']['wavelengths']
n_final = 4

# Plot lam vs dlam/lam for each order
debug_all=False
for indx, iorder in enumerate(order_vec):
    if np.any(populated[indx, :]):
        nsolns = np.sum( populated[indx, :])
        this_ech = ech_angle[indx, populated[indx, :]]
        this_xd_angle = xd_angle[indx, populated[indx, :]]
        this_lambda_cen = lambda_cen[indx, populated[indx, :]]
        this_wave = wave[indx, populated[indx, :], :]
        this_arc = arcspec[indx, populated[indx, :], :]

        this_dwave = np.zeros_like(this_wave)
        for ii, iwave in enumerate(this_wave):
            this_dwave[ii, :] = wvutils.get_delta_wave(iwave, (iwave > 0.0))

        # Now try a fit
        med_dlam = np.median(this_dwave[this_wave > 1.0])
        fit = robust_fit(this_wave.flatten(), this_dwave.flatten(), 3, maxiter=25, maxdev = 0.10*med_dlam, groupbadpix=True)
        wave_grid_fit, wave_grid_fit_mid, dsamp = wvutils.get_wave_grid(this_wave.T,wave_method='log10')
        dwave_fit = fit.eval(wave_grid_fit)
        gpm = fit.bool_gpm.copy()
        gpm.resize(this_wave.shape)

        for ii, iwave in enumerate(this_wave):
            this_color=next(colors)
            this_gpm = gpm[ii, :]
            plt.plot(iwave[this_gpm], this_dwave[ii, this_gpm], marker='o', markersize=1.0, mfc=this_color,
                     fillstyle='full',  linestyle='None', zorder=1)
            plt.plot(iwave[np.logical_not(this_gpm)], this_dwave[ii, np.logical_not(this_gpm)], marker='o',
                     markersize=2.0, mfc='red', fillstyle='full', zorder=3, linestyle='None')

        plt.plot(wave_grid_fit, dwave_fit, color='black', label='fit', zorder=10)
        plt.title(f'order={iorder}', fontsize=14)
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

        for slit, iwave in enumerate(this_wave):
            print('Working on soln={:d}'.format(slit))
            # Trim the template to the relevant range. Hack this for now
            itmpl = (wave_grid_mid >= 0.999*iwave.min()) & (wave_grid_mid <= 1.001*iwave.max())
            arcspec_tmpl_trim = arcspec_tmpl[itmpl]
            wave_grid_mid_trim = wave_grid_mid[itmpl]
            arc_in_pad = np.zeros_like(arcspec_tmpl_trim)
            in_gpm = this_arc[slit, :] != 0.0
            npix = np.sum(in_gpm)
            arc_in_pad[:npix] = this_arc[slit, in_gpm]
            detections[str(slit)], spec_cont_sub, all_patt_dict[str(slit)] = autoid.reidentify(
                arc_in_pad, arcspec_tmpl_trim, wave_grid_mid,  tot_line_list, par['nreid_min'],
                cc_thresh=par['cc_thresh'], match_toler=par['match_toler'], cc_local_thresh=par['cc_local_thresh'],
                nlocal_cc=par['nlocal_cc'], nonlinear_counts=1e10,
                sigdetect=par['sigdetect'], fwhm=par['fwhm'], debug_peaks=True, debug_xcorr=True, debug_reid=True)
#            detections[str(slit)], spec_cont_sub, all_patt_dict[str(slit)] = autoid.reidentify(
#                this_arc[slit,:], this_arc[slit+1,:], this_wave[slit+1,:],  tot_line_list, par['nreid_min'],
#                cc_thresh=par['cc_thresh'], match_toler=par['match_toler'], cc_local_thresh=par['cc_local_thresh'],
#                nlocal_cc=par['nlocal_cc'], nonlinear_counts=1e10,
#                sigdetect=par['sigdetect'], fwhm=par['fwhm'], debug_peaks=True, debug_xcorr=True, debug_reid=True)

            # Check if an acceptable reidentification solution was found
            if not all_patt_dict[str(slit)]['acceptable']:
                wv_calib[str(slit)] = None
                continue

            final_fit = wv_fitting.fit_slit(spec_cont_sub, all_patt_dict[str(slit)], detections[str(slit)],
                                         tot_line_list, match_toler=par['match_toler'],func=par['func'], n_first=par['n_first'],
                                         sigrej_first=par['sigrej_first'], n_final=n_final,sigrej_final=par['sigrej_final'])

            autoid.arc_fit_qa(final_fit, title='Silt: {}'.format(str(slit)))


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
for indx, iorder in enumerate(order_vec):
    if np.any(populated[indx, :]):
        this_ech = ech_angle[indx, populated[indx, :]]
        this_xd_angle = xd_angle[indx, populated[indx, :]]
        this_lambda_cen = lambda_cen[indx, populated[indx, :]]
        plt.plot(this_ech, this_lambda_cen, 'k.', label=f'order={iorder}')
        plt.legend()
        plt.show()


for xdisp in ['UV', 'RED']:
    for idet in [1,2,3]:
        indx = (XDISP_is_red_file == (xdisp == 'RED')) & (det_file == idet)
        plt.plot(xd_angle_file[indx], bluest_order[indx], 'k.', label=f'XDISP={xdisp}, det={idet}')
        plt.legend()
        plt.show()





