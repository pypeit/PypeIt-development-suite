
import argparse
import sys

import scipy
import yaml
from copy import deepcopy
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from functools import partial
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy import table
from astropy import stats

from pypeit.core import fitting, flux_calib
from pypeit.core.telluric import ZP_UNIT_CONST
from pypeit.core.wavecal import wvutils
from pypeit.spectrographs.util import load_spectrograph
from pypeit.sensfunc import SensFunc
from pypeit import msgs
from pypeit import io
from pypeit import utils

from IPython import embed

def get_std_dict(sensfuncs):
    """Get the standard star dictionary from the sensitivity functions.

    Args:
        sensfuncs (list):
            List of SensFunc objects.

    Returns:
        dict: Dictionary with the standard star information.
    """

    cal_names = []
    cal_ras = []
    cal_decs = []
    for sensobj in sensfuncs:
        if sensobj.std_name not in cal_names:
            cal_names.append(sensobj.std_name)
            cal_ras.append(sensobj.std_ra)
            cal_decs.append(sensobj.std_dec)

    all_std_dicts = {}
    for i, cal_name in enumerate(cal_names):
        std_dict = flux_calib.get_standard_spectrum(ra=cal_ras[i], dec=cal_decs[i])
        all_std_dicts[cal_name] = std_dict

    return all_std_dicts


# PLOTS RELATED FUNCTIONS #
def color_distance(c1, c2):
    return np.sqrt(sum((a - b) ** 2 for a, b in zip(mcolors.to_rgb(c1), mcolors.to_rgb(c2))))


def is_blackish_or_whiteish(color, bw_th=0.8):
    rgb = np.array(mcolors.to_rgb(color))
    return np.all(rgb > bw_th) or np.all(rgb < (1 - bw_th))


def get_unique_colors(n, color_dict):
    colors = list(color_dict.keys())
    unique_colors = []
    bw_th = 0.8
    while bw_th < 0.99 and len(unique_colors) < n:
        filtered_colors = [color for color in colors if not is_blackish_or_whiteish(color, bw_th)]
        threshold = 0.3
        while len(unique_colors) < n and threshold > 0:
            unique_colors = list(mcolors.TABLEAU_COLORS)
            for color in filtered_colors:
                if all(color_distance(color, uc) > threshold for uc in unique_colors):
                    unique_colors.append(color)
                if len(unique_colors) == n:
                    break
            threshold -= 0.01
        bw_th += 0.01

    if len(unique_colors) < n:
        raise ValueError("Not enough unique colors available")
    return unique_colors


# def on_pick(event, selected_lines, key, key_type='order'):
#     # Get the line that was picked
#     line = event.artist
#     label = line.get_label().split(' - ')[0]
#
#     # Add or remove the label to/from the dictionary under the specific key
#     if key not in selected_lines:
#         selected_lines[key] = []
#     if label in selected_lines[key]:
#         selected_lines[key].remove(label)
#         print(f"Removed: {label} from {key_type}: {key}")
#     else:
#         selected_lines[key].append(label)
#         print(f"Selected: {label} in {key_type}: {key}")

def on_pick(event, selected_lines, removed_lines, key, key_type='order'):
    # Get the line that was picked
    line = event.artist
    label = line.get_label().split(' - ')[0]
    ax = event.canvas.figure.gca()  # Get the current Axes instance

    # Remove the label from the dictionary under the specific key and from the plot
    if key in selected_lines.keys() and label in selected_lines[key]:
        # remove from plot
        line.set_visible(False)
        # update legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, fontsize=2)
        event.canvas.draw()
        # remove from dictionary
        removed_lines.append((label, line))
        selected_lines[key].remove(label)
        # print info
        print(f"\n-- Removed {label} from {key_type}: {key}")
        print(f"\n  Remaining selected lines:")
        for l,_lab in enumerate(selected_lines[key]):
            print(f"  {l+1}) {_lab}")
        print('')


def undo(event, selected_lines, removed_lines, key, key_type='order'):
    if event.key == 'z' and removed_lines:
        # recover removed lines and make them visible in the plot
        label, line = removed_lines.pop()
        line.set_visible(True)
        ax = event.canvas.figure.gca()  # Get the current Axes instance
        # update legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, fontsize=5)
        event.canvas.draw()
        # re-add the line in the dictionary
        if key in selected_lines.keys():
            selected_lines[key].append(label)
        else:
            selected_lines[key] = [label]
        # print info
        print(f"++ Re-Added {label} in {key_type}: {key}")

###########


def colormap(objlist):
    """Create a colormap for a list of objects.

    Args:
        objlist (list or `~numpy.ndarray`):
            List of objects.

    Returns:
        dict: Dictionary with the objects as keys and the colors as values.

    """
    unique_colors = get_unique_colors(len(objlist), mcolors.CSS4_COLORS)
    color_map = {obj: unique_colors[i] for i, obj in enumerate(objlist)}
    return color_map


def print_datasets(sensfuncs):
    """ ONLY FOR DEBUGGING. Print info on the sensfuncs.

    Args:
        sensfuncs (list):
            List of SensFunc objects.

    Returns:
        Table: Table with the information of the SensFunc objects.

    """
    # for debugging
    name = []
    airmass = []
    exptime = []
    order_min = []
    order_max = []
    no_gap_orders = []
    for sensobj in sensfuncs:
        name.append(Path(sensobj.spec1df).name.split('-')[0].split('_')[1]+'.fits')
        airmass.append(sensobj.airmass)
        exptime.append(sensobj.exptime)
        order_min.append(sensobj.sens['ECH_ORDERS'].data.min())
        order_max.append(sensobj.sens['ECH_ORDERS'].data.max())
        no_gap_orders.append(np.all(np.diff(sensobj.sens['ECH_ORDERS'].data) == -1))
    tab = Table()
    tab['name'] = name
    tab['airmass'] = airmass
    tab['exptime'] = exptime
    tab['order_min'] = order_min
    tab['order_max'] = order_max
    tab['no_gap_orders'] = no_gap_orders
    stars_tab = Table.read('../hires_std_stars - hires_std_star_files.csv', format='csv')
    stars_tab.keep_columns(
        ['koaid', 'ra', 'dec', 'progpi', 'airmass', 'waveblue', 'wavered', 'deckname', 'xdispers', 'fil1name',
         'guidfwhm', 'echangl', 'xdangl', 'binning', 'Std Stars', 'skyprobe extinction'])
    tab_all = table.join(tab, stars_tab, keys_left='name', keys_right='koaid')
    tab_all.sort(['airmass_1', 'skyprobe extinction'])
    tab_all.write('sensfunc_info.csv', format='csv', overwrite=True)
    tab_all.pprint_all()
    return tab_all


def get_sname_from_sensobj(sensobj):
    """Get the sensitivity function file name from the SensFunc object.

    Args:
        sensobj (:obj:`SensFunc`):
            SensFunc object.

    Returns:
        str: Sensitivity function file name.

    """
    return Path(sensobj.spec1df).name.replace('spec1d', 'sens').replace(' (1)', '')


def match_sensobj2sensfile(sensobjs, sensnames):
    """Match the SensFunc objects to the sensitivity function file names.

    Args:
        sensobjs (list):
            List of SensFunc objects.
        sensnames (list):
            List of sensitivity function file names.

    Returns:
        list: list of indexes to match the SensFunc objects to the sensitivity function file names.

    """
    mtch = []
    for n in sensnames:
        for s, sf in enumerate(sensobjs):
            if get_sname_from_sensobj(sf) == n.replace(' (1)', ''):
                mtch.append(s)
    return mtch

def create_sens_files(spec1d_files, spec1d_files_path, sens_files_path, boxcar=False,
                      use_flat=True, skip_existing=True):
    """Generate sensitivity functions from the spec1d files provided in the list.

    Args:
        spec1d_files (list):
            List of spec1d file names used to generate the sensitivity functions.
        spec1d_files_path (str or Path):
            Path to the spec1d files.
        sens_files_path (str or Path):
            Path to save the sensitivity functions.
        boxcar (bool):
            Use boxcar extraction for the sensitivity function computation.
        use_flat (bool):
            Use flat fielding for the sensitivity function computation.
        skip_existing (bool):
            Skip computing the sensitivity function if the output file already exists.
        parse (str):
            Parse the sensitivity functions. If None, no parsing will be done
            only the plotting will be done. If 'use_all', all the sensitivity functions
            will be parsed per each order. If 'select' the user will be able to select
            in a matplotlib plot the sensitivity functions to use per each order.
        plot_all (bool):
            Plot the sensitivity functions.
        ptype (str):
            Type of plot to generate. It can be 'zeropoint', 'throughput', or 'counts_per_angs'.

    Returns:
        tuple: Dictionary with the selected SensFunc objects and the selected sensitivity function file names
            to use in the combined sensitivity function per ech order.

    """

    # check if spec1d_files_path and sens_files_path are Path objects
    if not isinstance(spec1d_files_path, Path):
        spec1d_files_path = Path(spec1d_files_path).absolute()
    if not isinstance(sens_files_path, Path):
        sens_files_path = Path(sens_files_path).absolute()

    # check if the paths exist
    if not spec1d_files_path.exists():
        msgs.error(f"Spec1d files path {spec1d_files_path} does not exist.")
    if not sens_files_path.exists():
        msgs.error(f"Sensitivity functions path {sens_files_path} does not exist.")

    # load the spectrograph
    spectrograph = load_spectrograph("keck_hires")
    par = spectrograph.default_pypeit_par()
    par['sensfunc']['extrap_blu'] = 0.
    par['sensfunc']['extrap_red'] = 0.
    if use_flat:
        par['sensfunc']['use_flat'] = True
    if boxcar:
        par['sensfunc']['extr'] = 'BOX'
    # this is just for the QAplot
    par['fluxcalib']['extrap_sens'] = True

    par_outfile = sens_files_path / "sensfunc.par"
    print(f'Writing the sensfunc parameters to {par_outfile}')
    par['sensfunc'].to_config(par_outfile, section_name='sensfunc', include_descr=False)

    sensfuncs = []
    sensnames = []

    orders = np.array([])

    use_spec1dfiles_list = []
    all_orders = np.array([])
    for spec1d_file in spec1d_files:
        sens_file = sens_files_path / spec1d_file.replace('spec1d', 'sens')
        if sens_file.exists() and skip_existing:
            print(f'Sensitivity function {sens_file} already exists. Skipping.')
            continue
        try:
            # compute the sensitivity function
            spec1d = spec1d_files_path / spec1d_file
            sensobj = SensFunc.get_instance(str(spec1d), str(sens_file), par['sensfunc'],
                                            par_fluxcalib=par['fluxcalib'], chk_version=False)

            # Generate the sensfunc
            sensobj.run()
            msgs.info(f'Sensitivity function for {spec1d_file} computed.')
            sensfuncs.append(sensobj)
            orders = np.append(orders, sensobj.sens['ECH_ORDERS'].data)
            sensnames.append(sens_file.name)
            use_spec1dfiles_list.append(spec1d_file)
            # write the sensitivity function to a file
            save_sensobj_to_fits(sensobj, sens_file, fits.open(spec1d))

        except Exception as e:
            msgs.warn(f'Error computing sensitivity function for {spec1d_file}.\n{e}')

    if len(sensfuncs) == 0:
        msgs.error("No sensitivity functions loaded.")

    # save to file
    use_spec1dfiles_file = f'used_spec1ds_v{fits.getval(sensfuncs[0], "DMODVER", 1)}.txt'
    with open(use_spec1dfiles_file, 'w') as f:
        for fname in use_spec1dfiles_list:
            f.write(f"{fname}\n")


def load_sensfunc(sens_files_path, sens_fnames=None, sens_fnames_dict=None,
                  parse=None, plot_all=False, ptype='counts_per_angs'):
    """Load the sensitivity functions from a list of files.

    Args:
        sens_files_path (str or Path):
            Path to the sensitivity function files.
        sens_fnames (list, optional):
            List of sensitivity function file names.
        sens_fnames_dict (dict, optional):
            Dictionary with the orders as keys and the
            sensitivity function file names as values.
        parse (str):
            Parse the sensitivity functions. If None, no parsing will be done
            only the plotting will be done. If 'use_all', all the sensitivity functions
            will be parsed per each order. If 'select' the user will be able to select
            in a matplotlib plot the sensitivity functions to use per each order.
        plot_all (bool):
            Plot the sensitivity functions.
        ptype (str):
            Type of plot to generate. It can be 'zeropoint', 'throughput', or 'counts_per_angs'.

    Returns:
        tuple: Dictionary with the selected SensFunc objects and the selected sensitivity function file names
            to use in the combined sensitivity function per ech order.
    """

    # check if sens_files_path is a Path object
    if not isinstance(sens_files_path, Path):
        sens_files_path = Path(sens_files_path).absolute()
    # check if the path exists
    if not sens_files_path.exists():
        msgs.error(f"Sensitivity functions path {sens_files_path} does not exist.")

    if sens_fnames is None and sens_fnames_dict is None:
        msgs.error("No sensitivity function files provided. "
                   "Please provide a list of files or a dictionary with the orders and the files.")
    elif sens_fnames is not None and sens_fnames_dict is not None:
        msgs.warn("Both a list of files and a dictionary provided. "
                   "Using the dictionary and ignoring the list of files.")
        sens_fnames = np.unique([item for sublist in sens_fnames_dict.values() for item in sublist if sublist])
    elif sens_fnames is None and sens_fnames_dict is not None:
        sens_fnames = np.unique([item for sublist in sens_fnames_dict.values() for item in sublist if sublist])

    sensnames = []

    for sens_file in sens_fnames:
        if 'HZ44' in sens_file:
            continue
        sens_file = sens_files_path / sens_file
        if not sens_file.exists():
            msgs.warn(f'Sensitivity function file {sens_file} does not exist.')
        else:
            sensnames.append(sens_file.name)

    if len(sensnames) == 0:
        msgs.error("No sensitivity functions loaded.")

    zps, zp_scales, waves, ords, snames = plot_parse_loaded(sensnames, sens_files_path,
                                                            selected_lines=sens_fnames_dict,
                                                            parse=parse, plot_all=plot_all, ptype=ptype)

    if len(zps) == 0:
        msgs.error("No sensitivity functions loaded. Try setting the parameter --parse.")
    return zps, zp_scales, waves, ords, snames


def plot_parse_loaded(sensnames, sens_files_path, cut_left=100, cut_right=-100,
                      selected_lines=None, parse=None, plot_all=False,
                      savetofile=False, ptype='counts_per_angs'):
    """Plot all the loaded sensitivity functions.

    Args:
        sensnames (list):
            List of sensitivity function file names.
        sens_files_path (Path):
            Path to the sensitivity function files.
        selected_lines (dict):
            Dictionary with the selected sensitivity function file names per ech order.
        parse (str):
            Parse the sensitivity functions. If None, no parsing will be done
            only the plotting will be done. If 'use_all', all the sensitivity functions
            will be parsed per each order. If 'select' the user will be able to select
            in a matplotlib plot the sensitivity functions to use per each order.
        plot_all (bool):
            Plot all the sensitivity functions.
        savetofile (bool):
            Save the selected sensitivity functions to a yaml file.
        ptype (str):
            Type of plot to generate. It can be 'zeropoint', 'throughput', or 'counts_per_angs'.

    """

    # define the order vector. HARD CODED FOR NOW
    order_vec = np.arange(93, 34, -1)
    # get SensFunc objects
    sensfuncs = [SensFunc.from_file(sens_files_path / sn, chk_version=False) for sn in sensnames]

    # colors by SensFunc object
    color_map_sens = colormap(sensnames)
    # color by order
    color_map_ord = colormap(order_vec)
    y_label = 'Zeropoint (AB mag)' if ptype == 'zeropoint' else 'Throughput' if ptype == 'throughput' \
        else 'Counts /(Angstrom s)'
    outfile_name = f'collection_{ptype}_by_order.pdf'

    y_all = []
    y_smooth_all = []
    y_scale_all = []
    y_maxmax_all = []
    y_med_all = []
    wave_all = []
    order_all = []
    sname_all = []
    with PdfPages(outfile_name) as pdf:
        for iord in order_vec:
            # get the SensFunc objects for the current order
            if selected_lines is not None and iord in selected_lines.keys() and len(selected_lines[iord]) > 0:
                i_sensnames = selected_lines[iord]
                mtch = match_sensobj2sensfile(sensfuncs, i_sensnames)
                i_sensfuncs = [sensfuncs[m] for m in mtch]
            else:
                i_sensfuncs = sensfuncs
                i_sensnames = sensnames

            # reset selected_lines if it already exists, otherwise create it
            if selected_lines is None:
                selected_lines = {}
            selected_lines[iord] = []

            # plot
            fig = plt.figure(figsize=(23, 6.))
            plt.minorticks_on()
            plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
            y_iord = []
            y_iord_smooth = []
            y_iord_smooth_max = []
            wave_iord = []
            wave_iord_min = []
            wave_iord_max = []
            legend_iord = []
            color_iord = []
            for s, sensobj in enumerate(i_sensfuncs):
                _indx = np.where(sensobj.sens['ECH_ORDERS'].data == iord)[0]
                if _indx.size > 0:
                    indx = _indx[0]
                    sname = i_sensnames[s]
                    airmass = sensobj.airmass
                    exptime = sensobj.exptime
                    decker = fits.getval(sens_files_path / sname, 'DECKER')
                    filter = fits.getval(sens_files_path / sname, 'FILTER1')
                    echangle = fits.getval(sens_files_path / sname, 'ECHANGLE')
                    xdangle = fits.getval(sens_files_path / sname, 'XDANGLE')

                    color = color_map_sens[sname]
                    wave = sensobj.sens['SENS_WAVE'].data[indx]
                    wave_gmp = wave > 1.0
                    wmin = wave[wave_gmp].min()
                    wmax = wave[wave_gmp].max()
                    if ptype == 'zeropoint':
                        # zeropoint_data = sensobj.sens['SENS_ZEROPOINT'].data[indx]
                        # plt.plot(wave[wave_gmp], zeropoint_data[wave_gmp], alpha=0.4, color=color, lw=0.6,
                        #          zorder=-2)
                        y = sensobj.sens['SENS_ZEROPOINT_FIT'].data[indx]
                        zcut = y > 12.
                        y = y[zcut]
                        wave = wave[zcut]
                        wave_gmp = wave > 1.0
                        y_tresh = 20.
                    elif ptype == 'throughput':
                        wave = sensobj.wave[:, indx]
                        wcut = (wave >= wmin) & (wave <= wmax)
                        wave = wave[wcut]
                        wave_gmp = wave > 1.0
                        y = sensobj.throughput[:, indx][wcut]
                        y_tresh = 0.2
                    else:
                        # counts per Angstrom
                        y = sensobj.sens['SENS_COUNTS_PER_ANG'].data[indx] / exptime
                        y_tresh = 10000.

                    # smooth y
                    filt = 100
                    y_smooth = utils.fast_running_median(y[wave_gmp], filt)
                    y_smooth_max = np.nanmax(y_smooth[y_smooth < y_tresh])
                    # append by order
                    y_iord.append(y[wave_gmp])
                    y_iord_smooth.append(y_smooth)
                    y_iord_smooth_max.append(y_smooth_max)
                    wave_iord.append(wave[wave_gmp])
                    y_mintresh = 0.0
                    # if ptype == 'zeropoint':
                    #     y_mintresh = 13. if iord in [93, 92, 35] else 16.
                    _w = wave[wave_gmp][y_smooth > y_mintresh]
                    wave_iord_min.append(_w.min())
                    wave_iord_max.append(_w.max())
                    color_iord.append(color)
                    # legend_iord.append(
                    #     f'{sname} - {decker} - {filter} - ech: {echangle:.3f} - xd: {xdangle:.3f}\n'
                    #     f'- airmass: {airmass:.2f} - exptime: {exptime:.1f}')
                    legend_iord.append(f'{sname} - {decker} - {filter} - ech: {echangle:.3f} - xd: {xdangle:.3f}')

                    # append to all
                    y_all.append(y[wave_gmp])
                    y_smooth_all.append(y_smooth)
                    wave_all.append(wave[wave_gmp])
                    order_all.append(iord)
                    sname_all.append(sname)
                    selected_lines[iord].append(sname)
            if len(y_iord) == 0:
                plt.close(fig)
                continue
            # get the wavelength range for the current order that is common to all the sensitivity functions
            if iord == 93:
                wave_mid = 3840.
            else:
                wave_min = np.max(wave_iord_min)
                wave_max = np.min(wave_iord_max)
                wave_mid = (wave_min + wave_max) / 2
            # find the value of the smoothed y at wave_mid if it exists otherwise use the max value
            y_iord_smooth_ref = []
            for i in range(len(wave_iord)):
                # find the index of the wave_mid in the wave_iord array
                diff = np.abs(wave_iord[i] - wave_mid)
                if np.min(diff) < 5.:
                    ind = np.argmin(diff)
                    y_iord_smooth_ref.append(y_iord_smooth[i][ind])
                else:
                    y_iord_smooth_ref.append(y_iord_smooth_max[i])

            # get the max value of the smoothed y
            y_maxmax = np.nanmax(y_iord_smooth_ref)
            y_med = stats.sigma_clipped_stats(y_iord_smooth_ref, sigma=3, maxiters=5)[1]
            # append to all
            y_maxmax_all.append(y_maxmax)
            y_med_all.append(y_med)

            # plot by order
            for y_i, y_sm, y_ref, wave, color, legend in zip(y_iord, y_iord_smooth, y_iord_smooth_ref, wave_iord,
                                                             color_iord, legend_iord):
                s = y_med / y_ref
                legend += f' - scale: {s:.2f}'
                # if iord != 93:
                #     s = 1
                plt.plot(wave, y_i * s, color=color, alpha=0.3, lw=0.5, zorder=0)
                plt.plot(wave, y_sm * s, color=color, lw=1.5, zorder=1, label=legend, picker=True)
                y_scale_all.append(s)
            plt.axvline(wave_mid, color='k', ls='--', lw=0.5)
            plt.axhline(y_med, color='k', ls='--', lw=0.5)
            plt.title(f'Order {iord}')
            plt.xlabel('Wavelength (Angstroms)')
            plt.ylabel(y_label)
            if ptype == 'zeropoint':
                plt.ylim(y_med*0.8, y_maxmax * 1.1)
            elif ptype == 'throughput':
                plt.ylim(0.0,  y_maxmax * 1.2)
            else:
                plt.ylim(0.0, y_maxmax * 1.2)
            plt.legend(fontsize=2, loc='upper right')
            fig.tight_layout()
            if parse == 'select':
                removed_lines = []
                fig.canvas.mpl_connect('pick_event', partial(on_pick, selected_lines=selected_lines,
                                                             removed_lines=removed_lines, key=iord))
                fig.canvas.mpl_connect('key_press_event', partial(undo, selected_lines=selected_lines,
                                                                  removed_lines=removed_lines, key=iord))
                # sort the selected lines
                selected_lines[iord] = sorted(selected_lines[iord])

                print('')
                plt.show()
            if plot_all:
                pdf.savefig(dpi=50)
            plt.close(fig)
        if plot_all:
            # plot all orders
            fig = plt.figure(figsize=(23, 6.))
            plt.minorticks_on()
            plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
            plt.title(f'All orders')
            in_order = None
            for y_a, y_sm_a, s_a, wave_a, o_a, sn_a in zip(y_all, y_smooth_all, y_scale_all, wave_all, order_all, sname_all):
                if sn_a not in selected_lines[o_a]:
                    continue
                color = color_map_ord[o_a]
                leg = f'order {o_a}' if o_a != in_order else ''
                # plt.plot(wave_a, y_a * s_a, lw=0.5, color=color, alpha=0.3, zorder=0, label=leg)
                plt.plot(wave_a, y_sm_a * s_a, lw=1.5, color=color, zorder=1, label=leg)
                in_order = o_a

            plt.xlabel('Wavelength (Angstroms)')
            plt.ylabel(y_label)
            plt.legend(fontsize=5, loc='upper right')
            if ptype == 'zeropoint':
                plt.ylim(np.nanmin(y_med_all) * 0.95, np.nanmax(y_maxmax_all) * 1.05)
            elif ptype == 'throughput':
                plt.ylim(0.0, np.nanmax(y_maxmax_all) * 1.2)
            else:
                plt.ylim(0.0, np.nanmax(y_maxmax_all) * 1.2)
            fig.tight_layout()
            pdf.savefig(dpi=50)
            plt.close(fig)

    if parse == 'select' or savetofile:
        # save to file
        dmver = fits.getval(sens_files_path / sensnames[0], 'DMODVER', 2)
        parsed_sens_file = f'used_sensfuncs_v{dmver}.yaml'
        _parsed_sensnames = {str(key): value for key, value in selected_lines.items()}
        with open(parsed_sens_file, 'w') as f:
            yaml.dump(_parsed_sensnames, f)

    if ptype in ['throughput', 'counts_per_angs'] or parse == 'select':
        sys.exit('Exiting after parsing the sensitivity functions.' if parse == 'select' else 'Exiting after plotting.')
    else:
        return y_all, y_scale_all, wave_all, order_all, sname_all


def combine_sensfuncs(zps, zp_scales, waves, ords, snames, sens_file,
                      cut_left=50, cut_right=-100, poly_order=4, qa_plots=False):

    order_vec = np.unique(ords)[::-1]
    if np.any(np.diff(order_vec) != -1):
        msgs.warns("Orders are not contiguous. There might be missing orders.")

    # get the wavelength grid
    wave_grid, wave_grid_mid, dsamp = wvutils.get_wave_grid(waves, wave_method='velocity',
                                                            wave_grid_min=0.98*np.concatenate(waves).min(),
                                                            wave_grid_max=1.02*np.concatenate(waves).max())

    comb_zeropoints = []
    comb_waves = []
    comb_orders = []
    comb_coeff = []

    for iord in order_vec:
        # get the SensFunc objects that will be used for the current order
        i_zpoints = np.concatenate([zps[i] * zp_scales[i] for i, o in enumerate(ords) if o == iord])
        i_waves = np.concatenate([waves[i] for i, o in enumerate(ords) if o == iord])
        # fit the zeropoints
        func = 'legendre'
        pypeitFit = fitting.robust_fit(i_waves, i_zpoints, poly_order, function=func, minx=i_waves.min(),
                                       maxx=i_waves.max(), in_gpm=i_zpoints>9., lower=10, upper=10, use_mad=True)
        i_comb_zpoint = fitting.evaluate_fit(pypeitFit.fitc, func, wave_grid, minx=i_waves.min(), maxx=i_waves.max())
        i_comb_wave = wave_grid.copy()
        bad_pixs = (i_comb_wave < i_waves.min()) | (i_comb_wave > i_waves.max())
        i_comb_zpoint[bad_pixs] = 0.
        i_comb_wave[bad_pixs] = 0.

        # append the results to the lists
        comb_zeropoints.append(i_comb_zpoint)
        comb_waves.append(i_comb_wave)
        comb_orders.append(iord)
        comb_coeff.append(pypeitFit.fitc)

    comb_sensobj = fill_combSensObj(comb_zeropoints, comb_waves, comb_orders, poly_order, comb_coeff, sens_file)

    if qa_plots:
        plot_conbined_sensfunc(comb_sensobj, zps, zp_scales, waves, ords, snames)

    return comb_sensobj


def plot_conbined_sensfunc(comb_sensobj, zps, zp_scales, waves, ords, snames):
    """Plot the combined sensitivity function.

    Args:
        comb_sensobj (:obj:`SensFunc`):
            Combined sensitivity function object.
        zps (list):
            List of zeropoints.
        zp_scales (list):
            List of zeropoint scales.
        waves (list):
            List of wavelength arrays.
        ords (list):
            List of orders.
        snames (list):
            List of sensitivity function file names.

    """
    order_vec = comb_sensobj.sens['ECH_ORDERS'].data

    # colors by SensFunc object
    color_map_sens = colormap(list(np.unique(snames)))
    # color by order
    color_map_ord = colormap(list(order_vec))

    outfile_name = 'combined_sensfunc.pdf'
    factor = 10  # this is for spacing the data points in the plot
    with PdfPages(outfile_name) as pdf:
        for iord in order_vec:
            fig = plt.figure(figsize=(23, 6.))
            plt.minorticks_on()
            plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
            i_zpoints = [zps[i] * zp_scales[i] for i, o in enumerate(ords) if o == iord]
            i_waves = [waves[i] for i, o in enumerate(ords) if o == iord]
            i_snames = [snames[i] for i, o in enumerate(ords) if o == iord]
            if len(i_zpoints) == 0:
                plt.close(fig)
                continue
            for (i_zpoint, i_wave, i_sname) in zip(i_zpoints, i_waves, i_snames):
                color = color_map_sens[i_sname]
                plt.plot(i_wave[::factor], i_zpoint[::factor], color=color, ls='', marker='.',
                         ms=2, alpha=0.2, zorder=-2, label=i_sname)

            # plot the combined sensitivity function
            indx = np.where(comb_sensobj.sens['ECH_ORDERS'].data == iord)[0]
            if indx.size > 0:
                i_comb_zpoint = comb_sensobj.sens['SENS_ZEROPOINT'].data[indx]
                i_comb_wave = comb_sensobj.sens['SENS_WAVE'].data[indx]
                plt.plot(i_comb_wave[i_comb_wave>0], i_comb_zpoint[i_comb_wave>0],
                         color='k', lw=3., zorder=1, label='Combined')
            plt.title(f'Order {iord}')
            plt.xlabel('Wavelength (Angstroms)')
            plt.ylabel('Zeropoint (AB mag)')
            plt.ylim(i_comb_zpoint[i_comb_wave>0].min()*0.95, i_comb_zpoint[i_comb_wave>0].max()*1.05)
            #plt.legend(fontsize=4, loc='upper right')
            fig.tight_layout()
            pdf.savefig(dpi=50)
            plt.close(fig)
        # plot all orders
        fig = plt.figure(figsize=(23, 6.))
        plt.minorticks_on()
        plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
        plt.title(f'All orders - Zeropoint')
        in_order = None
        for i_zpoint, i_zp_scale, i_wave, i_sname, i_order in zip(zps, zp_scales, waves, snames, ords):
            color = color_map_ord[i_order]
            leg = f'order {i_order}' if i_order != in_order else ''
            plt.plot(i_wave, i_zpoint * i_zp_scale, lw=0.5, color=color, alpha=0.3, zorder=0, label=leg)
            in_order = i_order
        w_comb = comb_sensobj.wave.T
        z_comb = comb_sensobj.zeropoint.T
        for c, (cw, cz) in enumerate(zip(w_comb, z_comb)):
            plt.plot(cw[cw>1], cz[cw>1], lw=1.5, color='k', zorder=1, label='Combined' if c == 0 else '')
        plt.xlabel('Wavelength (Angstroms)')
        plt.ylabel('Zeropoint (AB mag)')
        plt.ylim(z_comb[w_comb > 0].min()*0.95, z_comb[w_comb > 0].max()*1.05)
        plt.legend(fontsize=4, loc='upper right')
        fig.tight_layout()
        pdf.savefig(dpi=50)
        plt.close(fig)
        # plot combined throughput
        fig = plt.figure(figsize=(23, 6.))
        plt.minorticks_on()
        plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
        plt.title(f'All orders - Throughput')
        w_tru = comb_sensobj.wave.T
        t_tru = comb_sensobj.throughput.T
        in_order = None
        for (o, wave, throughput) in zip(order_vec, w_tru, t_tru):
            color = color_map_ord[o]
            leg = f'order {o}' if o != in_order else ''
            plt.plot(wave[wave>1], throughput[wave>1], lw=1, color=color, zorder=0, label=leg)
            in_order = o
        plt.xlabel('Wavelength (Angstroms)')
        plt.ylabel('Throughput')
        plt.ylim(t_tru[w_tru > 0].min()*0.95, t_tru[w_tru > 0].max()*1.05)
        plt.legend(fontsize=4, loc='upper right')
        fig.tight_layout()
        pdf.savefig(dpi=50)
        plt.close(fig)


def fill_combSensObj(comb_zeropoint_fit_list, comb_wave_list, comb_orders_list, poly_order, comb_coeff_list, sensfile):
    """Fill the SensFunc object with the combined zeropoints.

    Args:
        comb_zeropoint_fit_list (list):
            List of combined zeropoints for each order.
        comb_wave_list (list):
            List of wavelength arrays for the combined zeropoints.
        comb_orders_list (list):
            List of orders.
        poly_order (int):
            Order of the polynomial to fit the zeropoints.
        comb_coeff_list (list):
            List of coefficients of the polynomial fits to the zeropoints.
        sensfile (Path):
            Path to a sensitivity function file. It will be used to initialize the SensFunc object.

    Returns:
        SensFunc: SensFunc object with the combined zeropoints.
    """

    # Find the maximum length of the arrays
    max_length = max(len(arr) for arr in comb_zeropoint_fit_list)

    # Pad the arrays with zeros to make them the same size
    comb_zeropoint_fit_list = [np.pad(arr, (0, max_length - len(arr)), constant_values=0) for arr in comb_zeropoint_fit_list]
    comb_wave_list = [np.pad(arr, (0, max_length - len(arr)), constant_values=0) for arr in comb_wave_list]

    # Convert lists to 2D arrays
    comb_zeropoint_fit_array = np.array(comb_zeropoint_fit_list)
    comb_wave_array = np.array(comb_wave_list)

    # initialize a SensFunc object with the same attributes as the input sensobj
    sensobj = SensFunc.from_file(sensfile)
    comb_sensobj = deepcopy(sensobj)
    # remove attributes that are not general to all the sensfuncs
    att_to_reinit = ['spec1df', 'std_name', 'std_cal', 'std_ra', 'std_dec', 'airmass', 'exptime', 'sens', 'wave',
                     'zeropoint', 'throughput']
    for att in att_to_reinit:
        setattr(comb_sensobj, att, None)
    tell_att_to_reinit = ['std_src', 'std_name', 'std_cal', 'airmass', 'exptime', 'std_ra', 'std_dec', 'model']
    for att in tell_att_to_reinit:
        setattr(comb_sensobj.telluric, att, None)

    # Save the combined zeropoints to a sensfunc file
    comb_sensobj.sens = comb_sensobj.empty_sensfunc_table(len(comb_orders_list), comb_zeropoint_fit_array.shape[0], 0)
    comb_sensobj.sens['ECH_ORDERS'] = comb_orders_list
    comb_sensobj.sens['SENS_COEFF'] = comb_coeff_list
    comb_sensobj.sens['SENS_WAVE'] = comb_wave_array
    comb_sensobj.sens['SENS_ZEROPOINT'] = comb_zeropoint_fit_array
    comb_sensobj.sens['SENS_ZEROPOINT_GPM'] = comb_wave_array > 1.0
    comb_sensobj.sens['SENS_ZEROPOINT_FIT'] = comb_zeropoint_fit_array
    comb_sensobj.sens['SENS_ZEROPOINT_FIT_GPM'] = comb_wave_array > 1.0
    comb_sensobj.sens['WAVE_MIN'] = np.array([np.min(w[w > 1]) if w[w > 1].size > 0 else 0 for w in comb_wave_array])
    comb_sensobj.sens['WAVE_MAX'] = np.array([np.max(w[w > 1]) if w[w > 1].size > 0 else 0 for w in comb_wave_array])
    comb_sensobj.sens['POLYORDER_VEC'] = poly_order

    comb_sensobj.wave = comb_wave_array.T
    comb_sensobj.zeropoint = comb_zeropoint_fit_array.T

    comb_sensobj.spectrograph = load_spectrograph("keck_hires")
    comb_sensobj.throughput = comb_sensobj.compute_throughput()[0]

    return comb_sensobj


def save_sensobj_to_fits(sensobj, outfile, hdul):
    """Save the SensFunc object to a fits file.

    Args:
        sensobj (SensFunc):
            SensFunc object.
        outfile (str or Path):
            Output file name.
        hdul (HDUList):
            HDUList object of a spec1d file.
    """

    # check if outfile is a Path object
    if not isinstance(outfile, Path):
        outfile = Path(outfile).absolute()

    spectrograph = load_spectrograph(hdul[0].header['PYP_SPEC'])
    spectrograph_config_par = spectrograph.config_specific_par(hdul)

    # Construct a primary FITS header that includes the spectrograph's
    #   config keys for inclusion in the output sensfunc file
    primary_hdr = io.initialize_header()
    add_keys = (['PYP_SPEC', 'DATE-OBS', 'TELESCOP', 'INSTRUME', 'DETECTOR']
                + spectrograph.configuration_keys() + spectrograph.raw_header_cards())
    for key in add_keys:
        if key.upper() in hdul[0].header.keys():
            primary_hdr[key.upper()] = hdul[0].header[key.upper()]

    # save the SensFunc object to a fits file
    sensobj.to_file(outfile, primary_hdr=primary_hdr, overwrite=True)


def parse_args():
    parser = argparse.ArgumentParser(description='Combine HIRES sensitivity functions to create a general purpose one.')
    parser.add_argument("outfile", type=str,
                        help="Path and name of the combined sensitivity function.")
    parser.add_argument("sensfuncs_path", type=str,
                        help="Path to place or read from the sensitivity functions.")
    parser.add_argument("--sens_fnames", type=str, default=None,
                        help="Text file with a list of sensitivity function file names to combine. If not provided, "
                             "all the sensitivity functions in the sensfuncs_path will be combined.")
    parser.add_argument("--spec1ds_path", type=str,
                        help="Path to the spec1d files used to generate the sensitivity functions.")
    parser.add_argument("--spec1d_fnames", type=str, default=None,
                        help="Text file with a list of spec1d file names used to generate the sensitivity functions. "
                             "If not provided, all the spec1d files in the spec1ds_path will be used.")
    parser.add_argument("--QA", action="store_true", default=False,
                        help="Generate and save the QA plots to a pdf file")
    parser.add_argument("--reuse", action="store_true", default=False,
                        help="Reuse existing sensfunc and spec1d files.")
    parser.add_argument("--parse", type=str, choices=['use_all', 'select'], default=None,
                        help="Parse the sensitivity functions. Options are 'use_all' to use all the sensitivity functions, "
                        'in each order, or "select" to interactively select only the one to use in each order.')
    parser.add_argument("--ptype", type=str, choices=['zeropoint', 'throughput', 'counts_per_angs'],
                        default='zeropoint', help="Type of plot to generate. It can be 'zeropoint', "
                                                  "'throughput', or 'counts_per_angs'.")
    parser.add_argument("--boxcar", action="store_true", default=False,
                        help="Use boxcar extraction for the sensitivity function computation.This is used only when "
                             "--reuse is not set.")
    parser.add_argument("--skip_existing", action="store_true", default=False,
                        help="Skip computing the sensitivity functions that already exists. This is used only when "
                             "--reuse is not set.")
    parser.add_argument("--use_flat", action="store_true", default=False,
                        help="Use flat fielding for the sensitivity functions computation. This is used only when "
                             "--reuse is not set.")

    return parser.parse_args()


def main(args):

    sensfuncs_path = Path(args.sensfuncs_path).absolute()
    if not sensfuncs_path.exists():
        msgs.error(f"Sensitivity functions path {sensfuncs_path} does not exist.")

    if args.reuse:
        sens_fnames_dict = None
        if args.sens_fnames is not None:
            # this needs to be a dictionary with order as key and sensfunc file name as value
            with open(args.sens_fnames, 'r') as file:
                sens_fnames_dict = yaml.safe_load(file)
            # convert the keys to integers
            sens_fnames_dict = {int(key): value for key, value in sens_fnames_dict.items()}
            sens_fnames = np.unique([item for sublist in sens_fnames_dict.values() for item in sublist if sublist])
        else:
            sens_fnames = [f.name for f in sensfuncs_path.glob('sens*.fits')]
        # load the sensitivity functions
        zps, zp_scales, waves, ords, snames = load_sensfunc(sensfuncs_path, sens_fnames=sens_fnames,
                                                            sens_fnames_dict=sens_fnames_dict,
                                                            parse=args.parse, plot_all=args.QA, ptype=args.ptype)

    else:
        spec1ds_path = Path(args.spec1ds_path).absolute()
        if not spec1ds_path.exists():
            msgs.error(f"Spec1d files path {spec1ds_path} does not exist.")
        if args.spec1d_fnames is not None:
            spec1d_fnames = np.loadtxt(args.spec1d_fnames, dtype=str)
        else:
            spec1d_fnames = [f.name for f in spec1ds_path.glob('spec1d*.fits')]
        if len(spec1d_fnames) == 0:
            msgs.error(f"No spec1d files found.")
        # create the sensitivity functions and save them to disk
        create_sens_files(spec1d_fnames, spec1ds_path, sensfuncs_path,
                          boxcar=args.boxcar, use_flat=args.use_flat, skip_existing=args.skip_existing)
        sens_fnames = [f.name for f in sensfuncs_path.glob('sens*.fits')]
        zps, zp_scales, waves, ords, snames = load_sensfunc(sensfuncs_path, sens_fnames=sens_fnames,
                                                         parse=args.parse, plot_all=args.QA, ptype=args.ptype)

    if len(zps) == 0:
        msgs.error(f"No sensitivity functions found.")


    # for debugging
    # tab_all = print_datasets(sensfuncs)

    # combine the sensitivity functions
    sensfunc = combine_sensfuncs(zps, zp_scales, waves, ords, snames, sensfuncs_path /  sens_fnames[0],
                                 poly_order=4, qa_plots=args.QA)

    # save the combined sensitivity function
    # Construct a primary FITS header
    with io.fits_open(sensfuncs_path / sens_fnames[0]) as hdul:
        save_sensobj_to_fits(sensfunc, args.outfile, hdul)


if __name__ == '__main__':
    args = parse_args()
    args.outfile = 'sens_keck_hires_RED_orders_93-35_sensfunc.fits'
    main(args)


