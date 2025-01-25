
import argparse
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


# PLOTS RELATED FUNCTIONS #
def color_distance(c1, c2):
    return np.sqrt(sum((a - b) ** 2 for a, b in zip(mcolors.to_rgb(c1), mcolors.to_rgb(c2))))


def is_blackish_or_whiteish(color):
    rgb = np.array(mcolors.to_rgb(color))
    return np.all(rgb > 0.8) or np.all(rgb < 0.1)


def get_unique_colors(n, color_dict):
    colors = [color for color in color_dict.keys() if not is_blackish_or_whiteish(color)]
    unique_colors = []
    threshold = 0.3
    while len(unique_colors) < n and threshold > 0:
        unique_colors = []
        for color in colors:
            if all(color_distance(color, uc) > threshold for uc in unique_colors):
                unique_colors.append(color)
            if len(unique_colors) == n:
                break
        threshold -= 0.01
    if len(unique_colors) < n:
        raise ValueError("Not enough unique colors available")
    return unique_colors


def on_pick(event, selected_lines, key, key_type='order'):
    # Get the line that was picked
    line = event.artist
    label = line.get_label().split(' - ')[0]

    # Add or remove the label to/from the dictionary under the specific key
    if key not in selected_lines:
        selected_lines[key] = []
    if label in selected_lines[key]:
        selected_lines[key].remove(label)
        print(f"Removed: {label} from {key_type}: {key}")
    else:
        selected_lines[key].append(label)
        print(f"Selected: {label} in {key_type}: {key}")
###########


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


def create_sens_files(spec1d_files, spec1d_files_path, sens_files_path, boxcar=False,
                      use_flat=True, skip_existing=True, only=None, parse=False):
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
        only (`numpy.ndarray`_):
            Array of files to use (raw filename used).
            If None, all the files in spec1d_files will be used.
        parse (bool):
            Parse the sensitivity functions. If True, all the sensitivity functions
            will be plotted per order and the user will be able to select the ones to use.

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
    deckers = []
    filters = []
    use_spec1dfiles_list = []
    all_orders = np.array([])
    for spec1d_file in spec1d_files:
        if only is not None and spec1d_file.split('-')[0].split('_')[1]+'.fits' not in only:
            continue
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
            sensnames.append(sens_file.name)
            deckers.append(fits.getval(sens_file, 'DECKER'))
            filters.append(fits.getval(sens_file, 'FILTER1'))
            use_spec1dfiles_list.append(spec1d_file)
            all_orders = np.append(all_orders, sensobj.sens['ECH_ORDERS'].data)

        except Exception as e:
            msgs.warn(f'Error computing sensitivity function for {spec1d_file}.\n{e}')

    if len(sensfuncs) == 0:
        msgs.error("No sensitivity functions loaded.")

    # generate the orders vector for the combined sensitivity function
    order_vec = np.arange(all_orders.min(), all_orders.max() + 1, dtype=int)
    # invert the order vector so that it goes from the highest order to the lowest (i.e., from blue to red)
    order_vec = order_vec[::-1]

    parsed_sensobjs, parsed_sensnames = parse_sensfunc(order_vec, sensfuncs, sensnames, deckers, filters,
                                                       use_all=(not parse))
    # save to file
    parsed_sens_file = f'used_sensfuncs_v{sensfuncs[0].version}.yaml'
    with open(parsed_sens_file, 'w') as f:
        yaml.dump(parsed_sensnames, f)

    # get the spec1d file names corresponding to the parsed sensobjs
    use_spec1dfiles_list = []
    for parsed in parsed_sensobjs.values():
        for sensobj in parsed:
            use_spec1dfiles_list.append(Path(sensobj.spec1df).name)
    # save to file
    use_spec1dfiles_file = f'used_spec1ds_v{fits.getval(sensfuncs[0], "DMODVER", 1)}.txt'
    with open(use_spec1dfiles_file, 'w') as f:
        for fname in use_spec1dfiles_list:
            f.write(f"{fname}\n")

    return parsed_sensobjs, parsed_sensnames


def load_sensfunc(sens_files_path, sens_fnames=None, sens_fnames_dict=None,
                  only=None, parse=False, plot_all=False, ptype='counts_per_angs'):
    """Load the sensitivity functions from a list of files.

    Args:
        sens_files_path (str or Path):
            Path to the sensitivity function files.
        sens_fnames (list, optional):
            List of sensitivity function file names.
        sens_fnames_dict (dict, optional):
            Dictionary with the orders as keys and the
            sensitivity function file names as values.
        only (`numpy.ndarray`_):
            Array of files to use (raw filename used).
            If None, all the files in sens_files will be used.
        parse (bool):
            Parse the sensitivity functions. If True, all the sensitivity functions
            will be plotted per order and the user will be able to select the ones to use.
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

    sensfuncs = []
    sensnames = []
    deckers = []
    filters = []
    echangles = []
    xdangles = []
    orders = np.array([])

    for sens_file in sens_fnames:
        if only is not None and sens_file.split('-')[0].split('_')[1]+'.fits' not in only:
            continue
        sens_file = sens_files_path / sens_file
        if not sens_file.exists():
            msgs.warn(f'Sensitivity function file {sens_file} does not exist.')
        else:
            sensobj = SensFunc.from_file(sens_file, chk_version=False)
            sensfuncs.append(sensobj)
            orders = np.append(orders, sensobj.sens['ECH_ORDERS'].data)
            sensnames.append(sens_file.name)
            deckers.append(fits.getval(sens_file, 'DECKER'))
            filters.append(fits.getval(sens_file, 'FILTER1'))
            echangles.append(fits.getval(sens_file, 'ECHANGLE'))
            xdangles.append(fits.getval(sens_file, 'XDANGLE'))

    if len(sensfuncs) == 0:
        msgs.error("No sensitivity functions loaded.")

    # order_vec = np.arange(orders.min(), orders.max() + 1, dtype=int)
    order_vec = np.arange(orders.min(), 93 + 1, dtype=int)
    order_vec = order_vec[::-1]

    if sens_fnames_dict is None and parse:
        parsed_sensobjs, parsed_sensnames = parse_sensfunc(order_vec, sensfuncs, sensnames, deckers, filters,
                                                           use_all=False)
    elif sens_fnames_dict is None and not parse:
        parsed_sensobjs, parsed_sensnames = parse_sensfunc(order_vec, sensfuncs, sensnames, deckers, filters,
                                                           use_all=True)
    else:
        parsed_sensobjs = {}
        for iord, sens_fnames in sens_fnames_dict.items():
            parsed_sensobjs[iord] = [sensobj for sensobj in sensfuncs if
                                     Path(sensobj.spec1df).name.replace('spec1d', 'sens') in sens_fnames]
        parsed_sensnames = sens_fnames_dict

    # save to file
    parsed_sens_file = f'used_sensfuncs_v{sensfuncs[0].version}.yaml'
    _parsed_sensnames = {str(key): value for key, value in parsed_sensnames.items()}
    with open(parsed_sens_file, 'w') as f:
        yaml.dump(_parsed_sensnames, f)

    if plot_all:
        # colors by SensFunc object
        unique_colors_sens = get_unique_colors(len(sensfuncs), mcolors.CSS4_COLORS)
        color_map_sens = {sensobj: unique_colors_sens[i] for i, sensobj in enumerate(sensfuncs)}
        # color by order
        unique_colors_ord = get_unique_colors(len(order_vec), mcolors.CSS4_COLORS)
        color_map_ord = {ord: unique_colors_ord[i] for i, ord in enumerate(order_vec)}
        y_label = 'Zeropoint (AB mag)' if ptype == 'zeropoint' else 'Throughput' if ptype == 'throughput' \
            else 'Counts /(Angstrom s)'
        outfile_name = f'collection_{ptype}_by_order.pdf'
        y_all = []
        y_smooth_all = []
        y_scale_all = []
        y_maxmax_all = []
        wave_all = []
        order_all = []
        with PdfPages(outfile_name) as pdf:
            for iord in order_vec:
                i_sensfuncs = parsed_sensobjs[iord]
                if len(i_sensfuncs) == 0:
                    continue
                fig = plt.figure(figsize=(23, 6.))
                plt.minorticks_on()
                plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
                y_iord = []
                y_iord_smooth = []
                y_iord_smooth_max = []
                wave_iord = []
                legend_iord = []
                color_iord = []
                for s,sensobj in enumerate(i_sensfuncs):
                    _indx = np.where(sensobj.sens['ECH_ORDERS'].data == iord)[0]
                    if _indx.size > 0:
                        indx = _indx[0]
                        name = Path(sensobj.spec1df).name.split('_')[1]
                        airmass = sensobj.airmass
                        exptime = sensobj.exptime
                        color = color_map_sens[sensobj]
                        wave = sensobj.sens['SENS_WAVE'].data[indx]
                        wave_gmp = wave > 1.0
                        wmin = wave[wave_gmp].min()
                        wmax = wave[wave_gmp].max()
                        if ptype == 'zeropoint':
                            zeropoint_data = sensobj.sens['SENS_ZEROPOINT'].data[indx]
                            plt.plot(wave[wave_gmp], zeropoint_data[wave_gmp], alpha=0.4, color=color, lw=0.6,
                                     zorder=-2)
                            y = sensobj.sens['SENS_ZEROPOINT_FIT'].data[indx]
                        elif ptype == 'throughput':
                            wave = sensobj.wave[:,indx]
                            wcut = (wave >= wmin) & (wave <= wmax)
                            wave = wave[wcut]
                            wave_gmp = wave > 1.0
                            y = sensobj.throughput[:, indx][wcut]
                        else:
                            # counts per Angstrom
                            y = sensobj.sens['SENS_COUNTS_PER_ANG'].data[indx]/exptime

                        # smooth y
                        filt = 100
                        y_smooth = utils.fast_running_median(y[wave_gmp], filt)
                        y_smooth_max = np.nanmax(y_smooth[y_smooth < 10000])
                        # append by order
                        y_iord.append(y[wave_gmp])
                        y_iord_smooth.append(y_smooth)
                        y_iord_smooth_max.append(y_smooth_max)
                        wave_iord.append(wave[wave_gmp])
                        color_iord.append(color)
                        legend_iord.append(f'{name} - {deckers[s]} - {filters[s]} - ech: {echangles[s]} - xd: {xdangles[s]}\n'
                                           f'- airmass: {airmass:.2f} - exptime: {exptime:.1f}')

                        # append to all
                        y_all.append(y[wave_gmp])
                        y_smooth_all.append(y_smooth)
                        wave_all.append(wave[wave_gmp])
                        order_all.append(iord)
                if len(y_iord) == 0:
                    plt.close(fig)
                    continue
                # get the max value of the smoothed y
                y_maxmax = np.nanmax(y_iord_smooth_max)
                # append to all
                y_maxmax_all.append(y_maxmax)

                # plot by order
                for y, y_sm, y_max, wave, color, legend in zip(y_iord, y_iord_smooth, y_iord_smooth_max, wave_iord, color_iord, legend_iord):
                    s = y_maxmax/y_max
                    plt.plot(wave, y*s, color=color, alpha=0.6, lw=0.6, zorder=0, label=legend)
                    plt.plot(wave, y_sm*s, color=color, lw=1.5, zorder=1)
                    y_scale_all.append(s)
                plt.title(f'Order {iord}')
                plt.xlabel('Wavelength (Angstroms)')
                plt.ylabel(y_label)
                if ptype == 'zeropoint':
                    plt.ylim(13.1, 20.9)
                elif ptype == 'throughput':
                    plt.ylim(0.0, 0.15)
                else:
                    plt.ylim(0.0, y_maxmax*1.2)
                plt.legend(fontsize=5)
                fig.tight_layout()
                pdf.savefig(dpi=60)
                plt.close(fig)
            # plot all orders
            fig = plt.figure(figsize=(23, 6.))
            plt.minorticks_on()
            plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
            plt.title(f'All orders')
            for y, y_sm, s, wave, o in zip(y_all, y_smooth_all, y_scale_all, wave_all, order_all):
                color = color_map_ord[o]
                leg = f'order {o}' if o not in [line.get_label() for line in plt.gca().get_lines()] else ''
                plt.plot(wave, y*s, lw=0.6, color=color, alpha=0.6, zorder=0, label=leg)
                plt.plot(wave, y_sm*s, lw=1.5, color=color, zorder=1)

            plt.xlabel('Wavelength (Angstroms)')
            plt.ylabel(y_label)
            plt.legend(fontsize=5)
            if ptype == 'zeropoint':
                plt.ylim(13.1, 20.9)
            elif ptype == 'throughput':
                plt.ylim(0.0, 0.15)
            else:
                plt.ylim(0.0, np.nanmax(y_maxmax_all)*1.2)
            fig.tight_layout()
            pdf.savefig(dpi=60)
            plt.close(fig)
    #embed()
    return parsed_sensobjs, parsed_sensnames


def parse_sensfunc(ordervec, sensfuncs, sensnames, deckers, filters, use_all=True, rescale=False):
    """Parse the sensitivity functions.

    Args:
        ordervec (list or `numpy.ndarray`_):
            List of orders.
        sensfuncs (list):
            List of SensFunc objects.
        sensnames (list):
            List of sensitivity function file names.
        deckers (list):
            List of deckers of the SensFunc objects.
        filters (list):
            List of filters of the SensFunc objects.
        use_all (bool):
            Use all the sensitivity functions. If False,
            the user will be able to select the ones to use in each order.
        rescale (bool):
            Rescale the zeropoints to the same median value. Only when use_all is False.

    Returns:
        dict: Dictionary with the selected SensFunc objects.

    """

    if rescale:
        zpoint_fit_list = []
        for sensobj in sensfuncs:
            wv_gpm = sensobj.sens['SENS_WAVE'].data > 1.0
            zpoint_fit_list.append(sensobj.sens['SENS_ZEROPOINT_FIT'].data[wv_gpm])

        # get the rescaling factor
        medians = [stats.sigma_clipped_stats(zp, mask_value=0., sigma=3)[1] for zp in zpoint_fit_list]
        max_median = np.nanmax(medians)
        scales = np.ones(len(medians)) * max_median / medians
    else:
        scales = np.ones(len(sensfuncs))

    # colors
    unique_colors = get_unique_colors(len(sensfuncs), mcolors.CSS4_COLORS)
    color_map = {sensobj: unique_colors[i] for i, sensobj in enumerate(sensfuncs)}

    selected_sensfuncs = {}
    selected_lines = {}
    for iord in ordervec:
        if use_all:
            selected_sensfuncs[iord] = sensfuncs
            selected_lines[iord] = sensnames
        else:
            selected_lines[iord] = []
            fig = plt.figure(figsize=(23, 6.))
            ax = plt.subplot()
            ax.set_title(f'Order {iord}')
            ax.set_xlabel('Wavelength (Angstroms)')
            ax.set_ylabel('Zeropoint (AB mag)')
            ax.minorticks_on()
            ax.tick_params(axis='both', direction='in', top=True, right=True, which='both')
            for sensobj, sname, scale, decker, filt in zip(sensfuncs, sensnames, scales, deckers, filters):
                _indx = np.where(sensobj.sens['ECH_ORDERS'].data == iord)[0]
                if _indx.size > 0:
                    indx = _indx[0]
                    airmass = sensobj.airmass
                    wave = sensobj.sens['SENS_WAVE'].data[indx]
                    wave_gmp = wave > 1.0
                    y = sensobj.sens['SENS_ZEROPOINT_FIT'][indx].data * scale
                    color = color_map[sensobj]
                    ax.plot(wave[wave_gmp], y[wave_gmp], color=color, lw=3., zorder=0,
                            label=f'{sname} - {decker} - {filt} -airmass: {airmass:.2f}', picker=True)
            ax.legend(fontsize=5)
            ax.set_ylim(13.1, 20.9)
            fig.tight_layout()
            fig.canvas.mpl_connect('pick_event', partial(on_pick, selected_lines=selected_lines, key=iord))
            print(' ')
            plt.show()
            selected_sensfuncs[iord] = [sensobj for sensobj in sensfuncs if
                                        Path(sensobj.spec1df).name.replace('spec1d', 'sens') in selected_lines[iord]]

    return selected_sensfuncs, selected_lines


# def get_std_dict(sensfuncs):
#     """Get the standard star dictionary from the sensitivity functions.
#
#     Args:
#         sensfuncs (list):
#             List of SensFunc objects.
#
#     Returns:
#         dict: Dictionary with the standard star information.
#     """
#
#     cal_names = []
#     cal_ras = []
#     cal_decs = []
#     for sensobj in sensfuncs:
#         if sensobj.std_name not in cal_names:
#             cal_names.append(sensobj.std_name)
#             cal_ras.append(sensobj.std_ra)
#             cal_decs.append(sensobj.std_dec)
#
#     all_std_dicts = {}
#     for i, cal_name in enumerate(cal_names):
#         std_dict = flux_calib.get_standard_spectrum(ra=cal_ras[i], dec=cal_decs[i])
#         all_std_dicts[cal_name] = std_dict
#
#     return all_std_dicts


def fill_combSensObj(comb_zeropoint_fit_list, comb_wave_list, comb_orders_list, poly_order, comb_coeff_list, sensobj):
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
        sensobj (SensFunc):
            SensFunc object of one of the sensitivity functions used to compute the combined zeropoints.
            This object will be used to initialize the combined SensFunc object.

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

    # initialize a SensFunc object
    comb_sensobj = sensobj
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


def combine_sensfuncs(sensfuncs_dict, sensnames_dict, cut_left=50, cut_right=-100, rescale=True,
                      poly_order=4, qa_plots=False):
    """Combine the sensitivity functions.

    Args:
        sensfuncs_dict (dict):
            Dictionary with the selected SensFunc objects per ech order.
        sensnames_dict (dict):
            Dictionary with the selected sensitivity function file names per ech order.
        cut_left (int):
            Number of pixels to cut from the left.
        cut_right (int):
            Number of pixels to cut from the right.
        rescale (bool):
            Rescale the zeropoints to the same median value.
        poly_order (int):
            Order of the polynomial to fit the zeropoints.
        qa_plots (bool):
            Generate QA plots.

    Returns:
        SensFunc: Combined sensitivity function object.
    """

    sensfuncs = list({s for sf in sensfuncs_dict.values() for s in sf})
    sens_orders = np.sort(list(sensfuncs_dict.keys()))[::-1]
    if np.any(np.diff(sens_orders) != -1):
        msgs.warns("Orders are not contiguous. There might be missing orders.")

    # if rescale:
    #     medians = []
    #     for sensobj in sensfuncs:
    #         zpoint_poly = np.array([])
    #         gpm = np.array([])
    #         for i,sens in enumerate(sensobj.sens):
    #             wave = sens['SENS_WAVE']
    #             wv_gpm = wave > 1.0
    #             wave = wave[wv_gpm][cut_left:cut_right]
    #             zpoint_data = sens['SENS_ZEROPOINT'].data[wv_gpm][cut_left:cut_right]
    #             zpoint_data_gpm = sens['SENS_ZEROPOINT_GPM'].data[wv_gpm][cut_left:cut_right]
    #             log10blaze = sens['SENS_LOG10_BLAZE_FUNCTION'].data[wv_gpm][cut_left:cut_right]
    #             # get zeropoint_poly
    #             _zpoint_poly = zpoint_data + 5.0*np.log10(wave) - ZP_UNIT_CONST
    #             if not np.all(log10blaze==0):
    #                 _zpoint_poly -= 2.5*log10blaze
    #             zpoint_poly = np.concatenate([zpoint_poly, _zpoint_poly])
    #             gpm = np.concatenate([gpm, zpoint_data_gpm])
    #
    #         medians.append(stats.sigma_clipped_stats(zpoint_poly,
    #                                                  mask=np.isinf(zpoint_poly) | np.isnan(zpoint_poly) | np.logical_not(gpm.astype(bool)),
    #                                                  mask_value=0., sigma=3)[1])
    #     max_median = np.nanmax(medians)
    #     scales = np.ones(len(medians)) * max_median / medians
    # else:
    #     scales = np.ones(len(sensfuncs))

    # initialize the collection of zeropoints as a dictionary
    zpoints_all = {}

    for iord in sens_orders:
        # get the SensFunc objects that will be used for the current order
        sensobjs = [sfs for sfs in sensfuncs_dict[iord]]
        sensnames = [sns for sns in sensnames_dict[iord]]
        # get the zeropoints for the current order
        zpoints_iord = {}
        zpoints_iord['name'] = []
        zpoints_iord['wave'] = []
        zpoints_iord['zpoint_data'] = []
        zpoints_iord['zpoint_data_gpm'] = []
        zpoints_iord['log10blaze'] = []
        zpoints_iord['zpoint_poly'] = []
        zpoints_iord['median'] = []
        zpoints_iord['max_med'] = []
        for sensobj, sensname in zip(sensobjs, sensnames):
            _indx = np.where(sensobj.sens['ECH_ORDERS'].data == iord)[0]
            if np.any(_indx):
                indx = _indx[0]
                wave = sensobj.sens['SENS_WAVE'].data[indx]
                wave_gpm = wave > 1.0
                wave = wave[wave_gpm][cut_left:cut_right]
                zpoint_data = sensobj.sens['SENS_ZEROPOINT'].data[indx][wave_gpm][cut_left:cut_right]
                zpoint_fit = sensobj.sens['SENS_ZEROPOINT_FIT'].data[indx][wave_gpm][cut_left:cut_right]
                zpoint_data_gpm = sensobj.sens['SENS_ZEROPOINT_GPM'].data[indx][wave_gpm][cut_left:cut_right]
                log10blaze = sensobj.sens['SENS_LOG10_BLAZE_FUNCTION'].data[indx][wave_gpm][cut_left:cut_right]
                zpoint_coeff = sensobj.sens['SENS_COEFF'].data[indx]
                # # get zeropoint_poly
                # zpoint_poly = zpoint_data + 5.0*np.log10(wave) - ZP_UNIT_CONST
                # if not np.all(log10blaze==0):
                #     zpoint_poly -= 2.5*log10blaze
                # this is actually the zeropoint
                # zpoint_poly = flux_calib.eval_zeropoint(zpoint_coeff, 'legendre', wave, wave.min(), wave.max())
                zpoint_poly = zpoint_fit
                med = stats.sigma_clipped_stats(zpoint_poly,
                                                mask=np.isinf(zpoint_poly) | np.isnan(zpoint_poly) | np.logical_not(zpoint_data_gpm.astype(bool)),
                                                mask_value=0., sigma=3)[1]
                zpoints_iord['name'].append(sensname)
                zpoints_iord['wave'].append(wave)
                zpoints_iord['zpoint_data'].append(zpoint_data)
                zpoints_iord['zpoint_data_gpm'].append(zpoint_data_gpm)
                zpoints_iord['log10blaze'].append(log10blaze)
                zpoints_iord['zpoint_poly'].append(zpoint_poly)
                zpoints_iord['median'].append(med)
        max_median = np.nanmax(zpoints_iord['median']) if len(zpoints_iord['median']) > 0 else []
        zpoints_iord['max_med'].append(max_median)

        zpoints_all[iord] = zpoints_iord

    # get the wavelength grid from the telluric grid
    hdul = io.load_telluric_grid(sensfuncs[0].telluric.telgrid)
    wave_grid_tell = hdul[1].data
    # keep only the wavelengths that are within the range of the zeropoints
    pad_frac = 0.1
    wmin_grid = (1.0 - pad_frac)*np.min([z.min() for zpoints in zpoints_all.values() if len(zpoints['wave']) > 0 for z in zpoints['wave']])
    wmax_grid = (1.0 + pad_frac)*np.max([z.max() for zpoints in zpoints_all.values() if len(zpoints['wave']) > 0 for z in zpoints['wave']])
    wave_grid = wave_grid_tell[(wave_grid_tell > wmin_grid) & (wave_grid_tell < wmax_grid)]

    # initialize the combined zeropoints
    comb_zeropoint_fit_list = []
    comb_wave_list = []
    comb_orders_list = []
    comb_coeff_list = []
    for iord in sens_orders:
        if len(zpoints_all[iord]['wave']) == 0:
            continue
        zpoints_iord = zpoints_all[iord]
        wave_iord = np.concatenate(zpoints_iord['wave'])
        medians_iord = np.array(zpoints_iord['median'])
        max_median_iord = np.array(zpoints_iord['max_med'])
        # rescale the zeropoints
        scale_iord = max_median_iord / medians_iord
        zpoint_poly_iord = np.concatenate([z * s for z,s in zip(zpoints_iord['zpoint_poly'], scale_iord)])
        gpm_iord = np.concatenate(zpoints_iord['zpoint_data_gpm'])
        func = 'legendre'
        pypeitFit = fitting.robust_fit(wave_iord, zpoint_poly_iord, poly_order,
                                       function=func, minx=wave_iord.min(), maxx=wave_iord.max(),
                                       in_gpm=gpm_iord, lower=10, upper=10, use_mad=True)
        # comb_zeropoint_fit_iord = flux_calib.eval_zeropoint(pypeitFit.fitc, func, wave_grid,
        #                                                     wave_iord.min(), wave_iord.max())
        comb_zeropoint_fit_iord = fitting.evaluate_fit(pypeitFit.fitc, func, wave_grid, minx=wave_iord.min(),
                                                        maxx=wave_iord.max())
        comb_wave = wave_grid.copy()
        outside_gpm = (comb_wave < wave_iord.min()) | (comb_wave > wave_iord.max())
        comb_zeropoint_fit_iord[outside_gpm] = 0.
        comb_wave[outside_gpm] = 0.

        # plt.plot(wave_iord[wave_iord > 1], zpoint_poly_iord[wave_iord > 1], marker='.', ms=4, ls='')
        # plt.plot(comb_wave[comb_wave > 0], comb_zeropoint_fit_iord[comb_wave > 0], '--r')
        # plt.show()
        #embed()
        # append the results to the lists
        comb_zeropoint_fit_list.append(comb_zeropoint_fit_iord)
        comb_wave_list.append(comb_wave)
        comb_orders_list.append(iord)
        comb_coeff_list.append(pypeitFit.fitc)
        # save also the scale
        zpoints_all[iord]['scale'] = scale_iord

    comb_sensobj = fill_combSensObj(comb_zeropoint_fit_list, comb_wave_list, comb_orders_list, poly_order,
                                    comb_coeff_list, deepcopy(sensfuncs[0]))

    # QA plots
    if qa_plots:
        plot_by_order(comb_sensobj, zpoints_all)


    if qa_plots:
        # plot the zeropoints for all orders
        plot_all_orders(comb_sensobj, comb_type='zeropoint')
        # plot the throughput for all orders
        plot_all_orders(comb_sensobj,comb_type='throughput')

    return comb_sensobj


def plot_by_order(comb_sensobj, zpoints_all):
    """Plot the zeropoints by order.

    Args:
        comb_sensobj (SensFunc):
            Combined sensitivity function object.
        zpoints_all
    """

    unique_sensfiles = np.unique([z for zpoints in zpoints_all.values() if len(zpoints['name']) > 0 for z in zpoints['name']])
    unique_colors = get_unique_colors(unique_sensfiles.size, mcolors.CSS4_COLORS)
    color_map = {sensobj: unique_colors[i] for i, sensobj in enumerate(unique_sensfiles)}

    outfile = 'collection_RESCALEDzeropoints_by_order_with_comb2.pdf'

    # get the values to plot
    orders = comb_sensobj.sens['ECH_ORDERS'].data
    waves = comb_sensobj.wave.T
    comb_zeropoint_fit = comb_sensobj.zeropoint.T # this is without the blaze correction
    comb_coeff = comb_sensobj.sens['SENS_COEFF'].data

    with PdfPages(outfile) as pdf:
        for i, iord in enumerate(orders):
            if len(zpoints_all[iord]['wave']) == 0:
                continue
            zpoints_iord = zpoints_all[iord]
            waves_iord = zpoints_iord['wave']
            zeropoints_data_iord = zpoints_iord['zpoint_data']
            zeropoints_fit_iord = zpoints_iord['zpoint_poly']
            scales_iord = zpoints_iord['scale']
            logblaze_iord = zpoints_iord['log10blaze']
            sensfiles_iord = zpoints_iord['name']
            # combined zpoints
            comb_zeropoint = comb_zeropoint_fit[i]
            wave = waves[i]
            coeff = comb_coeff[i]

            fig = plt.figure(figsize=(23, 6.))
            axis = plt.subplot()
            plt.minorticks_on()
            plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
            plt.title(f'Order {iord}')
            for i in range(len(sensfiles_iord)):
                label = f"{sensfiles_iord[i]}"
                color = color_map[sensfiles_iord[i]]

                axis.plot(waves_iord[i], zeropoints_data_iord[i]*scales_iord[i], drawstyle='steps-mid',
                          color=color, ls='', marker='.', ms=2, alpha=0.2, zorder=-2)
                # axis.plot(waves_iord[i], zeropoints_fit_iord[i], color=color, linewidth=1., alpha=1.,
                #           label=label, zorder=0)
                # recover the original zeropoints using the combined coefficients
                # poly_model = fitting.evaluate_fit(coeff, 'legendre', waves_iord[i],
                #                                   minx=waves_iord[i].min(), maxx=waves_iord[i].max())
                # # correct for scale
                # poly_model /= scales_iord[i]
                # this_comb = poly_model - 5.0 * np.log10(waves_iord[i]) + ZP_UNIT_CONST
                # if np.all(logblaze_iord[i]!=0):
                #     this_comb += 2.5 * logblaze_iord[i]
                # comb_zpoint = flux_calib.eval_zeropoint(coeff, 'legendre', waves_iord[i],
                #                                         waves_iord[i].min(), waves_iord[i].max(),
                #                                         log10_blaze_func_per_ang=logblaze_iord[i])
                if wave[wave>0].size > 0:
                    comb_zpoint = scipy.interpolate.interp1d(wave[wave>0], comb_zeropoint[wave>0],
                                                             bounds_error=False, fill_value=0.)(waves_iord[i])
                    axis.plot(waves_iord[i], comb_zpoint, color=color, linewidth=1.5, alpha=1.,
                              label=label, zorder=0)
                    axis.plot(waves_iord[i], comb_zpoint, color='k', linewidth=1.5, alpha=1.,zorder=1)
            if iord in [93, 35, 36, 37]:
                axis.set_ylim(13.1, 17.1)
            else:
                axis.set_ylim(16.5, 20.5)
            # axis.set_xlim(_wmin, _wmax)
            axis.legend(fontsize=5)
            axis.set_xlabel('Wavelength (Angstroms)')
            axis.set_ylabel('Zeropoint (AB mag)')
            fig.tight_layout()
            pdf.savefig(dpi=60)
            # plt.show()
            plt.close(fig)


def plot_all_orders(comb_sensobj, comb_type='zeropoint'):
    """Plot the combined zeropoint or throughput for all orders.

    Args:
        comb_sensobj (SensFunc):
            Combined sensitivity function object.
        comb_type (str):
            Type of data to plot. It can be 'zeropoint' or 'throughput'.

    """
    # get the values to plot
    orders = comb_sensobj.sens['ECH_ORDERS'].data
    waves = comb_sensobj.wave.T
    combined = comb_sensobj.zeropoint.T if comb_type == 'zeropoint' else comb_sensobj.throughput.T

    unique_colors = get_unique_colors(len(orders), mcolors.CSS4_COLORS)
    pcolors = {sensobj: unique_colors[i] for i, sensobj in enumerate(orders)}
    outfile = f'Combined_{comb_type}_all_orders.pdf'

    fig = plt.figure(figsize=(23, 6.))
    axis = plt.subplot()
    plt.minorticks_on()
    plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
    plt.title(f'PypeIt {comb_type} for keck_hires - orders {orders[0]} to {orders[-1]}')
    for (order, wave, comb) in zip(orders, waves, combined):
        color = pcolors[order]
        axis.plot(wave[wave>1], comb[wave>1], color=color, linewidth=1, alpha=1., label=f'order {order}', zorder=1)
    axis.set_ylim(13.1, 20.9) if comb_type == 'zeropoint' else axis.set_ylim(0.0, 0.15)
    # axis.set_ylim(10.1, 20.9)
    axis.set_xlim(3600., 10400.)
    axis.legend(loc='lower left', fontsize=4)
    axis.set_xlabel('Wavelength (Angstroms)')
    axis.set_ylabel('Zeropoint (AB mag)') if comb_type == 'zeropoint' else axis.set_ylabel('Throughput')
    fig.tight_layout()
    plt.savefig(outfile, dpi=60)
    # plt.show()
    plt.close(fig)


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
    parser.add_argument("--parse", action="store_true", default=False,
                        help="Parse the sensitivity functions to select the ones to use for each order.")
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
        sensfuncs_dict, sens_fnames_dict = load_sensfunc(sensfuncs_path, sens_fnames=sens_fnames,
                                                         sens_fnames_dict=sens_fnames_dict,
                                                         parse=args.parse, plot_all=args.QA)

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
        # create the sensitivity functions
        sensfuncs_dict, sens_fnames_dict = create_sens_files(spec1d_fnames, spec1ds_path, sensfuncs_path,
                                                             boxcar=args.boxcar, use_flat=args.use_flat,
                                                             skip_existing=args.skip_existing, parse=args.parse)
        sens_fnames = np.unique([item for sublist in sens_fnames_dict.values() for item in sublist if sublist])

    if len(sens_fnames) == 0:
        msgs.error(f"No sensitivity functions found.")


    # for debugging
    # tab_all = print_datasets(sensfuncs)

    # combine the sensitivity functions
    sensfunc = combine_sensfuncs(sensfuncs_dict, sens_fnames_dict, poly_order=4, qa_plots=args.QA)

    # save the combined sensitivity function
    # Construct a primary FITS header
    with io.fits_open(sensfuncs_path / sens_fnames[0]) as hdul:
        spectrograph = load_spectrograph(hdul[0].header['PYP_SPEC'])
        primary_hdr = io.initialize_header()
        add_keys = (['PYP_SPEC', 'DATE-OBS', 'TELESCOP', 'INSTRUME', 'DETECTOR'] +
                    spectrograph.configuration_keys() + spectrograph.raw_header_cards())
        for key in add_keys:
            if key.upper() in hdul[0].header.keys():
                primary_hdr[key.upper()] = hdul[0].header[key.upper()]
    sensfunc.to_file(args.outfile, primary_hdr=primary_hdr, overwrite=True)


if __name__ == '__main__':
    args = parse_args()
    args.outfile = 'sens_keck_hires_RED_orders_93-35_sensfunc.fits'
    main(args)


