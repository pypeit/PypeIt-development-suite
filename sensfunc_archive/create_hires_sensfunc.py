
import argparse
import scipy
from copy import deepcopy
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
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
from pypeit import dataPaths

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
###########

def select_datasets(sensfuncs):
    """ ONLY FOR DEBUGGING. GET info on the sensfuncs and select the ones to use for the combined sensfunc.

    Args:
        sensfuncs (list):
            List of SensFunc objects.

    Returns:
        tuple: Table with the information of the SensFunc objects and a boolean array to select the objects to use.

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
    keep = only_these_datasets()
    aa = np.isin(tab_all['name'].data, keep)
    tab_all[aa].pprint_all()
    return tab_all, aa


def create_sens_files(spec1d_files, spec1d_files_path, sens_files_path, boxcar=False,
                      use_flat=True, skip_existing=True, only=None):
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

    Returns:
        tuple: List of SensFunc objects and list of sensitivity function file names.

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
    use_sensfiles_list = []
    use_spec1dfiles_list = []
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
            use_spec1dfiles_list.append(spec1d_file)
            use_sensfiles_list.append(sens_file.name)
        except Exception as e:
            msgs.warn(f'Error computing sensitivity function for {spec1d_file}.\n{e}')

    if len(sensfuncs) > 0:
        # save the sensfunc and spec1d files that were used to compute the combined zeropoints to a text file and copy them
        # to a new directory
        sensfunc_basename = f'used_sensfuncs_v{sensobj.version}'
        spec1d_basename = f'used_spec1ds_v{fits.getval(sensfuncs[0], "DMODVER", 1)}'
        # Save the file paths to text files
        with open(f'{sensfunc_basename}.txt', 'w') as f:
            for fname in use_sensfiles_list:
                f.write(f"{fname}\n")

        with open(f'{spec1d_basename}.txt', 'w') as f:
            for fname in use_spec1dfiles_list:
                f.write(f"{fname}\n")

    return sensfuncs, use_sensfiles_list


def load_sensfunc(sens_files, sens_files_path, only=None, plot_all=False, ptype='zeropoint'):
    """Load the sensitivity functions from a list of files.

    Args:
        sens_files (list):
            List of sensitivity function file names.
        sens_files_path (str or Path):
            Path to the sensitivity function files.
        only (`numpy.ndarray`_):
            Array of files to use (raw filename used).
            If None, all the files in sens_files will be used.
        plot_all (bool):
            Plot the sensitivity functions.
        ptype (str):
            Type of plot to generate. It can be 'zeropoint', 'throughput', or 'counts_per_angs'.

    Returns:
        list: List of SensFunc objects.
    """

    # check if sens_files_path is a Path object
    if not isinstance(sens_files_path, Path):
        sens_files_path = Path(sens_files_path).absolute()
    # check if the path exists
    if not sens_files_path.exists():
        msgs.error(f"Sensitivity functions path {sens_files_path} does not exist.")

    sensfuncs = []
    deckers = []
    filters = []

    for sens_file in sens_files:
        if only is not None and sens_file.split('-')[0].split('_')[1]+'.fits' not in only:
            continue
        sens_file = sens_files_path / sens_file
        if not sens_file.exists():
            msgs.warn(f'Sensitivity function file {sens_file} does not exist.')
        else:
            sensfuncs.append(SensFunc.from_file(sens_file, chk_version=False))
            deckers.append(fits.getval(sens_file, 'DECKER'))
            filters.append(fits.getval(sens_file, 'FILTER1'))


    if len(sensfuncs) == 0:
        msgs.error("No sensitivity functions loaded.")

    if plot_all:
        unique_colors = get_unique_colors(len(sensfuncs), mcolors.CSS4_COLORS)
        color_map = {sensobj: unique_colors[i] for i, sensobj in enumerate(sensfuncs)}
        order_vec = np.arange(35, 103, dtype=int)
        outfile_name = f'collection_{ptype}_by_order.pdf'
        y_list = []
        wave_list = []
        with PdfPages(outfile_name) as pdf:
            for o in order_vec[::-1]:
                fig = plt.figure(figsize=(23, 6.))
                plt.minorticks_on()
                plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
                for s,sensobj in enumerate(sensfuncs):
                    indx = sensobj.sens['ECH_ORDERS'].data == o
                    if np.any(indx):
                        name = Path(sensobj.spec1df).name.split('-')[0].split('_')[1]
                        airmass = sensobj.airmass
                        exptime = sensobj.exptime
                        wave = sensobj.sens['SENS_WAVE'][indx]
                        wave_gmp = wave > 1.0
                        wmin = wave[wave_gmp].min()
                        wmax = wave[wave_gmp].max()
                        if ptype == 'zeropoint':
                            zeropoint_data = sensobj.sens['SENS_ZEROPOINT'][indx]
                            y = sensobj.sens['SENS_ZEROPOINT_FIT'][indx]
                            y_label = 'Zeropoint (AB mag)'
                        elif ptype == 'throughput':
                            y = sensobj.throughput[:,indx]
                            wave = sensobj.wave[:,indx]
                            wave_gmp = (wave >= wmin) & (wave <= wmax)
                            y_label = 'Throughput'
                        elif ptype == 'counts_per_angs':
                            y = sensobj.sens['SENS_COUNTS_PER_ANG'][indx]/exptime
                            y_label = 'Counts /(Angstrom s)'
                        color = color_map[sensobj]
                        if ptype == 'zeropoint':
                            plt.plot(wave[wave_gmp], zeropoint_data[wave_gmp], alpha=0.4, color=color, lw=0.6, zorder=-2)
                        plt.plot(wave[wave_gmp], y[wave_gmp], color=color, lw=0.6, zorder=0,
                                 label=f'{name} - {deckers[s]} - {filters[s]} -airmass: {airmass:.2f}')
                        y_list.append(y[wave_gmp])
                        wave_list.append(wave[wave_gmp])
                if np.sum(plt.axis()[:2]) < 10.:
                    plt.close(fig)
                    continue
                plt.title(f'Order {o}')
                plt.xlabel('Wavelength (Angstroms)')
                plt.ylabel(y_label)
                if ptype == 'zeropoint':
                    plt.ylim(13.1, 20.9)
                elif ptype == 'throughput':
                    plt.ylim(0.0, 0.15)
                plt.legend()
                fig.tight_layout()
                pdf.savefig(dpi=60)
                plt.close(fig)
            fig = plt.figure(figsize=(23, 6.))
            plt.minorticks_on()
            plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
            plt.title(f'All orders')
            for y, wave in zip(y_list, wave_list):
                plt.plot(wave, y, lw=0.6, color='tab:blue')
            plt.xlabel('Wavelength (Angstroms)')
            plt.ylabel(y_label)
            if ptype == 'zeropoint':
                plt.ylim(13.1, 20.9)
            elif ptype == 'throughput':
                plt.ylim(0.0, 0.15)
            fig.tight_layout()
            pdf.savefig(dpi=60)
            plt.close(fig)

    return sensfuncs

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
    comb_sensobj.sens['WAVE_MIN'] = np.array([np.min(w[w>1]) for w in comb_wave_array])
    comb_sensobj.sens['WAVE_MAX'] = np.array([np.max(w[w>1]) for w in comb_wave_array])
    comb_sensobj.sens['POLYORDER_VEC'] = poly_order

    comb_sensobj.wave = comb_wave_array.T
    comb_sensobj.zeropoint = comb_zeropoint_fit_array.T

    comb_sensobj.spectrograph = load_spectrograph("keck_hires")
    comb_sensobj.throughput = comb_sensobj.compute_throughput()[0]

    return comb_sensobj

def combine_sensfuncs(sensfuncs, cut_left=50, cut_right=-100, poly_order=4, qa_plots=False):
    """Combine the sensitivity functions.

    Args:
        sensfuncs (list):
            List of SensFunc objects.
        cut_left (int):
            Number of pixels to cut from the left.
        cut_right (int):
            Number of pixels to cut from the right.
        poly_order (int):
            Order of the polynomial to fit the zeropoints.
        qa_plots (bool):
            Generate QA plots.

    Returns:
        SensFunc: Combined sensitivity function object.
    """

    all_orders = np.array([])
    counts_per_angs_list = []
    for sensobj in sensfuncs:
        all_orders = np.append(all_orders, sensobj.sens['ECH_ORDERS'].data)
        wv_gpm = sensobj.sens['SENS_WAVE'] > 1.0
        counts_per_angs_list.append(sensobj.sens['SENS_COUNTS_PER_ANG'].data[wv_gpm][cut_left:cut_right]/sensobj.exptime)

    # generate the orders vector for the combined sensitivity function
    order_vec = np.arange(all_orders.min(), all_orders.max() + 1, dtype=int)
    # invert the order vector so that it goes from the highest order to the lowest (i.e., from blue to red)
    order_vec = order_vec[::-1]

    # rescale the counts_per_angs to the same median value
    medians = [stats.sigma_clipped_stats(cpa, mask_value=0., sigma=3)[1] for cpa in counts_per_angs_list]
    max_median = np.nanmax(medians)
    scales = np.ones(len(medians)) * max_median / medians

    # load all the calibration standard stars
    all_std_dicts = get_std_dict(sensfuncs)

    # initialize the combined zeropoints
    comb_zeropoint_fit_list = []
    comb_wave_list = []
    comb_orders_list = []
    comb_coeff_list = []
    scale_list = []
    # list of the zeropoints for each order
    zeropoints_fit_list = []
    zeropoints_data_list = []
    zeropoint_poly_list = []
    logblaze_list = []
    waves_list = []
    sensfiles_iord_list = []
    gpm_list = []

    # loop over the orders
    for iord in order_vec:
        zeropoint_medians = []
        waves_iord, zeropoints_data_iord, zeropoints_fit_iord, zeropoint_poly_iord, gpms_iord, sensfiles_iord, log10blaze_iord = \
            [], [], [], [], [], [], []
        for _s, sensobj in enumerate(sensfuncs):
            for sens in sensobj.sens:
                if sens['ECH_ORDERS'] == iord:
                    wave = sens['SENS_WAVE']
                    wv_gpm = wave > 1.0
                    wave = wave[wv_gpm][cut_left:cut_right]
                    counts_per_angs =scales[_s] * sens['SENS_COUNTS_PER_ANG'][wv_gpm][cut_left:cut_right]/sensobj.exptime
                    std_dict = all_std_dicts[sensobj.std_name]
                    flam_true = scipy.interpolate.interp1d(std_dict['wave'].value, std_dict['flux'].value, kind='linear',
                                                           bounds_error=False, fill_value=-1e20)(wave)
                    flam_true_gpm = (wave >= std_dict['wave'].value.min()) & (wave <= std_dict['wave'].value.max())

                    zeropoint_data, zeropoint_data_gpm = flux_calib.compute_zeropoint(wave, counts_per_angs,
                                                                                      flam_true_gpm, flam_true)
                    log10blaze = sens['SENS_LOG10_BLAZE_FUNCTION'][wv_gpm][cut_left:cut_right]
                    # get zeropoint_poly
                    zeropoint_poly = zeropoint_data + 5.0*np.log10(wave) - ZP_UNIT_CONST
                    if not np.all(log10blaze==0):
                        zeropoint_poly -= 2.5*log10blaze



                    # zeropoint_data2 = sens['SENS_ZEROPOINT'][wv_gpm][cut_left:cut_right]
                    # zeropoint_data_gpm = sens['SENS_ZEROPOINT_GPM'][wv_gpm][cut_left:cut_right]
                    zeropoint_fit = sens['SENS_ZEROPOINT_FIT'][wv_gpm][cut_left:cut_right]
                    # zeropoint_fit_gpm = sens['SENS_ZEROPOINT_FIT_GPM'][wv_gpm][cut_left:cut_right]


                    # compute median of the zeropoints
                    # zeropoint_medians.append(stats.sigma_clipped_stats(zeropoint_poly, mask_value=0., sigma=3)[1])
                    waves_iord.append(wave)
                    gpms_iord.append(zeropoint_data_gpm)
                    zeropoints_fit_iord.append(zeropoint_fit)
                    zeropoints_data_iord.append(zeropoint_data)
                    zeropoint_poly_iord.append(zeropoint_poly)
                    log10blaze_iord.append(log10blaze)
                    sensfiles_iord.append(f"{Path(sensobj.spec1df).name.replace('spec1d', 'sens')}")

        # # scale the zeropoints to the same median value (max of all the medians)
        # zeropoint_medians = np.array(zeropoint_medians)
        # max_median = np.nanmax(zeropoint_medians)
        # scales = np.ones(len(zeropoint_medians)) * max_median / zeropoint_medians
        # zeropoints_fit_scaled_iord = [z * s for z, s in zip(zeropoints_fit_iord, scales)]
        # zeropoints_data_scaled_iord = [z * s for z, s in zip(zeropoints_data_iord, scales)]
        # zeropoint_poly_scaled_iord = [z * s for z, s in zip(zeropoint_poly_iord, scales)]

        # discard the zeropoints that have values below a certain threshold. This allows to discard zeropoints
        # that go down because of the blocking filter or other reasons.
        # lower limit for the zeropoints
        low_thresh = 7.
        # if iord == 35:
        #     low_thresh = 13.
        # elif iord in [43,66,84]:
        #     low_thresh = 18.
        # elif iord == 88:
        #     low_thresh = 17.
        # elif iord == 76:
        #     low_thresh = 17.5
        # ind_to_remove = []
        # for i in range(len(zeropoints_fit_scaled_iord)):
        #     if np.any(zeropoints_fit_scaled_iord[i] < low_thresh) or np.any(zeropoints_fit_scaled_iord[i] > 20.):
        #         ind_to_remove.append(i)
        # if len(ind_to_remove) > 0:
        #     zeropoints_data_scaled_iord = [zeropoints_data_scaled_iord[k] for k in range(len(zeropoints_data_scaled_iord)) if k not in ind_to_remove]
        #     zeropoints_fit_scaled_iord = [zeropoints_fit_scaled_iord[k] for k in range(len(zeropoints_fit_scaled_iord)) if k not in ind_to_remove]
        #     zeropoint_poly_scaled_iord = [zeropoint_poly_scaled_iord[k] for k in range(len(zeropoint_poly_scaled_iord)) if k not in ind_to_remove]
        #     log10blaze_iord = [log10blaze_iord[k] for k in range(len(log10blaze_iord)) if k not in ind_to_remove]
        #     waves_iord = [waves_iord[k] for k in range(len(waves_iord)) if k not in ind_to_remove]
        #     gpms_iord = [gpms_iord[k] for k in range(len(gpms_iord)) if k not in ind_to_remove]
        #     sensfiles_iord = [sensfiles_iord[k] for k in range(len(sensfiles_iord)) if k not in ind_to_remove]

        if len(waves_iord) > 0:
            comb_orders_list.append(iord)
            zeropoints_fit_list.append(zeropoints_fit_iord)
            zeropoints_data_list.append(zeropoints_data_iord)
            zeropoint_poly_list.append(zeropoint_poly_iord)
            scale_list.append(scales)
            logblaze_list.append(log10blaze_iord)
            waves_list.append(waves_iord)
            gpm_list.append(gpms_iord)
            sensfiles_iord_list.append(sensfiles_iord)

    # get the wavelength grid from the telluric grid
    hdul = io.load_telluric_grid(sensfuncs[0].telluric.telgrid)
    wave_grid_tell = hdul[1].data
    # keep only the wavelengths that are within the range of the zeropoints
    wmin_grid = np.array([[w.min(), w.max()] for sublist in waves_list for w in sublist]).min()
    wmax_grid = np.array([[w.min(), w.max()] for sublist in waves_list for w in sublist]).max()
    wave_grid = wave_grid_tell[(wave_grid_tell > wmin_grid) & (wave_grid_tell < wmax_grid)]

    for iord, order in enumerate(comb_orders_list):
        # interpolate the zeropoints to the wave_grid
        zeropoint_poly_iord = np.concatenate(zeropoint_poly_list[iord])
        wave_iord = np.concatenate(waves_list[iord])
        gpm_iord = np.concatenate(gpm_list[iord])
        pypeitFit = fitting.robust_fit(wave_iord, zeropoint_poly_iord, poly_order,
                                       function='legendre', minx=wave_iord.min(), maxx=wave_iord.max(),
                                       in_gpm=gpm_iord, lower=3, upper=3, use_mad=True)
        comb_zeropoint_fit_iord = flux_calib.eval_zeropoint(pypeitFit.fitc, 'legendre', wave_grid,
                                                            wave_iord.min(), wave_iord.max())
        comb_wave = wave_grid.copy()
        outside_gpm = (comb_wave < wave_iord.min()) | (comb_wave > wave_iord.max())
        comb_zeropoint_fit_iord[outside_gpm] = 0.
        comb_wave[outside_gpm] = 0.
        # append the results to the lists
        comb_zeropoint_fit_list.append(comb_zeropoint_fit_iord)
        comb_wave_list.append(comb_wave)
        comb_coeff_list.append(pypeitFit.fitc)

    comb_sensobj = fill_combSensObj(comb_zeropoint_fit_list, comb_wave_list, comb_orders_list, poly_order,
                                    comb_coeff_list, deepcopy(sensfuncs[0]))

    # QA plots
    if qa_plots:
        plot_by_order(comb_sensobj, zeropoints_data_list, zeropoints_fit_list, waves_list, scale_list,
                      sensfiles_iord_list, logblaze=logblaze_list)


    if qa_plots:
        # plot the zeropoints for all orders
        plot_all_orders(comb_sensobj, comb_type='zeropoint')
        # plot the throughput for all orders
        plot_all_orders(comb_sensobj,comb_type='throughput')

    return comb_sensobj


# def combine_sensfuncs(sensfuncs, cut_left=50, cut_right=-100, poly_order=4, qa_plots=False):
#     """Combine the sensitivity functions.
#
#     Args:
#         sensfuncs (list):
#             List of SensFunc objects.
#         cut_left (int):
#             Number of pixels to cut from the left.
#         cut_right (int):
#             Number of pixels to cut from the right.
#         poly_order (int):
#             Order of the polynomial to fit the zeropoints.
#         qa_plots (bool):
#             Generate QA plots.
#
#     Returns:
#         SensFunc: Combined sensitivity function object.
#     """
#
#     all_orders = np.array([])
#     for sensobj in sensfuncs:
#         all_orders = np.append(all_orders, sensobj.sens['ECH_ORDERS'].data)
#     order_vec = np.arange(all_orders.min(), all_orders.max() + 1, dtype=int)
#     # invert the order vector so that it goes from the highest order to the lowest (i.e., from blue to red)
#     order_vec = order_vec[::-1]
#
#     # initialize the combined zeropoints
#     comb_zeropoint_fit_list = []
#     comb_wave_list = []
#     comb_orders_list = []
#     comb_coeff_list = []
#     scale_list = []
#     # list of the zeropoints for each order
#     zeropoints_fit_list = []
#     zeropoints_data_list = []
#     zeropoint_poly_list = []
#     logblaze_list = []
#     waves_list = []
#     sensfiles_iord_list = []
#     gpm_list = []
#     # loop over the orders
#     for iord in order_vec:
#         zeropoint_medians = []
#         waves_iord, zeropoints_data_iord, zeropoints_fit_iord, zeropoint_poly_iord, gpms_iord, sensfiles_iord, log10blaze_iord = \
#             [], [], [], [], [], [], []
#         for s, sensobj in enumerate(sensfuncs):
#             for sens in sensobj.sens:
#                 if sens['ECH_ORDERS'] == iord:
#                     wave = sens['SENS_WAVE']
#                     wv_gpm = wave > 1.0
#                     wave = wave[wv_gpm][cut_left:cut_right]
#                     zeropoint_data = sens['SENS_ZEROPOINT'][wv_gpm][cut_left:cut_right]
#                     zeropoint_data_gpm = sens['SENS_ZEROPOINT_GPM'][wv_gpm][cut_left:cut_right]
#                     zeropoint_fit = sens['SENS_ZEROPOINT_FIT'][wv_gpm][cut_left:cut_right]
#                     zeropoint_fit_gpm = sens['SENS_ZEROPOINT_FIT_GPM'][wv_gpm][cut_left:cut_right]
#                     log10blaze = sens['SENS_LOG10_BLAZE_FUNCTION'][wv_gpm][cut_left:cut_right]
#                     # get zeropoint_poly
#                     zeropoint_poly = zeropoint_data + 5.0*np.log10(wave) - ZP_UNIT_CONST
#                     if not np.all(log10blaze==0):
#                         zeropoint_poly -= 2.5*log10blaze
#
#                     # compute median of the zeropoints
#                     zeropoint_medians.append(stats.sigma_clipped_stats(zeropoint_poly, mask_value=0., sigma=3)[1])
#                     waves_iord.append(wave)
#                     gpms_iord.append(zeropoint_data_gpm)
#                     zeropoints_fit_iord.append(zeropoint_fit)
#                     zeropoints_data_iord.append(zeropoint_data)
#                     zeropoint_poly_iord.append(zeropoint_poly)
#                     log10blaze_iord.append(log10blaze)
#                     sensfiles_iord.append(f"{Path(sensobj.spec1df).name.replace('spec1d', 'sens')}")
#
#         # scale the zeropoints to the same median value (max of all the medians)
#         zeropoint_medians = np.array(zeropoint_medians)
#         max_median = np.nanmax(zeropoint_medians)
#         scales = np.ones(len(zeropoint_medians)) * max_median / zeropoint_medians
#         zeropoints_fit_scaled_iord = [z * s for z, s in zip(zeropoints_fit_iord, scales)]
#         zeropoints_data_scaled_iord = [z * s for z, s in zip(zeropoints_data_iord, scales)]
#         zeropoint_poly_scaled_iord = [z * s for z, s in zip(zeropoint_poly_iord, scales)]
#
#         # discard the zeropoints that have values below a certain threshold. This allows to discard zeropoints
#         # that go down because of the blocking filter or other reasons.
#         # lower limit for the zeropoints
#         low_thresh = 7.
#         # if iord == 35:
#         #     low_thresh = 13.
#         # elif iord in [43,66,84]:
#         #     low_thresh = 18.
#         # elif iord == 88:
#         #     low_thresh = 17.
#         # elif iord == 76:
#         #     low_thresh = 17.5
#         ind_to_remove = []
#         for i in range(len(zeropoints_fit_scaled_iord)):
#             if np.any(zeropoints_fit_scaled_iord[i] < low_thresh) or np.any(zeropoints_fit_scaled_iord[i] > 20.):
#                 ind_to_remove.append(i)
#         if len(ind_to_remove) > 0:
#             zeropoints_data_scaled_iord = [zeropoints_data_scaled_iord[k] for k in range(len(zeropoints_data_scaled_iord)) if k not in ind_to_remove]
#             zeropoints_fit_scaled_iord = [zeropoints_fit_scaled_iord[k] for k in range(len(zeropoints_fit_scaled_iord)) if k not in ind_to_remove]
#             zeropoint_poly_scaled_iord = [zeropoint_poly_scaled_iord[k] for k in range(len(zeropoint_poly_scaled_iord)) if k not in ind_to_remove]
#             log10blaze_iord = [log10blaze_iord[k] for k in range(len(log10blaze_iord)) if k not in ind_to_remove]
#             waves_iord = [waves_iord[k] for k in range(len(waves_iord)) if k not in ind_to_remove]
#             gpms_iord = [gpms_iord[k] for k in range(len(gpms_iord)) if k not in ind_to_remove]
#             sensfiles_iord = [sensfiles_iord[k] for k in range(len(sensfiles_iord)) if k not in ind_to_remove]
#
#         if len(waves_iord) > 0:
#             comb_orders_list.append(iord)
#             zeropoints_fit_list.append(zeropoints_fit_iord)
#             zeropoints_data_list.append(zeropoints_data_iord)
#             zeropoint_poly_list.append(zeropoint_poly_scaled_iord)
#             scale_list.append(scales)
#             logblaze_list.append(log10blaze_iord)
#             waves_list.append(waves_iord)
#             gpm_list.append(gpms_iord)
#             sensfiles_iord_list.append(sensfiles_iord)
#
#     # get the wavelength grid from the telluric grid
#     hdul = io.load_telluric_grid(sensfuncs[0].telluric.telgrid)
#     wave_grid_tell = hdul[1].data
#     # keep only the wavelengths that are within the range of the zeropoints
#     wmin_grid = np.array([[w.min(), w.max()] for sublist in waves_list for w in sublist]).min()
#     wmax_grid = np.array([[w.min(), w.max()] for sublist in waves_list for w in sublist]).max()
#     wave_grid = wave_grid_tell[(wave_grid_tell > wmin_grid) & (wave_grid_tell < wmax_grid)]
#
#     for iord, order in enumerate(comb_orders_list):
#         # interpolate the zeropoints to the wave_grid
#         zeropoint_poly_iord = np.concatenate(zeropoint_poly_list[iord])
#         wave_iord = np.concatenate(waves_list[iord])
#         gpm_iord = np.concatenate(gpm_list[iord])
#         pypeitFit = fitting.robust_fit(wave_iord, zeropoint_poly_iord, poly_order,
#                                        function='legendre', minx=wave_iord.min(), maxx=wave_iord.max(),
#                                        in_gpm=gpm_iord, lower=3, upper=3, use_mad=True)
#         comb_zeropoint_fit_iord = flux_calib.eval_zeropoint(pypeitFit.fitc, 'legendre', wave_grid,
#                                                             wave_iord.min(), wave_iord.max())
#         comb_wave = wave_grid.copy()
#         outside_gpm = (comb_wave < wave_iord.min()) | (comb_wave > wave_iord.max())
#         comb_zeropoint_fit_iord[outside_gpm] = 0.
#         comb_wave[outside_gpm] = 0.
#         # append the results to the lists
#         comb_zeropoint_fit_list.append(comb_zeropoint_fit_iord)
#         comb_wave_list.append(comb_wave)
#         comb_coeff_list.append(pypeitFit.fitc)
#
#     comb_sensobj = fill_combSensObj(comb_zeropoint_fit_list, comb_wave_list, comb_orders_list, poly_order,
#                                     comb_coeff_list, deepcopy(sensfuncs[0]))
#
#     # QA plots
#     if qa_plots:
#         plot_by_order(comb_sensobj, zeropoints_data_list, zeropoints_fit_list, waves_list, scale_list,
#                       sensfiles_iord_list, logblaze=logblaze_list)
#
#
#     if qa_plots:
#         # plot the zeropoints for all orders
#         plot_all_orders(comb_sensobj, comb_type='zeropoint')
#         # plot the throughput for all orders
#         plot_all_orders(comb_sensobj,comb_type='throughput')
#
#     return comb_sensobj


def plot_by_order(comb_sensobj, zeropoints_data, zeropoints_fit, waves_list, scale_list, sensfile_list, logblaze=None):
    """Plot the zeropoints by order.

    Args:
        comb_sensobj (SensFunc):
            Combined sensitivity function object.
        zeropoints_data (list):
            List of zeropoints data that were combined (several per order)
        zeropoints_fit (list):
            List of zeropoints fits that were combined (several per order)
        waves_list (list):
            List of wavelength arrays for the zeropoints (several per order)
        scale_list (list):
            List of scales used to rescale the zeropoints (several per order)
        sensfile_list (list):
            List of sensitivity function file names (several per order)
        logblaze (list, optional):
            List of log10 blaze functions (several per order). If None the blaze is not applied.
    """

    unique_sensfiles = np.unique([sensfile for sublist in sensfile_list for sensfile in sublist])
    unique_colors = get_unique_colors(unique_sensfiles.size, mcolors.CSS4_COLORS)
    color_map = {sensobj: unique_colors[i] for i, sensobj in enumerate(unique_sensfiles)}

    outfile = 'collection_RESCALEDzeropoints_by_order_with_comb.pdf'

    # get the values to plot
    orders = comb_sensobj.sens['ECH_ORDERS'].data
    waves = comb_sensobj.wave.T
    comb_zeropoint_fit = comb_sensobj.zeropoint.T # this is without the blaze correction
    comb_coeff = comb_sensobj.sens['SENS_COEFF'].data

    with PdfPages(outfile) as pdf:
        for (order, wave, comb_zeropoint, coeff, zeropoint_data, zeropoint_fits, waves_iord, scales, lblaze, sensfiles_iord) in \
                zip(orders, waves, comb_zeropoint_fit, comb_coeff, zeropoints_data, zeropoints_fit, waves_list, scale_list, logblaze, sensfile_list):

            fig = plt.figure(figsize=(23, 6.))
            axis = plt.subplot()
            plt.minorticks_on()
            plt.tick_params(axis='both', direction='in', top=True, right=True, which='both')
            plt.title(f'Order {order}')
            for i in range(len(sensfiles_iord)):
                label = f"{sensfiles_iord[i]}"
                color = color_map[sensfiles_iord[i]]

                axis.plot(waves_iord[i], zeropoint_data[i], drawstyle='steps-mid',
                          color=color, linewidth=0.6, alpha=0.4, zorder=-2)
                axis.plot(waves_iord[i], zeropoint_fits[i], color=color, linewidth=1., alpha=1.,
                          label=label, zorder=0)
                # recover the original zeropoints using the combined coefficients
                this_comb = flux_calib.eval_zeropoint(coeff, 'legendre', waves_iord[i],
                                                            waves_iord[i].min(), waves_iord[i].max(),
                                                            log10_blaze_func_per_ang=lblaze[i])
                axis.plot(waves_iord[i], this_comb, color='k', marker='x', markersize=5, linewidth=None,
                          alpha=0.4, zorder=2)

            axis.plot(wave[wave>1], comb_zeropoint[wave>1], color='k', linewidth=2.5, ls='--', alpha=1., zorder=1,
                      label='combined zeropoint - no blaze')

            #axis.set_ylim(13.1, 20.9)
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
    #axis.set_ylim(13.1, 20.9) if comb_type == 'zeropoint' else axis.set_ylim(0.0, 0.15)
    axis.set_xlim(3600., 10400.)
    axis.legend(loc='lower left', fontsize=4)
    axis.set_xlabel('Wavelength (Angstroms)')
    axis.set_ylabel('Zeropoint (AB mag)') if comb_type == 'zeropoint' else axis.set_ylabel('Throughput')
    fig.tight_layout()
    plt.savefig(outfile, dpi=60)
    # plt.show()
    plt.close(fig)

def only_these_datasets():
    return np.array(['HI.20191209.57374.fits', 'HI.20191209.16095.fits', 'HI.20151214.15526.fits', 'HI.20080704.52715.fits',
         'HI.20151215.16490.fits', 'HI.20151215.16075.fits', 'HI.20151214.16715.fits', 'HI.20151214.16343.fits'])


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
        if args.sens_fnames is not None:
            sens_fnames = np.loadtxt(args.sens_fnames, dtype=str)
        else:
            sens_fnames = [f.name for f in sensfuncs_path.glob('sens*.fits')]
        if len(sens_fnames) == 0:
            msgs.error(f"No sensitivity function files found.")
        sensfuncs = load_sensfunc(sens_fnames, sensfuncs_path)

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
        sensfuncs, sens_fnames = create_sens_files(spec1d_fnames, spec1ds_path, sensfuncs_path,
                                                   boxcar=args.boxcar, use_flat=args.use_flat,
                                                   skip_existing=args.skip_existing, only=only_these_datasets())

    # for debugging
    # tab_all, aa = select_datasets(sensfuncs)

    # combine the sensitivity functions
    sensfunc = combine_sensfuncs(sensfuncs, qa_plots=args.QA)

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


