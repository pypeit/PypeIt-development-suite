# after running the script hires_zeropoint_collection.py, run this script to plot the
# collection of zeropoints in a single figure.



import os
from pathlib import Path
from pypeit import sensfunc
from pypeit.core.wavecal import wvutils
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
import scipy
from astropy import stats
from IPython import embed

# for plots
plt.rc('font', family='serif')
plt.rcParams['mathtext.default'] = u'regular'
plt.rc('axes', labelsize=15)
plt.rc('legend', fontsize=7, borderpad=1., handlelength=1., handleheight=1., labelspacing=0.5, handletextpad=0.5,
       borderaxespad=0.5, columnspacing=0.)
plt.rc('xtick', labelsize=12, direction='in', top=True)
plt.rc('xtick.major', size=8, width=1.)
plt.rc('xtick.minor', size=4, width=1.)
plt.rc('ytick', labelsize=12, direction='in', right=True)
plt.rc('ytick.major', size=8, width=0.8)
plt.rc('ytick.minor', size=4, width=0.8)
plt.rc('axes', linewidth=1.5)

plot_by_order = True


# REDUX/keck_hires path
# hires_redux = Path(os.getenv('PYPEIT_DEV') + '/REDUX_OUT/keck_hires/')
hires_redux = Path('/Users/dpelliccia/Desktop/adap2020/')
# hires_redux = Path('/Volumes/GoogleDrive/Shared drives/PypeIt ADAP 2020/sensfuncs/')
# grab all the sens*.fits files
# sens_files = list(hires_redux.glob('*/sens*.fits'))
sens_files = list(hires_redux.glob('*/*/sens*.fits'))

# list of sens objects
sensobjs_list = []
wmins = []
wmaxs = []
datasets = []
all_orders = np.array([])

for sfile in sens_files:
    _sens = sensfunc.SensFunc.from_file(sfile, chk_version=False)
    sensobjs_list.append(_sens)
    wmins.append(_sens.sens['WAVE_MIN'].data.min())
    wmaxs.append(_sens.sens['WAVE_MAX'].data.max())
    datasets.append(sfile.parent.name)
    all_orders = np.append(all_orders, _sens.sens['ECH_ORDERS'].data)

order_vec = np.arange(all_orders.min(), all_orders.max()+1, dtype=int) if plot_by_order else [1]

wave_grids = []
if order_vec.size > 1:
    # create wave grid for each order to be used for the extrapolation of the zeropoints
    for iord in order_vec[::-1]:
        _wave_mins = []
        _wave_maxs = []
        dwaves = []
        for sensobj in sensobjs_list:
            sens_tab = sensobj.sens
            for sens in sens_tab:
                if sens['ECH_ORDERS'] == iord:
                    wave = sens['SENS_WAVE']
                    wv_gpm = (wave > 1.0) & sens['SENS_ZEROPOINT_GPM']
                    _wave_mins.append(wave[wv_gpm].min())
                    _wave_maxs.append(wave[wv_gpm].max())
                    dwave_data, _, _, _ = wvutils.get_sampling(wave[wv_gpm])
                    dwaves.append(dwave_data)
        wave_min = np.min(_wave_mins)
        wave_max = np.max(_wave_maxs)
        dwave = np.min(dwaves)
        nspec = int(np.ceil((wave_max - wave_min) / dwave))
        _wave_grid = np.arange(nspec)*((wave_max - wave_min) / (nspec - 1)) + (np.ones(nspec)*wave_min)
        wave_grids.append(_wave_grid)

# colors for the plots
colors_values = np.array(list(mcolors.CSS4_COLORS.values()))
# exclude light colors. Colors that have v (of hsv) values less than 0.85 are kept
hsv_v = np.array([mcolors.rgb_to_hsv(mcolors.to_rgba(c)[:3])[-1] for c in colors_values])
pcolors = colors_values[hsv_v < 0.9]


# wave ranges for the plots
plot_wmin = np.min(wmins)
plot_wmax = np.max(wmaxs)
if plot_by_order:
    wave_ranges = [(plot_wmin, plot_wmax)]
else:
    wv_break = (plot_wmax - plot_wmin)/6.
    wave_ranges = [(plot_wmin, plot_wmax - 5*wv_break), (plot_wmin + wv_break, plot_wmax - 4*wv_break),
                   (plot_wmin + 2*wv_break, plot_wmax - 3*wv_break), (plot_wmin + 3*wv_break, plot_wmax - 2*wv_break),
                   (plot_wmin + 4*wv_break, plot_wmax - wv_break), (plot_wmin + 5*wv_break, plot_wmax)]

# plot the zeropoints
plotname = 'collection_zeropoints.pdf' if plot_by_order is False else 'collection_zeropoints_by_order_with_comb.pdf'

with PdfPages(plotname) as pdf:
    for wr in wave_ranges:
        _wmin, _wmax = wr
        for o,iord in enumerate(order_vec[::-1]):
            fig = plt.figure(figsize=(23, 6.))
            axis = plt.subplot()
            plt.minorticks_on()
            zeropoint_fit_iord = np.zeros((len(sensobjs_list),len(wave_grids[o])))
            zeropoint_medians = np.ones((len(sensobjs_list)))
            for i,sensobj in enumerate(sensobjs_list):
                nplot = 0
                stdname = sensobj.std_name
                sens_tab = sensobj.sens
                for sens in sens_tab:
                    if (sens['ECH_ORDERS'] == iord) if plot_by_order else (sens['WAVE_MIN'] >= _wmin and sens['WAVE_MAX'] <= _wmax):
                        wave = sens['SENS_WAVE']
                        zeropoint_data = sens['SENS_ZEROPOINT']
                        zeropoint_data_gpm = sens['SENS_ZEROPOINT_GPM']
                        zeropoint_fit = sens['SENS_ZEROPOINT_FIT']
                        zeropoint_fit_gpm = sens['SENS_ZEROPOINT_FIT_GPM']
                        wv_gpm = (wave > 1.0) & zeropoint_data_gpm

                        zeropoint_fit_extr = scipy.interpolate.interp1d(wave[wv_gpm], zeropoint_fit[wv_gpm],
                                                                        kind='linear', bounds_error=False,
                                                                        fill_value='extrapolate')(wave_grids[o])
                        # actually mask the extrapolated values. We are extrapolating
                        # so that the zeropoints are on the same wavegrid and we can combine them
                        mask = (wave_grids[o] <= wave[wv_gpm].min()) | (wave_grids[o] >= wave[wv_gpm].max())
                        zeropoint_fit_extr[mask] = 0
                        zeropoint_fit_iord[i] = zeropoint_fit_extr
                        zeropoint_medians[i] = np.nanmax(zeropoint_fit[wv_gpm])
                        # rejmask = zeropoint_data_gpm[wv_gpm] & np.logical_not(zeropoint_fit_gpm[wv_gpm])
                        axis.plot(wave[wv_gpm], zeropoint_data[wv_gpm], drawstyle='steps-mid',
                                  color=pcolors[i], linewidth=0.6, alpha=0.4, zorder=-2)
                        if nplot == 0:
                            label = f'{stdname} - {datasets[i]}' if plot_by_order is False \
                                else f'{stdname} - {datasets[i]} - order {iord}'
                        else:
                            label = None
                        axis.plot(wave[wv_gpm], zeropoint_fit[wv_gpm], color=pcolors[i], linewidth=2.0, alpha=1.,
                                  label=label,zorder=0)

                        # axis.plot(wave_grids[o], zeropoint_fit_extr, color='k', ls='--', linewidth=2.0, alpha=1., zorder=0)

                        # axis.plot(wave[wv_gpm][rejmask], zeropoint_data[wv_gpm][rejmask], 's', zorder=10, mfc='None',
                        #           mec='blue', mew=0.7, label='rejected pixels from fit')
                        # axis.plot(wave[wv_gpm][np.logical_not(zeropoint_data_gpm[wv_gpm])],
                        #           zeropoint_data[wv_gpm][np.logical_not(zeropoint_data_gpm[wv_gpm])], 'v',
                        #           zorder=11, mfc='None', mec='orange', mew=0.7, label='originally masked')
                        nplot += 1
            if len(list(axis.get_lines())) == 0:
                plt.close()
                continue
            med_median = np.median(zeropoint_medians[zeropoint_medians!=1])
            weights = np.ones(len(zeropoint_medians)) * med_median/zeropoint_medians
            weights = np.tile(weights, (len(wave_grids[o]), 1)).T

            comb_zeropoint = stats.sigma_clipped_stats(weights*zeropoint_fit_iord, mask_value=0., axis=0, sigma=3)[1]
            axis.plot(wave_grids[o], comb_zeropoint, color='r', linewidth=2.0, alpha=1., zorder=1, label='combined zeropoint')
            axis.set_ylim(9,21)
            # axis.set_xlim(_wmin, _wmax)
            axis.legend()
            axis.set_xlabel('Wavelength')
            axis.set_ylabel('Zeropoint (AB mag)')
            fig.tight_layout()
            pdf.savefig(dpi=60)
            # plt.show()

