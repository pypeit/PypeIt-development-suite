# after running the script hires_zeropoint_collection.py, run this script to plot the
# collection of zeropoints in a single figure.



import os
from pathlib import Path
from pypeit import sensfunc
from pypeit.core.wavecal import wvutils
from pypeit.core import fitting
import numpy as np
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
plt.rc('legend', fontsize=4, borderpad=1., handlelength=1., handleheight=1., labelspacing=0.5, handletextpad=0.5,
       borderaxespad=0.5, columnspacing=0.)
plt.rc('xtick', labelsize=12, direction='in', top=True)
plt.rc('xtick.major', size=8, width=1.)
plt.rc('xtick.minor', size=4, width=1.)
plt.rc('ytick', labelsize=12, direction='in', right=True)
plt.rc('ytick.major', size=8, width=0.8)
plt.rc('ytick.minor', size=4, width=0.8)
plt.rc('axes', linewidth=1.5)

# colors for the plots
colors_values = np.array(list(mcolors.XKCD_COLORS.values()))
# exclude light and dark colors. Colors that have v (of hsv) values between 0.15 and 0.85 are kept
hsv_v = np.array([mcolors.rgb_to_hsv(mcolors.to_rgba(c)[:3])[-1] for c in colors_values])
pcolors = colors_values[(hsv_v >= 0.2) & (hsv_v <= 0.85)]

plot_by_order = True


# REDUX/keck_hires path
hires_redux = Path('/Users/dpelliccia/Desktop/adap2020/')
# hires_redux = Path('/Volumes/GoogleDrive/Shared drives/PypeIt ADAP 2020/sensfuncs/')
# grab all the sens*.fits files
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

# trim the sensfuncs 50 pixels from left and 100 pixels from right. This is to avoid the edges drops
cut_left, cut_right = 50, -100

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
plotname = 'collection_RESCALEDzeropoints.pdf' if plot_by_order is False else 'collection_RESCALEDzeropoints_by_order_with_comb.pdf'

with PdfPages(plotname) as pdf:
    for wr in wave_ranges:
        _wmin, _wmax = wr
        for o,iord in enumerate(order_vec[::-1]):
            zeropoint_medians = []
            datasets_iord, stdnames_iord = [], []
            waves_iord, zeropoints_data_iord, zeropoints_fit_iord, gpms_iord = [], [], [], []
            for i,sensobj in enumerate(sensobjs_list):
                stdname = sensobj.std_name
                sens_tab = sensobj.sens
                for sens in sens_tab:
                    if (sens['ECH_ORDERS'] == iord) if plot_by_order else (sens['WAVE_MIN'] >= _wmin and sens['WAVE_MAX'] <= _wmax):
                        wave = sens['SENS_WAVE']
                        wv_gpm = wave > 1.0
                        wave = wave[wv_gpm][cut_left:cut_right]
                        zeropoint_data = sens['SENS_ZEROPOINT'][wv_gpm][cut_left:cut_right]
                        zeropoint_data_gpm = sens['SENS_ZEROPOINT_GPM'][wv_gpm][cut_left:cut_right]
                        zeropoint_fit = sens['SENS_ZEROPOINT_FIT'][wv_gpm][cut_left:cut_right]
                        zeropoint_fit_gpm = sens['SENS_ZEROPOINT_FIT_GPM'][wv_gpm][cut_left:cut_right]
                        # compute median of the zeropoints
                        zeropoint_medians.append(stats.sigma_clipped_stats(zeropoint_fit, mask_value=0., sigma=3)[1])
                        waves_iord.append(wave)
                        gpms_iord.append(zeropoint_data_gpm & zeropoint_fit_gpm)
                        zeropoints_fit_iord.append(zeropoint_fit)
                        zeropoints_data_iord.append(zeropoint_data)
                        datasets_iord.append(datasets[i])
                        stdnames_iord.append(stdname)

            # scale the zeropoints to the same median value (max of all the medians)
            zeropoint_medians = np.array(zeropoint_medians)
            max_median = np.nanmax(zeropoint_medians)
            scales = np.ones(len(zeropoint_medians)) * max_median/zeropoint_medians
            zeropoints_fit_scaled = [z*s for z,s in zip(zeropoints_fit_iord, scales)]
            zeropoints_data_scaled = [z*s for z,s in zip(zeropoints_data_iord, scales)]
            # discard the zeropoints that have values below a certain threshold. This allows to discard zeropoints
            # that go down because of the blocking filter or other reasons. Then plot the remaining zeropoints
            fig = plt.figure(figsize=(23, 6.))
            axis = plt.subplot()
            plt.minorticks_on()
            ok_zeropoints = []
            ok_waves = []
            ok_gpms = []
            # lower limit for the zeropoints
            low_thresh = 16.
            if iord == 35:
                low_thresh = 13.
            elif iord in [43,66,84]:
                low_thresh = 18.
            elif iord == 88:
                low_thresh = 17.
            for i in range(len(zeropoints_fit_scaled)):
                if np.any(zeropoints_fit_scaled[i] < low_thresh) or np.any(zeropoints_fit_scaled[i] > 20.):
                    continue
                ok_zeropoints.append(zeropoints_fit_scaled[i])
                ok_waves.append(waves_iord[i])
                ok_gpms.append(gpms_iord[i])
                # plot the zeropoints
                label = f'{stdnames_iord[i]} - {datasets_iord[i]}' if plot_by_order is False \
                    else f'{stdnames_iord[i]} - {datasets_iord[i]} - order {iord}'

                axis.plot(waves_iord[i], zeropoints_data_scaled[i], drawstyle='steps-mid',
                          color=pcolors[i], linewidth=0.6, alpha=0.4, zorder=-2)
                axis.plot(waves_iord[i], zeropoints_fit_scaled[i], color=pcolors[i], linewidth=1., alpha=1.,
                          label=label, zorder=0)

            if len(ok_waves) > 0:
                if plot_by_order:
                    # create wave grid
                    wave_grid, wave_grid_mid, _ = wvutils.get_wave_grid(waves=ok_waves)
                    # fit a polynomial to all the zeropoints
                    poly_order = 4  # Define the order of the polynomial
                    all_zeropoints = np.concatenate(ok_zeropoints)
                    all_waves = np.concatenate(ok_waves)
                    all_gpms = np.concatenate(ok_gpms)
                    wave_min, wave_max = all_waves.min(), all_waves.max()
                    pypeitFit = fitting.robust_fit(all_waves, all_zeropoints, poly_order, minx=wave_min, maxx=wave_max,
                                                   in_gpm=all_gpms, lower=3, upper=3, use_mad=True)
                    # evaluate the polynomial
                    combined_zeropoints = pypeitFit.eval(wave_grid)
                    axis.plot(wave_grid, combined_zeropoints, color='k', linewidth=2.5, ls='--', alpha=1., zorder=1,
                              label='combined zeropoint')

                axis.set_ylim(13.1, 20.9)
                # axis.set_xlim(_wmin, _wmax)
                axis.legend()
                axis.set_xlabel('Wavelength')
                axis.set_ylabel('Zeropoint (AB mag)')
                fig.tight_layout()
                pdf.savefig(dpi=60)
                # plt.show()
            plt.close(fig)




