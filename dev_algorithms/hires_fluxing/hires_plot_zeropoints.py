# after running the script hires_zeropoint_collection.py, run this script to plot the
# collection of zeropoints in a single figure.



import os
from pathlib import Path
from pypeit import sensfunc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
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


# REDUX/keck_hires path
# hires_redux = Path(os.getenv('PYPEIT_DEV') + '/REDUX_OUT/keck_hires/')
hires_redux = Path('/Users/dpelliccia/Desktop/adap2020/')
# grab all the sens*.fits files
# sens_files = list(hires_redux.glob('*/sens*.fits'))
sens_files = list(hires_redux.glob('*/*/sens*.fits'))

# list of sens objects
sensobjs_list = []
wmins = []
wmaxs = []
datasets = []

for sfile in sens_files:
    _sens = sensfunc.SensFunc.from_file(sfile, chk_version=False)
    sensobjs_list.append(_sens)
    wmins.append(_sens.sens['WAVE_MIN'].data.min())
    wmaxs.append(_sens.sens['WAVE_MAX'].data.max())
    datasets.append(sfile.parent.name)

# colors for the plots
colors_values = np.array(list(mcolors.CSS4_COLORS.values()))
# exclude light colors. Colors that have v (of hsv) values less than 0.85 are kept
hsv_v = np.array([mcolors.rgb_to_hsv(mcolors.to_rgba(c)[:3])[-1] for c in colors_values])
pcolors = colors_values[hsv_v < 0.9]


# wave ranges for the plots
plot_wmin = np.min(wmins) #- 50.
plot_wmax = np.max(wmaxs) #+ 50.
wv_break = (plot_wmax - plot_wmin)/6.
wave_ranges = [(plot_wmin, plot_wmin + wv_break), (plot_wmin + wv_break, plot_wmin + 2*wv_break),
               (plot_wmin + 2*wv_break, plot_wmin + 3*wv_break), (plot_wmin + 3*wv_break, plot_wmin + 4*wv_break),
               (plot_wmin + 4*wv_break, plot_wmin + 5*wv_break), (plot_wmin + 5*wv_break, plot_wmax)]

# plot the zeropoints
plotname = 'collection_zeropoints.pdf'
with PdfPages(plotname) as pdf:
    for wr in wave_ranges:
        _wmin, _wmax = wr
        fig = plt.figure(figsize=(23, 6.))
        axis = plt.subplot()
        plt.minorticks_on()
        for i,sensobj in enumerate(sensobjs_list):
            nplot = 0
            stdname = sensobj.std_name
            sens_tab = sensobj.sens
            for sens in sens_tab:
                if sens['WAVE_MIN'] > _wmin and sens['WAVE_MAX'] <= _wmax:
                    wave = sens['SENS_WAVE']
                    zeropoint_data = sens['SENS_ZEROPOINT']
                    zeropoint_data_gpm = sens['SENS_ZEROPOINT_GPM']
                    zeropoint_fit = sens['SENS_ZEROPOINT_FIT']
                    zeropoint_fit_gpm = sens['SENS_ZEROPOINT_FIT_GPM']

                    wv_gpm = wave > 1.0
                    # rejmask = zeropoint_data_gpm[wv_gpm] & np.logical_not(zeropoint_fit_gpm[wv_gpm])
                    axis.plot(wave[wv_gpm], zeropoint_data[wv_gpm], drawstyle='steps-mid',
                              color=pcolors[i], linewidth=0.6, alpha=0.4, zorder=5)
                    if nplot == 0:
                        axis.plot(wave[wv_gpm], zeropoint_fit[wv_gpm], label=f'{stdname} - {datasets[i]}',
                                  color=pcolors[i], linewidth=1.3, alpha=1., zorder=0)
                    else:
                        axis.plot(wave[wv_gpm], zeropoint_fit[wv_gpm], color=pcolors[i], linewidth=2.0, alpha=1., zorder=0)
                    # axis.plot(wave[wv_gpm][rejmask], zeropoint_data[wv_gpm][rejmask], 's', zorder=10, mfc='None',
                    #           mec='blue', mew=0.7, label='rejected pixels from fit')
                    # axis.plot(wave[wv_gpm][np.logical_not(zeropoint_data_gpm[wv_gpm])],
                    #           zeropoint_data[wv_gpm][np.logical_not(zeropoint_data_gpm[wv_gpm])], 'v',
                    #           zorder=11, mfc='None', mec='orange', mew=0.7, label='originally masked')
                    nplot += 1
        axis.set_ylim(9,21)
        axis.legend(loc='lower left')
        axis.set_xlabel('Wavelength')
        axis.set_ylabel('Zeropoint (AB mag)')
        fig.tight_layout()
        pdf.savefig(dpi=60)
        # plt.show()
