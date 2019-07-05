
import numpy as np
import os
from pypeit.traceslits import TraceSlits
from pypeit.wavetilts import WaveTilts
from pypeit.spectrographs import util
from pypeit.core import pixels
from astropy.io import fits
from pypeit import ginga
from astropy import stats
import scipy
from pypeit import ginga
from pypeit.core import arc
from pypeit import utils
from matplotlib import pyplot as plt
from pypeit.core import trace_slits
from pypeit import msgs


dev_path = os.getenv('PYPEIT_DEV')
lris_path =  os.path.join(dev_path, 'REDUX_OUT/Keck_LRIS_blue/multi_600_4000_d560')
trc_file = os.path.join(lris_path, 'Masters', 'MasterTrace_A_1_01.fits')
tilts_file = os.path.join(lris_path, 'Masters', 'MasterTilts_A_1_01.fits')
spec2dfile = os.path.join(lris_path, 'Science', 'spec2d_b170320_2083-c17_60L._LRISb_2017Mar20T055336.211.fits')

spectrograph = util.load_spectrograph('keck_lris_blue')
par = spectrograph.default_pypeit_par()
tslits_dict = TraceSlits.load_from_file(trc_file)[0]


tilts_dict = WaveTilts.load_from_file(tilts_file, return_header=False)

hdu = fits.open(spec2dfile)
sciimg = hdu[1].data
nspec, nspat = sciimg.shape

slitmask = pixels.tslits2mask(tslits_dict)
onslits = (slitmask > -1)
corr_slits = (onslits.astype(float)).flatten()
#corr_roll = np.roll(corr_slits, 10, axis=1).astype(float)

# Compute
(mean_sci, med_sci, stddev_sci) = stats.sigma_clipped_stats(sciimg[onslits])
thresh =  med_sci + 5.0*stddev_sci
corr_sci = np.fmin(sciimg.flatten(), thresh)


maxlag = 20
lags, xcorr = utils.cross_correlate(corr_sci, corr_slits, maxlag)
xcorr_denom = np.sqrt(np.sum(corr_sci*corr_sci)*np.sum(corr_slits*corr_slits))
xcorr_norm = xcorr / xcorr_denom
tampl_true, tampl, pix_max, twid, centerr, ww, arc_cont, nsig = arc.detect_lines(xcorr_norm, sigdetect=3.0,
                                                                                 fit_frac_fwhm=1.5, fwhm=5.0,
                                                                                 cont_frac_fwhm=1.0, cont_samp=30,
                                                                                 nfind=1, debug=True)
xcorr_max = np.interp(pix_max, np.arange(lags.shape[0]), xcorr_norm)
lag_max = np.interp(pix_max, np.arange(lags.shape[0]), lags)
msgs.info('Flexure compensation ')
debug=True
if debug:
    # Interpolate for bad lines since the fitting code often returns nan
    plt.figure(figsize=(14, 6))
    plt.plot(lags, xcorr_norm, color='black', drawstyle='steps-mid', lw=3, label='x-corr', linewidth=1.0)
    plt.plot(lag_max[0], xcorr_max[0], 'g+', markersize=6.0, label='peak')
    plt.title('Best shift = {:5.3f}'.format(lag_max[0]) + ',  corr_max = {:5.3f}'.format(xcorr_max[0]))
    plt.legend()
    plt.show()


# Now translate the slits in the tslits_dict
tslits_shift = trace_slits.shift_slits(tslits_dict, lag_max)
# Now translate the tilts

slitmask_shift = pixels.tslits2mask(tslits_shift)
viewer, ch = ginga.show_image (sciimg)
ginga.show_slits(viewer, ch, tslits_shift['slit_left'], tslits_shift['slit_righ'])
ginga.show_slits(viewer, ch, tslits_dict['slit_left'], tslits_dict['slit_righ'])



sys.exit(-1)




#lags = np.arange(-nspat + 1, nspat)
#corr = scipy.signal.correlate(corr_sci, corr_slits, mode='full')





if debug:
    # Interpolate for bad lines since the fitting code often returns nan
    plt.figure(figsize=(14, 6))
    plt.plot(lags, corr_norm, color='black', drawstyle='steps-mid', lw=3, label='x-corr', linewidth=1.0)
    plt.plot(lag_max[0], corr_max[0], 'g+', markersize=6.0, label='peak')
    plt.title('Best shift = {:5.3f}'.format(lag_max[0]) + ',  corr_max = {:5.3f}'.format(corr_max[0]))
    plt.legend()
    plt.show()




