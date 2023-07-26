import os
from IPython import embed

import numpy as np

from astropy.io import fits
from astropy.time import Time
from astropy.table import Table

from pypeit import msgs
from pypeit.core import meta
from pypeit import utils
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core import coadd, flux_calib
from matplotlib import pyplot as plt
from pypeit.utils import fast_running_median


spectrograph = load_spectrograph('vlt_xshooter_vis')
par = spectrograph.default_pypeit_par()
# redux_path = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/keck_hires/J0100+2802_RED_C1_ECH_-0.91_XD_1.46_1x2/'
redux_path = '/Users/silviaonorato/Projects/highz_qso_redux/highz_qso_redux/PypeIt_data/REDUX_OUT/vlt_xshooter/test_dev_suite/'
star = 'Feige110'
other_star = 'LTT3218'
#star = 'BD+28'
#other_star = 'G191'
if star == 'Feige110':
    coadd1dfile = os.path.join(redux_path, 'Feige110_coadd_fluxed_with_LTT3218.fits')
elif star == 'LTT3218':
    coadd1dfile = os.path.join(redux_path, 'LTT3218_coadd_fluxed_with_Feige110.fits')

# Read in the coadd1d
std_star_table = Table.read(coadd1dfile, format='fits')
wave = std_star_table['wave_grid_mid'].data
flux = std_star_table['flux'].data
ivar = std_star_table['ivar'].data
gpm  = std_star_table['mask'].data
sigma = np.sqrt(utils.inverse(ivar))

# Read in the header
header = fits.getheader(coadd1dfile)
star_ra, star_dec = header['RA'], header['DEC']

# Convert to decimal deg, as needed
star_ra, star_dec = meta.convert_radec(star_ra, star_dec)

# Read in standard star dictionary
std_dict = flux_calib.get_standard_spectrum(star_type=None, star_mag=None, ra=star_ra, dec=star_dec)

nsmooth = 50
flux_sm = fast_running_median(flux, nsmooth)
model_sm = fast_running_median(std_dict['flux'].value, nsmooth)
ymin_max = np.min(flux_sm), np.max(flux_sm)
wave_gpm = (wave > 1.0)

# Renormalize the data to match the model at 7500A #12000A
renorm = np.interp(7500., std_dict['wave'].value, std_dict['flux'].value)/np.interp(7500., wave[wave_gpm], flux[wave_gpm])

plt.plot(wave[wave_gpm], renorm*flux[wave_gpm], drawstyle='steps-mid', linewidth=0.5, color='black',
         label='{} fluxed with {}'.format(star, other_star))
plt.plot(std_dict['wave'], std_dict['flux'], 'green', label='True spectrum of {}'.format(star))

wave_min = wave[wave_gpm].min()
wave_max = wave[wave_gpm].max()
wave_std_gpm = (std_dict['wave'].value >= wave_min) & (std_dict['wave'].value <= wave_max)
ymax = 1.2*std_dict['flux'][wave_std_gpm].value.max()
plt.ylim(-0.1, ymax)
plt.xlim(0.98*wave_min, 1.02*wave_max)
plt.legend()
plt.show()