import os

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy.time import Time

from pypeit import msgs
from pypeit import coadd1d
from pypeit import specobjs
from pypeit import inputfiles
from pypeit.core import meta
from pypeit.par import pypeitpar
from pypeit.scripts import scriptbase
from pypeit.spectrographs.util import load_spectrograph
from pypeit import sensfunc
from pypeit.core import coadd, flux_calib
from matplotlib import pyplot as plt
from pypeit.utils import fast_running_median


spectrograph = load_spectrograph('keck_hires')
par = spectrograph.default_pypeit_par()
redux_path = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/keck_hires/J0100+2802_RED_C1_ECH_-0.91_XD_1.46_1x2/'

star = 'BD+28'
if star == 'BD+28':
    sensfile = os.path.join(redux_path, 'sens_HI.20151102.54663-G191-B2B_HIRES_20151102T151104.762.fits')
    scifiles = ['spec1d_HI.20151102.15862-BD+28d4211_HIRES_20151102T042424.163.fits']
    objids = ['OBJ0479-MSC01']
elif star == 'G191':
    sensfile = os.path.join(redux_path, 'sens_HI.20151102.15862-BD+28d4211_HIRES_20151102T042424.163.fits')
    scifiles = ['spec1d_HI.20151102.54663-G191-B2B_HIRES_20151102T151104.762.fits']
    objids = ['OBJ0472-MSC01']

filenames = [os.path.join(redux_path, 'Science', f) for f in scifiles]


# Instantiate the coadd1d object
coadd1d_hires = coadd1d.CoAdd1D.get_instance(filenames, objids, spectrograph=spectrograph, par=par['coadd1d'],
                                           sensfile=sensfile, debug=True, show=True)
iexp = 0
sobjs = specobjs.SpecObjs.from_fitsfile(filenames[iexp])

indx = sobjs.name_indices(objids[iexp])
if not np.any(indx):
    msgs.error("No matching objects for {:s}.  Odds are you input the wrong OBJID".format(objids[iexp]))
wave_iexp, flux_iexp, ivar_iexp, gpm_iexp, trace_spec, trace_spat, meta_spec, header = \
    sobjs[indx].unpack_object(ret_flam=par['coadd1d']['flux_value'], extract_type=par['coadd1d']['ex_value'])
#

# If the user provided RA and DEC use those instead of what is in meta
star_ra = meta_spec['RA']
star_dec = meta_spec['DEC']
# Convert to decimal deg, as needed
star_ra, star_dec = meta.convert_radec(star_ra, star_dec)

# Read in standard star dictionary
std_dict = flux_calib.get_standard_spectrum(star_type=None, star_mag=None, ra=star_ra, dec=star_dec)

waves, fluxes, ivars, gpms, header_out = coadd1d_hires.load_arrays()

iexp = 0
norders = waves.shape[1]
nspec = waves.shape[0]
nsmooth = 101
colors = plt.get_cmap('tab10').colors
orders = meta_spec['ECH_ORDERS'] #61 - np.arange(norders)

ymin_max = np.zeros((norders, 2))
for ii, iord in enumerate(orders):
    rr = (np.max(iord) - iord)/np.maximum(np.max(orders) - np.min(orders), 1)
    gg = 0.0
    bb = (iord- np.min(orders))/np.maximum(np.max(orders) - np.min(orders), 1)

    wave_gpm =  waves[:, ii, iexp] > 1.0
    flux_sm = fast_running_median(fluxes[wave_gpm, ii, iexp], nsmooth)
    ymin_max[ii, 0] = np.min(flux_sm)
    ymin_max[ii, 1] = np.max(flux_sm)
    plt.plot(waves[wave_gpm, ii, iexp], fluxes[wave_gpm, ii, iexp], drawstyle='steps-mid',
             linewidth=0.5, color=(rr, gg, bb),
             label='Order {:d}'.format(iord))

plt.plot(std_dict['wave'], std_dict['flux'], 'green', label='Standard')
#ymin = np.min(ymin_max[:, 0])
#ymax = np.max(ymin_max[:, 1])
wave_min = waves[waves > 1.0].min()
wave_max = waves[waves > 1.0].max()
wave_std_gpm = (std_dict['wave'].value >= wave_min) & (std_dict['wave'].value <= wave_max)
ymax = 1.2*std_dict['flux'][wave_std_gpm].value.max()
plt.ylim(-0.1, ymax)
plt.xlim(0.98*wave_min, 1.02*wave_max)
plt.show()