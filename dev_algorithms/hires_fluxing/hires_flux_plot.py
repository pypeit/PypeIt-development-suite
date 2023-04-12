import os

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy.time import Time

from pypeit import msgs
from pypeit import coadd1d
from pypeit import inputfiles
from pypeit.par import pypeitpar
from pypeit.scripts import scriptbase
from pypeit.spectrographs.util import load_spectrograph
from pypeit import sensfunc
from pypeit.core import coadd, flux_calib
from matplotlib import pyplot as plt
from pypeit.utils import fast_running_median


spectrograph = load_spectrograph('keck_hires')
par = spectrograph.default_pypeit_par()
redux_path = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/keck_hires/RED_C1_ECH_-0.82_XD_1.62/'
sensfile = os.path.join(redux_path, 'sens_HI.20151214.16715-Feige110_HIRES_20151214T043836.845.fits')
scifiles = ['spec1d_HI.20151214.17593-J0100+2802_HIRES_20151214T045314.323.fits',
            'spec1d_HI.20151214.20581-J0100+2802_HIRES_20151214T054302.726.fits']
filenames = [os.path.join(redux_path, 'Science/', f) for f in scifiles]

objids = ['OBJ0450-MSC01', 'OBJ0454-MSC01']


# Instantiate the coadd1d object
coadd1d_hires = coadd1d.CoAdd1D.get_instance(filenames, objids, spectrograph=spectrograph, par=par['coadd1d'],
                                           sensfile=sensfile, debug=True, show=True)

waves, fluxes, ivars, gpms, header_out = coadd1d_hires.load_arrays()

iexp = 0
norders = waves.shape[1]
nspec = waves.shape[0]
nsmooth = 101
colors = plt.get_cmap('tab10').colors
orders = 61 - np.arange(norders)

ymin_max = np.zeros((norders, 2))
for ii, iord in enumerate(orders):
    rr = (np.max(iord) - iord)/np.maximum(np.max(orders) - np.min(orders), 1)
    gg = 0.0
    bb = (iord- np.min(orders))/np.maximum(np.max(orders) - np.min(orders), 1)

    flux_sm = fast_running_median(fluxes[:, ii, iexp], nsmooth)
    ymin_max[ii, 0] = np.min(flux_sm)
    ymin_max[ii, 1] = np.max(flux_sm)
    plt.plot(waves[:, ii, iexp], fluxes[:, ii, iexp], drawstyle='steps-mid', color=(rr, gg, bb),
             label='Order {:d}'.format(iord))

plt.show()