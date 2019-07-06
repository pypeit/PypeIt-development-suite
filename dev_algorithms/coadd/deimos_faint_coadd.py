import os
#from pypeit.core import coadd1d
from pypeit.core import coadd1d
from pypeit.core import load
import numpy as np
import os
from pypeit import utils
from astropy.table import Table
from astropy.io import fits
from matplotlib import pyplot as plt
import sys
from IPython import embed



def deimos_fnames_1451():
    reduxpath_6700 = '/Users/joe/DEIMOS_redux/July19/1630A/6700/Science/'
    fnames_6700 = [ 'spec1d_d0705_0048-J1630A_DEIMOS_2019Jul05T065811.251.fits',
                    'spec1d_d0705_0049-J1630A_DEIMOS_2019Jul05T072921.898.fits',
                    'spec1d_d0705_0050-J1630A_DEIMOS_2019Jul05T080052.502.fits']
    objids_6700 = [ 'SPAT1453-SLIT0007-DET03',  'SPAT1452-SLIT0007-DET03', 'SPAT1451-SLIT0007-DET03']

    reduxpath_6720 = '/Users/joe/DEIMOS_redux/July19/1630A/6720/Science/'
    fnames_6720 = ['spec1d_d0705_0052-J1630A_DEIMOS_2019Jul05T090350.083.fits',
                    'spec1d_d0705_0053-J1630A_DEIMOS_2019Jul05T093500.730.fits',
                    'spec1d_d0705_0054-J1630A_DEIMOS_2019Jul05T100613.450.fits']
    objids_6720 = ['SPAT1451-SLIT0007-DET03', 'SPAT1451-SLIT0007-DET03', 'SPAT1452-SLIT0007-DET03']

    fnames = np.core.defchararray.add(reduxpath_6700, fnames_6700).tolist() + \
             np.core.defchararray.add(reduxpath_6720, fnames_6720).tolist()
    objids = objids_6700 + objids_6720


    return fnames, objids



def deimos_fnames_1433():
    reduxpath_6700 = '/Users/joe/DEIMOS_redux/July19/1630A/6700/Science/'
    fnames_6700 = [ 'spec1d_d0705_0048-J1630A_DEIMOS_2019Jul05T065811.251.fits',
                    'spec1d_d0705_0049-J1630A_DEIMOS_2019Jul05T072921.898.fits',
                    'spec1d_d0705_0050-J1630A_DEIMOS_2019Jul05T080052.502.fits']
    objids_6700 = [ 'SPAT1433-SLIT0007-DET02', 'SPAT1432-SLIT0007-DET02', 'SPAT1431-SLIT0007-DET02']

    reduxpath_6720 = '/Users/joe/DEIMOS_redux/July19/1630A/6720/Science/'
    fnames_6720 = ['spec1d_d0705_0051-J1630A_DEIMOS_2019Jul05T083238.832.fits',
                   'spec1d_d0705_0052-J1630A_DEIMOS_2019Jul05T090350.083.fits',
                   'spec1d_d0705_0053-J1630A_DEIMOS_2019Jul05T093500.730.fits',
                   'spec1d_d0705_0054-J1630A_DEIMOS_2019Jul05T100613.450.fits']
    objids_6720 = ['SPAT1432-SLIT0007-DET02']*4

    fnames = np.core.defchararray.add(reduxpath_6700, fnames_6700).tolist() + \
             np.core.defchararray.add(reduxpath_6720, fnames_6720).tolist()
    objids = objids_6700 + objids_6720


    return fnames, objids



fnames, objids = deimos_fnames_1433()
wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.multi_combspec(fnames, objids, flux_value=False,
                                                                        show=True, debug=True,
                                                                        debug_scale=True, show_scale=True,
                                                                        outfile='SLIT14_SPAT1433_24.42_coadd.fits')


