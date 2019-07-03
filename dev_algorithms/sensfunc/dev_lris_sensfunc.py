import os
import numpy as np

import telluric
from flux1d import apply_sensfunc
from pypeit.core import coadd1d
from pypeit import msgs
from pypeit.core import load

debug = False
show = True
show_exp=True
do_sens = False


datapath = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT/Keck_LRIS_red/long_600_7500_d560/Science/')
objid ='SPAT0924-SLIT0000-DET01'
std1dfile = os.path.join(datapath, 'spec1d_LR.20160216.17613-G1910B2B_LRISr_2016Feb16T045333.331.fits')

sensfile = 'G191B2B_sens_tell_nires.fits'
telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/TelFit_MaunaKea_3100_26100_R20000.fits')
TelSens = telluric.sensfunc_telluric(std1dfile, telgridfile, sensfile, mask_abs_lines=True, debug=True)


