


import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from sklearn import mixture
from astropy.io import fits
from pypeit.core import pydl
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from pypeit.core import load
from pypeit.core import save
from pypeit.core import coadd2d
from pypeit.core import coadd1d
from pypeit.spectrographs import util
from pypeit import utils
from pypeit import msgs
import pickle
PYPEIT_FLUX_SCALE = 1e-17
from astropy.io import fits
import copy
import telluric
import IPython


dev_path = os.getenv('PYPEIT_DEV')
do_sens = True
do_qso = False


##################################################
# Example of sensitivity function determination  #
##################################################
if do_sens:
    spec1dfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/spec1d_XSHOO.2017-11-23T08:25:54.754-LTT3218_XShooter_NIR_2017Nov23T082554.754.fits')
    telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
    outfile = 'LTT3218_sens_tell.fits'
    TelSens = telluric.sensfunc_telluric(spec1dfile, telgridfile, outfile, mask_abs_lines=True,debug=False)


###############################
# Example of quasar model fit #
###############################
if do_qso:
    spec1dfluxfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/NIR_Stack/spec1d_stack_J0224-4711.fits')
    pca_file = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux//qso_pca_1200_3100.pckl')
    telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
    telloutfile = 'J0224-4711_tell.fits'
    outfile = 'spec1d_stack_tell_J0224-4711.fits'
    z_qso = 6.51
    TelQSO = telluric.qso_telluric(spec1dfluxfile, telgridfile, pca_file, z_qso, telloutfile, outfile,
                                   create_bal_mask=None, debug=True, show=True)



