


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
from pypeit.core import flux
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
do_sens = False
do_qso = True


##################################################
# Example of sensitivity function determination  #
##################################################
if do_sens:
    spec1dfile = '/d2/Feige/Dropbox/OBS_DATA/NIRES/ut190519/Science/spec1d_s190519_0059-GD153_NIRES_2019May19T083811.995.fits'
    telgridfile = '/d2/Feige/Dropbox/TellGrid/TelFit_MaunaKea_3100_26100_R20000.fits'
    outfile = 'GD153_sens_tell_nires.fits'
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



