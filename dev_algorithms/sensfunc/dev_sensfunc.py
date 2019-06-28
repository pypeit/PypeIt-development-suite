


import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from sklearn import mixture
from astropy.io import fits
from pypeit.core import pydl
import astropy.units as u
from astropy.io import fits
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

spec1dfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0020-3653/NIR/Science/spec1d_STD,FLUX_XShooter_NIR_2017Dec17T082243.751.fits')
telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
outfile = 'LTT3218_sens_tell.fits'
TelObj = telluric.sensfunc_telluric(spec1dfile, telgridfile, outfile, mask_abs_lines=True, only_orders=[5], debug=True)

spec1dfluxfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/NIR_Stack/spec1d_stack_J0224-4711.fits')
sys.exit(-1)


