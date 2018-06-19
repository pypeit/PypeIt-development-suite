
""" For the development and testing of Extraction
"""


import copy
from astropy.table import Table
from pydl.pydlutils.trace import TraceSet
from pydl.pydlutils.image import djs_maskinterp

#from pypit.idl_stats import djs_iterstat
from pypit import ginga
#from pypit.extract_boxcar import extract_boxcar
from matplotlib import pyplot as plt
from astropy.stats import sigma_clipped_stats

from pydl.pydlutils.math import djs_median

import numpy as np
from astropy.io import fits
import glob

from pypit import msgs
from pypit import ardebug as debugger
from pypit import ginga
#from pypit import arload
#from pypit import arproc
#from pypit import arcomb
#from pypit import ardeimos
#from pypit import arlris
#from pypit import arpixels
#from pypit import arsave
#from pypit import traceslits

debug = debugger.init()
debug['develop'] = True
msgs.reset(debug=debug, verbosity=2)
import sys
from scipy.io import readsav
from scipy.interpolate import RectBivariateSpline



from astropy.table import Table
from pydl.pydlutils.trace import TraceSet
from pydl.pydlutils.image import djs_maskinterp
from scipy.special import erfcinv

#from pypit.idl_stats import djs_iterstat
from pypit import ginga
#from pypit.extract_boxcar import extract_boxcar
from matplotlib import pyplot as plt
from astropy.stats import sigma_clipped_stats

from pydl.pydlutils.math import djs_median
from pydl.pydlutils.bspline import iterfit as bspline_iterfit
from IPython import embed
from pydl.pydlutils.bspline import bspline

from pydl.goddard.math import flegendre
from scipy.special import ndtr

#


savefile='/Users/joe/gprofile_develop/objfind.sav'
idl_dict=readsav(savefile)
sciimg = idl_dict['sciimg']
sciivar = idl_dict['sciivar']
skyimage = idl_dict['skyimage']
#piximg = idl_dict['piximg']
#waveimg = idl_dict['waveimg']
ximg = idl_dict['ximg']
#objstruct = idl_dict['objstruct']
slitmask = idl_dict['slitmask']
xx1 = idl_dict['xx1']
xx2 = idl_dict['xx2']

slitid = 3 # work on slit number 3

slit_left = xx1[2,:]
slit_righ = xx2[2,:]

thismask = (slitmask == slitid)

# View things
idims = sciimg.shape
nspat = idims[1]
nspec = idims[0]
viewer, ch = ginga.show_image((sciimg - skyimage) * (slitmask == slitid))
ginga.show_slits(viewer, ch, np.reshape(slit_left,(nspec,1)), np.reshape(slit_righ, (nspec,1)),[slitid])


# Call to objfind
# def objfind(sciimg, thismask, slit_left, slit_righ,  sciivar = None, fwhm=3.0 )

#sciimg = science image
# thismask is a mask which is true on the slit and false outside it
# slit_left and slit_right are nspec ndarrays for the left and right slit edge
# return (specobjs, skymask, objmask)
ncoeff = 5 # Number of coefficients to fit (i.e. order of polynomial)
nperslit = 10 # Number of objects to find
peakthresh = 0.0 #Flux threshhold for finding objects; the flux must be at least this fraction of the brightest object in each slit. [default = 0].
absthresh = None #Absolute flux threshold for finding objects; the peakflux must be at least htis large; default to 0. If both peakthresh and
# absthresh are set, absthresh overrides peakthresh.
PEAK_SMTH = 5.0 #  Let's rename this to NFWHM and think about whether it could be deprecated
SIG_THRESH = 5.0  #  Sigma threshold for objects  [default=5.0]
fwhm = 3.0 # FWHM in pixels for convolving flux along the slit before peak-finding. [Default= 3.0]
OBJTHRESH = 0.5 # threshold for object masking
HAND_DICT = None # dictionary with entires HAND_SPAT, HAND_SPEC, and HAND_FWHM which are np.arrays for hand apertures
# Deprecated sky_fwhm, crude, sigma_sky, simple_sub is only used in the nirspec pipeline?

# Approximately compute the noise image if one is not given
if sciivar is None:
    sciivar = 1.0/(np.abs(sciimg) + 100.0)

idims = sciimg.shape
nspat = idims[1]
nspec = idims[0]
specmid = nspec/2

skymask = thismask
objmask = np.zeros_like(thismask)

if HAND_DICT is not None:
    



