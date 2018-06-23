
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
from astropy.stats import sigma_clipped_stats
import glob
from pypit.arpixels import core_slit_pixels, ximg_and_edgemask
from pypit.arutils import find_nminima
from pypit.core.arextract import extract_asymbox2, extract_boxcar

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
import scipy
#from scipy.interpolate import RectBivariateSpline
#from scipy.signal import medfilt


from pydl.pydlutils.trace import TraceSet
from pydl.pydlutils.image import djs_maskinterp


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
#ximg = idl_dict['ximg']
#objstruct = idl_dict['objstruct']
slitmask = idl_dict['slitmask']
xx1 = idl_dict['xx1']
xx2 = idl_dict['xx2']

slitid = 3 # work on slit number 3

slit_left = xx1[2,:]
slit_righ = xx2[2,:]

#thismask = (slitmask == slitid)

# View things
idims = sciimg.shape
nspat = idims[1]
nspec = idims[0]
viewer, ch = ginga.show_image((sciimg - skyimage) * (slitmask == slitid))
ginga.show_slits(viewer, ch, np.reshape(slit_left,(nspec,1)), np.reshape(slit_righ, (nspec,1)),[slitid])

# testing boxcar extraction
#savefile='/Users/joe/gprofile_develop/local_skysub_dev.sav'
#idl_skysub=readsav(savefile)
#objstruct = idl_skysub['objstruct']
#trace = np.zeros((nspec,2))
#trace[:,0] = objstruct[0]['xpos']
#trace[:,1] = objstruct[1]['xpos']

#box_rad=7.0
#flux  = extract_boxcar(sciimg-skyimage,trace,box_rad)

#

# Call to objfind, objfind operates on a single slit. It will be called by a code which loops over the slits, or that
# will occur in the driver routine



    #handslits = thismask[int(np.rint(HAND_SPEC)), int(np.rint(HAND_SPAT))]

# def objfind(image, slit_left, slit_righ, mask = None, thismask=None, ximg = None, fwhm=3.0, HAND_DICT = None)

mask = (sciivar > 0.0)
image = sciimg - skyimage
ximg = None
thismask = None
HAND_DICT={'HAND_SPEC':1024.0,'HAND_SPAT': 750.0, 'HAND_FWHM': 4.0}

#sciimg = science image
# thismask is a mask which is true on the slit and false outside it
# slit_left and slit_right are nspec ndarrays for the left and right slit edge
# return (specobjs, skymask, objmask)
ncoeff = 5 # Number of coefficients to fit (i.e. order of polynomial)
nperslit = 10 # Number of objects to find
peakthresh = 0.0 #Flux threshhold for finding objects; the flux must be at least this fraction of the brightest object in each slit. [default = 0].
absthresh = None #Absolute flux threshold for finding objects; the peakflux must be at least htis large; default to 0. If both peakthresh and
# absthresh are set, absthresh overrides peakthresh.
BG_SMTH = 5.0 #  Smoothing parameter for b/g subtracting the smashed object peak image. Could this be omitted as an input parameter?
PKWDTH = 3.0 # Width of peaks for find_nminima
SIG_THRESH = 5.0  #  Sigma threshold for objects  [default=5.0]
FWHM = 3.0 # FWHM in pixels for convolving flux along the slit before peak-finding. [Default= 3.0]
OBJTHRESH = 0.5 # threshold for object masking
HAND_DICT = None # dictionary with entires HAND_SPAT, HAND_SPEC, and HAND_FWHM which are np.arrays for hand apertures
# Deprecated sky_fwhm, crude, sigma_sky, simple_sub is only used in the nirspec pipeline?

frameshape = sciimg.shape
nspec = frameshape[0]
nspat = frameshape[1]

# Synthesize necessary inputs if user did not provide them as optional arguments
if thismask is None:
    pad =0
    slitpix = core_slit_pixels(slit_left, slit_righ, frameshape, pad)
    thismask = (slitpix > 0)

if ximg is None:
    ximg, edgmask = ximg_and_edgemask(slit_left, slit_righ, thismask)

if mask is None:
    mask = np.ones_like(thismask)


idims = sciimg.shape
nspat = idims[1]
nspec = idims[0]
specmid = nspec/2

skymask = thismask & (edgmask == False)
objmask = np.zeros_like(thismask)

thisimg =image*thismask*mask
#  Smash the image (for this slit) into a single flux vector.  How many pixels wide is the slit at each Y?
xsize = slit_righ - slit_left
nsamp = np.ceil(np.median(xsize))
# Mask skypixels with 2 FWHM of edge
left_asym = np.outer(slit_left,np.ones(int(nsamp))) + np.outer(xsize/nsamp, np.arange(nsamp))
righ_asym = left_asym + np.outer(xsize/nsamp, np.ones(int(nsamp)))
flux_spec = extract_asymbox2(thisimg, left_asym, righ_asym)
flux_mean, flux_median, flux_sig = sigma_clipped_stats(flux_spec,axis=0, sigma = 4.0)
if (nsamp < 9.0*FWHM):
    fluxsub = flux_mean.data - np.median(flux_mean.data)
else:
    kernel_size= int(np.ceil(BG_SMTH*FWHM) // 2 * 2 + 1) # This ensure kernel_size is odd
    fluxsub = flux_mean.data - scipy.signal.medfilt(flux_mean.data, kernel_size=kernel_size)
    # This little bit below deals with degenerate cases for which the slit gets brighter toward the edge, i.e. when
    # alignment stars saturate and bleed over into other slits. In this case the median smoothed profile is the nearly
    # everywhere the same as the profile itself, and fluxsub is full of zeros (bad!). If 90% or more of fluxsub is zero,
    # default to use the unfiltered case
    isub_bad = (fluxsub == 0.0)
    frac_bad = np.sum(isub_bad)/nsamp
    if frac_bad > 0.9:
        fluxsub = flux_mean.data - np.median(flux_mean.data)

    fluxconv = scipy.ndimage.filters.gaussian_filter1d(fluxsub, FWHM/2.3548,mode='reflect')
    xcen, sigma, ledg, redg = find_nminima(-fluxconv,nfind=nperslit, width = PKWDTH, minsep = np.fmax(FWHM, PKWDTH))
    ypeak = np.interp(xcen,np.arange(nsamp),fluxconv)

# HAND_DICT parsing
if HAND_DICT is not None:
    if ('HAND_SPEC' not in HAND_DICT.keys() | 'HAND_SPAT' not in HAND_DICT.keys()):
        raise ValueError('HAND_SPEC and HAND_SPAT must be set in the HAND_DICT')

    HAND_SPEC=np.asarray(HAND_DICT['HAND_SPEC'])
    HAND_SPAT=np.asarray(HAND_DICT['HAND_SPAT'])
    if(HAND_SPEC.size != HAND_SPAT.size):
        raise ValueError('HAND_SPEC and HAND_SPAT must have the same size in the HAND_DICT')
    nhand = HAND_SPEC.size
    HAND_FLAG = np.ones_like(HAND_SPEC,dtype=bool)

    HAND_FWHM = HAND_DICT.get('HAND_FWHM')
    if HAND_FWHM is not None:
        HAND_FWHM = np.asarray(HAND_FWHM)
        if(HAND_FWHM.size==HAND_SPEC.size):
            pass
        elif (HAND_FWHM.size == 1):
            HAND_FWHM = np.full(nhand, HAND_FWHM)
        else:
            raise ValueError('HAND_FWHM must either be a number of have the same size as HAND_SPEC and HAND_SPAT')


