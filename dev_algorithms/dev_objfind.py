
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

from pypit.core.artraceslits import trace_fweight, trace_gweight

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


from pydl.pydlutils.trace import traceset2xy, xy2traceset
from pydl.pydlutils.image import djs_maskinterp

from specobj import SpecObj
#from pypit.idl_stats import djs_iterstat
from pypit import ginga
#from pypit.extract_boxcar import extract_boxcar
from matplotlib import pyplot as plt
from astropy.stats import sigma_clipped_stats

from pydl.pydlutils.math import djs_median
from pydl.pydlutils.bspline import iterfit as bspline_iterfit
from IPython import embed
from pydl.pydlutils.trace import

from pydl.goddard.math import flegendre






def parse_hand_dict(HAND_DICT):
    """ Utility routine to parase the HAND_DICT dictionary for hand selected apertures

    Parameters
    ----------
    HAND_DICT:   dictionary

    Returns
    -------
    hand_spec:  spectral pixel location, numpy float 1-d array with size equal to number of hand aperatures requested
    hand_spat:  spatial pixel location, numpy float 1-d array with size equal to number of hand aperatures requested
    hand_fwhm:  hand aperture fwhm, numpy float 1-d array with size equal to number of hand aperatures requested

    Revision History
    ----------------
    23-June-2018  Written by J. Hennawi
    """


    if ('HAND_SPEC' not in HAND_DICT.keys() | 'HAND_SPAT' not in HAND_DICT.keys()):
        raise ValueError('HAND_SPEC and HAND_SPAT must be set in the HAND_DICT')

    HAND_SPEC=np.asarray(HAND_DICT['HAND_SPEC'])
    HAND_SPAT=np.asarray(HAND_DICT['HAND_SPAT'])
    HAND_DET = np.asarray(HAND_DICT['HAND_DET'])
    if(HAND_SPEC.size == HAND_SPAT.size == HAND_DET.size) == False:
        raise ValueError('HAND_SPEC, HAND_SPAT, and HAND_DET must have the same size in the HAND_DICT')
    nhand = HAND_SPEC.size

    HAND_FWHM = HAND_DICT.get('HAND_FWHM')
    if HAND_FWHM is not None:
        HAND_FWHM = np.asarray(HAND_FWHM)
        if(HAND_FWHM.size==HAND_SPEC.size):
            pass
        elif (HAND_FWHM.size == 1):
            HAND_FWHM = np.full(nhand, HAND_FWHM)
        else:
            raise ValueError('HAND_FWHM must either be a number of have the same size as HAND_SPEC and HAND_SPAT')
    else:
        HAND_FWHM = np.full(nhand, None)

    return (HAND_SPEC, HAND_SPAT, HAND_DET, HAND_FWHM)


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

# def objfind(image, slit_left, slit_righ, mask = None, thismask=None, ximg = None, edgmask = None, fwhm=3.0, HAND_DICT = None,
# config = None, slitid = None, scidx = 1, det = 1, objtype = 'science', std_trace = None, CHK = True)

# Things that will presumably be in some settings import which are currently arguments to specobj
config = 'lris_b1200'
scidx = 1
det = 1
objtype= 'science'
std_trace = None
slitid = None
CHK = True

# arguments provided to the code
mask = (sciivar > 0.0)
image = sciimg - skyimage
ximg = None
thismask = None
HAND_DICT={'HAND_SPEC':1024.0,'HAND_SPAT': 750.0, 'HAND_DET':1, 'HAND_FWHM': 4.0}

#sciimg = science image
# thismask is a mask which is true on the slit and false outside it
# slit_left and slit_right are nspec ndarrays for the left and right slit edge
# return (specobjs, skymask, objmask)
ncoeff = 5 # Number of coefficients to fit (i.e. order of polynomial)
nperslit = 10 # Number of objects to find
SIG_THRESH = 5.0  #  Sigma threshold for objects  [default=5.0].  # This is the significance relative to fluctuations
# in the sprectral-direction smashed image used for object finding
PEAK_THRESH = 0.0 #Add itional flux threshhold criter for finding objects; the flux must be at least this fraction of the brightest object in each slit. [default = 0].
# Must be between 0.0 and 1.0
ABS_THRESH = 0.0 #Absolute flux threshold for finding objects; the peakflux must be at least htis large; default to 0. If both peakthresh and
# abs_thresh are set, abs_thresh overrides peakthresh.
BG_SMTH = 5.0 #  Smoothing parameter for b/g subtracting the smashed object peak image. Could this be omitted as an input parameter?
PKWDTH = 3.0 # Width of peaks for find_nminima
FWHM = 3.0 # FWHM in pixels for convolving flux along the slit before peak-finding. [Default= 3.0]
OBJTHRESH = 0.5 # threshold for object masking
#HAND_DICT = None # dictionary with entires HAND_SPAT, HAND_SPEC, and HAND_FWHM which are np.arrays for hand apertures
# Deprecated sky_fwhm, crude, sigma_sky, simple_sub is only used in the nirspec pipeline?


if ((PEAK_THRESH >=0.0) & (PEAK_THRESH <=1.0)) == False:
    raise ValueError('Invalid value of PEAK_THRESH. It must be between 0.0 and 1.0')

frameshape = sciimg.shape
nspec = frameshape[0]
nspat = frameshape[1]

# Some information about this slit we need for later when we instantiate specobj objects
spec_vec = np.arange(nspec)
spat_vec = np.arange(nspat)
slit_spec_pos = nspec/2.0

slit_spat_pos = (np.interp(slit_spec_pos, spec_vec, slit_left), np.interp(slit_spec_pos, spec_vec, slit_righ))

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

#skymask = thismask & (edgmask == False)
objmask = np.zeros_like(thismask)

thisimg =image*thismask*mask
#  Smash the image (for this slit) into a single flux vector.  How many pixels wide is the slit at each Y?
xsize = slit_righ - slit_left
nsamp = np.ceil(np.median(xsize))
# Mask skypixels with 2 FWHM of edge
left_asym = np.outer(slit_left,np.ones(int(nsamp))) + np.outer(xsize/nsamp, np.arange(nsamp))
righ_asym = left_asym + np.outer(xsize/nsamp, np.ones(int(nsamp)))
# This extract_asymbox2 call smashes the image in the spectral direction along the curved object traces
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
# Create a mask for pixels to use for a background flucutation level estimate. Mask spatial pixels that hit an object
imask = np.ones(int(nsamp), dtype=bool)
xvec = np.arange(nsamp)
for zz in range(len(xcen)):
    ibad = (np.abs(xvec - xcen[zz]) <= 2.0*FWHM)
    imask[ibad] = False

# Good pixels for flucutation level estimate. Omit edge pixels and pixels within a FWHM of a candidate object
igd = imask & (xvec > 3.0) & (xvec <= (nsamp-3.0))
if np.any(igd) == False:
    igd = np.ones(int(nsamp),dtype=bool) # if all pixels are masked for some reason, don't mask

(mean, med_sn2, skythresh) = sigma_clipped_stats(fluxconv[igd], sigma=1.5)
(mean, med_sn2, sigma)     = sigma_clipped_stats(fluxconv[igd], sigma=2.5)
if(skythresh == 0.0) & (sigma != 0.0):
    skythresh = sigma
elif(skythresh == 0.0) & (sigma==0.0):  # if both SKYTHRESH and sigma are zero mask out the zero pixels and reavaluate
    good = fluxconv > 0.0
    if np.any(good) == True:
        (mean, med_sn2, skythresh) = sigma_clipped_stats(fluxconv[good], sigma=1.5)
        (mean, med_sn2, sigma) = sigma_clipped_stats(fluxconv[good], sigma=2.5)
    else:
        raise ValueError('Object finding failed. All the elements of the fluxconv spatial profile array are zero')
else:
    pass

# Get rid of peaks within 3% of slit edge which are almost always spurious
not_near_edge = (xcen > nsamp*0.03) & (xcen < nsamp*0.97)
xcen = xcen[not_near_edge]
ypeak = ypeak[not_near_edge]
npeak = len(xcen)

specobjs =[]

# Choose which ones to keep and discard based on threshold params. Create SpecObj objects
if npeak > 0:
    # Possible thresholds    [significance,  fraction of brightest, absolute]
    threshvec = np.array([SIG_THRESH*sigma, PEAK_THRESH*ypeak.max(), ABS_THRESH])
    threshold = threshvec.max()
    msgs.info('Using object finding threshold of: {:5.2f}'.format(threshold))
    # Trim to only objects above this threshold
    ikeep = (ypeak >= threshold)
    xcen = xcen[ikeep]
    ypeak = ypeak[ikeep]
    nobj_reg = len(xcen)
    # Now create SpecObj objects for all of these
    for iobj in range(nobj_reg):
        specobj = SpecObj(frameshape, slit_spat_pos, slit_spec_pos, det = det, config = config, slitid = slitid, scidx = scidx, objtype=objtype)
        specobj.spat_fracpos = xcen[iobj]/nsamp
        specobj.smash_peakflux = ypeak[iobj]
        specobjs.append(specobj)


# Now deal with the hand apertures if a HAND_DICT was passed in. Add these to the SpecObj objects
if HAND_DICT is not None:
    # First Parse the hand_dict
    HAND_SPEC, HAND_SPAT, HAND_DET, HAND_FWHM = parse_hand_dict(HAND_DICT)
    # Determine if these hand apertures land on the slit in question
    hand_on_slit = thismask[int(np.rint(HAND_SPEC)),int(np.rint(HAND_SPAT))]
    HAND_SPEC = HAND_SPEC[hand_on_slit]
    HAND_SPAT = HAND_SPAT[hand_on_slit]
    HAND_DET  = HAND_DET[hand_on_slit]
    HAND_FWHM = HAND_FWHM[hand_on_slit]
    nobj_hand = len(HAND_SPEC)
    for iobj in range(nobj_hand):
        specobj = SpecObj(frameshape, slit_spat_pos, slit_spec_pos, det = det, config = config, slitid = slitid, scidx = scidx, objtype=objtype)
        specobj.HAND_SPEC = HAND_SPEC[iobj]
        specobj.HAND_SPAT = HAND_SPAT[iobj]
        specobj.HAND_DET = HAND_DET[iobj]
        specobj.HAND_FWHM = HAND_FWHM[iobj]
        specobj.HAND_FLAG = True
        f_ximg = scipy.interpolate.RectBivariateSpline(spec_vec, spat_vec, ximg)
        specobj.spat_fracpos = f_ximg(specobj.HAND_SPEC, specobj.HAND_SPAT, grid=False) # interpolate from ximg
        specobj.smash_peakflux = np.interp(specobj.spat_fracpos*nsamp,np.arange(nsamp),fluxconv) # interpolate from fluxconv
        specobjs.append(specobj)


nobj = nobj_reg + nobj_hand
# Now loop over all apertures and assign preliminary traces to them. For the regular objects these will be refined below.
# For the hand objects, we will also update these again later with brightest objects on the slit
for iobj in range(nobj):
    # Was a standard trace provided? If so, use that as a crutch.
    if std_trace is not None:
        msgs.info('Using input STANDARD star trace as crutch for object tracing'.format(threshold))
        x_trace = np.interp(specmid, spec_vec, std_trace)
        shift = slit_left + xsize*specobjs[iobj].spat_fracpos - x_trace
        specobjs[iobj].trace_spat = std_trace + shift
    else:    # If no standard is provided shift left slit boundary over to be initial trace
        specobjs[iobj].trace_spat = slit_left  + xsize*specobjs[iobj].spat_fracpos
    # Determine the FWHM max for the regular apertures
    if specobjs[iobj].HAND_FLAG == False:
        yhalf = 0.5*specobjs[iobj].smash_peakflux
        xpk = specobjs[iobj].spat_fracpos*nsamp
        x0 = int(np.rint(xpk))
        # TODO It seems we have two codes that do similar things, i.e. findfwhm in arextract.py. Could imagine having one
        # Find right location where smash profile croses yhalf
        if x0 < (int(nsamp)-1):
            ind_righ, = np.where(fluxconv[x0:] < yhalf)
            if len(ind_righ) > 0:
                i2 = ind_righ[0]
                if i2 == 0:
                    xrigh = None
                else:
                    xrigh_int = scipy.interpolate.interp1d(fluxconv[x0 + i2-1:x0 + i2 + 1], x0 + np.array([i2-1,i2],dtype=float),assume_sorted=False)
                    xrigh = xrigh_int([yhalf])[0]
            else:
                xrigh = None
        else:
            xrigh = None
        # Find left location where smash profile crosses yhalf
        if x0 > 0:
            ind_left, = np.where(fluxconv[0:np.fmin(x0+1,int(nsamp)-1)] < yhalf)
            if len(ind_left) > 0:
                i1 = (ind_left[::-1])[0]
                if i1 == (int(nsamp)-1):
                    xleft = None
                else:
                    xleft_int = scipy.interpolate.interp1d(fluxconv[i1:i1+2],np.array([i1,i1+1],dtype=float))
                    xleft = xleft_int([yhalf])[0]
            else:
                xleft = None
        else:
            xleft = None

        if (xleft is None) & (xrigh is None):
            fwhm_measure = None
        elif xrigh is None:
            fwhm_measure = 2.0*(xpk- xleft)
        elif xleft is None:
            fwhm_measure = 2.0*(xrigh - xpk)
        else:
            fwhm_measure = (xrigh - xleft)

        if fwhm_measure is not None:
            specobjs[iobj].fwhm = np.sqrt(np.fmax(fwhm_measure**2 - FWHM**2, (FWHM/2.0)**2)) # Set a floor of FWHM/2 on FWHM
        else:
            specobjs[iobj].fwhm = FWHM

if len(specobjs) == 0:
    print('Put a return statement here')
    # print_objects(
    # return None, None, None, None

msgs.info('Tweaking the object traces')


niter = 12
fwhm_vec = np.zeros(niter)
fwhm_vec[0:niter//3] = 1.3*FWHM
fwhm_vec[niter//3:2*niter//3] = 1.1*FWHM
fwhm_vec[2*niter//3:] = FWHM

xpos0 = np.stack([spec.trace_spat for spec in specobjs], axis=1)
# This xinvvar is used to handle truncated slits/orders so that they don't break the fits
xfit1 = xpos0
ypos = np.outer(spec_vec, np.ones(nobj))
yind = np.arange(nspec,dtype=int)
for iiter in range(niter):
    xpos1, xerr1 = trace_fweight(image*mask,xfit1, radius = fwhm_vec[iiter])
    # Get the indices of the current trace
    xinvvar = np.ones_like(xpos0) # invvar used only for masking, not for weighted fitting
    # Loop over centroids and mask out anything that has left the slit/order.
    xind = (np.fmax(np.fmin(np.rint(xpos1),nspat-1),0)).astype(int)
    for iobj in range(nobj):
        xinvvar[:,iobj] = thismask[yind, xind[:,iobj]].astype(float)
    # Mask out anything that left the image.
    off_image = (xpos1 < -0.2*nspat) | (xpos1 > 1.2*nspat)
    xinvvar[off_image] = 0.0
    pos_set1 = xy2traceset(ypos,xpos1, ncoeff=ncoeff,maxdev = 5.0,invvar = xinvvar)






# For hand apertures we presume there is not enough flux to measure a FWHM.  So we use:
# 1) the specified HAND_FWHM, if it was provided
# 2) if no HAND_FWHM was specified in the dict, or it was set to None we will give the hand apertures the median value of all the other objects
# on this slit.
# 3) If there are no other objects on the slit, we will use the value FWHM for the hand apertures

# Sort the apertures

# Make the objmask

# Make the skymask

CHK == False
# If requested display the model fits for this grouping
if CHK == True:
    viewer, ch = ginga.show_image(image*(thismask==True))
    ginga.show_slits(viewer, ch, slit_left.T, slit_righ.T, [slitid])
    # TODO figure out a way to overplot the pixels that were masked in red like as a scatter plot
    for iobj in range(nobj):
        ginga.show_trace(viewer, ch, specobjs[iobj].trace_spat, specobjs[iobj].idx, color='green')


