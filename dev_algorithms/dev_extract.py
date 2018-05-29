
""" For the development and testing of Extraction
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)



from astropy.table import Table
from pydl.pydlutils.trace import TraceSet
from pydl.pydlutils.image import djs_maskinterp

#from pypit.idl_stats import djs_iterstat
from pypit import ginga
#from pypit.extract_boxcar import extract_boxcar
from matplotlib import pyplot as plt
from astropy.stats import sigma_clipped_stats

from pydl.pydlutils.math import djs_median
from pydl.pydlutils.bspline import iterfit as bspline_iterfit

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





def extract_asymbox2(image,left,right,ycen = None,weight_image = None):
    """ Extract the total flux within a boxcar window at many positions.
    This routine will accept an asymmetric/variable window
    Traces are expected to run vertically to be consistent with other
    extract_  routines. Based on idlspec2d/spec2d/extract_asymbox2.pro

    Parameters
    ----------
    image :   numpy float 2-d array [nspec, nspat]
    left  :   Lower boundary of boxcar window (given as floating pt pixels) [nTrace,nspec]
    right     - Upper boundary of boxcar window (given as floating pt pixels) [nTrace,nspec]

    Optional Parameters
    -------------------
    ycen :    Y positions corresponding to "left" (expected as integers) [nTrace, nspec]
    weight_image:  Weights to be applied to image before boxcar [nspec, nspat]

    Returns
    -------
    fextract:   Extracted flux at positions specified by (left<-->right, ycen) [nTrace, nspec]

    Revision History
    ----------------
    24-Mar-1999  Written by David Schlegel, Princeton.
    17-Feb-2003  Written with slow IDL routine, S. Burles, MIT
    27-Jan-2004  Adopted to do asymmetric/varying boxcar
    22-Apr-2018  Ported to python by Joe Hennawi
    """

    dim = left.shape
    ndim = len(dim)
    npix = dim[1]
    if (ndim == 1):
        nTrace = 1
    else:
        nTrace = dim[0]

    if ycen == None:
        if ndim == 1:
            ycen = np.arange(npix, dtype='int')
        elif ndim == 2:
            ycen = np.outer(np.ones(nTrace, dtype='int'), np.arange(npix, dtype='int'), )
        else:
            raise ValueError('left is not 1 or 2 dimensional')

    if np.size(left) != np.size(ycen):
        raise ValueError('Number of elements in left and ycen must be equal')

    idims = image.shape
    nspat = idims[1]
    nspec = idims[0]

    maxwindow = np.max(right - left)
    tempx = np.int(maxwindow + 3.0)

    bigleft = np.outer(left[:], np.ones(tempx))
    bigright = np.outer(right[:], np.ones(tempx))
    spot = np.outer(np.ones(npix * nTrace), np.arange(tempx)) + bigleft - 1
    bigy = np.outer(ycen[:], np.ones(tempx, dtype='int'))

    fullspot = np.array(np.fmin(np.fmax(np.round(spot + 1) - 1, 0), nspat - 1), int)
    fracleft = np.fmax(np.fmin(fullspot - bigleft, 0.5), -0.5)
    fracright = np.fmax(np.fmin(bigright - fullspot, 0.5), -0.5)
    del bigleft
    del bigright
    bool_mask1 = (spot >= -0.5) & (spot < (nspat - 0.5))
    bool_mask2 = (bigy >= 0) & (bigy <= (nspec - 1))
    weight = (np.fmin(np.fmax(fracleft + fracright, 0), 1)) * bool_mask1 * bool_mask2
    del spot
    del fracleft
    del fracright
    bigy = np.fmin(np.fmax(bigy, 0), nspec - 1)

    if weight_image != None:
        temp = np.array([weight_image[x1, y1] * image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2 = np.reshape(weight.flatten() * temp, (nTrace, npix, tempx))
        fextract = np.sum(temp2, axis=2)
        temp_wi = np.array([weight_image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2_wi = np.reshape(weight.flatten() * temp_wi, (nTrace, npix, tempx))
        f_ivar = np.sum(temp2_wi, axis=2)
        fextract = fextract / (f_ivar + (f_ivar == 0)) * (f_ivar > 0)
    else:
        # Might be more pythonic way to code this. I needed to switch the flattening order in order to get
        # this to work
        temp = np.array([image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2 = np.reshape(weight.flatten() * temp, (nTrace, npix, tempx))
        fextract = np.sum(temp2, axis=2)

    # IDL version model functionality not implemented yet
    # At the moment I'm not reutnring the f_ivar for the weight_image mode. I'm not sure that this functionality is even
    # ever used

    return fextract

def extract_boxcar(image,trace, radius, ycen = None):
    """ Extract the total flux within a boxcar window at many positions. Based on idlspec2d/spec2d/extract_boxcar.pro

    Parameters
    ----------
    image :   numpy float 2-d array [nspec, nspat]
    trace :   Lower boundary of boxcar window (given as floating pt pixels) [nTrace,nspec]
    radius :  boxcar radius (given as floating pt pixels)

    Optional Parameters
    -------------------
    ycen :    Y positions corresponding to "trace" (expected as integers) [nTrace, nspec]

    Returns
    -------
    fextract:   Extracted flux at positions within (trace +- radius, ycen) [nTrace, nspec]

    Revision History
    ----------------
    24-Mar-1999  Written by David Schlegel, Princeton.
    22-Apr-2018  Ported to python by Joe Hennawi
    """

    dim = trace.shape
    ndim = len(dim)
    npix = dim[1]
    if (ndim == 1):
        nTrace = 1
    else:
        nTrace = dim[0]

    if ycen == None:
        if ndim == 1:
            ycen = np.arange(npix, dtype='int')
        elif ndim == 2:
            ycen = np.outer(np.ones(nTrace, dtype='int'), np.arange(npix, dtype='int'), )
        else:
            raise ValueError('trace is not 1 or 2 dimensional')

    if np.size(trace) != np.size(ycen):
        raise ValueError('Number of elements in trace and ycen must be equal')

    left = trace - radius
    right = trace + radius
    fextract = extract_asymbox2(image, left, right, ycen)

    return fextract


def findfwhm(model, sig_x):
    """ Calculate the spatial FWHM from an object profile.

    Parameters
    ----------
    model :   numpy float 2-d array [nspec, nspat]
    x :

    Returns
    -------
    peak :  Peak value of the profile model
    peak_x:  sig_x location where the peak value is obtained
    lwhm:   Value of sig_x at the left width at half maximum
    rwhm:   Value of sig_x at the right width at half maximum

    Revision History
    ----------------
    11-Mar-2005  Written by J. Hennawi and S. Burles David Schlegel, Princeton.
    28-May-2018  Ported to python by J. Hennawi
    """


    peak = (model*(np.abs(sig_x) < 1.)).max()
    peak_x = sig_x[(model*(np.abs(sig_x) < 1.)).argmax()]

    lrev = ((sig_x < peak_x) & (model < 0.5*peak))[::-1]
    lind, = np.where(lrev)
    if(lind.size > 0):
        lh = lind.min()
        lwhm = (sig_x[::-1])[lh]
    else:
        lwhm = -0.5*2.3548

    rind, = np.where((sig_x > peak_x) & (model < 0.5*peak))
    if(rind.size > 0):
        rh = rind.min()
        rwhm = sig_x[rh]
    else:
        rwhm = 0.5 * 2.3548

    return (peak, peak_x, lwhm, rwhm)






## Directory for IDL tests is /Users/joe/gprofile_develop/
## Run long_reduce, islit=3

path = '/Users/joe/REDUX/lris_redux/May_2008/0938+5317/Blue1200/'

# slitfile
slitfile = path + 'slits-lblue0059.fits.gz'
# science image
scifile = path + 'Science/sci-lblue0084.fits.gz'
wavefile = path + 'wave-lblue0028.fits.gz'



# Open up science output
hdu_sci = fits.open(scifile)
sciimg = hdu_sci[0].data
sciivar = hdu_sci[1].data
sky_model = hdu_sci[2].data
objstruct = Table.read(scifile,hdu=5)
# Grab the object trace
nobj = len(objstruct)
trace_in = objstruct['XPOS']
slit_str = np.array(objstruct['SLITID'],dtype=str)
obj_str = np.array(objstruct['OBJID'],dtype=str)
id_str = [x + '-' + y for x,y in zip(slit_str,obj_str)]
title_string = np.core.defchararray.add(np.array(objstruct['SLITID'],dtype=str),np.array(objstruct['OBJID'],dtype=str))
# Open  up the wavelength image
hdu_wave = fits.open(wavefile)
waveimg = hdu_wave[0].data


# Open up slitfile and read in tset_slits
tset_fits = Table.read(slitfile)
hdu_slit = fits.open(slitfile)
slitmask = hdu_slit[0].data
left_dum = Table(tset_fits[0])
righ_dum = Table(tset_fits[1])
left_dum.write('left_dum.fits',overwrite=True)
righ_dum.write('righ_dum.fits',overwrite=True)


hdu_left = fits.open('left_dum.fits')
fits_rec_left = hdu_left[1].data
hdu_righ = fits.open('righ_dum.fits')
fits_rec_righ = hdu_righ[1].data
#tset_fits = hdu[1].data
#tset_left = Table(ttset_fits[0],'dum_file.fits')
# Convert to python objects
tset_left = TraceSet(fits_rec_left)
tset_righ = TraceSet(fits_rec_righ)
# Do the initialization here
#
rows,left_edge = tset_left.xy()
rows,righ_edge = tset_righ.xy()
nslits = tset_left.nTrace

# Let's take a look at the image and the slits
#viewer, ch = ginga.show_image(sciimg-sky_model)
#ginga.show_slits(viewer, ch, left_edge.T, righ_edge.T,np.arange(nslits+1))
#for iobj in range(nobj):
#    ginga.show_trace(viewer, ch, trace_in[iobj,:],id_str[iobj], color='green')

# parameters for extract_asymbox2
box_rad = 7
trace_in = trace_in[5,:]
#right = trace_in[5:,:] + box_rad
image = sciimg - sky_model
ivar = sciivar*(slitmask == 3)
wave = objstruct['WAVE_BOX'][5]
flux = objstruct['FLUX_BOX'][5]
fluxivar = objstruct['IVAR_BOX'][5]
hwidth = objstruct['MASKWIDTH'][5]
thisfwhm = objstruct['FWHM'][5]

# This is the argument list
SN_GAUSS = None
MAX_TRACE_CORR = None
wvmnx = None
GAUSS = None # Do a Gaussian extraction
PROF_NSIGMA = None

if SN_GAUSS is None: SN_GAUSS = 3.0
if thisfwhm is None: thisfwhm = 4.0
if hwidth is None: 3.0*(np.max(thisfwhm) + 1.0)
if MAX_TRACE_CORR is None:  MAX_TRACE_CORR = 2.0
if wvmnx is None: wvmnx = [2900., 30000]

thisfwhm = np.fmax(thisfwhm,1.0) # require the FWHM to be greater than 1 pixel

xnew = trace_in
dims = image.shape
nspat = dims[1]
nspec = dims[0]

# Right now I'm not dealing with sub-images, but I may do so in the future
#top = np.fmin(np.max(np.where(np.sum(ivar == 0,0) < nspec)),nspat)
#bot = np.fmax(np.min(np.where(np.sum(ivar == 0,0) < nspec)),0)
#min_column = np.fmax(np.min(trace_in - hwidth)),bot)
#max_column = long(max(trace_in + hwidth)) <  top

# create some images we will need
profile_model = np.zeros((nspec,nspat))
sub_obj = image
sub_ivar = ivar
sub_wave = waveimg
sub_trace = trace_in
sub_x = np.arange(nspat)
sn2_sub = np.zeros((nspec,nspat))
spline_sub = np.zeros((nspec,nspat))



flux_sm = djs_median(flux, width = 5, boundary = 'reflect')
fluxivar_sm =  djs_median(fluxivar, width = 5, boundary = 'reflect')
fluxivar_sm = fluxivar_sm*(fluxivar > 0.0)

indsp = (wave > wvmnx[0]) & (wave < wvmnx[1]) & \
         np.isfinite(flux_sm) & (flux_sm < 5.0e5) &  \
         (flux_sm > -1000.0) & (fluxivar_sm > 0.0)
nsp = np.sum(indsp)

# Not all the djs_reject keywords args are implemented yet
b_answer, bmask   = bspline_iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp],everyn = 1.5)
b_answer, bmask2  = bspline_iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp]*bmask,everyn = 1.5)
c_answer, cmask   = bspline_iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp]*bmask2,everyn = 30)
spline_flux, _ = b_answer.value(wave[indsp])
cont_flux, _ = c_answer.value(wave[indsp])

sn2 = (np.fmax(spline_flux*(np.sqrt(np.fmax(fluxivar_sm[indsp], 0))*bmask2),0))**2
ind_nonzero = (sn2 > 0)
nonzero = np.sum(ind_nonzero)
if(nonzero >0):
    (mean, med_sn2, stddev) = sigma_clipped_stats(sn2)
else: med_sn2 = 0.0
sn2_med = djs_median(sn2, width = 9, boundary = 'reflect')
igood = (ivar > 0.0)
ngd = np.sum(igood)
if(ngd > 0):
    isrt = np.argsort(wave[indsp])
    sn2_sub[igood] = np.interp(sub_wave[igood],(wave[indsp])[isrt],sn2_med[isrt])
msgs.info('sqrt(med(S/N)^2) = ' + "{:5.2f}".format(np.sqrt(med_sn2)))

min_wave = np.min(wave[indsp])
max_wave = np.max(wave[indsp])
spline_flux1 = np.zeros(nspec)
cont_flux1 = np.zeros(nspec)
sn2_1 = np.zeros(nspec)
ispline = (wave >= min_wave) & (wave <= max_wave)
spline_tmp, _ = b_answer.value(wave[ispline])
spline_flux1[ispline] = spline_tmp
cont_tmp, _ = c_answer.value(wave[ispline])
cont_flux1[ispline] = cont_tmp
sn2_1[ispline] = np.interp(wave[ispline], (wave[indsp])[isrt], sn2[isrt])
bmask = np.zeros(nspec,dtype='bool')
bmask[indsp] = bmask2
spline_flux1 = djs_maskinterp(spline_flux1,(bmask == False))
cmask2 = np.zeros(nspec,dtype='bool')
cmask2[indsp] = cmask
cont_flux1 = djs_maskinterp(cont_flux1,(cmask2 == False))

if(med_sn2 <= 2.0):
    (_, _, sigma1) = sigma_clipped_stats(flux[indsp])
    spline_sub[igood]= np.fmax(sigma1,0)
else:
    if((med_sn2 <=5.0) and (med_sn2 > 2.0)):
        spline_flux1 = cont_flux1
    # Interp over points <= 0 in boxcar flux or masked points using cont model
    badpix = (spline_flux1 <= 0.5) | (bmask == False)
    goodval = (cont_flux1 > 0.0) & (cont_flux1 < 5e5)
    indbad1 = badpix & goodval
    nbad1 = np.sum(indbad1)
    if(nbad1 > 0):
        spline_flux1[indbad1] = cont_flux1[indbad1]
    indbad2 = badpix & ~goodval
    nbad2 = np.sum(indbad2)
    ngood0 = np.sum(~badpix)
    if((nbad2 > 0) or (ngood0 > 0)):
        spline_flux1[indbad2] = djs_median(spline_flux1[~badpix])
    # take a 5-pixel median to filter out some hot pixels
    spline_flux1 = djs_median(spline_flux1,width=5,boundary ='reflect')
    # Create the normalized object image
    if(ngd > 0):
        isrt = np.argsort(wave)
        spline_sub[igood] = np.interp(sub_wave[igood],wave[isrt],spline_flux1[isrt])


norm_obj = (spline_sub != 0.0)*sub_obj/(spline_sub + (spline_sub == 0.0))
norm_ivar = sub_ivar*spline_sub**2

# Cap very large inverse variances
ivar_mask = (norm_obj > -0.2) & (norm_obj < 0.7) & (sub_ivar > 0.0) & np.isfinite(norm_obj) & np.isfinite(norm_ivar)
norm_ivar = norm_ivar*ivar_mask
good = (norm_ivar > 0.0)
ngood = np.sum(good)


xtemp = np.cumsum(np.outer(4.0 + np.sqrt(np.fmax(sn2_1, 0.0)),np.ones(nspat)))
xtemp = xtemp/xtemp.max()

# norm_x is the x position along the image centered on the object  trace
norm_x = np.outer(np.ones(nspec), sub_x) - np.outer(sub_trace,np.ones(nspat))

sigma = np.full(nspat, thisfwhm/2.3548)
fwhmfit = sigma*2.3548
trace_corr = np.zeros(nspat)

# If we have too few pixels to fit a profile or S/N is too low, just use a Gaussian profile
# TODO Clean up the logic below. It is formally correct but
# redundant since no trace correction has been created or applied yet.  I'm only postponing doing it
# to preserve this syntax for later when it is needed
if((ngood < 10) or (med_sn2 < SN_GAUSS) or (GAUSS != None)):
    msgs.info("Too few good pixels or S/N <" + "{:5.1f}".format(SN_GAUSS) + " or GAUSS flag set")
    msgs.info("Returning Gaussian profile")
    sigma_x = norm_x/(np.outer(sigma, np.ones(nspat)) - np.outer(trace_corr,np.ones(nspat)))
    profile_model = np.exp(-0.5*sigma_x**2)/np.sqrt(2.0*np.pi)*(sigma_x**2 < 25.)
    msgs.info("FWHM="  + "{:6.2f}".format(thisfwhm) + ", S/N=" + "{:8.3f}".format(np.sqrt(med_sn2)))
    nxinf = np.sum(np.isfinite(xnew) == False)
    if(nxinf != 0):
        msgs.warn("Nan pixel values in trace correction")
        msgs.warn("Returning original trace....")
        xnew = trace_in
    inf = np.isfinite(profile_model) == False
    ninf = np.sum(inf)
    if (ninf != 0):
        msgs.warn("Nan pixel values in object profile... setting them to zero")
        profile_model[inf] = 0.0
    # Normalize profile
    norm = np.outer(np.sum(profile_model,1),np.ones(nspat))
    if(np.sum(norm) > 0.0):
        profile_model = profile_model/norm
    #TODO Put in call to QA here

    # TODO Put in a return statement here? Code returns updated trace and profile
    # return (xnew, profile_model)

msgs.info("Gaussian vs b-spline of width " + "{:6.2f}".format(thisfwhm) + " pixels")
area = 1.0
sigma_x = norm_x / (np.outer(sigma, np.ones(nspat)) - np.outer(trace_corr, np.ones(nspat)))

mask = np.full(nspec*nspat, False, dtype=bool)
skymask = np.full(nspec*nspat, True, dtype=bool)

# The following lines set the limits for the b-spline fit
limit = erfcinv(0.1/np.sqrt(med_sn2))*np.sqrt(2.0)
if(PROF_NSIGMA==None):
    sinh_space = 0.25*np.log10(np.fmax((1000./np.sqrt(med_sn2)),10.))
    abs_sigma = np.fmin((np.abs(sigma_x[good])).max(),2.0*limit)
    min_sigma = np.fmax(sigma_x[good].min(), (-abs_sigma))
    max_sigma = np.fmin(sigma_x[good].max(), (abs_sigma))
    nb = (np.arcsinh(abs_sigma)/sinh_space).astype(int) + 1
else:
    msgs.info("Using PROF_NSIGMA= " + "{:6.2f}".format(PROF_NSIGMA) + " for extended/bright objects")
    nb = np.round(PROF_NSIGMA > 10)
    max_sigma = PROF_NSIGMA
    min_sigma = -1*PROF_NSIGMA
    sinh_space = np.arcsinh(PROF_NSIGMA)/nb

rb = np.sinh((np.arange(nb) + 0.5) * sinh_space)
bkpt = np.concatenate([(-rb)[::-1], rb])
keep = ((bkpt >= min_sigma) & (bkpt <= max_sigma))
bkpt = bkpt[keep]

# Attempt B-spline first
GOOD_PIX = (sn2_sub > SN_GAUSS) & (norm_ivar > 0)
IN_PIX   = (sigma_x >= min_sigma) & (sigma_x <= max_sigma) & (norm_ivar > 0)
ngoodpix = np.sum(GOOD_PIX)
ninpix     = np.sum(IN_PIX)

if (ngoodpix >= 0.2*ninpix):
    inside, = np.where((GOOD_PIX & IN_PIX).flatten())
else:
    inside, = np.where(IN_PIX.flatten())


sigma_x_flat = sigma_x.flatten()
norm_obj_flat = norm_obj.flatten()
norm_ivar_flat = norm_ivar.flatten()


si = inside[np.argsort(sigma_x_flat[inside])]
sr = si[::-1]

bset, bmask = bspline_iterfit(sigma_x_flat[si],norm_obj_flat[si], invvar = norm_ivar_flat[si]
                       , nord = 4, bkpt = bkpt, maxiter = 15, upper = 1, lower = 1)
mode_fit, _ = bset.value(sigma_x_flat[si])
median_fit = np.median(norm_obj[norm_ivar > 0.0])

# TODO I don't follow the logic behind this statement but I'm leaving it for now. If the median is large it is used, otherwise we  user zero???
if (np.abs(median_fit) > 0.01):
    msgs.info("Median flux level in profile is not zero: median = " + "{:7.4f}".format(median_fit))
else:
    median_fit = 0.0

# Find the peak and FWHM based this profile fit
(peak, peak_x, lwhm, rwhm) = findfwhm(mode_fit - median_fit, sigma_x_flat[si])
trace_corr = np.full(nspec, peak_x)
min_level = peak*np.exp(-0.5*limit**2)

bspline_fwhm = (rwhm - lwhm)*thisfwhm/2.3548
msgs.info("Bspline FWHM: " + "{:7.4f}".format(bspline_fwhm) + ", compared to initial object finding FWHM: " + "{:7.4f}".format(thisfwhm) )
sigma = sigma * (rwhm-lwhm)/2.3548

limit = limit * (rwhm-lwhm)/2.3548

rev_fit = mode_fit[::-1]
lind, = np.where(((rev_fit < (min_level+median_fit)) & (sigma_x_flat[sr] < peak_x)) | (sigma_x_flat[sr] < (peak_x-limit)))
if (lind.size > 0):
    lp = lind.min()
    l_limit = sigma_x_flat[sr[lp]]
else:
    l_limit = min_sigma

rind, = np.where(((mode_fit < (min_level+median_fit)) & (sigma_x_flat[si] > peak_x)) | (sigma_x_flat[si] > (peak_x+limit)))
if (rind.size > 0):
    rp = rind.min()
    r_limit = sigma_x_flat[si[rp]]
else:
    r_limit = max_sigma

msgs.info("Trace limits: limit = " + "{:7.4f}".format(limit) + ", min_level = " + "{:7.4f}".format(min_level) +
          ", l_limit = " + "{:7.4f}".format(l_limit) + ", r_limit = " + "{:7.4f}".format(r_limit))

# Just grab the data points within the limits
mask[si]=((norm_ivar_flat[si] > 0) & (np.abs(norm_obj_flat[si] - mode_fit) < 0.1))
inside, = np.where((sigma_x_flat[si] > l_limit) & (sigma_x_flat[si] < r_limit) & mask[si])
ninside = inside.size


# If we have too few pixels after this step, then again just use a Gaussian profile and return. Note that
# we are following the original IDL code here and not using the improved sigma and trace correction for the
# profile at this stage since ninside is so small.
if(ninside < 10):
    msgs.info("Too few pixels inside l_limit and r_limit")
    msgs.info("Returning Gaussian profile")
    profile_model = np.exp(-0.5*sigma_x**2)/np.sqrt(2.0*np.pi)*(sigma_x**2 < 25.)
    msgs.info("FWHM="  + "{:6.2f}".format(thisfwhm) + ", S/N=" + "{:8.3f}".format(np.sqrt(med_sn2)))
    nxinf = np.sum(np.isfinite(xnew) == False)
    if(nxinf != 0):
        msgs.warn("Nan pixel values in trace correction")
        msgs.warn("Returning original trace....")
        xnew = trace_in
    inf = np.isfinite(profile_model) == False
    ninf = np.sum(inf)
    if (ninf != 0):
        msgs.warn("Nan pixel values in object profile... setting them to zero")
        profile_model[inf] = 0.0
    # Normalize profile
    norm = np.outer(np.sum(profile_model,1),np.ones(nspat))
    if(np.sum(norm) > 0.0):
        profile_model = profile_model/norm
    #TODO Put in call to QA here

    # TODO Put in a return statement here? Code returns updated trace and profile
    # return (xnew, profile_model)

sigma_iter = 3
isort =  (xtemp[si[inside]]).argsort()
inside = si[inside[isort]]
pb =np.ones(inside.size)

# Add the loop later
#for iiter in range(1,sigma_iter + 1):
mode_zero, _ = bset.value(sigma_x_flat[inside])
mode_zero = mode_zero*pb

mode_min05, _ = bset.value(sigma_x_flat[inside]-0.5)
mode_plu05, _ = bset.value(sigma_x_flat[inside]+0.5)
mode_shift = (mode_min05  - mode_plu05)*pb*((sigma_x_flat[inside] > (l_limit + 0.5)) &
                                            (sigma_x_flat[inside] < (r_limit - 0.5)))

mode_by13, _ = bset.value(sigma_x_flat[inside]/1.3)
mode_stretch = mode_by13*pb/1.3 - mode_zero

## Everything below here is bspline_longslit development
profile_basis = np.column_stack((mode_zero,mode_shift))
xdata = xtemp[inside]
ydata = norm_obj_flat[inside]
invvar = norm_ivar_flat[inside]
import pydl.pydlutils.bspline as bspline

maxiter =1
#def bspline_longslit(xdata, ydata, invvar, profile_basis, upper=5, lower=5,maxiter=10, **kwargs):
from pydl.pydlutils.math import djs_reject

nx = xdata.size
if ydata.size != nx:
    raise ValueError('Dimensions of xdata and ydata do not agree.')
if invvar.size != nx:
    raise ValueError('Dimensions of xdata and invvar do not agree.')

# Omit the model functionality right now
#if model != None:
#    model = image*0
#    modelwt = image*0



#ncol =


# flux  = extract_boxcar(image,trace,box_rad)
