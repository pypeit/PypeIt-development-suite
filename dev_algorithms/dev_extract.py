
""" For the development and testing of Extraction
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)



from astropy.io import fits
from astropy.table import Table
from pydl.pydlutils.trace import TraceSet
from pydl.pydlutils.image import djs_maskinterp

#from pypit.idl_stats import djs_iterstat
from pypit import ginga
from pypit.extract_boxcar import extract_boxcar
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
from pypit import arload
from pypit import arproc
from pypit import arcomb
from pypit import ardeimos
from pypit import arlris
from pypit import arpixels
from pypit import arsave
from pypit import traceslits

debug = debugger.init()
debug['develop'] = True
msgs.reset(debug=debug, verbosity=2)


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



## Directory for IDL tests is /Users/joe/gprofile_develop/

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
trace_in = trace_in[5:,:]
#right = trace_in[5:,:] + box_rad
image = sciimg - sky_model
ivar = sciivar*(slitmask == 3)
wave = objstruct['WAVE_BOX'][5]
flux = objstruct['FLUX_BOX'][5]
fluxivar = objstruct['IVAR_BOX'][5]
hwidth = objstruct['MASKWIDTH'][5]
thisfwhm = objstruct['FWHM'][5]

SN_GAUSS = None
MAX_TRACE_CORR = None
wvmnx = None

if SN_GAUSS == None: SN_GAUSS = 3.0
if thisfwhm == None: thisfwhm = 4.0
if hwidth == None: 3.0*(np.max(thisfwhm) + 1.0)
if MAX_TRACE_CORR == None:  MAX_TRACE_CORR = 2.0
if wvmnx == None: wvmnx = [2900., 30000]

thisfwhm = np.fmax(thisfwhm,1.0) # require the FWHM to be greater than 1 pixel

xnew = trace_in
dims = image.shape
nspat = dims[1]
nspec = dims[0]

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
print('sqrt(med(S/N)^2) = ', np.sqrt(med_sn2))

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


#flux  = extract_boxcar(image,trace,box_rad)





# Omit the model functionality right now
#if model != None:
#    model = image*0
#    modelwt = image*0



#ncol =



