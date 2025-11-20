
""" For the development and testing of Extraction
"""

from __future__ import (print_function, absolute_import, division, unicode_literals)



from astropy.table import Table
from pypit import ginga
from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
import sys
from scipy.io import readsav
from specobj import SpecObj
from scipy.special import ndtr
from IPython import embed


# PYPIT imports
from pypit import ginga
from pypit import msgs
from pypit.core.arextract import extract_boxcar, fit_profile, extract_optimal
from pypit.core.arskysub import skyoptimal
from pypit.core.pydl import bspline

from pypit import ardebug as debugger
debug = debugger.init()
debug['develop'] = True
msgs.reset(debug=debug, verbosity=2)




# Hack in IDL LOW_REDUX outputs. To be run on PYPIT something later
savefile='/Users/joe/gprofile_develop/local_skysub_dev.sav'
idl_dict=readsav(savefile)
sciimg = idl_dict['sciimg']
sciivar = idl_dict['sciivar']
skyimage = idl_dict['skyimage']
piximg = idl_dict['piximg']
waveimg = idl_dict['waveimg']
ximg = idl_dict['ximg']
objstruct = idl_dict['objstruct']
slitmask = idl_dict['slitmask']
slitid = idl_dict['slitid']
xx1 = idl_dict['xx1']
xx2 = idl_dict['xx2']
edgmask = np.array(idl_dict['edgmask'],dtype=bool)
bsp = idl_dict['bsp']
rn_img = idl_dict['rn_img']
outmask = idl_dict['outmaskt']
modelivar = idl_dict['modelivart']

nspat =sciimg.shape[1]
nspec =sciimg.shape[0]
slitid = idl_dict['slitid']
nobj = idl_dict['nobj']

slit_left = xx1[2,:]
slit_righ = xx2[2,:]

thismask = (slitmask == slitid)



# Create the specobj object, which is an argument to localskysub
from specobj import SpecObj

# Some information about this slit we need for later when we instantiate specobj objects
config = 'lris_b1200'
scidx = 1
det = 1
objtype= 'science'
CHK = True
spec_vec = np.arange(nspec)
spat_vec = np.arange(nspat)
slit_spec_pos = nspec/2.0
slit_spat_pos = (np.interp(slit_spec_pos, spec_vec, slit_left), np.interp(slit_spec_pos, spec_vec, slit_righ))


specobjs =[]
for ii in range(nobj):
    specobj = SpecObj(sciimg.shape, slit_spat_pos, slit_spec_pos, det=det, config=config, slitid=slitid, scidx=scidx,
                      objtype=objtype)
    # Add some attributes that I will need
    specobj.trace_spat = objstruct[ii]['xpos']
    specobj.trace_spec = objstruct[ii]['ypos']
    specobj.maskwidth = objstruct[ii]['maskwidth']
    specobj.fwhm = objstruct[ii]['fwhm']
    specobj.fwhmfit = np.zeros(nspec)
    specobj.slitid = objstruct[ii]['slitid']
    specobj.objid = objstruct[ii]['objid']
    specobj.spat_medpos = np.median(objstruct[ii]['xpos'])
    specobjs.append(specobj)


#
#TODO Remove djs_median dependency and replace with scipy.medfilt

# This is the argument list
#def localskysub(sciimg, sciivar, skyimage, rn_img, piximg, waveimg, ximg, thismask, edgmask, slit_left, slit_righ, bsp, outmask, modelivar, specobjs,
# PROF_NSIGMA = None, niter=4, box_rad = 7, sigrej = 3.5, skysample = False, FULLWELL = 5e5,MINWELL = -1000.0, SN_GAUSS = 3.0):

# outmask modelivar are created in the return and returned
# specobjs are modified "in place"

## ximg and edgmask should be created in the routine

## rn_img = Read noise image is created upstram in PYPIT from arprocimg

## slit_left and slit_right could be the trace slits object? Or maybe it easier to not have an object here

# Optional arguments
SN_GAUSS = 3.0
niter = 4
box_rad = 7
sigrej = 3.5
skysample = False
NOLOCAL = False
PROF_NSIGMA = None
FULLWELL = 5e5 # pixels above saturation level are masked
MINWELL = -1000.0 # Extremely negative pixels are also masked
SKYSAMPLE = False
COADD_2D=False
STD = False # is this a standard star?
CHK = True # Inspect the model fits with ginga?
#modelivar = None

# local skysub starts here

# Copy the specobjs that will be the output
nobj = len(specobjs)
#specobjs = copy.deepcopy(specobjs_in)

if(PROF_NSIGMA is None):
    prof_nsigma1 = np.full(len(specobjs),None)
elif len(PROF_NSIGMA) == 1:
    prof_nsigma1 = np.full(nobj, PROF_NSIGMA)
elif len(PROF_NSIGMA) == nobj:
    prof_nsigma1 = PROF_NSIGMA
else:
    raise ValueError('Invalid size for PROF_NSIGMA.')

# Set some rejection parameters based on whether this is a STD or not. Only reject extreme outliers for standards
# since super high S/N and low order profile models imply we will always have large outliers
if STD is True:
    chi2_sigrej = 100.0
    sigrej_ceil = 1e10
else:
    chi2_sigrej = 6.0
    sigrej_ceil = 10.0
# We will use this number later
gauss_prob = 1.0 - 2.0 * ndtr(-sigrej)


for iobj in range(nobj):
    specobjs[iobj].prof_nsigma = prof_nsigma1[iobj]

nspat =sciimg.shape[1]
nspec =sciimg.shape[0]

# Initialize the output mask
outmask = (sciivar > 0.0) & thismask & np.isfinite(sciimg) & (sciimg < FULLWELL) & (sciimg > MINWELL)
# Initiatlize modelivar
modelivar = sciivar
objimage = np.zeros_like(sciimg)


#sciimg_this = sciimg[thismask]
#outmask[thismask] = thismask[thismask] & (sciivar[thismask] > 0.0)



varnoobj = np.abs(skyimage - np.sqrt(2.0) * rn_img) + rn_img ** 2

xarr = np.outer(np.ones(nspec),np.arange(nspat))
yarr = np.outer(np.arange(nspec),np.ones(nspat))

xa_min = xarr[thismask].min()
xa_max = xarr[thismask].max()
ya_min = yarr[thismask].min()
ya_max = yarr[thismask].max()

xsize = slit_righ - slit_left
spatial_img = thismask*ximg*(np.outer(xsize,np.ones(nspat)))

# Loop over objects and group them
i1 = 0
while i1 < nobj:
    group = []
    group.append(i1)
    # The default value of maskwidth = 3.0 * FWHM = 7.05 * sigma in long_objfind with a log(S/N) correction for bright objects
    mincols = np.maximum(specobjs[i1].trace_spat - specobjs[i1].maskwidth - 1,slit_left)
    maxcols = np.minimum(specobjs[i1].trace_spat + specobjs[i1].maskwidth + 1,slit_righ)
    for i2 in range(i1+1,nobj):
        left_edge = specobjs[i2].trace_spat - specobjs[i2].maskwidth - 1
        righ_edge = specobjs[i2].trace_spat + specobjs[i2].maskwidth + 1
        touch = (left_edge < maxcols) & (specobjs[i2].trace_spat > slit_left) & (righ_edge > mincols)
        if touch.any():
            maxcols = np.minimum(np.maximum(righ_edge, maxcols),slit_righ)
            mincols = np.maximum(np.minimum(left_edge, mincols),slit_left)
            group.append(i2)
    # Keep for next iteration
    i1 = max(group) + 1
    # Some bookeeping to define the sub-image and make sure it does not land off the mask
    objwork = len(group)
    scope = np.sum(thismask,axis=0)
    iscp, = np.where(scope)
    imin = min(iscp)
    imax = max(iscp)
    mincol = np.fmax(np.floor(min(mincols)),imin)
    maxcol = np.fmin(np.ceil(max(maxcols)),imax)
    nc = int(maxcol - mincol + 1)
    rows = np.arange(nspec, dtype=np.intp)
    columns = np.arange(mincol, mincol + nc, dtype=np.intp)
    ipix = np.ix_(rows,columns)
    skymask = outmask & ~edgmask
    if nc > 100:
        npoly = 3
    elif nc > 40:
        npoly = 2
    else:
        npoly = 1
    obj_profiles = np.zeros((nspec,nspat,objwork), dtype=float)
    sigrej_eff = sigrej
    for iiter in range(1,niter):
        msgs.info("Iteration # " + "{:2d}".format(iiter) + " of " + "{:2d}".format(niter))
        img_minsky = sciimg - skyimage
        for ii in range(objwork):
            iobj = group[ii]
            if iiter == 1:
                # If this is the first iteration, print status message. Initiate profile fitting with a simple
                # boxcar extraction.
                msgs.info("-------------------REDUCING-------------------")
                msgs.info("Fitting profile for obj #: " + "{:d}".format(specobjs[iobj].objid) + " of {:d}".format(nobj))
                msgs.info("At x = {:5.2f}".format(specobjs[iobj].spat_medpos) + " on slit # {:d}".format(specobjs[iobj].slitid))
                msgs.info("----------------------------------------------")
                flux = extract_boxcar(img_minsky*outmask, specobjs[iobj].trace_spat, box_rad, ycen = specobjs[iobj].trace_spec)
                mvarimg = 1.0/(modelivar + (modelivar == 0))
                mvar_box = extract_boxcar(mvarimg*outmask, specobjs[iobj].trace_spat, box_rad, ycen = specobjs[iobj].trace_spec)
                pixtot = extract_boxcar(0*mvarimg + 1.0, specobjs[iobj].trace_spat, box_rad, ycen = specobjs[iobj].trace_spec)
                mask_box = (extract_boxcar(~outmask, specobjs[iobj].trace_spat, box_rad, ycen=specobjs[iobj].trace_spec) != pixtot)
                box_denom = extract_boxcar(waveimg > 0.0, specobjs[iobj].trace_spat, box_rad, ycen = specobjs[iobj].trace_spec)
                wave = extract_boxcar(waveimg, specobjs[iobj].trace_spat, box_rad, ycen = specobjs[iobj].trace_spec)/(box_denom + (box_denom == 0.0))
                fluxivar = mask_box/(mvar_box + (mvar_box == 0.0))
            else:
                # For later iterations, profile fitting is based on an optimal extraction
                last_profile = obj_profiles[:,:,ii]
                trace = np.outer(specobjs[iobj].trace_spat, np.ones(nspat))
                objmask = ((xarr >= (trace - 2.0*box_rad)) & (xarr <= (trace + 2.0*box_rad)))
                extract_optimal(waveimg,img_minsky,modelivar, (outmask & objmask), last_profile, skyimage,rn_img,box_rad, specobjs[iobj])
                # If the extraction is bad do not update
                if specobjs[iobj].optimal['MASK_OPT'].any():
                    flux = specobjs[iobj].optimal['FLUX_OPT']
                    fluxivar = specobjs[iobj].optimal['IVAR_OPT']
                    wave = specobjs[iobj].optimal['WAVE_OPT']

            if wave.any():
                (profile_model, xnew, fwhmfit, med_sn2) = fit_profile(img_minsky[ipix], (modelivar*outmask)[ipix],
                                                                      waveimg[ipix],specobjs[iobj].trace_spat - mincol,
                                                                      wave, flux, fluxivar,
                                                                      thisfwhm = specobjs[iobj].fwhm,
                                                                      hwidth = specobjs[iobj].maskwidth,
                                                                      PROF_NSIGMA = specobjs[iobj].prof_nsigma,
                                                                      SN_GAUSS =SN_GAUSS)
                # Update the object profile and the fwhm and mask parameters
                obj_profiles[ipix[0], ipix[1], ii] = profile_model
                specobjs[iobj].trace_spat = xnew + mincol
                specobjs[iobj].fwhmfit = fwhmfit
                specobjs[iobj].fwhm = np.median(fwhmfit)
                mask_fact = 1.0 + 0.5*np.log10(np.fmax(np.sqrt(np.fmax(med_sn2,0.0)),1.0))
                maskwidth = 3.0*np.median(fwhmfit)*mask_fact
                if specobjs[iobj].prof_nsigma is None:
                    specobjs[iobj].maskwidth = maskwidth
                else:
                    specobjs[iobj].maskwidth = specobjs[iobj].prof_nsigma*(specobjs[iobj].fwhm/2.3548)

            else:
                msgs.warn("Bad extracted wavelengths in local_skysub")
                msgs.warn("Skipping this profile fit and continuing.....")

        sky_bmodel = np.array(0.0)
        iterbsp = 0
        while (sky_bmodel.any()==False) & (iterbsp <=5):
            bsp_now = (1.2**iterbsp)*bsp
            # if skysample is set, determine optimal break-point spacing
            # directly measuring how well we are sampling of the sky. The
            # bsp in this case correspons to the minimum distance between
            # breakpoints which we allow.
            if SKYSAMPLE:
                sampmask = (waveimg > 0.0) & (thismask == True)
                # fullbkpt = skybkpts()
                # TODO Port long_skybkpts.pro code and put it here.
            else:
                pixvec = piximg[skymask]
                srt = pixvec.flatten().argsort()
                bset0 = bspline(pixvec.flat[srt],nord=4, bkspace=bsp_now)
                fullbkpt = bset0.breakpoints
            # check to see if only a subset of the image is used.
            # if so truncate input pixels since this can result in singular matrices
            ibool = (yarr >= ya_min) & (yarr <= ya_max) & (xarr >= xa_min) & (xarr <= xa_max) & (xarr >= mincol) & (xarr <= maxcol) & thismask
            isub, = np.where(ibool.flatten())
            sortpix = (piximg.flat[isub]).argsort()
            ithis, = np.where(thismask.flat[isub])
            keep = (fullbkpt >= piximg.flat[isub[ithis]].min()) & (fullbkpt <= piximg.flat[isub[ithis]].max())
            fullbkpt = fullbkpt[keep]
            obj_profiles_flat = obj_profiles.reshape(nspec*nspat, objwork)
            (sky_bmodel, obj_bmodel, outmask_opt) = skyoptimal(piximg.flat[isub], sciimg.flat[isub],
                                                               (modelivar*skymask).flat[isub],
                                                               obj_profiles_flat[isub, :], sortpix, spatial=spatial_img.flat[isub],
                                                               fullbkpt=fullbkpt, sigrej=sigrej_eff, npoly=npoly)
            iterbsp = iterbsp + 1
            if (sky_bmodel.any() is False) & (iterbsp <= 4):
                msgs.warn('***************************************')
                msgs.warn('WARNING: bspline sky-subtraction failed')
                msgs.warn('Increasing bkpt spacing by 20%. Retry')
                msgs.warn('Old bsp = {:5.2f}'.format(bsp_now) + '; New bsp = {:5.2f}'.format(1.2**(iterbsp)*bsp))
                msgs.warn('***************************************')

        if(sky_bmodel.any() == True):
            skyimage.flat[isub] = sky_bmodel
            objimage.flat[isub] = obj_bmodel
            img_minsky.flat[isub]=sciimg.flat[isub] - sky_bmodel
            var = np.abs(sky_bmodel + obj_bmodel - np.sqrt(2.0)*rn_img.flat[isub]) + rn_img.flat[isub]**2
            var_no = np.abs(sky_bmodel - np.sqrt(2.0)*rn_img.flat[isub]) + rn_img.flat[isub]**2
            igood1 = skymask.flat[isub]
            #  update the outmask for only those pixels that were fit. This prevents masking of slit edges in outmask
            outmask.flat[isub[igood1]]=outmask_opt[igood1]
            #  For weighted co-adds, the variance of the image is no longer equal to the image, and so the modelivar
            #  eqn. below is not valid. However, co-adds already have the model noise propagated correctly in sciivar,
            #  so no need to re-model the variance
            if COADD_2D is False:
                modelivar.flat[isub] = (var > 0.0)/(var + (var == 0.0))
                varnoobj.flat[isub]  = var_no
            # Now do some masking based on this round of model fits
            chi2 = (img_minsky.flat[isub] - obj_bmodel) ** 2 * modelivar.flat[isub]
            igood = (skymask.flat[isub]) & (chi2 <= chi2_sigrej ** 2)
            ngd = np.sum(igood)
            if ngd > 0:
                chi2_good = chi2[igood]
                chi2_srt = np.sort(chi2_good)
                sigind = np.fmin(int(np.rint(gauss_prob * float(ngd))), ngd - 1)
                chi2_sigrej = chi2_srt[sigind]
                sigrej_eff = np.fmax(np.sqrt(chi2_sigrej), sigrej)
                #  Maximum sigrej is sigrej_ceil (unless this is a standard)
                sigrej_eff = np.fmin(sigrej_eff, sigrej_ceil)
                msgs.info('Measured effective rejection from distribution of chi^2')
                msgs.info('Instead of rejecting sigrej = {:5.2f}'.format(sigrej) +
                ', use threshold sigrej_eff = {:5.2f}'.format(sigrej_eff))
                # Explicitly mask > sigrej outliers using the distribution of chi2 but only in the region that was actually fit.
                # This prevents e.g. excessive masking of slit edges
                outmask.flat[isub[igood1]] = outmask.flat[isub[igood1]] & (chi2[igood1] < chi2_sigrej) & (sciivar.flat[isub[igood1]] > 0.0)
                nrej = outmask.flat[isub[igood1]].sum()
                msgs.info('Iteration = {:d}'.format(iiter) + ', rejected {:d}'.format(nrej) + ' of ' + '{:d}'.format(igood1.sum()) + 'fit pixels')

        else:
            msgs.warn('ERROR: Bspline sky subtraction failed after 4 iterations of bkpt spacing')
            msgs.warn('       Moving on......')
            obj_profiles= np.zeros_like(obj_profiles)


    # Now that the iterations of profile fitting and sky subtraction are completed,
    # loop over the objwork objects in this grouping and perform the final extractions.
    for ii in range(objwork):
        iobj = group[ii]
        msgs.info('Extracting for obj # {:d}'.format(iobj+1) + ' of {:d}'.format(nobj) +
                  ' on slit # {:d}'.format(specobjs[iobj].slitid) + ' at x = {:5.2f}'.format(np.median(specobjs[iobj].trace_spat)))
        this_profile = obj_profiles[:,:,ii]
        trace = np.outer(specobjs[iobj].trace_spat, np.ones(nspat))
        objmask = ((xarr >= (trace - 2.0 * box_rad)) & (xarr <= (trace + 2.0 * box_rad)))
        extract_optimal(waveimg, img_minsky, modelivar*thismask, (outmask & objmask), this_profile, skyimage, rn_img, box_rad, specobjs[iobj])
        specobjs[iobj].mincol = mincol
        specobjs[iobj].maxcol = maxcol

    # If requested display the model fits for this grouping
    if CHK == True:
        viewer, ch = ginga.show_image((sciimg - skyimage)*np.sqrt(modelivar))
        # TODO figure out a way to overplot the pixels that were masked in red like as a scatter plot
        for ii in range(objwork):
            iobj = group[ii]
            ginga.show_trace(viewer, ch, specobjs[iobj].trace_spat, specobjs[iobj].idx, color='green')

#    return None, need to think about what we return versus what is changed in place

'''    #Alternative QA with hack to show traces
    if CHK == True:
        qaimg =sciimg - skyimage)*np.sqrt(modelivar)
        samp = (np.rint(np.arange(nspec/4)*4)).astype(int)
        for ii in range(objwork):
            iobj = group[ii]
            ispec = (specobjs[iobj].trace_spec[samp]).astype(int)
            ispat = (specobjs[iobj].trace_spat[samp]).astype(int)
            qaimg[ispec, ispat] = -10000
        ginga.show_image(qaimg)
'''
# return something?



'''
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
objind = 6
trace_in = trace_in[objind,:]
#right = trace_in[5:,:] + box_rad
image = sciimg - sky_model
ivar = sciivar*(slitmask == 3)
wave = objstruct['WAVE_BOX'][objind]
flux = objstruct['FLUX_BOX'][objind]
fluxivar = objstruct['IVAR_BOX'][objind]
hwidth = objstruct['MASKWIDTH'][objind]
thisfwhm = objstruct['FWHM'][objind]

# This is the argument list
SN_GAUSS = None # S/N threshold for giving up on profile fitting and using a Gaussian with the measured FWHM
MAX_TRACE_CORR = None # Maximum trace correction that can be applied
wvmnx = None # minimum and maximum wavelength to use for fits
GAUSS = None # Do a Gaussian extraction
PROF_NSIGMA = None # Width of profile specified by hand for extended objects
NO_DERIV = False # disable profile derivative computation



# Check for infinities in trace correction




#


#(mode_shift_set, mode_shift_mask, mode_shift_fit, red_chi) = mode_shift_out


# Omit the model functionality right now
#if model != None:
#    model = image*0
#    modelwt = image*0



#ncol =


# flux  = extract_boxcar(image,trace,box_rad)

# Debugging bspline_longslit
#from scipy.io import readsav
#savefile='/Users/joe/gprofile_develop/bspline_out.sav'
#idl_dict=readsav(savefile)
#xdata = idl_dict['xdata']
#ydata = idl_dict['ydata']
#prof_basis = idl_dict['profile_basis'].T
#invvar = idl_dict['invvar']
#fullbkpt = idl_dict['fullbkpt']
#mode_shift_out = bspline_longslit(xdata,ydata,invvar,prof_basis,maxiter=1,fullbkpt=fullbkpt)
#sys.exit(-1)
'''