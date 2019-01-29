


import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from astropy import time
from pypeit import msgs
from pypeit.core import arc
from pypeit import utils
from pypeit.core import parse
from pypeit.core import pca
from pypeit.core import qa
from pypeit import traceslits
from pypeit.core import trace_slits
from pypeit.core import extract
from pypeit.spectrographs import util
from astropy.table import Table
from astropy import stats
from astropy import constants as const
from astropy import units as u
from pypeit.core import extract
from pypeit.core import skysub
from pypeit.core import procimg
from pypeit import ginga
from pypeit import masterframe
from pypeit.core import load
from pypeit.core import coadd2d
from pypeit.core import save
from pypeit.core import pixels
from pypeit import reduce
from collections import OrderedDict
from linetools import utils as ltu
import scipy

from pypeit import processimages
import glob


# NIRES
redux_path = '/Users/joe/Dropbox/PypeIt_Redux/NIRES/J0252/ut181001_new/'
objprefix = 'J0252-0503'
spec2d_files = glob.glob(redux_path + 'Science/spec2d_' + objprefix + '*')
#objid=np.full(len(spec2d_files),1)

#objid=1
# GNIRS
#redux_path = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Gemini_GNIRS/'
#objprefix ='pisco'
#slitid = 2
#spec2d_files = glob.glob(redux_path + 'Science/spec2d_' + objprefix + '*')
#objid=np.full(len(spec2d_files),1)

#objid = 1


# XSHOOTER
#redux_path = '/Users/joe/Dropbox/PypeIt_Redux/XSHOOTER/J0020m3653/NIR/'
#objprefix ='VHSJ0020'
#slitid=11
#objid = [1,1,1,3]
# read in the stacks


def extract_one_coadd2d(spec2d_files, ir_redux=False, par=None, show=False):

    # Read in the stack, grab some meta data we will need
    stack_dict = coadd2d.load_coadd2d_stacks(spec2d_files)
    head1d = stack_dict['head1d_list'][0]
    head2d = stack_dict['head2d_list'][0]
    try:
        mjd = head1d['mjd']  # recorded as 'mjd' in fitstbl
    except KeyError:
        mjd = head1d['MJD-OBS']
    obstime = time.Time(mjd, format='mjd')
    filename = os.path.basename(spec2d_files[0])
    basename =

    # Find the objid of the brighest object, and the average snr across all orders
    nslits = stack_dict['tslits_dict']['slit_left'].shape[1]
    objid, snr_bar = coadd2d.get_brightest_obj(stack_dict['specobjs_list'], echelle=True)
    # TODO Print out a report here on the image stack

    spectrograph = util.load_spectrograph(stack_dict['spectrograph'])
    par = spectrograph.default_pypeit_par() if par is None else par

    # Grab the wavelength grid that we will rectify onto
    wave_grid = spectrograph.wavegrid()

    ## Loop over the slit and create these stacked images for each order.
    ## the rest of our routines will run on them, like ech_objfind etc. Generalize the Scienceimage class to handle
    # this and/or strip those methods out to functions.

    coadd_list = []
    nspec_vec = np.zeros(nslits,dtype=int)
    nspat_vec = np.zeros(nslits,dtype=int)
    # ToDO Generalize this to be a loop over detectors, such tha the coadd_list is an ordered dict (perhaps) with
    # all the slits on all detectors
    for islit in range(nslits):
        msgs.info('Performing 2d coadd for slit: {:d}/{:d}'.format(islit,nslits-1))
        # Determine the wavelength dependent optimal weights and grab the reference trace
        rms_sn, weights, trace_stack, wave_stack = coadd2d.optimal_weights(stack_dict['specobjs_list'], islit, objid)
        thismask_stack = stack_dict['slitmask_stack'] == islit
        # Perform the 2d coadd
        coadd_dict = coadd2d.coadd2d(trace_stack, stack_dict['sciimg_stack'], stack_dict['sciivar_stack'],
                                     stack_dict['skymodel_stack'], (stack_dict['mask_stack'] == 0),
                                     stack_dict['tilts_stack'], stack_dict['waveimg_stack'], thismask_stack,
                                     weights = weights, wave_grid=wave_grid)
        coadd_list.append(coadd_dict)
        nspec_vec[islit]=coadd_dict['nspec']
        nspat_vec[islit]=coadd_dict['nspat']

    # Determine the size of the psuedo image
    nspat_pad = 10
    nspec_psuedo = nspec_vec.max()
    nspat_psuedo = np.sum(nspat_vec) + (nslits + 1)*nspat_pad
    spec_vec_psuedo = np.arange(nspec_psuedo)
    shape_psuedo = (nspec_psuedo, nspat_psuedo)
    imgminsky_psuedo = np.zeros(shape_psuedo)
    sciivar_psuedo = np.zeros(shape_psuedo)
    waveimg_psuedo = np.zeros(shape_psuedo)
    tilts_psuedo = np.zeros(shape_psuedo)
    spat_psuedo = np.zeros(shape_psuedo)
    nused_psuedo = np.zeros(shape_psuedo, dtype=int)
    inmask_psuedo = np.zeros(shape_psuedo, dtype=bool)
    wave_mid = np.zeros((nspec_psuedo, nslits))
    dspat_mid = np.zeros((nspat_psuedo, nslits))

    spat_left = nspat_pad
    slit_left = np.zeros((nspec_psuedo, nslits))
    slit_righ = np.zeros((nspec_psuedo, nslits))
    spec_min1 = np.zeros(nslits)
    spec_max1 = np.zeros(nslits)
    for islit, coadd_dict in enumerate(coadd_list):
        spat_righ = spat_left + nspat_vec[islit]
        ispec = slice(0,nspec_vec[islit])
        ispat = slice(spat_left,spat_righ)
        imgminsky_psuedo[ispec, ispat] = coadd_dict['imgminsky']
        sciivar_psuedo[ispec, ispat] = coadd_dict['sciivar']
        waveimg_psuedo[ispec, ispat] = coadd_dict['waveimg']
        tilts_psuedo[ispec, ispat] = coadd_dict['tilts']
        # spat_psuedo is the sub-pixel image position on the rebinned psuedo image
        inmask_psuedo[ispec, ispat] = coadd_dict['outmask']
        image_temp = (coadd_dict['dspat'] -  coadd_dict['dspat_mid'][0] + spat_left)*coadd_dict['outmask']
        spat_psuedo[ispec, ispat] = image_temp
        nused_psuedo[ispec, ispat] = coadd_dict['nused']
        wave_mid[ispec, islit] = coadd_dict['wave_mid']
        dspat_mid[ispat, islit] = coadd_dict['dspat_mid']
        slit_left[:,islit] = np.full(nspec_psuedo, spat_left)
        slit_righ[:,islit] = np.full(nspec_psuedo, spat_righ)
        spec_max1[islit] = nspec_vec[islit]-1
        spat_left = spat_righ + nspat_pad

    tslits_dict_psuedo = dict(slit_left=slit_left, slit_righ=slit_righ, nspec=nspec_psuedo, nspat=nspat_psuedo, pad=0,
                              nslits = nslits, binspectral=1, binspatial=1, spectrograph=spectrograph.spectrograph,
                              spec_min=spec_min1, spec_max=spec_max1)

    slitmask_psuedo = pixels.tslits2mask(tslits_dict_psuedo)
    # This is a kludge to deal with cases where bad wavelengths result in large regions where the slit is poorly sampled,
    # which wreaks havoc on the local sky-subtraction
    min_slit_frac = 0.70
    spec_min = np.zeros(nslits)
    spec_max = np.zeros(nslits)
    for islit in range(nslits):
        slit_width = np.sum(inmask_psuedo*(slitmask_psuedo == islit),axis=1)
        slit_width_img = np.outer(slit_width, np.ones(nspat_psuedo))
        med_slit_width = np.median(slit_width_img[slitmask_psuedo==islit])
        nspec_eff = np.sum(slit_width > min_slit_frac*med_slit_width)
        nsmooth = int(np.fmax(np.ceil(nspec_eff*0.02),10))
        slit_width_sm = scipy.ndimage.filters.median_filter(slit_width, size=nsmooth, mode='reflect')
        igood = (slit_width_sm > min_slit_frac*med_slit_width)
        spec_min[islit] = spec_vec_psuedo[igood].min()
        spec_max[islit] = spec_vec_psuedo[igood].max()
        bad_pix = (slit_width_img < min_slit_frac*med_slit_width) & (slitmask_psuedo == islit)
        inmask_psuedo[bad_pix] = False

    # Update with tslits_dict_psuedo
    tslits_dict_psuedo['spec_min'] = spec_min
    tslits_dict_psuedo['spec_max'] = spec_max
    slitmask_psuedo = pixels.tslits2mask(tslits_dict_psuedo)

    # Make a fake bitmask from the outmask. We are kludging the crmask to be the outmask_psuedo here, and setting the bpm to
    # be good everywhere
    mask = processimages.ProcessImages.build_mask(imgminsky_psuedo, sciivar_psuedo, np.invert(inmask_psuedo),
                                                  np.zeros_like(inmask_psuedo), slitmask=slitmask_psuedo)

    redux = reduce.instantiate_me(spectrograph, tslits_dict_psuedo, mask, ir_redux=ir_redux,
                                  par=par['scienceimage'], frame_par=par['scienceframe'],
                                  objtype = 'science')

    # Object finding
    sobjs_obj, nobj, skymask_init = redux.find_objects(imgminsky_psuedo, sciivar_psuedo, ir_redux=ir_redux,
                                                       show_peaks=False, show=show)
    # Local sky-subtraction
    global_sky_psuedo = np.zeros_like(imgminsky_psuedo) # No global sky for co-adds since we go straight to local
    rn2img_psuedo = global_sky_psuedo # No rn2img for co-adds since we go do not model noise
    skymodel_psuedo, objmodel_psuedo, ivarmodel_psuedo, outmask_psuedo, sobjs = \
        redux.local_skysub_extract(imgminsky_psuedo, sciivar_psuedo, tilts_psuedo, waveimg_psuedo, global_sky_psuedo,
                                   rn2img_psuedo, sobjs_obj, spat_pix=spat_psuedo,
                                   model_noise=False, show_profile=False, show=show)

    from IPython import embed
    embed()

    if ir_redux:
        sobjs.purge_neg()

    # Flexure correction
    redux.flexure_correct(sobjs, basename)

    # Grab coord
    radec = ltu.radec_to_coord(head1d['ra'], head1d['dec'])
    vel_corr = redux.helio_correct(sobjs, radec, self.obstime)

    # TODO Add flexure and heliocentric correction here
    vel_corr=None

    sci_dict = {}
    sci_dict['meta'] = {}
    sci_dict['meta']['vel_corr'] = 0.

    sci_dict['sciimg'] = imgminsky_psuedo
    sci_dict['sciivar'] = sciivar_psuedo
    sci_dict['skymodel']= skymodel_psuedo
    sci_dict['objmodel']=objmodel_psuedo
    sci_dict['ivarmodel']=ivarmodel_psuedo
    sci_dict['outmask'] = outmask_psuedo
    sci_dict['specobjs'] = sobjs
    if vel_corr is not None:
        sci_dict['meta']['vel_corr'] = vel_corr

    # Create a master key dict
    master_key_dict={}
    # Create a new master dir?
    master_dir=''

    # headers

    save.save_all(sci_dict, master_key_dict, master_dir, spectrograph, head1d, head2d, scipath, basename)


extract_one_coadd2d(spec2d_files, ir_redux=True, par=None, show=True)

igd = sobjs_out[0].optimal['WAVE'] > 1.0
wave = sobjs_out[0].optimal['WAVE'][igd]
flux = sobjs_out[0].optimal['COUNTS'][igd]
sig = sobjs_out[0].optimal['COUNTS_SIG'][igd]

plt.plot(wave, flux,drawstyle='steps-mid')
plt.plot(wave, sig,drawstyle='steps-mid')
plt.show()
sys.exit(-1)
# This is the wavelength grid of the rectified images. It differs from the 'true' wavelengths by 2% of a pixel
loglam_bins = np.log10(wave_bins)
loglam_mid = ((loglam_bins + np.roll(loglam_bins,1))/2.0)[1:]
#wave_mid = np.power(10.0,loglam_mid)
plt.plot(loglam_mid[igd],(loglam_mid[igd] - np.log10(wave))/loglam_mid[igd]/dloglam)


dspat_mid = ((dspat_bins + np.roll(dspat_bins,1))/2.0)[1:]
loglam_img = np.outer(loglam_mid,np.ones(nspat_rect))
# Denom is computed in case the trace goes off the edge of the image
loglam_extract_exp = np.log10(waveimg_rect_stack[0, :, :] + (waveimg_rect_stack[0,:,:] == 0.0))
# We centered everything with respect to the 0 position, so the trace runs vertically
zero_dspat = np.interp(0.0, dspat_mid, np.arange(dspat_mid.shape[0]))
trace = np.full(loglam_extract_exp.shape[0],zero_dspat)


# Which image do we use for creating the mask?
#data = np.ma.MaskedArray(imgminsky_rect_stack, (norm_rect_stack == 0))
#sigclip = stats.SigmaClip(sigma=sigrej, maxiters=maxiters, cenfunc='median')
#data_clipped = sigclip(data, axis=0, masked=True)
#mask_rect_stack = np.invert(data_clipped.mask)  # mask_rect_stack = True are good values
#nused1 = np.sum(mask_rect_stack,axis=0)
#weights_stack = np.einsum('i,ijk->ijk', weights, mask_rect_stack)
#weights_sum = np.sum(weights_stack, axis=0)
#sciimg1 = np.sum(sciimg_rect_stack*weights_stack,axis=0)/(weights_sum + (weights_sum == 0.0))
#waveimg1 = np.sum(waveimg_rect_stack*weights_stack,axis=0)/(weights_sum + (weights_sum == 0.0))
#tilts1 = np.sum(tilts_rect_stack*weights_stack,axis=0)/(weights_sum + (weights_sum == 0.0))
#imgminsky1 = np.sum(imgminsky_rect_stack*weights_stack,axis=0)/(weights_sum + (weights_sum == 0.0))
#varfinal1 = np.sum(var_rect_stack * weights_stack ** 2, axis=0) / (weights_sum + (weights_sum == 0.0)) ** 2
#sciivar1 = utils.calc_ivar(varfinal1)
#outmask1 = np.sum(mask_rect_stack,axis=0) > 0


#pad_spat = 5
#spec_ind, spat_ind = np.where(thismask)
#min_spat = spat_ind.min() - pad_spat
#max_spat = spat_ind.max() + pad_spat
#num_spat = max_spat - min_spat + 1

# Need to determine the limits of the spat img
#spat_lin = np.linspace(min_spat,max_spat,num = int(np.round(num_spat)))
#spat_img = np.outer(np.ones(nwave), spat_lin)
#spat_img, loglam_img = np.meshgrid(spat_lin, loglam_slit)
#    slit_cen_lin = (interpolate.interp1d(wave_stack[iexp,:],trace_stack[iexp,:],bounds_error=False ,fill_value='extrapolate'))(wave_slit)



# Single exposure boxcar extraction
#box_radius=7.0
#box_denom_exp = extract.extract_boxcar(loglam_extract_exp*norm_rect_stack[0,:,:] > 0.0, trace, box_radius)
#loglam_box_exp = extract.extract_boxcar(loglam_extract_exp*(norm_rect_stack[0,:,:] > 0.0),trace, box_radius)/(box_denom_exp + (box_denom_exp == 0.0))
#delta_exp =((loglam_box_exp-loglam_mid)/dloglam)[box_denom_exp > 0]
#plt.hist(delta_exp,bins=40)
#plt.ylabel('Number of Pixels')
#plt.xlabel('Offset from fixed grid (pixels)')

# Coadd boxcar extraction
#loglam_extract = np.log10(waveimg + (waveimg == 0.0))
#box_radius=7.0
#box_denom = extract.extract_boxcar(loglam_extract*outmask > 0.0, trace, box_radius)
#loglam_box = extract.extract_boxcar(loglam_extract*outmask,trace, box_radius)/(box_denom + (box_denom == 0.0))
#delta =((loglam_box-loglam_mid)/dloglam)[box_denom > 0]
#plt.hist(delta,bins=40)
#plt.ylabel('Number of Pixels')
#plt.xlabel('Offset from fixed grid (pixels)')

# Now do a co-add

# Extract this as a waveimg
#delta_loglam = (np.log10(waveimg_rect_stack[0, :, :]) - loglam_img)/dloglam


#ginga.show_image((np.log10(waveimg_rect_stack[0, :, :]) - loglam_img)/dloglam)


#shape = (nexp, nspec_rect, nspat_rect)
#sciimg_rect_stack = np.zeros(shape)
#imgminsky_rect_stack = np.zeros(shape)
#tilts_rect_stack = np.zeros(shape)
#waveimg_rect_stack = np.zeros(shape)
#dspat_rect_stack = np.zeros(shape)
#norm_rect_stack = np.zeros(shape)
#nsamp_rect_stack = np.zeros(shape,dtype=int)
#var_rect_stack = np.zeros(shape)
#
#
#
# for iexp in range(nexp):
#     # This fist image is purely for bookeeping purposes to determine the number of times each pixel
#     # could have been sampled
#     thismask = thismask_stack[iexp,:,:]
#     spec_rebin_this = waveimg_stack[iexp,:,:][thismask]
#     spat_rebin_this = dspat_stack[iexp,:,:][thismask]
#
#     norm_img_this, spec_edges, spat_edges = np.histogram2d(spec_rebin_this, spat_rebin_this,
#                                                   bins=[wave_bins, dspat_bins], density=False)
#     nsamp_rect_stack[iexp,:,:] = norm_img_this
#
#     finmask = thismask & (mask_stack[iexp,:,:] == 0)
#     spec_rebin = waveimg_stack[iexp,:,:][finmask]
#     spat_rebin = dspat_stack[iexp,:,:][finmask]
#     norm_img, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
#                                                   bins=[wave_bins, dspat_bins], density=False)
#     norm_rect_stack[iexp,:,:] = norm_img
#     # Rebin the science image (we are not actually using this at the moment)
#     img_now = sciimg_stack[iexp,:,:]
#     weigh_img, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
#                                                    bins=[wave_bins, dspat_bins], density=False,
#                                                    weights = img_now[finmask])
#     sciimg_rect_stack[iexp,:,:] =(norm_img > 0.0)*weigh_img/(norm_img + (norm_img == 0.0))
#     # Rebin the sky subtracted image
#     imgminsky_now = (sciimg_stack[iexp,:,:] - skymodel_stack[iexp,:,:])
#     weigh_imgminsky, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
#                                                    bins=[wave_bins, dspat_bins], density=False,
#                                                    weights = imgminsky_now[finmask])
#     imgminsky_rect_stack[iexp,:,:] =(norm_img > 0.0)*weigh_imgminsky/(norm_img + (norm_img == 0.0))
#
#     # Rebin the tilts
#     tilts_now = tilts_stack[iexp,:,:]
#     weigh_tilts, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
#                                                    bins=[wave_bins, dspat_bins], density=False,
#                                                    weights = tilts_now[finmask])
#     tilts_rect_stack[iexp,:,:] =(norm_img > 0.0)*weigh_tilts/(norm_img + (norm_img == 0.0))
#     # Rebin the wavelength image
#     waveimg_now = waveimg_stack[iexp,:,:]
#     weigh_waveimg, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
#                                                    bins=[wave_bins, dspat_bins], density=False,
#                                                    weights = waveimg_now[finmask])
#     waveimg_rect_stack[iexp,:,:] =(norm_img > 0.0)*weigh_waveimg/(norm_img + (norm_img == 0.0))
#
#     # Rebin the sptial offset from the trace in pixels
#     dspat_now = dspat_stack[iexp,:,:]
#     weigh_dspat, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
#                                                            bins=[wave_bins, dspat_bins], density=False,
#                                                            weights=dspat_now[finmask])
#     dspat_rect_stack[iexp, :, :] = (norm_img > 0.0) * weigh_dspat/(norm_img + (norm_img == 0.0))
#
#     # Rebin the variance
#     var_now = utils.calc_ivar(sciivar_stack[iexp,:,:])
#     weigh_var, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
#                                                    bins=[wave_bins, dspat_bins], density=False,
#                                                    weights = var_now[finmask])
#     var_rect_stack[iexp,:,:] =(norm_img > 0.0)*weigh_var/(norm_img + (norm_img == 0.0))**2
#

#dloglam_order = np.zeros(spectrograph.norders)
#sobjs_now = specobjs_list[0]
#for iorder in range(spectrograph.norders):
#    ithis = (sobjs_now.objid == objid[0]) & (sobjs_now.ech_orderindx == iorder)
#    wave_order = sobjs_now[ithis][0].optimal['WAVE'].value
#    igd = wave_order > 1.0
#    wave_order = wave_order[igd]
#    loglam = np.log10(wave_order)
#    dloglam = (loglam - np.roll(loglam,1))[1:]
#    dloglam_order[iorder] = np.median(dloglam)

#dloglam_mean = np.mean(dloglam_order)
#delta_loglam = (dloglam_order - dloglam_mean)/dloglam_mean
#wave_order = np.zeros(norder)
# median wavelength separation

