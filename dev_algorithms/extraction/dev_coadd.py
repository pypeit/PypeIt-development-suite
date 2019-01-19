


import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
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
import glob


# NIRES
redux_path = '/Users/joe/Dropbox/PypeIt_Redux/NIRES/J0252/ut181001/'
objprefix = 'J0252-0503'
spec2d_files = glob.glob(redux_path + 'Science/spec2d_' + objprefix + '*')
slitid = 3
objid=np.full(len(spec2d_files),1)
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
specobjs_list, tslits_dict_orig, slitmask_stack, sciimg_stack, sciivar_stack, skymodel_stack, mask_stack, tilts_stack, waveimg_stack = \
    coadd2d.load_coadd2d_stacks(spec2d_files)



nexp, nspec, nspat  = sciimg_stack.shape
# Create the thismask_stack for this slit.
thismask_stack = slitmask_stack == slitid
thismask = thismask_stack[0,:,:]
slitmask = slitmask_stack[0,:,:]


# Grab the traces for this slit and objid.
trace_stack = np.zeros((nexp, nspec),dtype=float)
for iexp, sobjs in enumerate(specobjs_list):
    ithis = (sobjs.slitid == slitid) & (sobjs.objid == objid[iexp])
    trace_stack[iexp,:] = sobjs[ithis].trace_spat


spectrograph = util.load_spectrograph(tslits_dict_orig['spectrograph'])

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
wave_grid = spectrograph.wavegrid()


wave_bins, dspat_bins, sciimg_rect, sciivar_rect, imgminsky_rect, outmask_rect, nused, tilts_rect, waveimg_rect, \
    dspat_rect, thismask_rect, tslits_dict = \
    coadd2d.coadd2d(trace_stack, sciimg_stack, sciivar_stack, skymodel_stack, (mask_stack == 0), tilts_stack, waveimg_stack,
    thismask_stack, wave_grid=wave_grid)

sobjs, _ = extract.objfind(imgminsky_rect, thismask_rect, tslits_dict['lcen'], tslits_dict['rcen'],
                           inmask=outmask_rect, show_peaks=True,show_fits=True, show_trace=True)
sobjs_neg, _ = extract.objfind(-imgminsky_rect, thismask_rect, tslits_dict['lcen'], tslits_dict['rcen'],
                               inmask=outmask_rect, show_peaks=True,show_fits=True, show_trace=True)
sobjs.append_neg(sobjs_neg)
# Local sky subtraction and extraction
skymodel_rect = np.zeros_like(imgminsky_rect)
global_sky_rect = np.zeros_like(imgminsky_rect)
rn2img_rect = np.zeros_like(imgminsky_rect)
objmodel_rect = np.zeros_like(imgminsky_rect)
ivarmodel_rect = np.zeros_like(imgminsky_rect)
extractmask_rect = np.zeros_like(thismask_rect)
par = spectrograph.default_pypeit_par()
# TODO Modify profile fitting so that we can pass in a position image which will allow us to improve the spatial sampling
# in this final extraction step
skymodel_rect[thismask_rect], objmodel_rect[thismask_rect], ivarmodel_rect[thismask_rect], extractmask_rect[thismask_rect] = \
    skysub.local_skysub_extract(imgminsky_rect, sciivar_rect, tilts_rect, waveimg_rect, global_sky_rect, rn2img_rect, thismask_rect,
    tslits_dict['lcen'], tslits_dict['rcen'], sobjs, model_noise=False,bsp=par['scienceimage']['bspline_spacing'],
    sn_gauss=par['scienceimage']['sn_gauss'],inmask=outmask_rect, show_profile=True)
resids_rect = ((imgminsky_rect - skymodel_rect - objmodel_rect)*np.sqrt(ivarmodel_rect))
resids_pix = resids_rect[extractmask_rect]
std_dev = np.std(resids_pix)
mean_resid = np.mean(resids_pix)
sobjs_out = sobjs.copy()
sobjs_out.purge_neg()


viewer, ch = ginga.show_image(imgminsky_rect,'imgminsky')
viewer, ch = ginga.show_image((imgminsky_rect - skymodel_rect)*np.sqrt(ivarmodel_rect),'sky_resids')
viewer, ch = ginga.show_image((imgminsky_rect - skymodel_rect - objmodel_rect)*np.sqrt(ivarmodel_rect),'resids')

# After displaying all the images sync up the images with WCS_MATCH
shell = viewer.shell()
out = shell.start_global_plugin('WCSMatch')
out = shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', ['imgminsky'], {})


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
