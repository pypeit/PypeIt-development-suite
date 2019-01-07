


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
from pypeit.spectrographs.util import load_spectrograph
from astropy.table import Table
from astropy import stats
from pypeit.core import extract
from pypeit.core import skysub
from pypeit.core import procimg
from pypeit import ginga

def get_wave_ind(wave_grid, wave_min, wave_max):
    diff = wave_grid - wave_min
    diff[diff > 0] = np.inf
    ind_lower = np.argmin(np.abs(diff))
    diff = wave_max - wave_grid
    diff[diff > 0] = np.inf
    ind_upper = np.argmin(np.abs(diff))
    return ind_lower, ind_upper



# Read in the masters
path = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT_old/Gemini_GNIRS/GNIRS/'
master_dir = path + 'MF_gemini_gnirs/'
scidir = path + 'Science/'
waveimgfiles = [os.path.join(master_dir, ifile) for ifile in ['MasterWave_A_1_01.fits','MasterWave_A_2_01.fits', 'MasterWave_A_4_01.fits',
                                                              'MasterWave_A_8_01.fits']]
tiltfiles = [os.path.join(master_dir, ifile) for ifile in ['MasterTilts_A_1_01.fits','MasterTilts_A_2_01.fits', 'MasterTilts_A_4_01.fits',
                                                           'MasterTilts_A_8_01.fits']]
tracefile = os.path.join(master_dir,'MasterTrace_A_15_01')
spec2d_files = [os.path.join(scidir, ifile) for ifile in ['spec2d_pisco_GNIRS_2017Mar31T085412.181.fits',
                                                              'spec2d_pisco_GNIRS_2017Mar31T085933.097.fits',
                                                              'spec2d_pisco_GNIRS_2017Mar31T091538.731.fits',
                                                              'spec2d_pisco_GNIRS_2017Mar31T092059.597.fits']]
spec1d_files = [os.path.join(scidir, ifile) for ifile in ['spec1d_pisco_GNIRS_2017Mar31T085412.181.fits',
                                                              'spec1d_pisco_GNIRS_2017Mar31T085933.097.fits',
                                                              'spec1d_pisco_GNIRS_2017Mar31T091538.731.fits',
                                                              'spec1d_pisco_GNIRS_2017Mar31T092059.597.fits']]

slit = 4

nexp = len(waveimgfiles)

for ifile in range(nexp):
    hdu_wave = fits.open(waveimgfiles[ifile])
    waveimg = hdu_wave[0].data
    hdu_tilts = fits.open(tiltfiles[ifile])
    tilts = hdu_tilts[0].data
    hdu_sci = fits.open(spec2d_files[ifile])
    sciimg = hdu_sci[1].data
    sky = hdu_sci[3].data
    sciivar = hdu_sci[5].data
    mask = hdu_sci[6].data
    hdu_sci1d = fits.open(spec1d_files[ifile])
    obj_table = Table(hdu_sci1d[slit+1].data)
    trace = obj_table['TRACE']
    wave_trace = obj_table['BOX_WAVE']
    if ifile == 0:
        waveimg_stack = np.zeros((nexp,waveimg.shape[0],waveimg.shape[1]),dtype=float)
        tilts_stack = np.zeros((nexp,waveimg.shape[0],waveimg.shape[1]),dtype=float)
        sciimg_stack = np.zeros((nexp,sciimg.shape[0],sciimg.shape[1]),dtype=float)
        skymodel_stack = np.zeros((nexp,sciimg.shape[0],sciimg.shape[1]),dtype=float)
        sciivar_stack = np.zeros((nexp,sciimg.shape[0],sciimg.shape[1]),dtype=float)
        mask_stack = np.zeros((nexp,sciimg.shape[0],sciimg.shape[1]),dtype=float)
        trace_stack = np.zeros((nexp,trace.shape[0]))
        wave_stack = np.zeros((nexp,trace.shape[0]))

    waveimg_stack[ifile,:,:] = waveimg
    tilts_stack[ifile,:,:] = tilts
    sciimg_stack[ifile,:,:] = sciimg
    sciivar_stack[ifile,:,:] = sciivar
    mask_stack[ifile,:,:] = mask
    skymodel_stack[ifile,:,:] = sky
    trace_stack[ifile,:] = trace
    wave_stack[ifile,:] = wave_trace

# Read in the tslits_dict
tslits_dict = traceslits.load_tslits_dict(tracefile)
#Tslits = traceslits.TraceSlits.from_master_files(tracefile)
#tslits_dict = Tslits._fill_tslits_dict()
spectrograph = load_spectrograph('gemini_gnirs')
slitmask = spectrograph.slitmask(tslits_dict)
nslits = tslits_dict['lcen'].shape[1]
slitmask_stack = np.einsum('i,jk->ijk', np.ones(nexp), slitmask)
thismask_stack = slitmask_stack == slit

# Define the new wavelength grid
ngrid = 5000
dloglam = 0.000127888 # this is the average of the median dispersions
logmin = 3.777
#loglam = logmin + dloglam*np.arange(ngrid)
# Define the oversampled grid
osamp = 1.0
loglam_grid = logmin + (dloglam/osamp)*np.arange(int(np.ceil(osamp*ngrid)))
wave_grid = None
# def coadd_2d(sciimg_stack, sciivar_stack, inmask_stack, skymodel_stack, tilts_stack, waveimg_stack, trace_stack, thismask_stack,
# loglam_grid = None, wave_grid = None


nspec, nspat = sciimg_stack[0,:,:].shape

# Determine the wavelength grid that we will use for the current slit/order
wave_min = waveimg_stack[thismask_stack].min()
wave_max = waveimg_stack[thismask_stack].max()
if loglam_grid is not None:
    ind_lower, ind_upper = get_wave_ind(loglam_grid, np.log10(wave_min), np.log10(wave_max))
    loglam_bins = loglam_grid[ind_lower:ind_upper+1]
    wave_bins = np.power(10.0, loglam_bins)
elif wave_grid is not None:
    ind_lower, ind_upper = get_wave_ind(wave_grid, wave_min, wave_max)
    wave_bins = wave_grid[ind_lower:ind_upper+1]
    loglam_bins = np.log10(wave_bins)
else:
    msgs.error('You must input either a uniformly space loglam grid or wave grid')

nbins = wave_bins.shape[0]

# Create the slit_cen_stack and determine the minimum and maximum spatial offsets that we need to cover to determine
# the spatial bins
spat_img = np.outer(np.ones(nspec), np.arange(nspat))
dspat_stack = np.zeros_like(sciimg_stack)
spat_min = np.inf
spat_max = -np.inf
for iexp in range(nexp):
    # center of the slit replicated spatially
    slit_cen_img = np.outer(trace_stack[iexp,:], np.ones(nspat))
    dspat_iexp = (spat_img - slit_cen_img)
    dspat_stack[iexp,:,:] = dspat_iexp
    thismask = thismask_stack[iexp,:,:]
    spat_min = np.fmin(spat_min,dspat_iexp[thismask].min())
    spat_max = np.fmax(spat_max,dspat_iexp[thismask].max())

spat_min_int = int(np.floor(spat_min))
spat_max_int = int(np.ceil(spat_max))
dspat_bins = np.arange(spat_min_int,spat_max_int +1,1)

# Allocate the rect_stack images
nspec_rect = nbins - 1
nspat_rect = dspat_bins.shape[0]-1

sci_list = [sciimg_stack, sciimg_stack - skymodel_stack, tilts_stack, waveimg_stack, dspat_stack]
var_list = [utils.calc_ivar(sciivar_stack)]

sci_list_rebin, var_list_rebin, norm_rebin_stack, nsmp_rebin_stack = procimg.rebin2d(
    wave_bins, dspat_bins, waveimg_stack, dspat_stack, thismask_stack, (mask_stack == 0), sci_list, var_list)

# Now compute the final stack with sigma clipping
sigrej = 3.0
maxiters = 10
# Replace this later with S/N weights
weights = np.ones(nexp)/float(nexp)
sci_list_out, var_list_out, outmask, nused = procimg.weighted_combine(
    weights, sci_list_rebin, var_list_rebin, (norm_rebin_stack != 0),
    sigma_clip=True, sigma_clip_stack=sci_list_rebin[1], sigrej=sigrej, maxiters=maxiters)
sciimg, imgminsky, tilts, waveimg, dspat = sci_list_out
sciivar = utils.calc_ivar(var_list_out[0])

sys.exit(-1)

# Now let's try to do an extraction
nsamp = np.sum(nsamp_rect_stack,axis=0) # This image indicates the number of times each pixel could have been sampled
thismask_rect = np.ones_like(nsamp,dtype=bool)
#skymask = np.zeros_like(thismask_rect)
slit_left = np.full(nspec_rect,0.0)
slit_righ = np.full(nspec_rect,nspat_rect-1)
tslits_dict = dict(lcen=slit_left,rcen=slit_righ)
sobjs, _ = extract.objfind(imgminsky, thismask_rect, tslits_dict['lcen'], tslits_dict['rcen'],
                                                inmask=outmask, show_peaks=True,show_fits=True, show_trace=True)
sobjs_neg, _ = extract.objfind(-imgminsky, thismask_rect, tslits_dict['lcen'], tslits_dict['rcen'],
                               inmask=outmask, show_peaks=True,show_fits=True, show_trace=True)
sobjs.append_neg(sobjs_neg)
# Local sky subtraction and extraction
skymodel_rect = np.zeros_like(imgminsky)
global_sky = np.zeros_like(imgminsky)
rn2img = np.zeros_like(imgminsky)
objmodel_rect = np.zeros_like(imgminsky)
ivarmodel_rect = np.zeros_like(imgminsky)
extractmask_rect = np.zeros_like(thismask_rect)
par = spectrograph.default_pypeit_par()
# TODO Modify profile fitting so that we can pass in a position image which will allow us to improve the spatial sampling
# in this final extraction step
skymodel_rect[thismask_rect], objmodel_rect[thismask_rect], ivarmodel_rect[thismask_rect], extractmask_rect[thismask_rect] = \
    skysub.local_skysub_extract(imgminsky, sciivar, tilts, waveimg, global_sky, rn2img, thismask_rect,
    tslits_dict['lcen'], tslits_dict['rcen'], sobjs, model_noise=False,bsp=par['scienceimage']['bspline_spacing'],
    sn_gauss=par['scienceimage']['sn_gauss'],inmask=outmask, show_profile=True)
resids_rect = ((imgminsky - skymodel_rect - objmodel_rect)*np.sqrt(ivarmodel_rect))
resids_pix = resids_rect[extractmask_rect]
std_dev = np.std(resids_pix)
mean_resid = np.mean(resids_pix)
sobjs_out = sobjs.copy()
sobjs_out.purge_neg()

sys.exit(-1)

# Need to deal with the read noise image!
loglam_mid = ((loglam_bins + np.roll(loglam_bins,1))/2.0)[1:]
wave_mid = np.power(10.0,loglam_mid)
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
