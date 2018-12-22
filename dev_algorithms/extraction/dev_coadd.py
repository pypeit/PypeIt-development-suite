
#

import inspect
import copy

import numpy as np

from scipy import interpolate

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


from pypeit import ginga

# Read in the masters
path = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT_old/Gemini_GNIRS/GNIRS/'
master_dir = path + 'MF_gemini_gnirs/'
scidir = path + 'Science/'
waveimgfiles = [os.path.join(master_dir, ifile) for ifile in ['MasterWave_A_1_01.fits','MasterWave_A_2_01.fits', 'MasterWave_A_4_01.fits',
                                                              'MasterWave_A_8_01.fits']]
tiltfiles = [os.path.join(master_dir, ifile) for ifile in ['MasterTilts_A_1_01.fits','MasterTilts_A_2_01.fits', 'MasterTilts_A_4_01.fits',
                                                           'MasterWave_A_8_01.fits']]
tracefile = os.path.join(master_dir,'MasterTrace_A_15_01')
spec2d_files = [os.path.join(scidir, ifile) for ifile in ['spec2d_pisco_GNIRS_2017Mar31T085412.181.fits',
                                                              'spec2d_pisco_GNIRS_2017Mar31T085933.097.fits',
                                                              'spec2d_pisco_GNIRS_2017Mar31T091538.731.fits',
                                                              'spec2d_pisco_GNIRS_2017Mar31T092059.597.fits']]
spec1d_files = [os.path.join(scidir, ifile) for ifile in ['spec1d_pisco_GNIRS_2017Mar31T085412.181.fits',
                                                              'spec1d_pisco_GNIRS_2017Mar31T085933.097.fits',
                                                              'spec1d_pisco_GNIRS_2017Mar31T091538.731.fits',
                                                              'spec1d_pisco_GNIRS_2017Mar31T092059.597.fits']]

# Define the new wavelength grid
ngrid = 5000
dloglam = 0.000127888 # this is the average of the median dispersions
logmin = 3.777
loglam = logmin + dloglam*np.arange(ngrid)
# Define the oversampled grid
nosamp = 1
loglam_osamp = logmin + (dloglam/nosamp)*np.arange(nosamp*ngrid)

nexp = len(waveimgfiles)
# Read in the tslits_dict
Tslits = traceslits.TraceSlits.from_master_files(tracefile)
tslits_dict = Tslits._fill_tslits_dict()
spectrograph = load_spectrograph('gemini_gnirs')

slitmask = spectrograph.slitmask(tslits_dict)
nslits = tslits_dict['lcen'].shape[1]
slit = 5
slitmask_stack = np.einsum('i,jk->ijk', np.ones(nexp), slitmask)

for ifile in range(nexp):
    hdu = fits.open(waveimgfiles[ifile])
    waveimg = hdu[0].data
    hdu_tilts = fits.open(tiltfiles[ifile])
    tilts = hdu[0].data
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
        sky_stack = np.zeros((nexp,sciimg.shape[0],sciimg.shape[1]),dtype=float)
        sciivar_stack = np.zeros((nexp,sciimg.shape[0],sciimg.shape[1]),dtype=float)
        mask_stack = np.zeros((nexp,sciimg.shape[0],sciimg.shape[1]),dtype=float)
        trace_stack = np.zeros((nexp,trace.shape[0]))
        wave_stack = np.zeros((nexp,trace.shape[0]))

    waveimg_stack[ifile,:,:] = waveimg
    tilts_stack[ifile,:,:] = tilts
    sciimg_stack[ifile,:,:] = sciimg
    sciivar_stack[ifile,:,:] = sciivar
    mask_stack[ifile,:,:] = mask
    sky_stack[ifile,:,:] = sky
    trace_stack[ifile,:] = trace
    wave_stack[ifile,:] = wave_trace

nspec, nspat = sciimg_stack[0,:,:].shape

# Determine the wavelength grid that we will use for the current slit/order
thismask_stack = slitmask_stack == slit
thismask = thismask_stack[0,:,:]
wave_min = waveimg_stack[thismask_stack].min()
wave_max = waveimg_stack[thismask_stack].max()
diff = loglam_osamp - np.log10(wave_min)
diff[diff > 0] = 1e10
ind_lower = np.argmin(np.abs(diff))
diff = np.log10(wave_max) - loglam_osamp
diff[diff > 0] = 1e10
ind_upper = np.argmin(np.abs(diff))
loglam_bins = loglam_osamp[ind_lower] + (dloglam/nosamp)*np.arange(ind_upper-ind_lower + 1)
wave_bins = np.power(10.0,loglam_bins)
nwave = wave_bins.shape[0]

spat_img = np.outer(np.ones(nspec), np.arange(nspat))
#pad_spat = 5
#spec_ind, spat_ind = np.where(thismask)
#min_spat = spat_ind.min() - pad_spat
#max_spat = spat_ind.max() + pad_spat
#num_spat = max_spat - min_spat + 1

# Need to determine the limits of the spat img
#spat_lin = np.linspace(min_spat,max_spat,num = int(np.round(num_spat)))
#spat_img = np.outer(np.ones(nwave), spat_lin)
#spat_img, loglam_img = np.meshgrid(spat_lin, loglam_slit)

# Normalized spatial offset image (from object trace)
spat_min = 1e10
spat_max = -1e10
for iexp in range(nexp):
#    slit_cen_lin = (interpolate.interp1d(wave_stack[iexp,:],trace_stack[iexp,:],bounds_error=False ,fill_value='extrapolate'))(wave_slit)
    slit_cen_img = np.outer(trace_stack[iexp,:], np.ones(nspat))  # center of the slit replicated spatially
    dspat_img = (spat_img - slit_cen_img)
    spat_min = np.fmin(spat_min,dspat_img[thismask].min())
    spat_max = np.fmax(spat_max,dspat_img[thismask].max())

spat_min_int = int(np.floor(spat_min))
spat_max_int = int(np.ceil(spat_max))
dspat_bins = np.arange(spat_min_int,spat_max_int +1,1)

nspec_out = nwave -1
nspat_out = dspat_bins.shape[0]-1
shape = (nexp, nspec_out, nspat_out)
sciimg_out_stack = np.zeros(shape)
imgminsky_out_stack = np.zeros(shape)
var_out_stack = np.zeros(shape)
tilts_out_stack = np.zeros(shape)
waveimg_out_stack = np.zeros(shape)
norm_out_stack = np.zeros(shape)


for iexp in range(nexp):
    slit_cen_img = np.outer(trace_stack[iexp,:], np.ones(nspat))  # center of the slit replicated spatially
    dspat_img = (spat_img - slit_cen_img)
    finmask = thismask & (mask_stack[iexp,:,:] == 0)
    spec_rebin = np.log10(waveimg_stack[iexp,:,:][finmask])
    spat_rebin = dspat_img[finmask]
    norm_img, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                  bins=[loglam_bins, dspat_bins], density=False)
    norm_out_stack[iexp,:,:] = norm_img
    # Rebin the science image
    img_now = sciimg_stack[iexp,:,:]
    weigh_img, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                   bins=[loglam_bins, dspat_bins], density=False,
                                                   weights = img_now[finmask])
    sciimg_out_stack[iexp,:,:] =(norm_img > 0.0)*weigh_img/(norm_img + (norm_img == 0.0))
    # Rebin the sky subtracted image
    imgminsky_now = (sciimg_stack[iexp,:,:] - sky_stack[iexp,:,:])
    weigh_imgminsky, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                   bins=[loglam_bins, dspat_bins], density=False,
                                                   weights = imgminsky_now[finmask])
    imgminsky_out_stack[iexp,:,:] =(norm_img > 0.0)*weigh_imgminsky/(norm_img + (norm_img == 0.0))
    # Propagate the variance
    var_now = utils.calc_ivar(sciivar_stack[iexp,:,:])
    weigh_var, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                   bins=[loglam_bins, dspat_bins], density=False,
                                                   weights = var_now[finmask])
    var_out_stack[iexp,:,:] =(norm_img > 0.0)*weigh_var/(norm_img + (norm_img == 0.0))**2
    # Rebin the tilts
    tilts_now = tilts_stack[iexp,:,:]
    weigh_tilts, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                   bins=[loglam_bins, dspat_bins], density=False,
                                                   weights = tilts_now[finmask])
    tilts_out_stack[iexp,:,:] =(norm_img > 0.0)*weigh_tilts/(norm_img + (norm_img == 0.0))
    # Rebin the wavelength image
    waveimg_now = tilts_stack[iexp,:,:]
    weigh_waveimg, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                   bins=[loglam_bins, dspat_bins], density=False,
                                                   weights = waveimg_now[finmask])
    waveimg_out_stack[iexp,:,:] =(norm_img > 0.0)*weigh_waveimg/(norm_img + (norm_img == 0.0))

# Now compute the final stack with sigma clipping
sigrej = 3.0
maxiters = 10
data = np.ma.MaskedArray(sciimg_out_stack, (norm_out_stack > 0))
sigclip = stats.SigmaClip(sigma=sigrej, maxiters=maxiters, cenfunc='median')
data_clipped = sigclip(data, axis=0, masked=True)
mask_out_stack = np.invert(data_clipped.mask)  # outmask = True are good values
# Replace this later with S/N weights
weights = np.ones(nexp)/float(nexp)
weights_stack = np.einsum('i,ijk->ijk', weights, mask_out_stack)
weights_sum = np.sum(weights_stack, axis=0)
sciimg = np.sum(sciimg_out_stack*weights_stack,axis=0)/(weights_sum + (weights_sum == 0.0))
varfinal = np.sum(var_out_stack * weights_stack ** 2, axis=0) / (weights_sum + (weights_sum == 0.0)) ** 2
sciivar = utils.calc_ivar(varfinal)
# Need to deal with the read noise image!