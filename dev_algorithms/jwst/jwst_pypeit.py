
import os
import numpy as np
import scipy
from matplotlib import pyplot as plt
from astropy.stats import sigma_clipped_stats

from IPython import embed
# set environment variables
os.environ['CRDS_PATH'] = '/Users/joe/crds_cache/jwst_pub'
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds-pub.stsci.edu'
from matplotlib import pyplot as plt
from astropy.io import fits
from gwcs import wcstools


## JWST imports
# The calwebb_spec and spec3 pipelines
from jwst.pipeline import Spec2Pipeline
from jwst.pipeline import Spec3Pipeline
# individual steps
from jwst.assign_wcs import AssignWcsStep
from jwst.background import BackgroundStep
from jwst.imprint import ImprintStep
from jwst.msaflagopen import MSAFlagOpenStep
from jwst.extract_2d import Extract2dStep
from jwst.srctype import SourceTypeStep
from jwst.wavecorr import WavecorrStep
from jwst.flatfield import FlatFieldStep
from jwst.pathloss import PathLossStep
from jwst.barshadow import BarShadowStep
from jwst.photom import PhotomStep
from jwst.resample import ResampleSpecStep
from jwst.extract_1d import Extract1dStep
from jwst import datamodels
DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']


# PypeIt imports
from jwst_utils import compute_diff, get_cuts, show_2dspec, get_jwst_slits
from pypeit.display import display
from pypeit.utils import inverse
from pypeit.core import findobj_skymask
from pypeit.core import skysub, coadd
from pypeit.core import procimg

det = 'nrs1'
#disperser = 'G395M'
disperser='PRISM'
if 'PRISM' in disperser:
    # PRISM data
    rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/Raw'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/output'

    # NIRSPEC 3-point dither
    # dither center
    scifile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + det + '_rate.fits')
    bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + det + '_rate.fits')
    bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + det + '_rate.fits')

    # dither offset
    #scifile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + det + '_rate.fits')
    #bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + det + '_rate.fits')
    #bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + det + '_rate.fits')
elif 'G395M' in disperser:
    # Use islit = 37 for nrs1
    # G395M data
    rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/Raw'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/output'

    # NIRSPEC 3-point dither
    scifile = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + det + '_rate.fits')
    bkgfile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + det + '_rate.fits')
    bkgfile2 = os.path.join(rawpath_level2, 'jw02736007001_03103_00003_' + det + '_rate.fits')


# Plot the 2d differnence image
rawscience, diff = compute_diff(scifile, bkgfile1, bkgfile2, )
#sci_rate = datamodels.open(scifile)


#viewer_diff, ch_diff = display.show_image(diff.T, cuts=get_cuts(diff), chname='diff2d')
#viewer_sci,  ch_sci = display.show_image(sci.T, cuts=get_cuts(sci), chname='raw', wcs_match=True)
basename = os.path.basename(scifile).replace('rate.fits', '')

# TODO Should we flat field. The flat field and flat field error are wonky and probably nonsense
param_dict = {
    'extract_2d': {'save_results': True},
    'bkg_subtract': {'skip': True},
    'master_background_mos': {'skip': True},
    'srctype': {'source_type':'EXTENDED'},
    'flat_field': {'skip': True},
    'resample_spec': {'skip': True},
    'extract_1d': {'skip': True},
    #    'flat_field': {'save_interpolated_flat': True},
}


runflag = True
if runflag:
    spec2 = Spec2Pipeline(steps=param_dict)
    spec2.save_results = True
    spec2.output_dir = output_dir
    result = spec2(scifile)

# Read in the files
#intflat_output_file = os.path.join(output_dir, basename + 'interpolatedflat.fits')
e2d_output_file = os.path.join(output_dir, basename + 'extract_2d.fits')
cal_output_file = os.path.join(output_dir, basename + 'cal.fits')
s2d_output_file = os.path.join(output_dir, basename + 's2d.fits')
# TESTING
#final2d = datamodels.open(s2d_output_file)
#intflat = None
e2d = datamodels.open(e2d_output_file)
final2d = datamodels.open(cal_output_file)
#intflat = datamodels.open(intflat_output_file)
intflat = None
islit = 10
#islit = 37
#islit =18

show_2dspec(rawscience, final2d, islit, intflat=intflat, emb=False, clear=True)

# Try to reverse engineer all the things they multiply into the data
slit_name = final2d.slits[islit].name
pathloss = np.array(final2d.slits[islit].pathloss_uniform.T, dtype=float) \
    if final2d.slits[islit].source_type == 'EXTENDED' else np.array(final2d.slits[islit].pathloss_point.T, dtype=float)
#flat = np.array(intflat.slits[islit].data.T, dtype=float)
flat = np.ones_like(pathloss)
barshadow = np.array(final2d.slits[islit].barshadow.T, dtype=float)
photom_conversion = final2d.slits[islit].meta.photometry.conversion_megajanskys

# This is the conversion between final2d and e2d, i.e. final2d = jwst_scale*e2d
jwst_scale = photom_conversion/flat/pathloss/barshadow
# I think there is a logical inconsistency here associated with the flat, since it works great without flat fielding
count_scale = flat*pathloss*barshadow # These are the things that the raw rate data were divided by
t_eff = final2d.slits[islit].meta.exposure.effective_exposure_time # TODO I don't konw what this means! Find out

# Let's get these images into counts so that PypeIt will make sense
flux_to_counts = t_eff/photom_conversion
science = np.array(final2d.slits[islit].data.T, dtype=float)*flux_to_counts
# TESTING!!  kludge the error by multiplying by a small number
#kludge_err = 0.1667
#kludge_err = 0.28
#kludge_err =0.36
kludge_err = 1.0
# TODO Currently the var_flat is nonsense I think and so I'm just going to use the var_poisson and var_rnoise to get
# the noise.
#err = kludge_err*np.array(final2d.slits[islit].err.T, dtype=float)*flux_to_counts
var_poisson = final2d.slits[islit].var_poisson.T*flux_to_counts**2
var_rnoise = final2d.slits[islit].var_rnoise.T*flux_to_counts**2
var = kludge_err*np.array(var_poisson + var_rnoise, dtype=float)
# This needs to be multiplied by count_scale to get it into units of counts which is what pypeit requires. I checked
# that this base_var is equal to e2d.var_rnoise if you remove the flux_to_counts factor.
base_var = np.array(final2d.slits[islit].var_rnoise.T, dtype=float)*flux_to_counts**2*count_scale**2
dq = np.array(final2d.slits[islit].dq.T, dtype=int)
waveimg = np.array(final2d.slits[islit].wavelength.T, dtype=float)

slit_wcs = final2d.slits[islit].meta.wcs
x, y = wcstools.grid_from_bounding_box(slit_wcs.bounding_box, step=(1, 1))
calra, caldec, calwave = slit_wcs(x, y)
ra = calra.T

gpm = np.logical_not(dq & DO_NOT_USE)
thismask = np.isfinite(science)
nanmask = np.logical_not(thismask)
science[nanmask]= 0.0
#err[nanmask] = 0.0
var[nanmask] = 0.0
sciivar = inverse(var)*gpm
base_var[nanmask] = 0.0
count_scale[nanmask] = 0.0

# Wave nanmask is different from data nanmask
nanmask_wave = np.logical_not(np.isfinite(waveimg))
wave_min = np.min(waveimg[np.logical_not(nanmask_wave)])
wave_max = np.max(waveimg[np.logical_not(nanmask_wave)])
nanmask_ra = np.logical_not(np.isfinite(ra))
ra_min = np.min(ra[np.logical_not(nanmask_ra)])
ra_max = np.max(ra[np.logical_not(nanmask_ra)])
waveimg[nanmask_wave] = 0.0
ra[nanmask_ra]=0.0

nspec, nspat = science.shape


slit_left, slit_righ = get_jwst_slits(thismask)
# Generate some tilts and a spatial image
tilts = np.zeros_like(waveimg)
tilts[np.isfinite(waveimg)] = (waveimg[np.isfinite(waveimg)] - wave_min)/(wave_max-wave_min)
# TODO Fix this spat_pix to make it increasing with pixel. For now don't use
spat_pix = (ra - ra_min)/(ra_max - ra_min)*(nspat-1)
spat_pix[nanmask_ra] = 0.0


trim_edg = (0,0)
boxcar_rad_pix = 8.0
fwhm = 2.0
bsp = 2.5
sn_gauss=3.0
nperslit = 2
no_poly=False
ncoeff = 5
snr_thresh = 5.0
pos_mask = False
model_noise = True

# First pass sky-subtraction and object finding
initial_sky0 = np.zeros_like(science)
initial_sky0[thismask] = skysub.global_skysub(science, sciivar, tilts, thismask, slit_left, slit_righ,
                                              inmask = gpm, bsp=bsp, pos_mask=pos_mask, no_poly=no_poly, show_fit=True,
                                              trim_edg=trim_edg)
sobjs_slit0 = findobj_skymask.objs_in_slit(science-initial_sky0, sciivar, thismask, slit_left, slit_righ, inmask=gpm, ncoeff=ncoeff,
                                          snr_thresh=snr_thresh, show_peaks=True, show_trace=True,
                                          trim_edg=trim_edg,  fwhm=fwhm, boxcar_rad=boxcar_rad_pix, maxdev = 2.0, find_min_max=None,
                                          qa_title='objfind_QA_' + slit_name, nperslit=nperslit,
                                          objfindQA_filename=None, debug_all=False)
# Create skymask and perfrom second pass sky-subtraction and object finding
skymask = np.ones_like(thismask)
skymask[thismask] = findobj_skymask.create_skymask(sobjs_slit0, thismask,
                                                   slit_left, slit_righ,
                                                   trim_edg=trim_edg) #, box_rad_pix=boxcar_rad_pix,)

# TODO On the second pass update the variance model like in pypeit
initial_sky = np.zeros_like(science)
initial_sky[thismask] = skysub.global_skysub(science, sciivar, tilts, thismask, slit_left, slit_righ,
                                             inmask = (gpm & skymask), bsp=bsp, pos_mask=pos_mask, no_poly=no_poly, show_fit=True,
                                             trim_edg=trim_edg)
# TESTING Let's try to update the variance model and see if it makes sense


sobjs_slit = findobj_skymask.objs_in_slit(science-initial_sky, sciivar, thismask, slit_left, slit_righ, inmask=gpm, ncoeff=ncoeff,
                                          snr_thresh=snr_thresh, show_peaks=True, show_trace=True,
                                          trim_edg=trim_edg,  fwhm=fwhm, boxcar_rad=boxcar_rad_pix, maxdev = 2.0, find_min_max=None,
                                          qa_title='objfind_QA_' + slit_name, nperslit=nperslit,
                                          objfindQA_filename=None, debug_all=False)

viewer, ch = display.show_image(science-initial_sky, cuts=get_cuts(science), chname='science', wcs_match=True)
display.show_slits(viewer, ch, slit_left, slit_righ, pstep=1, slit_ids=np.array([int(slit_name)]))
for spec in sobjs_slit:
    if spec.hand_extract_flag == False:
        color = 'orange'
    else:
        color = 'blue'
    display.show_trace(viewer, ch, spec.TRACE_SPAT, trc_name=spec.NAME, color=color)

#initial_sky = initial_sky*0.0
# Local sky subtraction and extraction
skymodel = initial_sky.copy()
objmodel = np.zeros_like(science)
ivarmodel = sciivar.copy()
extractmask = gpm.copy()
sobjs = sobjs_slit.copy()
# TODO Need to figure out what the count_scale is order to update the noise in the low-background regime.
skymodel[thismask], objmodel[thismask], ivarmodel[thismask], extractmask[thismask] = skysub.local_skysub_extract(
    science, sciivar, tilts, waveimg, initial_sky, thismask, slit_left, slit_righ, sobjs, ingpm = gpm,
    base_var = base_var, count_scale=count_scale, bsp = bsp, sn_gauss = sn_gauss, trim_edg=trim_edg, spat_pix=None, model_full_slit=True,
    model_noise=model_noise,  show_profile=True, show_resids=True, debug_bkpts=True)

# This checks out
#adderr = 0.01
#var_check = procimg.variance_model(base_var, counts=(skymodel + objmodel), count_scale=count_scale, noise_floor=adderr)

n_bins = 50
sig_range = 7.0
binsize = 2.0 * sig_range / n_bins
bins_histo = -sig_range + np.arange(n_bins) * binsize + binsize / 2.0

xvals = np.arange(-10.0, 10, 0.02)
gauss = scipy.stats.norm(loc=0.0, scale=1.0)
gauss_corr = scipy.stats.norm(loc=0.0, scale=1.0)

chi = (science-objmodel-skymodel)*np.sqrt(ivarmodel)
chi = chi.flatten()
maskchi = extractmask.flatten()
sigma_corr, maskchi = coadd.renormalize_errors(chi, maskchi, max_corr = 20.0, title='jwst_sigma_corr', debug=True)

plt.plot(sobjs[0].OPT_WAVE,sobjs[0].OPT_COUNTS, drawstyle='steps-mid', color='black', label='optimal')
plt.plot(sobjs[0].OPT_WAVE,sobjs[0].OPT_COUNTS_SIG, drawstyle='steps-mid', color='red', label='optimal')
plt.plot(sobjs[0].BOX_WAVE,sobjs[0].BOX_COUNTS, drawstyle='steps-mid', color='green', label='boxcar')
plt.plot(sobjs[0].BOX_WAVE,sobjs[0].BOX_COUNTS_SIG, drawstyle='steps-mid', color='blue', label='boxcar')
plt.legend()
plt.show()


# Now play around with a PypeIt extraction



