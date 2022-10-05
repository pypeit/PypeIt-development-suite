import os
import numpy as np
import scipy
from astropy.io import fits
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
from jwst_utils import compute_diff, get_cuts, jwst_nircam_proc, jwst_nircam_subimgs
from pypeit.display import display
from pypeit import specobjs
from pypeit import slittrace
from pypeit.utils import inverse, fast_running_median
from pypeit.core import findobj_skymask
from pypeit.core import skysub, coadd
from pypeit.core import procimg
from pypeit.core import flat

from pypeit.spectrographs.util import load_spectrograph
from pypeit.images import pypeitimage
from pypeit import calibrations
from pypeit import find_objects
from pypeit import extraction
from pypeit import msgs
from pypeit import spec2dobj
from pypeit import coadd2d
from pypeit.images.pypeitimage import PypeItImage
from pypeit.scripts.show_2dspec import show_trace


DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']

# NIRCAM 2 flat field exposures
rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRCAM_WFSS/GTO/'
output_dir = '/Users/joe/jwst_redux/redux/NIRCAM_WFSS/GTO/calwebb/output/'
pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRCAM_WFSS/GTO/calwebb/pypeit/'

scifile1 = os.path.join(rawpath_level2, 'jw01243001001_02101_00001_nrcalong_rate.fits')
scifile2 = os.path.join(rawpath_level2, 'jw01243001001_02101_00002_nrcalong_rate.fits')
scifile3 = os.path.join(rawpath_level2, 'jw01243001001_02101_00003_nrcalong_rate.fits')
scifiles = [scifile1, scifile2, scifile3]


# TODO: This needs to be defined by the user
scipath = os.path.join(pypeit_output_dir, 'Science')
if not os.path.isdir(scipath):
    msgs.info('Creating directory for Science output: {0}'.format(scipath))
    os.makedirs(scipath)

# NIRCAM Mirage configuration file
condir = '/Users/joe/jwst_redux/Raw/NIRCAM_WFSS/GTO/sens/'
filt = 'F356W'
module='a'
grism = 'R'
configfile = os.path.join(condir, 'NIRCAM_' + filt + '_mod' + module.upper() + '_' + grism + '.conf')
detname ='a'
# print(configfile)


# TODO Should we flat field. The flat field and flat field error are wonky and probably nonsense
param_dict = {
    'extract_2d': {'skip': True},
    'bkg_subtract': {'skip': True},
    'flat_field': {'save_interpolated_flat': True},
    'wfss_contam': {'skip': True},
    'photom': {'skip': True},
    'extract_1d': {'skip': True},
}

basenames = []
for sci in scifiles:
    basenames.append(os.path.basename(sci).replace('_rate.fits', ''))

# Run the spec2 pipeline
runflag = False
if runflag:
    for sci in scifiles:
        step = AssignWcsStep()
        step.save_results = True
        step.output_dir = output_dir
        result = step.run(sci)

        step = FlatFieldStep()
        step.output_dir = output_dir
        # step.override_flat = new_flat_file
        step.save_results = True
        result = step.run(result)

#        spec2 = Spec2Pipeline(steps=param_dict)
#        spec2.save_results = True
#        spec2.output_dir = output_dir
#        result = spec2(sci)


# Some pypeit things
spectrograph = load_spectrograph('jwst_nircam')
par = spectrograph.default_pypeit_par()
det_container = spectrograph.get_detector_par(1)

pypeline = 'MultiSlit'
par['rdx']['redux_path'] = pypeit_output_dir
qa_dir = os.path.join(pypeit_output_dir, 'QA')
par['rdx']['qadir'] = 'QA'
png_dir = os.path.join(qa_dir, 'PNGs')
if not os.path.isdir(qa_dir):
    msgs.info('Creating directory for QA output: {0}'.format(qa_dir))
    os.makedirs(qa_dir)
if not os.path.isdir(png_dir):
    os.makedirs(png_dir)

# Set some parameters for difference imaging
#if diff_redux:
#    par['reduce']['findobj']['skip_skysub'] = True # Do not sky-subtract when object finding
 #   par['reduce']['extraction']['skip_optimal'] = True # Skip local_skysubtraction and profile fitting



# TODO Fix this, currently does not work if target names have - or _
filename_first = os.path.basename(scifiles[0])
filename_last = os.path.basename(scifiles[-1])
split_first = filename_first.split('_')
split_last = filename_last.split('_')
prefix_first = ''
for ii in range(len(split_first) - 3):
    prefix_first += split_first[ii] + "_"
out_filename = prefix_first + split_first[-3] + "-" + split_last[-3]

show = True
show_peaks=True
show_skysub_fit=True

flat_output_files = []
for base in basenames:
    flat_output_files.append(os.path.join(output_dir, base + '_flatfieldstep.fits'))

RA, DEC =  15.05428, 28.040513
rate_obj, science, sciivar, gpm, dq_gpm, base_var, count_scale, finitemask, tilts, waveimg, spat_img = jwst_nircam_proc(
    flat_output_files[0], configfile, RA, DEC, kludge_err=1.15, noise_floor=0.01, saturation=65000)


# Create PypeIt sciimg object
sciImg = pypeitimage.PypeItImage(image=science,ivar=sciivar, base_var=base_var, img_scale=count_scale,
                                 detector=det_container, bpm=np.logical_not(gpm).astype(int))
# Create bogus rectilinear slit boundaries -- replace with a tracing model when that finally works
nspec, nspat = science.shape
slit_left = np.full(nspec, 0.0, dtype=float)
slit_righ = np.full(nspec, nspat, dtype=float)

# Create some bogus slit boundaries
slits = slittrace.SlitTraceSet(slit_left, slit_righ, pypeline, detname=detname, nspat=nspat,
                               PYP_SPEC=spectrograph.name)
slitmask = slits.slit_img()

sciImg.build_mask(slitmask=slitmask)

if show:
    display.clear_all()

if show:
    display.show_image(rate_obj.data, chname='raw rate')
    viewer_data, ch_data = display.show_image(sciImg.image, waveimg=waveimg, cuts=get_cuts(science), chname='pypeit image')
    display.show_slits(viewer_data, ch_data, slit_left, slit_righ, pstep=1, slit_ids=np.array(['this_slit']))
    viewer_wave, ch_wave = display.show_image(waveimg, chname='wave')
    viewer_tilts, ch_tilts = display.show_image(tilts, waveimg=waveimg, chname='tilts')
    #viewer_flat, ch_flat = display.show_image(flatfield, waveimg=waveimg, chname='flat_{:s}'.format(detectors[idet]))
    viewer_spat, ch_spat = display.show_image(spat_img, waveimg=waveimg, chname='spat_img')


# Insantiate FindObjects object
objFind = find_objects.FindObjects.get_instance(sciImg, slits, spectrograph, par, 'science_coadd2d', tilts=tilts,
                                                basename=basenames[0], clear_ginga=False, show=show)

#if show:
#    gpm = sciImg.select_flag(invert=True)
#    objFind.show('image', image=sciImg.image*gpm.astype(float), chname='science', slits=True)

global_sky0, sobjs_obj = objFind.run(show_peaks=show or show_peaks, show_skysub_fit=show_skysub_fit)

# TODO Add the option to update global sky as an option to FindObjects
skymask = objFind.create_skymask(sobjs_obj)
global_sky = objFind.global_skysub(previous_sky=global_sky0, skymask=skymask, show=True)

# Initiate Extract object
extract = extraction.Extract.get_instance(sciImg, slits, sobjs_obj, spectrograph, par, 'science_coadd2d',
                                          tilts=tilts, waveimg=waveimg, basename=basenames[0], show=show)


if not par['reduce']['extraction']['skip_extraction']:
    skymodel, objmodel, ivarmodel, outmask, sobjs, waveimg, tilts = extract.run(global_sky, sobjs_obj)
else:
    # Although exrtaction is not performed, still need to prepare some masks and the tilts
    # self.exTract.prepare_extraction()
    # Since the extraction was not performed, fill the arrays with the best available information
    skymodel = global_sky
    objmodel = np.zeros_like(extract.sciImg.image)
    ivarmodel = np.copy(extract.sciImg.ivar)
    outmask = extract.sciImg.fullmask
    sobjs = sobjs_obj

# TODO -- Do this upstream
# Tack on detector
for sobj in sobjs:
    sobj.DETECTOR = sciImg.detector


# Construct the Spec2DObj with the positive image
spec2DObj = spec2dobj.Spec2DObj(sciimg=sciImg.image,
                                ivarraw=sciImg.ivar,
                                skymodel=skymodel,
                                objmodel=objmodel,
                                ivarmodel=ivarmodel,
                                scaleimg=None,
                                waveimg=waveimg,
                                bpmmask=outmask,
                                detector=sciImg.detector,
                                sci_spat_flexure=sciImg.spat_flexure,
                                sci_spec_flexure=None,
                                vel_corr=None,
                                vel_type=None,
                                tilts=tilts,
                                slits=slits,
                                wavesol=None,
                                maskdef_designtab=None)
spec2DObj.process_steps = sciImg.process_steps

# QA
spec2DObj.gen_qa()

# Make a plot of the residuals for a random slit
chi = (science-objmodel-skymodel)*np.sqrt(ivarmodel)*outmask

islit = 0
slitmask = slits.slit_img(initial=False, flexure=None, exclude_flag=None)
slitord_id = slits.slitord_id[islit]
thismask = slitmask == slitord_id

n_bins = 50
sig_range = 7.0
binsize = 2.0 * sig_range / n_bins
bins_histo = -sig_range + np.arange(n_bins) * binsize + binsize / 2.0

xvals = np.arange(-10.0, 10, 0.02)
gauss = scipy.stats.norm(loc=0.0, scale=1.0)
gauss_corr = scipy.stats.norm(loc=0.0, scale=1.0)

maskchi = thismask & outmask
sigma_corr, maskchi = coadd.renormalize_errors(chi, maskchi, max_corr = 20.0, title='jwst_sigma_corr', debug=True)


#
#
#
# rate_obj = datamodels.open(flat_output_files[0])
#
#
# #TODO Fix all these numbers to be correct
# kludge_err = 1.0
# ronoise=5.17
# saturation=65000
# noise_floor=0.01
#
# RA, DEC =  15.05428, 28.040513
# t_eff = rate_obj.meta.exposure.effective_exposure_time
#
# # TODO Use the spat_img which PypeIt can take, but we are not using!!
# raw_sub, var_tot_sub, var_poisson_sub, var_rnoise_sub, dq_sub, waveimg_sub, spat_img_sub  = jwst_nircam_subimgs(configfile, RA, DEC, flat_output_files[0])
#
#
# # slit_slice, slit_left, slit_righ, slit_left_orig, slit_righ_orig, spec_vals_orig, src_trace_ra, src_trace_dec, dq, \
# # ra, dec, waveimg, tilts, flatfield, pathloss, barshadow, photom_conversion, final = jwst_extract_subimgs(
# #    final_slit, intflat_slit)
#
# # Now deal with the image processing
# # if not np.any(finitemask):
# #    return (None,)*21
#
# # Read in the output after msa_flagging. Extract the sub-images, rotate to PypeIt format.
# rate = raw_sub.T
# rate_var_tot = var_tot_sub.T
# rate_var_poisson = var_poisson_sub.T
# rate_var_rnoise = var_rnoise_sub.T
# # TODO Check that the errors don't have nonsense from the flat field error budget
# dq = dq_sub.T
# # Now perform the image processing
# raw_counts = rate * t_eff
# raw_var_poisson = kludge_err ** 2 * rate_var_poisson * t_eff ** 2
# raw_var_rnoise = kludge_err ** 2 * rate_var_rnoise * t_eff ** 2
# # Is this correct? I'm not sure I should be using their poisson variance for the noise floor
# raw_var = procimg.variance_model(raw_var_rnoise, counts=raw_var_poisson, noise_floor=noise_floor)
# # TODO This  is a hack until I can understand how to get rid of the hot pixels in the JWST variance arrays using DQ flags.
# # I don't know what the value of this parameter currently set to 20 should be?? Look into this via a github issue.
# # raw_gpm = (raw_var_rnoise < 20.0*ronoise**2) & (raw_var_poisson < saturation)
# raw_gpm = (raw_var_rnoise < saturation) & (raw_var_poisson < saturation)
# # raw_var_poisson + raw_var_rnoise # TODO Leaving out problematic flat field term from pipeline
#
#
# # This is the conversion between final2d and e2d, i.e. final2d = jwst_scale*e2d
# # total_flat = flatfield*pathloss*barshadow
# # flux_to_counts = t_eff / photom_conversion  # This converts s2d outputs of flux to counts.
# # jwst_scale = photom_conversion/flatfield/pathloss/barshadow
#
# # TODO Perform the flat field correction yourself so as to update the noise model. Setting this to unit
# total_flat = np.ones_like(rate)
# finitemask = np.isfinite(rate, dtype=bool) # This was from NIRSPEC and may not be necessary
# total_flat_square = np.square(total_flat)
#
# count_scale = inverse(total_flat)  # This is the quantity that goes into PypeIt for var modeling
# science, flat_bpm = flat.flatfield(raw_counts, total_flat)
# var_poisson, _ = flat.flatfield(raw_var_poisson, total_flat_square)
# base_var, _ = flat.flatfield(raw_var_rnoise, total_flat_square)
# var, _ = flat.flatfield(raw_var, total_flat_square)
# sciivar = inverse(var)
# dq_gpm = np.logical_not(dq & DO_NOT_USE)
# gpm = finitemask & dq_gpm & np.logical_not(flat_bpm) & (sciivar > 0.0) & raw_gpm
#
# nanmask = np.logical_not(finitemask)
# count_scale[nanmask] = 0.0
# science[nanmask] = 0.0
# var_poisson[nanmask] = 0.0
# base_var[nanmask] = 0.0
# var[nanmask] = 0.0
# sciivar[nanmask] = 0.0