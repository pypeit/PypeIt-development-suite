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
from jwst_utils import compute_diff, get_cuts, jwst_show_spec2, jwst_show_msa, jwst_proc, jwst_extract_subimgs
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

# detname = 'nrs1'
# detector = 1 if 'nrs1' in detname else 2

disperser = 'G395M_Maseda'
#disperser = 'G395M'
#disperser = 'PRISM_01117'
# disperser = 'G235M'
#disperser='PRISM_01133'
# detectors = ['nrs1', 'nrs2']
# disperser='PRISM_01117'
# disperser='PRISM_FS'
mode = 'MSA'
# mode ='FS'
detectors = ['nrs1', 'nrs2']
exp_list = []
diff_redux = True
# If diff_redux is False, the code will model the sky and the object profile and perform optimal extraction.
# If diff_redux is True, the code will difference image and simply boxcar extract (optimal not implemented yet)
for detname in detectors:
    # TODO add the kendrew FS SN data to this.
    if 'PRISM_01133' == disperser:
        ## Prorgram for Slit Loss Characterization for MSA shutters
        # PRISM data
        rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/Raw'
        output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/output'
        pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/pypeit'

        # NIRSPEC 3-point dither
        # dither center
        scifile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_rate.fits')
        scifile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_rate.fits')
        scifile3 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_rate.fits')

        # dither offset
        # scifile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_rate.fits')
        # bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_rate.fits')
        # bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_rate.fits')
    if 'PRISM_FS' == disperser:
        ## Prorgram for Slit Loss Characterization for MSA shutters
        # PRISM data
        rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/2072/level_12'
        output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/02027_PRISM/calwebb/output'
        pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/02027_PRISM/calwebb/pypeit'

        # NIRSPEC 3-point dither
        # dither center
        scifile1 = os.path.join(rawpath_level2, 'jw02072002001_05101_00001_' + detname + '_rate.fits')
        scifile2 = os.path.join(rawpath_level2, 'jw02072002001_05101_00002_' + detname + '_rate.fits')
        scifile3 = os.path.join(rawpath_level2, 'jw02072002001_05101_00003_' + detname + '_rate.fits')

    elif 'PRISM_01117' in disperser:
        # PRISM data
        rawpath_level2 = '//Users/joe/jwst_redux/Raw/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/level_12/01117'
        output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb/output'
        pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb/pypeit'

        # NIRSPEC 3-point dither
        # dither center
        scifile1 = os.path.join(rawpath_level2, 'jw01117007001_03101_00002_' + detname + '_rate.fits')
        scifile2 = os.path.join(rawpath_level2, 'jw01117007001_03101_00003_' + detname + '_rate.fits')
        scifile3 = os.path.join(rawpath_level2, 'jw01117007001_03101_00004_' + detname + '_rate.fits')

    elif 'G395M' == disperser:
        # Use islit = 37 for nrs1
        # G395M data
        rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/Raw'
        output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/output'
        pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/pypeit'

        # NIRSPEC 3-point dither
        scifile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + detname + '_rate.fits')
        scifile2 = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + detname + '_rate.fits')
        scifile3 = os.path.join(rawpath_level2, 'jw02736007001_03103_00003_' + detname + '_rate.fits')
    elif 'G395M_Maseda' == disperser:
        # Use islit = 37 for nrs1
        # G395M data
        rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/Maseda/'
        output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/Maseda/395M/calwebb/output'
        pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/Maseda/395M/calwebb/pypeit'

        # NIRSPEC 3-point dither
        scifile1 = os.path.join(rawpath_level2, 'jw01671001001_03101_00002_' + detname + '_rate.fits')
        scifile2 = os.path.join(rawpath_level2, 'jw01671001001_03101_00003_' + detname + '_rate.fits')
        scifile3 = os.path.join(rawpath_level2, 'jw01671001001_03101_00004_' + detname + '_rate.fits')

    elif 'G235M' == disperser:
        # Use islit = 38 for nrs1
        # G235M data
        rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/level_2/'
        output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G235M/calwebb/output'
        pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G235M/calwebb/pypeit/'

        # NIRSPEC 3-point dither
        scifile1 = os.path.join(rawpath_level2, 'jw02736007001_03101_00002_' + detname + '_rate.fits')
        scifile2 = os.path.join(rawpath_level2, 'jw02736007001_03101_00003_' + detname + '_rate.fits')
        scifile3 = os.path.join(rawpath_level2, 'jw02736007001_03101_00004_' + detname + '_rate.fits')
    exp_list.append([scifile1, scifile2, scifile3])

if 'MSA' in mode:
    offsets_pixels_list = [[0, 5.0, -5.0], [0, 5.0, -5.0]]
elif 'FS' in mode:
    offsets_pixels_list = [[0, 8.0, 18.0], [0, 8.0, 18.0]]

scifiles_1 = exp_list[0]
scifiles_2 = exp_list[1] if len(exp_list) > 1 else []
scifiles = [scifiles_1, scifiles_2]
scifiles_all = scifiles_1 + scifiles_2
nexp = len(scifiles_1)
# Make the new Science dir
# TODO: This needs to be defined by the user
scipath = os.path.join(pypeit_output_dir, 'Science')
if not os.path.isdir(scipath):
    msgs.info('Creating directory for Science output: {0}'.format(scipath))
    os.makedirs(scipath)

# TODO Should we flat field. The flat field and flat field error are wonky and probably nonsense
param_dict = {
    'extract_2d': {'save_results': True},
    'bkg_subtract': {'skip': True},
    'imprint_subtract': {'save_results': True},
    'msa_flagging': {'save_results': True},
    'master_background_mos': {'skip': True},
    'srctype': {'source_type': 'EXTENDED'},
    #    'flat_field': {'skip': True},
    'resample_spec': {'skip': True},
    'extract_1d': {'skip': True},
    'flat_field': {'save_interpolated_flat': True},  # Flats appear to be just nonsense. So skip for now.
}

basenames_1 = []
basenames_2 = []
for sci1, sci2 in zip(scifiles_1, scifiles_2):
    basenames_1.append(os.path.basename(sci1).replace('_rate.fits', ''))
    basenames_2.append(os.path.basename(sci2).replace('_rate.fits', ''))

# Run the spec2 pipeline
runflag = False
if runflag:
    for sci in scifiles_all:
        Spec2Pipeline.call(sci, save_results=True, output_dir=output_dir, steps=param_dict)
        #spec2 = Spec2Pipeline(steps=param_dict)
        #spec2.save_results = True
        #spec2.output_dir = output_dir
        #result = spec2(sci)

# Output file names
intflat_output_files_1 = []
#e2d_output_files_1 = []
cal_output_files_1 = []
msa_output_files_1 = []

intflat_output_files_2 = []
#e2d_output_files_2 = []
cal_output_files_2 = []
msa_output_files_2 = []

for base1, base2 in zip(basenames_1, basenames_2):
    #e2d_output_files_1.append(os.path.join(output_dir, base1 + '_extract_2d.fits'))
    intflat_output_files_1.append(os.path.join(output_dir, base1 + '_interpolatedflat.fits'))
    cal_output_files_1.append(os.path.join(output_dir, base1 + '_cal.fits'))
    msa_output_files_1.append(os.path.join(output_dir, base1 + '_msa_flagging.fits'))

    #e2d_output_files_2.append(os.path.join(output_dir, base2 + '_extract_2d.fits'))
    intflat_output_files_2.append(os.path.join(output_dir, base2 + '_interpolatedflat.fits'))
    cal_output_files_2.append(os.path.join(output_dir, base2 + '_cal.fits'))
    msa_output_files_2.append(os.path.join(output_dir, base2 + '_msa_flagging.fits'))

# Read in calwebb outputs for everytihng

# Read in multi exposure calwebb outputs
msa_multi_list_1 = []
#e2d_multi_list_1 = []
intflat_multi_list_1 = []
final_multi_list_1 = []
msa_multi_list_2 = []
#e2d_multi_list_2 = []
intflat_multi_list_2 = []
final_multi_list_2 = []
nslits_1 = np.zeros(nexp, dtype=int)
nslits_2 = np.zeros(nexp, dtype=int)
t_eff = np.zeros(nexp, dtype=float)

# TODO Figure out why this is so damn slow! I suspect it is calwebb
for iexp in range(nexp):
    # Open some JWST data models
    #e2d_multi_list_1.append(datamodels.open(e2d_output_files_1[iexp]))
    msa_multi_list_1.append(datamodels.open(msa_output_files_1[iexp]))
    intflat_multi_list_1.append(datamodels.open(intflat_output_files_1[iexp]))
    final_multi_list_1.append(datamodels.open(cal_output_files_1[iexp]))

    msa_multi_list_2.append(datamodels.open(msa_output_files_2[iexp]))
    intflat_multi_list_2.append(datamodels.open(intflat_output_files_2[iexp]))
    final_multi_list_2.append(datamodels.open(cal_output_files_2[iexp]))

    t_eff[iexp] = final_multi_list_1[iexp].meta.exposure.effective_exposure_time
    nslits_1[iexp] = len(final_multi_list_1[iexp].slits)
    nslits_2[iexp] = len(final_multi_list_2[iexp].slits)



# Some pypeit things
spectrograph = load_spectrograph('jwst_nirspec')
par = spectrograph.default_pypeit_par()
det_container_list = [spectrograph.get_detector_par(1), spectrograph.get_detector_par(2)]

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
if diff_redux:
    par['reduce']['findobj']['skip_skysub'] = True # Do not sky-subtract when object finding
    par['reduce']['extraction']['skip_optimal'] = True # Skip local_skysubtraction and profile fitting



# TODO Fix this, currently does not work if target names have - or _
filename_first = os.path.basename(scifiles_1[0])
filename_last = os.path.basename(scifiles_1[-1])
split_first = filename_first.split('_')
split_last = filename_last.split('_')
prefix_first = ''
for ii in range(len(split_first) - 3):
    prefix_first += split_first[ii] + "_"
diff_str = 'diff_' if diff_redux else ''
out_filename = diff_str + prefix_first + split_first[-3] + "-" + split_last[-3]

show = True

spec_samp_fact = 1.0
spat_samp_fact = 1.0

# Use the first exposure to se the slit names
slit_names_1 = [slit.name for slit in final_multi_list_1[0].slits]
slit_names_2 = [slit.name for slit in final_multi_list_2[0].slits]
slit_names_uni = np.unique(np.hstack([slit_names_1, slit_names_2]))

# Loop over slits
# islit = '10'
# islit = 'S200A1'
islit = '83'
#islit = None
gdslits = slit_names_uni[::-1] if islit is None else [islit]
bad_slits = []

# First index is detector, second index is exposure
#e2d_multi_list = [e2d_multi_list_1, e2d_multi_list_2]
msa_multi_list = [msa_multi_list_1, msa_multi_list_2]
intflat_multi_list = [intflat_multi_list_1, intflat_multi_list_2]
final_multi_list = [final_multi_list_1, final_multi_list_2]
slit_names_list = [slit_names_1, slit_names_2]
kludge_err = 1.5

# Is this correct?
bkg_indices = [1, 2, 1]

# TODO Are the tilts working correctly for the mosaics? I'm not so sure??

# Loop over all slits, create a list of spec2d objects and run 2d coadd
for ii, islit in enumerate(gdslits):
    spec2d_list = []
    offsets_pixels = []
    if show:
        display.clear_all()
    # Loop over detectors, and exposures
    for idet in range(len(detectors)):
        for iexp in range(nexp):
            indx = np.where(np.array(slit_names_list[idet]) == islit)[0]
            if len(indx) > 0:
                jj = indx[0]
                ibkg = bkg_indices[iexp]
                # Grab the calibration subimages for the first exposure, use for all subsequent exposures since
                # sometimes (because of a bug) the sizes of these calibraition (and science) sub-images can change.
                if iexp == 0:
                    slit_slice, slit_left, slit_righ, slit_left_orig, slit_righ_orig, spec_vals_orig, \
                    src_trace_ra, src_trace_dec, dq_sub, ra, dec, finitemask, waveimg, tilts, flatfield, pathloss, barshadow, \
                    photom_conversion, calwebb = jwst_extract_subimgs(
                        final_multi_list[idet][iexp].slits[jj], intflat_multi_list[idet][iexp].slits[jj])

                # Process the science image for this exposure
                science, sciivar, gpm, base_var, count_scale = jwst_proc(
                    msa_multi_list[idet][iexp], t_eff[iexp], slit_slice, finitemask, pathloss, barshadow,
                    noise_floor=par['scienceframe']['process']['noise_floor'],
                    kludge_err=kludge_err, ronoise=det_container_list[idet].ronoise)

                # If there are not good pixels continue
                if not np.any(gpm):
                    bad_slits.append(islit)
                    continue

                # Instantiate
                sciImg = PypeItImage(image=science, ivar=sciivar, bpm=np.logical_not(gpm).astype(int),
                                     base_var=base_var, img_scale=count_scale,
                                     rn2img=np.full_like(science, det_container_list[idet].ronoise[0]**2),
                                     detector=det_container_list[idet])

                if diff_redux:
                    # Process the background image using the same calibrations as the science
                    bkg, bkgivar, bkg_gpm, bkg_base_var, bkg_count_scale = jwst_proc(
                        msa_multi_list[idet][ibkg], t_eff[ibkg], slit_slice, finitemask, pathloss, barshadow,
                        noise_floor=par['scienceframe']['process']['noise_floor'],
                        kludge_err=kludge_err, ronoise=det_container_list[idet].ronoise)

                    # If there are not good pixels continue
                    if not np.any(bkg_gpm):
                        bad_slits.append(islit)
                        continue

                    # Instantiate
                    bkgImg = PypeItImage(image=bkg, ivar=bkgivar, bpm=np.logical_not(bkg_gpm).astype(int),
                                         base_var=bkg_base_var, img_scale=bkg_count_scale,
                                         rn2img=np.full_like(bkg, det_container_list[idet].ronoise[0]**2),
                                         detector=det_container_list[idet])

                    # Perform the difference imaging, propagate the error and masking
                    sciImg = sciImg.sub(bkgImg, par['scienceframe']['process'])

                nspec, nspat = waveimg.shape
                slits = slittrace.SlitTraceSet(slit_left, slit_righ, pypeline, detname=det_container_list[idet].name,
                                               nspat=nspat,
                                               PYP_SPEC=spectrograph.name)
                slits.maskdef_id = np.array([ii])


                # Construct the Spec2DObj with the positive image
                spec2DObj = spec2dobj.Spec2DObj(sciimg=sciImg.image,
                                                ivarraw=sciImg.ivar,
                                                skymodel=np.zeros_like(sciImg.image),
                                                objmodel=np.zeros_like(sciImg.image),
                                                ivarmodel=sciImg.ivar,
                                                scaleimg=None,
                                                waveimg=waveimg,
                                                bpmmask=sciImg.select_flag().astype(int),
                                                detector=sciImg.detector,
                                                sci_spat_flexure=None,
                                                sci_spec_flexure=None,
                                                vel_corr=None,
                                                vel_type=None,
                                                tilts=tilts,
                                                slits=slits,
                                                wavesol=None,
                                                maskdef_designtab=None)

                spec2d_list.append(spec2DObj)
                offsets_pixels.append(offsets_pixels_list[idet][iexp])

                if show and (iexp == 0):
                    sci_rate = datamodels.open(msa_multi_list[idet][iexp])
                    sci_data = np.array(sci_rate.data.T, dtype=float)

                    bkg_rate = datamodels.open(msa_multi_list[idet][ibkg])
                    bkg_data = np.array(bkg_rate.data.T, dtype=float)

                    viewer_raw, ch_raw = display.show_image(sci_data, cuts=get_cuts(sci_data),
                                                            chname='raw_rate_iexp_{:d}_{:s}'.format(iexp, detectors[idet]))
                    viewer_sci, ch_sci = display.show_image(sci_data-bkg_data, cuts=get_cuts(sci_data-bkg_data),
                                                            chname='diff_rate_iexp_{:d}_{:s}'.format(iexp, detectors[idet]))
                    display.show_slits(viewer_sci, ch_sci, slit_left_orig, slit_righ_orig, spec_vals=spec_vals_orig,
                                       pstep=1, slit_ids=np.array([islit]))
                    display.show_slits(viewer_raw, ch_raw, slit_left_orig, slit_righ_orig, spec_vals=spec_vals_orig,
                                       pstep=1, slit_ids=np.array([islit]))
                    # Show individual calibrations
                    if idet == 0:
                        viewer_data, ch_data = display.show_image(calwebb, waveimg=waveimg, cuts=get_cuts(calwebb),
                                                                  chname='calwebb_{:s}'.format(detectors[idet]))
                        display.show_trace(viewer_data, ch_data, src_trace_ra, 'trace-RA_{:s}'.format(detectors[idet]), color='#f0e442', pstep=1)
                        display.show_trace(viewer_data, ch_data, src_trace_dec, 'trace-DEC_{:s}'.format(detectors[idet]), color='#f0e442', pstep=1)
                        viewer_pypeit, ch_pypeit = display.show_image(sciImg.image, waveimg=waveimg, cuts=get_cuts(sciImg.image),
                                                                      chname='pypeit_{:s}'.format(detectors[idet]))
                        display.show_slits(viewer_pypeit, ch_pypeit, slit_left, slit_righ, pstep=1, slit_ids=np.array([islit]))
                        viewer_wave, ch_wave = display.show_image(waveimg, chname='wave_{:s}'.format(detectors[idet]))
                        viewer_tilts, ch_tilts = display.show_image(tilts, waveimg=waveimg, chname='tilts_{:s}'.format(detectors[idet]))
                        viewer_flat, ch_flat = display.show_image(flatfield, waveimg=waveimg, chname='flat_{:s}'.format(detectors[idet]))
                        viewer_path, ch_path = display.show_image(pathloss, waveimg=waveimg, chname='pathloss_{:s}'.format(detectors[idet]))
                        viewer_bar, ch_bar = display.show_image(barshadow, waveimg=waveimg, chname='barshadow_{:s}'.format(detectors[idet]),
                                                                cuts=(0.0, 1.0))


            else:
                continue

    if len(spec2d_list) > 0:
        basename = '{:s}_{:s}'.format(out_filename, 'slit' + islit)

        # TODO implement an option to extract everything onto the same wavelength grid optimized to match the
        # coverage of JWST. For example for the prism things are quite nonlinear

        # TODO Not sure what to do with the detector container here
        # Instantiate Coadd2d
        coAdd = coadd2d.CoAdd2D.get_instance(spec2d_list, spectrograph, par, det=det_container_list[0].det,
                                             offsets=offsets_pixels, weights='uniform',
                                             spec_samp_fact=spec_samp_fact,
                                             spat_samp_fact=spat_samp_fact,
                                             bkg_redux=False, debug=show)

        coadd_dict_list = coAdd.coadd(only_slits=None, interp_dspat=True)

        # Create the pseudo images
        pseudo_dict = coAdd.create_pseudo_image(coadd_dict_list)

        sciimg_coadd, sciivar_coadd, skymodel_coadd, objmodel_coadd, ivarmodel_coadd, \
        outmask_coadd, sobjs_coadd, detector_coadd, slits_coadd, tilts_coadd, waveimg_coadd = coAdd.reduce(
            pseudo_dict, show=True, clear_ginga=False, show_peaks=show,
            basename=basename)

        # Tack on detector (similarly to pypeit.extract_one)
        for sobj in sobjs_coadd:
            sobj.DETECTOR = det_container_list[0]

        # TODO not currently using counts_scale and base_var. Need to rework coadd2d to operate on sciimgs

        # Construct the Spec2DObj with the positive image
        spec2DObj_coadd = spec2dobj.Spec2DObj(sciimg=sciimg_coadd,
                                              ivarraw=sciivar_coadd,
                                              skymodel=skymodel_coadd,
                                              objmodel=objmodel_coadd,
                                              ivarmodel=ivarmodel_coadd,
                                              scaleimg=None,
                                              waveimg=waveimg_coadd,
                                              bpmmask=outmask_coadd,
                                              detector=det_container_list[0],
                                              sci_spat_flexure=None,
                                              sci_spec_flexure=None,
                                              vel_corr=None,
                                              vel_type=None,
                                              tilts=tilts_coadd,
                                              slits=slits_coadd,
                                              wavesol=None,
                                              maskdef_designtab=None)

        # QA
        if show:
            # Plot the 2d
            nobj = len(sobjs_coadd)
            left, right, mask = spec2DObj_coadd.slits.select_edges()
            slid_IDs = spec2DObj_coadd.slits.slitord_id
            image = spec2DObj_coadd.sciimg  # Processed science image
            mean, med, sigma = sigma_clipped_stats(image[spec2DObj_coadd.bpmmask == 0], sigma_lower=5.0, sigma_upper=5.0)
            cut_min = mean - 4.0 * sigma
            cut_max = mean + 4.0 * sigma
            chname_sci = 'slit-science-' + islit
            # Clear all channels at the beginning
            viewer, ch_sci = display.show_image(image, chname=chname_sci,
                                                waveimg=spec2DObj_coadd.waveimg,
                                                clear=False,
                                                cuts=(cut_min, cut_max))

            if nobj > 0:
                show_trace(sobjs_coadd, 'DET01', viewer, ch_sci)
            display.show_slits(viewer, ch_sci, left, right, slit_ids=slid_IDs)

            image = np.sqrt(inverse(spec2DObj_coadd.ivarmodel))  # Processed science image
            mean, med, sigma = sigma_clipped_stats(image[spec2DObj_coadd.bpmmask == 0], sigma_lower=5.0, sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            chname_sig = 'slit-sigma-' + islit
            # Clear all channels at the beginning
            viewer, ch_sig = display.show_image(image, chname=chname_sig,
                                                waveimg=spec2DObj_coadd.waveimg,
                                                clear=False,
                                                cuts=(cut_min, cut_max))

            if nobj > 0:
                show_trace(sobjs_coadd, 'DET01', viewer, ch_sig)
            display.show_slits(viewer, ch_sig, left, right, slit_ids=slid_IDs)

            channel_names = [chname_sci, chname_sig]

            # After displaying all the images sync up the images with WCS_MATCH
            shell = viewer.shell()
            shell.start_global_plugin('WCSMatch')
            shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [channel_names[-1]],
                                            {})
            if not diff_redux:
                slitmask_coadd = slits_coadd.slit_img(initial=False, flexure=None, exclude_flag=None)
                slitord_id = slits_coadd.slitord_id[0]
                thismask = slitmask_coadd == slitord_id

                gpm_extract = spec2DObj_coadd.bpmmask == 0
                # Make a plot of the residuals for a random slit
                chi = (spec2DObj_coadd.sciimg - spec2DObj_coadd.objmodel - spec2DObj_coadd.skymodel) * np.sqrt(
                    spec2DObj_coadd.ivarmodel) * gpm_extract

                maskchi = thismask & gpm_extract

                n_bins = 50
                sig_range = 7.0
                binsize = 2.0 * sig_range / n_bins
                bins_histo = -sig_range + np.arange(n_bins) * binsize + binsize / 2.0

                xvals = np.arange(-10.0, 10, 0.02)
                gauss = scipy.stats.norm(loc=0.0, scale=1.0)
                gauss_corr = scipy.stats.norm(loc=0.0, scale=1.0)

                sigma_corr, maskchi = coadd.renormalize_errors(chi, maskchi, max_corr=20.0, title='jwst_sigma_corr',
                                                               debug=True)


            if nobj > 0:
                wv_gpm = sobjs_coadd[0].BOX_WAVE > 1.0
                plt.plot(sobjs_coadd[0].BOX_WAVE[wv_gpm], sobjs_coadd[0].BOX_COUNTS[wv_gpm]*sobjs_coadd[0].BOX_MASK[wv_gpm],
                color='black', drawstyle='steps-mid', label='Counts')
                plt.plot(sobjs_coadd[0].BOX_WAVE[wv_gpm], sobjs_coadd[0].BOX_COUNTS_SIG[wv_gpm]*sobjs_coadd[0].BOX_MASK[wv_gpm],
                color='red', drawstyle='steps-mid', label='Counts Error')
                plt.legend()
                plt.show()




        # container for specobjs and Spec2d
        all_specobjs = specobjs.SpecObjs()

        all_spec2d = spec2dobj.AllSpec2DObj()
        # set some meta
        all_spec2d['meta']['bkg_redux'] = False
        all_spec2d['meta']['find_negative'] = False

        # fill the specobjs container
        all_specobjs.add_sobj(sobjs_coadd)
        all_spec2d[det_container_list[0].name] = spec2DObj_coadd

        # THE FOLLOWING MIMICS THE CODE IN pypeit.save_exposure()
        scipath = os.path.join(pypeit_output_dir, 'Science')
        if not os.path.isdir(scipath):
            msgs.info('Creating directory for Science output: {0}'.format(scipath))

        # Write out specobjs
        # Build header for spec2d
        head2d = fits.getheader(cal_output_files_1[0])
        subheader = spectrograph.subheader_for_spec(head2d, head2d, allow_missing=True)
        if all_specobjs.nobj > 0:
            outfile1d = os.path.join(scipath, 'spec1d_{:s}.fits'.format(basename))
            all_specobjs.write_to_fits(subheader, outfile1d)

            # Info
            outfiletxt = os.path.join(scipath, 'spec1d_{:s}.txt'.format(basename))
            all_specobjs.write_info(outfiletxt, spectrograph.pypeline)

        # Build header for spec2d
        outfile2d = os.path.join(scipath, 'spec2d_{:s}.fits'.format(basename))
        # TODO For the moment hack so that we can write this out
        pri_hdr = all_spec2d.build_primary_hdr(head2d, spectrograph, subheader=subheader,
                                               redux_path=None, master_key_dict=None, master_dir=None)
        # Write spec2d
        all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr, overwrite=True)



