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
from jwst_utils import compute_diff, get_cuts, jwst_show_spec2, jwst_show_msa, jwst_proc, jwst_extract_subimgs, jwst_get_slits
from jwst_utils import NIRSpecSlitCalibrations, jwst_mosaic, jwst_reduce
from pypeit.metadata import PypeItMetaData
from pypeit.display import display
from pypeit import specobjs
from pypeit import slittrace
from pypeit.utils import inverse, fast_running_median
from pypeit.core import findobj_skymask
from pypeit.core import skysub, coadd
from pypeit.core import procimg
from pypeit.core import flat

from pypeit.images.mosaic import Mosaic


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

# JFH: This is the main up to date routine. Ignore the others.


DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']

# detname = 'nrs1'
# detector = 1 if 'nrs1' in detname else 2

#disperser = 'J0313_G235M'
#disperser = 'G395M_Maseda'
#disperser = 'G395M'
#disperser = 'PRISM_01117'
# disperser = 'G235M'
#disperser='PRISM_01133'
disperser = 'PRISM_02756'
# detectors = ['nrs1', 'nrs2']
# disperser='PRISM_01117'
# disperser='PRISM_FS'

detectors = ['nrs1', 'nrs2']
exp_list = []

bkg_redux = True
runflag = False
mode = 'MSA'
#mode ='FS'
#slit = 'S200A1'
#islit = 'S200A2'
#islit = '37'
#reduce_sources = ['2756_10025']
reduce_sources = ['2756_202']
reduce_slits = None

# If bkg_redux is False, the code will model the sky and the object profile and perform optimal extraction.
# If bkg_redux is True, the code will difference image and simply boxcar extract (optimal not implemented yet)
for detname in detectors:
    # TODO add the kendrew FS SN data to this.
    if 'PRISM_02756' == disperser:
        ## Prorgram for Slit Loss Characterization for MSA shutters
        # PRISM data
        rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/NIRSPEC_2756/level_12/02756/'
        redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/02756_CLEAR_PRISM/calwebb'
        output_dir = os.path.join(redux_dir, 'output')
        pypeit_output_dir = os.path.join(redux_dir, 'pypeit')

        # NIRSPEC 3-point dither
        # dither center
        scifile1 = os.path.join(rawpath_level2, 'jw02756001001_03101_00001_' + detname + '_rate.fits')
        scifile2 = os.path.join(rawpath_level2, 'jw02756001001_03101_00002_' + detname + '_rate.fits')
        scifile3 = os.path.join(rawpath_level2, 'jw02756001001_03101_00003_' + detname + '_rate.fits')

        # dither offset
        # scifile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_rate.fits')
        # bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_rate.fits')
        # bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_rate.fits')

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
    elif 'PRISM_FS' == disperser:
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

    elif 'J0313_G235M' == disperser:
        ## Prorgram for Slit Loss Characterization for MSA shutters
        # PRISM data
        rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1764/level_12/01764/'
        output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J0313_G235M/calwebb/output'
        pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J0313_G235M/calwebb/pypeit'

        # NIRSPEC 3-point dither
        # dither center

        if reduce_slits[0] == 'S200A1':
            scifile1 = os.path.join(rawpath_level2, 'jw01764014001_03102_00001_' + detname + '_rate.fits')
            scifile2 = os.path.join(rawpath_level2, 'jw01764014001_03102_00002_' + detname + '_rate.fits')
            scifile3 = os.path.join(rawpath_level2, 'jw01764014001_03102_00003_' + detname + '_rate.fits')
        elif reduce_slits[0] == 'S200A2':
            scifile1 = os.path.join(rawpath_level2, 'jw01764014001_03104_00001_' + detname + '_rate.fits')
            scifile2 = os.path.join(rawpath_level2, 'jw01764014001_03104_00002_' + detname + '_rate.fits')
            scifile3 = os.path.join(rawpath_level2, 'jw01764014001_03104_00003_' + detname + '_rate.fits')


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

#if 'MSA' in mode:
#    offsets_pixels_list = [[0, 5.0, -5.0], [0, 5.0, -5.0]]
#elif 'FS' in mode:
#    offsets_pixels_list = [[0, 8.0, 18.0], [0, 8.0, 18.0]]


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


# Some pypeit things
spectrograph = load_spectrograph('jwst_nirspec')
par = spectrograph.default_pypeit_par()
det_container_list = [spectrograph.get_detector_par(1), spectrograph.get_detector_par(2)]
fitstbl_1 = PypeItMetaData(spectrograph, par=par,files=scifiles_1, strict=True)
fitstbl_2 = PypeItMetaData(spectrograph, par=par,files=scifiles_2, strict=True)
fitstbls = [fitstbl_1, fitstbl_2]

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
if bkg_redux:
    par['reduce']['findobj']['skip_skysub'] = True # Do not sky-subtract when object finding
    par['reduce']['extraction']['skip_optimal'] = True # Skip local_skysubtraction and profile fitting


# TODO Should we flat field. The flat field and flat field error are wonky and probably nonsense
param_dict = {
    'extract_2d': {'save_results': True},
    'bkg_subtract': {'skip': True},
    'imprint_subtract': {'save_results': True}, # TODO Check up on whether imprint subtraction is being done by us???
    'master_background_mos': {'skip': True},
    'srctype': {'source_type': 'EXTENDED'},
    #    'flat_field': {'skip': True},
    'resample_spec': {'skip': True},
    'extract_1d': {'skip': True},
    'flat_field': {'save_interpolated_flat': True},  # Flats appear to be just nonsense. So skip for now.
}

# For MSA data, we need to run the MSA flagging step and this generates the full 2d frames we operate on, for
# FS data, the MSA flagging is not performed and so the last step is the assign_wcs step
if mode =='MSA':
    param_dict['msa_flagging'] = {'save_results': True}
elif mode == 'FS':
    param_dict['assign_wcs'] = {'save_results': True}

basenames = []
basenames_1 = []
basenames_2 = []
for sci1, sci2 in zip(scifiles_1, scifiles_2):
    b1 = os.path.basename(sci1).replace('_rate.fits', '')
    b2 = os.path.basename(sci2).replace('_rate.fits', '')
    basenames_1.append(b1)
    basenames_2.append(b2)
    basenames.append(b1.replace('_nrs1', ''))

# Run the spec2 pipeline
if runflag:
    for sci in scifiles_all:
        Spec2Pipeline.call(sci, save_results=True, output_dir=output_dir, steps=param_dict)
        #spec2 = Spec2Pipeline(steps=param_dict)
        #spec2.save_results = True
        #spec2.output_dir = output_dir
        #result = spec2(sci)

# Output file names
intflat_output_files_1 = []
msa_output_files_1 = []
cal_output_files_1 = []

intflat_output_files_2 = []
msa_output_files_2 = []
cal_output_files_2 = []

for base1, base2 in zip(basenames_1, basenames_2):
    if mode == 'MSA':
        msa_output_files_1.append(os.path.join(output_dir, base1 + '_msa_flagging.fits'))
        msa_output_files_2.append(os.path.join(output_dir, base2 + '_msa_flagging.fits'))
    else:
        msa_output_files_1.append(os.path.join(output_dir, base1 + '_assign_wcs.fits'))
        msa_output_files_2.append(os.path.join(output_dir, base2 + '_assign_wcs.fits'))

    intflat_output_files_1.append(os.path.join(output_dir, base1 + '_interpolatedflat.fits'))
    cal_output_files_1.append(os.path.join(output_dir, base1 + '_cal.fits'))
    intflat_output_files_2.append(os.path.join(output_dir, base2 + '_interpolatedflat.fits'))
    cal_output_files_2.append(os.path.join(output_dir, base2 + '_cal.fits'))

# Read in calwebb outputs for everytihng

# Read in multi exposure calwebb outputs
msa_multi_list_1 = []
intflat_multi_list_1 = []
final_multi_list_1 = []
msa_multi_list_2 = []
intflat_multi_list_2 = []
final_multi_list_2 = []
#nslits_1 = np.zeros(nexp, dtype=int)
#nslits_2 = np.zeros(nexp, dtype=int)
#t_eff = np.zeros(nexp, dtype=float)

ndetectors = 2
# Create arrays to hold JWST spec2, but only load the files when they're needed
msa_data = np.empty((ndetectors, nexp), dtype=object)
flat_data = np.empty((ndetectors, nexp), dtype=object)
cal_data = np.empty((ndetectors, nexp), dtype=object)

dither_offsets = np.zeros((ndetectors,nexp), dtype=float)
# TODO: This probably isn't correct.  I.e., need to know offsets and slit
# position angle.
for iexp in range(nexp):
    with fits.open(scifiles_1[iexp]) as hdu:
        dither_offsets[0,iexp] = hdu[0].header['YOFFSET']
for idet in range(1,ndetectors):
    dither_offsets[idet] = dither_offsets[0]
dither_offsets_pixels = dither_offsets.copy()
for idet in range(ndetectors):
    dither_offsets_pixels[idet] /= det_container_list[idet].platescale
# NOTE: Sign convention requires this calculation of the offset
dither_offsets_pixels = dither_offsets_pixels[:,0,None] - dither_offsets_pixels
print(dither_offsets_pixels)



# TODO Figure out why this is so damn slow! I suspect it is calwebb1
for iexp in range(nexp):
    # Open some JWST data models
    #e2d_multi_list_1.append(datamodels.open(e2d_output_files_1[iexp]))
    msa_data[0, iexp] = datamodels.open(msa_output_files_1[iexp])
    flat_data[0, iexp] = datamodels.open(intflat_output_files_1[iexp])
    cal_data[0, iexp] = datamodels.open(cal_output_files_1[iexp])

    msa_data[1, iexp] = datamodels.open(msa_output_files_2[iexp])
    flat_data[1, iexp] = datamodels.open(intflat_output_files_2[iexp])
    cal_data[1, iexp] = datamodels.open(cal_output_files_2[iexp])


show = True

# Create a set of aligned slit and source names for both detectors


# Use the first exposure to se the slit names
# (ndet, nslit)
slit_names_1 = [slit.name for slit in cal_data[0,0].slits]
slit_names_2 = [slit.name for slit in cal_data[1,0].slits]
slit_names_tot = np.hstack([slit_names_1, slit_names_2])
source_names_1 = [slit.source_name for slit in cal_data[0,0].slits]
source_names_2 = [slit.source_name for slit in cal_data[1,0].slits]
source_names_tot = np.hstack([source_names_1, source_names_2])

# Find the unique slit names and the unique sources aligned with those slits
slit_names_uni, uni_indx = np.unique(slit_names_tot, return_index=True)
source_names_uni = source_names_tot[uni_indx]
slit_sources_uni = [(slit, source) for slit, source in zip(slit_names_uni, source_names_uni)]


# Loop over slits
#islit = '10'
#islit = 'S200A1'
#islit = '83'
#islit = None
#islit = '63'

#bad_slits = []
#gdsources = source_names_uni[::-1] if source is None else [source]



# First index is detector, second index is exposure
#msa_multi_list = [msa_multi_list_1, msa_multi_list_2]
#msa_multi_list = [msa_multi_list_1, msa_multi_list_2]
#intflat_multi_list = [intflat_multi_list_1, intflat_multi_list_2]
#final_multi_list = [final_multi_list_1, final_multi_list_2]
#slit_names_list = [slit_names_1, slit_names_2]
kludge_err = 1.5

# Is this correct?
bkg_indices = [1, 2, 1]
#detector_gap = 180 # 18" divided by 0.1" pixels from the JWST website


iexp_ref = 0
if not os.path.isdir(scipath):
    msgs.info('Creating directory for Science output: {0}'.format(scipath))

# TODO Fix this, currently does not work if target names have - or _
out_filenames = basenames
#diff_str = 'diff_' if bkg_redux else ''
#out_filenames = [diff_str + base for base in basenames]


if reduce_slits is not None:
    gd_slits_sources = [(slt, src) for slt, src in slit_sources_uni for slit in reduce_slits if slt == slit]
elif reduce_sources is not None:
    gd_slits_sources = [(slt, src) for slt, src in slit_sources_uni for source in reduce_sources if src == source]
else:
    gd_slits_sources = slit_sources_uni

# Loop over all exposures, loop over all slits, create a mosaic and save them to individual PypeIt spec2d files.
for iexp in range(nexp):
    for ii, (islit, isource) in enumerate(gd_slits_sources):
        # Container for all the Spec2DObj, different spec2dobj and specobjs for each slit
        all_spec2d = spec2dobj.AllSpec2DObj()
        all_spec2d['meta']['bkg_redux'] = bkg_redux
        all_spec2d['meta']['find_negative'] = bkg_redux
        # Container for the specobjs
        all_specobjs = specobjs.SpecObjs()

        # TODO this step is being executed repeatedly for each new exposure?
        # Generate the calibrations from the reference exposure
        # TODO This step is only performed with a reference exposure because calwebb has an annoying property that
        # it does not always extract the same subimage spectral pixels for the different dithers in the dither pattern.
        # This seems to be a bug in calwebb, since it is unclear why the subimage calibrations should change.
        # This is problem for 2d coadding, since then the offsets in the detector frame will be not allow one to register
        # the frames. It is possible to fix this by using the RA/DEC images provided by calwebb to determine the
        # actual locations on the sky, which would be preferable. However, this does not appear to be working correctly
        # in calwebb. So for now, we just use the first exposure as the reference exposure for the calibrations.
        CalibrationsNRS1 = NIRSpecSlitCalibrations(
            det_container_list[0], cal_data[0, iexp_ref], flat_data[0, iexp_ref], islit)
        CalibrationsNRS2 = NIRSpecSlitCalibrations(
            det_container_list[1], cal_data[1, iexp_ref], flat_data[1, iexp_ref], islit)

        ibkg = bkg_indices[iexp]
        # Create the image mosaic
        sciImg, slits, waveimg, tilts, ndet = jwst_mosaic(msa_data[:, iexp], [CalibrationsNRS1, CalibrationsNRS2], kludge_err=kludge_err,
            noise_floor=par['scienceframe']['process']['noise_floor'], bkg_image_model_tuple=msa_data[:, ibkg],
            show=show & (iexp == iexp_ref))

        # If this is a bkg_redux, perform background subtraction
        if bkg_redux:
            bkgImg, _, _, _, _= jwst_mosaic(msa_data[:, ibkg], [CalibrationsNRS1, CalibrationsNRS2], kludge_err=kludge_err,
                                          noise_floor=par['scienceframe']['process']['noise_floor'])
            sciImg = sciImg.sub(bkgImg)


        # Run the reduction
        all_spec2d[sciImg.detector.name], tmp_sobjs = jwst_reduce(sciImg, slits, waveimg, tilts, spectrograph, par,
                                                     show=show, find_negative=bkg_redux, bkg_redux=bkg_redux,
                                                     clear_ginga=False, show_peaks=show, show_skysub_fit=True,
                                                     basename=basenames[iexp])
        # Hold em
        if tmp_sobjs.nobj > 0:
            all_specobjs.add_sobj(tmp_sobjs)

            if show:
                # Plot boxcar
                wv_gpm_box = tmp_sobjs[0].BOX_WAVE > 1.0
                plt.plot(tmp_sobjs[0].BOX_WAVE[wv_gpm_box],tmp_sobjs[0].BOX_COUNTS[wv_gpm_box] * tmp_sobjs[0].BOX_MASK[wv_gpm_box],
                color='green', drawstyle='steps-mid', label='Boxcar Counts')
                plt.plot(tmp_sobjs[0].BOX_WAVE[wv_gpm_box], tmp_sobjs[0].BOX_COUNTS_SIG[wv_gpm_box] * tmp_sobjs[0].BOX_MASK[wv_gpm_box],
                     color='cyan', drawstyle='steps-mid', label='Boxcar Counts Error')

                # plot optimal
                if tmp_sobjs[0].OPT_WAVE is not None:
                    wv_gpm_opt = tmp_sobjs[0].OPT_WAVE > 1.0
                    plt.plot(tmp_sobjs[0].OPT_WAVE[wv_gpm_opt], tmp_sobjs[0].OPT_COUNTS[wv_gpm_opt] * tmp_sobjs[0].OPT_MASK[wv_gpm_opt],
                             color='black', drawstyle='steps-mid', label='Optimal Counts')
                    plt.plot(tmp_sobjs[0].OPT_WAVE[wv_gpm_opt], tmp_sobjs[0].OPT_COUNTS_SIG[wv_gpm_opt] * tmp_sobjs[0].OPT_MASK[wv_gpm_opt],
                             color='red', drawstyle='steps-mid', label='Optimal Counts Error')

                plt.legend()
                plt.show()

        # THE FOLLOWING MIMICS THE CODE IN pypeit.save_exposure()
        basename = '{:s}_{:s}'.format(out_filenames[iexp], 'source_' + isource)

        # TODO Populate the header with metadata relevant to this source?

        # Write out specobjs
        # Build header for spec2d
        head2d = fits.getheader(scifiles_1[iexp])
        subheader = spectrograph.subheader_for_spec(fitstbl_1[iexp], head2d, allow_missing=False)
        # Overload the target name with the source name
        subheader['target'] = isource
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
                                               redux_path=None, calib_dir=None)
        # Write spec2d
        all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr, overwrite=True)

