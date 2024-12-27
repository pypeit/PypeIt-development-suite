
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
DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']


detname = 'nrs1'
detector = 1 if 'nrs1' in detname else 2
disperser = 'G395M'
#disperser = 'G235M'
#disperser='PRISM_01133'
#disperser='PRISM_01117'
if 'PRISM_01133' in disperser:
    # PRISM data
    rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/Raw'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/output'
    pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/pypeit'


    # NIRSPEC 3-point dither
    # dither center
    scifile1  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_rate.fits')
    scifile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_rate.fits')
    scifile3 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_rate.fits')

    # dither offset
    #scifile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_rate.fits')
    #bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_rate.fits')
    #bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_rate.fits')
elif 'PRISM_01117' in disperser:
    # PRISM data
    rawpath_level2 = '//Users/joe/jwst_redux/Raw/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/level_12/01117'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb/output'
    pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb/pypeit'


    # NIRSPEC 3-point dither
    # dither center
    scifile1  = os.path.join(rawpath_level2, 'jw01117007001_03101_00002_' + detname + '_rate.fits')
    scifile2 = os.path.join(rawpath_level2, 'jw01117007001_03101_00003_' + detname + '_rate.fits')
    scifile3 = os.path.join(rawpath_level2, 'jw01117007001_03101_00004_' + detname + '_rate.fits')

elif 'G395M' in disperser:
    # Use islit = 37 for nrs1
    # G395M data
    rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/Raw'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/output'
    pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/pypeit'

    # NIRSPEC 3-point dither
    scifile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + detname + '_rate.fits')
    scifile2 = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + detname + '_rate.fits')
    scifile3 = os.path.join(rawpath_level2, 'jw02736007001_03103_00003_' + detname + '_rate.fits')
elif 'G235M' in disperser:
    # Use islit = 38 for nrs1
    # G235M data
    rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/level_2/'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/output'
    pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/pypeit/'

    # NIRSPEC 3-point dither
    scifile1  = os.path.join(rawpath_level2, 'jw02736007001_03101_00002_' + detname + '_rate.fits')
    scifile2 = os.path.join(rawpath_level2, 'jw02736007001_03101_00003_' + detname + '_rate.fits')
    scifile3 = os.path.join(rawpath_level2, 'jw02736007001_03101_00004_' + detname + '_rate.fits')


# Make the new Science dir
# TODO: This needs to be defined by the user
scipath = os.path.join(pypeit_output_dir, 'Science')
if not os.path.isdir(scipath):
    msgs.info('Creating directory for Science output: {0}'.format(scipath))
    os.makedirs(scipath)

scifiles = [scifile1, scifile2, scifile3]

# TODO Should we flat field. The flat field and flat field error are wonky and probably nonsense
param_dict = {
    'extract_2d': {'save_results': True},
    'bkg_subtract': {'skip': True},
    'imprint_subtract': {'save_results': True},
    'msa_flagging': {'save_results': True},
    'master_background_mos': {'skip': True},
    'srctype': {'source_type':'EXTENDED'},
#    'flat_field': {'skip': True},
    'resample_spec': {'skip': True},
    'extract_1d': {'skip': True},
    'flat_field': {'save_interpolated_flat': True}, # Flats appear to be just nonsense. So skip for now.
}



basenames = []
for sci in scifiles:
    basenames.append(os.path.basename(sci).replace('_rate.fits', ''))

# Run the spec2 pipeline
runflag = False
if runflag:
    for sci in scifiles:
        spec2 = Spec2Pipeline(steps=param_dict)
        spec2.save_results = True
        spec2.output_dir = output_dir
        result = spec2(sci)

# Output file names
intflat_output_files = []
e2d_output_files =[]
cal_output_files =[]

for ii, basename in enumerate(basenames):
    e2d_output_files.append(os.path.join(output_dir, basename + '_extract_2d.fits'))
    intflat_output_files.append(os.path.join(output_dir, basename + '_interpolatedflat.fits'))
    cal_output_files.append(os.path.join(output_dir, basename + '_cal.fits'))


# Some pypeit things
spectrograph = load_spectrograph('jwst_nirspec')
par = spectrograph.default_pypeit_par()
det_container = spectrograph.get_detector_par(detector)
pypeline = 'MultiSlit'
par['rdx']['redux_path'] = pypeit_output_dir
qa_dir = os.path.join(pypeit_output_dir, 'QA')
par['rdx']['qadir'] = 'QA'
png_dir = os.path.join(qa_dir,'PNGs')
if not os.path.isdir(qa_dir):
    msgs.info('Creating directory for QA output: {0}'.format(qa_dir))
    os.makedirs(qa_dir)
if not os.path.isdir(png_dir):
    os.makedirs(png_dir)

# TODO Fix this, currently does not work if target names have - or _
filename_first = os.path.basename(scifiles[0])
filename_last = os.path.basename(scifiles[-1])
split_first = filename_first.split('_')
split_last = filename_last.split('_')
prefix_first = ''
for ii in range(len(split_first)-3):
    prefix_first += split_first[ii] + "_"
out_filename = prefix_first + split_first[-3] + "-" + split_last[-3]

nexp = len(scifiles)
nslits = np.zeros(nexp, dtype=int)
t_eff = np.zeros(nexp, dtype=float)

#islit = 37
show=True

offsets_pixels = [0, 5.0, -5.0]
spec_samp_fact = 1.0
spat_samp_fact = 1.0

# Read in multi exposure calwebb outputs
e2d_multi_list = []
intflat_multi_list = []
final_multi_list = []
for iexp in range(nexp):
    # Open some JWST data models
    e2d_multi_list.append(datamodels.open(e2d_output_files[iexp]))
    intflat_multi_list.append(datamodels.open(intflat_output_files[iexp]))
    final_multi_list.append(datamodels.open(cal_output_files[iexp]))

    t_eff[iexp] = e2d_multi_list[iexp].meta.exposure.effective_exposure_time
    nslits[iexp] = len(final_multi_list[iexp].slits)

islit = None
gdslits = np.arange(nslits[0]) if islit is None else [islit]
bad_slits = []

islit = 8

ra_list = []
waveimg_list = []
slit_left_list = []
slit_righ_list = []
slit_left_orig_list = []
slit_righ_orig_list = []

for iexp in range(nexp):
    slit_left, slit_righ, slit_left_orig, slit_righ_orig, spec_vals_orig, src_trace_ra, src_trace_dec, \
    rate, rate_var_rnoise, rate_var_poisson, rate_var_tot, dq, ra, dec, waveimg, tilts, flatfield, pathloss, \
    barshadow, photom_conversion, final = jwst_extract_subimgs(
        e2d_multi_list[iexp].slits[islit], final_multi_list[iexp].slits[islit], intflat_multi_list[iexp].slits[islit])
    ra_list.append(ra)
    waveimg_list.append(waveimg)
    slit_left_list.append(slit_left)
    slit_righ_list.append(slit_righ)
    slit_left_orig_list.append(slit_left_orig)
    slit_righ_orig_list.append(slit_righ_orig)

