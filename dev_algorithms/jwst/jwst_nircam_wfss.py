

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
runflag = True
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

flat_output_files = []
for base in basenames:
    flat_output_files.append(os.path.join(output_dir, base + '_flatfieldstep.fits'))

rate_obj = datamodels.open(flat_output_files[0])

if show:
    display.show_image(rate_obj.data, chname='raw rate')

