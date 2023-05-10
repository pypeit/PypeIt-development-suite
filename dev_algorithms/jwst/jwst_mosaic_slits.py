import os
import numpy as np
import scipy
from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.stats import sigma_clipped_stats

from IPython import embed

# set environment variables
os.environ['CRDS_PATH'] = '/Users/suksientie/crds_cache/' #'/Users/joe/crds_cache/jwst_pub'
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu' #'https://jwst-crds-pub.stsci.edu'
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

def get_one_calib_obj(calfile, flatfile):

    spectrograph = load_spectrograph('jwst_nirspec')
    par = spectrograph.default_pypeit_par()

    if 'nrs1' in os.path.basename(calfile):
        det_container_list = [spectrograph.get_detector_par(1)]
    elif 'nrs2' in os.path.basename(calfile):
        det_container_list = [spectrograph.get_detector_par(2)]

    ndetectors = 1
    nexp = 1
    cal_data = np.empty((ndetectors, nexp), dtype=object)
    flat_data = np.empty((ndetectors, nexp), dtype=object)

    cal_data[0, 0] = datamodels.open(calfile)
    flat_data[0, 0] = datamodels.open(flatfile)

    slit_names_1 = [slit.name for slit in cal_data[0, 0].slits]
    slit_names_uni = slit_names_1 #np.unique(np.hstack([slit_names_1, slit_names_2]))

    islit = None
    gdslits = slit_names_uni[::-1] if islit is None else [islit]
    iexp_ref = 0

    nrs_calib = []

    for ii, islit in enumerate(gdslits):
        # Container for all the Spec2DObj, different spec2dobj and specobjs for each slit
        # all_spec2d = spec2dobj.AllSpec2DObj()
        # all_spec2d['meta']['bkg_redux'] = bkg_redux
        # all_spec2d['meta']['find_negative'] = bkg_redux
        # Container for the specobjs
        # all_specobjs = specobjs.SpecObjs()

        # TODO this step is being executed repeatedly for each new exposure?
        # Generate the calibrations from the reference exposure
        # TODO This step is only performed with a reference exposure because calwebb has an annoying property that
        # it does not always extract the same subimage spectral pixels for the different dithers in the dither pattern.
        # This seems to be a bug in calwebb, since it is unclear why the subimage calibrations should change.
        # This is problem for 2d coadding, since then the offsets in the detector frame will be not allow one to register
        # the frames. It is possible to fix this by using the RA/DEC images provided by calwebb to determine the
        # actual locations on the sky, which would be preferable. However, this does not appear to be working correctly
        # in calwebb. So for now, we just use the first exposure as the reference exposure for the calibrations.
        CalibrationsNRS_i = NIRSpecSlitCalibrations(
            det_container_list[0], cal_data[0, iexp_ref], flat_data[0, iexp_ref], islit)

        nrs_calib.append(CalibrationsNRS_i)

    return nrs_calib

def get_calib_obj(nrs1_calfile, nrs2_calfile, nrs1_flatfile, nrs2_flatfile):
    spectrograph = load_spectrograph('jwst_nirspec')
    par = spectrograph.default_pypeit_par()
    det_container_list = [spectrograph.get_detector_par(1), spectrograph.get_detector_par(2)]

    ndetectors = 2
    nexp = 1
    cal_data = np.empty((ndetectors, nexp), dtype=object)
    flat_data = np.empty((ndetectors, nexp), dtype=object)

    cal_data[0, 0] = datamodels.open(nrs1_calfile)
    cal_data[1, 0] = datamodels.open(nrs2_calfile)
    flat_data[0, 0] = datamodels.open(nrs1_flatfile)
    flat_data[1, 0] = datamodels.open(nrs2_flatfile)

    slit_names_1 = [slit.name for slit in cal_data[0, 0].slits]
    slit_names_uni = slit_names_1 #np.unique(np.hstack([slit_names_1, slit_names_2]))

    islit = None
    gdslits = slit_names_uni[::-1] if islit is None else [islit]
    iexp_ref = 0

    nrs1_calib = []
    nrs2_calib = []

    for ii, islit in enumerate(gdslits):
        # Container for all the Spec2DObj, different spec2dobj and specobjs for each slit
        # all_spec2d = spec2dobj.AllSpec2DObj()
        # all_spec2d['meta']['bkg_redux'] = bkg_redux
        # all_spec2d['meta']['find_negative'] = bkg_redux
        # Container for the specobjs
        # all_specobjs = specobjs.SpecObjs()

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

        nrs1_calib.append(CalibrationsNRS1)
        nrs2_calib.append(CalibrationsNRS2)

    return nrs1_calib, nrs1_calib


