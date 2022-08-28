
import os
import numpy as np
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

# PypeIt imports
from jwst_utils import show_diff, get_cuts
from pypeit.display import display


rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/Raw'
output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/output'


# NIRSPEC 3-point dither
det = 'nrs1'
scifile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + det + '_rate.fits')
bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + det + '_rate.fits')
bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + det + '_rate.fits')
# Plot the 2d differnence image
sci, diff = show_diff(scifile, bkgfile1, bkgfile2)

#viewer_diff, ch_diff = display.show_image(diff.T, cuts=get_cuts(diff), chname='diff2d')
viewer_sci,  ch_sci = display.show_image(sci.T, cuts=get_cuts(sci), chname='raw', wcs_match=True)
basename = os.path.basename(scifile).replace('rate.fits', '')

param_dict = {
    'bkg_subtract': {'skip': True},
    'srctype': {'source_type':'EXTENDED'},
    'flat_field': {'save_interpolated_flat':True},
}


runflag = True
if runflag:
    spec2 = Spec2Pipeline()
    spec2.save_results = True
    spec2.output_dir = output_dir
    result = spec2(scifile)


