
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


#detname = 'nrs1'
#detector = 1 if 'nrs1' in detname else 2

disperser = 'G395M'
#disperser = 'G235M'
#disperser='PRISM_01133'
detectors = ['nrs1', 'nrs2']
#disperser='PRISM_01117'
scifiles = []
for detname in detectors:
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
    scifiles += [scifile1, scifile2, scifile3]

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


e2d_nrs1 = datamodels.open(e2d_output_files[0])
e2d_nrs2 = datamodels.open(e2d_output_files[3])

final_nrs1 = datamodels.open(cal_output_files[0])
final_nrs2 = datamodels.open(cal_output_files[3])

intflat_nrs1 = datamodels.open(intflat_output_files[0])
intflat_nrs2 = datamodels.open(intflat_output_files[3])

islit_1 = 37
slit_names_1 = [slit.name for slit in e2d_nrs1.slits]
slit_names_2 = [slit.name for slit in e2d_nrs2.slits]
slit_name = slit_names_1[islit_1]
islit_2 = np.where(np.array(slit_names_2) == slit_name)[0][0]


slit_left_1, slit_righ_1, slit_left_orig_1, slit_righ_orig_1, spec_vals_orig_1, src_trace_ra_1, src_trace_dec_1, \
rate_1, rate_var_rnoise_1, rate_var_poisson_1, rate_var_tot_1, dq_1, ra_1, dec_1, waveimg_1, tilts_1, flatfield_1, pathloss_1, \
barshadow_1, photom_conversion_1, final_1 = jwst_extract_subimgs(e2d_nrs1.slits[islit_1], final_nrs1.slits[islit_1], intflat_nrs1.slits[islit_1])

slit_left_2, slit_righ_2, slit_left_orig_2, slit_righ_orig_2, spec_vals_orig_2, src_trace_ra_2, src_trace_dec_2, \
rate_2, rate_var_rnoise_2, rate_var_poisson_2, rate_var_tot_2, dq_2, ra_2, dec_2, waveimg_2, tilts_2, flatfield_2, pathloss_2, \
barshadow_2, photom_conversion_2, final_2 = jwst_extract_subimgs(e2d_nrs2.slits[islit_2], final_nrs2.slits[islit_2], intflat_nrs2.slits[islit_2])


# Find common slit to both
nspec_1 = rate_1.shape[0]
nspec_2 = rate_2.shape[0]
min_slit_ra_1 = np.zeros(nspec_1)
max_slit_ra_1 = np.zeros(nspec_1)
min_slit_ra_2 = np.zeros(nspec_2)
max_slit_ra_2 = np.zeros(nspec_2)

for ii in range(nspec_1):
    ra_seg = ra_1[ii, :]
    min_slit_ra_1[ii] = np.min(ra_seg[np.isfinite(ra_seg)])
    max_slit_ra_1[ii] = np.max(ra_seg[np.isfinite(ra_seg)])

for ii in range(nspec_2):
    ra_seg = ra_2[ii, :]
    min_slit_ra_2[ii] = np.min(ra_seg[np.isfinite(ra_seg)])
    max_slit_ra_2[ii] = np.max(ra_seg[np.isfinite(ra_seg)])



sci_nrs1 = datamodels.open(scifiles[0])
sci_nrs2 = datamodels.open(scifiles[3])

nspec, nspat = sci_nrs1.data.shape
gap = 180
nspec_mosaic = 2*nspec + gap

mosaic = np.zeros((nspec_mosaic, nspat))
mosaic[0:nspec, :] = sci_nrs1.data.T
mosaic[nspec+gap:, :] = sci_nrs2.data.T

display.show_image(mosaic, chname='Mosaic')


