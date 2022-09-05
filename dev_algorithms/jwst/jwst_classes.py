
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
from jwst_utils import compute_diff, get_cuts, jwst_show_spec2, jwst_show_msa, jwst_proc
from pypeit.display import display
from pypeit import slittrace
from pypeit.utils import inverse, fast_running_median
from pypeit.core import findobj_skymask
from pypeit.core import skysub, coadd
from pypeit.core import procimg
from pypeit.spectrographs.util import load_spectrograph


detname = 'nrs1'
#disperser = 'G395M'
#disperser = 'G235M'
#disperser='PRISM_01133'
disperser='PRISM_01117'
if 'PRISM_01133' in disperser:
    # PRISM data
    rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/Raw'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/output'

    # NIRSPEC 3-point dither
    # dither center
    bkgfile1  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_rate.fits')
    scifile = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_rate.fits')
    bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_rate.fits')

    # dither offset
    #scifile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_rate.fits')
    #bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_rate.fits')
    #bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_rate.fits')
elif 'PRISM_01117' in disperser:
    # PRISM data
    rawpath_level2 = '//Users/joe/jwst_redux/Raw/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/level_12/01117'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb/output'

    # NIRSPEC 3-point dither
    # dither center
    bkgfile1  = os.path.join(rawpath_level2, 'jw01117007001_03101_00002_' + detname + '_rate.fits')
    scifile = os.path.join(rawpath_level2, 'jw01117007001_03101_00003_' + detname + '_rate.fits')
    bkgfile2 = os.path.join(rawpath_level2, 'jw01117007001_03101_00004_' + detname + '_rate.fits')

elif 'G395M' in disperser:
    # Use islit = 37 for nrs1
    # G395M data
    rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/Raw'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/output'

    # NIRSPEC 3-point dither
    bkgfile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + detname + '_rate.fits')
    scifile = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + detname + '_rate.fits')
    bkgfile2 = os.path.join(rawpath_level2, 'jw02736007001_03103_00003_' + detname + '_rate.fits')
elif 'G235M' in disperser:
    # Use islit = 38 for nrs1
    # G235M data
    rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/level_2/'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/output'

    # NIRSPEC 3-point dither
    bkgfile2  = os.path.join(rawpath_level2, 'jw02736007001_03101_00002_' + detname + '_rate.fits')
    bkgfile1 = os.path.join(rawpath_level2, 'jw02736007001_03101_00003_' + detname + '_rate.fits')
    scifile = os.path.join(rawpath_level2, 'jw02736007001_03101_00004_' + detname + '_rate.fits')




# Plot the 2d differnence image

#rawscience, diff = compute_diff(scifile, bkgfile1, bkgfile2, )
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
    #'flat_field': {'skip': True},
    'resample_spec': {'skip': True},
    'extract_1d': {'skip': True},
    'flat_field': {'save_interpolated_flat': True},
}


runflag = False
if runflag:
    spec2 = Spec2Pipeline(steps=param_dict)
    spec2.save_results = True
    spec2.output_dir = output_dir
    result = spec2(scifile)

# Read in the files
intflat_output_file = os.path.join(output_dir, basename + 'interpolatedflat.fits')
e2d_output_file = os.path.join(output_dir, basename + 'extract_2d.fits')
cal_output_file = os.path.join(output_dir, basename + 'cal.fits')
s2d_output_file = os.path.join(output_dir, basename + 's2d.fits')
# TESTING
#final2d = datamodels.open(s2d_output_file)
#intflat = None
#e2d = datamodels.open(e2d_output_file)
final2d = datamodels.open(cal_output_file)
intflat = datamodels.open(intflat_output_file)
nslits = len(final2d.slits)

reduce_gpm = np.ones(nslits, dtype=bool)
for islit, slit in enumerate(final2d.slits):
    if np.all(np.isnan(slit.data)):
        reduce_gpm[islit] = False

gdslits =np.where(reduce_gpm)[0]
nslits_good = len(gdslits)

sci_rate = datamodels.open(scifile)
#jwst_show_msa(sci_rate, final2d, clear=True)

#jwst_show_spec2(final2d.slits[islit], intflat_slit=intflat.slits[islit], emb=False, clear=False)

spectrograph = load_spectrograph('jwst_nirspec')
sci_data = np.array(sci_rate.data.T, dtype=float)
nspec, nspat = sci_data.shape
viewer_sci, ch_sci = display.show_image(sci_data, cuts=get_cuts(sci_data), chname='raw rate', clear=True)

science = np.zeros((nspec, nspat))
sciivar = np.zeros_like(science)
gpm = np.zeros_like(science, dtype=bool)
tilts = np.zeros_like(science)
waveimg = np.zeros_like(science)
count_scale = np.zeros_like(science)
base_var = np.zeros_like(science)
slit_left = np.zeros((nspec, nslits))
slit_righ = np.zeros((nspec, nslits))
spec_min = np.zeros(nslits)
spec_max = np.zeros(nslits)

for islit in gdslits:
    # Read in data print out slit name
    slit_name = final2d.slits[islit].name
    print('Slit={:s}'.format(slit_name))
    nspec_sub, nspat_sub = final2d.slits[islit].data.T.shape

    ########################
    # Plot the image segment being used for each slit
    spec_lo = final2d.slits[islit].xstart - 1
    spec_hi = spec_lo + final2d.slits[islit].xsize
    spat_lo = final2d.slits[islit].ystart - 1
    spat_hi = spat_lo + final2d.slits[islit].ysize
    # This is the segment of the 2d image
    slit_slice = np.s_[spec_lo: spec_hi, spat_lo: spat_hi]

    seg_left = np.full(nspec_sub, spat_lo)
    seg_righ = np.full(nspec_sub, spat_hi)
    spec_val = spec_lo + np.arange(spec_hi - spec_lo)


    sub_science, sub_sciivar, sub_gpm, sub_base_var, sub_count_scale, sub_tilts, sub_waveimg, sub_thismask, \
    sub_slit_left, sub_slit_righ, t_eff = jwst_proc(final2d.slits[islit], intflat_slit=intflat.slits[islit])

    science[slit_slice][sub_thismask]     = sub_science[sub_thismask]
    sciivar[slit_slice][sub_thismask]     = sub_sciivar[sub_thismask]
    gpm[slit_slice][sub_thismask]         = sub_gpm[sub_thismask]
    base_var[slit_slice][sub_thismask]    = sub_base_var[sub_thismask]
    count_scale[slit_slice][sub_thismask] = sub_count_scale[sub_thismask]
    tilts[slit_slice][sub_thismask]       = sub_tilts[sub_thismask]
    waveimg[slit_slice][sub_thismask]     = sub_waveimg[sub_thismask]
    slit_left[spec_lo: spec_hi, islit]    = spat_lo + sub_slit_left
    slit_righ[spec_lo: spec_hi, islit]    = spat_lo + sub_slit_righ
    spec_min[islit] = spec_lo
    spec_max[islit] = spec_hi

    #display.show_slits(viewer_sci, ch_sci, seg_left, seg_righ, spec_vals=spec_val, pstep=1,
    #                   slit_ids=np.array([int(slit_name)]))
    display.show_slits(viewer_sci, ch_sci, slit_left[spec_lo:spec_hi, islit], slit_righ[spec_lo:spec_hi, islit], spec_vals=spec_val, pstep=1,
                       slit_ids=np.array([int(slit_name)]))

pypeline='MultiSlit'

slits = slittrace.SlitTraceSet(slit_left, slit_righ, pypeline, detname=detname, nspat=nspat,
                             PYP_SPEC=spectrograph.name, specmin=spec_min, specmax=spec_max)
# master_key=self.stack_dict['master_key_dict']['trace'],
# master_dir=self.master_dir)



