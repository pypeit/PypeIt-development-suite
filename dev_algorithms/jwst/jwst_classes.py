
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
from jwst_utils import compute_diff, get_cuts, jwst_show_spec2, jwst_show_msa, jwst_proc, jwst_populate_calibs
from pypeit.display import display
from pypeit import slittrace
from pypeit.utils import inverse, fast_running_median
from pypeit.core import findobj_skymask
from pypeit.core import skysub, coadd
from pypeit.core import procimg
from pypeit.spectrographs.util import load_spectrograph
from pypeit.images import pypeitimage
from pypeit import calibrations
from pypeit import find_objects
from pypeit import extraction
from pypeit import spec2dobj


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

    # NIRSPEC 3-point dither
    # dither center
    scifile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_rate.fits')
    bkgfile = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_rate.fits')
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
    scifile  = os.path.join(rawpath_level2, 'jw01117007001_03101_00002_' + detname + '_rate.fits')
    bkgfile1 = os.path.join(rawpath_level2, 'jw01117007001_03101_00003_' + detname + '_rate.fits')
    bkgfile2 = os.path.join(rawpath_level2, 'jw01117007001_03101_00004_' + detname + '_rate.fits')

elif 'G395M' in disperser:
    # Use islit = 37 for nrs1
    # G395M data
    rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/Raw'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/output'

    # NIRSPEC 3-point dither
    scifile = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + detname + '_rate.fits')
    bkgfile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + detname + '_rate.fits')
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
    'imprint_subtract': {'save_results': True},
    'msa_flagging': {'save_results': True},
    'master_background_mos': {'skip': True},
    'srctype': {'source_type':'EXTENDED'},
#    'flat_field': {'skip': True},
    'resample_spec': {'skip': True},
    'extract_1d': {'skip': True},
    'flat_field': {'save_interpolated_flat': True}, # Flats appear to be just nonsense. So skip for now.
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
# This is the output after the MSA flagging step, which is a 2d image that has been background subtraced (disabled),
# imprint subtracted (currently this does nothing since files are missing), and msa flagging
msa_flag_output_file = os.path.join(output_dir, basename +  'msa_flagging.fits')


#s2d_output_file = os.path.join(output_dir, basename + 's2d.fits')
# TESTING
#final2d = datamodels.open(s2d_output_file)
#intflat = None
e2d = datamodels.open(e2d_output_file)
final_multi = datamodels.open(cal_output_file)
intflat_multi = datamodels.open(intflat_output_file)

sci_rate = datamodels.open(scifile)
msa_flagged_rate = datamodels.open(msa_flag_output_file)

sci_data = np.array(msa_flagged_rate.data.T, dtype=float)
nspec, nspat = sci_data.shape
viewer_sci, ch_sci = display.show_image(sci_data, cuts=get_cuts(sci_data), chname='raw rate', clear=True)

spectrograph = load_spectrograph('jwst_nirspec')
reduce_gpm, ra, dec, waveimg, tilts, flatfield, pathloss, barshadow, photom_conversion, calwebb_final, subimg_count, \
slit_left, slit_righ, spec_min, spec_max, meta_list = jwst_populate_calibs(nspec, nspat, final_multi, intflat_multi)
sys.exit(-1)

#intflat = None
#nslits_final2d = len(final2d.slits)

#final2d_reduce_gpm = np.ones(nslits_final2d, dtype=bool)
#for islit, slit in enumerate(final2d.slits):
#    if np.all(np.isnan(slit.data)):
#        final2d_reduce_gpm[islit] = False

#gdslits = np.where(final2d_reduce_gpm)[0]
#nslits = len(gdslits)


#jwst_show_msa(sci_rate, final2d, clear=True)

#jwst_show_spec2(final2d.slits[islit], intflat_slit=intflat.slits[islit], emb=False, clear=False)




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


for ii, islit in enumerate(gdslits):
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
    #intflat_slit = intflat.slits[islit] if intflat is not None else None
    sub_science, sub_sciivar, sub_gpm, sub_base_var, sub_count_scale, sub_tilts, sub_waveimg, sub_thismask, \
    sub_slit_left, sub_slit_righ, t_eff = jwst_proc(e2d.slits[islit], final2d.slits[islit], intflat_slit=intflat.slits[islit])


    science[slit_slice][sub_thismask]     = sub_science[sub_thismask]
    sciivar[slit_slice][sub_thismask]     = sub_sciivar[sub_thismask]
    gpm[slit_slice][sub_thismask]         = sub_gpm[sub_thismask]
    base_var[slit_slice][sub_thismask]    = sub_base_var[sub_thismask]
    count_scale[slit_slice][sub_thismask] = sub_count_scale[sub_thismask]
    tilts[slit_slice][sub_thismask]       = sub_tilts[sub_thismask]
    waveimg[slit_slice][sub_thismask]     = sub_waveimg[sub_thismask]
    slit_left[spec_lo: spec_hi, ii]    = spat_lo + sub_slit_left
    slit_righ[spec_lo: spec_hi, ii]    = spat_lo + sub_slit_righ
    # This is a hack for now until we figure out how to deal with spec_min and spec_max slits
    slit_left[:spec_lo, ii] = spat_lo + sub_slit_left[0]
    slit_left[spec_hi:, ii] = spat_lo + sub_slit_left[-1]
    slit_righ[:spec_lo, ii] = spat_lo + sub_slit_righ[0]
    slit_righ[spec_hi:, ii] = spat_lo + sub_slit_righ[-1]

    spec_min[ii] = spec_lo
    spec_max[ii] = spec_hi

    #display.show_slits(viewer_sci, ch_sci, seg_left, seg_righ, spec_vals=spec_val, pstep=1,
    #                   slit_ids=np.array([int(slit_name)]))
    display.show_slits(viewer_sci, ch_sci, slit_left[spec_lo:spec_hi, ii], slit_righ[spec_lo:spec_hi, ii], spec_vals=spec_val, pstep=1,
                       slit_ids=np.array([int(slit_name)]))


show=True
pypeline='MultiSlit'

slits = slittrace.SlitTraceSet(slit_left, slit_righ, pypeline, detname=detname, nspat=nspat,
                             PYP_SPEC=spectrograph.name, specmin=spec_min, specmax=spec_max)


det_container = spectrograph.get_detector_par(detector)

sciImg = pypeitimage.PypeItImage(image=science,
                                   ivar=sciivar,
                                   base_var=base_var,
                                   img_scale=count_scale,
                                   detector=det_container,
                                   crmask=np.invert(gpm))

slitmask = slits.slit_img()
par = spectrograph.default_pypeit_par()


# TODO add hooks here for manual extraction using the location of the trace as predicted by the JWST meta data. Right now
# this does not work since the WCS and hence trace is garbabe. See coadd2d.py for an example


# Initiate FindObjects object
objFind = find_objects.FindObjects.get_instance(sciImg, slits, spectrograph, par, 'science_coadd2d', tilts=tilts,
                                                manual=None, show=True)
global_sky0, sobjs_obj = objFind.run(show_peaks=True)

# TODO add this as optional to objFind.run()
skymask = objFind.create_skymask(sobjs_obj)
global_sky = objFind.global_skysub(previous_sky=global_sky0, skymask=skymask, show=True)

# Initiate Extract object
extract = extraction.Extract.get_instance(sciImg, slits, sobjs_obj, spectrograph, par, 'science_coadd2d',
                                          tilts=tilts, waveimg=waveimg, basename=basename, show=True)



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
    waveImg = waveimg
    tilts = tilts
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

islit = 10
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

