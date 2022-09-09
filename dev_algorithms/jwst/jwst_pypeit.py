
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
DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']


detname = 'nrs1'
detector = 1 if 'nrs1' in detname else 2
#disperser = 'G395M'
#disperser = 'G235M'
#disperser='PRISM_01133'
disperser='PRISM_01117'
if 'PRISM_01133' in disperser:
    # PRISM data
    rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/Raw'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/output'
    pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/pypeit'


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
    pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb/pypeit'


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
    pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/pypeit'

    # NIRSPEC 3-point dither
    scifile = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + detname + '_rate.fits')
    bkgfile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + detname + '_rate.fits')
    bkgfile2 = os.path.join(rawpath_level2, 'jw02736007001_03103_00003_' + detname + '_rate.fits')
elif 'G235M' in disperser:
    # Use islit = 38 for nrs1
    # G235M data
    rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/level_2/'
    output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/output'
    pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/calwebb/pypeit/'

    # NIRSPEC 3-point dither
    bkgfile2  = os.path.join(rawpath_level2, 'jw02736007001_03101_00002_' + detname + '_rate.fits')
    bkgfile1 = os.path.join(rawpath_level2, 'jw02736007001_03101_00003_' + detname + '_rate.fits')
    scifile = os.path.join(rawpath_level2, 'jw02736007001_03101_00004_' + detname + '_rate.fits')

# Plot the 2d differnence image
#sci, diff = compute_diff(scifile, bkgfile1, bkgfile2)
#viewer_diff, ch_diff = display.show_image(diff.T, cuts=get_cuts(diff), chname='diff2d')
#viewer_sci,  ch_sci = display.show_image(sci.T, cuts=get_cuts(sci), chname='science')
#sys.exit(-1)


# Plot the 2d differnence image

#rawscience, diff = compute_diff(scifile, bkgfile1, bkgfile2, )
#sci_rate = datamodels.open(scifile)
#viewer_diff, ch_diff = display.show_image(diff.T, cuts=get_cuts(diff), chname='diff2d')
#viewer_sci,  ch_sci = display.show_image(sci.T, cuts=get_cuts(sci), chname='raw', wcs_match=True)
basename = os.path.basename(scifile).replace('_rate.fits', '')

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

# Run the pipeline
runflag = False
if runflag:
    spec2 = Spec2Pipeline(steps=param_dict)
    spec2.save_results = True
    spec2.output_dir = output_dir
    result = spec2(scifile)

# The output filen ames
intflat_output_file = os.path.join(output_dir, basename + '_interpolatedflat.fits')
e2d_output_file = os.path.join(output_dir, basename + '_extract_2d.fits')
cal_output_file = os.path.join(output_dir, basename + '_cal.fits')
head2d = fits.getheader(cal_output_file)

# Some pypeit things
spectrograph = load_spectrograph('jwst_nirspec')
par = spectrograph.default_pypeit_par()
det_container = spectrograph.get_detector_par(detector)
pypeline = 'MultiSlit'


# Open some JWST data models
e2d_multi = datamodels.open(e2d_output_file)
final_multi = datamodels.open(cal_output_file)
intflat_multi = datamodels.open(intflat_output_file)
nslits = len(final_multi.slits)

sci_rate = datamodels.open(scifile)
sci_data = np.array(sci_rate.data.T, dtype=float)
nspec_raw, nspat_raw = sci_data.shape

kludge_err = 1.0
t_eff = e2d_multi.meta.exposure.effective_exposure_time
# TODO I don't know how the t_eff quantity is defined. Better would be some proxy for the exposure time per pixel
# The science data is divided by (flatfield*pathloss*barshadow) and then multiplied by photom_conversion. Since
# we work in units of counts, we divide by the photom conversion and multiply by t_eff.

show=True
if show:
    display.connect_to_ginga(raise_err=True, allow_new=True)

#islit = 37
islit = None
gdslits = np.arange(nslits) if islit is None else [islit]
bad_slits = []

# container for specobjs
all_specobjs = specobjs.SpecObjs()
# container for spec2dobj
all_spec2d = spec2dobj.AllSpec2DObj()
# set some meta
all_spec2d['meta']['bkg_redux'] = False
all_spec2d['meta']['find_negative'] = False

for ii, islit in enumerate(gdslits):
    # Read in data print out slit name
    slit_name = final_multi.slits[islit].name

    slit_left, slit_righ, slit_left_orig, slit_righ_orig, spec_vals_orig, src_trace_ra, src_trace_dec, \
    rate, rate_var_rnoise, rate_var_poisson, rate_var_tot, dq, \
    ra, dec, waveimg, tilts, flatfield, pathloss, barshadow, photom_conversion, final = jwst_extract_subimgs(
        e2d_multi.slits[islit], final_multi.slits[islit], intflat_multi.slits[islit])

    # Now deal with the image processing
    nspec, nspat = rate.shape
    finitemask = np.isfinite(waveimg)
    if not np.any(finitemask):
        bad_slits.append(islit)
        continue

    # Now perform the image processing
    raw_counts = rate*t_eff
    raw_var_poisson = kludge_err**2*rate_var_poisson*t_eff**2
    raw_var_rnoise = kludge_err**2*rate_var_rnoise*t_eff**2
    # Is this correct? I'm not sure I should be using their poisson variance for the noise floor
    raw_var = procimg.variance_model(raw_var_rnoise, counts = raw_var_poisson, noise_floor=par['scienceframe']['process']['noise_floor'])
    # TODO This  is a hack until I can understand how to get rid of the hot pixels in the JWST variance arrays using DQ flags
    raw_gpm = (raw_var_rnoise < 2*det_container.ronoise[0]**2) & (raw_var_poisson < det_container.saturation)
    #raw_var_poisson + raw_var_rnoise # TODO Leaving out problematic flat field term from pipeline

    # This is the conversion between final2d and e2d, i.e. final2d = jwst_scale*e2d
    # total_flat = flatfield*pathloss*barshadow
    #flux_to_counts = t_eff / photom_conversion  # This converts s2d outputs of flux to counts.
    #jwst_scale = photom_conversion/flatfield/pathloss/barshadow
    total_flat = pathloss*barshadow
    total_flat_square = np.square(total_flat)

    count_scale = inverse(total_flat)  # This is the quantity that goes into PypeIt for var modeling
    science, flat_bpm = flat.flatfield(raw_counts, total_flat)
    var_poisson, _ = flat.flatfield(raw_var_poisson, total_flat**2)
    base_var, _ = flat.flatfield(raw_var_rnoise, total_flat**2)
    var, _ = flat.flatfield(raw_var, total_flat**2)
    sciivar = inverse(var)
    dq_gpm = np.logical_not(dq & DO_NOT_USE)
    gpm = finitemask & dq_gpm & np.logical_not(flat_bpm) & (sciivar > 0.0) & raw_gpm

    if show:
        viewer_sci, ch_sci = display.show_image(sci_data, cuts=get_cuts(sci_data), chname='raw rate', clear=True)
        display.show_slits(viewer_sci, ch_sci, slit_left_orig, slit_righ_orig, spec_vals=spec_vals_orig, pstep=1,
                           slit_ids=np.array([int(slit_name)]))
        viewer_data, ch_data = display.show_image(science, waveimg = waveimg, cuts = get_cuts(science), chname='science')
        viewer_wave, ch_wave = display.show_image(waveimg, chname='wave')
        viewer_tilts, ch_tilts = display.show_image(tilts, waveimg=waveimg, chname='tilts')
        viewer_flat, ch_flat = display.show_image(flatfield, waveimg=waveimg, chname='flat')
        viewer_path, ch_path = display.show_image(pathloss, waveimg=waveimg, chname='pathloss')
        viewer_bar , ch_bar  = display.show_image(barshadow, waveimg=waveimg, chname='barshadow')
        display.show_trace(viewer_data, ch_data, src_trace_ra, 'trace-RA', color='#f0e442')
        display.show_trace(viewer_data, ch_data, src_trace_dec, 'trace-DEC', color='#f0e442')


    slits = slittrace.SlitTraceSet(slit_left, slit_righ, pypeline, detname=detname, nspat=nspat,
                                   PYP_SPEC=spectrograph.name)
    slits.maskdef_id = np.array([int(slit_name)])
    # TODO assign the slitord_id to the JWST meta data slit name?

    # I cannot figure out how to instantiate this class from an existing mask. Same problem exists in 2d coadd as well
    # as QA. I think it should be done with fullmask but cannot figure out the bits.
    sciImg = pypeitimage.PypeItImage(image=science,
                                     ivar=sciivar*gpm,
                                     base_var=base_var,
                                     img_scale=count_scale,
                                     detector=det_container,
                                     bpm=np.logical_not(gpm).astype(int), # TODO Why is this an integer??
                                     crmask=np.logical_not(gpm))

    slitmask = slits.slit_img()
    par = spectrograph.default_pypeit_par()

    # TODO add hooks here for manual extraction using the location of the trace as predicted by the JWST meta data. Right now
    # this does not work since the WCS and hence trace is garbabe. See coadd2d.py for an example

    # Initiate FindObjects object
    objFind = find_objects.FindObjects.get_instance(sciImg, slits, spectrograph, par, 'science_coadd2d', tilts=tilts,
                                                    manual=None, show=True, clear_ginga=False)
    global_sky0, sobjs_obj = objFind.run(show_peaks=True)

    # TODO add this as optional to objFind.run()
    skymask = objFind.create_skymask(sobjs_obj)
    global_sky = objFind.global_skysub(previous_sky=global_sky0, skymask=skymask, show=True)

    # Initiate Extract object
    extract = extraction.Extract.get_instance(sciImg, slits, sobjs_obj, spectrograph, par, 'science_coadd2d',
                                              tilts=tilts, waveimg=waveimg, basename=basename, show=True)

    if not par['reduce']['extraction']['skip_extraction']:
        skymodel, objmodel, ivarmodel, outmask, sobjs, waveimg, tilts = extract.run(global_sky*0.0, sobjs_obj)
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

    # Now append
    all_specobjs.add_sobj(sobjs)

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

    all_spec2d[slit_name] = spec2DObj



    if show:
        slitmask = slits.slit_img(initial=False, flexure=None, exclude_flag=None)
        slitord_id = slits.slitord_id[0]
        thismask = slitmask == slitord_id

        gpm_extract = spec2DObj.bpmmask == 0
        # Make a plot of the residuals for a random slit
        chi = (spec2DObj.sciimg - spec2DObj.objmodel - spec2DObj.skymodel) * np.sqrt(spec2DObj.ivarmodel) * gpm_extract


        maskchi = thismask & gpm_extract

        n_bins = 50
        sig_range = 7.0
        binsize = 2.0 * sig_range / n_bins
        bins_histo = -sig_range + np.arange(n_bins) * binsize + binsize / 2.0

        xvals = np.arange(-10.0, 10, 0.02)
        gauss = scipy.stats.norm(loc=0.0, scale=1.0)
        gauss_corr = scipy.stats.norm(loc=0.0, scale=1.0)


        sigma_corr, maskchi = coadd.renormalize_errors(chi, maskchi, max_corr=20.0, title='jwst_sigma_corr', debug=True)



# Make the new Science dir
# TODO: This needs to be defined by the user
scipath = os.path.join(pypeit_output_dir, 'Science')
if not os.path.isdir(scipath):
    msgs.info('Creating directory for Science output: {0}'.format(scipath))
    os.makedirs(scipath)

# THE FOLLOWING MIMICS THE CODE IN pypeit.save_exposure()
subheader = spectrograph.subheader_for_spec(head2d, head2d, allow_missing=True)
# Write spec1D
if all_specobjs.nobj > 0:
    outfile1d = os.path.join(scipath, 'spec1d_{:s}.fits'.format(basename))
    all_specobjs.write_to_fits(subheader, outfile1d)

    # Info
    outfiletxt = os.path.join(scipath, 'spec1d_{:s}.txt'.format(basename))
    sobjs = specobjs.SpecObjs.from_fitsfile(outfile1d, chk_version=False)
    sobjs.write_info(outfiletxt, spectrograph.pypeline)

    # Build header for spec2d
    outfile2d = os.path.join(scipath, 'spec2d_{:s}.fits'.format(basename))
    # TODO For the moment hack so that we can write this out
    pri_hdr = all_spec2d.build_primary_hdr(head2d, spectrograph, subheader=subheader,
                                           redux_path=None, master_key_dict=None, master_dir=None)

    # Write spec2d
    all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr)
