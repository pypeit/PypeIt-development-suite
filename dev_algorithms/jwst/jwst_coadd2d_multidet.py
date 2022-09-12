
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
exp_list = []
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
        rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/Raw'
        output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/output'
        pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/pypeit'

        # NIRSPEC 3-point dither
        scifile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + detname + '_rate.fits')
        scifile2 = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + detname + '_rate.fits')
        scifile3 = os.path.join(rawpath_level2, 'jw02736007001_03103_00003_' + detname + '_rate.fits')
    elif 'G235M' in disperser:
        # Use islit = 38 for nrs1
        # G235M data
        rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/level_2/'
        output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G235M/calwebb/output'
        pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G235M/calwebb/pypeit/'

        # NIRSPEC 3-point dither
        scifile1  = os.path.join(rawpath_level2, 'jw02736007001_03101_00002_' + detname + '_rate.fits')
        scifile2 = os.path.join(rawpath_level2, 'jw02736007001_03101_00003_' + detname + '_rate.fits')
        scifile3 = os.path.join(rawpath_level2, 'jw02736007001_03101_00004_' + detname + '_rate.fits')
    exp_list.append([scifile1, scifile2, scifile3])

scifiles_1 = exp_list[0]
scifiles_2 = exp_list[1]
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
    'srctype': {'source_type':'EXTENDED'},
#    'flat_field': {'skip': True},
    'resample_spec': {'skip': True},
    'extract_1d': {'skip': True},
    'flat_field': {'save_interpolated_flat': True}, # Flats appear to be just nonsense. So skip for now.
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
        spec2 = Spec2Pipeline(steps=param_dict)
        spec2.save_results = True
        spec2.output_dir = output_dir
        result = spec2(sci)



# Output file names
intflat_output_files_1 = []
e2d_output_files_1 =[]
cal_output_files_1 =[]

intflat_output_files_2 = []
e2d_output_files_2 =[]
cal_output_files_2 =[]
for base1, base2 in zip(basenames_1, basenames_2):
    e2d_output_files_1.append(os.path.join(output_dir, base1 + '_extract_2d.fits'))
    intflat_output_files_1.append(os.path.join(output_dir, base1 + '_interpolatedflat.fits'))
    cal_output_files_1.append(os.path.join(output_dir, base1 + '_cal.fits'))

    e2d_output_files_2.append(os.path.join(output_dir, base2 + '_extract_2d.fits'))
    intflat_output_files_2.append(os.path.join(output_dir, base2 + '_interpolatedflat.fits'))
    cal_output_files_2.append(os.path.join(output_dir, base2 + '_cal.fits'))

# Read in calwebb outputs for everytihng

# Read in multi exposure calwebb outputs
e2d_multi_list_1 = []
intflat_multi_list_1 = []
final_multi_list_1 = []
e2d_multi_list_2 = []
intflat_multi_list_2 = []
final_multi_list_2 = []
nslits_1 = np.zeros(nexp, dtype=int)
nslits_2 = np.zeros(nexp, dtype=int)
t_eff = np.zeros(nexp, dtype=float)

for iexp in range(nexp):
    # Open some JWST data models
    e2d_multi_list_1.append(datamodels.open(e2d_output_files_1[iexp]))
    intflat_multi_list_1.append(datamodels.open(intflat_output_files_1[iexp]))
    final_multi_list_1.append(datamodels.open(cal_output_files_1[iexp]))

    e2d_multi_list_2.append(datamodels.open(e2d_output_files_2[iexp]))
    intflat_multi_list_2.append(datamodels.open(intflat_output_files_2[iexp]))
    final_multi_list_2.append(datamodels.open(cal_output_files_2[iexp]))

    t_eff[iexp] = e2d_multi_list_1[iexp].meta.exposure.effective_exposure_time
    nslits_1[iexp] = len(final_multi_list_1[iexp].slits)
    nslits_2[iexp] = len(final_multi_list_2[iexp].slits)

#
# the multi-exposure calwebb outputs


#slit_name = slit_names_1[islit_1]
#islit_2 = np.where(np.array(slit_names_2) == slit_name)[0][0]



# Some pypeit things
spectrograph = load_spectrograph('jwst_nirspec')
par = spectrograph.default_pypeit_par()
det_container_list = [spectrograph.get_detector_par(1), spectrograph.get_detector_par(2)]

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
filename_first = os.path.basename(scifiles_1[0])
filename_last = os.path.basename(scifiles_1[-1])
split_first = filename_first.split('_')
split_last = filename_last.split('_')
prefix_first = ''
for ii in range(len(split_first)-3):
    prefix_first += split_first[ii] + "_"
out_filename = prefix_first + split_first[-3] + "-" + split_last[-3]


show=True

offsets_pixels_list = [[0, 5.0, -5.0], [0, 5.0, -5.0]]


spec_samp_fact = 1.0
spat_samp_fact = 1.0


slit_names_1 = [int(slit.name) for slit in e2d_multi_list_1[0].slits]
slit_names_2 = [int(slit.name) for slit in e2d_multi_list_2[0].slits]
slit_names_uni = np.unique(np.hstack([slit_names_1, slit_names_2]))

# Loop over slits
islit = 10
#islit = 64
#islit=None
gdslits = slit_names_uni[::-1] if islit is None else [islit]
bad_slits = []

# First index is detector, second index is exposure
e2d_multi_list = [e2d_multi_list_1, e2d_multi_list_2]
intflat_multi_list = [intflat_multi_list_1, intflat_multi_list_2]
final_multi_list = [final_multi_list_1, final_multi_list_2]
slit_names_list = [slit_names_1, slit_names_2]
kludge_err = 1.5



# Loop over all slits, create a list of spec2d objects and run 2d coadd
for islit in gdslits:
    slit_name_str = str(islit)
    spec2d_list = []
    offsets_pixels = []
    if show:
        display.clear_all()
    for idet in range(2):
        for iexp in range(nexp):
            indx = np.where(np.array(slit_names_list[idet]) == islit)[0]
            if len(indx) > 0:
                ii = indx[0]
                waveimg, tilts, slit_left, slit_righ, science, sciivar, gpm, base_var, count_scale, \
                slit_left_orig, slit_righ_orig, spec_vals_orig =jwst_proc(
                    t_eff[iexp], e2d_multi_list[idet][iexp].slits[ii], final_multi_list[idet][iexp].slits[ii],
                    intflat_multi_list[idet][iexp].slits[ii], noise_floor=par['scienceframe']['process']['noise_floor'],
                    kludge_err=kludge_err, show=(iexp == 0), ronoise=det_container_list[idet].ronoise)

                # If no finite pixels in the waveimg then skip this slit
                if not np.any(gpm):
                    bad_slits.append(islit)
                    continue

                if show and (iexp ==0):
                    sci_rate = datamodels.open(scifiles[idet][iexp])
                    sci_data = np.array(sci_rate.data.T, dtype=float)
                    viewer_sci, ch_sci = display.show_image(sci_data, cuts=get_cuts(sci_data),
                                                            chname='raw rate_iexp_{:d}_idet_{:d}'.format(iexp, idet), clear=False)
                    display.show_slits(viewer_sci, ch_sci, slit_left_orig, slit_righ_orig, spec_vals=spec_vals_orig, pstep=1,
                                       slit_ids=np.array([islit]))

                nspec, nspat = waveimg.shape
                slits = slittrace.SlitTraceSet(slit_left, slit_righ, pypeline, detname=det_container_list[idet].name, nspat=nspat,
                                               PYP_SPEC=spectrograph.name)
                slits.maskdef_id = np.array([islit])

                # Construct the Spec2DObj with the positive image
                spec2DObj = spec2dobj.Spec2DObj(sciimg=science,
                                                ivarraw=sciivar,
                                                skymodel=np.zeros_like(science),
                                                objmodel=np.zeros_like(science),
                                                ivarmodel=sciivar,
                                                scaleimg=None,
                                                waveimg=waveimg,
                                                bpmmask=np.logical_not(gpm).astype(int),
                                                detector=det_container_list[idet],
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
            else:
                continue

    if len(spec2d_list) > 0:
        basename = '{:s}_{:s}'.format(out_filename, 'slit' + slit_name_str)


        #TODO Not sure what to do with the detector container here
        # Instantiate Coadd2d
        coAdd = coadd2d.CoAdd2D.get_instance(spec2d_list, spectrograph, par, det=det_container_list[0].det,
                                             offsets=offsets_pixels, weights='uniform',
                                             spec_samp_fact=spec_samp_fact,
                                             spat_samp_fact=spat_samp_fact,
                                             bkg_redux=False, debug=show)

        coadd_dict_list = coAdd.coadd(only_slits=None, interp_dspat=True)


        # TODO THe rebin2d images have all kinds of holes in them, which happens before any kind of sigma clipping.
        # I need to get to the bottom of this!!
        #display.show_image(coadd_dict_list[0]['rebin_sciimg_stack'][3, :, :], chname='img4')

        # Create the pseudo images
        pseudo_dict = coAdd.create_pseudo_image(coadd_dict_list)

        sciimg_coadd, sciivar_coadd, skymodel_coadd, objmodel_coadd, ivarmodel_coadd, \
        outmask_coadd, sobjs_coadd, detector_coadd, slits_coadd, tilts_coadd, waveimg_coadd = coAdd.reduce(
            pseudo_dict, global_sky_subtract=True, show_skysub_fit=True, show=show, clear_ginga=False, show_peaks=show, basename=basename)

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
            spec2DObj_coadd.gen_qa()


            slitmask_coadd = slits_coadd.slit_img(initial=False, flexure=None, exclude_flag=None)
            slitord_id = slits_coadd.slitord_id[0]
            thismask = slitmask_coadd == slitord_id

            gpm_extract = spec2DObj_coadd.bpmmask == 0
            # Make a plot of the residuals for a random slit
            chi = (spec2DObj_coadd.sciimg - spec2DObj_coadd.objmodel - spec2DObj_coadd.skymodel) * np.sqrt(spec2DObj_coadd.ivarmodel) * gpm_extract

            maskchi = thismask & gpm_extract

            n_bins = 50
            sig_range = 7.0
            binsize = 2.0 * sig_range / n_bins
            bins_histo = -sig_range + np.arange(n_bins) * binsize + binsize / 2.0

            xvals = np.arange(-10.0, 10, 0.02)
            gauss = scipy.stats.norm(loc=0.0, scale=1.0)
            gauss_corr = scipy.stats.norm(loc=0.0, scale=1.0)

            sigma_corr, maskchi = coadd.renormalize_errors(chi, maskchi, max_corr=20.0, title='jwst_sigma_corr', debug=True)


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



