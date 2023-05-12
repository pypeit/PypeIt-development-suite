# Copied from ../jwst_diff_coadd2d.py 10 Jan 2023

import os
from pathlib import Path

from IPython import embed

import numpy as np
import scipy
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

# JWST imports
from jwst import datamodels

# PypeIt imports
from pypeit.display import display
from pypeit import specobjs
from pypeit import slittrace
from pypeit.utils import inverse
from pypeit.core import coadd, combine

from pypeit.spectrographs.util import load_spectrograph
from pypeit import msgs
from pypeit import spec2dobj
from pypeit import coadd2d
from pypeit.images.pypeitimage import PypeItImage
from pypeit.scripts.show_2dspec import show_trace

from jwst_utils import get_cuts, jwst_show_spec2, jwst_show_msa, jwst_proc, jwst_extract_subimgs

# Raw directory
raw_dir = Path('/Users/westfall/Work/JWST/data/10Jan2023/2756_nirspec/raw').resolve()

# JWST spec2 outputs
jwst_rdx_dir = raw_dir.parent / 'jwst_rdx'
if not jwst_rdx_dir.exists():
    jwst_rdx_dir.mkdir(parents=True)

# PypeIt Directory setup
pypeit_rdx_dir = raw_dir.parent / 'pypeit_rdx_2'
pypeit_sci_dir = pypeit_rdx_dir / 'Science'
if not pypeit_sci_dir.exists():
    msgs.info(f'Creating directory for Science output: {pypeit_sci_dir}')
    pypeit_sci_dir.mkdir(parents=True)
pypeit_qa_dir = pypeit_rdx_dir / 'QA'
if not pypeit_qa_dir.exists():
    msgs.info(f'Creating directory for QA output: {pypeit_qa_dir}')
    pypeit_qa_dir.mkdir(parents=True)
pypeit_png_dir = pypeit_qa_dir / 'PNGs'
if not pypeit_png_dir.exists():
    msgs.info(f'Creating directory for QA output: {pypeit_qa_dir}')
    pypeit_png_dir.mkdir(parents=True)

# PypeIt spectrograph and parameters
spectrograph = load_spectrograph('jwst_nirspec')
# TODO: Why isn't this defined by the spectrograph instance
spectrograph.pypeline = 'MultiSlit'  
# TODO: config_specific_par()?
par = spectrograph.default_pypeit_par()

# If diff_redux is False, the code will model the sky and the object profile and
# perform optimal extraction.  If diff_redux is True, the code will difference
# image and simply boxcar extract (optimal not implemented yet).
diff_redux = True

# Parameter alterations
par['rdx']['redux_path'] = str(pypeit_rdx_dir)
par['rdx']['qadir'] = 'QA'
if diff_redux:
    # Do not sky-subtract when object finding
    par['reduce']['findobj']['skip_skysub'] = True
    # Skip local_skysubtraction and profile fitting
    par['reduce']['extraction']['skip_optimal'] = True

# Show stuff to the screen/ginga
show = True

# Sampling factors
spec_samp_fact = 2.0
spat_samp_fact = 1.5

# Kludge JWST spec2 error
kludge_err = 1.5
noise_floor = par['scienceframe']['process']['noise_floor']

# NOTE: Annoyingly, source_name is not defined in the flat field datamodel...
# Use source name instead of slit name because slit names change between
# exposure sequence 3101 and 3103, but the source name doesn't.
# Image B of source JD1
#exp_seq = 'jw02756001001_03101'
#exp_seq = 'jw02756001001_03103'
#sources = ['2756_202']
#par['reduce']['findobj']['snr_thresh'] = 2.0
# YD8
exp_seq = 'jw02756001001_03101'
sources = ['2756_10025']
par['reduce']['findobj']['snr_thresh'] = 7.
# YD7 or YD4 (unsure): Check against their Table 1
# TODO: Check the exp is correct.
#exp_seq = 'jw02756001001_03103'
#sources = ['2756_100002']
#par['reduce']['findobj']['snr_thresh'] = 6.0
# Bright source
# TODO: Check the exp is correct.
#exp_seq = 'jw02756001001_03103'
#sources = ['2756_320035']

# TODO:
#   - Use association files to setup which files to use as science and background

# Get the list of JWST rate files in the dither sequence
#detectors = ['nrs1', 'nrs2']
#det_par = [spectrograph.get_detector_par(1), spectrograph.get_detector_par(2)]
detectors = ['nrs2']
det_par = [spectrograph.get_detector_par(2)]
rate_files = []
for d in detectors:
    rate_files += [list(sorted(raw_dir.glob(f'{exp_seq}*{d}_rate.fits')))]
# Must be the same number of exposures per detector, so the length of each list
# in rate_files is the same. i.e., can treat rate_files as an array
rate_files = np.asarray(rate_files)
bkg_indices = [1, 0, 0]
#bkg_indices = [[1,2], [0,2], [0,1]]
#rate_files = np.asarray(rate_files)[:,1:]
#bkg_indices = [1, 0]
mode = 'MSA'

# The shape of the file list is (ndet,nexp), number of detectors (2) by the
# number of exposures
ndet, nexp = rate_files.shape

dither_offsets = np.zeros((ndet,nexp), dtype=float)
# TODO: This probably isn't correct.  I.e., need to know offsets and slit
# position angle.
for iexp in range(nexp):
    with fits.open(rate_files[0,iexp]) as hdu:
        dither_offsets[0,iexp] = hdu[0].header['YOFFSET']
for idet in range(1,ndet):
    dither_offsets[idet] = dither_offsets[0]
dither_offsets_pixels = dither_offsets.copy()
for idet in range(ndet):
    dither_offsets_pixels[idet] /= det_par[idet].platescale
# NOTE: Sign convention requires this calculation of the offset
dither_offsets_pixels = dither_offsets_pixels[:,0,None] - dither_offsets_pixels
print(dither_offsets_pixels)

# Hard-code the dither pattern.  Dither patterns are described here: 
#   https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-operations/nirspec-dithers-and-nods/nirspec-mos-dither-and-nod-patterns
# This is the "3-shutter Nod" pattern
# TODO: Can this be pulled from the header?
#       - Yes, use hdu[0].header['XOFFSET'] and hdu[0].header['YOFFSET']
#offsets_pixels_list = [[0, 5.29, -5.29], [0, 5.29, -5.29]]  # Offsets are in pixels; pixel scale is 0.1 "/pixel

# Set the exposure base name
basenames = np.array([f.name.replace('_rate.fits','') for f in rate_files.flat]).reshape(ndet,nexp)

# Dither number
dither = np.array([n.split('_')[2] for n in basenames.flat]).reshape(ndet,nexp)

# Output root name
out_root = f'{exp_seq}_{dither[0][0]}-{dither[0][-1]}'

# Output files used from JWST spec2
msa_files = np.array([jwst_rdx_dir / f'{b}_msa_flagging.fits' for b in basenames.flat]
                     ).reshape(ndet,nexp)
if not all([f.exists for f in msa_files.flat]):
    msgs.error('Missing msa_flagging files.')
flat_files = np.array([jwst_rdx_dir / f'{b}_interpolatedflat.fits' for b in basenames.flat]
                      ).reshape(ndet,nexp)
if not all([f.exists for f in flat_files.flat]):
    msgs.error('Missing interpolatedflat files.')
cal_files = np.array([jwst_rdx_dir / f'{b}_cal.fits' for b in basenames.flat]).reshape(ndet,nexp)
if not all([f.exists for f in cal_files.flat]):
    msgs.error('Missing cal files.')

# Create arrays to hold JWST spec2, but only load the files when they're needed
msa_data = np.empty((ndet, nexp), dtype=object)
flat_data = np.empty((ndet, nexp), dtype=object)
cal_data = np.empty((ndet, nexp), dtype=object)

# Connect to the ginga window
if show:
    msgs.info('Connecting to ginga')
    display.connect_to_ginga(raise_err=True, allow_new=True)

# Loop over all slits, create a list of spec2d objects and run 2d coadd
for ii, source in enumerate(sources):
    spec2d_list = []
    offsets_pixels = []
    if show:
        msgs.info('Clearing ginga')
        display.clear_all()
    # Loop over detectors, and exposures
    for idet in range(ndet):
        for iexp in range(nexp):
            if cal_data[idet,iexp] is None:
                msgs.info(f'Loading {cal_files[idet,iexp]}')
                cal_data[idet,iexp] = datamodels.open(cal_files[idet,iexp])
            indx = np.where([s.source_name == source for s in cal_data[idet,iexp].slits])[0]
            if indx.size == 0:
                msgs.info(f'Source {source} is not observed by {basenames[idet,iexp]}.')
                continue

            # Load the rest of the needed JWST spec2 data
            if msa_data[idet,iexp] is None:
                msgs.info(f'Loading {msa_files[idet,iexp]}')
                msa_data[idet,iexp] = datamodels.open(msa_files[idet,iexp])
            if flat_data[idet,iexp] is None:
                msgs.info(f'Loading {flat_files[idet,iexp]}')
                flat_data[idet,iexp] = datamodels.open(flat_files[idet,iexp])

            jj = indx[0]

            # Grab the calibration subimages for the first exposure, use for all
            # subsequent exposures since sometimes (because of a bug) the sizes
            # of these calibration (and science) sub-images can change.

#            if iexp == 0:
            slit_slice, slit_left, slit_righ, slit_left_orig, slit_righ_orig, spec_vals_orig, \
                    src_trace_ra, src_trace_dec, dq_sub, ra, dec, finitemask, waveimg, tilts, \
                    flatfield, pathloss, barshadow, photom_conversion, calwebb \
                        = jwst_extract_subimgs(cal_data[idet,iexp].slits[jj],
                                               flat_data[idet,iexp].slits[jj])

            # Process the science image for this exposure
            t_eff = cal_data[idet,iexp].meta.exposure.effective_exposure_time
            science, sciivar, gpm, base_var, count_scale \
                    = jwst_proc(msa_data[idet,iexp], t_eff, slit_slice, finitemask, flatfield, pathloss,
                                barshadow, photom_conversion, noise_floor=noise_floor, kludge_err=kludge_err,
                                ronoise=det_par[idet].ronoise)

            # If there are not good pixels continue
            if not np.any(gpm):
                msgs.info(f'No good pixels for {source} in exposure {basenames[idet,iexp]}')
                continue

            # Instantiate
            rn2img = np.full_like(science, det_par[idet].ronoise[0]**2)
            # TODO: If we have base_var, do we need rn2img?
            sciImg = PypeItImage(image=science, ivar=sciivar, base_var=base_var,
                                 img_scale=count_scale, rn2img=rn2img, detector=det_par[idet],
                                 bpm=np.logical_not(gpm))

            if diff_redux:
                ibkg = bkg_indices[iexp]
                ibkg = ibkg if isinstance(ibkg, list) else [ibkg]
                bkg = [None]*len(ibkg)
                bkgivar = [None]*len(ibkg)
                bkg_gpm = [None]*len(ibkg)
                bkg_base_var = [None]*len(ibkg)
                bkg_count_scale = [None]*len(ibkg)
                for i, _ibkg in enumerate(ibkg):
                    if cal_data[idet,_ibkg] is None:
                        msgs.info(f'Loading {cal_files[idet,_ibkg]}')
                        cal_data[idet,_ibkg] = datamodels.open(cal_files[idet,_ibkg])
                    if msa_data[idet,_ibkg] is None:
                        msgs.info(f'Loading {msa_files[idet,_ibkg]}')
                        msa_data[idet,_ibkg] = datamodels.open(msa_files[idet,_ibkg])
                    # Process the background image using the same calibrations as the science
                    t_eff = cal_data[idet,_ibkg].meta.exposure.effective_exposure_time
                    print(f'procing {i}')
                    bkg[i], bkgivar[i], bkg_gpm[i], bkg_base_var[i], bkg_count_scale[i] \
                            = jwst_proc(msa_data[idet,_ibkg], t_eff, slit_slice, finitemask, flatfield, pathloss,
                                        barshadow, photom_conversion, noise_floor=noise_floor, kludge_err=kludge_err,
                                        ronoise=det_par[idet].ronoise)
                if len(ibkg) > 1:
                    _bkg = np.array(bkg)
                    if _bkg.ndim != 3:
                        msgs.error('bad shape')
                    weights = np.ones(_bkg.shape[0], dtype=float)/_bkg.shape[0]
                    img_list_out, var_list_out, gpm, nstack \
                        = combine.weighted_combine(weights,
                                               [_bkg, np.array(bkg_count_scale)],
                                               [inverse(np.array(bkgivar)), np.array(bkg_base_var)],
                                               np.array(bkg_gpm), sigma_clip=True,
                                               sigma_clip_stack=_bkg, sigrej=3., maxiters=5)
                    bkg, bkg_count_scale = img_list_out
                    bkgivar, bkg_base_var = var_list_out
                    bkgivar = inverse(bkgivar)
                    bkgivar[gpm] /= bkg_count_scale[gpm]**2
                    bkg_base_var[gpm] /= bkg_count_scale[gpm]**2
                    bkg_gpm = gpm
                else:
                    bkg = bkg[0]
                    bkgivar = bkgivar[0]
                    bkg_gpm = bkg_gpm[0]
                    bkg_base_var = bkg_base_var[0]
                    bkg_count_scale = bkg_count_scale[0]

                # If there are not good pixels continue
                if not np.any(bkg_gpm):
                    msgs.info(f'No good pixels for {source} in background exposure '
                              f'{basenames[idet,_ibkg]}')
                    continue

                # Instantiate
                rn2img = np.full_like(bkg, det_par[idet].ronoise[0]**2)
                bkgImg = PypeItImage(image=bkg, ivar=bkgivar, base_var=bkg_base_var,
                                     img_scale=bkg_count_scale, rn2img=rn2img,
                                     detector=det_par[idet], bpm=np.logical_not(bkg_gpm))

                # Perform the difference imaging, propagate the error and masking
                sciImg = sciImg.sub(bkgImg)

            # Create the SlitTraceSet
            nspec, nspat = waveimg.shape
            slits = slittrace.SlitTraceSet(slit_left, slit_righ, spectrograph.pypeline,
                                           detname=det_par[idet].name, nspat=nspat,
                                           PYP_SPEC=spectrograph.name)
            slits.maskdef_id = np.array([ii])

            # Construct the Spec2DObj with the positive image
            spec2DObj = spec2dobj.Spec2DObj(sciimg=sciImg.image,
                                            ivarraw=sciImg.ivar,
                                            skymodel=np.zeros_like(sciImg.image),
                                            objmodel=np.zeros_like(sciImg.image),
                                            ivarmodel=sciImg.ivar,
                                            scaleimg=None,
                                            waveimg=waveimg,
                                            bpmmask=sciImg.fullmask,
                                            detector=sciImg.detector,
                                            sci_spat_flexure=None,
                                            sci_spec_flexure=None,
                                            vel_corr=None,
                                            vel_type=None,
                                            tilts=tilts,
                                            slits=slits,
                                            wavesol=None,
                                            maskdef_designtab=None)

            spec2d_list.append(spec2DObj)
            offsets_pixels.append(dither_offsets_pixels[idet,iexp].tolist())

            if show: # and iexp == 0:
                sci_rate = msa_data[idet,iexp]
                sci_data = np.array(sci_rate.data.T, dtype=float)

                viewer_raw, ch_raw \
                        = display.show_image(sci_data, cuts=get_cuts(sci_data),
                                             chname=f'raw_rate_iexp_{iexp}_{detectors[idet]}')
                display.show_slits(viewer_raw, ch_raw, slit_left_orig, slit_righ_orig,
                                   spec_vals=spec_vals_orig, pstep=1, slit_ids=[source])
                if diff_redux:
                    ibkg = bkg_indices[iexp]
                    ibkg = ibkg if isinstance(ibkg, list) else [ibkg]
                    bkg_rate = msa_data[idet,ibkg[0]]
                    bkg_data = np.array(bkg_rate.data.T, dtype=float)
                    bg_sub = science - bkg
                    viewer_sci, ch_sci \
                            = display.show_image(sci_data-bkg_data,
                                                 cuts=get_cuts(sci_data-bkg_data), 
                                                 chname=f'diff_rate_iexp_{iexp}_{detectors[idet]}')
                    display.show_slits(viewer_sci, ch_sci, slit_left_orig, slit_righ_orig,
                                       spec_vals=spec_vals_orig, pstep=1, slit_ids=[source])

                # Show individual calibrations
                if idet == 0:
                    viewer_data, ch_data \
                            = display.show_image(calwebb, waveimg=waveimg, cuts=get_cuts(calwebb),
                                                 chname=f'calwebb_{iexp}_{detectors[idet]}')
                    display.show_trace(viewer_data, ch_data, src_trace_ra,
                                       f'trace-RA_{detectors[idet]}', color='#f0e442', pstep=1)
                    display.show_trace(viewer_data, ch_data, src_trace_dec,
                                       f'trace-DEC_{detectors[idet]}', color='#f0e442', pstep=1)
                    viewer_pypeit, ch_pypeit \
                            = display.show_image(sciImg.image, waveimg=waveimg,
                                                 cuts=get_cuts(sciImg.image),
                                                 chname=f'pypeit_{detectors[idet]}')
                    display.show_slits(viewer_pypeit, ch_pypeit, slit_left, slit_righ, pstep=1,
                                       slit_ids=[source])
                    viewer_wave, ch_wave \
                            = display.show_image(waveimg, chname=f'wave_{detectors[idet]}')
                    viewer_tilts, ch_tilts = display.show_image(tilts, waveimg=waveimg,
                                                                chname=f'tilts_{detectors[idet]}')
                    viewer_flat, ch_flat = display.show_image(flatfield, waveimg=waveimg,
                                                              chname=f'flat_{detectors[idet]}')
                    viewer_path, ch_path \
                            = display.show_image(pathloss, waveimg=waveimg,
                                                 chname=f'pathloss_{detectors[idet]}')
                    viewer_bar, ch_bar \
                            = display.show_image(barshadow, waveimg=waveimg, cuts=(0.0, 1.0),
                                                 chname=f'barshadow_{detectors[idet]}')

    if len(spec2d_list) == 0:
        msgs.info(f'No data found for {source}.')
        continue

    basename = f'{out_root}_{source}'

    # TODO: Implement an option to extract everything onto the same
    # wavelength grid optimized to match the coverage of JWST. For example,
    # for the prism things are quite nonlinear.

    # TODO: Not sure what to do with the detector container here
    # Instantiate Coadd2d
    # COADD2D doesn't yet support input offsets?
    coAdd = coadd2d.CoAdd2D.get_instance(spec2d_list, spectrograph, par, det=det_par[0].det,
                                         offsets=offsets_pixels, weights='uniform',
                                         spec_samp_fact=spec_samp_fact,
                                         spat_samp_fact=spat_samp_fact,
                                         #bkg_redux=diff_redux, debug=show)
                                         bkg_redux=False, debug=show)

    coadd_dict_list = coAdd.coadd(only_slits=None, interp_dspat=True)

    # Create the pseudo images
    pseudo_dict = coAdd.create_pseudo_image(coadd_dict_list)

    sciimg_coadd, sciivar_coadd, skymodel_coadd, objmodel_coadd, ivarmodel_coadd, \
        outmask_coadd, sobjs_coadd, detector_coadd, slits_coadd, tilts_coadd, waveimg_coadd \
            = coAdd.reduce(pseudo_dict, show=True, clear_ginga=False, show_peaks=show,
                           basename=basename)

    # Tack on detector (similarly to pypeit.extract_one)
    for sobj in sobjs_coadd:
        sobj.DETECTOR = det_par[0]

    # TODO: not currently using counts_scale and base_var. Need to rework
    # coadd2d to operate on sciimgs

    # Construct the Spec2DObj with the positive image
    spec2DObj_coadd = spec2dobj.Spec2DObj(sciimg=sciimg_coadd,
                                          ivarraw=sciivar_coadd,
                                          skymodel=skymodel_coadd,
                                          objmodel=objmodel_coadd,
                                          ivarmodel=ivarmodel_coadd,
                                          scaleimg=None,
                                          waveimg=waveimg_coadd,
                                          bpmmask=outmask_coadd,
                                          detector=det_par[0],
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
        # Plot the 2d
        nobj = len(sobjs_coadd)
        left, right, mask = spec2DObj_coadd.slits.select_edges()
        slid_IDs = spec2DObj_coadd.slits.slitord_id
        image = spec2DObj_coadd.sciimg  # Processed science image
        gpm = spec2DObj_coadd.bpmmask.mask == 0
        mean, med, sigma = sigma_clipped_stats(image[gpm], sigma_lower=5.0, sigma_upper=5.0)
        cut_min = mean - 4.0 * sigma
        cut_max = mean + 4.0 * sigma
        chname_sci = f'slit-science-{source}'
        # Clear all channels at the beginning
        viewer, ch_sci = display.show_image(image, chname=chname_sci,
                                            waveimg=spec2DObj_coadd.waveimg, clear=False,
                                            cuts=(cut_min, cut_max))

        if nobj > 0:
            show_trace(sobjs_coadd, 'DET01', viewer, ch_sci)
        display.show_slits(viewer, ch_sci, left, right, slit_ids=slid_IDs)

        image = np.sqrt(inverse(spec2DObj_coadd.ivarmodel))  # Processed science image
        mean, med, sigma = sigma_clipped_stats(image[gpm], sigma_lower=5.0, sigma_upper=5.0)
        cut_min = mean - 1.0 * sigma
        cut_max = mean + 4.0 * sigma
        chname_sig = f'slit-sigma-{source}'
        # Clear all channels at the beginning
        viewer, ch_sig = display.show_image(image, chname=chname_sig,
                                            waveimg=spec2DObj_coadd.waveimg, clear=False,
                                            cuts=(cut_min, cut_max))

        if nobj > 0:
            show_trace(sobjs_coadd, 'DET01', viewer, ch_sig)
        display.show_slits(viewer, ch_sig, left, right, slit_ids=slid_IDs)

        channel_names = [chname_sci, chname_sig]

        # After displaying all the images sync up the images with WCS_MATCH
        shell = viewer.shell()
        shell.start_global_plugin('WCSMatch')
        shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [channel_names[-1]],
                                        {})
        if not diff_redux:
            slitmask_coadd = slits_coadd.slit_img(initial=False, flexure=None, exclude_flag=None)
            slitord_id = slits_coadd.slitord_id[0]
            thismask = slitmask_coadd == slitord_id

            gpm_extract = spec2DObj_coadd.bpmmask.mask == 0
            # Make a plot of the residuals for a random slit
            chi = (spec2DObj_coadd.sciimg - spec2DObj_coadd.objmodel - spec2DObj_coadd.skymodel) \
                    * np.sqrt(spec2DObj_coadd.ivarmodel) * gpm_extract.astype(float)

            maskchi = thismask & gpm_extract

            n_bins = 50
            sig_range = 7.0
            binsize = 2.0 * sig_range / n_bins
            bins_histo = -sig_range + np.arange(n_bins) * binsize + binsize / 2.0

            xvals = np.arange(-10.0, 10, 0.02)
            gauss = scipy.stats.norm(loc=0.0, scale=1.0)
            gauss_corr = scipy.stats.norm(loc=0.0, scale=1.0)

            sigma_corr, maskchi = coadd.renormalize_errors(chi, maskchi, max_corr=20.0,
                                                           title='jwst_sigma_corr', debug=True)


        if nobj > 0:
            wv_gpm = sobjs_coadd[0].BOX_WAVE > 1.0
            plt.plot(sobjs_coadd[0].BOX_WAVE[wv_gpm],
                     sobjs_coadd[0].BOX_COUNTS[wv_gpm]*sobjs_coadd[0].BOX_MASK[wv_gpm],
            color='black', drawstyle='steps-mid', label='Counts')
            plt.plot(sobjs_coadd[0].BOX_WAVE[wv_gpm],
                     sobjs_coadd[0].BOX_COUNTS_SIG[wv_gpm]*sobjs_coadd[0].BOX_MASK[wv_gpm],
            color='red', drawstyle='steps-mid', label='Counts Error')
            plt.legend()
            plt.show()

    # container for specobjs and Spec2d
    all_specobjs = specobjs.SpecObjs()

    all_spec2d = spec2dobj.AllSpec2DObj()
    # set some meta
    all_spec2d['meta']['bkg_redux'] = diff_redux #False
    all_spec2d['meta']['find_negative'] = False

    # fill the specobjs container
    all_specobjs.add_sobj(sobjs_coadd)
    all_spec2d[det_par[0].name] = spec2DObj_coadd

    # Write out specobjs
    # Build header for spec2d
    # TODO: Get header some other way...
    head2d = fits.getheader(cal_files[0,0])
    subheader = spectrograph.subheader_for_spec(head2d, head2d, allow_missing=True)
    if all_specobjs.nobj > 0:
        outfile1d = pypeit_sci_dir / f'spec1d_{basename}.fits'
        all_specobjs.write_to_fits(subheader, str(outfile1d))

        # Info
        outfiletxt = pypeit_sci_dir / f'spec1d_{basename}.txt'
        all_specobjs.write_info(str(outfiletxt), spectrograph.pypeline)

    # Build header for spec2d
    outfile2d = pypeit_sci_dir / f'spec2d_{basename}.fits'
    # TODO For the moment hack so that we can write this out
    pri_hdr = all_spec2d.build_primary_hdr(head2d, spectrograph, subheader=subheader,
                                           redux_path=None, master_key_dict=None, master_dir=None)
    # Write spec2d
    all_spec2d.write_to_fits(str(outfile2d), pri_hdr=pri_hdr, overwrite=True)



