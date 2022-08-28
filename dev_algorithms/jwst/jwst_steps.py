
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


# 1) Assign WCS
step = AssignWcsStep()
step.save_results = True
step.output_dir = output_dir
result = step(scifile)


# The output file has the suffix _assignwcsstep appended:
awcs_output_file = os.path.join(output_dir, basename +  'assignwcsstep.fits')
print(awcs_output_file)

# load the output into a data model container
awcs = datamodels.open(awcs_output_file)
awcs

# get the WCS information populated by the algorithms of the assign_wcs step
wcsobj = awcs.meta.wcs

# list the frames of reference available for transformation
wcsobj.available_frames

# examples of the units for some of these frames
print(wcsobj.detector.unit) # xy pixel indices
print(wcsobj.slit_frame.unit)  # relative xy position in the slit aperture
print(wcsobj.msa_frame.unit)  # absolute xy position in the slit aperture
print(wcsobj.world.unit)  # RA, Dec sky position

# 2) Subtract the background
bg_subtract = False
if bg_subtract:
    backgrounds = [bkgfile1, bkgfile2]
    # these exposures are for the middle spectral dither positions at primary dither positions 1 and 3

    step = BackgroundStep()
    step.save_results = True
    step.output_dir = output_dir
    result = step(awcs_output_file, backgrounds)

    # The output file has the suffix _backgroundstep appended:
    bsub_output_file = os.path.join(output_dir, basename + 'backgroundstep.fits')

    # load the output into a data model container
    bsub = datamodels.open(bsub_output_file)
    bsub_img = np.array(bsub.data, dtype=float)
    # plot the image
    viewer_diff, ch_diff = display.show_image(bsub_img.T, cuts=get_cuts(bsub_img.T), chname='bg_subtracted', wcs_match=True)


# 3) Extract 2d
ex2d_input = awcs_output_file if not bg_subtract else bsub_output_file
step = Extract2dStep()
step.save_results = True
step.output_dir = output_dir
result = step(ex2d_input)

# The output file has the suffix _extrac2dstep appended:
e2d_output_file = os.path.join(output_dir, basename + 'extract2dstep.fits')
# load the output into a data model container
e2d = datamodels.open(e2d_output_file)
nslits = len(e2d.slits)
show_regions = True
if show_regions:
    for islit in np.arange(nslits):
        # Read in data print out slit name
        slit_name = e2d.slits[islit].name
        print(slit_name)
        calsci = np.array(e2d.slits[islit].data, dtype=float)  # contains the pixel data from the cal file (SCI extension)
        nspat, nspec = calsci.shape

        ########################
        # Plot the image segment being used for each slit
        xlo = e2d.slits[islit].xstart - 1
        xhi = xlo + e2d.slits[islit].xsize
        ylo = e2d.slits[islit].ystart - 1
        yhi = ylo + e2d.slits[islit].ysize
        # This is the segment of the 2d image
        slit_slice = np.s_[ylo: yhi, xlo: xhi]
        #xvals = xlo + np.arange(xhi - xlo)
        #yvals = ylo + np.arange(yhi - ylo)
        slit_left = np.full(nspec, ylo)
        slit_righ = np.full(nspec, yhi)
        spec_val = xlo + np.arange(xhi - xlo)
        display.show_slits(viewer_sci, ch_sci, slit_left, slit_righ, spec_vals = spec_val, pstep=1, slit_ids=np.array([int(slit_name)]))


# 4) Assign source type.
# TODO: So far this influences the wavecorr step (slit centering), the pathloss step (point vs uniform),
# as well as whether the barshadow correction is applied.
step = SourceTypeStep(source_type='EXTENDED')
step.save_results = True
step.output_dir = output_dir
result = step(e2d_output_file)

# The output file has the suffix _sourcetypestep appended:
stype_output_file = os.path.join(output_dir, basename + 'sourcetypestep.fits')
# load the output into a data model container
stype = datamodels.open(stype_output_file)
stype_list = [slit.source_type for slit in stype.slits]
assert np.all(np.array(stype_list) == 'EXTENDED')

# 5) Apply wevecorr correction
# According to calwebbb, the wavelength values per pixel need to be corrected when a point source is not located at
# the center of a slit aperture in the dispersion direction.
# TODO: I don't understand this, this correction needs to be applied irrespective of whether the source is point or
#  extended, although it is clearly more complicated if the the source is extended and has a light centroid favoring
#  one side of the slit realtive to the other. Maybe the thing to do here is perform this step after setting all sources
#  to point, but then change their source type so that the barshadow etc. correction will be performed. Anyway, I'm running
#  this below but I think it does nothing

step = WavecorrStep() # TODO: Can you specify the source type manually here?
step.save_results = True
step.output_dir = output_dir
result = step(stype_output_file)

# The output file has the suffix _wavecorr appended:
wavecorr_output_file = os.path.join(output_dir, basename + 'sourcetypestep.fits')
#wavecorr = datamodels.open(wavecorr_output_file)

# 6) Apply the flat field correction
step = FlatFieldStep()
step.save_interpolated_flat = True  # this will save the on-the-fly flat field correction values as an image
step.save_results = True
step.output_dir = output_dir
result = step(wavecorr_output_file)

# The output file has the suffix _flatfieldstep appended:
flat_output_file = os.path.join(output_dir, basename + 'flatfieldstep.fits')
# The optional saved flat correction image has the suffix _interpolatedflat appended:
intflat_output_file = os.path.join(output_dir, basename + 'interpolatedflat.fits')
# load the output into a data model container
#flat_data = datamodels.open(flat_output_file)
intflat = datamodels.open(intflat_output_file)


# 7) Apply the pathloss correction
# The pathloss correction scales the 2D spectrum as a function of wavelength to account for geometric and diffraction
# losses of a non-centered point source incurred by the slit aperture. This is only relative to a centered source, as
# the absolute pathloss in that case has already been corrected by the F flat component of the flat field correction.
# More information about this can be found in this ESA technical note (just type into Google): ESA-JWST-SCI-NRS-TN-2016-004
# TODO: I don't undersatnd how this interacts with source type, note you can pass source type as an argument here.
step = PathLossStep()
step.save_results = True
step.output_dir = output_dir
result = step(flat_output_file)
# The output file has the suffix _flatfieldstep appended:
ploss_output_file = os.path.join(output_dir, basename + 'pathlossstep.fits')

# 8) Apply the barshadow correction
step = BarShadowStep() # You can force this to do a specifc source_type
step.save_results = True
step.output_dir = output_dir
result = step(ploss_output_file)
# The output file has the suffix _barshadowstep appended:
bshadow_output_file = os.path.join(output_dir, basename + 'barshadowstep.fits')

# 9) Apply the photometric correction converting the image to physical units
# The photom step applies a scalar conversion factor to convert to physical units.  Once on-orbit observations of
# spectrophotometric standards are obtained, this may also include a wavelength-dependent vector, if necessary, to
# account for any small discrepancies found in the throughput corrections.
# TODO: You can specify a source type here via the source_type argument
step = PhotomStep()
step.save_results = True
step.output_dir = output_dir
result = step(bshadow_output_file)
# The output file has the suffix _photomstep appended:
phot_output_file = os.path.join(output_dir, basename + 'photomstep.fits')
# load the output into a data model container
phot = datamodels.open(phot_output_file)

# what are the flux calibration-related keywords?
phot.find_fits_keyword('PHOTUJA2')
#print('PHOTMJSR=', phot.slits[0].meta.photometry.conversion_megajanskys)
# note that despite the keyword name, when SRCTYPE=POINT, the units are actually MJy per pixel
# the pre-launch reference file has a placeholder value, so the output fluxes will appear unphysically large
print('PHOTUJA2=', phot.slits[0].meta.photometry.conversion_microjanskys)
# this is only applicable for extended sources
# TODO I cannot figure out how to plot the sensitivity function that is applied here. The code is confusing, but only
# the units appear to be recorded in the MultiSlitModel object

final2d = datamodels.open(phot_output_file)
show_cutouts = True
only_slit = 10
slit_indx = np.arange(nslits) if only_slit is None else [only_slit]

if show_cutouts:
    for islit in slit_indx:
        # Read in data print out slit name
        slit_name = phot.slits[islit].name
        print(slit_name)
        calsci = np.array(final2d.slits[islit].data, dtype=float)  # contains the pixel data from the cal file (SCI extension)
        nspat, nspec = calsci.shape

        # get the source RA and Dec coordinates from the metadata (also located in the header of the fits SCI extension)
        source_ra = final2d.slits[islit].meta.target.ra
        source_dec = final2d.slits[islit].meta.target.dec
        print('catalog RA,DEC:', source_ra, source_dec)
        # determine the wavelength scale of the cal data for plotting purposes
        # get the data model WCS object. This example is from the fixed slit notebook
        slit_wcs = final2d.slits[islit].meta.wcs
        x, y = wcstools.grid_from_bounding_box(slit_wcs.bounding_box, step=(1, 1))
        calra, caldec, calwave = slit_wcs(x, y)

        ## Old way from fixed slit notebook
        #y1, x1 = np.mgrid[:nspat,:nspec]  # grid of pixel x,y indices
        #det2sky = slit_wcs.get_transform('detector','world')  # the coordinate transform from detector space (pixels) to sky (RA, DEC in degrees)
        #calra, caldec, calwave = det2sky(x1, y1)  # RA, Dec, wavelength (microns) for each pixel
        cal_spec = np.arange(nspec)  # spectral position
        cal_spat = np.arange(nspat)  # spatial position
        cal_src_from_ra_spat = np.zeros(nspec) # Array to hold the source_RA as a function of spectral position
        cal_src_from_dec_spat = np.zeros(nspec) # Array to hold the source_DEC as a function of spectral position
        for ispec in range(nspec):
            ra_vs_spat = calra[:, ispec] #
            # Interpolate y-pixel as a functio of RA onto the source RA
            cal_src_from_ra_spat[ispec] = np.interp(source_ra, ra_vs_spat[np.isfinite(ra_vs_spat)], cal_spat[np.isfinite(ra_vs_spat)])
            dec_vs_spat = caldec[:, ispec]
            cal_src_from_dec_spat[ispec] = np.interp(source_dec, dec_vs_spat[np.isfinite(dec_vs_spat)], cal_spat[np.isfinite(dec_vs_spat)])

        # Now transpose everything to PypeIt convention for viewing.

        # plot the unrectified calibrated 2D spectrum
        waveimg = np.array(final2d.slits[islit].wavelength.T,dtype=float)
        pathloss = np.array(final2d.slits[islit].pathloss_uniform.T,dtype=float) \
            if final2d.slits[islit].source_type == 'EXTENDED' else np.array(final2d.slits[islit].pathloss_point.T,dtype=float)
        flat = np.array(intflat.slits[islit].data.T,dtype=float)
        barshadow = np.array(final2d.slits[islit].barshadow.T,dtype=float)
        viewer_data, ch_data = display.show_image(calsci.T, waveimg = waveimg, cuts = get_cuts(calsci.T), chname=slit_name + '_data')
        viewer_wave, ch_wave = display.show_image(waveimg, waveimg=waveimg, chname=slit_name + '_wave')
        viewer_flat, ch_flat = display.show_image(flat, waveimg=waveimg, chname=slit_name + '_flat')
        viewer_path, ch_path = display.show_image(pathloss, waveimg=waveimg, chname=slit_name + '_pathloss')
        viewer_path, ch_path = display.show_image(barshadow, waveimg=waveimg, chname=slit_name + '_barshadow')
        display.show_trace(viewer_data, ch_data, cal_src_from_ra_spat, 'RA', color='#f0e442')
        display.show_trace(viewer_data, ch_data, cal_src_from_dec_spat, 'DEC', color='#f0e442')

