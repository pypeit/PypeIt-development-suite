

import os
import numpy as np
from astropy.stats import sigma_clipped_stats

from IPython import embed
# set environment variables
os.environ['CRDS_PATH'] = '/Users/joe/crds_cache/jwst_pub'
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds-pub.stsci.edu'
from matplotlib import pyplot as plt
from astropy.io import fits

from pypeit.display import display
from jwst.pipeline import Spec2Pipeline
# data models
from jwst import datamodels
from jwst_utils import show_diff, get_cuts



rawpath_level2 = '/Users/joe/jwst_redux/NIRSPEC/Raw/02736_ERO_SMACS0723_G395MG235M/level_2'
output_dir = '/Users/joe/jwst_redux/NIRSPEC/redux/calwebb/output'

det = 'nrs1'
if 'nrs1' in det:
    # nrs1
    asn_file = os.path.join(rawpath_level2, 'jw02736-o007_20220712t161855_spec2_006_asn.json')
else:
    # nrs2
    asn_file =  os.path.join(rawpath_level2, 'jw02736-o007_20220712t161855_spec2_002_asn.json')


# NIRSPEC 3-point dither
scifile = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + det + '_rate.fits')
bkgfile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + det + '_rate.fits')
bkgfile2 = os.path.join(rawpath_level2, 'jw02736007001_03103_00003_' + det + '_rate.fits')
# Plot the 2d differnence image
sci, diff = show_diff(scifile, bkgfile1, bkgfile2)

viewer_diff, ch_diff = display.show_image(diff.T, cuts=get_cuts(diff), chname='diff2d')
viewer_sci,  ch_sci = display.show_image(sci.T, cuts=get_cuts(sci), chname='science')

runflag = False
if runflag:
    spec2 = Spec2Pipeline()
    spec2.save_results = True
    spec2.output_dir = output_dir
    result = spec2(asn_file)

# take a look at the results - open the level 2b files
calfile = os.path.join(output_dir, 'jw02736007001_03103_00001_' + det + '_cal.fits')
s2dfile = calfile.replace('_cal.fits', '_s2d.fits')
x1dfile = calfile.replace('_cal.fits', '_x1d.fits')
cal = datamodels.open(calfile)  # this contains the calibrated unrectified 2D spectra
s2d = datamodels.open(s2dfile)  # this contains the calibrated *rectified* 2D spectra
x1d = datamodels.open(x1dfile)  # # this contains the aperture-extracted 1D spectra
# Loop over slits to tryto visualize what the pipeline is doing
for i, slit in enumerate(cal.slits):
    # Read in data print out slit name
    print(slit.name)
    calsci = np.array(slit.data, dtype=float) # contains the pixel data from the cal file (SCI extension)
    s2dsci = np.array(s2d.slits[i].data, dtype=float) # contains the pixel data from the s2d file
    nspat, nspec = calsci.shape


    ########################
    # Plot the image segment being used for each slit
    xlo = slit.xstart - 1
    xhi = xlo + slit.xsize
    ylo = slit.ystart - 1
    yhi = ylo + slit.ysize
    # This is the segment of the 2d image
    slit_slice = np.s_[ylo: yhi, xlo: xhi]
    #xvals = xlo + np.arange(xhi - xlo)
    #yvals = ylo + np.arange(yhi - ylo)
    slit_left = np.full(nspec, ylo)
    slit_righ = np.full(nspec, yhi)
    spec_val = xlo + np.arange(xhi - xlo)
    display.show_slits(viewer_sci, ch_sci, slit_left, slit_righ, spec_vals = spec_val, pstep=1, slit_ids=np.array([int(slit.name)]))


    ########################
    # find the trace of the source coordinates

    # get the source RA and Dec coordinates from the metadata (also located in the header of the fits SCI extension)
    source_ra = cal.slits[i].meta.target.ra
    source_dec = cal.slits[i].meta.target.dec
    print('catalog RA,DEC:', source_ra, source_dec)
    # determine the wavelength scale of the cal data for plotting purposes
    # get the data model WCS object
    wcsobj = cal.slits[i].meta.wcs
    y, x = np.mgrid[:nspat,:nspec]  # grid of pixel x,y indices
    det2sky = wcsobj.get_transform('detector','world')  # the coordinate transform from detector space (pixels) to sky (RA, DEC in degrees)
    calra, caldec, calwave = det2sky(x, y)  # RA, Dec, wavelength (microns) for each pixel

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

    #s2dwaves = s2dwave[0, :]  # only need a single row of values since this is the rectified spectrum
    #xtint = np.arange(100, s2dsci.shape[1], 100)
    #xtlab = np.round(s2dwaves[xtint], 2)  # wavelength labels for the x-axis

    # get wavelength & flux from the x1d data model
    #x1dwave = x1d.spec[i].spec_table.WAVELENGTH
    #x1dflux = x1d.spec[i].spec_table.FLUX

    # Now transpose everything to PypeIt convention for viewing.

    # plot the unrectified calibrated 2D spectrum
    viewer, ch = display.show_image(calsci.T, waveimg = np.array(slit.wavelength.T,dtype=float),
                                    cuts = get_cuts(calsci.T), chname='unrectified')
    viewer_wave, ch_wave = display.show_image(np.array(slit.wavelength.T,dtype=float), chname='wavelength')
    display.show_trace(viewer, ch, cal_src_from_ra_spat, slit.name + '-RA', color='#f0e442')
    display.show_trace(viewer, ch, cal_src_from_dec_spat, slit.name + '-DEC', color='#f0e442')

    embed(header='slit={:s}'.format(slit.name))
    #input("Press Enter to continue...")

    #show_image(calsci, -0.01, 0.01, aspect=5., scale='linear', units='MJy')

    # plot the rectified 2D spectrum
    #display.show_image(s2dsci,  waveimg = np.array(s2dwave, dtype=float), cuts = get_cuts(s2dsci), chname='rectified')
    #show_image(s2dsci, -0.01, 0.01, aspect=5., scale='linear', units='MJy')
    #plt.xticks(xtint, xtlab)
    #plt.xlabel('wavelength (microns)')

    # plot the 1D extracted spectrum
    #fig = plt.figure(figsize=(19, 8))
    #plt.plot(x1dwave, x1dflux)
    #plt.xlabel('wavelength (microns)')
    #plt.ylabel('flux (Jy)')
    #plt.show()

