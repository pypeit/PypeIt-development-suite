
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from pypeit.display import display
from gwcs import wcstools

def show_diff(scifile, bkgfile1, bkgfile2):
    hdu_sci = fits.open(scifile)
    sci,  sci_err = hdu_sci[1].data, hdu_sci[2].data
    hdu_bkg1 = fits.open(bkgfile1)
    bkg1, bkg1_err = hdu_bkg1[1].data, hdu_bkg1[2].data
    hdu_bkg2 = fits.open(bkgfile2)
    bkg2, bkg2_err = hdu_bkg2[1].data, hdu_bkg2[2].data

    diff = sci - (bkg1 + bkg2)/2.0

    return sci, diff

def get_cuts(image):
    mean, med, sigma = sigma_clipped_stats(image, sigma_lower=5.0, sigma_upper=5.0)
    cut_min = mean - 1.0 * sigma
    cut_max = mean + 4.0 * sigma
    return (cut_min, cut_max)


def show_2dspec(final2d, islit, intflat=None):

    # Read in data print out slit name
    slit_name = final2d.slits[islit].name
    print('Slit={:s}'.format(slit_name))
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
    barshadow = np.array(final2d.slits[islit].barshadow.T,dtype=float)
    viewer_data, ch_data = display.show_image(calsci.T, waveimg = waveimg, cuts = get_cuts(calsci.T), chname=slit_name + '_data')
    viewer_wave, ch_wave = display.show_image(waveimg, waveimg=waveimg, chname=slit_name + '_wave')
    viewer_ra, ch_ra = display.show_image(calra.T, waveimg=waveimg, chname=slit_name + '_RA')
    if intflat is not None:
        flat = np.array(intflat.slits[islit].data.T,dtype=float)
        viewer_flat, ch_flat = display.show_image(flat, waveimg=waveimg, chname=slit_name + '_flat')
    viewer_path, ch_path = display.show_image(pathloss, waveimg=waveimg, chname=slit_name + '_pathloss')
    viewer_path, ch_path = display.show_image(barshadow, waveimg=waveimg, chname=slit_name + '_barshadow')
    display.show_trace(viewer_data, ch_data, cal_src_from_ra_spat, 'RA', color='#f0e442')
    display.show_trace(viewer_data, ch_data, cal_src_from_dec_spat, 'DEC', color='#f0e442')
