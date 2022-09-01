
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from pypeit.display import display
from pypeit.core import fitting
from gwcs import wcstools
from matplotlib import pyplot as plt
from jwst import datamodels


from IPython import embed

def compute_diff(scifile, bkgfile1, bkgfile2):
    sci_rate = datamodels.open(scifile)
    bkg1_rate = datamodels.open(bkgfile1)
    bkg2_rate = datamodels.open(bkgfile2)

    sci = sci_rate.data
    diff = sci_rate.data - (bkg1_rate.data + bkg2_rate.data)/2.0

    return sci, diff

def get_cuts(image):
    mean, med, sigma = sigma_clipped_stats(image, sigma_lower=5.0, sigma_upper=5.0)
    cut_min = mean - 1.0 * sigma
    cut_max = mean + 4.0 * sigma
    return (cut_min, cut_max)

def fit_slit(thismask, left_or_right, polyorder=2, function='legendre', debug=False):

    slit_width = np.sum(thismask, axis=1)
    med_slit_width = np.median(slit_width[slit_width > 0])
    nspec, nspat = thismask.shape
    spec_vec = np.arange(nspec, dtype=float)
    spat_vec = np.arange(nspat, dtype=float)
    spat_img, spec_img = np.meshgrid(spat_vec, spec_vec)

    dummy_spat_img = spat_img.copy()
    bad_value = +np.inf if 'left' in left_or_right else -np.inf
    dummy_spat_img[np.logical_not(thismask)] = bad_value
    slit_mask = np.min(dummy_spat_img, axis=1) if 'left' in left_or_right else np.max(dummy_spat_img, axis=1)
    good_for_slit = (slit_width > 0.5 * med_slit_width) & (slit_mask != bad_value)
    bad_for_slit = np.logical_not(good_for_slit)

    pypeitFit = fitting.robust_fit(spec_vec[good_for_slit], slit_mask[good_for_slit], polyorder, function=function,
                                   maxiter=25,
                                   lower=3.0, upper=2.0, maxrej=1, sticky=True,
                                   minx=0.0, maxx=float(nspec - 1))
    slit = pypeitFit.eval(spec_vec)
    if debug:
        plt.plot(spec_vec[good_for_slit], slit_mask[good_for_slit], 'k.')
        plt.plot(spec_vec[bad_for_slit], slit_mask[bad_for_slit], 'r.')
        plt.plot(spec_vec, slit, 'b')
        plt.show()


    return slit

def get_jwst_slits(thismask, polyorder=2, function='legendre', debug=False):

    slit_left = fit_slit(thismask, 'left', polyorder=polyorder, function=function, debug=debug)
    slit_righ = fit_slit(thismask, 'righ', polyorder=polyorder, function=function, debug=debug)
    return slit_left, slit_righ



def show_2dspec(rawscience, final2d, islit, intflat=None, clear=True, emb=None):


    # Read in data print out slit name
    slit_name = final2d.slits[islit].name
    print('Slit={:s}'.format(slit_name))
    calsci = np.array(final2d.slits[islit].data, dtype=float)  # contains the pixel data from the cal file (SCI extension)
    nspat, nspec = calsci.shape


    ########################
    # Plot the image segment being used for each slit
    xlo = final2d.slits[islit].xstart - 1
    xhi = xlo + final2d.slits[islit].xsize
    ylo = final2d.slits[islit].ystart - 1
    yhi = ylo + final2d.slits[islit].ysize
    # This is the segment of the 2d image
    slit_slice = np.s_[ylo: yhi, xlo: xhi]
    # xvals = xlo + np.arange(xhi - xlo)
    # yvals = ylo + np.arange(yhi - ylo)
    slit_left = np.full(nspec, ylo)
    slit_righ = np.full(nspec, yhi)
    spec_val = xlo + np.arange(xhi - xlo)
    viewer_sci, ch_sci = display.show_image(rawscience.T, cuts=get_cuts(rawscience), chname='raw', clear=clear)
    display.show_slits(viewer_sci, ch_sci, slit_left, slit_righ, spec_vals=spec_val, pstep=1,
                       slit_ids=np.array([int(slit_name)]))

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
    waveimg = calwave.T if (final2d.slits[islit].wavelength.shape == (0,0)) \
            else np.array(final2d.slits[islit].wavelength.T,dtype=float)
    pathloss = np.array(final2d.slits[islit].pathloss_uniform.T,dtype=float) \
        if final2d.slits[islit].source_type == 'EXTENDED' else np.array(final2d.slits[islit].pathloss_point.T,dtype=float)
    barshadow = np.array(final2d.slits[islit].barshadow.T,dtype=float)
    viewer_data, ch_data = display.show_image(calsci.T, waveimg = waveimg, cuts = get_cuts(calsci.T), chname=slit_name + '_data')
    viewer_wave, ch_wave = display.show_image(waveimg, waveimg=waveimg, chname=slit_name + '_wave')
    viewer_ra, ch_ra = display.show_image(calra.T, waveimg=waveimg, chname=slit_name + '_RA')
    if intflat is not None:
        flat = np.array(intflat.slits[islit].data.T,dtype=float)
        viewer_flat, ch_flat = display.show_image(flat, waveimg=waveimg, chname=slit_name + '_flat')
    display.show_trace(viewer_data, ch_data, cal_src_from_ra_spat, trc_name='RA', pstep=1, color='#f0e442')
    display.show_trace(viewer_data, ch_data, cal_src_from_dec_spat, trc_name='DEC', pstep=1,  color='#f0e442')
    if pathloss.shape != (0,0):
        viewer_path, ch_path = display.show_image(pathloss, waveimg=waveimg, chname=slit_name + '_pathloss')
    if barshadow.shape != (0,0):
        viewer_bar, ch_bar = display.show_image(barshadow, waveimg=waveimg, chname=slit_name + '_barshadow')

    if emb:
        embed(header='Slit={:s}'.format(slit_name))


