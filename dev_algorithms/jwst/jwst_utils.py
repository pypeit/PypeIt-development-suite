
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from pypeit.display import display
from pypeit.core import fitting
from gwcs import wcstools
from matplotlib import pyplot as plt
from jwst import datamodels
from pypeit.utils import inverse
DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']


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

def jwst_get_slits(thismask, polyorder=2, function='legendre', debug=False):

    slit_left = fit_slit(thismask, 'left', polyorder=polyorder, function=function, debug=debug)
    slit_righ = fit_slit(thismask, 'righ', polyorder=polyorder, function=function, debug=debug)
    return slit_left, slit_righ

def jwst_proc(slit, intflat_slit=None, kludge_err=1.0):


    # Try to reverse engineer all the things they multiply into the data
    slit_name = slit.name
    pathloss = np.array(slit.pathloss_uniform.T, dtype=float) if slit.source_type == 'EXTENDED' else \
        np.array(slit.pathloss_point.T, dtype=float)
    flat = np.array(intflat_slit.data.T, dtype=float)
    # flat = np.ones_like(pathloss)
    barshadow = np.array(slit.barshadow.T, dtype=float)
    photom_conversion = slit.meta.photometry.conversion_megajanskys

    t_eff = slit.meta.exposure.effective_exposure_time
    # TODO I don't know how the t_eff quantity is defined. Better would be some proxy for the exposure time per pixel
    # This is the conversion between final2d and e2d, i.e. final2d = jwst_scale*e2d
    jwst_scale = photom_conversion / flat / pathloss / barshadow
    # The science data is divided by (flat*pathloss*barshadow) and then multiplied by photom_conversion. Since
    # we work in units of counts, we divide by the photom conversion and multiply by t_eff.
    flux_to_counts = t_eff / photom_conversion # This converts s2d outputs of flux to counts.
    count_scale = inverse(flat * pathloss * barshadow)   # This is the quantity that goes into PypeIt

    science = np.array(slit.data.T, dtype=float) * flux_to_counts

    # TODO Currently the var_flat is nonsense I think and so I'm just going to use the var_poisson and var_rnoise to get
    # the noise. If this gets fixed use the line below which includes the var_flat.
    # err = kludge_err*np.array(slit.err.T, dtype=float)*flux_to_counts
    var_poisson = slit.var_poisson.T * flux_to_counts ** 2
    var_rnoise = slit.var_rnoise.T * flux_to_counts ** 2
    var = kludge_err**2*np.array(var_poisson + var_rnoise, dtype=float)
    # This needs to be multiplied by count_scale to get it into units of counts which is what pypeit requires. I checked
    # that this base_var is equal to e2d.var_rnoise if you remove the flux_to_counts factor.
    # base_var = np.array(final2d.slits[islit].var_rnoise.T, dtype=float)*flux_to_counts**2*count_scale**2
    base_var = np.array(slit.var_rnoise.T, dtype=float) * flux_to_counts ** 2
    dq = np.array(slit.dq.T, dtype=int)
    waveimg = np.array(slit.wavelength.T, dtype=float)

    gpm = np.logical_not(dq & DO_NOT_USE)
    thismask = np.isfinite(science)
    nanmask = np.logical_not(thismask)
    science[nanmask] = 0.0
    # err[nanmask] = 0.0
    var[nanmask] = 0.0
    sciivar = inverse(var) * gpm
    base_var[nanmask] = 0.0
    count_scale[nanmask] = 0.0
    # Wave nanmask is different from data nanmask
    slit_wcs = slit.meta.wcs
    x, y = wcstools.grid_from_bounding_box(slit_wcs.bounding_box, step=(1, 1))
    calra, caldec, calwave = slit_wcs(x, y)
    ra = calra.T
    nanmask_wave = np.logical_not(np.isfinite(waveimg))
    wave_min = np.min(waveimg[np.logical_not(nanmask_wave)])
    wave_max = np.max(waveimg[np.logical_not(nanmask_wave)])
    nanmask_ra = np.logical_not(np.isfinite(ra))
    ra_min = np.min(ra[np.logical_not(nanmask_ra)])
    ra_max = np.max(ra[np.logical_not(nanmask_ra)])
    waveimg[nanmask_wave] = 0.0
    ra[nanmask_ra] = 0.0


    # TODO Figure out a way to get the slit boundaries from the WCS itself instead of this kludge with the nan values
    slit_left, slit_righ = jwst_get_slits(thismask)
    # Generate some tilts and a spatial image
    tilts = np.zeros_like(waveimg)
    tilts[np.isfinite(waveimg)] = (waveimg[np.isfinite(waveimg)] - wave_min) / (wave_max - wave_min)

    # TODO Fix this spat_pix to make it increasing with pixel. For now don't use it
    nspec, nspat = science.shape
    spat_pix = (ra - ra_min) / (ra_max - ra_min) * (nspat - 1)
    spat_pix[nanmask_ra] = 0.0

    return science, sciivar, gpm, base_var, count_scale, tilts, waveimg, thismask, slit_left, slit_righ, t_eff



def jwst_show_msa(sci_rate, final2d, clear=True):

    sci_data = sci_rate.data.T
    viewer_sci, ch_sci = display.show_image(sci_data, cuts=get_cuts(sci_data), chname='raw rate', clear=clear)

    for islit, slit in enumerate(final2d.slits):
        # Read in data print out slit name
        slit_name = final2d.slits[islit].name
        calsci = np.array(final2d.slits[islit].data, dtype=float)  # contains the pixel data from the cal file (SCI extension)
        print('Slit={:s}'.format(slit_name))
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
        display.show_slits(viewer_sci, ch_sci, slit_left, slit_righ, spec_vals=spec_val, pstep=1,
                           slit_ids=np.array([int(slit_name)]))


def jwst_show_spec2(slit, intflat_slit=None, clear=True, emb=False):


    # Read in data print out slit name
    slit_name = slit.name
    print('Slit={:s}'.format(slit_name))
    calsci = np.array(slit.data, dtype=float)  # contains the pixel data from the cal file (SCI extension)
    nspat, nspec = calsci.shape


    ########################
    # Plot the image segment being used for each slit
    #xlo = final2d.slits[islit].xstart - 1
    #xhi = xlo + final2d.slits[islit].xsize
    #ylo = final2d.slits[islit].ystart - 1
    #yhi = ylo + final2d.slits[islit].ysize
    # This is the segment of the 2d image
    #slit_slice = np.s_[ylo: yhi, xlo: xhi]
    # xvals = xlo + np.arange(xhi - xlo)
    # yvals = ylo + np.arange(yhi - ylo)
    #slit_left = np.full(nspec, ylo)
    #slit_righ = np.full(nspec, yhi)
    #spec_val = xlo + np.arange(xhi - xlo)
    #viewer_sci, ch_sci = display.show_image(rawscience.T, cuts=get_cuts(rawscience), chname='raw', clear=clear)
    #display.show_slits(viewer_sci, ch_sci, slit_left, slit_righ, spec_vals=spec_val, pstep=1,
    #                   slit_ids=np.array([int(slit_name)]))

    # get the source RA and Dec coordinates from the metadata (also located in the header of the fits SCI extension)
    source_ra = slit.meta.target.ra
    source_dec = slit.meta.target.dec
    print('catalog RA,DEC:', source_ra, source_dec)
    # determine the wavelength scale of the cal data for plotting purposes
    # get the data model WCS object. This example is from the fixed slit notebook
    slit_wcs = slit.meta.wcs
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
    waveimg = calwave.T if (slit.wavelength.shape == (0,0)) else np.array(slit.wavelength.T,dtype=float)
    pathloss = np.array(slit.pathloss_uniform.T,dtype=float) if slit.source_type == 'EXTENDED' else \
        np.array(slit.pathloss_point.T,dtype=float)
    barshadow = np.array(slit.barshadow.T,dtype=float)
    viewer_data, ch_data = display.show_image(calsci.T, waveimg = waveimg, cuts = get_cuts(calsci.T),
                                              chname=slit_name + '_data', clear=clear)
    viewer_wave, ch_wave = display.show_image(waveimg, waveimg=waveimg, chname=slit_name + '_wave')
    viewer_ra, ch_ra = display.show_image(calra.T, waveimg=waveimg, chname=slit_name + '_RA')
    if intflat_slit is not None:
        flat = np.array(intflat_slit.data.T,dtype=float)
        viewer_flat, ch_flat = display.show_image(flat, waveimg=waveimg, chname=slit_name + '_flat')
    display.show_trace(viewer_data, ch_data, cal_src_from_ra_spat, trc_name='RA', pstep=1, color='#f0e442')
    display.show_trace(viewer_data, ch_data, cal_src_from_dec_spat, trc_name='DEC', pstep=1,  color='#f0e442')
    if pathloss.shape != (0,0):
        viewer_path, ch_path = display.show_image(pathloss, waveimg=waveimg, chname=slit_name + '_pathloss')
    if barshadow.shape != (0,0):
        viewer_bar, ch_bar = display.show_image(barshadow, waveimg=waveimg, chname=slit_name + '_barshadow')

    if emb:
        embed(header='Slit={:s}'.format(slit_name))


