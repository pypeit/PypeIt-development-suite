
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
from pypeit import msgs
from pypeit.core import flat
from pypeit.core import procimg


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
                                   maxiter=25, lower=3.0, upper=3.0, maxrej=1, sticky=True, verbose=False,
                                   minx=0.0, maxx=float(nspec - 1))
    slit = pypeitFit.eval(spec_vec)
    if debug:
        plt.plot(spec_vec[good_for_slit], slit_mask[good_for_slit], 'k.')
        plt.plot(spec_vec[bad_for_slit], slit_mask[bad_for_slit], 'r.')
        plt.plot(spec_vec, slit, 'b')
        plt.show()


    return slit

def jwst_get_slits(finitemask, polyorder=5, function='legendre', debug=False):

    slit_left = fit_slit(finitemask, 'left', polyorder=polyorder, function=function, debug=debug)
    slit_righ = fit_slit(finitemask, 'righ', polyorder=polyorder, function=function, debug=debug)
    return slit_left, slit_righ


def jwst_proc(t_eff, e2d_slit, final_slit, intflat_slit, kludge_err=1.0, ronoise=5.17, saturation=65000, noise_floor=0.01,
              show=False):

    # Read in data print out slit name
    slit_name = final_slit.name

    slit_left, slit_righ, slit_left_orig, slit_righ_orig, spec_vals_orig, src_trace_ra, src_trace_dec, \
    rate, rate_var_rnoise, rate_var_poisson, rate_var_tot, dq, ra, dec, waveimg, tilts, flatfield, pathloss, \
    barshadow, photom_conversion, final = jwst_extract_subimgs(e2d_slit, final_slit, intflat_slit)

    # Now deal with the image processing
    finitemask = np.isfinite(waveimg)
    if not np.any(finitemask):
        return (None,)*21


    # Now perform the image processing
    raw_counts = rate*t_eff
    raw_var_poisson = kludge_err**2*rate_var_poisson*t_eff**2
    raw_var_rnoise = kludge_err**2*rate_var_rnoise*t_eff**2
    # Is this correct? I'm not sure I should be using their poisson variance for the noise floor
    raw_var = procimg.variance_model(raw_var_rnoise, counts = raw_var_poisson, noise_floor=noise_floor)
    # TODO This  is a hack until I can understand how to get rid of the hot pixels in the JWST variance arrays using DQ flags
    raw_gpm = (raw_var_rnoise < 2*ronoise**2) & (raw_var_poisson < saturation)
    #raw_var_poisson + raw_var_rnoise # TODO Leaving out problematic flat field term from pipeline

    # This is the conversion between final2d and e2d, i.e. final2d = jwst_scale*e2d
    # total_flat = flatfield*pathloss*barshadow
    #flux_to_counts = t_eff / photom_conversion  # This converts s2d outputs of flux to counts.
    #jwst_scale = photom_conversion/flatfield/pathloss/barshadow
    total_flat = pathloss*barshadow
    total_flat_square = np.square(total_flat)

    count_scale = inverse(total_flat)  # This is the quantity that goes into PypeIt for var modeling
    science, flat_bpm = flat.flatfield(raw_counts, total_flat)
    var_poisson, _ = flat.flatfield(raw_var_poisson, total_flat_square)
    base_var, _ = flat.flatfield(raw_var_rnoise, total_flat_square)
    var, _ = flat.flatfield(raw_var, total_flat_square)
    sciivar = inverse(var)
    dq_gpm = np.logical_not(dq & DO_NOT_USE)
    gpm = finitemask & dq_gpm & np.logical_not(flat_bpm) & (sciivar > 0.0) & raw_gpm


    if show:
        viewer_data, ch_data = display.show_image(science, waveimg = waveimg, cuts = get_cuts(science), chname='science')
        viewer_wave, ch_wave = display.show_image(waveimg, chname='wave')
        viewer_tilts, ch_tilts = display.show_image(tilts, waveimg=waveimg, chname='tilts')
        viewer_flat, ch_flat = display.show_image(flatfield, waveimg=waveimg, chname='flat')
        viewer_path, ch_path = display.show_image(pathloss, waveimg=waveimg, chname='pathloss')
        viewer_bar , ch_bar  = display.show_image(barshadow, waveimg=waveimg, chname='barshadow')
        display.show_trace(viewer_data, ch_data, src_trace_ra, 'trace-RA', color='#f0e442')
        display.show_trace(viewer_data, ch_data, src_trace_dec, 'trace-DEC', color='#f0e442')

    nanmask = np.logical_not(finitemask)
    count_scale[nanmask] = 0.0
    science[nanmask] = 0.0
    var_poisson[nanmask] = 0.0
    base_var[nanmask] = 0.0
    var[nanmask] = 0.0
    sciivar[nanmask] = 0.0
    tilts[nanmask] = 0.0
    waveimg[nanmask] = 0.0

    waveimg = 1e4*waveimg  # Convert to angstroms for pypeit


    return waveimg, tilts, slit_left, slit_righ, science, sciivar, gpm, base_var, count_scale, \
           slit_left_orig, slit_righ_orig, spec_vals_orig

def jwst_extract_subimgs(e2d_slit, final_slit, intflat_slit):

    # The various multiplicative calibrations we need.
    slit_name = final_slit.name
    waveimg = np.array(final_slit.wavelength.T, dtype=float)
    slit_wcs = final_slit.meta.wcs
    x, y = wcstools.grid_from_bounding_box(slit_wcs.bounding_box, step=(1, 1))
    calra, caldec, calwave = slit_wcs(x, y)
    ra = calra.T
    dec = caldec.T

    # get the source RA and Dec coordinates from the metadata (also located in the header of the fits SCI extension)
    nspec, nspat = ra.shape
    src_ra, src_dec= final_slit.meta.target.ra, final_slit.meta.target.dec

    cal_spat = np.arange(nspat)  # spatial position
    src_trace_ra = np.zeros(nspec)  # Array to hold the source_RA as a function of spectral position
    src_trace_dec = np.zeros(nspec)  # Array to hold the source_DEC as a function of spectral position
    for ispec in range(nspec):
        ra_vs_spat = calra[:, ispec]  #
        # Interpolate y-pixel as a functio of RA onto the source RA
        src_trace_ra[ispec] = np.interp(src_ra, ra_vs_spat[np.isfinite(ra_vs_spat)],
                                                cal_spat[np.isfinite(ra_vs_spat)])
        dec_vs_spat = caldec[:, ispec]
        src_trace_dec[ispec] = np.interp(src_dec, dec_vs_spat[np.isfinite(dec_vs_spat)], cal_spat[np.isfinite(dec_vs_spat)])



    waveimg_from_wcs = calwave.T
    assert np.allclose(waveimg, waveimg_from_wcs, equal_nan=True)

    flatfield = np.array(intflat_slit.data.T, dtype=float) #if intflat_slit is not None else np.ones_like(pathloss)
    pathloss = np.array(final_slit.pathloss_uniform.T, dtype=float) if final_slit.source_type == 'EXTENDED' else \
        np.array(final_slit.pathloss_point.T, dtype=float)
    if pathloss.shape == (0,0):
        msgs.warn('No pathloss for slit {0}'.format(slit_name) + ', setting to 1.0')
        pathloss = np.ones_like(flatfield)

    barshadow = np.array(final_slit.barshadow.T, dtype=float)
    photom_conversion = final_slit.meta.photometry.conversion_megajanskys
    final = np.array(final_slit.data.T, dtype=float)
    rate = np.array(e2d_slit.data.T, dtype=float)
    rate_var_rnoise = np.array(e2d_slit.var_rnoise.T, dtype=float)
    rate_var_poisson = np.array(e2d_slit.var_poisson.T, dtype=float)
    # This is currently buggy as it includes flat field error
    rate_var_tot = np.square(np.array(e2d_slit.err.T, dtype=float))
    dq = np.array(final_slit.dq.T, dtype=int)


    # Generate some tilts and a spatial image
    finitemask = np.isfinite(waveimg)
    wave_min, wave_max = np.min(waveimg[finitemask]), np.max(waveimg[finitemask])

    tilts = np.zeros_like(waveimg)
    tilts[finitemask] = (waveimg[finitemask] - wave_min) / (wave_max - wave_min)

    # TODO Fix this spat_pix to make it increasing with pixel. For now don't use it
    # This currnetly depends on poisition angle which I need to hack to fix
    # ra_min, ra_max = np.min(ra_sub[finitemask_sub]),  np.max(ra_sub[finitemask_sub])
    # spat_pix_sub = np.zeros_like(ra_sub)
    # spat_pix_sub[finitemask_sub] = spat_lo + (ra[finitemask_sub] - ra_min) / (ra_max - ra_min) * (nspat_sub - 1)

    # Get slit bounadries
    slit_left, slit_righ = jwst_get_slits(finitemask)

    ########################
    # The image segment being used for each slit
    spec_lo = final_slit.xstart - 1
    spec_hi = spec_lo + final_slit.xsize
    spat_lo = final_slit.ystart - 1
    spat_hi = spat_lo + final_slit.ysize
    # slice object for the segment
    slit_slice = np.s_[spec_lo: spec_hi, spat_lo: spat_hi]

    slit_left_orig = spat_lo + slit_left
    slit_righ_orig = spat_lo + slit_righ
    spec_vals_orig = spec_lo + np.arange(spec_hi - spec_lo)


    return slit_left, slit_righ, slit_left_orig, slit_righ_orig, spec_vals_orig, src_trace_ra, src_trace_dec, \
           rate, rate_var_rnoise, rate_var_poisson, rate_var_tot, dq, \
           ra, dec, waveimg, tilts, flatfield, pathloss, barshadow, photom_conversion, final


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



# TODO Deprecated. This won't work now that I realize slits can overlap. So just use it as a visualization tool or something??
def jwst_populate_calibs(nspec, nspat, e2d_multi, final_multi, intflat_multi):

    ra = np.zeros((nspec, nspat))
    dec = np.zeros((nspec, nspat))
    waveimg = np.zeros((nspec, nspat))
    tilts = np.zeros((nspec, nspat))
    # The product of these three are the total flat so we instantiate all to 1.0
    flatfield = np.ones((nspec, nspat))
    pathloss = np.ones((nspec, nspat))
    barshadow = np.ones((nspec, nspat))
    photom_conversion = np.zeros((nspec, nspat)) # Currently a constant, but build the possiibility to have it be an image
    calwebb_final = np.zeros((nspec, nspat))
    subimg_count = np.zeros((nspec, nspat), dtype=int)
    slit_name_mask = np.zeros((nspec, nspat), dtype=int)

    # slit boundary stuff
    nslits = len(final_multi.slits)
    slit_left = np.zeros((nspec, nslits))
    slit_righ = np.zeros((nspec, nslits))
    spec_min = np.zeros(nslits)
    spec_max = np.zeros(nslits)

    # TODO print out a warning message here about the slits that are bad, i.e. nan everwhere
    reduce_gpm = np.ones(nslits, dtype=bool)
    meta_list = []

    for islit in range(nslits):
        # Read in data print out slit name
        slit_name = final_multi.slits[islit].name
        meta_list.append(final_multi.slits[islit].meta)

        slit_left_sub, slit_right_sub, rate_sub, rate_var_rnoise_sub, rate_var_poisson_sub, rate_var_tot_sub, \
        ra_sub, dec_sub, waveimg_sub, tilts_sub, flatfield_sub, pathloss_sub, barshadow_sub, photom_conversion_sub, \
        final_sub = jwst_extract_subimgs(e2d_multi.slits[islit], final_multi.slits[islit], intflat_multi.slits[islit])

        # Determine the slit boundaries using the waveimg
        finitemask_sub = np.isfinite(waveimg_sub)
        if not np.any(finitemask_sub):
            reduce_gpm[islit] = False
            msgs.warn('All nan wavelengths for Slit={:s}. Not extracting calibrations'.format(slit_name))
        else:
            #msgs.info('Extracting calibrations for Slit={:s}'.format(slit_name))

            ########################
            # The image segment being used for each slit
            spec_lo = final_multi.slits[islit].xstart - 1
            spec_hi = spec_lo + final_multi.slits[islit].xsize
            spat_lo = final_multi.slits[islit].ystart - 1
            spat_hi = spat_lo + final_multi.slits[islit].ysize
            # slice object for the segment
            slit_slice = np.s_[spec_lo: spec_hi, spat_lo: spat_hi]

            # Get slit bounadries
            nspec_sub, nspat_sub = ra_sub.shape
            sub_slit_left, sub_slit_righ = jwst_get_slits(finitemask_sub)



            # TODO Fix this spat_pix to make it increasing with pixel. For now don't use it
            # This currnetly depends on poisition angle which I need to hack to fix
            #ra_min, ra_max = np.min(ra_sub[finitemask_sub]),  np.max(ra_sub[finitemask_sub])
            #spat_pix_sub = np.zeros_like(ra_sub)
            #spat_pix_sub[finitemask_sub] = spat_lo + (ra[finitemask_sub] - ra_min) / (ra_max - ra_min) * (nspat_sub - 1)

            slit_left[spec_lo: spec_hi, islit] = spat_lo + sub_slit_left
            slit_righ[spec_lo: spec_hi, islit] = spat_lo + sub_slit_righ
            # This is a hack for now until we figure out how to deal with spec_min and spec_max slits.
            # I'm just setting the boundaries to be everywhere the last defined boundary location in the sub-image
            slit_left[:spec_lo, islit] = spat_lo + sub_slit_left[0]
            slit_left[spec_hi:, islit] = spat_lo + sub_slit_left[-1]
            slit_righ[:spec_lo, islit] = spat_lo + sub_slit_righ[0]
            slit_righ[spec_hi:, islit] = spat_lo + sub_slit_righ[-1]

            spec_min[islit] = spec_lo
            spec_max[islit] = spec_hi

            # Populate the 2d images in the regions where the JWST calibrations are finite
            ra[slit_slice][finitemask_sub] = ra_sub[finitemask_sub]
            dec[slit_slice][finitemask_sub] = dec_sub[finitemask_sub]
            waveimg[slit_slice][finitemask_sub] = waveimg_sub[finitemask_sub]
            tilts[slit_slice][finitemask_sub] = tilts_sub[finitemask_sub]
            flatfield[slit_slice][finitemask_sub] = flatfield_sub[finitemask_sub]
            pathloss[slit_slice][finitemask_sub] = pathloss_sub[finitemask_sub]
            barshadow[slit_slice][finitemask_sub] = barshadow_sub[finitemask_sub]
            photom_conversion[slit_slice][finitemask_sub] = photom_conversion_sub # Currently just a float but may be an image in the future
            calwebb_final[slit_slice][finitemask_sub] = final_sub[finitemask_sub]
            subimg_count[slit_slice][finitemask_sub] += 1

    return reduce_gpm, ra, dec, waveimg, tilts, flatfield, pathloss, barshadow, photom_conversion, calwebb_final, subimg_count, \
           slit_left, slit_righ, spec_min, spec_max, meta_list




#TODO Deprecated
def jwst_proc_old(e2d_slit, final_slit, intflat_slit=None, kludge_err=1.0):


    # Try to reverse engineer all the things they multiply into the data
    slit_name = e2d_slit.name

    t_eff = e2d_slit.meta.exposure.effective_exposure_time
    # TODO I don't know how the t_eff quantity is defined. Better would be some proxy for the exposure time per pixel
    # The science data is divided by (flat*pathloss*barshadow) and then multiplied by photom_conversion. Since
    # we work in units of counts, we divide by the photom conversion and multiply by t_eff.

    # This is the raw e2d data before the pipeline does idiotic things
    raw_data_counts = np.array(e2d_slit.data.T, dtype=float)*t_eff
    raw_var_poisson = kludge_err**2*np.array(e2d_slit.var_poisson.T, dtype=float)*t_eff**2
    raw_var_rnoise = kludge_err**2*np.array(e2d_slit.var_rnoise.T, dtype=float)*t_eff**2
    raw_var = kludge_err**2*np.square(np.array(e2d_slit.err.T, dtype=float))*t_eff**2

    photom_conversion = final_slit.meta.photometry.conversion_megajanskys
    pathloss = np.array(final_slit.pathloss_uniform.T, dtype=float) if final_slit.source_type == 'EXTENDED' else \
        np.array(final_slit.pathloss_point.T, dtype=float)
    if pathloss.shape == (0,0):
        msgs.warn('No pathloss for slit {0}'.format(slit_name) + ', setting to 1.0')
        pathloss = np.ones_like(raw_data_counts)
    flatfield = np.array(intflat_slit.data.T, dtype=float) if intflat_slit is not None else np.ones_like(raw_data_counts)
    barshadow = np.array(final_slit.barshadow.T, dtype=float)

    # This is the conversion between final2d and e2d, i.e. final2d = jwst_scale*e2d
    jwst_scale = photom_conversion / flatfield / pathloss / barshadow
    flux_to_counts = t_eff / photom_conversion # This converts s2d outputs of flux to counts.

    #science = np.array(e2d_slit.data.T, dtype=float) * flux_to_counts

    total_flat = flatfield * pathloss * barshadow
    total_flat_square = np.square(total_flat)
    count_scale = inverse(total_flat)  # This is the quantity that goes into PypeIt

    science, flat_bpm = flat.flatfield(raw_data_counts, total_flat)
    var_poisson, _ = flat.flatfield(raw_var_poisson, total_flat**2)
    base_var, _ = flat.flatfield(raw_var_rnoise, total_flat**2)
    var, _ = flat.flatfield(raw_var, total_flat**2)
    sciivar = inverse(var)

    # TODO Currently the var_flat is nonsense I think and so I'm just going to use the var_poisson and var_rnoise to get
    # the noise. If this gets fixed use the line below which includes the var_flat.
    # err = kludge_err*np.array(slit.err.T, dtype=float)*flux_to_counts
    #var_poisson = slit.var_poisson.T * flux_to_counts ** 2
    #var_rnoise = slit.var_rnoise.T * flux_to_counts ** 2
    #var = kludge_err**2*np.array(var_poisson + var_rnoise, dtype=float)
    # This needs to be multiplied by count_scale to get it into units of counts which is what pypeit requires. I checked
    # that this base_var is equal to e2d.var_rnoise if you remove the flux_to_counts factor.
    # base_var = np.array(final2d.slits[islit].var_rnoise.T, dtype=float)*flux_to_counts**2*count_scale**2
    #base_var = np.array(slit.var_rnoise.T, dtype=float) * flux_to_counts ** 2

    # TODO I'm unsure about these
    dq = np.array(final_slit.dq.T, dtype=int)
    waveimg = np.array(final_slit.wavelength.T, dtype=float)

    gpm = np.logical_not(dq & DO_NOT_USE)

    finite_mask = np.isfinite(science)
    nanmask = np.logical_not(finite_mask)
    science[nanmask] = 0.0
    # err[nanmask] = 0.0
    var[nanmask] = 0.0
    sciivar = inverse(var) * gpm
    base_var[nanmask] = 0.0
    count_scale[nanmask] = 0.0
    # Wave nanmask is different from data nanmask
    slit_wcs = final_slit.meta.wcs
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
    slit_left, slit_righ = jwst_get_slits(finite_mask)
    # Generate some tilts and a spatial image
    tilts = np.zeros_like(waveimg)
    tilts[np.isfinite(waveimg)] = (waveimg[np.isfinite(waveimg)] - wave_min) / (wave_max - wave_min)

    # TODO Fix this spat_pix to make it increasing with pixel. For now don't use it
    nspec, nspat = science.shape
    spat_pix = (ra - ra_min) / (ra_max - ra_min) * (nspat - 1)
    spat_pix[nanmask_ra] = 0.0

    embed()

    return science, sciivar, gpm, base_var, count_scale, tilts, waveimg, finite_mask, slit_left, slit_righ, t_eff


