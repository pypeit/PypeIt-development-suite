
import copy
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from pypeit.display import display
from pypeit.core import fitting
from gwcs import wcstools
from matplotlib import pyplot as plt
from jwst import datamodels
from pypeit.utils import inverse, zero_not_finite
DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']
from pypeit import msgs
from pypeit import datamodel
from pypeit.core import flat
from pypeit.core import procimg
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.pypeitimage import PypeItImage
from pypeit.images.mosaic import Mosaic
from pypeit import slittrace, spec2dobj
from pypeit import find_objects, extraction

from astropy.io import fits
from astropy.wcs import WCS
import grismconf



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


def jwst_nircam_proc(rate_file, configfile, RA, DEC, kludge_err=1.0, noise_floor=0.01, saturation=65000):


    # TODO Use the spat_img which PypeIt can take, but we are not using!!
    rate_obj, raw_sub, var_tot_sub, var_poisson_sub, var_rnoise_sub, dq_sub, waveimg_sub, spat_img_sub = jwst_nircam_subimgs(
        configfile, RA, DEC, rate_file)
    t_eff = rate_obj.meta.exposure.effective_exposure_time

    # Read in the output after msa_flagging. Extract the sub-images, rotate to PypeIt format.
    rate = raw_sub.T
    rate_var_tot = var_tot_sub.T
    rate_var_poisson = var_poisson_sub.T
    rate_var_rnoise = var_rnoise_sub.T
    # TODO Check that the errors don't have nonsense from the flat field error budget
    dq = dq_sub.T
    # Now perform the image processing
    raw_counts = rate * t_eff
    raw_var_poisson = kludge_err ** 2 * rate_var_poisson * t_eff ** 2
    raw_var_rnoise = kludge_err ** 2 * rate_var_rnoise * t_eff ** 2
    # Is this correct? I'm not sure I should be using their poisson variance for the noise floor
    raw_var = procimg.variance_model(raw_var_rnoise, counts=raw_var_poisson, noise_floor=noise_floor)
    # TODO This  is a hack until I can understand how to get rid of the hot pixels in the JWST variance arrays using DQ flags.
    # I don't know what the value of this parameter currently set to 20 should be?? Look into this via a github issue.
    # raw_gpm = (raw_var_rnoise < 20.0*ronoise**2) & (raw_var_poisson < saturation)
    raw_gpm = (raw_var_rnoise < saturation) & (raw_var_poisson < saturation)
    # raw_var_poisson + raw_var_rnoise # TODO Leaving out problematic flat field term from pipeline

    # This is the conversion between final2d and e2d, i.e. final2d = jwst_scale*e2d
    # total_flat = flatfield*pathloss*barshadow
    # flux_to_counts = t_eff / photom_conversion  # This converts s2d outputs of flux to counts.
    # jwst_scale = photom_conversion/flatfield/pathloss/barshadow

    # TODO Perform the flat field correction yourself so as to update the noise model. Setting this to unit
    total_flat = np.ones_like(rate)
    finitemask = np.isfinite(rate, dtype=bool)  # This was from NIRSPEC and may not be necessary
    total_flat_square = np.square(total_flat)

    count_scale = inverse(total_flat)  # This is the quantity that goes into PypeIt for var modeling
    science, flat_bpm = flat.flatfield(raw_counts, total_flat)
    var_poisson, _ = flat.flatfield(raw_var_poisson, total_flat_square)
    base_var, _ = flat.flatfield(raw_var_rnoise, total_flat_square)
    var, _ = flat.flatfield(raw_var, total_flat_square)
    sciivar = inverse(var)
    dq_gpm = np.logical_not(dq & DO_NOT_USE)
    gpm = finitemask & dq_gpm & np.logical_not(flat_bpm) & (sciivar > 0.0) & raw_gpm

    nanmask = np.logical_not(finitemask)
    count_scale[nanmask] = 0.0
    science[nanmask] = 0.0
    var_poisson[nanmask] = 0.0
    base_var[nanmask] = 0.0
    var[nanmask] = 0.0
    sciivar[nanmask] = 0.0

    # TODO This is kludge city!!
    waveimg = 1e4 * waveimg_sub.T
    wave_min, wave_max = np.min(waveimg), np.max(waveimg)
    tilts = (waveimg - wave_min) / (wave_max - wave_min)
    spat_img = spat_img_sub.T

    return rate_obj, science, sciivar, gpm, dq_gpm, base_var, count_scale, finitemask, tilts, waveimg, spat_img

def jwst_nircam_subimgs(configfile, RA, DEC, rate_file, senscorrect=False, h5name=False, yoffset=0., yhsize=20., use_wcs=False):


    hdu = fits.open(rate_file)
    wcs = WCS(hdu[1].header)

    # Information about observing mode of this rate file
    h = hdu[0].header
    filt = h["FILTER"]  # Filter name, e.g. F410M
    grism = h["PUPIL"][-1]  # R or C
    module = h["MODULE"]  # Which NIRCAM module, A or B

    C = grismconf.Config(configfile)

    # Compute the position of the source in the image in pixel coordinates
    grism_with_wcs = datamodels.open(rate_file)

    world_to_pix = grism_with_wcs.meta.wcs.get_transform('world', 'detector')
    x0, y0, foo, foo2 = world_to_pix(RA, DEC, 0, 0)
    print('Target is at ', x0, y0)

    y0 = y0 + yoffset
    xref = x0  # 2048.
    yref = y0  # 2048.

    t = C.INVDISPL('+1', xref, yref, np.array([2.99, 4.21]))
    dx = C.DISPX('+1', xref, yref, t)
    dy = C.DISPY('+1', xref, yref, t)
    x_line = x0 + dx
    y_line = y0 + dy
    minx0 = np.max([0, np.int32(np.min(x_line))])
    maxx0 = np.min([2047, np.int32(np.max(x_line))])  # where do these numbers come from??
    miny0 = np.max([0, np.int32(np.min(y_line - yhsize))])
    maxy0 = np.min([2047, np.int32(np.max(y_line + yhsize))])

    data = grism_with_wcs.data
    var_tot = grism_with_wcs.err ** 2
    var_poisson = grism_with_wcs.var_poisson
    var_rnoise = grism_with_wcs.var_rnoise
    dq = grism_with_wcs.dq
    #data = hdu["SCI"].data
    #err = hdu["ERR"].data
    #dq = hdu["DQ"].data

    # We trim our data to be the stamp containing the spectrum we want to extract
    data = data[miny0:maxy0 + 1, minx0:maxx0 + 1]
    var_tot = var_tot[miny0:maxy0 + 1, minx0:maxx0 + 1]
    var_poisson = var_poisson[miny0:maxy0 + 1, minx0:maxx0 + 1]
    var_rnoise = var_rnoise[miny0:maxy0 + 1, minx0:maxx0 + 1]
    dq = dq[miny0:maxy0 + 1, minx0:maxx0 + 1]
    print(data.shape, var_tot.shape, dq.shape)  # , model0.shape)

    # These are the coordinates of all the pixels in our 2D stamp, but in the full image (wrt calibration is known)
    ys, xs = np.indices((maxy0 - miny0 + 1, maxx0 - minx0 + 1))
    # xs and ys are now the relative dx and dy offsets from the position of our source. They are both 2D arrays of x and y coordinates.
    xs = xs + minx0 - x0
    ys = ys + miny0 - y0

    # Depending on whether the grism disperse in the x or y direction, we use the INVDISPX or INVDISPY functions
    # to compute the value for t for every pixel in our 2D stamps
    if grism == "R":
        ts = C.INVDISPX("+1", x0, y0, xs)
        dys = C.DISPY("+1", x0, y0, ts) + ys - 2 * C.DISPY("+1", x0, y0, ts)

    if grism == "C":
        ts = C.INVDISPY("+1", x0, y0, ys)
        dys = C.DISPX("+1", x0, y0, ts) + xs

    # Now compute the wavelength of every pixel in our 2D stamp
    ws = C.DISPL("+1", x0, y0, ts)

    # Now, depending of whether things are in the row or col, we transpose things so that we can look at them properly (i.e. row direction)
    if grism == "C":
        # m = np.transpose(model0) # The model counts in each pixel
        l = np.transpose(ws)  # The wavelength of each pixel
        d = np.transpose(data)  # The data counts in each pixel
        var_tot_out = np.transpose(var_tot)  # THe data error estimates in each pixel
        var_poisson_out = np.transpose(var_poisson)  # THe data error estimates in each pixel
        var_rnoise_out = np.transpose(var_rnoise)  # THe data error estimates in each pixel
        q = np.transpose(dq)  # The data DQ in each pixel
        y = np.transpose(dys)  # The cross-dispersion distance of each pixel from the trace

    if grism == "R":
        # m = model0
        l = ws
        d = data
        var_tot_out = var_tot
        var_poisson_out = var_poisson
        var_rnoise_out = var_rnoise
        q = dq
        y = dys

    # correct for the sensitivity function of the filter
    if senscorrect:
        lam = np.nanmean(l, axis=0)
        sens = 1E-18 * C.SENS["+1"](lam)  # always use the one for module A because Module B data has been rescaled???

        d = d / sens
        var_tot_out = var_tot / sens
        var_poisson_out = var_poisson / sens
        var_rnoise_out= var_rnoise / sens

    # TODO change variable names to be informative. Annoying that we have to recast
    return grism_with_wcs, d.astype(float), var_tot_out.astype(float), var_poisson_out.astype(float), var_rnoise_out.astype(float), \
           q.astype(bool), l.astype(float), y.astype(float)  # , m


def jwst_get_slits(finitemask, polyorder=5, function='legendre', debug=False):

    slit_left = fit_slit(finitemask, 'left', polyorder=polyorder, function=function, debug=debug)
    slit_righ = fit_slit(finitemask, 'righ', polyorder=polyorder, function=function, debug=debug)
    return slit_left, slit_righ


def jwst_proc(msa_data, slit_slice, finitemask, flatfield, pathloss, barshadow, photom_conversion, ronoise,
              kludge_err=1.0, saturation=65000, noise_floor=0.01, use_flat=True):


    #slit_slice, slit_left, slit_righ, slit_left_orig, slit_righ_orig, spec_vals_orig, src_trace_ra, src_trace_dec, dq, \
    #ra, dec, waveimg, tilts, flatfield, pathloss, barshadow, photom_conversion, final = jwst_extract_subimgs(
    #    final_slit, intflat_slit)

    # Now deal with the image processing
    #if not np.any(finitemask):
    #    return (None,)*21

    # Read in the output after msa_flagging. Extract the sub-images, rotate to PypeIt format.
    rate = np.array(msa_data.data.T[slit_slice], dtype=float)
    rate_var_rnoise = np.array(msa_data.var_rnoise.T[slit_slice], dtype=float)
    rate_var_poisson = np.array(msa_data.var_poisson.T[slit_slice], dtype=float)
    # This is currently buggy as it includes flat field error
    # rate_var_tot = np.square(np.array(e2d_slit.err.T, dtype=float))
    dq = np.array(msa_data.dq.T[slit_slice], dtype=int)
    t_eff = msa_data.meta.exposure.effective_exposure_time

    # Now perform the image processing
    raw_counts = rate*t_eff
    raw_var_poisson = kludge_err**2*rate_var_poisson*t_eff**2
    raw_var_rnoise = kludge_err**2*rate_var_rnoise*t_eff**2
    # Is this correct? I'm not sure I should be using their poisson variance for the noise floor
    raw_var = procimg.variance_model(raw_var_rnoise, counts = raw_var_poisson, noise_floor=noise_floor)
    # TODO This  is a hack until I can understand how to get rid of the hot pixels in the JWST variance arrays using DQ flags.
    # I don't know what the value of this parameter currently set to 20 should be?? Look into this via a github issue.
    #raw_gpm = (raw_var_rnoise < 20.0*ronoise**2) & (raw_var_poisson < saturation)
    raw_gpm = (raw_var_rnoise < saturation) & (raw_var_poisson < saturation)
    #raw_var_poisson + raw_var_rnoise # TODO Leaving out problematic flat field term from pipeline

    # This is the conversion between final2d and e2d, i.e. final2d = jwst_scale*e2d
    # total_flat = flatfield*pathloss*barshadow
    #flux_to_counts = t_eff / photom_conversion  # This converts s2d outputs of flux to counts.
    #jwst_scale = photom_conversion/flatfield/pathloss/barshadow

    # This is what Kyle is using in his code
    if use_flat:
        total_flat = flatfield*pathloss*barshadow/photom_conversion
    else:
        total_flat = pathloss*barshadow

    total_flat_square = np.square(total_flat)

    count_scale = inverse(total_flat)  # This is the quantity that goes into PypeIt for var modeling
    science, flat_bpm = flat.flatfield(raw_counts, total_flat)
    var_poisson, _ = flat.flatfield(raw_var_poisson, total_flat_square)
    base_var, _ = flat.flatfield(raw_var_rnoise, total_flat_square)
    var, _ = flat.flatfield(raw_var, total_flat_square)
    sciivar = inverse(var)
    dq_gpm = np.logical_not(dq & DO_NOT_USE)
    gpm = finitemask & dq_gpm & np.logical_not(flat_bpm) & (sciivar > 0.0) & raw_gpm & np.isfinite(science) & np.isfinite(sciivar)

    # This finitemask is based on where the waveimg is defined. Note however, the JWST images have nans in other
    # places which is exposure specific.
    nanmask = np.logical_not(finitemask)
    count_scale[nanmask] = 0.0
    science[nanmask] = 0.0
    var_poisson[nanmask] = 0.0
    base_var[nanmask] = 0.0
    var[nanmask] = 0.0
    sciivar[nanmask] = 0.0
    rn2_img = np.zeros_like(science)
    rn2_img[finitemask] = ronoise**2

    # Make sure that we zero out any nan pixels since this causes problems in PypeIt
    return zero_not_finite(science), zero_not_finite(sciivar), gpm, zero_not_finite(base_var), \
        zero_not_finite(count_scale), rn2_img

class NIRSpecSlitCalibrations(datamodel.DataContainer):
    version = '1.0.0'
    """Datamodel version."""

    datamodel = {'slit_name': dict(otype=str, descr='Name of slit'),
                 'source_name': dict(otype=str, descr='Name of source'),
                 'on_detector': dict(otype=bool, descr='True if the slit is on the detector, otherwise False'),
                 'det_name': dict(otype=str, descr='Name of NIRSpec detector, i.e. either NRS1 or NRS2'),
                 'detector': dict(otype=DetectorContainer,
                                  descr='The detector (see :class:`~pypeit.images.detector_container.DetectorContainer`) '
                                        'parameters'),
                 'slit_indx': dict(otype=int, descr='Index in JWST datamodel if slit is on nrs1, else -1'),
                 'slit_slice': dict(otype=tuple, descr='Slice for nrs1'),
                 'slit_left': dict(otype=np.ndarray, atype=float, descr='Left slit edge for nrs1'),
                 'slit_righ': dict(otype=np.ndarray, atype=float, descr='Right slit edge for nrs1'),
                 'slit_left_orig': dict(otype=np.ndarray, atype=float, descr='Original left slit edge for nrs1'),
                 'slit_righ_orig': dict(otype=np.ndarray, atype=float, descr='Original right slit edge for nrs1'),
                 'spec_vals_orig': dict(otype=np.ndarray, atype=np.int64, descr='Original spectral values for nrs1'),
                 'src_trace_ra': dict(otype=np.ndarray, atype=float, descr='Source trace RA for nrs1'),
                 'src_trace_dec': dict(otype=np.ndarray, atype=float, descr='Source trace DEC for nrs1'),
                 'dq_sub': dict(otype=np.ndarray, atype=np.int64, descr='DQ flags images for nrs1'),
                 'ra': dict(otype=np.ndarray, atype=float, descr='RA image map for nrs1'),
                 'dec': dict(otype=np.ndarray, atype=float, descr='DEC image map for nrs1'),
                 'finitemask': dict(otype=np.ndarray, atype=np.bool_, descr='Mask indicating where wavelengths are defined for nrs1'),
                 'waveimg': dict(otype=np.ndarray, atype=float, descr='Wavelength image for nrs1'),
                 'tilts': dict(otype=np.ndarray, atype=float, descr='Tilt image for nrs1'),
                 'flatfield': dict(otype=np.ndarray, atype=float, descr='Flat field image for nrs1'),
                 'pathloss': dict(otype=np.ndarray, atype=float, descr='Pathloss image for nrs1'),
                 'barshadow': dict(otype=np.ndarray, atype=float, descr='Barshadow image for nrs1'),
                 'photom_conversion': dict(otype=float, descr='Photom conversion in MJy for nrs1'),
                 'calwebb_proc': dict(otype=np.ndarray, atype=float, descr='Calwebb processed image for nrs1'),
                 }

    """DataContainer datamodel."""

    internals = ['_indx', 'slit_names', 'intflat_slit_names',]
    def __init__(self, detector, ms_model, ms_model_flat, slit_name):

        # Instantiate as an empty DataContainer
        super().__init__()
        self.slit_name = slit_name
        self.det_name = ms_model.meta.instrument.detector
        self.detector = detector
        # Is this slit on nrs1?
        slit_names = np.array([slit.name for slit in ms_model.slits])
        intflat_slit_names = np.array([slit.name for slit in ms_model_flat.slits])
        _indx = np.where((slit_names == slit_name) & (intflat_slit_names == slit_name))[0]
        self.on_detector = _indx.size != 0
        self.slit_indx = int(_indx[0]) if self.on_detector else -1
        self.source_name = ms_model.slits[self.slit_indx].source_name

        # Assign calibrations
        if self.on_detector:
            self.slit_slice, self.slit_left, self.slit_righ, self.slit_left_orig, self.slit_righ_orig, self.spec_vals_orig, \
                self.src_trace_ra, self.src_trace_dec, self.dq_sub, self.ra, self.dec, \
                self.finitemask, self.waveimg, self.tilts, self.flatfield, self.pathloss, \
                self.barshadow, self.photom_conversion, self.calwebb_proc = jwst_extract_subimgs(
                ms_model.slits[self.slit_indx], ms_model_flat.slits[self.slit_indx])

    def show(self):
        # Connect to an open ginga window, or open a new one
        display.connect_to_ginga(raise_err=True, allow_new=True)
        ch_list = []
        image_list = [self.waveimg, self.tilts, self.flatfield, self.pathloss, self.barshadow, self.calwebb_proc]
        waveimg_list = [None] + [self.waveimg] * 5
        chname_list = ['wave_{:s}'.format(self.det_name), 'tilts_{:s}'.format(self.det_name), 'flat_{:s}'.format(self.det_name),
                       'pathloss_{:s}'.format(self.det_name), 'barshadow_{:s}'.format(self.det_name), 'calwebb_{:s}'.format(self.det_name)]
        cuts_list = [None] * 6
        cuts_list[-2]= (0.0, 1.0)
        for image, waveimg, chname, cuts in zip(image_list, waveimg_list, chname_list, cuts_list):
            viewer, ch_tmp = display.show_image(image, waveimg=waveimg, chname=chname, cuts=cuts)
            ch_list.append(ch_tmp)
        display.show_slits(viewer, ch_list[0], self.slit_left, self.slit_righ, pstep=1, slit_ids=np.array([self.slit_name]))
        display.show_trace(viewer, ch_list[-1], self.src_trace_ra, 'trace-RA_{:s}'.format(self.det_name),color='#f0e442', pstep=1)
        display.show_trace(viewer, ch_list[-1], self.src_trace_dec, 'trace-DEC_{:s}'.format(self.det_name), color='#f0e442', pstep=1)

        # After displaying all the images sync up the images with WCS_MATCH
        shell = viewer.shell()
        shell.start_global_plugin('WCSMatch')
        shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [ch_list[0]], {})

        return viewer, ch_list



def jwst_mosaic(image_model_tuple, Calibrations_tuple, kludge_err=1.0,
                noise_floor=0.01, show=False, bkg_image_model_tuple=(None,None)):

    dets = [1,2]
    sciimg_list, sciivar_list, gpm_list, base_var_list, count_scale_list, rn2_img_list = [], [], [], [], [], []
    waveimg_list, tilts_list, slit_slice_list, slit_left_list, slit_righ_list, det_list, calib_list = [], [], [], [], [], [], []
    for det, image_model, bkg_image_model, Calib in zip(dets, image_model_tuple, bkg_image_model_tuple, Calibrations_tuple):

        if Calib.on_detector:
            sciimg, sciivar, gpm, base_var, count_scale, rn2_img = jwst_proc(
                image_model, Calib.slit_slice, Calib.finitemask, Calib.flatfield, Calib.pathloss, Calib.barshadow,
                Calib.photom_conversion, Calib.detector.ronoise, noise_floor=noise_floor, kludge_err=kludge_err)
            sciimg_list.append(sciimg)
            sciivar_list.append(sciivar)
            gpm_list.append(gpm)
            base_var_list.append(base_var)
            count_scale_list.append(count_scale)
            rn2_img_list.append(rn2_img)
            waveimg_list.append(Calib.waveimg)
            tilts_list.append(Calib.tilts)
            slit_slice_list.append(Calib.slit_slice)
            slit_left_list.append(Calib.slit_left)
            slit_righ_list.append(Calib.slit_righ)
            det_list.append(Calib.detector)
            calib_list.append(Calib)
            if show:
                display.connect_to_ginga(raise_err=True, allow_new=True)
                # Show the raw rate image
                rate_image = np.array(image_model.data.T)
                viewer, ch_raw = display.show_image(rate_image, cuts=get_cuts(rate_image),
                                                        chname='raw_rate_{:s}'.format(Calib.det_name))
                display.show_slits(viewer, ch_raw, Calib.slit_left_orig, Calib.slit_righ_orig,
                                   spec_vals=Calib.spec_vals_orig, pstep=1, slit_ids=np.array([Calib.slit_name]))
                # If bkg_image models were provided, show the difference image
                if bkg_image_model is not None:
                    bkg_rate_image = np.array(bkg_image_model.data.T)
                    viewer, ch_bkg = display.show_image(rate_image - bkg_rate_image,
                                                            cuts=get_cuts(rate_image - bkg_rate_image),
                                                            chname='diff_rate_{:s}'.format(Calib.det_name))
                    display.show_slits(viewer, ch_bkg, Calib.slit_left_orig, Calib.slit_righ_orig,
                                       spec_vals=Calib.spec_vals_orig, pstep=1, slit_ids=np.array([Calib.slit_name]))

                # Show the calibrations
                viewer, ch_list = Calib.show()
                # Show the pypeit processed image
                viewer_pypeit, ch_pypeit = display.show_image(sciimg, waveimg=Calib.waveimg, cuts=get_cuts(sciimg),
                                                              chname='pypeit_{:s}'.format(Calib.det_name))
                display.show_slits(viewer_pypeit, ch_pypeit, Calib.slit_left, Calib.slit_righ, pstep=1,
                                   slit_ids=np.array([Calib.slit_name]))
                # Sync up the images that have the same shape, namely calibration image and pypeit processed image
                shell = viewer.shell()
                shell.start_global_plugin('WCSMatch')
                shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [ch_list[-1]], {})

    ndet = len(calib_list)

    # TODO I would like to create an image indicating which detector contributed to which pixels
    if ndet == 1:
        # Assign the data frames
        sciimg_tot, sciivar_tot, gpm_tot, base_var_tot, count_scale_tot = \
            sciimg_list[0], sciivar_list[0], gpm_list[0], base_var_list[0], count_scale_list[0]
        # Assign the calibraitons
        waveimg_tot, tilts_tot, rn2_img_tot = waveimg_list[0], tilts_list[0], rn2_img_list[0]
        slit_left_tot, slit_righ_tot = slit_left_list[0], slit_righ_list[0]
        det_or_mosaic = det_list[0]
        shape = sciimg_tot.shape
    elif ndet == 2:
        # Spatial offset is relative to NRS1
        detector_gap = int(calib_list[0].detector.xgap)
        spat_offset = (slit_slice_list[0][1].start - slit_slice_list[1][1].start)
        spec_lo1, spec_hi1 = 0, sciimg_list[0].shape[0]
        spec_lo2, spec_hi2 = sciimg_list[0].shape[0] + detector_gap, \
                             sciimg_list[0].shape[0] + detector_gap + sciimg_list[1].shape[0]
        shape = (sciimg_list[0].shape[0] + sciimg_list[1].shape[0] + detector_gap,
                 np.max([sciimg_list[0].shape[1], sciimg_list[1].shape[1]]) + np.abs(spat_offset))
        if spat_offset < 0:
            spat_lo1, spat_hi1 = -spat_offset, sciimg_list[0].shape[1] - spat_offset
            spat_lo2, spat_hi2 = 0, sciimg_list[1].shape[1]
        else:
            spat_lo1, spat_hi1 = 0, sciimg_list[0].shape[1]
            spat_lo2, spat_hi2 = spat_offset, spat_offset + sciimg_list[1].shape[1]

        nrs1_slice = np.s_[spec_lo1: spec_hi1, spat_lo1: spat_hi1]
        nrs2_slice = np.s_[spec_lo2: spec_hi2, spat_lo2: spat_hi2]

        sciimg_tot = np.zeros(shape)
        sciivar_tot = np.zeros(shape)
        gpm_tot = np.zeros(shape, dtype=bool)
        base_var_tot = np.zeros(shape)
        count_scale_tot = np.zeros(shape)
        rn2_img_tot = np.zeros(shape)
        sciimg_tot[nrs1_slice] = sciimg_list[0]
        sciimg_tot[nrs2_slice] = sciimg_list[1]
        sciivar_tot[nrs1_slice] = sciivar_list[0]
        sciivar_tot[nrs2_slice] = sciivar_list[1]
        gpm_tot[nrs1_slice] = gpm_list[0]
        gpm_tot[nrs2_slice] = gpm_list[1]
        base_var_tot[nrs1_slice] = base_var_list[0]
        base_var_tot[nrs2_slice] = base_var_list[1]
        count_scale_tot[nrs1_slice] = count_scale_list[0]
        count_scale_tot[nrs2_slice] = count_scale_list[1]
        rn2_img_tot[nrs1_slice] = rn2_img_list[0]
        rn2_img_tot[nrs2_slice] = rn2_img_list[1]

        waveimg_tot = np.full(shape, np.nan)
        waveimg_tot[nrs1_slice] = waveimg_list[0]
        waveimg_tot[nrs2_slice] = waveimg_list[1]
        finitemask_tot = np.isfinite(waveimg_tot)
        wave_min, wave_max = np.min(waveimg_tot[finitemask_tot]), np.max(waveimg_tot[finitemask_tot])
        tilts_tot = np.zeros(shape)
        tilts_tot[finitemask_tot] = (waveimg_tot[finitemask_tot] - wave_min) / (wave_max - wave_min)
        slit_left_tot, slit_righ_tot = jwst_get_slits(finitemask_tot)
        waveimg_tot[np.logical_not(finitemask_tot)] = 0.0

        det_or_mosaic = Mosaic(1, np.array(det_list), shape, None, None, None, None)
    else:
        msgs.error('Invalid number of detectors. There is a problem with this slit')


    # Instantiate
    sciImg = PypeItImage(image=zero_not_finite(sciimg_tot), ivar=zero_not_finite(sciivar_tot),
                         base_var=zero_not_finite(base_var_tot),
                         img_scale=zero_not_finite(count_scale_tot),
                         rn2img=zero_not_finite(rn2_img_tot),
                         detector=det_or_mosaic, bpm=np.logical_not(gpm_tot))
    slits = slittrace.SlitTraceSet(slit_left_tot, slit_righ_tot, 'MultiSlit', detname=det_or_mosaic.name,
                                   nspat=int(shape[1]),PYP_SPEC='jwst_nirspec')



    return sciImg, slits, zero_not_finite(waveimg_tot), zero_not_finite(tilts_tot), ndet


def jwst_reduce(sciImg, slits, waveimg, tilts, spectrograph, par, show=False, find_negative=False, bkg_redux=False,
                  clear_ginga=True, show_peaks=False, show_skysub_fit=False, basename=None):
    """
    Method to run the reduction on coadd2d pseudo images

    Args:
        pseudo_dict (dict):
           Dictionary containing coadd2d pseudo images
        show (bool):
           If True, show the outputs to ginga and the screen analogous to run_pypeit with the -s option
        show_peaks (bool):
           If True, plot the object finding QA to the screen.
        basename (str):
           The basename for the spec2d output files.

    Returns:

    """


    #
    slitmask_pseudo = slits.slit_img()
    sciImg.build_mask(slitmask=slitmask_pseudo)


    # TODO implement manual extraction.
    manual_obj = None
    # Get bpm mask. There should not be any masked slits because we excluded those already
    # before the coadd, but we need to pass a bpm to FindObjects and Extract

    # Initiate FindObjects object
    objFind = find_objects.FindObjects.get_instance(sciImg, slits, spectrograph, par,
                                                    'science_coadd2d', tilts=tilts,
                                                    bkg_redux=bkg_redux, manual=manual_obj,
                                                    find_negative=find_negative, basename=basename,
                                                    clear_ginga=clear_ginga, show=show)
    global_sky0, sobjs_obj = objFind.run(show_peaks=show or show_peaks, show_skysub_fit=show_skysub_fit)

    # TODO add this as optional to objFind.run()?

    # For JWST bkg_redux we don't perform any global skysubtraction whatsoever. We could in principle check
    # bkg_redux here, but that is equiavlent since skip_skysub is set to True in the bkg_redux case.
    if par['reduce']['findobj']['skip_skysub']:
        final_global_sky=global_sky0
    else:
        skymask = objFind.create_skymask(sobjs_obj)
        final_global_sky = objFind.global_skysub(previous_sky=global_sky0, skymask=skymask, show=show)


    # Initiate Extract object
    exTract = extraction.Extract.get_instance(sciImg, slits, sobjs_obj, spectrograph, par,
                                              'science_coadd2d', global_sky=final_global_sky, tilts=tilts,
                                              waveimg=waveimg, bkg_redux=bkg_redux,
                                              basename=basename, show=show)

    skymodel, objmodel, ivarmodel, outmask, sobjs, _, _ = exTract.run()

    # TODO -- Do this upstream
    # Tack on detector and wavelength RMS
    for sobj in sobjs:
        sobj.DETECTOR = sciImg.detector

    # Construct the Spec2DObj
    spec2DObj = spec2dobj.Spec2DObj(sciimg=sciImg.image,
                                    ivarraw=sciImg.ivar,
                                    skymodel=skymodel,
                                    objmodel=objmodel,
                                    ivarmodel=ivarmodel,
                                    scaleimg=objFind.scaleimg,
                                    waveimg=waveimg,
                                    tilts=tilts,
                                    bpmmask=outmask,
                                    detector=sciImg.detector,
                                    slits=slits,
                                    wavesol=None,
                                    maskdef_designtab=None,
                                    sci_spat_flexure=sciImg.spat_flexure,
                                    sci_spec_flexure=None,
                                    vel_corr=None,
                                    vel_type=None)

    spec2DObj.process_steps = sciImg.process_steps

    # QA
    spec2DObj.gen_qa()

    return spec2DObj, sobjs




def jwst_extract_subimgs(final_slit, intflat_slit):

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
    # Sometimes this fails at the 1e-4 level and disagreess about nans???
    #assert np.allclose(waveimg, waveimg_from_wcs, rtol=1e-3, atol=1e-3, equal_nan=True)


    flatfield = np.array(intflat_slit.data.T, dtype=float) #if intflat_slit is not None else np.ones_like(pathloss)
    pathloss = np.array(final_slit.pathloss_uniform.T, dtype=float) if final_slit.source_type == 'EXTENDED' else \
        np.array(final_slit.pathloss_point.T, dtype=float)
    if pathloss.shape == (0,0):
        msgs.warn('No pathloss for slit {0}'.format(slit_name) + ', setting to 1.0')
        pathloss = np.ones_like(flatfield)

    barshadow = np.array(final_slit.barshadow.T, dtype=float)
    if barshadow.shape == (0,0):
        msgs.warn('No barshadow for slit {0}'.format(slit_name) + ', setting to 1.0')
        barshadow = np.ones_like(flatfield)

    photom_conversion = final_slit.meta.photometry.conversion_megajanskys
    final = np.array(final_slit.data.T, dtype=float)


    # Generate some tilts and a spatial image
    finitemask = np.isfinite(waveimg)
    # Get slit bounadries
    slit_left, slit_righ = jwst_get_slits(finitemask)

    waveimg = 1e4*waveimg
    #waveimg[np.logical_not(finitemask)] = 0.0
    wave_min, wave_max = np.min(waveimg[finitemask]), np.max(waveimg[finitemask])

    tilts = np.zeros_like(waveimg)
    tilts[finitemask] = (waveimg[finitemask] - wave_min) / (wave_max - wave_min)

    # TODO Fix this spat_pix to make it increasing with pixel. For now don't use it
    # This currnetly depends on poisition angle which I need to hack to fix
    # ra_min, ra_max = np.min(ra_sub[finitemask_sub]),  np.max(ra_sub[finitemask_sub])
    # spat_pix_sub = np.zeros_like(ra_sub)
    # spat_pix_sub[finitemask_sub] = spat_lo + (ra[finitemask_sub] - ra_min) / (ra_max - ra_min) * (nspat_sub - 1)



    ########################
    # The image segment being used for each slit
    spec_lo = final_slit.xstart - 1
    spec_hi = spec_lo + final_slit.xsize
    spat_lo = final_slit.ystart - 1
    spat_hi = spat_lo + final_slit.ysize
    # slice object for the segment
    slit_slice = np.s_[spec_lo: spec_hi, spat_lo: spat_hi]

    #embed()
    #rate = np.array(e2d_slit.data.T, dtype=float)
    #rate_var_rnoise = np.array(e2d_slit.var_rnoise.T, dtype=float)
    #rate_var_poisson = np.array(e2d_slit.var_poisson.T, dtype=float)
    # This is currently buggy as it includes flat field error
    #rate_var_tot = np.square(np.array(e2d_slit.err.T, dtype=float))
    dq = np.array(final_slit.dq.T, dtype=int)

    slit_left_orig = spat_lo + slit_left
    slit_righ_orig = spat_lo + slit_righ
    spec_vals_orig = spec_lo + np.arange(spec_hi - spec_lo)


    return slit_slice, slit_left, slit_righ, slit_left_orig, slit_righ_orig, spec_vals_orig, src_trace_ra, src_trace_dec, dq, \
           ra, dec, finitemask, waveimg, tilts, flatfield, pathloss, barshadow, photom_conversion, final


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


    return science, sciivar, gpm, base_var, count_scale, tilts, waveimg, finite_mask, slit_left, slit_righ, t_eff




def jwst_extract_subimgs_old(e2d_slit, final_slit, intflat_slit):

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
    # Sometimes this fails at the 1e-4 level and disagreess about nans???
    #assert np.allclose(waveimg, waveimg_from_wcs, rtol=1e-3, atol=1e-3, equal_nan=True)


    flatfield = np.array(intflat_slit.data.T, dtype=float) #if intflat_slit is not None else np.ones_like(pathloss)
    pathloss = np.array(final_slit.pathloss_uniform.T, dtype=float) if final_slit.source_type == 'EXTENDED' else \
        np.array(final_slit.pathloss_point.T, dtype=float)
    if pathloss.shape == (0,0):
        msgs.warn('No pathloss for slit {0}'.format(slit_name) + ', setting to 1.0')
        pathloss = np.ones_like(flatfield)

    barshadow = np.array(final_slit.barshadow.T, dtype=float)
    if barshadow.shape == (0,0):
        msgs.warn('No barshadow for slit {0}'.format(slit_name) + ', setting to 1.0')
        barshadow = np.ones_like(flatfield)

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

