from __future__ import (print_function, absolute_import)
from __future__ import (division, unicode_literals)

# importing

import numpy as np
from astropy.io import fits
from glob import glob
import os
import argparse

from pypeit import ginga
from pypeit.core import pixels
from pypeit import traceslits
from pypeit import processimages
from pypeit import scienceimage
from pypeit import arcimage
from pypeit.core import arc
from pypeit import wavecalib
#from pypit import wavecalib
from pypeit import wavetilts
from pypeit import waveimage
from pypeit import flatfield
from pypeit import traceimage
from pypeit import specobjs
from pypeit import utils

# Spectrgraph and Settings
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import pypeitpar
from pypeit import ginga
from pypeit.core.skysub import global_skysub

from linetools import utils as ltu

from IPython import embed
from matplotlib import pyplot as plt

# Load Spectrograph
spectro_name = 'keck_nires'
spectrograph = load_spectrograph(spectrograph=spectro_name)

### Load settings

## Detector settings
par = spectrograph.default_pypeit_par()

## Trace settings
par['calibrations']['traceframe']['process']['overscan'] = 'median'
par['calibrations']['slits']['polyorder'] = 5
par['calibrations']['slits']['maxshift'] = 3.
par['calibrations']['slits']['pcatype'] = 'order'

# Arc settings
par['calibrations']['arcframe']['process']['overscan'] = 'median'

# Flat settings
par['calibrations']['pixelflatframe']['process']['overscan'] = 'median'

# Tilts settings
par['calibrations']['tilts']['tracethresh'] = [50, 50, 60, 60, 2000]

# Science settings
par['scienceframe']['process']['overscan'] = 'median'


def get_tslits_nires(flat_files,
                     user_settings=par,
                     gingashow=True,
                     tilt_root='tilt_nires'):
    """Precess flat files and get titlts for NIRES
    """

    # Process flat images
    tImage = traceimage.TraceImage(spectrograph,
                                   file_list=flat_files,
                                   par=par['calibrations']['traceframe'])

    tflat = tImage.process(bias_subtract='overscan',
                           trim=False)

    mstrace = tflat.copy()

    # Define pixlocn and bpm
    pixlocn = pixels.gen_pixloc(tImage.stack.shape)
    bpm = spectrograph.bpm(shape=tflat.shape, det=1)

    # Instantiate Trace
    tSlits = traceslits.TraceSlits(mstrace,
                                   pixlocn,
                                   par=par['calibrations']['slits'],
                                   binbpx=bpm)
    tslits_dict = tSlits.run(plate_scale = 0.123)

    if gingashow:
        # Look at what TraceSlits was actually trying to trace
        viewer, ch = ginga.show_image(tSlits.edgearr)
        # Look at the sawtooth convolved image
        viewer, ch = ginga.show_image(tSlits.siglev)

        tmp = tSlits.edgearr * 100.
        tmp[np.where(tmp == 0.)] = 1.
        ginga.show_image(tSlits.mstrace * tmp)
        ginga.show_slits(viewer,
                         ch,
                         tSlits.lcen,
                         tSlits.rcen,
                         slit_ids=np.arange(tSlits.lcen.shape[1]) + 1,
                         pstep=50)

    if tilt_root is not None:
        # Write dict on a json file
        jdict = ltu.jsonify(tslits_dict.copy())
        ltu.savejson(tilt_root + '.json', jdict, overwrite=True, indent=None, easy_to_read=True)
        print("Wrote: {:s}".format(tilt_root + '.json'))

    return tslits_dict


EXAMPLES = """

    dev_nires_quicklook.py -f s180604_0023.fits.gz -a s180604_0107.fits.gz -b s180604_0108.fits.gz

    dev_nires_quicklook.py -a s180604_0107.fits.gz -b s180604_0108.fits.gz

    Calling sequence from the dev suite: 

    python dev_nires_quicklook.py -f 
    /Users/joe/python/PypeIt-development-suite/RAW_DATA/Keck_NIRES/s180604_0023.fits.gz 
    /Users/joe/python/PypeIt-development-suite/RAW_DATA/Keck_NIRES/s180604_0024.fits.gz 
    /Users/joe/python/PypeIt-development-suite/RAW_DATA/Keck_NIRES/s180604_0025.fits.gz 
    -a /Users/joe/python/PypeIt-development-suite/RAW_DATA/Keck_NIRES/s180604_0089.fits.gz 
    -b /Users/joe/python/PypeIt-development-suite/RAW_DATA/Keck_NIRES/s180604_0090.fits.gz
    """


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='''
            Quick look extraction of NIRES spectrum.
            TO DO:
                Currently only works for single AB sequence!
            TWO STEPS:
                1, Generate tslit if flatfiles is not None
                2, Reduce a sequence of science files in A and B positions.
                      ''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('-f', '--flatfiles', nargs='+', type=str, default=None,
                        help='A list of flat images')

    parser.add_argument('-a', '--sciAfiles', nargs='+', type=str, default=None,
                        required=True, help='A list of science images at A position')

    parser.add_argument('-b', '--sciBfiles', nargs='+', type=str, default=None,
                        required=True, help='A list of science images at B position')

    parser.add_argument('--show', default=True, type=bool,
                        help='Need to open ginga with ginga --module=RC first')
    parser.add_argument('--version', action='version', version='1.1')

    return parser.parse_args()


############################################

if __name__ == '__main__':
    args = parse_arguments()

    if os.path.exists('mkdir -p QA/PNGs/'):
        pass
    else:
        os.system('mkdir -p QA/PNGs/')

    # ToDo -- fix ifs
    if args.show:
        gingashow = True
    ## Get tslits_dict either from flatfiles or load from a json file
    if args.flatfiles is None:
        jdict = ltu.loadjson('tilt_nires.json')
        tslits_dict = jdict.copy()
        for tkey in tslits_dict.keys():
            tslits_dict[tkey] = np.array(tslits_dict[tkey])
    else:
        tslits_dict = get_tslits_nires(args.flatfiles, user_settings=par, gingashow=gingashow)

    # Get Tilts from scienceB image
    aImage = arcimage.ArcImage(spectrograph,
                               file_list=args.sciBfiles,
                               par=par['calibrations']['arcframe'])
    msarc = aImage.process(bias_subtract='overscan',
                           trim=False)
    pixlocn = pixels.gen_pixloc(aImage.stack.shape)
    bpm = spectrograph.bpm(shape=msarc.shape, det=1)
    # Extract spectrum at the center
    arccen, mask, _ = arc.get_censpec(tslits_dict['lcen'],
                                      tslits_dict['rcen'],
                                      pixlocn,
                                      msarc,
                                      1)
    waveTilts = wavetilts.WaveTilts(msarc, spectrograph=spectrograph,
                                    par=par['calibrations']['tilts'], det=1,
                                    tslits_dict=tslits_dict, pixlocn=pixlocn)
    nslits = tslits_dict['lcen'].shape[1]
    maskslits = np.zeros(nslits, dtype=bool)
    # QA is not working here
    mstilts, wt_maskslits = waveTilts.run(maskslits=maskslits,
                                          wv_calib=None,
                                          doqa=False)

    # Wavelength calibration from scienceB image
    arcparam = {}
    spectrograph.setup_arcparam(arcparam)
    setup = 'NIRES'
    waveCalib = wavecalib.WaveCalib(msarc,
                                    spectrograph=spectrograph,
                                    #par=par['calibrations']['wavelengths'],
                                    det=1,
                                    setup=setup,
                                    arcparam=arcparam)

    wv_calib, _ = waveCalib.run(tslits_dict['lcen'],
                                tslits_dict['rcen'],
                                pixlocn,
                                nonlinear=spectrograph.detector[0]['nonlinear'])
    # Get 2-D wavelength solution
    waveImage = waveimage.WaveImage(tslits_dict['slitpix'],
                                    mstilts,
                                    wv_calib,
                                    setup=setup,
                                    maskslits=maskslits)
    waveimg = waveImage._build_wave()

    # Process science images
    datasec_img = np.zeros_like(msarc) + 1

    # Hack, create the fitstbl from the header of one of the science frames
    # ToDo -- this step should not be necessary
    #from astropy.io import fits
    #from astropy.table import Table
    #hdul = fits.open(args.sciAfiles[0])
    #fitstbl = Table(rows=[(hdul[0].header['ITIME'], '1x1')], names=('exptime', 'binning',))

    # Processing A
    sciIMG_A = scienceimage.ScienceImage(spectrograph,
                                         file_list=args.sciAfiles,
                                         frame_par=par['scienceframe'],
                                         #tslits_dict=tslits_dict,
                                         #tilts=mstilts,
                                         det=1,
                                         #datasec_img=datasec_img,
                                         #bpm=bpm,
                                         #pixlocn=pixlocn,
                                         #fitstbl=fitstbl)
                                         )
    sciframe_A, sciivar_A, rawvarframe_A, crmask_A = sciIMG_A.process(bias_subtract=None,
                                                           pixel_flat=None,
                                                           bpm=bpm,
                                                           apply_gain=True,
                                                           trim=False)
    # Processing B
    sciIMG_B = scienceimage.ScienceImage(spectrograph,
                                         file_list=args.sciBfiles,
                                         frame_par=par['scienceframe'],
                                         #tslits_dict=tslits_dict,
                                         #tilts=mstilts,
                                         det=1,
                                         #datasec_img=datasec_img,
                                         #bpm=bpm,
                                         #pixlocn=pixlocn,
                                         #fitstbl=fitstbl)
                                         )
    sciframe_B, sciivar_B, rawvarframe_B, crmask_B = sciIMG_B.process(bias_subtract=None,
                                                           pixel_flat=None,
                                                           bpm=bpm,
                                                           apply_gain=True,
                                                           trim=False)
    # This addition of a tiny bit of effor imposes a maximum S/N ratio which stabilizes the fits
    adderr = 0.01  # Additional error to add to the formal errors, as a fraction of the flux; default to 0.01 (1 per cent).
    ivar_A = utils.calc_ivar(rawvarframe_A) * (crmask_A == False)
    gmask_A = (ivar_A > 0.0).astype(int)  # = 1 for good points, 0 for bad points
    ivar_A = gmask_A / (1.0 / (ivar_A + (1.0 - gmask_A)) + (adderr ** 2) * (np.abs(sciframe_A)) ** 2)
    ivar_B = utils.calc_ivar(rawvarframe_B) * (crmask_B == False)
    gmask_B = (ivar_B > 0.0).astype(int)  # = 1 for good points, 0 for bad points
    ivar_B = gmask_B / (1.0 / (ivar_B + (1.0 - gmask_B)) + (adderr ** 2) * (np.abs(sciframe_B)) ** 2)
    # Diff the images
    diff_AB = sciframe_A - sciframe_B
    var_A = utils.calc_ivar(ivar_A)
    var_B = utils.calc_ivar(ivar_B)
    var_AB = var_A + var_B
    # I think processimages already masks saturation but am not sure
    mask_AB = (sciframe_A > 0.0) & (sciframe_B > 0.0) & \
              (sciframe_A < spectrograph.detector[0]['saturation']) & \
              (sciframe_B < spectrograph.detector[0]['saturation']) & \
              (crmask_A == False) & (crmask_B == False)
    ivar_AB = utils.calc_ivar(var_AB) * (mask_AB == True)

    # Sky subtraction
    from pypeit.core import skysub
    from pypeit.core import extract

    slitpix = tslits_dict['slitpix']
    lcen = tslits_dict['lcen']
    rcen = tslits_dict['rcen']
    ximg = tslits_dict['ximg']
    edgmask = tslits_dict['edge_mask']
    FWHM = 5.0
    bsp = 0.8
    residual_img = np.zeros_like(sciframe_A)
    specobjs_pos = specobjs.SpecObjs()
    specobjs_neg = specobjs.SpecObjs()
    skymask = (slitpix > 0)
    for islit in range(1, nslits + 1):
        thismask = (slitpix == islit)
        residual_img[thismask] = skysub.global_skysub(diff_AB,
                                                      ivar_AB,
                                                      mstilts,
                                                      thismask,
                                                      lcen[:,islit-1],
                                                      rcen[:,islit-1],
                                                      inmask=((edgmask == False) & (mask_AB == True)),
                                                      bsp=bsp,
                                                      pos_mask=False,
                                                      show_fit=False)
        image = diff_AB - residual_img
        # Extract negative trace
        specobj_slit_neg, skymask_neg, objmask_neg = extract.objfind(-image,
                                                                     #ivar_AB,
                                                                     thismask,
                                                                     lcen[:,islit-1],
                                                                     rcen[:,islit-1],
                                                                     sig_thresh=3.0,
                                                                     inmask=mask_AB,
                                                                     FWHM=FWHM,
                                                                     nperslit=1,
                                                                     trim_edg=(3, 3),
                                                                     show_trace=False,
                                                                     show_peaks=False,
                                                                     show_fits =False)
        if specobj_slit_neg is not None:
            specobjs_neg.add_sobj(specobj_slit_neg.specobjs.tolist())
        # Extract positive trace
        specobj_slit_pos, skymask_pos, objmask_pos = extract.objfind(image,
                                                                     #ivar_AB,
                                                                     thismask,
                                                                     lcen[:,islit-1],
                                                                     rcen[:,islit-1],
                                                                     sig_thresh=3.0,
                                                                     inmask=mask_AB,
                                                                     FWHM=FWHM,
                                                                     nperslit=1,
                                                                     trim_edg=(3, 3),
                                                                     show_trace=False,
                                                                     show_peaks=False,
                                                                     show_fits =False)
        if specobj_slit_pos is not None:
            specobjs_pos.add_sobj(specobj_slit_pos.specobjs.tolist())
        skymask[thismask] = (skymask_pos & skymask_neg)
    # Show results on ginga
    if gingashow:
        # Plot the chi image
        chi = (diff_AB - residual_img) * np.sqrt(ivar_AB) * (slitpix > 0) * ((edgmask == False) & (mask_AB == True))
        viewer, ch = ginga.show_image(chi)
        ginga.show_slits(viewer, ch, lcen, rcen, slit_ids=None)
        for islit in range(0, nslits):
            ginga.show_trace(viewer, ch, specobjs_pos[islit].trace_spat, trc_name=specobjs_pos[islit].idx, color='blue')
            ginga.show_trace(viewer, ch, specobjs_neg[islit].trace_spat, trc_name=specobjs_neg[islit].idx,
                             color='orange')

    # Boxcar extraction
    from pypeit.core.extract import extract_boxcar

    outmask = (slitpix > 0) * ((edgmask == False) & (mask_AB == True))
    box_rad = 8.0
    # ToDo -- Check for indexes in islit [0-based or 1-based?]

    for islit in range(1, nslits + 1):
        # Positive trace
        flux = extract_boxcar(image,
                              specobjs_pos[islit - 1].trace_spat,
                              box_rad)
        mvarimg = 1.0 / (ivar_AB + (ivar_AB == 0))
        mvar_box = extract.extract_boxcar(mvarimg,
                                          specobjs_pos[islit - 1].trace_spat,
                                          box_rad,
                                          ycen=specobjs_pos[islit - 1].trace_spec)
        pixtot = extract_boxcar(0 * mvarimg + 1.0,
                                specobjs_pos[islit - 1].trace_spat, box_rad, ycen=specobjs_pos[islit - 1].trace_spec)
        mask_box = (extract_boxcar(~outmask, specobjs_pos[islit - 1].trace_spat, box_rad,
                                   ycen=specobjs_pos[islit - 1].trace_spec) != pixtot)
        box_denom = extract_boxcar(waveimg > 0.0, specobjs_pos[islit - 1].trace_spat, box_rad,
                                   ycen=specobjs_pos[islit - 1].trace_spec)
        wave = extract_boxcar(waveimg, specobjs_pos[islit - 1].trace_spat, box_rad,
                              ycen=specobjs_pos[islit - 1].trace_spec) / (
                       box_denom + (box_denom == 0.0))
        fluxivar = mask_box / (mvar_box + (mvar_box == 0.0))
        specobjs_pos[islit - 1].boxcar['wave'] = wave
        specobjs_pos[islit - 1].boxcar['flux'] = flux
        specobjs_pos[islit - 1].boxcar['var'] = 1. / fluxivar
        # Negative trace
        flux = extract_boxcar(-image,
                              specobjs_neg[islit - 1].trace_spat,
                              box_rad)
        mvarimg = 1.0 / (ivar_AB + (ivar_AB == 0))
        mvar_box = extract.extract_boxcar(mvarimg,
                                          specobjs_neg[islit - 1].trace_spat,
                                          box_rad,
                                          ycen=specobjs_neg[islit - 1].trace_spec)
        pixtot = extract_boxcar(0 * mvarimg + 1.0,
                                specobjs_neg[islit - 1].trace_spat, box_rad, ycen=specobjs_neg[islit - 1].trace_spec)
        mask_box = (extract_boxcar(~outmask, specobjs_neg[islit - 1].trace_spat, box_rad,
                                   ycen=specobjs_neg[islit - 1].trace_spec) != pixtot)
        box_denom = extract_boxcar(waveimg > 0.0, specobjs_pos[islit - 1].trace_spat, box_rad,
                                   ycen=specobjs_pos[islit - 1].trace_spec)
        wave = extract_boxcar(waveimg, specobjs_neg[islit - 1].trace_spat, box_rad,
                              ycen=specobjs_neg[islit - 1].trace_spec) / (
                       box_denom + (box_denom == 0.0))
        fluxivar = mask_box / (mvar_box + (mvar_box == 0.0))
        specobjs_neg[islit - 1].boxcar['wave'] = wave
        specobjs_neg[islit - 1].boxcar['flux'] = flux
        specobjs_neg[islit - 1].boxcar['var'] = 1. / fluxivar

    # Show spectrum
    import matplotlib.pyplot as plt

    plt.close()
    for islit in range(1, nslits + 1):
        plt.plot(specobjs_neg[islit - 1].boxcar['wave'], color='b')
        plt.plot(specobjs_pos[islit - 1].boxcar['wave'], color='red')
    plt.legend()
    plt.title('Check Wavelength Solution')
    plt.xlabel('Pixels')
    plt.ylabel('Wavelength')
    plt.savefig('Wave.png', bbox_inches='tight')

    import matplotlib.pyplot as plt

    plt.close()
    for islit in range(1, nslits + 1):
        plt.plot(specobjs_neg[islit - 1].boxcar['wave'], specobjs_neg[islit - 1].boxcar['flux'], color='b')
        #plt.plot(specobjs_neg[islit - 1].boxcar['wave'], np.sqrt(specobjs_neg[islit - 1].boxcar['var']), color='b',
        #         alpha=0.5)
        plt.plot(specobjs_pos[islit - 1].boxcar['wave'], specobjs_pos[islit - 1].boxcar['flux'], color='red')
        #plt.plot(specobjs_pos[islit - 1].boxcar['wave'], np.sqrt(specobjs_pos[islit - 1].boxcar['var']), color='red',
        #         alpha=0.5)
    plt.legend()
    plt.title('Boxcar extraction')
    plt.xlabel('Wavelength')
    plt.ylabel('Counts')
    plt.savefig('Extract.png', bbox_inches='tight')
    plt.show()

    """
    # Write dict on a json file
    jdict_pos = ltu.jsonify(specobjs_pos.tolist())
    jdict_neg = ltu.jsonify(specobjs_neg.tolist())
    ltu.savejson('Spectrum_pos.json', jdict_pos, overwrite=True, indent=None, easy_to_read=True)
    ltu.savejson('Spectrum_neg.json', jdict_neg, overwrite=True, indent=None, easy_to_read=True)
    print("Wrote: {:s}".format('Spectrum_pos.json'))
    print("Wrote: {:s}".format('Spectrum_neg.json'))
    """
