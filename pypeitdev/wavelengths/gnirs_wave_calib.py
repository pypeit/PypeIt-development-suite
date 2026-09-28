

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.table import Table
from pypeit.core import wavecal
from pypeit.core import qa
from pypeit import wavecalib
import scipy
import os

def wv_calib_from_extern(wave_soln, arc, lamps, outfile = None, sigdetect=5.0, fwhm=4.0, nonlinear_counts=1e10, outroot='./', debug=False):
    """

    Args:
        wave_soln:
        arc:
        lamps:
        outfile:
        sigdetect:
        fwhm:
        nonlinear_counts:
        outroot:
        debug:

    Returns:

    """

    # TODO add array size checking etc.
    nslits = wave_soln.shape[1]
    nspec  = wave_soln.shape[0]
    line_lists = wavecal.waveio.load_line_lists(lamps)
    wv_calib = {}
    spec_vec = np.arange(nspec)
    for islit in range(nslits):
        print(str(islit))
        # Find peaks for this slit
        tcent, ecent, cut_tcent, icut, spec_cont_sub = wavecal.wvutils.arc_lines_from_spec(arc[:,islit], sigdetect=sigdetect,
                                                                                   nonlinear_counts=nonlinear_counts,
                                                                                   fwhm=fwhm, debug=debug)
        detections = tcent[icut]
        wave_det = (scipy.interpolate.interp1d(spec_vec, wave_soln[:, islit], kind='cubic'))(detections)
        patt_dict = {}
        patt_dict['mask'] = np.ones_like(detections, dtype=bool)
        patt_dict['IDs'] = wave_det
        patt_dict['bdisp'] = (np.max(wave_soln[:,islit]) - np.min(wave_soln[:,islit]))/nspec
        final_fit = wavecal.fitting.fit_slit(spec_cont_sub, patt_dict, detections, line_lists, vel_tol=300.0, outroot=outroot, verbose=True)
        wv_calib[str(islit)] = final_fit
        qa_file = qa.set_qa_filename('GNIRS', 'arc_fit_qa', slit=islit, out_dir=outroot)
        wavecal.qa.arc_fit_qa(wv_calib[str(islit)],outfile=qa_file)
        if debug:
            # Show the QA
            wavecal.qa.arc_fit_qa(wv_calib[str(islit)])

    waveCalib = wavecalib.WaveCalib(None,None)
    if outfile is not None:
        waveCalib.save_master(wv_calib, outfile=outfile)

    return wv_calib



#scifile ='/Users/joe/GN-LP-7/redux/J1335+3533/Science/J1335+3533_10/sci-cN20160127S0399-402.fits'
scifile ='/Users/joe/REDUX/gnirs_redux/1420+0227/Science/1420+0227_1/sci-N20110622S0304-307.fits'

obj = Table.read(scifile,hdu=4)

nspec = obj[0]['WAVE_OPT'].size
norders = int(len(obj)/2)
wave = np.zeros((nspec, norders))
arcspec = np.zeros((nspec, norders))
for iord, igem in zip(range(norders),range(2*norders-2,-2,-2)):
    wave[:,iord] = (obj[igem]['WAVE_OPT'])[::-1]
    arcspec[:,iord]  = (obj[igem]['SKY_OPT'])[::-1]

lamps = ['OH_GNIRS']
fwhm = 4.0

wv_calib = wv_calib_from_extern(wave,arcspec, lamps,debug=False, outfile='./gemini_gnirs_idl.json')