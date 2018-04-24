from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from astropy.io import fits
import glob

from pypit import msgs
from pypit import ardebug as debugger
from pypit import ginga
from pypit import arload
from pypit import arproc
from pypit import arcomb

debug = debugger.init()
debug['develop'] = True
msgs.reset(debug=debug, verbosity=2)

from pypit import artrace

def_settings=dict(trace={'slits': {'single': [],
                                   'function': 'legendre',
                                   'polyorder': 3,
                                   'diffpolyorder': 2,
                                   'fracignore': 0.01,
                                   'number': -1,
                                   'maxgap': None,
                                   'sigdetect': 20.,
                                   'pca': {'params': [3,2,1,0,0,0], 'type': 'pixel', 'extrapolate': {'pos': 0, 'neg':0}},
                                   'sobel': {'mode': 'nearest'}},
                         'combine': {'match': -1.,
                                     'satpix': 'reject',
                                     'method': 'weightmean',
                                     'reject': {'cosmics': 20., 'lowhigh': [0,0], 'level': [3.,3.], 'replace': 'maxnonsat'}}})

def pypit_trace_slits():
    pass


def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Developing/testing/checking trace_slits [v1.1]')
    parser.add_argument("instr", type=str, help="Instrument [deimos]")
    parser.add_argument("--det", default=1, type=int, help="Detector")
    parser.add_argument("--show", default=False, action="store_true", help="Show the image with traces")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):

    xgap = 0.0   # Gap between the square detector pixels (expressed as a fraction of the x pixel size)
    ygap = 0.0   # Gap between the square detector pixels (expressed as a fraction of the x pixel size)
    ysize = 1.0  # The size of a pixel in the y-direction as a multiple of the x pixel size

    settings = def_settings.copy()
    if pargs.instr == 'deimos':
        spectrograph = 'keck_deimos'
        saturation = 65535.0              # The detector Saturation level
        files = glob.glob('data/DE*')
        #files = ['../RAW_DATA/Keck_DEIMOS/830G_L/'+ifile for ifile in ['d0914_0014.fits', 'd0914_0015.fits']]
        numamplifiers=1

        # Grab data info
        datasec, oscansec, naxis0, naxis1 = arproc.get_datasec(spectrograph, files[0],
                                                               numamplifiers=numamplifiers, det=pargs.det)

        # Load + process images
        frames = []
        for ifile in files:
            rawframe, head0 = arload.load_raw_frame(spectrograph, ifile, pargs.det, disp_dir=0)
            # Bias subtract
            newframe = arproc.sub_overscan(rawframe.copy(), numamplifiers, datasec, oscansec)
            # Trim
            frame = arproc.trim(newframe, numamplifiers, datasec)
            # Append
            frames.append(frame)

        # Convert to array
        frames_arr = np.zeros((frames[0].shape[0], frames[0].shape[1], len(frames)))
        for ii in range(len(frames)):
            frames_arr[:,:,ii] = frames[ii]

        # Combine
        mstrace = arcomb.comb_frames(frames_arr, frametype='Unknown',
                                     method=settings['trace']['combine']['method'],
                                     reject=settings['trace']['combine']['reject'],
                                     satpix=None, saturation=saturation)
        # binpx
        binbpx = np.zeros_like(mstrace)

        #hdul = fits.open('trace_slit.fits')
        settings['trace']['slits']['sigdetect'] = 50.0
        settings['trace']['slits']['pca']['params'] = [3,2,1,0]
    else:
        debugger.set_trace()

    # pixlocn
    pixlocn = artrace.gen_pixloc(mstrace, xgap, ygap, ysize)

    #binbpx = hdul[0].data
    #pixlocn = hdul[1].data
    #mstrace = hdul[2].data

    lcenint, rcenint, extrapord = artrace.refactor_trace_slits(pargs.det, mstrace, binbpx, pixlocn,
                                         settings=settings, pcadesc="", maskBadRows=False, min_sqm=30.)
    if pargs.show:
        viewer, ch = ginga.show_image(mstrace)
        nslit = lcenint.shape[1]
        ginga.show_slits(viewer, ch, lcenint, rcenint, np.arange(nslit) + 1, pstep=50)

    debugger.set_trace()

if __name__ == '__main__':
    pargs = parser()
    main(pargs)
