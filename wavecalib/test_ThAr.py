""" Run the current, favored holy grail routine on a set of
standard input spectra.  Generate summary output
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import time
import glob
import astropy.io.fits as pyfits
import numpy as np
import warnings
import pdb

from arclines.holy import grail

#from xastropy.xutils import xdebug as xdb

import arclines
test_arc_path = arclines.__path__[0]+'/data/test_arcs/'
outdir = 'TEST_SUITE_OUTPUT/'


def tst_holy(name, spec, wav_id, pix_id, test='semi_brute', toler=0.001):
    """
    toler : float
      Tolerance for scoring correct matches
    """
    # Favored parameters (should match those in the defaults)
    siglev=20.
    min_ampl=1000.
    min_match = 10
    lines = ['ThAr']

    # Run
    outroot = outdir+name
    #print(wav_id)
    if test == 'general':
        best_dict, final_fit = grail.general(spec, lines, siglev=siglev,
                                             min_ampl=min_ampl, min_nmatch=min_match, outroot=outroot)
        if best_dict is None:
            return "FAILED", None, None
    else:
        pdb.set_trace()

    # Score
    grade = 'PASSED'
    nfail = 0
    for pp in range(pix_id.size):
        wid = np.argmin(np.abs(pix_id[pp]-final_fit['xfit']))
        if np.abs(wav_id[pp]-final_fit['yfit'][wid]) > toler:
            nfail += 1
    if nfail != 0:
        grade = "FAILED"
    # Warn
    warnings.warn("Solution for {:s} failed {:d}/{:d} lines".format(name, nfail, pix_id.size))
    return grade, best_dict, final_fit


def main(flg_tst):

    if flg_tst in [1]:
        test = 'semi_brute'
    elif flg_tst in [2]:
        test = 'general'

    # Initialize
    names = []
    specs = []
    wavid = []
    pixid = []
    if True:
        # Test all HIREDUX solutions
        files = glob.glob(test_arc_path + "HIREDUX/*aspec.fits.gz")
        for fn in files:
            names += [fn.split("/")[-1]]
            filename = pyfits.open(fn)
            fx = filename[0].data
            ids = pyfits.open(fn.replace("aspec.fits.gz", "lines.fits.gz"))
            pxs = ids[1].data["PIX"]
            wvs = ids[1].data["WV"]
            ordspec, ordwavs, ordpixs = [], [], []
            for ord in range(fx.shape[0]):
                # Store the spectrum
                ordspec += [fx[ord, :]]
                # Store the 'true' line IDs
                ww = np.where(wvs[ord, :] != 0.0)
                ordwavs += [wvs[ord, ww[0]]]
                ordpixs += [pxs[ord, ww[0]]]
            specs += [ordspec]
            wavid += [ordwavs]
            pixid += [ordpixs]

    # Run it
    sv_grade = []  # for the end, just in case
    for name, spec, wvid, pxid in zip(names, specs, wavid, pixid):
        grades = np.zeros(len(spec), dtype=np.int)
        bwave, bdisp = np.zeros(len(spec)), np.zeros(len(spec))
        timstart = time.time()
        for ord in range(len(spec)):
            print("Analyzing order {0:d}/{1:d}".format(ord, fx.shape[0]))
            grade, best_dict, final_fit = tst_holy(name, spec[ord], wvid[ord], pxid[ord], test=test)
            if best_dict is None:
                grades = grades[:ord-1]
                break
            bwave[ord] = best_dict['bwv']
            bdisp[ord] = best_dict['bdisp']
            if grade == 'PASSED':
                grades[ord] = 1
        timfinish = time.time()
        sv_grade.append(grades.copy())
        pdb.set_trace()
        print("Completion time (minutes):", (timfinish-timstart)/60.0)
        from matplotlib import pyplot as plt
        plt.plot(bwave, bdisp, 'bo')
        plt.show()

    # Report it
    print('==============================================================')
    for name, grade in zip(names, sv_grade):
        print("{:s} :: PASSED {:d}/{:d}".format(name, np.sum(grade), grade.size))


# Test
if __name__ == '__main__':
    # flg_tst = 1   # Run em all with semi-brute
    flg_tst = 2   # Run em all with general

    main(flg_tst)
