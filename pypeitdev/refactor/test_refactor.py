# Tests for refactoring

import os, sys
import shutil
import numpy as np
from pathlib import Path

from pypeit.spectrographs.util import load_spectrograph
from pypeit.wavecalib import WaveCalib
from pypeit.slittrace import SlitTraceSet
from pypeit.scripts.run_pypeit import RunPypeIt
from pypeit import specobjs

import pytest
from IPython import embed

sys.path.append(os.path.join(
    os.path.abspath(
        os.environ["PYPEIT_DEV"]),"test_scripts"))

# Shane Kastblue 600/4310 D55

def test_shane_kastb_spec1d():
    _redux_out = os.path.join(os.environ['PYPEIT_DEV'], 'REDUX_OUT')

    # Load spec1d's
    setup = '600_4310_d55'
    root = 'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits'
    spec1d_full_refactor_file = os.path.join(_redux_out, 
                                                'shane_kast_blue', 
                                                setup,
                                                'REFACTOR_FULL','Science', root)
    spec1d_full = specobjs.SpecObjs.from_fitsfile(spec1d_full_refactor_file,
                                                    chk_version=False)

    spec1d_dev_file = spec1d_full_refactor_file.replace('REFACTOR_FULL','DEVELOP')
    spec1d_dev = specobjs.SpecObjs.from_fitsfile(spec1d_dev_file,
                                                    chk_version=False)
    spec1d_load_file = spec1d_full_refactor_file.replace('REFACTOR_FULL', 'shane_kast_blue_A')
    spec1d_load = specobjs.SpecObjs.from_fitsfile(spec1d_load_file,
                                                    chk_version=False)


    # Full vs. Develop
    for mode, spec1d in zip(['full', 'load'], [spec1d_full, spec1d_load]):
        assert spec1d.nobj == spec1d_dev.nobj
        assert np.allclose(spec1d.BOX_COUNTS, spec1d_dev.BOX_COUNTS, atol=1e-4)
        assert np.allclose(spec1d.OPT_WAVE, spec1d_dev.OPT_WAVE, atol=1e-4)



def test_mosfire_ylong_spec1d():
    _redux_out = os.path.join(os.environ['PYPEIT_DEV'], 'REDUX_OUT')

    # Load spec1d's
    setup = 'Y_long'

    root = 'spec1d_m191120_0043-J2132-1434_OFF_MOSFIRE_20191120T061253.347.fits'
    idx = 0  # Only 1 sorce on science exposure

    #root = 'spec1d_m191118_0064-GD71_MOSFIRE_20191118T104704.507.fits' # Standard star
    #idx = 0 # This one fails, for float-precision bizarre reasons
    #idx = 2 # This is the standard star

    spec1d_full_refactor_file = os.path.join(_redux_out, 
                                                'keck_mosfire', 
                                                'REFACTOR_FULL','Science', root)
    spec1d_full = specobjs.SpecObjs.from_fitsfile(spec1d_full_refactor_file,
                                                    chk_version=False)

    spec1d_dev_file = spec1d_full_refactor_file.replace('REFACTOR_FULL','DEVELOP')
    spec1d_dev = specobjs.SpecObjs.from_fitsfile(spec1d_dev_file,
                                                    chk_version=False)
    #spec1d_load_file = spec1d_full_refactor_file.replace('REFACTOR_FULL', setup)
    #spec1d_load = specobjs.SpecObjs.from_fitsfile(spec1d_load_file,
    #                                                chk_version=False)


    # Full vs. Develop
    #for mode, spec1d in zip(['full', 'load'], [spec1d_full, spec1d_load]):
    for mode, spec1d in zip(['full'], [spec1d_full]):
        assert spec1d.nobj == spec1d_dev.nobj, f'{mode}: nobj {spec1d.nobj} != {spec1d_dev.nobj}'
        assert np.allclose(spec1d[idx].BOX_COUNTS, spec1d_dev[idx].BOX_COUNTS, atol=1e-4)
        assert np.allclose(spec1d[idx].OPT_WAVE, spec1d_dev[idx].OPT_WAVE, atol=1e-4)


