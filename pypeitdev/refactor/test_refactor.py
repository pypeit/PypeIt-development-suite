# Tests for refactoring

import os, sys
import shutil
import numpy as np
from pathlib import Path

from pypeit.spectrographs.util import load_spectrograph
from pypeit.wavecalib import WaveCalib
from pypeit.slittrace import SlitTraceSet
from pypeit.flatfield import FlatImages
from pypeit.scripts.reduce_by_step import ReducebyStep
from pypeit import specobjs
from pypeit import spec2dobj

import pytest
from IPython import embed

sys.path.append(os.path.join(
    os.path.abspath(
        os.environ["PYPEIT_DEV"]),"test_scripts"))

# Shane Kastblue 600/4310 D55

def run_the_steps():
    # Parser
    #args = ReducebyStep.parse_args(['shane_kast_blue_A.pypeit', 'b28.fits.gz', 'process', '--det', '1'])
    pass
    return
    '''
    # Need to do the Standard first to get the trace (not sure it is even used)
    pypeit_reduce_by_step shane_kast_blue_A.pypeit b24.fits.gz process --det 1
    pypeit_reduce_by_step shane_kast_blue_A.pypeit b24.fits.gz findobj --det 1
    pypeit_reduce_by_step shane_kast_blue_A.pypeit b24.fits.gz extract --det 1

    pypeit_reduce_by_step shane_kast_blue_A.pypeit b27.fits.gz process --det 1
    #pypeit_view_fits shane_kast_blue Intermediate/sciImg_b28-J1217p3905_KASTb_20150520T051801.470_DET01.fits --inter

    pypeit_reduce_by_step shane_kast_blue_A.pypeit b27.fits.gz findobj --det 1
    #pypeit_view_fits shane_kast_blue Intermediate/initSky_b28-J1217p3905_KASTb_20150520T051801.470_DET01.fits --inter

    pypeit_reduce_by_step shane_kast_blue_A.pypeit b27.fits.gz extract --det 1
    #pypeit_show_2dspec Science/spec2d_b28-J1217p3905_KASTb_20150520T051801.470.fits
    '''

def test_shane_kastb_calibs():
    setup = '600_4310_d55'
    _redux_out = os.path.join(os.environ['PYPEIT_DEV'], 'REDUX_OUT')

    # Slits
    root = 'Slits_A_0_DET01.fits.gz'
    slits_full_refactor_file = os.path.join(_redux_out, 'shane_kast_blue', setup,
                                                'REFACTOR_FULL','Calibrations', root)
    slits_dev_file = slits_full_refactor_file.replace('REFACTOR_FULL','DEVELOP')

    # Load
    slits_full = SlitTraceSet.from_file(slits_full_refactor_file)
    slits_dev = SlitTraceSet.from_file(slits_dev_file)
    assert np.allclose(slits_full.left_tweak, slits_dev.left_tweak, rtol=1e-8)

    # Flats
    root = 'Flat_A_0_DET01.fits'
    flats_full_refactor_file = os.path.join(_redux_out, 'shane_kast_blue', setup,
                                                'REFACTOR_FULL','Calibrations', root)
    flats_dev_file = flats_full_refactor_file.replace('REFACTOR_FULL','DEVELOP')

    # Load
    flats_full = FlatImages.from_file(flats_full_refactor_file)
    flats_dev = FlatImages.from_file(flats_dev_file)
    assert np.allclose(flats_full.pixelflat_norm, 
                    flats_dev.pixelflat_norm, rtol=1e-8)

def test_shane_kastb_sciimg():
    _redux_out = os.path.join(os.environ['PYPEIT_DEV'], 'REDUX_OUT')
    setup = '600_4310_d55'
    root = 'spec2d_b27-J1217p3905_KASTb_20150520T045733.560.fits'
    spec2d_full_refactor_file = os.path.join(_redux_out, 
                                                'shane_kast_blue', 
                                                setup,
                                                'REFACTOR_FULL','Science', root)

    # Load
    spec2d_full = spec2dobj.Spec2DObj.from_file(spec2d_full_refactor_file, 'DET01',
                                                    chk_version=False)
    spec2d_dev_file = spec2d_full_refactor_file.replace('REFACTOR_FULL','DEVELOP')
    spec2d_dev = spec2dobj.Spec2DObj.from_file(spec2d_dev_file, 'DET01',
                                                        chk_version=False)
    # Compare
    #pytest.set_trace()
    assert np.allclose(spec2d_full.sciimg, spec2d_dev.sciimg, rtol=1e-7)


'''
# FUSSING WITH TRACES
#sobjs_obj.write_to_fits({}, 'initial_sobjs.fits')
#sobjs.write_to_fits({}, 'extract_sobjs.fits')

_redux_out = os.path.join(os.environ['PYPEIT_DEV'], 'REDUX_OUT')

# Load spec1d's
setup = '600_4310_d55'
root_initial = 'initial_sobjs.fits'
root_ex = 'extract_sobjs.fits'

init_spec1d_full_refactor_file = os.path.join(_redux_out, 
                                            'shane_kast_blue', 
                                            setup,
                                            'REFACTOR_FULL', root_initial)
ex_spec1d_full_refactor_file = os.path.join(_redux_out, 
                                            'shane_kast_blue', 
                                            setup,
                                            'REFACTOR_FULL', root_ex)

init_spec1d_full = specobjs.SpecObjs.from_fitsfile(init_spec1d_full_refactor_file,
                                                chk_version=False)
init_spec1d_dev_file = init_spec1d_full_refactor_file.replace('REFACTOR_FULL','DEVELOP')
init_spec1d_dev = specobjs.SpecObjs.from_fitsfile(init_spec1d_dev_file,
                                                chk_version=False)
assert np.allclose(init_spec1d_full.TRACE_SPAT, init_spec1d_dev.TRACE_SPAT, rtol=1e-8)

ex_spec1d_full = specobjs.SpecObjs.from_fitsfile(ex_spec1d_full_refactor_file,
                                                chk_version=False)
ex_spec1d_dev_file = ex_spec1d_full_refactor_file.replace('REFACTOR_FULL','DEVELOP')
ex_spec1d_dev = specobjs.SpecObjs.from_fitsfile(ex_spec1d_dev_file,
                                                chk_version=False)

embed(header='End of loading spec1d')
assert np.allclose(ex_spec1d_full.TRACE_SPAT, ex_spec1d_dev.TRACE_SPAT)#, atol=1e-4)
'''

# FWHM (DEVELOP)
#In [23]: sobjs.FWHM
#Out[23]: array([71.78102448,  5.16074157, 13.061549  ])

#In [1]: sobjs.FWHM
#Out[1]: array([72.08452694,  5.15945478, 23.80638717])

# Simple RERUN from an embed

#In [1]: sobjs.FWHM
#Out[1]: array([71.58973139,  5.1598249 , 24.41008157])

#In [3]: sobjs.FWHM
#Out[3]: array([71.68468681,  5.160167  , 22.69068257])

#In [5]: sobjs.FWHM
#Out[5]: array([71.64438716,  5.15979802, 12.22337752])

# Calls to local_skysub_extract
#In [3]: self.sobjs[thisobj].FWHM
#Out[3]: array([71.8780582 ,  5.15966312, 15.98746125])

#In [8]: self.sobjs[thisobj].FWHM
#Out[8]: array([71.68468681,  5.160167  , 22.69068257])

# From within the subroutine
#In [6]: sobjs.FWHM
#Out[6]: array([71.67431597,  5.16014457, 22.56220688])

#In [8]: sobjs.FWHM
#Out[8]: array([71.74012149,  5.16000288, 22.66154206])




#_redux_out = os.path.join(os.environ['PYPEIT_DEV'], 'REDUX_OUT')
#setup = '600_4310_d55'
#root_ex = 'spec1d_b24-Feige66_KASTb_20150520T041246.960.fits'
#ex1_spec1d_full_refactor_file = os.path.join(_redux_out, 
#                                            'shane_kast_blue', 
#                                            setup,
#                                            'DEVELOP', 'Sci_1', root_ex)
#ex2_spec1d_full_refactor_file = ex1_spec1d_full_refactor_file.replace('Sci_1','Sci_2')
#ex3_spec1d_full_refactor_file = ex1_spec1d_full_refactor_file.replace('Sci_1','Sci_3')
#
## Load em
#ex1_spec1d_full = specobjs.SpecObjs.from_fitsfile(ex1_spec1d_full_refactor_file,
#                                                chk_version=False)
#ex2_spec1d_full = specobjs.SpecObjs.from_fitsfile(ex2_spec1d_full_refactor_file,
#                                                chk_version=False)
#ex3_spec1d_full = specobjs.SpecObjs.from_fitsfile(ex3_spec1d_full_refactor_file,
#                                                chk_version=False)
#
#embed(header='End of loading spec1d')
## Compare FWHM
#assert np.allclose(ex1_spec1d_full.FWHM, ex2_spec1d_full.FWHM, rtol=1e-8)
#assert np.allclose(ex1_spec1d_full.FWHM, ex3_spec1d_full.FWHM, rtol=1e-8)


def test_shane_kastb_spec1d():
    # Initial sky, traces look identical

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

    spec1d_steps_file = spec1d_full_refactor_file.replace('REFACTOR_FULL','STEP_BY_STEP')
    spec1d_steps = specobjs.SpecObjs.from_fitsfile(spec1d_steps_file,
                                                    chk_version=False)


    #spec1d_load_file = spec1d_full_refactor_file.replace('REFACTOR_FULL', 'shane_kast_blue_A')
    #spec1d_load = specobjs.SpecObjs.from_fitsfile(spec1d_load_file,
    #                                                chk_version=False)


    # Full vs. Develop
    #for mode, spec1d in zip(['full', 'load'], [spec1d_full, spec1d_load]):

    #pytest.set_trace()
    #assert np.allclose(spec1d_full.TRACE_SPAT, spec1d_dev.TRACE_SPAT)#, atol=1e-4)
    for mode, spec1d in zip(['Full', 'Steps'], [spec1d_full, spec1d_steps]):
        atol = 1e-7 if mode == 'Full' else 1e-5
        assert spec1d.nobj == spec1d_dev.nobj
        assert np.allclose(spec1d.TRACE_SPAT, spec1d_dev.TRACE_SPAT, atol=atol), f'{mode}: TRACE_SPAT'
        assert np.allclose(spec1d.BOX_COUNTS, spec1d_dev.BOX_COUNTS, atol=atol), f'Failed {mode} BOX_COUNTS'
        assert np.allclose(spec1d.OPT_WAVE, spec1d_dev.OPT_WAVE, atol=atol), f'Failed {mode} OPT_WAVE'



def test_mosfire_ylong_spec1d():
    _redux_out = os.path.join(os.environ['PYPEIT_DEV'], 'REDUX_OUT')

    # Load spec1d's
    setup = 'Y_long'

    root = 'spec1d_m191120_0043-J2132-1434_OFF_MOSFIRE_20191120T061253.347.fits'
    idx = 0  # Only 1 source on science exposure

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
    #pytest.set_trace()
    for mode, spec1d in zip(['full'], [spec1d_full]):
        #atol = 1e-4 if mode == 'full' else 1e-3
        assert spec1d.nobj == spec1d_dev.nobj, f'{mode}: nobj {spec1d.nobj} != {spec1d_dev.nobj}'
        assert np.allclose(spec1d[idx].BOX_COUNTS, spec1d_dev[idx].BOX_COUNTS, atol=1e-4), f'Failed {mode} BOX_COUNTS'
        assert np.allclose(spec1d[idx].OPT_WAVE, spec1d_dev[idx].OPT_WAVE, atol=1e-6), f'Failed {mode} OPT_WAVE'


# DEIMOS
#    pypeit_reduce_by_step keck_deimos_600zd_m_6500.pypeit d1010_0056.fits.gz process --det 1,5
#    pypeit_reduce_by_step keck_deimos_600zd_m_6500.pypeit d1010_0056.fits.gz process --det 2,6
#    pypeit_reduce_by_step keck_deimos_600zd_m_6500.pypeit d1010_0056.fits.gz process --det 3,7
#    pypeit_reduce_by_step keck_deimos_600zd_m_6500.pypeit d1010_0056.fits.gz process --det 4,8