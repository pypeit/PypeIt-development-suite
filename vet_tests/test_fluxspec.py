"""
Module to run tests on SensFunc and FluxCalibrate classes
Requires files in Development suite (REDUX_OUT) and an Environmental variable
"""
import os

import pytest

from IPython import embed

import numpy as np

from pypeit import fluxcalibrate
from pypeit import sensfunc
from pypeit.tests.tstutils import data_output_path
from pypeit.spectrographs.util import load_spectrograph
from pypeit import specobjs


@pytest.fixture
def kast_blue_files(request):
    redux_out = request.config.getoption("--redux_out")
    std_file = os.path.join(redux_out,
                            'shane_kast_blue', 
                            '600_4310_d55', 
                            'shane_kast_blue_A', 'Science',
                            'spec1d_b24-Feige66_KASTb_20150520T041246.960.fits')
    sci_file = os.path.join(redux_out,
                            'shane_kast_blue', 
                            '600_4310_d55', 
                            'shane_kast_blue_A', 'Science',
                            'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits')
    return [std_file, sci_file]


def test_sensfunc(kast_blue_files, request):

    sens_file = data_output_path('sensfunc.fits')
    redux_out = request.config.getoption("--redux_out")
    kast_blue_out = os.path.join(redux_out, 'shane_kast_blue', '600_4310_d55', 'shane_kast_blue_A')

    # Test the meta_spec data of the SensFunc
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    std_file, sci_file = kast_blue_files
    # Instantiate
    sensFunc = sensfunc.UVISSensFunc(std_file, sens_file, par['sensfunc'])
    # Test the standard loaded
    assert sensFunc.meta_spec['BINNING'] == '1,1'
    assert sensFunc.meta_spec['TARGET'] == 'Feige 66'

    # Vet the sensitivity function generated by the devsuite.

    # Load
    sensFunc = sensfunc.SensFunc.from_file(os.path.join(kast_blue_out, 'sens_b24-Feige66_KASTb_20150520T041246.960.fits'))

    # Validate some sens table columns
    sens_colnames = ['SENS_ZEROPOINT', 'SENS_ZEROPOINT_FIT', 'SENS_ZEROPOINT_FIT_GPM', 'SENS_WAVE']
    for col in sens_colnames:
        assert col in sensFunc.sens.keys(), f'Missing column "{col}".'

    # Validate the SensFunc attributes
    assert os.path.basename(sensFunc.std_cal) == 'feige66_002.fits.gz'
    assert len(sensFunc.wave) > 0
    assert len(sensFunc.zeropoint) > 0
    assert np.all(np.isfinite(sensFunc.wave))
    assert np.all(np.isfinite(sensFunc.zeropoint))
    assert not np.any(sensFunc.wave < 0)

def test_flux(kast_blue_files):

    # Validate fluxing information
    sobjs = specobjs.SpecObjs.from_fitsfile(kast_blue_files[1])
    sobj = sobjs[0]
    assert 'OPT_FLAM' in sobj.keys(), f'OPT_FLAM missing from {kast_blue_files[1]}'
    assert sobj.OPT_FLAM is not None, f'OPT_FLAM from {kast_blue_files[1]} is None'
    assert len(sobj.OPT_FLAM) > 0, f'OPT_FLAM from {kast_blue_files[1]} is empty'
    assert np.all(np.isfinite(sobj.OPT_FLAM)), f'OPT_FLAM from {kast_blue_files[1]} has infinite or NaN values'
    assert 'WAVE_RMS' in sobj.keys(), f'WAVE_RMS missing from {kast_blue_files[1]}'
    assert sobj.WAVE_RMS is not None, f'WAVE_RMS from {kast_blue_files[1]} is None'
    assert sobj.WAVE_RMS < 0.08 , f'WAVE_RMS from {kast_blue_files[1]} above threshold'
