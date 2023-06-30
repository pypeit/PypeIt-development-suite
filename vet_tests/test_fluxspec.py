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
from pypeit.tests.tstutils import data_path
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


def test_gen_sensfunc(kast_blue_files):

    sens_file = data_path('sensfunc.fits')
    if os.path.isfile(sens_file):
        os.remove(sens_file)

    # Get it started
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    std_file, sci_file = kast_blue_files
    # Instantiate
    sensFunc = sensfunc.UVISSensFunc(std_file, sens_file, par['sensfunc'])
    # Test the standard loaded
    assert sensFunc.meta_spec['BINNING'] == '1,1'
    assert sensFunc.meta_spec['TARGET'] == 'Feige 66'

    # Generate the sensitivity function
    sensFunc.run()
    # Test
    assert os.path.basename(sensFunc.std_cal) == 'feige66_002.fits'
    # TODO: @jhennawi, please check this edit
    assert 'SENS_ZEROPOINT' in sensFunc.sens.keys(), 'Bad column names'
    # Write
    sensFunc.to_file(sens_file)

    os.remove(sens_file)


def test_from_sens_func(kast_blue_files):

    sens_file = data_path('sensfunc.fits')
    if os.path.isfile(sens_file):
        os.remove(sens_file)

    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    std_file, sci_file = kast_blue_files

    # Build the sensitivity function
    sensFunc = sensfunc.UVISSensFunc(std_file, sens_file, par['sensfunc'])
    sensFunc.run()
    sensFunc.to_file(sens_file)

    # Instantiate and run
    outfile = data_path(os.path.basename(sci_file))
    fluxCalibrate = fluxcalibrate.FluxCalibrate([sci_file], [sens_file], par=par['fluxcalib'],
                                              outfiles=[outfile])
    # Test
    sobjs = specobjs.SpecObjs.from_fitsfile(outfile)
    assert 'OPT_FLAM' in sobjs[0].keys()

    os.remove(sens_file)
    os.remove(outfile)


