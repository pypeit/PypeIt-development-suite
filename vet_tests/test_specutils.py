
from pathlib import Path

from IPython import embed

import pytest

from pypeit.specutils import pypeit_loaders
from pypeit.specutils import Spectrum1D, SpectrumList

import os
redux_out = Path(os.getenv('PYPEIT_DEV')).resolve() / 'REDUX_OUT'

def test_identify_as_pypeit_file(redux_out):
    rdx = Path(redux_out).resolve()

    # Spec2D file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'Science' \
                    / 'spec2d_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify spec2d file as a pypeit file'

    # Spec1D file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'Science' \
                    / 'spec1d_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify spec1d file as a pypeit file'
    
    # sensfunc file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' \
                    / 'sens_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify sensfunc file as a pypeit file'

    # coadd2d Spec2D file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'Science_coadd' \
                    / 'spec2d_cN20170331S0216-cN20170331S0221-pisco.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify coadd spec2d file as a pypeit file'

    # coadd2d Spec1D file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'Science_coadd' \
                    / 'spec1d_cN20170331S0216-cN20170331S0221-pisco.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify coadd spec1d file as a pypeit file'

    # coadd1d file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'pisco_coadd.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify coadd file as a pypeit file'

    # telluric-corrected coadd1d file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'pisco_coadd_tellcorr.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify telluric-corrected coadd file as a pypeit file'

    # coadd1d telluric model file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'pisco_coadd_tellmodel.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify telluric-corrected coadd file as a pypeit file'


def test_identify_as_spec1d_file(redux_out):
    rdx = Path(redux_out).resolve()

    # Spec1D file
    test_file = rdx / 'shane_kast_blue' / '600_4310_d55' / 'shane_kast_blue_A' \
                    / 'Science' / 'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_spec1d(None, test_file), \
                'Did not identify spec1d file as a pypeit file'

    # Spec1D file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'Science' \
                    / 'spec1d_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_spec1d(None, test_file), \
                'Did not identify spec1d file as a pypeit file'

    # coadd2d Spec1D file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'Science_coadd' \
                    / 'spec1d_cN20170331S0216-cN20170331S0221-pisco.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_spec1d(None, test_file), \
                'Did not identify coadd spec1d file as a pypeit file'


def test_identify_as_onespec_file(redux_out):
    rdx = Path(redux_out).resolve()

    # coadd1d file
    test_file = rdx / 'shane_kast_blue' / '600_4310_d55' / 'shane_kast_blue_A' \
                    / 'J1217p3905_coadd.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_onespec(None, test_file), \
                'Did not identify coadd file as a pypeit file'

    # coadd1d file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'pisco_coadd.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_onespec(None, test_file), \
                'Did not identify coadd file as a pypeit file'

    # telluric-corrected coadd1d file
    test_file = rdx / 'gemini_gnirs' / '32_SB_SXD' / 'pisco_coadd_tellcorr.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_onespec(None, test_file), \
                'Did not identify telluric-corrected coadd file as a pypeit file'


def test_load_spec1d(redux_out):
    rdx = Path(redux_out).resolve()

    # Spec1D file
    test_file = rdx / 'shane_kast_blue' / '600_4310_d55' / 'shane_kast_blue_A' \
                    / 'Science' / 'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    spec = SpectrumList.read(str(test_file))
    #spec_list = SpectrumList.read(str(test_file))

def test_load_onespec(redux_out):
    rdx = Path(redux_out).resolve()

    # coadd1d file
    test_file = rdx / 'shane_kast_blue' / '600_4310_d55' / 'shane_kast_blue_A' \
                    / 'J1217p3905_coadd.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    spec1 = Spectrum1D.read(str(test_file))

    spec2 = Spectrum1D.read(str(test_file), grid=True)

    embed()
    exit()

test_load_onespec(redux_out)

