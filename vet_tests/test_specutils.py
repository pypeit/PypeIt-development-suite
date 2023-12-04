
from pathlib import Path

from IPython import embed

import numpy as np

import pytest

try:
    from pypeit.specutils import pypeit_loaders
except ModuleNotFoundError:
    pypeit_loaders = None
from pypeit.specutils import Spectrum1D, SpectrumList

from pypeit.pypmsgs import PypeItError

specutils_required = pytest.mark.skipif(Spectrum1D is None or SpectrumList is None 
                                            or pypeit_loaders is None,
                                        reason='specutils not installed')


@specutils_required
def test_identify_as_pypeit_file(redux_out):
    rdx = Path(redux_out).resolve()

    # Spec2D file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'Science' \
                    / 'spec2d_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify spec2d file as a pypeit file'

    # Spec1D file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'Science' \
                    / 'spec1d_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify spec1d file as a pypeit file'
    
    # sensfunc file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' \
                    / 'sens_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify sensfunc file as a pypeit file'

    # coadd2d Spec2D file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'Science_coadd' \
                    / 'spec2d_cN20170331S0216-cN20170331S0221-pisco.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify coadd spec2d file as a pypeit file'

    # coadd2d Spec1D file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'Science_coadd' \
                    / 'spec1d_cN20170331S0216-cN20170331S0221-pisco.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify coadd spec1d file as a pypeit file'

    # coadd1d file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'pisco_coadd.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify coadd file as a pypeit file'

    # telluric-corrected coadd1d file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'pisco_coadd_tellcorr.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify telluric-corrected coadd file as a pypeit file'

    # coadd1d telluric model file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'pisco_coadd_tellmodel.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders._identify_pypeit(test_file), \
                'Did not identify telluric-corrected coadd file as a pypeit file'


@specutils_required
def test_identify_as_spec1d_file(redux_out):
    rdx = Path(redux_out).resolve()

    # Spec1D file
    test_file = rdx / 'shane_kast_blue' / '600_4310_d55' / 'shane_kast_blue_A' \
                    / 'Science' / 'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_spec1d(None, test_file), \
                'Did not identify spec1d file as a pypeit file'

    # Spec1D file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'Science' \
                    / 'spec1d_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_spec1d(None, test_file), \
                'Did not identify spec1d file as a pypeit file'

    # coadd2d Spec1D file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'Science_coadd' \
                    / 'spec1d_cN20170331S0216-cN20170331S0221-pisco.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_spec1d(None, test_file), \
                'Did not identify coadd spec1d file as a pypeit file'


@specutils_required
def test_identify_as_onespec_file(redux_out):
    rdx = Path(redux_out).resolve()

    # coadd1d file
    test_file = rdx / 'shane_kast_blue' / '600_4310_d55' / 'shane_kast_blue_A' \
                    / 'J1217p3905_coadd.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_onespec(None, test_file), \
                'Did not identify coadd file as a pypeit file'

    # coadd1d file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'pisco_coadd.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_onespec(None, test_file), \
                'Did not identify coadd file as a pypeit file'

    # telluric-corrected coadd1d file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'pisco_coadd_tellcorr.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    assert pypeit_loaders.identify_pypeit_onespec(None, test_file), \
                'Did not identify telluric-corrected coadd file as a pypeit file'


# TODO: Break some of these out into separate tests?
@specutils_required
def test_load_spec1d(redux_out):
    rdx = Path(redux_out).resolve()

    # Spec1D file
    test_file = rdx / 'shane_kast_blue' / '600_4310_d55' / 'shane_kast_blue_A' \
                    / 'Science' / 'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'

    # Test failure using a Spectrum1D read
    with pytest.raises(PypeItError):
        spec = Spectrum1D.read(str(test_file))

    # Correctly use SpectrumList to read the file
    spec = SpectrumList.read(str(test_file))
    assert len(spec) == 1, 'Should only be one spectrum'
    assert spec[0].meta['extract'] == 'OPT', 'Optimal extraction should have been performed'

    # Whether or not the spectrum has been flux-calibrated depends on which
    # tests were executed.  If the default read is fluxed, try reading the unfluxed version
    if spec[0].meta['fluxed']:
        unfluxed_spec = SpectrumList.read(str(test_file), fluxed=False)
        assert not unfluxed_spec[0].meta['fluxed'], 'Attempt to read unfluxed data failed'
        assert not np.array_equal(spec[0].flux.data, unfluxed_spec[0].flux.data), \
                'Fluxed spectrum should be different from unfluxed spectrum'

    # Try reading the unfluxed, box-extracted spectrum
    unfluxed_box_spec = SpectrumList.read(str(test_file), extract='BOX', fluxed=False)
    assert unfluxed_box_spec[0].meta['extract'] == 'BOX', 'Should have read boxcar extraction'

    # Try reading a file with multiple spectra
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'Science' \
                    / 'spec1d_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    spec = SpectrumList.read(str(test_file))
    assert len(spec) == 6, 'Expected 1 spectrum per order'


@specutils_required
def test_load_onespec(redux_out):
    rdx = Path(redux_out).resolve()

    # coadd1d file
    test_file = rdx / 'shane_kast_blue' / '600_4310_d55' / 'shane_kast_blue_A' \
                    / 'J1217p3905_coadd.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'

    # Try reading it as a list.  BEWARE: This works without bespoke pypeit code,
    # so there must be something within specutils that enables this.
    spec = SpectrumList.read(str(test_file))
    assert len(spec) == 1, 'Should have only read one spectrum'

    # Try reading as a single Spectrum1D
    spec = Spectrum1D.read(str(test_file))
    assert not spec.meta['grid'], 'Default read should use the contribution-weighted wavelengths'

    # Try reading the grid wavelength vector
    grid_spec = Spectrum1D.read(str(test_file), grid=True)
    assert grid_spec.meta['grid'], 'Did not set grid flag correctly.'
    assert not np.array_equal(spec.spectral_axis.data, grid_spec.spectral_axis.data), \
            'Wavelength vector did not change between grid and contribution-weighted read'

    # coadd1d file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'pisco_coadd.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    spec = Spectrum1D.read(str(test_file))
    assert spec.meta['extract'] == 'OPT', 'Expected optimal extraction'

    # telluric-corrected coadd1d file
    test_file = rdx / 'gemini_gnirs_echelle' / '32_SB_SXD' / 'pisco_coadd_tellcorr.fits'
    assert test_file.exists(), 'Output file does not exist or the name changed'
    spec = Spectrum1D.read(str(test_file))
    assert spec.meta['extract'] == 'OPT', 'Expected optimal extraction'


