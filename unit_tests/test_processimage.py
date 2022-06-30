"""
Module to run tests on ProcessImages class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob

from IPython import embed

import numpy as np

from pypeit.images import rawimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core import procimg
from pypeit.par.pypeitpar import ProcessImagesPar
from pypeit import utils


@pytest.fixture
def deimos_flat_files():
    # Longslit in dets 3,7
    deimos_flat_files = [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos',
                                      '830G_L_8400', ifile) 
                            for ifile in ['d0914_0014.fits.gz', 'd0914_0015.fits.gz']]
    assert len(deimos_flat_files) == 2
    return deimos_flat_files

@pytest.fixture
def kast_blue_bias_files():
    return glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'shane_kast_blue',
                                  '600_4310_d55', 'b1?.fits*'))

@pytest.fixture
def kast_blue_arc_file():
    return glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'shane_kast_blue',
                                  '600_4310_d55', 'b1.fits*'))

def test_instantiate(deimos_flat_files, kast_blue_bias_files):
    one_file = deimos_flat_files[0]
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    det = 3
    rawImage = rawimage.RawImage(one_file, spectograph, det)
    # Test
    assert rawImage.datasec_img.shape == (1, 4096, 2128), 'Wrong shape'

    # Kast blue
    det2 = 1
    one_file = kast_blue_bias_files[0]
    spectograph2 = load_spectrograph('shane_kast_blue')
    rawImage2 = rawimage.RawImage(one_file, spectograph2, det2)
    assert rawImage2.image.shape == (1, 350, 2112), 'Wrong shape'


def test_overscan_subtract(deimos_flat_files):
    one_file = deimos_flat_files[0]
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    det = 3
    rawImage = rawimage.RawImage(one_file, spectograph, det)
    rawImage.par = spectograph.default_pypeit_par()['scienceframe']['process']
    # Bias subtract
    pre_sub = rawImage.image.copy()
    rawImage.subtract_overscan()
    oscan = np.median(pre_sub-rawImage.image)
    assert np.isclose(oscan, 1001.2, rtol=0.01)
    # Trim
    rawImage.trim()
    # Test
    assert rawImage.steps['subtract_overscan']
    assert rawImage.steps['trim']
    assert rawImage.image.shape == (1,4096,2048)


def test_continuum_subtraction(kast_blue_arc_file):
    one_file = kast_blue_arc_file[0]
    spectograph = load_spectrograph('shane_kast_blue')
    # Kast
    det = 1
    rawImage = rawimage.RawImage(one_file, spectograph, det)
    defpar = spectograph.default_pypeit_par()['calibrations']['arcframe']['process']
    defpar['subtract_continuum'] = True
    rawImage.par = defpar
    # Subtract continuum
    rawImage.subtract_continuum(force=True)
    # Test
    assert rawImage.steps['subtract_continuum']

def test_lacosmic():
    spec = load_spectrograph('keck_deimos')
    file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos', '1200G_M_5500',
                        'd0315_45929.fits')
    par = ProcessImagesPar(use_biasimage=False, use_pixelflat=False, use_illumflat=False)
    img = rawimage.RawImage(file, spec, 1)
    pimg = img.process(par)
    test_img = pimg.image[:500,:500]
    test_var = utils.inverse(pimg.ivar[:500,:500])
    crmask = procimg.lacosmic(test_img, varframe=test_var, maxiter=1)
    assert np.sum(crmask) == 1240, 'L.A.Cosmic changed'

    _crmask = procimg.lacosmic(test_img, varframe=test_var, maxiter=2)
    assert np.sum(_crmask) > np.sum(crmask), '2nd iteration should find more cosmics'

    _crmask = procimg.lacosmic(test_img, saturation=6000., varframe=test_var, maxiter=1)
    assert np.sum(_crmask) < np.sum(crmask), 'Should have flagged some pixels as saturated'

    __crmask = procimg.lacosmic(test_img, saturation=np.full(test_img.shape, 6000.),
                                    varframe=test_var, maxiter=1)
    assert np.array_equal(__crmask, _crmask), 'Saturation array failed.'


