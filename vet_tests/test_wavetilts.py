"""
Module to run tests on WaveTilts and BuildWaveTilts classes
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import numpy as np


from pypeit import wavetilts
from pypeit import slittrace
from pypeit import edgetrace
from pypeit.images import buildimage 

from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph



kastb_dir = os.path.join('shane_kast_blue', '600_4310_d55','shane_kast_blue_A')

instant_dict = dict(coeffs=np.ones((6,4,1)),
                    nslit=1,
                    spat_order=np.array([3]),
                    spec_order=np.array([5]),
                    spat_id=np.array([150]),
                    func2d='legendre2d')


def test_instantiate_from_master(redux_out):
    master_file = os.path.join(redux_out, kastb_dir, 'Calibrations',
                               'Tilts_A_0_DET01.fits')
    slit_master_file = os.path.join(redux_out, kastb_dir, 'Calibrations',
                                    'Slits_A_0_DET01.fits.gz')
    slits = slittrace.SlitTraceSet.from_file(slit_master_file)
    waveTilts = wavetilts.WaveTilts.from_file(master_file)
    tilts = waveTilts.fit2tiltimg(slits.slit_img())
    assert isinstance(tilts, np.ndarray)


# Test rebuild tilts with a flexure offset
def test_flexure(redux_out):
    flexure = 1.
    master_file = os.path.join(redux_out, kastb_dir, 'Calibrations',
                               'Tilts_A_0_DET01.fits')
    waveTilts = wavetilts.WaveTilts.from_file(master_file)
    # Need slitmask
    slit_file = os.path.join(redux_out, kastb_dir, 'Calibrations',
                             'Slits_A_0_DET01.fits.gz')
    slits = slittrace.SlitTraceSet.from_file(slit_file)
    slitmask = slits.slit_img(flexure=flexure)
    # Do it
    new_tilts = waveTilts.fit2tiltimg(slitmask, flexure=flexure)
    # Test?

def test_run(redux_out):
    # Masters
    spectrograph = load_spectrograph('shane_kast_blue')
    master_file = os.path.join(redux_out, kastb_dir, 'Calibrations',
                               'Tiltimg_A_0_DET01.fits')
    mstilt = buildimage.TiltImage.from_file(master_file)
    # Slits
    slit_file = os.path.join(redux_out, kastb_dir, 'Calibrations',
                             'Slits_A_0_DET01.fits.gz')
    # Instantiate
    #spectrograph.detector[0]['saturation'] = 60000.
    #spectrograph.detector[0]['nonlinear'] = 0.9
    par = pypeitpar.WaveTiltsPar()
    wavepar = pypeitpar.WavelengthSolutionPar()
    slits = slittrace.SlitTraceSet.from_file(slit_file)
    buildwaveTilts = wavetilts.BuildWaveTilts(mstilt, slits, spectrograph, par, wavepar, det=1)
    # Run
    waveTilts = buildwaveTilts.run(doqa=False)
    assert isinstance(waveTilts.fit2tiltimg(slits.slit_img()), np.ndarray)


