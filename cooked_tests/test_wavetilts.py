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
from pypeit.tests.tstutils import data_path



kastb_dir = os.path.join(os.getenv('PYPEIT_DEV'), 
                            'REDUX_OUT',
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A')

@pytest.fixture
def master_dir():
    return os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'shane_kast_blue')

instant_dict = dict(coeffs=np.ones((6,4,1)),
                    nslit=1,
                    spat_order=np.array([3]),
                    spec_order=np.array([5]),
                    spat_id=np.array([150]),
                    func2d='legendre2d')


def test_instantiate_from_master(master_dir):
    master_file = os.path.join(kastb_dir, 'Masters',
                               'MasterTilts_A_1_DET01.fits')
    slit_master_file = os.path.join(kastb_dir, 'Masters',
                               'MasterSlits_A_1_DET01.fits.gz')
    slits = slittrace.SlitTraceSet.from_file(slit_master_file)
    waveTilts = wavetilts.WaveTilts.from_file(master_file)
    tilts = waveTilts.fit2tiltimg(slits.slit_img())
    assert isinstance(tilts, np.ndarray)


# Test rebuild tilts with a flexure offset
def test_flexure(master_dir):
    flexure = 1.
    master_file = os.path.join(kastb_dir, 'Masters',
                               'MasterTilts_A_1_DET01.fits')
    waveTilts = wavetilts.WaveTilts.from_file(master_file)
    # Need slitmask
    slit_file = os.path.join(kastb_dir, 'Masters',
                             'MasterSlits_A_1_DET01.fits.gz')
    slits = slittrace.SlitTraceSet.from_file(slit_file)
    slitmask = slits.slit_img(flexure=flexure)
    # Do it
    new_tilts = waveTilts.fit2tiltimg(slitmask, flexure=flexure)
    # Test?

def test_run(master_dir):
    # Masters
    spectrograph = load_spectrograph('shane_kast_blue')
    master_file = os.path.join(kastb_dir, 'Masters',
                               'MasterTiltimg_A_1_DET01.fits')
    mstilt = buildimage.TiltImage.from_file(master_file)
    trace_file = os.path.join(kastb_dir, 'Masters',
                               'MasterEdges_A_1_DET01.fits.gz')
    edges = edgetrace.EdgeTraceSet.from_file(trace_file)
    # Instantiate
    #spectrograph.detector[0]['saturation'] = 60000.
    #spectrograph.detector[0]['nonlinear'] = 0.9
    par = pypeitpar.WaveTiltsPar()
    wavepar = pypeitpar.WavelengthSolutionPar()
    slits = edges.get_slits()
    buildwaveTilts = wavetilts.BuildWaveTilts(mstilt, slits, spectrograph, par, wavepar, det=1)
    # Run
    waveTilts = buildwaveTilts.run(doqa=False)
    assert isinstance(waveTilts.fit2tiltimg(slits.slit_img()), np.ndarray)


