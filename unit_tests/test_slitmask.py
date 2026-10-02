"""
Simple tests that load files and parse the slit-mask design data.
"""

import os
from pathlib import Path

from IPython import embed
import numpy as np

from pypeit.spectrographs.util import load_spectrograph

# TODO: Add additional tests for other spectrographs with slit-mask design data

def test_get_slitmask_mmt_binospec():
    ifile = (
        Path(os.getenv('PYPEIT_DEV')).absolute() / 'RAW_DATA' / 'mmt_binospec' / 'Multislit_G270'
        / 'BOSS1441_350.Science.6785.fits'
    )
    spec = load_spectrograph('mmt_binospec')
    slits = spec.get_slitmask(ifile, det=1)

    assert slits.nslits == 72, 'Incorrect number of slits'
    assert np.array_equal(slits.slitid, np.arange(slits.nslits)+1), 'slit IDs are incorrect'

    assert slits.corners.shape == (72, 4, 2), 'Corner shape is incorrect'

    slits = spec.get_slitmask(ifile, det=2)

    assert slits.nslits == 73, 'Incorrect number of slits'
    assert np.array_equal(slits.slitid, np.arange(slits.nslits)+1), 'slit IDs are incorrect'

    assert slits.corners.shape == (73, 4, 2), 'Corner shape is incorrect'

    # Test that the slit mask code works for longslit data
    ifile = (
        Path(os.getenv('PYPEIT_DEV')).absolute() / 'RAW_DATA' / 'mmt_binospec' / 'Longslit_G600'
        / 'sci_img_2018.0209.035515.fits'
    )
    slits = spec.get_slitmask(ifile, det=1)
    assert slits.nslits == 1, 'Incorrect number of slits'

    # Test that the slit mask code works for longslit data
    ifile = (
        Path(os.getenv('PYPEIT_DEV')).absolute() / 'RAW_DATA' / 'mmt_binospec' / 'Longslit_G1000'
        / 'sci_img_2022.0327.063501.fits'
    )
    slits = spec.get_slitmask(ifile, det=1)
    assert slits.nslits == 1, 'Incorrect number of slits'
