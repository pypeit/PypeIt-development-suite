"""
Module to run tests on FlatField class
Requires files in Development suite and an Environmental variable
"""
from pathlib import Path
import os

import pytest
import glob
import shutil

import numpy as np

from pypeit import calibrations
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit import wavecalib
from IPython import embed

from pypeit.tests.tstutils import dummy_fitstbl, data_path


@pytest.fixture
def fitstbl():
    if os.getenv('PYPEIT_DEV') is None:
        fitstbl = dummy_fitstbl(directory=data_path(''))
        fitstbl['framebit'][0] = fitstbl.type_bitmask.turn_off(fitstbl['framebit'][0], flag='bias')
        fitstbl['filename'][1] = 'b1.fits.gz'
        fitstbl['filename'][5] = 'b27.fits.gz'
        return fitstbl

    fitstbl = dummy_fitstbl(directory=os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA',
                                                   'shane_kast_blue', '600_4310_d55'))
    # Set the Bias to known
    fitstbl['framebit'][0] = fitstbl.type_bitmask.turn_off(fitstbl['framebit'][0], flag='bias')
    fitstbl['filename'][1] = 'b1.fits.gz'
    for ii in range(2,5):
        fitstbl['filename'][ii] = 'b{0}.fits.gz'.format(ii)
    fitstbl['filename'][5] = 'b27.fits.gz'

    return fitstbl


@pytest.fixture
def multi_caliBrate(fitstbl):
    # Grab a science file for configuration specific parameters
    for idx, row in enumerate(fitstbl):
        if 'science' in row['frametype']:
            sci_file = os.path.join(row['directory'], row['filename'])
            break
    # Par
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.config_specific_par(sci_file)
    turn_off = dict(use_biasimage=False)
    par.reset_all_processimages_par(**turn_off)
    #
    calib_par = par['calibrations']
    calib_par['bpm_usebias'] = False
    #calib_par['biasframe']['useframe'] = 'none' # Only use overscan
    calib_par['slitedges']['sync_predict'] = 'nearest'

    multi_caliBrate = calibrations.MultiSlitCalibrations(fitstbl, calib_par, spectrograph,
                                                         data_path('Calibrations'))
    return reset_calib(multi_caliBrate)


def reset_calib(calib):
    # Find the first science row
    frame = calib.fitstbl.find_frames('science', index=True)[0]
    # Set
    det = 1
    calib.set_config(frame, det)
    return calib


###################################################
# TESTS BEGIN HERE


def test_it_all(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    #multi_caliBrate.get_pixlocn()
    multi_caliBrate.get_bpm()
    multi_caliBrate.get_arc()
    multi_caliBrate.get_tiltimg()
    slits = multi_caliBrate.get_slits()
    assert slits.PYP_SPEC == 'shane_kast_blue', 'Wrong spectrograph'
    assert (slits.nspec, slits.nspat) == multi_caliBrate.shape, 'Wrong image shape'
    assert slits.nslits == 1, 'Incorrect number of slits'
    assert slits.left_init.shape == (2048,1), 'Incorrect shape for left'
    assert slits.left_tweak is None, 'Tweaks should not exist'

    wv_calib = multi_caliBrate.get_wv_calib()
    assert isinstance(wv_calib, wavecalib.WaveCalib)
    assert 175 in wv_calib.spat_ids
    assert wv_calib.wv_fits[0]['rms'] < 0.2

    waveTilts = multi_caliBrate.get_tilts()
    assert waveTilts.nslit == 1

    multi_caliBrate.get_flats()
    flatImages = multi_caliBrate.get_flats()
    assert flatImages.pixelflat_norm.shape == (2048,350)
    assert flatImages.fit2illumflat(slits).shape == (2048,350)

    # Wave image
    slitmask = slits.slit_img()
    tilts = waveTilts.fit2tiltimg(slitmask)

    #
    mswave = wv_calib.build_waveimg(tilts, slits)
    assert mswave.shape == (2048,350)

def test_reuse(multi_caliBrate, fitstbl):
    """
    Test that Calibrations appropriately reuses existing calibrations frames.
    """
    # In case of previous data or failures
    calib_dir = Path(multi_caliBrate.calib_dir).resolve()
    if calib_dir.exists():
        shutil.rmtree(calib_dir)

    calib_dir.mkdir(parents=True)

    # Perform the calibrations and check that the data are correctly
    # stored in memory
    multi_caliBrate.shape = (2048,350)
    multi_caliBrate.get_bpm()

#   Make them all
    msarc = multi_caliBrate.get_arc()
    msarc = multi_caliBrate.get_tiltimg()
    multi_caliBrate.get_slits()
    multi_caliBrate.get_wv_calib()
    multi_caliBrate.get_tilts()
    multi_caliBrate.get_flats()

    # Reset
    #reset_calib(multi_caliBrate_reuse)
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    multi_caliBrate_reuse = calibrations.MultiSlitCalibrations(fitstbl, par['calibrations'],
                                                               spectrograph, str(calib_dir))
    multi_caliBrate_reuse.reuse = True
    reset_calib(multi_caliBrate_reuse)

    # Read the calibrations
    #   - These don't source a calibration file
    multi_caliBrate_reuse.shape = (2048,350)
    multi_caliBrate_reuse.get_bpm()
    _msarc = multi_caliBrate_reuse.get_arc()
    assert multi_caliBrate_reuse.msarc is not None, 'Should find cached data.'
    _msarc = multi_caliBrate_reuse.get_tiltimg()
    assert multi_caliBrate_reuse.mstilt is not None, 'Should find cached data.'
    multi_caliBrate_reuse.get_slits()
    multi_caliBrate_reuse.get_wv_calib()
    multi_caliBrate_reuse.get_tilts()
    multi_caliBrate_reuse.get_flats()

    # Clean-up
    shutil.rmtree(multi_caliBrate_reuse.calib_dir)

