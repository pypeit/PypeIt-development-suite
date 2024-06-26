"""
Module to run tests on FlatField class
Requires files in Development suite and an Environmental variable
"""
from pathlib import Path
import os
import shutil

from IPython import embed

import pytest

from pypeit import calibrations
from pypeit.spectrographs.util import load_spectrograph
from pypeit import wavecalib
from pypeit import pypeitsetup

from pypeit.tests.tstutils import dummy_fitstbl, data_output_path

@pytest.fixture
def fitstbl():
    # Check for files
    root = Path(os.getenv('PYPEIT_DEV'), 'RAW_DATA/shane_kast_blue/600_4310_d55').resolve()
    files = [ 
        root / 'b1.fits.gz',    # arc
        root / 'b11.fits.gz',   # trace
        root / 'b21.fits.gz',   # bias
        root / 'b24.fits.gz',   # standard
        root / 'b27.fits.gz'    # science
    ]

    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue')
    setupc.build_fitstbl(files)
    setupc.fitstbl.finalize_usr_build(None, 'A')
    return setupc.fitstbl


@pytest.fixture
def multi_caliBrate(fitstbl):
    # Grab a science file for configuration specific parameters
    indx = fitstbl.find_frames('science', index=True)[0]
    sci_file = fitstbl.frame_paths(indx)
    # Par
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.config_specific_par(sci_file)
    par.reset_all_processimages_par(use_biasimage=False)
    #
    calib_par = par['calibrations']
    calib_par['bpm_usebias'] = False
    calib_par['slitedges']['sync_predict'] = 'nearest'

    multi_caliBrate = calibrations.MultiSlitCalibrations(fitstbl, calib_par, spectrograph,
                                                         data_output_path('Calibrations'))
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

    # Remove any pre-existing directory
    calib_dir = Path(multi_caliBrate.calib_dir).resolve()
    if calib_dir.exists():
        shutil.rmtree(calib_dir)

    # Setup
    multi_caliBrate.shape = (2048,350)
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

    flatImages = multi_caliBrate.get_flats()
    assert flatImages.pixelflat_norm.shape == (2048,350)
    assert flatImages.fit2illumflat(slits).shape == (2048,350)

    # Wave image
    slitmask = slits.slit_img()
    tilts = waveTilts.fit2tiltimg(slitmask)

    #
    mswave = wv_calib.build_waveimg(tilts, slits)
    assert mswave.shape == (2048,350)

    # Clean-up
    if calib_dir.exists():
        shutil.rmtree(calib_dir)


def test_reuse(multi_caliBrate, fitstbl):
    """
    Test that Calibrations appropriately reuses existing calibrations frames.
    """
    # In case of previous data or failures
    calib_dir = Path(multi_caliBrate.calib_dir).resolve()
    if calib_dir.exists():
        shutil.rmtree(calib_dir)

    # Perform the calibrations and check that the data are correctly
    # stored in memory
    multi_caliBrate.shape = (2048,350)
    multi_caliBrate.get_bpm()

#   Make them all
    multi_caliBrate.get_arc()
    multi_caliBrate.get_tiltimg()
    multi_caliBrate.get_slits()
    multi_caliBrate.get_wv_calib()
    multi_caliBrate.get_tilts()
    multi_caliBrate.get_flats()

    # Reset
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    multi_caliBrate_reuse = calibrations.MultiSlitCalibrations(fitstbl, par['calibrations'],
                                                               spectrograph, str(calib_dir))
    multi_caliBrate_reuse.reuse_calibs = True
    reset_calib(multi_caliBrate_reuse)

    # Read the calibrations
    #   - These don't source a calibration file
    multi_caliBrate_reuse.shape = (2048,350)
    multi_caliBrate_reuse.get_bpm()
    multi_caliBrate_reuse.get_arc()
    assert multi_caliBrate_reuse.msarc is not None, 'Should find cached data.'
    multi_caliBrate_reuse.get_tiltimg()
    assert multi_caliBrate_reuse.mstilt is not None, 'Should find cached data.'
    multi_caliBrate_reuse.get_slits()
    multi_caliBrate_reuse.get_wv_calib()
    multi_caliBrate_reuse.get_tilts()
    multi_caliBrate_reuse.get_flats()

    # Clean-up
    if calib_dir.exists():
        shutil.rmtree(calib_dir)

