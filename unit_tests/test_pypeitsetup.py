"""
Module to run tests on PypeItSetup class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from astropy.table import Table


from pypeit import pypeitsetup
from pypeit.par import pypeitpar
from pypeit.metadata import PypeItMetaData
from pypeit.tests.tstutils import data_path


def get_files():
    # Check for files
    file_root = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA/shane_kast_blue/600_4310_d55/b')
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    return files


def get_lrisr_files():
    # Check for files
    file_root = os.path.join(os.getenv('PYPEIT_DEV'), 
                             'RAW_DATA/keck_lris_red/long_600_7500_d560/LR')
    files = glob.glob(file_root+'*')
    files.sort()
    assert len(files) > 0
    return files


def test_init():
    # Init
    files = get_files()
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue')
    assert setupc.nfiles == 0


def test_build_fitstbl():
    # Check for files
    file_root = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA/shane_kast_blue/600_4310_d55/b')
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    # Init
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue')
    #
    fitstbl = setupc.build_fitstbl(files)
    assert isinstance(fitstbl, Table)
    assert setupc.nfiles == 26


def test_image_type():
    # Check for files
    files = get_files()
    # Init
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue')
    fitstbl = setupc.build_fitstbl(files)
    # Type
    setupc.get_frame_types(flag_unknown=True)
    assert 'framebit' in setupc.fitstbl.keys()
    assert np.sum(setupc.fitstbl.find_frames('arc')) == 1
    assert np.sum(setupc.fitstbl.find_frames('None')) == 0

    assert np.sum(setupc.fitstbl.find_frames('pixelflat')
                    & setupc.fitstbl.find_frames('trace')) == 12


def test_type():
    # Check for files
    files = get_files()
    # Init
    cfg_lines = ['[rdx]',
                 'spectrograph = shane_kast_blue',
                 '[calibrations]',
                 '[[standardframe]]',
                 'exprng = None,60',
                 '[scienceframe]',
                 'exprng = 60,None']
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue',
                                     cfg_lines=cfg_lines)
    setupc.build_fitstbl(files)
    setupc.get_frame_types(flag_unknown=True)
    assert np.sum(setupc.fitstbl.find_frames('science')) == 2


def test_run():
    # Check for files
    files = get_files()
    # Init
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue')
    # Run
    par, spectrograph, fitstbl = setupc.run()
    # Test
    assert isinstance(par, pypeitpar.PypeItPar)
    assert isinstance(fitstbl, PypeItMetaData)


def test_run_setup():
    files = get_files()
    # Init
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue')
    # Run
    par, spectrograph, fitstbl = setupc.run(setup_only=True)


def test_run_on_bad_headers():
    files = get_lrisr_files()
    # Init
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_lris_red')
    # Run
    par, spectrograph, fitstbl = setupc.run(setup_only=True)
    # Test
    idx = np.where(setupc.fitstbl['filename'] == 'LR.20160216.05709.fits.gz')[0]
    assert setupc.fitstbl['ra'][idx][0] is None
    assert len(setupc.fitstbl) == 23


