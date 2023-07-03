import os
import numpy as np

from pypeit.spectrographs.util import load_spectrograph
from pypeit.wavecalib import WaveCalib

import pytest

def test_shane_kast_red(redux_out):

    for setup, setupID, index, rms in zip(
        ['300_7500_Ne', '600_7500_d57', '1200_5000_d57'],
        ['B_0', 'A_0', 'A_0'],
        [1, 0, 0],
        [0.05, 0.05, 0.055],
        ):
        # Check that spatial flexure shift was set!
        file_path = os.path.join(redux_out,
                             'shane_kast_red',
                             setup,
                             'Calibrations',
                             f'WaveCalib_{setupID}_DET01.fits')
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        assert waveCalib.wv_fits[index].rms < rms, f'RMS of shane_kast_red {setup} is too high!'

def test_not_alfosc(redux_out):

    for setup, rms in zip(
        ['grism3', 'grism4_nobin', 'grism5', 'grism7', 'grism10', 'grism11', 'grism17', 'grism18', 'grism19', 'grism20'],
        [0.45, 0.19, 0.15, 0.11, 0.20, 0.17, 0.15, 0.13, 0.05, 0.05],
        ):
        setupID = 'A_0'
        index = 0
        # Check that spatial flexure shift was set!
        file_path = os.path.join(redux_out,
                             'not_alfosc',
                             setup,
                             'Calibrations',
                             f'WaveCalib_{setupID}_DET01.fits')
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        assert waveCalib.wv_fits[index].rms < rms, f'RMS of not_alfosc {setup} is too high!'

def test_deimos(redux_out):

    for setup, index, rms, mosaic in zip(
        ['1200B_LVM_5200', '600ZD_M_6500', '900ZD_LVM_5500'],
        [3, 1, 1],
        [0.15, 0.35, 0.1],
        ['MSC03', 'MSC03', 'MSC03']
        ):
        setupID = 'A_0'
        # Check that spatial flexure shift was set!
        file_path = os.path.join(redux_out,
                             'keck_deimos',
                             setup,
                             'Calibrations',
                             f'WaveCalib_{setupID}_{mosaic}.fits')
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        assert waveCalib.wv_fits[index].rms < rms, f'RMS of keck_deimos {setup} is too high!'

def test_mdm_modspec(redux_out):

    for setup, rms in zip(
        ['Echelle'],
        [0.05],
        ):
        setupID = 'A_0'
        index = 0
        # Check that spatial flexure shift was set!
        file_path = os.path.join(redux_out,
                             'mdm_modspec',
                             setup,
                             'Calibrations',
                             f'WaveCalib_{setupID}_DET01.fits')
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        assert waveCalib.wv_fits[index].rms < rms, f'RMS of mdm_modspec {setup} is too high!'