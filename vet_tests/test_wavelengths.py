import os, sys
import shutil
import numpy as np
from pathlib import Path

from pypeit.spectrographs.util import load_spectrograph
from pypeit.wavecalib import WaveCalib
from pypeit.slittrace import SlitTraceSet
from pypeit.scripts.run_pypeit import RunPypeIt

import pytest
from IPython import embed

sys.path.append(os.path.join(
    os.path.abspath(
        os.environ["PYPEIT_DEV"]),"test_scripts"))
import pypeit_tests


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
        [0.45, 0.19, 0.15, 0.11, 0.20, 0.17, 0.15, 0.14, 0.05, 0.05],
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

def test_redoslits_kastr(redux_out):
    """ Test the redo_slits option using shane_kast_red

    Args:
        redux_out (str): path to REDUX_OUT
    """

    setup = '600_5000_d46'

    rdx_dir = os.path.join(redux_out,
                             'shane_kast_red',
                             setup)
    # Artificially make the slit bad
    slit_file = os.path.join(rdx_dir,
                             'Calibrations',
                             'Slits_A_0_DET01.fits.gz')
    # Copy the original
    orig_slit_file = slit_file.replace('.fits.gz', '_orig.fits.gz')
    shutil.copyfile(slit_file, orig_slit_file)

    # Modify
    slits = SlitTraceSet.from_file(slit_file)
    slits.mask[0] = slits.bitmask.turn_on(slits.mask[0], 'BADWVCALIB')
    slits.to_file(slit_file, overwrite=True)

    # Copy the pypeit file
    root_redoslit_file = 'shane_kast_red_redoslit_600_5000_d46.pypeit'
    pyp_file = os.path.join(os.path.abspath(
        os.environ["PYPEIT_DEV"]), "vet_tests", "files", root_redoslit_file)

    new_pyp_file = os.path.join(redux_out,
                             'shane_kast_red',
                             setup, root_redoslit_file)
                             
    raw_data_path =   os.path.join(os.path.abspath(
            os.environ["PYPEIT_DEV"]),"RAW_DATA",
            'shane_kast_red', setup)

    pypeit_tests.fix_pypeit_file_directory(
        pyp_file, None, raw_data_path,
        None, None, None, outfile=new_pyp_file)
        
    # Run 
    sv_cd = os.getcwd()
    os.chdir(rdx_dir)
    pargs = RunPypeIt.parse_args([str(new_pyp_file), '-o'])
    RunPypeIt.main(pargs)

    # Check
    slits2 = SlitTraceSet.from_file(slit_file)
    assert slits2.mask[0] == 0, 'Slit was not fixed!'

    # Copy back
    shutil.copyfile(orig_slit_file, slit_file)

    os.chdir(sv_cd)


def test_keck_lris_blue(redux_out):

    _redux_out = Path(redux_out).resolve()

    for setup, det, rms in zip(['multi_300_5000_d680', 'long_400_3400_d560', 'long_600_4000_d560'],
                               [1, 1, 1],
                               [0.32, 0.21, 0.08]):

        # get WaveCalib file
        this_calib = _redux_out / 'keck_lris_blue' / setup / 'Calibrations'
        assert this_calib.is_dir(), f'Calibration directory {this_calib} does not exist!'
        file_path = list(this_calib.glob(f'WaveCalib_*_DET0{det}.fits'))[0]
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        # get all the rms values
        rms_vals = np.array([ww.rms for ww in waveCalib.wv_fits if ww.rms is not None])
        # check the wavelength solution rms
        assert np.all(rms_vals <= rms), f'wave RMS for setup {setup} is too high!'


def test_keck_lris_blue_orig(redux_out):

    _redux_out = Path(redux_out).resolve()

    for setup, det, rms in zip(['long_600_4000_d500', 'multi_1200_3400_d460'],
                               [1, 1],
                               [0.25, 0.2]):

        # get WaveCalib file
        this_calib = _redux_out / 'keck_lris_blue_orig' / setup / 'Calibrations'
        assert this_calib.is_dir(), f'Calibration directory {this_calib} does not exist!'
        file_path = list(this_calib.glob(f'WaveCalib_*_DET0{det}.fits'))[0]
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        # get all the rms values
        rms_vals = np.array([ww.rms for ww in waveCalib.wv_fits if ww.rms is not None])
        # check the wavelength solution rms
        assert np.all(rms_vals <= rms), f'wave RMS for setup {setup} is too high!'


def test_keck_lris_red(redux_out):

    _redux_out = Path(redux_out).resolve()

    for setup, det, rms in zip(['long_150_7500_d560', 'long_300_5000_d560', 'long_400_8500_longread',
                                'multi_600_5000_d560', 'long_600_7500_d560', 'long_600_10000_d680',
                                'mulit_831_8200_d560', 'multi_900_5500_d560', 'long_1200_7500_d560',
                                'multi_1200_9000_d680'],
                               [2, 2, 2, 1, 1, 1, 2, 1, 1, 2],
                               [0.5, 0.1, 0.1, 0.23, 0.05, 0.1, 0.15, 0.1, 0.03, 0.07]):

        # get WaveCalib file
        this_calib = _redux_out / 'keck_lris_red' / setup / 'Calibrations'
        assert this_calib.is_dir(), f'Calibration directory {this_calib} does not exist!'
        file_path = list(this_calib.glob(f'WaveCalib_*_DET0{det}.fits'))[0]
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        # get all the rms values
        rms_vals = np.array([ww.rms for ww in waveCalib.wv_fits if ww.rms is not None])
        # check the wavelength solution rms
        assert np.all(rms_vals <= rms), f'wave RMS for setup {setup} is too high!'


def test_keck_lris_red_orig(redux_out):

    _redux_out = Path(redux_out).resolve()

    for setup, rms in zip(['long_150_7500_d500', 'long_300_5000', 'long_400_8500_d560', 'multi_600_5000_d500',
                           'long_600_7500_d680', 'long_600_10000_d460', 'long_831_8200_d460',
                           'long_900_5500_d560', 'long_1200_7500_d560'],
                          [0.15, 0.05, 0.08, 0.13, 0.18, 0.05, 0.05, 0.05, 0.08]):

        # get WaveCalib file
        this_calib = _redux_out / 'keck_lris_red_orig' / setup / 'Calibrations'
        assert this_calib.is_dir(), f'Calibration directory {this_calib} does not exist!'
        file_path = list(this_calib.glob(f'WaveCalib_*.fits'))[0]
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        # get all the rms values
        rms_vals = np.array([ww.rms for ww in waveCalib.wv_fits if ww.rms is not None])
        # check the wavelength solution rms
        assert np.all(rms_vals <= rms), f'wave RMS for setup {setup} is too high!'


def test_keck_lris_red_mark4(redux_out):

    _redux_out = Path(redux_out).resolve()

    for setup, rms in zip(['long_400_8500_d560', 'long_600_10000_d680'],
                          [0.07, 0.08]):

        # get WaveCalib file
        this_calib = _redux_out / 'keck_lris_red_mark4' / setup / 'Calibrations'
        assert this_calib.is_dir(), f'Calibration directory {this_calib} does not exist!'
        file_path = list(this_calib.glob(f'WaveCalib_*.fits'))[0]
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        # get all the rms values
        rms_vals = np.array([ww.rms for ww in waveCalib.wv_fits if ww.rms is not None])
        # check the wavelength solution rms
        assert np.all(rms_vals <= rms), f'wave RMS for setup {setup} is too high!'


def test_keck_hires(redux_out):

    _redux_out = Path(redux_out).resolve()
    for setup, rms in zip(['J0100+2802_H204Hr_RED_C1_ECH_-0.82_XD_1.62_1x2',
                           'J0100+2802_H204Hr_RED_C1_ECH_0.75_XD_1.69_1x2',
                           'J0100+2802_H237Hr_RED_C1_ECH_-0.91_XD_1.46_1x2',
                           'J0100+2802_H237Hr_RED_C1_ECH_0.88_XD_1.46_1x2',
                           'J0100+2802_N255Hr_RED_C2_ECH_0.74_XD_1.39_1x3',
                           'J0306+1853_U074_RED_C2_ECH_-0.86_XD_1.31_1x3',
                           'J0306+1853_U074_RED_C2_ECH_0.72_XD_1.42_1x3',
                           'J1723+2243_W241_RED_C5_ECH_-0.15_XD_0.90_2x2',
                           'Q1009+2956_G10H_BLUE_C5_ECH_-0.00_XD_1.02_1x3'],
                          [0.37,
                           0.29,
                           0.37,
                           0.25,
                           0.45,
                           0.32,
                           0.30,
                           0.21,
                           0.35]):

        # get WaveCalib file
        this_calib = _redux_out / 'keck_hires' / setup / 'Calibrations'
        assert this_calib.is_dir(), f'Calibration directory {this_calib} does not exist!'
        file_path = list(this_calib.glob(f'WaveCalib_*.fits'))[0]
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        # get all the rms values
        rms_vals = np.array([ww.rms for ww in waveCalib.wv_fits if ww.rms is not None])
        # check the wavelength solution rms
        assert np.all(rms_vals <= rms), f'wave RMS for setup {setup} is too high!'


def test_gmos(redux_out):

    for setup, index, rms, mosaic in zip(
        ['GS_HAM_B480_550'],
        [1],
        [0.30],
        ['MSC01'],
        ):
        setupID = 'A_0'
        # Check that spatial flexure shift was set!
        file_path = os.path.join(redux_out,
                             'gemini_gmos',
                             setup,
                             'Calibrations',
                             f'WaveCalib_{setupID}_{mosaic}.fits')
        # Load
        waveCalib = WaveCalib.from_file(file_path)
        assert waveCalib.wv_fits[index].rms < rms, f'RMS of gemini_gmos {setup} is too high!'

