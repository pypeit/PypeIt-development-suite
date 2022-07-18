import os
import glob
import yaml

from IPython import embed

import pytest

import numpy as np

from pypeit.pypeitsetup import PypeItSetup
from pypeit.tests.tstutils import data_path
from pypeit.metadata import PypeItMetaData
from pypeit.spectrographs.util import load_spectrograph
from pypeit import inputfiles


def test_lris_red_multi_400():
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_red',
                          'multi_400_8500_d560', '*.fits.gz'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.build_fitstbl()
    ps.get_frame_types(flag_unknown=True)
    cfgs = ps.fitstbl.unique_configurations()
    ps.fitstbl.set_configurations(cfgs)
    ps.fitstbl.set_calibration_groups() #global_frames=['bias', 'dark'])
    # Test
    assert np.all(ps.fitstbl['setup'] == 'A')


def test_lris_red_multi():
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_red',
                          'multi*', '*.fits*'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.build_fitstbl()
    ps.get_frame_types(flag_unknown=True)
    cfgs = ps.fitstbl.unique_configurations()
    ps.fitstbl.set_configurations(cfgs)
    ps.fitstbl.set_calibration_groups() #global_frames=['bias', 'dark'])


def test_lris_red_multi_calib():
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_red',
                          'multi_400_8500_d560', '*.fits.gz'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.build_fitstbl()
    ps.get_frame_types(flag_unknown=True)
    cfgs = ps.fitstbl.unique_configurations()
    ps.fitstbl.set_configurations(cfgs)
    ps.fitstbl.set_calibration_groups() #global_frames=['bias', 'dark'])

    cfile = data_path('test.calib') 
    ps.fitstbl.write_calib(cfile)
    with open(cfile, 'r') as f:
        calib = yaml.load(f, Loader=yaml.FullLoader)

    assert np.array_equal(list(calib['A'].keys()), ['--', 1]), \
            'Calibrations dictionary read incorrectly.'

    os.remove(cfile)


def test_lris_red_multi_run():
    # Perform the setup
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_red',
                          'multi*', '*.fits*'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.run(setup_only=True)

    # Test
    #assert len(ps.setup_dict) == 2, 'Should find two setups'
    assert len(ps.fitstbl) >= 40, 'Should find 40+ files'
    arcs = ps.fitstbl['filename'][ps.fitstbl.find_frames('arc')]
    assert len(arcs) >= 2, 'Should find two or more arcs'
    assert 'r170320_2017.fits.gz' in arcs, \
            'Should have identified r170320_2017.fits.gz as an arc'
    assert 'r170816_0057.fits' in ps.fitstbl['filename'][ps.fitstbl.find_frames('science')], \
            'Should have identified r170816_0057.fits as a science frame'

    # Clean-up
    #os.remove('keck_lris_red.lst')
    #os.remove('keck_lris_red.setups')
    os.remove('keck_lris_red.sorted')


def test_lris_blue_pypeit_overwrite():
    f = os.path.join(os.environ['PYPEIT_DEV'],
                     'pypeit_files/keck_lris_blue_long_400_3400_d560.pypeit')
    assert os.path.isfile(f), 'Could not find pypeit file.'
        
    pypeitFile = inputfiles.PypeItFile.from_file(f)

    # Reset path
    istr = pypeitFile.file_paths[0].find('RAW_DATA')
    pypeitFile.file_paths = [os.path.join(os.environ['PYPEIT_DEV'], 
                                          pypeitFile.file_paths[0][istr:])]
    data_files = pypeitFile.filenames

    # Read the fits table with and without the user data
    spectrograph = load_spectrograph('keck_lris_blue')
    par = spectrograph.default_pypeit_par()
    fitstbl = PypeItMetaData(spectrograph, par, files=data_files)
    fitstbl_usr = PypeItMetaData(spectrograph, par, files=data_files, 
                                 usrdata=pypeitFile.data)

    assert fitstbl['target'][0] == 'unknown', 'Grating name changed in file header'
    assert fitstbl_usr['target'][0] == 'test', 'Grating name changed in pypeit file'
    assert fitstbl['target'][0] != fitstbl_usr['target'][0], \
            'Fits header value and input pypeit file value expected to be different.'

