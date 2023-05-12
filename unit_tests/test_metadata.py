from pathlib import Path
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
from pypeit.calibrations import Calibrations


def get_kast_blue_files():
    # Check for files
    file_root = Path(os.getenv('PYPEIT_DEV')) / 'RAW_DATA' / 'shane_kast_blue' / '600_4310_d55'
    assert file_root.exists(), 'Bad raw data directory'
    files = sorted(file_root.glob('b*'))
    assert len(files) > 0
    return files


def test_rm_rows():
    # Check for files
    files = get_kast_blue_files()
    # Init
    ps = PypeItSetup(files, spectrograph_name='shane_kast_blue')
    # Run
    ps.run()
    assert len(ps.fitstbl) == 26, 'Missing files'

    # Find the biases
    is_bias = ps.fitstbl.find_frames('bias', index=True)
    assert len(is_bias) == 10, 'Should find 10 bias frames'

    # Remove them
    ps.fitstbl.remove_rows(is_bias)
    assert len(ps.fitstbl) == 16, 'Should have removed 10 frames'
    assert not any(ps.fitstbl.find_frames('bias')), 'Should not find any bias frames'

    # TODO: Leaving this inconsistency for now.  BEWARE!!
    assert len(ps.file_list) == 26, \
            'Direct manipulation of fitstbl means PypeItSetup is not self-consistent'

    # Rebuild
    ps = PypeItSetup(files, spectrograph_name='shane_kast_blue')
    ps.run()

    # Remove the standard
    is_std = ps.fitstbl.find_frames('standard', index=True)
    assert is_std.size == 1, 'Should find one standard frame'
    is_sci = ps.fitstbl.find_frames('science', index=True)
    assert is_sci.size == 2, 'Should find two science frames'
    assert np.array_equal(ps.fitstbl['comb_id'][is_std].data, [1]), 'Combination groups changed'
    assert np.array_equal(ps.fitstbl['comb_id'][is_sci].data, [2,3]), 'Combination groups changed'

    # Remove the standard, but don't regroup
    ps.fitstbl.remove_rows(is_std)
    assert len(ps.fitstbl) == 25, 'Should have removed one standard'
    is_sci = ps.fitstbl.find_frames('science', index=True)
    assert np.array_equal(ps.fitstbl['comb_id'][is_sci].data, [2,3]), 'Combination groups changed'

    # Rebuild
    ps = PypeItSetup(files, spectrograph_name='shane_kast_blue')
    ps.run()
    is_std = ps.fitstbl.find_frames('standard', index=True)
    # Remove the standard and regroup
    ps.fitstbl.remove_rows(is_std, regroup=True)
    is_sci = ps.fitstbl.find_frames('science', index=True)
    assert np.array_equal(ps.fitstbl['comb_id'][is_sci].data, [1,2]), \
            'Combination group numbers not updated'

    # Rebuild
    ps = PypeItSetup(files, spectrograph_name='shane_kast_blue')
    ps.run()

    # Use the PypeItSetup wrapper
    is_std = ps.fitstbl.find_frames('standard', index=True)
    ps.remove_table_rows(is_std, regroup=True)
    assert len(ps.fitstbl) == 25, 'Should have removed one standard'
    assert len(ps.fitstbl) == len(ps.file_list), 'Length of file list should match metadata table'
    # TODO: Need tests that test propagation of removed rows to frametype and usrdata


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
    ps.fitstbl.set_calibration_groups()

    calib_file = data_path('test.calib')
    caldir = Path(data_path('')).resolve() / ps.par['calibrations']['calib_dir']
    Calibrations.association_summary(calib_file, ps.fitstbl, ps.spectrograph, caldir,
                                     overwrite=True)
    with open(calib_file, 'r') as f:
        calib = yaml.load(f, Loader=yaml.FullLoader)

    assert np.array_equal(list(calib['A'].keys()), ['--', 0]), \
            'Calibrations dictionary read incorrectly.'

    os.remove(calib_file)


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

