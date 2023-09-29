
from pathlib import Path
import os
import glob
import shutil

from IPython import embed

import numpy as np

import pytest

from pypeit.pypeitsetup import PypeItSetup
from pypeit.inputfiles import PypeItFile

def test_deimos():
    # Raw DEIMOS directory
    raw_dir = os.path.join(os.getenv('PYPEIT_DEV'), 
                           'RAW_DATA', 'keck_deimos')

    # Get the list of setup directories
    setups = glob.glob(os.path.join(raw_dir, '*'))

    # Set the output path and *remove if* if it already exists
    output_path = os.path.join(os.getcwd(), 'output')
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)

    # Iterate through the setups
    for setup in setups:
 
        # Find the relevant pypeit file constructed by hand.
        by_hand_pypeit = os.path.join(os.getenv('PYPEIT_DEV'), 'pypeit_files',
                                      'keck_deimos_{0}.pypeit'.format(
                                        os.path.split(setup)[1].lower()))

        if not os.path.isfile(by_hand_pypeit):
            # It doesn't exist, so assume there is no by-hand pypeit
            # file to compare to
            continue

        # Run pypeit_setup
        ps = PypeItSetup.from_file_root(setup, 'keck_deimos')
        ps.run(setup_only=True)
        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        # Read the frame types from both the by-hand and automated
        # pypeit files
        byhand_pypeItFile = PypeItFile.from_file(by_hand_pypeit)
        auto_pypeItFile = PypeItFile.from_file(pypeit_files[0])

        # For each file in the by-hand list, check that the frame types
        # in the automatically generated pypeit file are identical
        for f in byhand_pypeItFile.frametypes.keys():
            type_list = np.sort(byhand_pypeItFile.frametypes[f].split(','))
            if 'science' in type_list or 'standard' in type_list:
                # Only ensuring that calibrations are correctly typed
                continue
            assert f in auto_pypeItFile.frametypes.keys(), \
                'Frame {0} not automatically parsed for setup {1}.'.format(f, setup)
            assert np.array_equal(type_list, np.sort(auto_pypeItFile.frametypes[f].split(','))), \
                'Frame types differ for file {0} in setup {1}\n'.format(f, setup) \
                 + '    By-hand types: {0}'.format(byhand_pypeItFile.frametypes[f]) \
                 + '    Automated types: {0}'.format(auto_pypeItFile.frametypes[f])

        # Clean up after every setup
        shutil.rmtree(output_path)


def test_mosfire():
    # Raw MOSFIRE directory
    raw_dir = os.path.join(os.getenv('PYPEIT_DEV'), 
                           'RAW_DATA', 'keck_mosfire')

    # Get the list of setup directories
    setups = [os.path.join(raw_dir, isetup)
                for isetup in ['Y_multi', 'mask1_K_with_continuum', 'mask1_J_with_continuum']]

    # Set the output path and *remove if* if it already exists
    output_path = os.path.join(os.getcwd(), 'output')
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)

    # Iterate through the setups
    for setup in setups:

        # Find the relevant pypeit file constructed by hand.
        by_hand_pypeit = os.path.join(os.getenv('PYPEIT_DEV'), 'pypeit_files',
                                      'keck_mosfire_{0}.pypeit'.format(
                                          os.path.split(setup)[1].lower()))

        if not os.path.isfile(by_hand_pypeit):
            # It doesn't exist, so assume there is no by-hand pypeit
            # file to compare to
            continue

        # Run pypeit_setup
        ps = PypeItSetup.from_file_root(setup, 'keck_mosfire')
        ps.run(setup_only=True)
        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        # Read the frame types from both the by-hand and automated
        # pypeit files
        byhand_pypeItFile = PypeItFile.from_file(by_hand_pypeit)
        auto_pypeItFile = PypeItFile.from_file(pypeit_files[0])

        # For each file in the by-hand list, check that the frame types
        # in the automatically generated pypeit file are identical
        for f in byhand_pypeItFile.frametypes.keys():
            type_list = np.sort(byhand_pypeItFile.frametypes[f].split(','))
            assert f in auto_pypeItFile.frametypes.keys(), \
                'Frame {0} not automatically parsed for setup {1}.'.format(f, setup)
            assert np.array_equal(type_list, np.sort(auto_pypeItFile.frametypes[f].split(','))), \
                'Frame types differ for file {0} in setup {1}\n'.format(f, setup) \
                + '    By-hand types: {0}'.format(byhand_pypeItFile.frametypes[f]) \
                + '    Automated types: {0}'.format(auto_pypeItFile.frametypes[f])

        # Clean up after every setup
        shutil.rmtree(output_path)


def test_nires():

    # Raw NIRES directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'

    raw_dir = dev_root / 'RAW_DATA' / 'keck_nires'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the list of setup directories
    setups = list(raw_dir.glob('*'))

    # Set the output path and *remove if* if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Iterate through the setups
    for setup in setups:
 
        # Find the relevant pypeit file constructed by hand.
        by_hand_pypeit = dev_root / 'pypeit_files' / f'keck_nires_{setup.name.lower()}.pypeit'

        if not by_hand_pypeit.exists():
            # It doesn't exist, so assume there is no by-hand pypeit
            # file to compare to
            continue

        # Run pypeit_setup
        ps = PypeItSetup.from_file_root(str(setup), 'keck_nires')
        ps.run(setup_only=True)

        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        # Read the frame types from both the by-hand and automated
        # pypeit files
        byhand_pypeItFile = PypeItFile.from_file(by_hand_pypeit)
        auto_pypeItFile = PypeItFile.from_file(pypeit_files[0])

        # For each file in the by-hand list, check that the frame types
        # in the automatically generated pypeit file are identical
        for f in byhand_pypeItFile.frametypes.keys():
            type_list = np.sort(byhand_pypeItFile.frametypes[f].split(','))
            if 'science' in type_list or 'standard' in type_list:
                # Only ensuring that calibrations are correctly typed
                continue
            assert f in auto_pypeItFile.frametypes.keys(), \
                'Frame {0} not automatically parsed for setup {1}.'.format(f, setup)
            assert np.array_equal(type_list, np.sort(auto_pypeItFile.frametypes[f].split(','))), \
                'Frame types differ for file {0} in setup {1}\n'.format(f, setup) \
                 + '    By-hand types: {0}'.format(byhand_pypeItFile.frametypes[f]) \
                 + '    Automated types: {0}'.format(auto_pypeItFile.frametypes[f])

        # Clean up after every setup
        shutil.rmtree(output_path)


def test_lris_red():
    # Raw LRIS RED directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'

    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_red'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the list of setup directories
    setups_names = ['long_600_7500_d560', 'multi_1200_9000_d680', 'multi_1200_9000_d680_1x2',
                    'long_400_8500_longread', 'long_600_10000_d680', 'long_150_7500_d560',
                    'long_300_5000_d560', 'long_1200_7500_d560', 'mulit_831_8200_d560', 'multi_900_5500_d560']
    setups = np.concatenate([sorted(raw_dir.glob(s)) for s in setups_names]).tolist()

    # Set the output path and *remove if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Iterate through the setups
    for setup in setups:
        # Find the relevant pypeit file constructed by hand.
        by_hand_pypeit = dev_root / 'pypeit_files' / f'keck_lris_red_{setup.name.lower()}.pypeit'

        if not by_hand_pypeit.exists():
            # It doesn't exist, so assume there is no by-hand pypeit
            # file to compare to
            continue

        # Run pypeit_setup
        ps = PypeItSetup.from_file_root(str(setup), 'keck_lris_red')
        ps.run(setup_only=True)

        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        # Read the frame types from both the by-hand and automated
        # pypeit files
        byhand_pypeItFile = PypeItFile.from_file(by_hand_pypeit)
        auto_pypeItFile = PypeItFile.from_file(pypeit_files[0])

        # For each file in the by-hand list, check that the frame types
        # in the automatically generated pypeit file are identical
        for f in byhand_pypeItFile.frametypes.keys():
            type_list = np.sort(byhand_pypeItFile.frametypes[f].split(','))
            if 'science' in type_list or 'standard' in type_list:
                # Only ensuring that calibrations are correctly typed
                continue
            assert f in auto_pypeItFile.frametypes.keys(), \
                'Frame {0} not automatically parsed for setup {1}.'.format(f, setup)
            assert np.array_equal(type_list, np.sort(auto_pypeItFile.frametypes[f].split(','))), \
                'Frame types differ for file {0} in setup {1}\n'.format(f, setup) \
                + '    By-hand types: {0}'.format(byhand_pypeItFile.frametypes[f]) \
                + '    Automated types: {0}'.format(auto_pypeItFile.frametypes[f])

        # Clean up after every setup
        shutil.rmtree(output_path)


def test_lris_red_orig():
    # Raw LRIS RED ORIG directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'

    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_red_orig'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the list of setup directories
    setups = list(raw_dir.glob('*'))

    # Set the output path and *remove if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Iterate through the setups
    for setup in setups:
        # Find the relevant pypeit file constructed by hand.
        by_hand_pypeit = dev_root / 'pypeit_files' / f'keck_lris_red_orig_{setup.name.lower()}.pypeit'

        if not by_hand_pypeit.exists():
            # It doesn't exist, so assume there is no by-hand pypeit
            # file to compare to
            continue

        # Run pypeit_setup
        ps = PypeItSetup.from_file_root(str(setup), 'keck_lris_red_orig')
        ps.run(setup_only=True)

        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        # Read the frame types from both the by-hand and automated
        # pypeit files
        byhand_pypeItFile = PypeItFile.from_file(by_hand_pypeit)
        auto_pypeItFile = PypeItFile.from_file(pypeit_files[0])

        # For each file in the by-hand list, check that the frame types
        # in the automatically generated pypeit file are identical
        for f in byhand_pypeItFile.frametypes.keys():
            type_list = np.sort(byhand_pypeItFile.frametypes[f].split(','))
            if 'science' in type_list or 'standard' in type_list:
                # Only ensuring that calibrations are correctly typed
                continue
            assert f in auto_pypeItFile.frametypes.keys(), \
                'Frame {0} not automatically parsed for setup {1}.'.format(f, setup)
            assert np.array_equal(type_list, np.sort(auto_pypeItFile.frametypes[f].split(','))), \
                'Frame types differ for file {0} in setup {1}\n'.format(f, setup) \
                + '    By-hand types: {0}'.format(byhand_pypeItFile.frametypes[f]) \
                + '    Automated types: {0}'.format(auto_pypeItFile.frametypes[f])

        # Clean up after every setup
        shutil.rmtree(output_path)


def test_lris_red_mark4():
    # Raw LRIS RED MARK 4 directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'

    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_red_mark4'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the list of setup directories
    setups_names = ['multi_600_10000_slitmask']
    setups = np.concatenate([sorted(raw_dir.glob(s)) for s in setups_names]).tolist()

    # Set the output path and *remove if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Iterate through the setups
    for setup in setups:
        # Find the relevant pypeit file constructed by hand.
        by_hand_pypeit = dev_root / 'pypeit_files' / f'keck_lris_red_mark4_{setup.name.lower()}.pypeit'

        if not by_hand_pypeit.exists():
            # It doesn't exist, so assume there is no by-hand pypeit
            # file to compare to
            continue

        # Run pypeit_setup
        ps = PypeItSetup.from_file_root(str(setup), 'keck_lris_red_mark4')
        ps.run(setup_only=True)

        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        # Read the frame types from both the by-hand and automated
        # pypeit files
        byhand_pypeItFile = PypeItFile.from_file(by_hand_pypeit)
        auto_pypeItFile = PypeItFile.from_file(pypeit_files[0])

        # For each file in the by-hand list, check that the frame types
        # in the automatically generated pypeit file are identical
        for f in byhand_pypeItFile.frametypes.keys():
            type_list = np.sort(byhand_pypeItFile.frametypes[f].split(','))
            if 'science' in type_list or 'standard' in type_list:
                # Only ensuring that calibrations are correctly typed
                continue
            assert f in auto_pypeItFile.frametypes.keys(), \
                'Frame {0} not automatically parsed for setup {1}.'.format(f, setup)
            assert np.array_equal(type_list, np.sort(auto_pypeItFile.frametypes[f].split(','))), \
                'Frame types differ for file {0} in setup {1}\n'.format(f, setup) \
                + '    By-hand types: {0}'.format(byhand_pypeItFile.frametypes[f]) \
                + '    Automated types: {0}'.format(auto_pypeItFile.frametypes[f])

        # Clean up after every setup
        shutil.rmtree(output_path)


def test_lris_blue():
    # Raw LRIS BLUE directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'

    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_blue'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the list of setup directories
    setups_names = ['multi_300_5000_d680', 'multi_600_4000_slitmask']
    setups = np.concatenate([sorted(raw_dir.glob(s)) for s in setups_names]).tolist()

    # Set the output path and *remove if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Iterate through the setups
    for setup in setups:
        # Find the relevant pypeit file constructed by hand.
        by_hand_pypeit = dev_root / 'pypeit_files' / f'keck_lris_blue_{setup.name.lower()}.pypeit'

        if not by_hand_pypeit.exists():
            # It doesn't exist, so assume there is no by-hand pypeit
            # file to compare to
            continue

        # Run pypeit_setup
        ps = PypeItSetup.from_file_root(str(setup), 'keck_lris_blue')
        ps.run(setup_only=True)

        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        # Read the frame types from both the by-hand and automated
        # pypeit files
        byhand_pypeItFile = PypeItFile.from_file(by_hand_pypeit)
        auto_pypeItFile = PypeItFile.from_file(pypeit_files[0])

        # For each file in the by-hand list, check that the frame types
        # in the automatically generated pypeit file are identical
        for f in byhand_pypeItFile.frametypes.keys():
            type_list = np.sort(byhand_pypeItFile.frametypes[f].split(','))
            if 'science' in type_list or 'standard' in type_list or 'pixelflat' in type_list:
                # Only ensuring that calibrations are correctly typed
                continue
            assert f in auto_pypeItFile.frametypes.keys(), \
                'Frame {0} not automatically parsed for setup {1}.'.format(f, setup)
            assert np.array_equal(type_list, np.sort(auto_pypeItFile.frametypes[f].split(','))), \
                'Frame types differ for file {0} in setup {1}\n'.format(f, setup) \
                + '    By-hand types: {0}'.format(byhand_pypeItFile.frametypes[f]) \
                + '    Automated types: {0}'.format(auto_pypeItFile.frametypes[f])

        # Clean up after every setup
        shutil.rmtree(output_path)


def test_lris_blue_orig():
    # Raw LRIS BLUE ORIG directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'

    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_blue_orig'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the list of setup directories
    setups = list(raw_dir.glob('*'))

    # Set the output path and *remove if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Iterate through the setups
    for setup in setups:
        # Find the relevant pypeit file constructed by hand.
        by_hand_pypeit = dev_root / 'pypeit_files' / f'keck_lris_blue_orig_{setup.name.lower()}.pypeit'

        if not by_hand_pypeit.exists():
            # It doesn't exist, so assume there is no by-hand pypeit
            # file to compare to
            continue

        # Run pypeit_setup
        ps = PypeItSetup.from_file_root(str(setup), 'keck_lris_blue_orig')
        ps.run(setup_only=True)

        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        # Read the frame types from both the by-hand and automated
        # pypeit files
        byhand_pypeItFile = PypeItFile.from_file(by_hand_pypeit)
        auto_pypeItFile = PypeItFile.from_file(pypeit_files[0])

        # For each file in the by-hand list, check that the frame types
        # in the automatically generated pypeit file are identical
        for f in byhand_pypeItFile.frametypes.keys():
            type_list = np.sort(byhand_pypeItFile.frametypes[f].split(','))
            if 'science' in type_list or 'standard' in type_list or 'pixelflat' in type_list:
                # Only ensuring that calibrations are correctly typed
                continue
            assert f in auto_pypeItFile.frametypes.keys(), \
                'Frame {0} not automatically parsed for setup {1}.'.format(f, setup)
            assert np.array_equal(type_list, np.sort(auto_pypeItFile.frametypes[f].split(','))), \
                'Frame types differ for file {0} in setup {1}\n'.format(f, setup) \
                + '    By-hand types: {0}'.format(byhand_pypeItFile.frametypes[f]) \
                + '    Automated types: {0}'.format(auto_pypeItFile.frametypes[f])

        # Clean up after every setup
        shutil.rmtree(output_path)


