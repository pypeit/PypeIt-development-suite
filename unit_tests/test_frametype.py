
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


def test_vet_assigned_ftypes():

    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'

    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_blue' / 'long_600_4000_d560_slitless'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    raw_files = list(raw_dir.glob('*'))
    assert len(raw_files) > 0, f'No raw files found in {raw_dir}'

    ps = PypeItSetup.from_file_root(str(raw_dir), 'keck_lris_blue')
    ps.run(setup_only=True)

    # CHECK N.1
    # check that if standard and science frames are both assigned to a frame, vet_assigned_ftypes
    # will keep the right frame type
    # find standard star
    std_indx = np.where(ps.fitstbl['target'] == 'GD71')[0]
    # check that the standard star is assigned the correct frame type (standard and not science)
    assert ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][std_indx], flag='standard') and \
        not ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][std_indx], flag='science'), \
        'Standard star not assigned the correct frame type.'

    # modify the standard star frame bits to add a science frame type
    ps.fitstbl['framebit'][std_indx] = ps.fitstbl.type_bitmask.turn_on(ps.fitstbl['framebit'][std_indx], flag='science')

    assert ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][std_indx], flag='standard') and \
        ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][std_indx], flag='science'), \
        'Science framebit not modified correctly.'

    # run the vet_assigned_ftypes method
    ps.spectrograph.vet_assigned_ftypes(ps.fitstbl['framebit'], ps.fitstbl)

    # recheck that the standard star is assigned the correct frame type (standard and not science)
    assert ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][std_indx], flag='standard') and \
        not ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][std_indx], flag='science'), \
        'Standard star not assigned the correct frame type.'

    # CHECK N.2
    # let's do the opposite now, add a standard frame type to a science frame and check that the frame type is corrected
    # find science frame
    sci_indx = np.where(ps.fitstbl['target'] == 'DLA_J1054_off')[0]
    # check that the science frame is assigned the correct frame type (science and not standard)
    assert ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][sci_indx], flag='science') and \
        not ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][sci_indx], flag='standard'), \
        'Science frame not assigned the correct frame type.'

    # modify the science frame bits to add a standard frame type
    ps.fitstbl['framebit'][sci_indx] = ps.fitstbl.type_bitmask.turn_on(ps.fitstbl['framebit'][sci_indx], flag='standard')

    assert ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][sci_indx], flag='science') and \
        ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][sci_indx], flag='standard'), \
        'Science framebit not modified correctly.'

    # run the vet_assigned_ftypes method
    ps.spectrograph.vet_assigned_ftypes(ps.fitstbl['framebit'], ps.fitstbl)

    # recheck that the science frame is assigned the correct frame type (science and not standard)
    assert ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][sci_indx], flag='science') and \
        not ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][sci_indx], flag='standard'), \
        'Science frame not assigned the correct frame type.'

    # CHECK N.3
    # check if this dataset has both a pixelflat and slitless_pixflat frame assigned, the vet_assigned_ftypes
    # will keep the slitless frame type (currently specific for LRIS)

    # do we have slitless_pixflat frames?
    slitless_frames = ps.fitstbl.find_frames('slitless_pixflat')
    assert np.any(slitless_frames), 'No slitless_pixflat frames found.'

    # add the pixelflat frame type to a illumflat frame
    ill_indx = np.where(ps.fitstbl.find_frames('illumflat'))[0]
    ps.fitstbl['framebit'][ill_indx] = ps.fitstbl.type_bitmask.turn_on(ps.fitstbl['framebit'][ill_indx], flag='pixelflat')

    # check
    assert np.all(ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][ill_indx], flag='pixelflat')), \
        'Pixelflat framebit not added correctly.'

    # run the vet_assigned_ftypes method
    ps.spectrograph.vet_assigned_ftypes(ps.fitstbl['framebit'], ps.fitstbl)

    # recheck that the pixelflat frame type is removed from the illumflat frame
    assert not np.all(ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][ill_indx], flag='pixelflat')), \
        'Pixelflat framebit not removed correctly.'

    # CHECK N.4
    # what happen if ra and/or dec are None in the case of both standard and science type assigned to a frame?

    # find standard star
    std = np.where(ps.fitstbl.find_frames('standard'))[0]
    # check that ra and dec are not None
    assert np.all(ps.fitstbl['ra'][std] != None) and np.all(ps.fitstbl['dec'][std] != None), \
        'RA and/or DEC are None for standard star.'

    # modify the standard star frame bits to add a science frame type
    ps.fitstbl['framebit'][std] = ps.fitstbl.type_bitmask.turn_on(ps.fitstbl['framebit'][std], flag='science')

    # set ra and dec to None
    ps.fitstbl['ra'][std] = None
    ps.fitstbl['dec'][std] = None

    # run the vet_assigned_ftypes method
    ps.spectrograph.vet_assigned_ftypes(ps.fitstbl['framebit'], ps.fitstbl)

    # if ra and dec are None, PypeIt cannot determine if the frame is a standard, so it assumes it is a science frame
    assert np.all(ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][std], flag='science')), \
        'Frames with None ra and dec are not assigned the correct frame type.'

    # CHECK N.5
    # what happen if ra and dec are not part of ps.fitstbl.keys()?
    # restore the standard star frame bits
    ps.fitstbl['framebit'][std] = ps.fitstbl.type_bitmask.turn_on(ps.fitstbl['framebit'][std], flag='standard')

    # remove ra and dec from the table
    ps.fitstbl.table.remove_columns(['ra', 'dec'])

    # run the vet_assigned_ftypes method
    ps.spectrograph.vet_assigned_ftypes(ps.fitstbl['framebit'], ps.fitstbl)

    # if ra and dec are not part of the table, PypeIt cannot determine if the frame is a standard,
    # so it assumes it is a science frame
    assert np.all(ps.fitstbl.type_bitmask.flagged(ps.fitstbl['framebit'][std], flag='science')), \
        'Frames without ra and dec are not assigned the correct frame type.'
