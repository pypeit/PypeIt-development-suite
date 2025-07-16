"""
Module to run tests on scripts
"""
from pathlib import Path
import os
import glob
import shutil

from IPython import embed

import numpy as np

from configobj import ConfigObj

from pypeit.metadata import PypeItMetaData
from pypeit.par import PypeItPar
from pypeit.scripts.setup import Setup
from pypeit.scripts.chk_for_calibs import ChkForCalibs
from pypeit.spectrographs.util import load_spectrograph
from pypeit import pypeitsetup
from pypeit import inputfiles


def expected_file_extensions():
    return ['.sorted', '.obslog', '.pypeit', '.calib']


def files_are_expected(path):
    files = sorted(list(path.glob('*')))
    expected = expected_file_extensions()
    assert np.all([f.suffix in expected for f in files]), \
            'File produced that does not have an expected extension'


def generic_setup_test(spec, setup, cfg=None, prefix=None, extension=None):
    # Define the path with the raw data
    # TODO: Make the path structure of these instruments the same as the rest!!
    if 'vlt_xshooter' in spec:
        spec_dir = 'vlt_xshooter'
    else:
        spec_dir = spec
    data_root = Path(os.environ['PYPEIT_DEV']).resolve() / 'RAW_DATA' / spec_dir / setup
    assert data_root.exists(), 'TEST ERROR: Raw data path does not exist'
    if prefix is not None:
        data_root /= prefix

    # Define the output directory and remove it if it already exist
    setup_path = Path().resolve() / (f'{spec}_A' if cfg else 'setup_files')
    if setup_path.exists():
        shutil.rmtree(setup_path)

    args = ['-r', str(data_root), '-s', spec]
    if extension is not None:
        args += ['--extension', extension]
    if cfg is not None:
        args += ['-c', cfg]
    pargs = Setup.parse_args(args)
    Setup.main(pargs)

    assert setup_path.exists(), 'No setup_files directory created'
    files_are_expected(setup_path)

    # Clean-up
    shutil.rmtree(setup_path)


def test_setup_keck_lris_red_mark4():
    spec = 'keck_lris_red_mark4'
    setup = 'long_400_8500_d560'
    generic_setup_test(spec, setup)


def test_setup_keck_lris_red_mark4_multiconfig():
    # Dev suite directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'
    # Raw data directory
    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_red_mark4'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the fits files
    datasets = ['long_400_8500_d560']
    files = np.concatenate([sorted(raw_dir.glob(f'{s}/*fits*')) for s in datasets]).tolist()

    # Set the output path and remove if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Run pypeit_setup
    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_lris_red_mark4')
    ps.run(setup_only=True)

    # Write the automatically generated pypeit data
    pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

    assert len(pypeit_files) == 2, 'Should have created two pypeit files'

    # Test the pypeit files for the correct configuration and
    # calibration group results
    for f, s, c in zip(pypeit_files, ['A', 'B'], ['0', '1']):
        # Read the pypeit file
        pypeitFile = inputfiles.PypeItFile.from_file(f)
        # Check setup name
        assert pypeitFile.setup_name == s, 'Setup is wrong'
        # check calibration group
        assert np.all(pypeitFile.data['calib'].astype(str) == c), 'Calibration group is wrong'
    # Clean-up
    shutil.rmtree(output_path)


def test_setup_keck_lris_red():
    spec = 'keck_lris_red'
    setup = 'multi_400_8500_d560'
    generic_setup_test(spec, setup)


def test_setup_keck_lris_red_multiconfig():
    # Dev suite directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'
    # Raw data directory
    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_red'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the fits files
    datasets = ['long_150_7500_d560', 'long_1200_7500_d560']
    files = np.concatenate([sorted(raw_dir.glob(f'{s}/*fits*')) for s in datasets]).tolist()

    # Set the output path and remove if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Run pypeit_setup
    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_lris_red')
    ps.run(setup_only=True)

    # Write the automatically generated pypeit data
    pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

    assert len(pypeit_files) == 2, 'Should have created two pypeit files'

    # Test the pypeit files for the correct configuration and
    # calibration group results
    for d, f, s, c in zip(datasets, pypeit_files, ['A', 'B'], ['0', '1']):
        # Read the pypeit file
        pypeitFile = inputfiles.PypeItFile.from_file(f)
        # Check setup name
        assert pypeitFile.setup_name == s, 'Setup is wrong'
        # check calibration group
        assert np.all(pypeitFile.data['calib'].astype(str) == c), 'Calibration group is wrong'
        # check that this is the right dataset
        assert Path(pypeitFile.file_paths[0]).name == d, 'Wrong dataset'
    # Clean-up
    shutil.rmtree(output_path)

def test_setup_keck_lris_red_orig():
    spec = 'keck_lris_red_orig'
    setup = 'long_300_5000'
    generic_setup_test(spec, setup)


def test_setup_keck_lris_red_orig_multiconfig():
    # Dev suite directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'
    # Raw data directory
    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_red_orig'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the fits files
    datasets = ['long_900_5500_d560', 'long_831_8200_d460']
    files = np.concatenate([sorted(raw_dir.glob(f'{s}/*fits*')) for s in datasets]).tolist()

    # Set the output path and remove if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Run pypeit_setup
    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_lris_red_orig')
    ps.run(setup_only=True)

    # Write the automatically generated pypeit data
    pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

    assert len(pypeit_files) == 2, 'Should have created two pypeit files'

    # Test the pypeit files for the correct configuration and
    # calibration group results
    for d, f, s, c in zip(datasets, pypeit_files, ['A', 'B'], ['0', '1']):
        # Read the pypeit file
        pypeitFile = inputfiles.PypeItFile.from_file(f)
        # Check setup name
        assert pypeitFile.setup_name == s, 'Setup is wrong'
        # check calibration group
        assert np.all(pypeitFile.data['calib'].astype(str) == c), 'Calibration group is wrong'
        # check that this is the right dataset
        assert Path(pypeitFile.file_paths[0]).name == d, 'Wrong dataset'
    # Clean-up
    shutil.rmtree(output_path)


def test_setup_keck_lris_blue():
    spec = 'keck_lris_blue'
    setup = 'multi_600_4000_d560'
    generic_setup_test(spec, setup)


def test_setup_keck_lris_blue_multiconfig():
    # Dev suite directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'
    # Raw data directory
    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_blue'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the fits files
    datasets = ['multi_300_5000_d680', 'multi_600_4000_slitmask']
    files = np.concatenate([sorted(raw_dir.glob(f'{s}/*fits*')) for s in datasets]).tolist()

    # Set the output path and remove if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Run pypeit_setup
    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_lris_blue')
    ps.run(setup_only=True)

    # Write the automatically generated pypeit data
    pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

    assert len(pypeit_files) == 2, 'Should have created two pypeit files'

    # Test the pypeit files for the correct configuration and
    # calibration group results
    for d, f, s, c in zip(datasets, pypeit_files, ['A', 'B'], ['0', '1']):
        # Read the pypeit file
        pypeitFile = inputfiles.PypeItFile.from_file(f)
        # Check setup name
        assert pypeitFile.setup_name == s, 'Setup is wrong'
        # check calibration group
        assert np.all(pypeitFile.data['calib'].astype(str) == c), 'Calibration group is wrong'
        # check that this is the right dataset
        assert Path(pypeitFile.file_paths[0]).name == d, 'Wrong dataset'
    # Clean-up
    shutil.rmtree(output_path)


def test_setup_keck_lris_blue_orig():
    spec = 'keck_lris_blue_orig'
    setup = 'long_600_4000_d500'
    generic_setup_test(spec, setup)


def test_setup_keck_lris_blue_orig_multiconfig():
    # Dev suite directory
    dev_root = Path(os.getenv('PYPEIT_DEV')).resolve()
    assert dev_root.exists(), f'PypeIt development suite directory does not exist: {dev_root}'
    # Raw data directory
    raw_dir = dev_root / 'RAW_DATA' / 'keck_lris_blue_orig'
    assert raw_dir.exists(), f'Raw data directory does not exist: {raw_dir}'

    # Get the fits files
    datasets = ['long_600_4000_d500', 'multi_1200_3400_d460']
    files = np.concatenate([sorted(raw_dir.glob(f'{s}/*fits*')) for s in datasets]).tolist()

    # Set the output path and remove if it already exists
    output_path = Path('.').resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)

    # Run pypeit_setup
    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_lris_blue_orig')
    ps.run(setup_only=True)

    # Write the automatically generated pypeit data
    pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

    assert len(pypeit_files) == 2, 'Should have created two pypeit files'

    # Test the pypeit files for the correct configuration and
    # calibration group results
    for d, f, s, c in zip(datasets, pypeit_files, ['A', 'B'], ['0', '1']):
        # Read the pypeit file
        pypeitFile = inputfiles.PypeItFile.from_file(f)
        # Check setup name
        assert pypeitFile.setup_name == s, 'Setup is wrong'
        # check calibration group
        assert np.all(pypeitFile.data['calib'].astype(str) == c), 'Calibration group is wrong'
        # check that this is the right dataset
        assert Path(pypeitFile.file_paths[0]).name == d, 'Wrong dataset'
    # Clean-up
    shutil.rmtree(output_path)


def test_setup_shane_kast_blue():
    spec = 'shane_kast_blue'
    setup = '600_4310_d55'
    generic_setup_test(spec, setup)
    generic_setup_test(spec, setup, cfg='all')


def test_setup_shane_kast_red():
    spec = 'shane_kast_red'
    setup = '600_7500_d55_ret'
    generic_setup_test(spec, setup)


def test_setup_keck_deimos():
    spec = 'keck_deimos'
    setup = '830G_M_8600'
    generic_setup_test(spec, setup)


def test_setup_keck_deimos_multiconfig():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos')
    files = glob.glob(os.path.join(root, '830G_L_8100', '*fits*'))
    files += glob.glob(os.path.join(root, '830G_L_8400', '*fits*'))

    output_path = os.path.join(os.getcwd(), 'output')
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)

    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_deimos')
    ps.run(setup_only=True)

    # Write the automatically generated pypeit data
    pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg,
                                           write_bkg_pairs=True)

    assert len(pypeit_files) == 2, 'Should have created two pypeit files'

    # Test the pypeit files for the correct configuration and
    # calibration group results
    for f, s, c in zip(pypeit_files, ['A', 'B'], ['0', '1']):

        # TODO: All of this front-end stuff, pulled from pypeit.py, should
        # be put into a function.

        # Read the pypeit file
        pypeitFile = inputfiles.PypeItFile.from_file(f)
        # Spectrograph
        cfg = ConfigObj(pypeitFile.cfg_lines)
        spectrograph = load_spectrograph(cfg['rdx']['spectrograph'])
        # Configuration-specific parameters
        for idx, row in enumerate(pypeitFile.data):
            if 'science' in row['frametype'] or 'standard' in row['frametype']:
                break
        spectrograph_cfg_lines = spectrograph.config_specific_par(
            pypeitFile.filenames[idx]).to_config()
        #  PypeIt parameters
        par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, 
                                       merge_with=pypeitFile.cfg_lines)
        #  Metadata
        fitstbl = PypeItMetaData(spectrograph, par, 
                                 files=pypeitFile.filenames, 
                                 usrdata=pypeitFile.data, 
                                 strict=True)
        fitstbl.finalize_usr_build(pypeitFile.frametypes, pypeitFile.setup_name)

        assert np.all(fitstbl['setup'] == s), 'Setup is wrong'
        assert np.all(fitstbl['calib'].astype(str) == c), 'Calibration group is wrong'

    # Clean-up
    shutil.rmtree(output_path)


def test_setup_keck_deimos_multiconfig_clean():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos')
    files = glob.glob(os.path.join(root, '830G_L_8100', '*fits*'))
    files += glob.glob(os.path.join(root, '830G_L_8400', '*fits*'))

    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_deimos')
    ps.build_fitstbl(strict=False)
    ps.get_frame_types(flag_unknown=True)

    # Test that the correct number of configurations are found
    cfgs = ps.fitstbl.unique_configurations()
    assert len(cfgs) == 2, 'Should find 2 configurations'

    # Test that the bias is assigned to the correct configuration
    ps.fitstbl.set_configurations(cfgs)
    biases = np.where(ps.fitstbl.find_frames('bias'))[0]
    assert biases.size == 1, 'Should only be 1 bias'
    assert ps.fitstbl['setup'][biases[0]] == 'B', 'Bias should be in configuration group B'

    # Table should have 25 rows
    assert len(ps.fitstbl) == 25, 'Incorrect number of table rows.'

    # All frames should be from valid configurations
    ps.fitstbl.clean_configurations()
    assert len(ps.fitstbl) == 25, 'Incorrect number of table rows.'

    # Artificially set the amplifier and mode of two frames to be
    # invalid
    ps.fitstbl['amp'][0] = 'DUAL:A+B'
    ps.fitstbl['mode'][1] = 'Direct'
    ps.fitstbl.clean_configurations()
    # Those two frames should have been removed
    assert len(ps.fitstbl) == 23, 'Incorrect number of table rows.'


def test_setup_keck_mosfire():
    spec = 'keck_mosfire'
    setup = 'J_multi'
    generic_setup_test(spec, setup)


def test_setup_keck_mosfire_multiconfig():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_mosfire')
    files = glob.glob(os.path.join(root, 'K_long', '*fits*'))
    files += glob.glob(os.path.join(root, 'long2pos1_H', '*fits*'))
    files += glob.glob(os.path.join(root, 'mask1_K_with_continuum', '*fits*'))
    files += glob.glob(os.path.join(root, 'Y_multi', '*fits*'))

    output_path = os.path.join(os.getcwd(), 'output')
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)

    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_mosfire')
    ps.run(setup_only=True) 
    # Write the automatically generated pypeit data
    pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg,
                                           write_bkg_pairs=True)

    assert len(pypeit_files) == 4, 'Should have created four pypeit files'

    # Test the pypeit files for the correct configuration,
    # calibration group and combination group results
    # expected values
    setups = ['A', 'B', 'C', 'D']
    calib_ids = ['0', '1', '2', '3']
    abba_dpat = [1,2,2,1], [2,1,1,2]
    long2pos_dpat = [5,6,7,8], [6,5,8,7]
    abab_dpat = [9,10,11,12],[10,9,12,11]   # tihs is ABA'B'
    masknod_dpat = [13,14,15,16], [14,13,16,15]
    for f, s, c, comb, bkg in zip(pypeit_files,setups, calib_ids,
                                  [abba_dpat[0], long2pos_dpat[0], abab_dpat[0], masknod_dpat[0]],
                                  [abba_dpat[1], long2pos_dpat[1], abab_dpat[1], masknod_dpat[1]]):

        # TODO: All of this front-end stuff, pulled from pypeit.py, should
        # be put into a function.

        # Read the pypeit file
        pypeitFile = inputfiles.PypeItFile.from_file(f)
        # Spectrograph
        cfg = ConfigObj(pypeitFile.cfg_lines)
        spectrograph = load_spectrograph(cfg['rdx']['spectrograph'])
        # Configuration-specific parameters
        for idx, row in enumerate(pypeitFile.data):
            if 'science' in row['frametype'] or 'standard' in row['frametype']:
                # assume there is always a science/standard for this test
                config_specific_file = pypeitFile.filenames[idx]
        spectrograph_cfg_lines = spectrograph.config_specific_par(config_specific_file).to_config()
        #  PypeIt parameters
        par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines,
                                       merge_with=pypeitFile.cfg_lines)
        #  Metadata
        fitstbl = PypeItMetaData(spectrograph, par,
                                 files=pypeitFile.filenames,
                                 usrdata=pypeitFile.data,
                                 strict=True)
        fitstbl.finalize_usr_build(pypeitFile.frametypes, pypeitFile.setup_name)

        # Check setup
        assert np.all(fitstbl['setup'] == s), 'Setup is wrong'
        # Check calibration group
        assert np.all(fitstbl['calib'].astype(str) == c), 'Calibration group is wrong'
        # Check combination and background group for only science/standard
        sci_std_idx = np.array(['science' in _tab or 'standard' in _tab
                                    for _tab in fitstbl['frametype']]) & (fitstbl['setup'] == s)
        assert np.all(fitstbl['comb_id'][sci_std_idx] == comb), 'Combination group is wrong'
        assert np.all(fitstbl['bkg_id'][sci_std_idx] == bkg), 'Background group is wrong'

    # Clean-up
    shutil.rmtree(output_path)


def test_setup_keck_hires_multiconfig():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_hires')
    files = glob.glob(os.path.join(root, 'J0306+1853_U074_RED_C2_ECH_0.72_XD_1.42_1x3', '*fits*'))
    files += glob.glob(os.path.join(root, 'J1218+2951_U116Hr_RED_C5_ECH_-0.22_XD_0.21_1x2', '*fits*'))

    output_path = os.path.join(os.getcwd(), 'output')
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)

    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_hires')
    ps.run(setup_only=True)

    # Write the automatically generated pypeit data
    pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg,
                                           write_bkg_pairs=True)

    assert len(pypeit_files) == 2, 'Should have created two pypeit files'

    # Test the pypeit files for the correct configuration and
    # calibration group results
    for f, s, c in zip(pypeit_files, ['A', 'B'], ['0', '1']):

        # TODO: All of this front-end stuff, pulled from pypeit.py, should
        # be put into a function.

        # Read the pypeit file
        pypeitFile = inputfiles.PypeItFile.from_file(f)
        # Spectrograph
        cfg = ConfigObj(pypeitFile.cfg_lines)
        spectrograph = load_spectrograph(cfg['rdx']['spectrograph'])
        # Configuration-specific parameters
        for idx, row in enumerate(pypeitFile.data):
            if 'science' in row['frametype'] or 'standard' in row['frametype']:
                break
        spectrograph_cfg_lines = spectrograph.config_specific_par(
            pypeitFile.filenames[idx]).to_config()
        #  PypeIt parameters
        par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines,
                                       merge_with=pypeitFile.cfg_lines)
        #  Metadata
        fitstbl = PypeItMetaData(spectrograph, par,
                                 files=pypeitFile.filenames,
                                 usrdata=pypeitFile.data,
                                 strict=True)
        fitstbl.finalize_usr_build(pypeitFile.frametypes, pypeitFile.setup_name)

        assert np.all(fitstbl['setup'] == s), 'Setup is wrong'
        assert np.all(fitstbl['calib'].astype(str) == c), 'Calibration group is wrong'

    # Clean-up
    shutil.rmtree(output_path)


def test_setup_keck_nires():
    spec = 'keck_nires'
    setup = 'ABBA_wstandard'
    generic_setup_test(spec, setup)



def test_setup_keck_nires_comb():

    dev_root = Path(os.environ['PYPEIT_DEV']).resolve()
    data_root = dev_root / 'RAW_DATA' / 'keck_nires'

    output_path = Path().resolve() / 'output'
    if output_path.exists():
        shutil.rmtree(output_path)
    output_path.mkdir()

    datasets = ['ABpat_wstandard', 'ABC_nostandard', 'ABBA_nostandard']
    for dset in datasets:
        files = sorted(list((data_root / dset).glob('*fits*')))
        correct_pypeit_file = dev_root / 'pypeit_files' / f'keck_nires_{dset.lower()}.pypeit'

        # run setup on raw files
        ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_nires')
        ps.run(setup_only=True)
        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg,
                                               write_bkg_pairs=True)


        # Test that the automatically generated configuration, calibration group, combination group
        # and background group for science frames are correct, i.e., are the same as in the pre-generated pypeit file.

        # check that only one setup is generated
        assert np.all(ps.fitstbl['setup'] == 'A'), 'Should have created one setup'

        # Read the correct pypeit file
        pypeitFile = inputfiles.PypeItFile.from_file(correct_pypeit_file)

        auto_science = ps.fitstbl['frametype'] == 'arc,science,tilt'
        correct_science = pypeitFile.data['frametype'] == 'arc,science,tilt'
        science_filenames = pypeitFile.data[correct_science]['filename'].data
        for i in range(science_filenames.size):
            where_this = ps.fitstbl['filename'] == science_filenames[i]
            # Check this file exists
            assert np.any(where_this), 'Science file does not exist in the correct pypeit file'
            # Check calibration group
            assert ps.fitstbl['calib'][where_this].data[0] == pypeitFile.data[correct_science]['calib'][i], \
                'Calibration group is wrong'
            # Check combination and background group
            assert ps.fitstbl['comb_id'][where_this].data[0] == int(pypeitFile.data[correct_science]['comb_id'][i]), \
                'Combination group is wrong'
            assert ps.fitstbl['bkg_id'][where_this].data[0] == int(pypeitFile.data[correct_science]['bkg_id'][i]), \
                'Background group is wrong'


def test_setup_keck_nirspec():
    spec = 'keck_nirspec_low'
    setup = 'LOW_NIRSPEC-1'
    generic_setup_test(spec, setup)


def test_setup_magellan_mage():
    spec = 'magellan_mage'
    setup = '1x1'
    generic_setup_test(spec, setup)

def test_setup_magellan_fire():
    spec = 'magellan_fire'
    setup = 'FIRE'
    generic_setup_test(spec, setup)


def test_setup_wht_isis_blue():
    spec = 'wht_isis_blue'
    setup = 'long_R300B_d5300'
    generic_setup_test(spec, setup)


def test_setup_vlt_xshooter_uvb():
    spec = 'vlt_xshooter_uvb'
    setup = 'UVB_1x1'
    prefix = 'XSHO'
    generic_setup_test(spec, setup, prefix=prefix)


def test_setup_vlt_xshooter_vis():
    spec = 'vlt_xshooter_vis'
    setup = 'VIS_1x1'
    prefix = 'XSHO'
    generic_setup_test(spec, setup, prefix=prefix)


def test_setup_vlt_xshooter_nir():
    spec = 'vlt_xshooter_nir'
    setup = 'NIR'
    prefix = 'XSHO'
    generic_setup_test(spec, setup, prefix=prefix)


def test_setup_gemini_gnirs():
    spec = 'gemini_gnirs_echelle'
    setup = '32_SB_SXD'
    prefix = 'cN'
    generic_setup_test(spec, setup, prefix=prefix)


def test_setup_not_alfosc():
    spec = 'not_alfosc'
    setup = 'grism4'
    prefix = 'ALD'
    generic_setup_test(spec, setup, prefix=prefix)
    generic_setup_test(spec, setup, cfg='A', prefix=prefix)


def test_setup_vlt_fors2():
    spec = 'vlt_fors2'
    setup = '300I'
    prefix = 'FORS2'
    generic_setup_test(spec, setup, prefix=prefix)

    # Now chk calib
    data_root = Path(os.environ['PYPEIT_DEV']).resolve() / 'RAW_DATA' / spec / setup / prefix
    pargs = ChkForCalibs.parse_args([str(data_root), '-s', 'vlt_fors2'])
    answers, ps = ChkForCalibs.main(pargs)
    assert answers['pass'][0], 'A must pass!'


def test_setup_ldt_deveny():
    spec = 'ldt_deveny'
    setup = 'DV1'
    generic_setup_test(spec, setup)
    # TODO: Think about how to test that all DeVeny setups are being detected


def test_setup_param_block():
    # Define the output directory and remove it if it already exist
    setup_path = Path().resolve() / 'ldt_deveny_A'
    if setup_path.exists():
        shutil.rmtree(setup_path)

    # Test this with LDT/DeVeny::DV6
    data_root = Path(os.environ['PYPEIT_DEV']).resolve() / 'RAW_DATA' / 'ldt_deveny' / 'DV6'
    assert data_root.exists(), 'TEST ERROR: Raw data path does not exist'

    # Test the ability to read in the extra parameters
    parblock_fn = Path(os.environ['PYPEIT_DEV']).resolve() / 'pypeit_files' / 'ldt_deveny_xtra_params.txt'
    args = ['-r', str(data_root), '-s', 'ldt_deveny', '-p', str(parblock_fn) , '-c' 'A']
    pargs = Setup.parse_args(args)
    Setup.main(pargs)

    # Read in the xtra_pars and the created PypeIt file
    with open(parblock_fn, 'r', encoding='utf-8') as par_fobj:
        xtra_pars = [l.rstrip() for l in par_fobj.readlines()]
    with open(setup_path / 'ldt_deveny_A.pypeit', 'r',encoding='utf-8') as pypeit_fobj:
        pypeit_contents = [l.rstrip() for l in pypeit_fobj.readlines()]

    # Check that each of the `xtra_pars` is in the created PypeIt file
    for par in xtra_pars:
        assert par in pypeit_contents

    # Clean-up
    shutil.rmtree(setup_path)


# TODO: Add other instruments!

