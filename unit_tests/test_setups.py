"""
Module to run tests on scripts
"""
import os
import glob
import shutil

from IPython import embed

import numpy as np

import pytest
from configobj import ConfigObj

from pypeit.pypmsgs import PypeItError
from pypeit.metadata import PypeItMetaData
from pypeit.par import PypeItPar
from pypeit.scripts.setup import Setup
from pypeit.scripts.chk_for_calibs import ChkForCalibs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests.tstutils import data_path
from pypeit import pypeit
from pypeit import pypeitsetup
from pypeit import inputfiles


def expected_file_extensions():
    return ['sorted']


def test_setup_keck_lris_red_mark4():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_lris_red_mark4/long_400_8500_d560')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'keck_lris_red_mark4'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_red_mark4*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_keck_lris_red():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_lris_red/multi_400_8500_d560')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'keck_lris_red'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_red*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
        'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_keck_lris_red_orig():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_lris_red_orig/long_300_5000')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'keck_lris_red_orig'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_red_orig*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_keck_lris_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_lris_blue/multi_600_4000_d560')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'keck_lris_blue'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_blue*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
        'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_keck_lris_blue_orig():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_lris_blue_orig/long_600_4000_d500')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'keck_lris_blue_orig'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_lris_blue_orig*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_shane_kast_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/shane_kast_blue/600_4310_d55')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'shane_kast_blue'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'shane_kast_blue*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_shane_kast_red():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/shane_kast_red/600_7500_d55_ret')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'shane_kast_red'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'shane_kast_red*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)

# TODO: We need a test data set for shane_kast_red_ret

def test_setup_keck_deimos():

    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_deimos/830G_M_8600')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'keck_deimos'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_deimos*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_keck_deimos_multiconfig():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos')
    files = glob.glob(os.path.join(root, '830G_L_8100', '*fits*'))
    files += glob.glob(os.path.join(root, '830G_L_8400', '*fits*'))

    output_path = os.path.join(os.getcwd(), 'output')
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)

    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_deimos')
    ps.run(setup_only=True, sort_dir=output_path)
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
        #cfg_lines, data_files, frametype, usrdata, setups, _ = parse_pypeit_file(f, runtime=True)
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
    ps.fitstbl['amp'][0] = 'SINGLE:A'
    ps.fitstbl['mode'][1] = 'Direct'
    ps.fitstbl.clean_configurations()
    # Those two frames should have been removed
    assert len(ps.fitstbl) == 23, 'Incorrect number of table rows.'


def test_setup_keck_mosfire():

    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_mosfire/J_multi')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'keck_mosfire'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_mosfire*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


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
    ps.run(setup_only=True, sort_dir=output_path, write_bkg_pairs=True)
    # Write the automatically generated pypeit data
    pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg,
                                           write_bkg_pairs=True)

    assert len(pypeit_files) == 4, 'Should have created two pypeit files'

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
        sci_std_idx = np.array(['science' in _tab or 'standard' in _tab for _tab in fitstbl['frametype']]) & \
                      (fitstbl['setup'] == s)
        assert np.all(fitstbl['comb_id'][sci_std_idx] == comb), 'Combination group is wrong'
        assert np.all(fitstbl['bkg_id'][sci_std_idx] == bkg), 'Background group is wrong'

    # Clean-up
    shutil.rmtree(output_path)


def test_setup_keck_nires():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_nires/NIRES/')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'keck_nires'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_nires*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_keck_nires_comb():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_nires')
    datasets = ['ABpat_wstandard', 'ABC_nostandard', 'ABBA_nostandard']
    for dset in datasets:
        files = glob.glob(os.path.join(root, dset, '*fits*'))
        correct_pypeit_file = os.path.join(os.environ['PYPEIT_DEV'], 'pypeit_files', f'keck_nires_{dset.lower()}.pypeit')

        # run setup on raw files
        ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_nires')
        ps.run(setup_only=True, write_bkg_pairs=True)

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
            assert ps.fitstbl['calib'][where_this].data[0] == int(pypeitFile.data[correct_science]['calib'][i]), \
                'Calibration group is wrong'
            # Check combination and background group
            assert ps.fitstbl['comb_id'][where_this].data[0] == int(pypeitFile.data[correct_science]['comb_id'][i]), \
                'Combination group is wrong'
            assert ps.fitstbl['bkg_id'][where_this].data[0] == int(pypeitFile.data[correct_science]['bkg_id'][i]), \
                'Background group is wrong'


def test_setup_keck_nirspec():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_nirspec/LOW_NIRSPEC-1')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'keck_nirspec_low'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'keck_nirspec*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_magellan_mage():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/magellan_mage/1x1')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'magellan_mage'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'magellan_mage*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_wht_isis_blue():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/wht_isis_blue/long_R300B_d5300')
    droot += '/'
    pargs = Setup.parse_args(['-r', droot, '-s', 'wht_isis_blue', '--extension', '.fit'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'wht_isis_blue*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_vlt_xshooter_uvb():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_xshooter/UVB_1x1')
    droot += '/XSHO'
    pargs = Setup.parse_args(['-r', droot, '-s', 'vlt_xshooter_uvb'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_uvb*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_vlt_xshooter_vis():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_xshooter/VIS_1x1')
    droot += '/XSHO'
    pargs = Setup.parse_args(['-r', droot, '-s', 'vlt_xshooter_vis'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_vis*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_vlt_xshooter_nir():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_xshooter/NIR')
    droot += '/XSHO'
    pargs = Setup.parse_args(['-r', droot, '-s', 'vlt_xshooter_nir'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_xshooter_nir*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_gemini_gnirs():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/gemini_gnirs/32_SB_SXD/')
    droot += '/cN'
    pargs = Setup.parse_args(['-r', droot, '-s', 'gemini_gnirs'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'gemini_gnirs*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
            'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)


def test_setup_not_alfosc():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/not_alfosc/grism4')
    droot += '/ALD'
    pargs = Setup.parse_args(['-r', droot, '-s', 'not_alfosc'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'not_alfosc*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
        'Did not find all setup file extensions: {0}'.format(expected)

    # Build a PypeIt file
    pargs = Setup.parse_args(['-r', droot, '-s', 'not_alfosc', '-c', 'A', '-d', data_path('')])
    Setup.main(pargs)
    pypeit_file = data_path('not_alfosc_A/not_alfosc_A.pypeit')
    # TODO: Why is this using pypeit.PypeIt and not pypeitsetup.PypeItSetup?
    pypeIt = pypeit.PypeIt(pypeit_file, calib_only=True)

    # Clean-up
    shutil.rmtree(setup_dir)
    shutil.rmtree(data_path('not_alfosc_A'))

def test_setup_vlt_fors2():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_fors2/300I/')
    droot += '/FORS2'
    pargs = Setup.parse_args(['-r', droot, '-s', 'vlt_fors2'])
    Setup.main(pargs)

    cwd = os.getcwd()
    setup_dir = os.path.join(cwd, 'setup_files')
    assert os.path.isdir(setup_dir), 'No setup_files directory created'

    files = glob.glob(os.path.join(setup_dir, 'vlt_fors2*'))
    ext = [f.split('.')[-1] for f in files]
    expected = expected_file_extensions()
    assert np.all([e in ext for e in expected]), \
        'Did not find all setup file extensions: {0}'.format(expected)

    # Clean-up
    shutil.rmtree(setup_dir)

    # Now chk calib
    pargs = ChkForCalibs.parse_args([droot, '-s', 'vlt_fors2'])
    answers, ps = ChkForCalibs.main(pargs)
    assert answers['pass'][0], 'A must pass!'

# TODO: Add other instruments!

