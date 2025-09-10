"""
Module to run tests on scripts
"""
import os
import pathlib

import matplotlib
import numpy as np

from pypeit.scripts.chk_for_calibs import ChkForCalibs
from pypeit.tests.tstutils import data_output_path

from IPython import embed

matplotlib.use('agg')  # For Travis


def test_chk_calibs_not():
    os.chdir(data_output_path(''))
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/not_alfosc/grism4')
    droot += '/ALD'

    pargs = ChkForCalibs.parse_args([droot, '-s', 'not_alfosc'])
    answers, _ = ChkForCalibs.main(pargs)

    assert answers['pass'][0], 'One or more failures!'


def test_chk_calibs_deimos():
    os.chdir(data_output_path(''))
    # 830G
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_deimos/830G_M_8600/')
    pargs = ChkForCalibs.parse_args([droot, '-s', 'keck_deimos'])
    answers, _ = ChkForCalibs.main(pargs)
    assert np.all(answers['pass'][:-1]), 'One or more failures!'

    # 600ZD
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_deimos/600ZD_M_7500/')
    pargs = ChkForCalibs.parse_args([droot, '-s', 'keck_deimos'])
    answers, _ = ChkForCalibs.main(pargs)
    assert np.all(answers['pass'][:-1]), 'One or more failures!'

def test_chk_calibs_save_setups():
    os.chdir(data_output_path(''))
    # 830G
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_deimos/830G_M_8600/')
    pargs = ChkForCalibs.parse_args([droot, '-s', 'keck_deimos', '--save_setups'])
    answers, ps = ChkForCalibs.main(pargs)
    assert np.all(answers['pass'][:-1]), 'One or more failures!'

    # Check for the created setup files
    output_path = pathlib.Path().absolute() / 'setup_files'
    assert output_path.is_dir()
    pypeit_file = output_path / ps.pypeit_file
    assert pypeit_file.is_file()
    assert pypeit_file.with_suffix('.sorted').is_file()
    assert pypeit_file.with_suffix('.calib').is_file()

    # Delete all created files
    pypeit_file.unlink()
    pypeit_file.with_suffix('.sorted').unlink()
    pypeit_file.with_suffix('.calib').unlink()
    output_path.unlink()

