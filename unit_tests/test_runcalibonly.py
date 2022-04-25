"""
Module to run calibrations only test with run_pypeit
"""
import os
import shutil
from IPython.terminal.embed import embed

from configobj import ConfigObj

import pytest

from pypeit.scripts.parse_calib_id import ParseCalibID
from pypeit.scripts.setup import Setup
from pypeit.scripts.run_pypeit import RunPypeIt


def test_run_pypeit_calib_only():
    # Get the directories
    rawdir = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'shane_kast_blue', '600_4310_d55')
    assert os.path.isdir(rawdir), 'Incorrect raw directory'

    master_key = 'A_1_DET01'

    # File list
    all_files = {
        'arcs': ['b1.fits.gz'],
        'flats': ['b11.fits.gz', 'b12.fits.gz', 'b13.fits.gz'],
        'bias': ['b21.fits.gz', 'b22.fits.gz', 'b23.fits.gz'],
    }
    all_masters = [f'MasterArc_{master_key}.fits',
                   f'MasterTiltimg_{master_key}.fits',
                   f'MasterBias_{master_key}.fits',
                   f'MasterTilts_{master_key}.fits',
                   f'MasterEdges_{master_key}.fits.gz',
                   f'MasterFlat_{master_key}.fits',
                   f'MasterWaveCalib_{master_key}.fits']

    # Just get a few files
    for ss, sub_files, masters in zip(range(3),
            [['arcs', 'flats', 'bias'],
             ['arcs', 'bias'],
             ['flats', 'bias']],
            [all_masters,
             [f'MasterArc_{master_key}.fits', f'MasterTiltimg_{master_key}.fits'],
             [f'MasterEdges_{master_key}.fits.gz']]):
        # Grab the subset
        files = []
        for sub_file in sub_files:
            files += all_files[sub_file]
        #
        testrawdir = os.path.join(rawdir, 'TEST')
        if os.path.isdir(testrawdir):
            shutil.rmtree(testrawdir)
        os.makedirs(testrawdir)
        for f in files:
            shutil.copy(os.path.join(rawdir, f), os.path.join(testrawdir, f))

        outdir = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT_TEST')

        # For previously failed tests
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)

        # Run the setup
        sargs = Setup.parse_args(['-r', testrawdir, '-s', 'shane_kast_blue', '-c all', '-o',
                                  '--output_path', outdir])
        Setup.main(sargs)

        # Change to the configuration directory and set the pypeit file
        configdir = os.path.join(outdir, 'shane_kast_blue_A')
        pyp_file = os.path.join(configdir, 'shane_kast_blue_A.pypeit')
        assert os.path.isfile(pyp_file), 'PypeIt file not written.'

        # Perform the calib-only reduction
        pargs = RunPypeIt.parse_args([pyp_file, '-c', '-r', configdir])
        RunPypeIt.main(pargs)

        # Test!
        for master_file in masters:
            assert os.path.isfile(os.path.join(configdir, 'Masters', master_file)
                                  ), 'Master File {:s} missing!'.format(master_file)

        # Now test parse_calib_id
        if ss == 0:
            pargs2 = ParseCalibID.parse_args([pyp_file])
            calib_dict = ParseCalibID.main(pargs2)
            assert isinstance(calib_dict, dict)
            assert len(calib_dict) > 0
            assert calib_dict['1'][master_key]['arc']['raw_files'][0] == 'b1.fits.gz'

        # Clean-up
        shutil.rmtree(outdir)
        shutil.rmtree(testrawdir)

