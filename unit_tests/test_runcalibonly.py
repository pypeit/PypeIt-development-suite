"""
Module to run calibrations only test with run_pypeit
"""
import os
import shutil
from IPython import embed

from pypeit.scripts.setup import Setup
from pypeit.scripts.run_pypeit import RunPypeIt


def test_run_pypeit_calib_only():
    # Get the directories
    rawdir = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'shane_kast_blue', '600_4310_d55')
    assert os.path.isdir(rawdir), 'Incorrect raw directory'

    calib_key = 'A_0_DET01'

    # File list
    all_files = {
        'arcs': ['b1.fits.gz'],
        'flats': ['b11.fits.gz', 'b12.fits.gz', 'b13.fits.gz'],
        'bias': ['b21.fits.gz', 'b22.fits.gz', 'b23.fits.gz'],
    }
    all_calibs = [f'Arc_{calib_key}.fits',
                   f'Tiltimg_{calib_key}.fits',
                   f'Bias_{calib_key}.fits',
                   f'Tilts_{calib_key}.fits',
                   f'Edges_{calib_key}.fits.gz',
                   f'Flat_{calib_key}.fits',
                   f'WaveCalib_{calib_key}.fits']

    # Just get a few files
    for ss, sub_files, calibs in zip(range(3),
            [['arcs', 'flats', 'bias'],
             ['arcs', 'bias'],
             ['flats', 'bias']],
            [all_calibs,
             [f'Arc_{calib_key}.fits', f'Tiltimg_{calib_key}.fits'],
             [f'Edges_{calib_key}.fits.gz']]):
        # Grab the subset
        files = []
        for sub_file in sub_files:
            files += all_files[sub_file]
        #
        testrawdir = os.path.join(os.getenv('PYPEIT_DEV'), 'TEST')
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
        for calib_file in calibs:
            assert os.path.isfile(os.path.join(configdir, 'Calibrations', calib_file)), \
                         f'Calibration file {calib_file} missing!'

        # Clean-up
        shutil.rmtree(outdir)
        shutil.rmtree(testrawdir)

