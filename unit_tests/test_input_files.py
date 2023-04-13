""" Tests Reading of PypeIt Input files """
from pathlib import Path
import os
import shutil
import glob

from IPython import embed

from pypeit import inputfiles
from pypeit.scripts.setup import Setup


def test_read_fluxing_files():
    """ Test reading fluxing files """
    # Grab em
    fluxing_files = glob.glob(os.path.join(
        os.getenv('PYPEIT_DEV'), 'fluxing_files', '*.flux'))
    # Loop
    for ifile in fluxing_files:
        fluxFile = inputfiles.FluxFile.from_file(ifile)
        assert len(fluxFile.data['filename']) > 0

def test_read_coadd1d_files():
    """ Test reading coadd1d files """
    # Grab em
    coadd1d_files = glob.glob(os.path.join(
        os.getenv('PYPEIT_DEV'), 'coadd1d_files', '*.coadd1d'))
    # Loop
    for ifile in coadd1d_files:
        coadd1dFile = inputfiles.Coadd1DFile.from_file(ifile)
        assert len(coadd1dFile.data['filename']) > 0

def test_read_coadd2d_files():
    """ Test reading coadd2d files """
    # Grab em
    coadd2d_files = glob.glob(os.path.join(
        os.getenv('PYPEIT_DEV'), 'coadd2d_files', '*.coadd2d'))
    # Loop
    for ifile in coadd2d_files:
        coadd2dFile = inputfiles.Coadd2DFile.from_file(ifile)
        assert len(coadd2dFile.data['filename']) > 0

def test_read_flexure_files():
    """ Test reading flexure files """
    # Grab em
    flexure_files = glob.glob(os.path.join(
        os.getenv('PYPEIT_DEV'), 'flexure_files', '*.flex'))
    # Loop
    for ifile in flexure_files:
        flexureFile = inputfiles.FlexureFile.from_file(ifile)
        assert len(flexureFile.data['filename']) > 0

def test_read_pypeit_files():
    """ Test reading PypeIt files """
    # Grab em
    pypeit_files = glob.glob(os.path.join(
        os.getenv('PYPEIT_DEV'), 'pypeit_files', '*.pypeit'))
    # Loop
    for ifile in pypeit_files:
        pypeitFile = inputfiles.PypeItFile.from_file(ifile)


def test_get_spectrograph():
    pypeit_file = Path(os.getenv('PYPEIT_DEV')) / 'pypeit_files' \
                    / 'shane_kast_blue_452_3306_d57.pypeit'
    assert pypeit_file.exists(), 'Missing test pypeit file'

    pypeitFile = inputfiles.PypeItFile.from_file(pypeit_file)
    spec = pypeitFile.get_spectrograph()
    assert spec.name == 'shane_kast_blue', 'Wrong spectrograph'


def test_get_pypeitpar():

    spec = 'shane_kast_blue'
    setup = '600_4310_d55'

    # Define the path with the raw data
    data_root = Path(os.environ['PYPEIT_DEV']).resolve() / 'RAW_DATA' / spec / setup
    assert data_root.exists(), 'TEST ERROR: Raw data path does not exist'

    # Define the output directory and remove it if it already exist
    setup_path = Path().resolve() / f'{spec}_A'
    if setup_path.exists():
        shutil.rmtree(setup_path)

    args = ['-r', str(data_root), '-s', spec, '-c', 'all']
    pargs = Setup.parse_args(args)
    Setup.main(pargs)

    assert setup_path.exists(), 'No setup_files directory created'
    pypeit_file = setup_path / f'{spec}_A.pypeit'
    assert pypeit_file.exists(), 'PypeIt file not written'

    pypeitFile = inputfiles.PypeItFile.from_file(pypeit_file)
    spectrograph, par, file = pypeitFile.get_pypeitpar()

    assert spectrograph.name == spec, 'Wrong spectrograph'
    assert par['rdx']['spectrograph'] == spec, 'Wrong spectrgraph'
    assert par['rdx']['scidir'] == 'Science', 'Wrong default science directory'
    assert Path(file).name == 'b27.fits.gz', 'File name changed'

    # Clean-up
    shutil.rmtree(setup_path)

