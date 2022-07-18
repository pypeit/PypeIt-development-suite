"""
Module to run tests on scripts
"""
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

import matplotlib
from IPython import embed
matplotlib.use('agg')  # For Travis

from astropy.io import fits

from pypeit import scripts
from pypeit.tests.tstutils import data_path
from pypeit import edgetrace
from pypeit import io
from pypeit import fluxcalibrate
from pypeit import onespec

from pypeit.pypeitsetup import PypeItSetup
from pypeit.pypmsgs import PypeItError


def test_quicklook():
    # The following needs the LRISb calibration files to be
    # found in a CALIBS/ folder in the DEV Suite.  You can get
    # that folder here:
    # https://drive.google.com/drive/folders/1NSg5Rmr8hD_1-fOchQc3WXjt59D6f9od?usp=sharing
    calib_dir = os.path.join(os.environ['PYPEIT_DEV'], 'CALIBS')
    if not os.path.isdir(calib_dir):
        raise IOError("You need to get the CALIBS folder as described above!!")

    # Define the output directories (HARDCODED!!)
    cdir = os.getcwd()
    os.chdir(data_path(''))
    outdir = data_path('keck_lris_blue_A')
    # Remove them if they already exist
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    # Raw path
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_blue',
                         'long_600_4000_d560')
    pixflat = os.path.join(calib_dir, 'PYPEIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz')
    scripts.ql_mos.QLMOS.main(
            scripts.ql_mos.QLMOS.parse_args(['keck_lris_blue', droot, 'b150910_2033.fits.gz',
                                             'b150910_2051.fits.gz', 'b150910_2070.fits.gz',
                                             '--det=2', f'--user_pixflat={pixflat}']))
    
    # Cleanup
    os.chdir(cdir)
    shutil.rmtree(outdir)

# THIS TEST HAS MAJOR PATH ISSUES
#def test_trace_edges():
#    # Define the output directories 
#    tmpdir = data_path('TMP')
#    if os.path.isdir(tmpdir):
#        shutil.rmtree(tmpdir)
#    os.mkdir(tmpdir) 
#    #
#    setupdir = os.path.join(tmpdir, 'setup_files')
#    outdir = os.path.join(tmpdir, 'shane_kast_blue_A')
#    masterdir = os.path.join(tmpdir, 'shane_kast_blue_A', 'Masters')
#    # Remove them if they already exist
#    if os.path.isdir(setupdir):
#        shutil.rmtree(setupdir)
#    if os.path.isdir(outdir):
#        shutil.rmtree(outdir)
#
#    # Perform the setup
#    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/shane_kast_blue/600_4310_d55')
#    droot += '/'
#    scripts.setup.Setup.main(scripts.setup.Setup.parse_args(['-r', droot, '-s', 'shane_kast_blue',
#                                                             '-c', 'all']))
#
#    # Generate the Masters folder
#    pytest.set_trace()
#    os.mkdir(masterdir)
#
#    # Define the pypeit file (HARDCODED!!)
#    pypeit_file = os.path.join(outdir, 'shane_kast_blue_A.pypeit')
#
#    # Run the tracing
#    scripts.trace_edges.TraceEdges.main(
#            scripts.trace_edges.TraceEdges.parse_args(['-f', pypeit_file]))
#
#    # Define the edges master file (HARDCODED!!)
#    trace_file = os.path.join(outdir, 'Masters', 'MasterEdges_A_1_DET01.fits.gz')
#
#    # Check that the correct number of traces were found
#    edges = edgetrace.EdgeTraceSet.from_file(trace_file)
#    assert edges.ntrace == 2, 'Did not find the expected number of traces.'
#
#    # Clean up
#    shutil.rmtree(setupdir)
#    shutil.rmtree(outdir)


def test_trace_add_rm():
    # Define the output directories (HARDCODED!!)
    setupdir = os.path.join(os.getcwd(), 'setup_files')
    outdir = os.path.join(os.getcwd(), 'shane_kast_blue_A')
    masterdir = os.path.join(os.getcwd(), 'shane_kast_blue_A', 'Masters')
    # Remove them if they already exist
    if os.path.isdir(setupdir):
        shutil.rmtree(setupdir)
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/shane_kast_blue/600_4310_d55')

    # Run the setup
    ps = PypeItSetup.from_file_root(droot, 'shane_kast_blue', output_path=setupdir)
    ps.run(setup_only=True, sort_dir=setupdir)

    # Add lines to remove and add slits. This removes the one slit that
    # is found and adds another.
    ps.user_cfg += ['[calibrations]', '[[slitedges]]', 'rm_slits = 1:1028:170',
                    'add_slits = 1:1028:30:300', 'add_predict = straight']

    # Use PypeItMetaData to write the complete PypeIt file
    pypeit_file = ps.fitstbl.write_pypeit(output_path=os.getcwd(), cfg_lines=ps.user_cfg,
                                          configs=['all'])[0]

    # Run the tracing
    scripts.trace_edges.TraceEdges.main(
            scripts.trace_edges.TraceEdges.parse_args(['-f', pypeit_file]))

    # Define the edges master file (HARDCODED!!)
    trace_file = os.path.join(outdir, 'Masters', 'MasterEdges_A_1_DET01.fits.gz')

    # Check that the correct number of traces were found
    edges = edgetrace.EdgeTraceSet.from_file(trace_file)
    assert edges.ntrace == 2, 'Did not find the expected number of traces.'

    # Clean up
    shutil.rmtree(setupdir)
    shutil.rmtree(outdir)



def test_view_fits_proc():
    """ Test that it works on a raw image
    """
    spec_file = Path(os.getenv('PYPEIT_DEV')).resolve() / 'RAW_DATA' / 'shane_kast_blue' \
                    / '830_3460_d46' / 'b100.fits.gz'
    pargs = scripts.view_fits.ViewFits.parse_args(['shane_kast_blue', str(spec_file), '--proc'])
    scripts.view_fits.ViewFits.main(pargs)


def test_view_fits_mosaic():
    """ Test the list option
    """
    file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'gemini_gmos', 'GN_HAM_R400_885',
                        'N20190205S0035.fits')
    pargs = scripts.view_fits.ViewFits.parse_args(['gemini_gmos_north_ham', file,
                                                   '--det', 'mosaic',
                                                   '--proc'])
    scripts.view_fits.ViewFits.main(pargs)



# TODO -- REMOVE if this test is being run in the DEV SUITE afterburn.
#   Am pretty sure it is..
def test_coadd1d_1(monkeypatch):
    """
    Test basic coadd using shane_kast_blue
    """
    dp = data_path('')
    # Change to the parent directory of the data path, so we can test that
    # coadding without a coadd output file specified places the output next
    # to the spec1ds. Using monkeypatch means the current working directory
    # will be restored after the test.
    monkeypatch.chdir(Path(dp).parent)

    # NOTE: flux_value is False
    parfile = 'files/coadd1d.par'
    if os.path.isfile(parfile):
        os.remove(parfile)
    coadd_ofile = data_path('coadd1d_J1217p3905_KASTb_20150520_20150520.fits')
    if os.path.isfile(coadd_ofile):
        os.remove(coadd_ofile)

    coadd_ifile = data_path('shane_kast_blue.coadd1d')
    scripts.coadd_1dspec.CoAdd1DSpec.main(
            scripts.coadd_1dspec.CoAdd1DSpec.parse_args([coadd_ifile, "--par_outfile", parfile]))
    hdu = io.fits_open(coadd_ofile)
    assert hdu[1].header['EXT_MODE'] == 'OPT'
    assert hdu[1].header['FLUXED'] is False
    # Test that the output file is kosher and contains the right quantities
    spec = onespec.OneSpec.from_file(coadd_ofile)
    assert spec.wave.shape == spec.wave_grid_mid.shape

    # Clean up
    hdu.close()
    os.remove(parfile)
    os.remove(coadd_ofile)



def test_obslog():
    # Define the output directories (HARDCODED!!)
    setupdir = os.path.join(os.getcwd(), 'setup_files')
    obslogfile = 'shane_kast_blue.obslog'
    # Remove the directory if it already exists
    if os.path.isdir(setupdir):
        shutil.rmtree(setupdir)

    # Perform the setup
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/shane_kast_blue/600_4310_d55')
    scripts.obslog.ObsLog.main(scripts.obslog.ObsLog.parse_args(['shane_kast_blue', '-r', droot,
                                                                 '-f', obslogfile, '-d', setupdir]))

    # Clean up
    shutil.rmtree(setupdir)



#  THIS TEST IS RUN IN THE DEV SUITE AS PART OF AFTERBURN
#   AND THE VERSION BELOW IS A MISH-MASH OF PYPEIT TEST + COOKED
#  TODO -- Could add the checking part as a VET test
#def test_coadd1d_2():
#    """
#    Test combining Echelle
#    """
#    # NOTE: flux_value is False
#    parfile = 'coadd1d.par'
#    if os.path.isfile(parfile):
#        os.remove(parfile)
#    coadd_ofile = data_path('pisco_coadd.fits')
#    if os.path.isfile(coadd_ofile):
#        os.remove(coadd_ofile)
#
#    coadd_ifile = data_path('gemini_gnirs_32_sb_sxd.coadd1d')
#    pytest.set_trace()
#    scripts.coadd_1dspec.CoAdd1DSpec.main(
#            scripts.coadd_1dspec.CoAdd1DSpec.parse_args([coadd_ifile, '--test_spec_path',
#                                                         data_path('')]))
#
#    hdu = io.fits_open(coadd_ofile)
#    assert hdu[1].header['EXT_MODE'] == 'OPT'
#    assert hdu[1].header['FLUXED'] is False
#
#    # Test that the output file is kosher and contains the right quantities
#    spec = onespec.OneSpec.from_file(coadd_ofile)
#    assert spec.wave.shape == spec.wave_grid_mid.shape
#
#    # Clean up
#    hdu.close()
#    os.remove(parfile)
#    os.remove(coadd_ofile)
