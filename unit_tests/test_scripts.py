"""
Module to run tests on scripts
"""
import os
import shutil
from pathlib import Path

from IPython import embed

import matplotlib
matplotlib.use('agg')  # For Travis

from pypeit import scripts
from pypeit import edgetrace

from pypeit.pypeitsetup import PypeItSetup


def test_trace_add_rm():
    # Define the output directories (HARDCODED!!)
    outdir = Path().resolve() / 'shane_kast_blue_A'
    calibdir = Path().resolve() / 'Calibrations'
    # Remove them if they already exist
    if outdir.exists():
        shutil.rmtree(outdir)
    if calibdir.exists():
        shutil.rmtree(calibdir)

    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/shane_kast_blue/600_4310_d55')

    # Run the setup
    ps = PypeItSetup.from_file_root(droot, 'shane_kast_blue')
    ps.run(setup_only=True)

    # Add lines to remove and add slits. This removes the one slit that
    # is found and adds another.
    ps.user_cfg += ['[calibrations]', '[[slitedges]]', 'rm_slits = 1:1028:170',
                    'add_slits = 1:1028:30:300', 'add_predict = straight']

    # Use PypeItMetaData to write the complete PypeIt file
    pypeit_file = ps.fitstbl.write_pypeit(output_path=Path().resolve(), cfg_lines=ps.user_cfg,
                                          configs=['all'])[0]

    # Run the tracing
    scripts.trace_edges.TraceEdges.main(
            scripts.trace_edges.TraceEdges.parse_args(['-f', pypeit_file]))

    # Define the edges master file (HARDCODED!!)
    trace_file = calibdir / 'Edges_A_0_DET01.fits.gz'

    # Check that the correct number of traces were found
    edges = edgetrace.EdgeTraceSet.from_file(trace_file)
    assert edges.ntrace == 2, 'Did not find the expected number of traces.'

    # Clean up
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


def test_show_arxiv():
    # Pargs
    pargs = scripts.show_arxiv.ShowArxiv.parse_args(['gemini_gmos_r831_ham.fits', '--det', '1',
                                                     '--test'])
    scripts.show_arxiv.ShowArxiv.main(pargs)


