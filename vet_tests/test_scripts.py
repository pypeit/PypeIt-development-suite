"""
Module to run tests on scripts
"""
import os
import numpy as np
import pytest

import matplotlib
from IPython import embed
matplotlib.use('agg')  # For Travis


from pypeit.scripts import parse_slits
from pypeit import scripts
from pypeit.tests.tstutils import data_path
from pypeit.display import display
from pypeit import wavecalib
from pypeit import coadd1d

from pypeit.pypmsgs import PypeItError


def test_parse_calibs(redux_out):
    pypeit_file = os.path.join(redux_out, 'keck_nires', 'NIRES', 'keck_nires_nires.pypeit')
    # Just list
    pargs = scripts.parse_calib_id.ParseCalibID.parse_args([pypeit_file])
    scripts.parse_calib_id.ParseCalibID.main(pargs)


def test_show_1dspec(redux_out):
    spec_file = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A', 'Science',
                             'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits')
    # Just list
    pargs = scripts.show_1dspec.Show1DSpec.parse_args([spec_file, '--list'])
    scripts.show_1dspec.Show1DSpec.main(pargs)


def test_show_2dspec(redux_out):
    droot = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A') 
    spec2d_file = os.path.join(droot, 'Science',
                             'spec2d_b27-J1217p3905_KASTb_20150520T045733.560.fits')
    # Ginga needs to be open in RC mode
    display.connect_to_ginga(raise_err=True, allow_new=True)
    # Save
    cdir = os.getcwd()
    os.chdir(droot)
    # List
    pargs = scripts.show_2dspec.Show2DSpec.parse_args([spec2d_file, '--list'])
    scripts.show_2dspec.Show2DSpec.main(pargs)
    # Show
    pargs = scripts.show_2dspec.Show2DSpec.parse_args([spec2d_file])
    scripts.show_2dspec.Show2DSpec.main(pargs)
    # Go back
    os.chdir(cdir)


def test_chk_edges(redux_out):
    mstrace_root = os.path.join(redux_out,
                                'keck_lris_red', 
                                'multi_400_8500_d560', 
                                'Masters', 
                                'MasterEdges_A_1_DET01.fits.gz')
    # Ginga needs to be open in RC mode
    display.connect_to_ginga(raise_err=True, allow_new=True)
    #
    pargs = scripts.chk_edges.ChkEdges.parse_args([mstrace_root])
    scripts.chk_edges.ChkEdges.main(pargs)


def test_view_fits_list(redux_out):
    """ Test the list option
    """
    spec_file = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A', 'Science',
                             'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits')
    #spec_file = data_path('spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    pargs = scripts.view_fits.ViewFits.parse_args(['shane_kast_blue', spec_file, '--list'])
    scripts.view_fits.ViewFits.main(pargs)


def test_view_fits_proc_fail(redux_out):
    """ Test that it fails when trying to proc an output pypeit image
    """
    droot = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A') 
    spec_file = os.path.join(droot, 'Science',
                             'spec2d_b27-J1217p3905_KASTb_20150520T045733.560.fits')
    #spec_file = data_path('spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    pargs = scripts.view_fits.ViewFits.parse_args(['shane_kast_blue', spec_file, '--proc'])
    with pytest.raises(PypeItError):
        scripts.view_fits.ViewFits.main(pargs)


def test_chk_flat(redux_out):
    droot = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A') 
    mstrace_root = os.path.join(droot,
                                'Masters',
                                'MasterFlat_A_1_DET01.fits')
    # Ginga needs to be open in RC mode
    display.connect_to_ginga(raise_err=True, allow_new=True)
    #
    pargs = scripts.chk_flats.ChkFlats.parse_args([mstrace_root])
    scripts.chk_flats.ChkFlats.main(pargs)


def test_chk_wavecalib(redux_out):
    droot = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A') 
    ms_root = os.path.join(droot,
                           'Masters',
                           'MasterWaveCalib_A_1_DET01.fits')
    #
    pargs = scripts.chk_wavecalib.ChkWaveCalib.parse_args([ms_root])
    scripts.chk_wavecalib.ChkWaveCalib.main(pargs)


def test_identify(redux_out):
    droot = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A') 
    arc_file = os.path.join(droot, 'Masters',
                             'MasterArc_A_1_DET01.fits')
    slits_file = os.path.join(droot, 'Masters',
                            'MasterSlits_A_1_DET01.fits.gz')
    # Just list
    pargs = scripts.identify.Identify.parse_args([arc_file, slits_file, '--test'])
    arcfitter = scripts.identify.Identify.main(pargs)

    # Load line list
    arcfitter.load_IDs(fname=data_path('waveid_tests.ascii'))
    assert arcfitter._detns.size == 31, 'Bad load'

    # Fit
    arcfitter._fitdict['polyorder'] = 3
    arcfitter.fitsol_fit()
    assert arcfitter._fitdict['fitc'].size == 4, 'Bad fit'

    # Auto
    arcfitter.auto_id()
    assert np.sum(arcfitter._lineflg < 3) > 10, 'Bad auto ID'
    arcfitter.fitsol_fit()

    # Write
    final_fit = arcfitter.get_results()

    waveCalib = wavecalib.WaveCalib(nslits=1, wv_fits=np.atleast_1d(arcfitter._fitdict['WaveFit']),
                              arc_spectra=np.atleast_2d(arcfitter.specdata).T,
                              spat_ids=np.atleast_1d(int(arcfitter._spatid)),
                              PYP_SPEC='shane_kast_blue',
                              )

    # If you touch the following line, you probably need to update the call in scripts/identify.py
    wvarxiv_fn = arcfitter.store_solution(final_fit, 1, rmstol=0.1, force_save=True, wvcalib=waveCalib)

    # Test we can read it
    tmp = wavecalib.WaveCalib.from_file('wvcalib.fits')

    # Clean up -- If these fail then the store solution failed
    os.remove('waveid.ascii')
    os.remove(wvarxiv_fn)
    os.remove('wvcalib.fits')


def test_compare_sky(redux_out):
    spec_file = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A', 'Science',
                             'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits')
    sky_file = 'sky_kastb_600.fits'

    # Running in `test` mode for boxcar extraction
    pargs = scripts.compare_sky.CompareSky.parse_args([spec_file, sky_file, '--test'])
    scripts.compare_sky.CompareSky.main(pargs)

    # Running in `test` mode for optimal extraction
    pargs = scripts.compare_sky.CompareSky.parse_args([spec_file, sky_file, '--test',
                                                       '--optimal'])
    scripts.compare_sky.CompareSky.main(pargs)


def test_collate_1d(tmp_path, monkeypatch, redux_out):
    kastb_dir = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A', 'Science')
    deimos_dir = os.path.join(redux_out, 
                              'keck_deimos','830G_M_8500', 
                              'Science')

    # Build up arguments for testing command line parsing
    args = ['--dry_run', '--ignore_flux', '--flux', '--outdir', '/outdir2', '--match', 'ra/dec', '--exclude_slit_bm', 'BOXSLIT', '--exclude_serendip', '--wv_rms_thresh', '0.2', '--refframe', 'heliocentric']
    spec1d_file = os.path.join(kastb_dir, 'spec1d_b27*fits')
    spec1d_args = ['--spec1d_files', spec1d_file]
    tol_args = ['--tolerance', '0.03d']
    alt_spec1d = os.path.join(deimos_dir, 'spec1d_DE.20100913.22358*fits')
    expanded_spec1d = os.path.join(kastb_dir,
                                   'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits')
    expanded_alt_spec1d = os.path.join(deimos_dir,
                                       'spec1d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits')
    spec1d_args = ['--spec1d_files', expanded_spec1d]

    # Create config files for testing config file parsing
    config_file_full = str(tmp_path / "test_collate1d_full.collate1d")

    with open(config_file_full, "w") as f:
        print("[coadd1d]", file=f)
        print("ex_value = BOX", file=f)
        print("[collate1d]", file=f)
        print("dry_run = False", file=f)
        print("flux = False", file=f)
        print("outdir = /outdir", file=f)
        print("ignore_flux = False", file=f)
        print("tolerance = 4.0", file=f)
        print("match_using = 'pixel'", file=f)
        print("exclude_slit_trace_bm = BADSKYSUB,BADEXTRACT", file=f)
        print("exclude_serendip = False", file=f)
        print("wv_rms_thresh = 0.1", file=f)
        print("refframe = 'observed'", file=f)
        print("spec1d read", file=f)
        print("filename", file=f)
        print(alt_spec1d, file=f)
        print("spec1d end", file=f)

    config_file_spec1d = str(tmp_path / "test_collate1d_spec1d_only.collate1d")
    with open(config_file_spec1d, "w") as f:
        print("[collate1d]", file=f)
        print("spec1d read", file=f)
        print("filename", file=f)
        print(spec1d_file, file=f)
        print("spec1d end", file=f)

    config_file_coadd1d = str(tmp_path / "test_collate1d_spec1d_only.coadd1d")
    with open(config_file_coadd1d, "w") as f:
        print("[coadd1d]", file=f)
        print("ex_value = BOX", file=f)
        # Need a data block
        print(" ", file=f)
        print("coadd1d read", file=f)
        print("filename | obj_id", file=f)
        print(alt_spec1d + ' | DUMMY', file=f)
        print("coadd1d end", file=f)

    # Args only, nospec1d files should exit with an errror
    with pytest.raises(SystemExit):
        parsed_args = scripts.collate_1d.Collate1D.parse_args(args + tol_args)
        params, spectrograph, expanded_spec1d_files \
                = scripts.collate_1d.build_parameters(parsed_args)

    # Everything passed via command line
    parsed_args = scripts.collate_1d.Collate1D.parse_args(args + tol_args + spec1d_args)
    params, spectrograph, expanded_spec1d_files = scripts.collate_1d.build_parameters(parsed_args)
    assert params['collate1d']['dry_run'] is True
    assert params['collate1d']['outdir'] == '/outdir2'
    assert params['collate1d']['match_using'] == 'ra/dec'
    assert params['collate1d']['tolerance'] == '0.03d'
    assert params['collate1d']['exclude_slit_trace_bm'] == ['BOXSLIT']
    assert params['collate1d']['exclude_serendip'] is True
    assert params['collate1d']['wv_rms_thresh'] == 0.2
    assert params['coadd1d']['ex_value'] == 'OPT'
    assert params['collate1d']['refframe'] == 'heliocentric'
    assert spectrograph.name == 'shane_kast_blue'
    assert len(expanded_spec1d_files) == 1 and expanded_spec1d_files[0] == expanded_spec1d

    # Full config file, should work
    parsed_args = scripts.collate_1d.Collate1D.parse_args([config_file_full])
    params, spectrograph, expanded_spec1d_files = scripts.collate_1d.build_parameters(parsed_args)
    assert params['collate1d']['dry_run'] is False
    assert params['collate1d']['outdir'] == '/outdir'
    assert params['collate1d']['ignore_flux'] == False
    assert params['collate1d']['flux'] == False
    assert params['collate1d']['tolerance'] == 4.0
    assert params['collate1d']['match_using'] == 'pixel'
    assert params['collate1d']['exclude_slit_trace_bm'] == ['BADSKYSUB','BADEXTRACT']
    assert params['collate1d']['exclude_serendip'] is False
    assert params['collate1d']['wv_rms_thresh'] == 0.1
    assert params['coadd1d']['ex_value'] == 'BOX'
    assert params['collate1d']['refframe'] == 'observed'
    assert spectrograph.name == 'keck_deimos'
    assert len(expanded_spec1d_files) == 1 and expanded_spec1d_files[0] == expanded_alt_spec1d

    # Test that a full command line overrides a config file
    parsed_args = scripts.collate_1d.Collate1D.parse_args(args + spec1d_args + tol_args
                                                          + [config_file_full])
    params, spectrograph, expanded_spec1d_files = scripts.collate_1d.build_parameters(parsed_args)
    assert params['collate1d']['dry_run'] is True
    assert params['collate1d']['outdir'] == '/outdir2'
    assert params['collate1d']['ignore_flux'] == True
    assert params['collate1d']['flux'] == True
    assert params['collate1d']['tolerance'] == '0.03d'
    assert params['collate1d']['match_using'] == 'ra/dec'
    assert params['collate1d']['exclude_slit_trace_bm'] == ['BOXSLIT']
    assert params['collate1d']['exclude_serendip'] is True
    assert params['collate1d']['wv_rms_thresh'] == 0.2
    assert params['collate1d']['refframe'] == 'heliocentric'
    assert spectrograph.name == 'shane_kast_blue'
    assert len(expanded_spec1d_files) == 1 and expanded_spec1d_files[0] == expanded_spec1d

    # Test that a config file with spec1d files. Test that default tolerance and match_using is used
    # Also test using an external coadd1d file with the same name
    parsed_args = scripts.collate_1d.Collate1D.parse_args([config_file_spec1d])
    params, spectrograph, expanded_spec1d_files = scripts.collate_1d.build_parameters(parsed_args)
    assert params['collate1d']['tolerance'] == 1.0
    assert params['collate1d']['match_using'] == 'ra/dec'
    assert params['coadd1d']['ex_value'] == 'BOX'
    assert spectrograph.name == 'shane_kast_blue'
    assert len(expanded_spec1d_files) == 1 and expanded_spec1d_files[0] == expanded_spec1d

    # Mocks for testing main
    class MockCoadd:
        def run(*args, **kwargs):
            pass

        def save(self, file):
            if os.path.basename(file) == "J232856.20-030325.90_DEIMOS_20100913.fits":
                raise ValueError("test exception")

    def mock_get_instance(*args, **kwargs):
        return MockCoadd()

    def mock_get_subdir(*args, **kwargs):
        return "subdir"

    with monkeypatch.context() as m:
        monkeypatch.setattr(coadd1d.CoAdd1D, "get_instance", mock_get_instance)

        os.chdir(tmp_path)
        par_file = str(tmp_path / 'collate1d.par')
        
        # Test:
        # * main
        # * creation of collate1d.par
        # * parsing of pixel tolerance
        # * detection of spec2d files and excluding by slit bitmask
        # * copying of spec1d file when doing refframe correction and spec1d_output is set.
        archive_dir = tmp_path / 'archive'

        parsed_args = scripts.collate_1d.Collate1D.parse_args(['--par_outfile', par_file, '--match',
                                                               'pixel', '--tolerance', '3',
                                                               '--spec1d_files', expanded_spec1d,
                                                               '--spec1d_outdir', str(tmp_path),
                                                               '--refframe', 'heliocentric',
                                                               '--exclude_slit_bm', 'BADSKYSUB,BADEXTRACT'])
        assert scripts.collate_1d.Collate1D.main(parsed_args) == 0
        assert os.path.exists(par_file)
        assert os.path.exists(os.path.join(str(tmp_path), os.path.basename(expanded_spec1d)))

        # Remove par_file to avoid a warning
        os.unlink(par_file)
        
        # Test:
        # * default units of arcsec for tolerance when match is ra/dec
        # * that a spec2d file isn't needed if exclude_slit_flags is empty.
        # * test exception handling when one file fails
        # The 240 arsec tolerance is to ensure there's only two outputs, one of which the mock 
        # coadd object will fail
        parsed_args = scripts.collate_1d.Collate1D.parse_args(['--par_outfile', par_file,
                                                               '--match', 'ra/dec', '--tolerance',
                                                               '240', '--spec1d_files',
                                                               expanded_alt_spec1d])                                                            
        assert scripts.collate_1d.Collate1D.main(parsed_args) == 0

        # Remove par_file to avoid a warning
        os.unlink(par_file)

        # Test parsing of units in ra/dec tolerance
        parsed_args = scripts.collate_1d.Collate1D.parse_args(['--par_outfile', par_file,
                                                               '--match', 'ra/dec', '--tolerance', '3d',
                                                               '--spec1d_files', expanded_alt_spec1d])
        assert scripts.collate_1d.Collate1D.main(parsed_args) == 0
        
# pypeit_parse_calib_id is tested in test_runpypeit

def test_parse_slits(redux_out):
    kastb_dir = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A')
    slits_file = os.path.join(kastb_dir, 'Masters',
                              'MasterSlits_A_1_DET01.fits.gz')
    spec2d_file = os.path.join(kastb_dir, 'Science',
                              'spec2d_b27-J1217p3905_KASTb_20150520T045733.560.fits')

    # Slits
    pargs = parse_slits.ParseSlits.parse_args([slits_file])
    parse_slits.ParseSlits.main(pargs)

    # Spec2d
    pargs = parse_slits.ParseSlits.parse_args([spec2d_file])
    parse_slits.ParseSlits.main(pargs)
    


# TODO: Include tests for coadd2d, sensfunc


