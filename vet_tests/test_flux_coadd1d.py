import os
from pathlib import Path
from configobj import ConfigObj

import numpy as np

from pypeit.scripts import flux_calib, coadd_1dspec, flux_setup
from pypeit import specobjs, onespec
from pypeit.inputfiles import FluxFile, Coadd1DFile, TelluricFile 

def flux_and_validate(flux_file, output_files):

    # Run pypeit_flux_calib on a fluxing file and validate the output files
    parser = flux_calib.FluxCalib.get_parser()
    args = parser.parse_args([str(flux_file)])
    assert flux_calib.FluxCalib.main(args) == 0

    for file in output_files:
        sobjs = specobjs.SpecObjs.from_fitsfile(str(file))
        sobj = sobjs[0]
        assert 'OPT_FLAM' in sobj.keys(), f'OPT_FLAM missing from {file}'
        assert sobj.OPT_FLAM is not None, f'OPT_FLAM from {file} is None'
        assert len(sobj.OPT_FLAM) > 0, f'OPT_FLAM from {file} is empty'
        assert np.all(np.isfinite(sobj.OPT_FLAM)), f'OPT_FLAM from {file} has infinite or NaN values'



def coadd(coadd_file, output_file):
    # 1d coadd the vlt_xshooter NIR and VIS data

    if not isinstance(output_file, Path):
        output_file = Path(output_file)

    if output_file.exists():
        # Remove results from previous run
        os.unlink(output_file)

    parser = coadd_1dspec.CoAdd1DSpec.get_parser()
    args = parser.parse_args([str(coadd_file)])
    coadd_1dspec.CoAdd1DSpec.main(args)

    assert output_file.exists(), f"coadd1d output {output_file} does not exist"

    coadd_spec = onespec.OneSpec.from_file(str(output_file))
    mask = coadd_spec.mask
    assert len(coadd_spec.wave[mask]) > 0, f"coadd1d output {output_file} has empty wavelength array."
    assert np.all(np.isfinite(coadd_spec.wave[mask])), f"coadd1d output {output_file} has infinite or NaN values in it's wavelength array."

    assert len(coadd_spec.flux[mask]) > 0, f"coadd1d output {output_file} has empty flux array."
    assert np.all(np.isfinite(coadd_spec.flux[mask])), f"coadd1d output {output_file} has infinite or NaN values in it's flux array."

    assert coadd_spec.fluxed is True, f"coadd1d output {output_file} was not fluxed."


def run_flux_setup(input_paths, output_name):
        
    parser = flux_setup.FluxSetup.get_parser()
    args = parser.parse_args(['--name', output_name] + [str(path) for path in input_paths])
    flux_setup.FluxSetup.main(args)

    flux_filename = f"{output_name}.flux"
    coadd1d_filename = f"{output_name}.coadd1d"
    telluric_filename = f"{output_name}.tell"

    # Validate the flux_setup output files  can be read
    flux_file = FluxFile.from_file(flux_filename)
    coadd1d_file =  Coadd1DFile.from_file(coadd1d_filename)
    telluric_file = TelluricFile.from_file(telluric_filename)

    return flux_file, coadd1d_file, telluric_file

def test_flux_setup_vlt_xshooter(redux_out, monkeypatch):

    redux_out_path = Path(redux_out)
    with monkeypatch.context() as m:
        monkeypatch.chdir(redux_out_path / 'vlt_xshooter')

        # Make sure spec1ds were created correctly
        spec1d_vis_patterns = ['spec1d_*2016-08-02T08:45:46.510*.fits','spec1d_*2016-08-02T09:16:51.980*.fits' ] 
        spec1d_nir_patterns = ['spec1d_*2016-08-02T09:00:56.860*.fits','spec1d_*2016-08-02T09:16:54.772*.fits'] 

        vis_sci_path = Path('VIS_2x1', 'Science')
        nir_sci_path = Path('NIR', 'Science')

        
        spec1d_files = [] # Gather the actual spec1d filenames for later
        for pattern in spec1d_vis_patterns:
             files = list(vis_sci_path.glob(pattern))
             assert len(files) == 1, f"Missing spec1d file matching {vis_sci_path / pattern}"
             spec1d_files += files

        for pattern in spec1d_nir_patterns:
             files = list(nir_sci_path.glob(pattern))
             assert len(files) == 1,  f"Missing spec1d file matching {nir_sci_path / pattern}"
             spec1d_files += files
        
        # Validate sensfiles were created correctly
        sensfunc_paths = [Path('VIS_1x1_Feige110'),
                          Path('NIR_Feige110'),
                          ]
                       
        for path in sensfunc_paths:
             files = list(path.glob("sens_*.fits"))
             assert len(files) == 1, f"No sens file found in {path}"

        # Run flux setup
        flux_setup_output = "vlt_xshooter_vis_nir"
        flux_file, coadd1d_file, telluric_file = run_flux_setup([vis_sci_path, nir_sci_path] + sensfunc_paths, flux_setup_output)

        # We don't do much with the telluric file
        assert telluric_file is not None

        # Now update the fluxing and coadd1d files
        updated_flux_filename = f'updated_{flux_setup_output}.flux'
        updated_coadd1d_filename = f"updated_{flux_setup_output}.coadd1d"

        
        # Flux setup, remove undesired spec1ds
        spec1d_names = [path.name for path in spec1d_files]
        files_to_delete = np.nonzero([file not in spec1d_names for file in flux_file.data['filename']])
        flux_file.data.remove_rows(files_to_delete)
        flux_file.write(updated_flux_filename)

        # Coadd 1D setup, filter out undesired spec1ds
        files_to_delete = np.nonzero([file not in spec1d_names for file in coadd1d_file.data['filename']])
        coadd1d_file.data.remove_rows(files_to_delete)

        # Set a specific output file, and be forgiving of PypeIt version
        coadd_output_file = "J0100p2802_XShooter_VIS_NIR_coadd.fits"
        coadd1d_file.config = ConfigObj({'rdx': {'chk_version': False},
                                         'coadd1d': {'coaddfile': coadd_output_file}})
        coadd1d_file.write(updated_coadd1d_filename)

        # Now test the fluxing
        flux_and_validate(updated_flux_filename, spec1d_files)

        # Now test the coadding
        coadd(updated_coadd1d_filename, coadd_output_file)

def test_feige110_ltt3218(redux_out, monkeypatch):

    redux_out_path = Path(redux_out)
    with monkeypatch.context() as m:
        monkeypatch.chdir(redux_out_path / 'vlt_xshooter')

        # Make sure spec1ds were created correctly
        feige110_spec1d_patterns = ['NIR_Feige110/Science/spec1d_*01.793*.fits',
                                    'NIR_Feige110/Science/spec1d_*01.304*.fits',
                                    'VIS_1x1_Feige110/Science/spec1d_*58.946*.fits',
                                    ]

        ltt3218_spec1d_patterns =  ['NIR_LTT3218/Science/spec1d_*13.970*.fits',
                                    'NIR_LTT3218/Science/spec1d_*04.807*.fits',
                                    'VIS_1x1_LTT3218/Science/spec1d_*11.182*.fits']

        
        feige_spec1d_files = [] # Gather the actual spec1d filenames for later
        for path_pattern in feige110_spec1d_patterns:
             p = Path(path_pattern)
             pattern = p.name
             files = list(p.parent.glob(pattern))
             assert len(files) == 1, f"Missing spec1d file matching {path_pattern}"
             feige_spec1d_files += files

        ltt3218_spec1d_files = [] # Gather the actual spec1d filenames for later
        for path_pattern in ltt3218_spec1d_patterns:
             p = Path(path_pattern)
             pattern = p.name
             files = list(p.parent.glob(pattern))
             assert len(files) == 1, f"Missing spec1d file matching {path_pattern}"
             ltt3218_spec1d_files += files
        
        # Validate sensfiles were created correctly
        feige110_setup_paths = [Path('VIS_1x1_Feige110'),
                                Path('NIR_Feige110')]
        ltt3218_setup_paths =  [Path('VIS_1x1_LTT3218'),
                                Path('NIR_LTT3218'),
                               ]
                       
        for path in feige110_setup_paths + ltt3218_setup_paths:
             files = list(path.glob("sens_*.fits"))
             assert len(files) == 1, f"No sens file found in {path}"


        # Generate fluxing files to flux each std star with the other

        # Run flux setup for fluxing feige110 with ltt3218
        flux_setup_feige110_with_ltt3218_output = "vlt_xshooter_feige110_with_ltt3218"
        flux_feige110_with_ltt3218_file, feige110_coadd1d_file, telluric_file = run_flux_setup([p / 'Science' for p in feige110_setup_paths] + ltt3218_setup_paths, flux_setup_feige110_with_ltt3218_output)

        # Run flux setup for fluxing ltt3218 with feige110
        flux_setup_ltt3218_with_feige110_output = "vlt_xshooter_ltt3218_with_feige110"
        flux_ltt3218_with_feige110_file, ltt3218_coadd1d_file, telluric_file = run_flux_setup([p / 'Science' for p in ltt3218_setup_paths] + feige110_setup_paths, flux_setup_ltt3218_with_feige110_output)

        # Now test the fluxing
        flux_and_validate(flux_setup_feige110_with_ltt3218_output + ".flux", feige_spec1d_files)
        flux_and_validate(flux_setup_ltt3218_with_feige110_output + ".flux", ltt3218_spec1d_files)

        # Now update the coadd1d files. We'll use the files from flux_setup updated with the parameters
        # from the $PYPEIT_DEV/coadd1d_files
        updated_feige110_filename = f'updated_{flux_setup_feige110_with_ltt3218_output}.coadd1d'
        feige110_premade_coadd1d_file = Path(os.environ["PYPEIT_DEV"], "coadd1d_files", 'vlt_xshooter_nir-vis_Feige110.coadd1d')
        coadd1d_file = Coadd1DFile.from_file(str(feige110_premade_coadd1d_file))
        feige110_coadd1d_file.config = coadd1d_file.config
        feige110_output = coadd1d_file.config['coadd1d']['coaddfile']
        feige110_coadd1d_file.write(str(updated_feige110_filename))

        updated_ltt3218_filename = f"updated_{flux_setup_ltt3218_with_feige110_output}.coadd1d"
        ltt3218_premade_coadd1d_file = Path(os.environ["PYPEIT_DEV"], "coadd1d_files", 'vlt_xshooter_nir-vis_LTT3218.coadd1d')
        coadd1d_file = Coadd1DFile.from_file(str(ltt3218_premade_coadd1d_file))
        ltt3218_coadd1d_file.config = coadd1d_file.config
        ltt3218_output = coadd1d_file.config['coadd1d']['coaddfile']
        ltt3218_coadd1d_file.write(str(updated_ltt3218_filename))

        # Now test the coadding
        coadd(updated_feige110_filename, feige110_output)
        coadd(updated_ltt3218_filename, ltt3218_output)
