"""
Module to run tests on arcoadd
"""
import os
from pathlib import Path
import shutil

from astropy.table import Table
from IPython import embed
import numpy as np
import pytest

from pypeit.coadd3d import CoAdd3D, DataCube
from pypeit.scripts.extract_datacube import ExtractDataCube
from pypeit.scripts.sensfunc import SensFunc
from pypeit.spectrographs.util import load_spectrograph
from pypeit import inputfiles, specobjs, utils
from pypeit.core import standard


def test_coadd_datacube(redux_out):
    """ Test the coaddition of spec2D files into datacubes """
    # Setup the dev path
    dev_path = Path(os.getenv('PYPEIT_DEV')).absolute()
    outdir = dev_path / 'REDUX_OUT_TEST'
    if outdir.is_dir():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    # Define the input files
    droot = Path(redux_out).absolute() / 'keck_kcwi' / 'small_bh2_4200' / 'Science'
    files = [
        'spec2d_KB.20191219.56886-BB1245p4238_KCWI_20191219T154806.538.fits',
        'spec2d_KB.20191219.57662-BB1245p4238_KCWI_20191219T160102.755.fits'
    ]
    config = [
        '[rdx]',
        '  spectrograph = keck_kcwi'
    ]
    output_filename = "BB1245p4238_KCWI_20191219.fits"
    # Fake data table
    tbl = Table()
    tbl['filename'] = files

    # Generate a mock coadd3dfile
    coadd3dfile = inputfiles.Coadd3DFile(
        config=config, file_paths=[droot], data_table=tbl, setup=None
    )
    # Grab the spectrograph and parset
    spec = load_spectrograph("keck_kcwi")
    parset = spec.default_pypeit_par()
    parset['reduce']['cube']['output_filename'] = output_filename
    parset['reduce']['cube']['align'] = True
    parset['reduce']['cube']['combine'] = True
    parset['reduce']['cube']['weight_method'] = 'relative'

    # Speed up the computation by reducing the number of subpixels
    parset['reduce']['cube']['spat_subpixel'] = 3
    parset['reduce']['cube']['spec_subpixel'] = 1
    parset['reduce']['cube']['slice_subpixel'] = 3
    parset['reduce']['cube']['wave_min'] = 3922.758514
    parset['reduce']['cube']['wave_max'] = 4469.062985
    parset['reduce']['cube']['wave_delta'] = 0.115005

    # Extract the options
    # TODO: Why are these here?  At least at the moment, all the options are
    # None.
    ra_offsets = coadd3dfile.options['ra_offset']
    dec_offsets = coadd3dfile.options['dec_offset']
    skysub_frame = coadd3dfile.options['skysub_frame']
    scale_corr = coadd3dfile.options['scale_corr']
    grating_corr = coadd3dfile.options['grating_corr']
    sensfuncfile = coadd3dfile.options['sensfile']

    # Instantiate CoAdd3d, and then coadd the frames
    coadd = CoAdd3D.get_instance(
        coadd3dfile.filenames, parset, output_dir=str(outdir), skysub_frame=skysub_frame,
        grating_corr=grating_corr, scale_corr=scale_corr, sensfile=sensfuncfile,
        ra_offsets=ra_offsets, dec_offsets=dec_offsets, spectrograph=spec, overwrite=True
    )
    coadd.run()
    # Check the file exists
    cube_file = coadd.scidir / output_filename
    assert cube_file.is_file(), 'output file not written'

    ######################################
    # Test the extraction of a 1D spectrum from the datacube
    # Prepare the output filename
    pargs = ExtractDataCube.parse_args([
        str(coadd.scidir / output_filename), "-o", "-s", output_filename
    ])
    ExtractDataCube.main(pargs)
    # Expected output files
    ext_spec1d_file = coadd.scidir / f'spec1d_{output_filename}'
    ext_spec2d_file = coadd.scidir / f'spec2d_{output_filename}'
    # Check the files exist
    assert ext_spec1d_file.is_file(), 'spec1d file not written'
    assert ext_spec2d_file.is_file(), 'spec2d file not written'
    # TODO: Check for the whitelight file?

    ######################################
    # Using the 1D spectrum, generate a sensitivity function to flux the datacube
    # Prepare the output filename
    outfile_sens = ext_spec1d_file.name.replace('.fits', '_sens.fits')
    input_senspar = dev_path / 'sensfunc_files' / 'keck_kcwi_small_bh2_4200.sens'
    pargs = SensFunc.parse_args([str(ext_spec1d_file), "-o", outfile_sens, "-s", str(input_senspar)])
    SensFunc.main(pargs)
    # Check the files exist
    assert Path(outfile_sens).is_file(), 'Sensitivity function not written'

    ######################################
    # Now test the fluxing of the datacube
    fluxed_output_file = "BB1245p4238_KCWI_20191219_fluxing.fits"
    parset['reduce']['cube']['output_filename'] = fluxed_output_file
    parset['reduce']['cube']['sensfile'] = outfile_sens

    # Extract the options
    ra_offsets = coadd3dfile.options['ra_offset']
    dec_offsets = coadd3dfile.options['dec_offset']
    skysub_frame = coadd3dfile.options['skysub_frame']
    scale_corr = coadd3dfile.options['scale_corr']

    # Instantiate CoAdd3d, and then coadd the frames
    coadd = CoAdd3D.get_instance(
        coadd3dfile.filenames, parset, skysub_frame=skysub_frame, grating_corr=grating_corr,
        scale_corr=scale_corr, sensfile=outfile_sens, ra_offsets=ra_offsets,
        dec_offsets=dec_offsets, spectrograph=spec, overwrite=True
    )

    coadd.run()

    # Check the files exist
    fluxed_cube_file = coadd.scidir / fluxed_output_file
    assert fluxed_cube_file.is_file(), 'output fluxed file not written'

    ######################################
    # Finally, test the extraction of a 1D fluxed spectrum from the datacube
    # Prepare the output filename

    pargs = ExtractDataCube.parse_args([
        str(coadd.scidir / fluxed_output_file), "-o", "-s", fluxed_output_file
    ])
    ExtractDataCube.main(pargs)
    # Expected output files
    fluxed_ext_spec1d_file = coadd.scidir / f'spec1d_{fluxed_output_file}'
    fluxed_ext_spec2d_file = coadd.scidir / f'spec2d_{fluxed_output_file}'
    # Check the files exist
    assert fluxed_ext_spec1d_file.is_file(), 'fluxed_spec1d file not written'
    assert fluxed_ext_spec2d_file.is_file(), 'fluxed_spec2d file not written'
    # TODO: Check for the whitelight file?

    ######################################
    # Load in the extracted spec1d file, and compare it to the expected values
    spec1d = specobjs.SpecObjs.from_fitsfile(fluxed_ext_spec1d_file)
    # Generate a spectrum of the standard star that was used to generate the sensitivity function
    # Load in the standard star spectrum
    ra, dec = 191.39844, 42.64016
    std_spec = standard.get_archive_standard(ra, dec, archives='blackbody')
    wave_std, flux_std = std_spec.wave, std_spec.flux
    # Test the optimal extraction
    # Interpolate the standard star spectrum to the same wavelength grid as the spec1d
    flux_std_interp = np.interp(spec1d[0].OPT_WAVE, wave_std, flux_std)
    # Compare the extracted spectrum to the standard star spectrum, and make sure that the residuals are small
    resid = (spec1d[0].OPT_FLAM-flux_std_interp)*utils.inverse(spec1d[0].OPT_FLAM_SIG)
    med, std = np.median(resid), 1.4826*np.median(np.abs(np.median(resid) - resid))
    assert(np.abs(med) < 0.2*std)
    # Test the boxcar extraction
    # Interpolate the standard star spectrum to the same wavelength grid as the spec1d
    flux_std_interp = np.interp(spec1d[0].BOX_WAVE, wave_std, flux_std)
    # Compare the extracted spectrum to the standard star spectrum, and make sure that the residuals are small
    resid = (spec1d[0].BOX_FLAM-flux_std_interp)*utils.inverse(spec1d[0].BOX_FLAM_SIG)
    med, std = np.median(resid), 1.4826*np.median(np.abs(np.median(resid) - resid))
    # The sensitivity function is based on the optimal extraction, so the optimal should be spot on.
    # The boxcar will be worse, so allow a larger tolerance
    assert(np.abs(med) < std)

    ######################################
    # Remove all of the created files
    shutil.rmtree(outdir)

#    # First remove the non-fluxed files
#    os.remove(output_filename)
#    os.remove(output1d_filename)
#    # Remove the sensitivity function files and the associated QA files
#    os.remove(outfile_sens)
#    os.remove('sensfunc.par')
#    os.remove(outfile_sens.replace('.fits', '_QA.pdf'))
#    os.remove(outfile_sens.replace('.fits', '_throughput.pdf'))
#    os.remove(outfile_sens.replace('.fits', '_fluxed_std.pdf'))
#    # Remove the fluxed files
#    os.remove(output_fileflux)
#    os.remove(output1d_fileflux)

def test_residuals(redux_out):
    """ Test the residuals of a spec2D DOMEFLAT file
    """
    # TODO: Fix as part of https://github.com/pypeit/PypeIt/issues/1951
    # Fixing the `assert` checks caused this test to fail.  See the issue link
    #   for a description of what needs to be dealt with in the main PypeIt repo
    #   before this test can be reinstated.
    return

    # Define the input files
    droot = os.path.join(redux_out,
                         'keck_kcwi',
                         'small_bh2_4200',
                         'Science')
    files = ['spec2d_KB.20191220.62342-DOMEPHLAT_KCWI_20191220T171902.438.fits']
    config = ['[rdx]',
              '  spectrograph = keck_kcwi']

    output_filename = "DOMEFLAT_BH2_333.fits"
    # Fake data table
    tbl = Table()
    tbl['filename'] = files

    # Generate a mock coadd3dfile
    coadd3dfile = inputfiles.Coadd3DFile(config=config,
                                         file_paths=[droot],
                                         data_table=tbl,
                                         setup=None)
    # Grab the spectrograph and parset
    spec = load_spectrograph("keck_kcwi")
    parset = spec.default_pypeit_par()
    parset['reduce']['cube']['output_filename'] = output_filename
    parset['reduce']['cube']['correct_dar'] = False
    parset['reduce']['cube']['combine'] = False
    parset['reduce']['cube']['align'] = False
    parset['reduce']['cube']['weight_method'] = 'relative'
    parset['reduce']['cube']['method'] = 'subpixel'
    parset['reduce']['cube']['spat_subpixel'] = 3
    parset['reduce']['cube']['spec_subpixel'] = 3
    parset['reduce']['cube']['slice_subpixel'] = 3
    parset['reduce']['cube']['wave_min'] = 3922.758514
    parset['reduce']['cube']['wave_max'] = 4469.062985
    parset['reduce']['cube']['wave_delta'] = 0.115005

    # Extract the options
    ra_offsets = coadd3dfile.options['ra_offset']
    dec_offsets = coadd3dfile.options['dec_offset']
    skysub_frame = coadd3dfile.options['skysub_frame']
    scale_corr = coadd3dfile.options['scale_corr']
    grating_corr = coadd3dfile.options['grating_corr']
    sensfuncfile = coadd3dfile.options['sensfile']

    # Instantiate CoAdd3d, and then coadd the frames
    coadd = CoAdd3D.get_instance(coadd3dfile.filenames, parset, skysub_frame=skysub_frame, grating_corr=grating_corr,
                                 scale_corr=scale_corr, sensfile=sensfuncfile,
                                 ra_offsets=ra_offsets, dec_offsets=dec_offsets, spectrograph=spec, overwrite=True)
    coadd.run()
    # Check the file exists
    assert(os.path.exists(output_filename))
    ######################################
    # Check the residuals are OK for method=subpixel
    cube = DataCube.from_file(output_filename)
    ww = np.where(cube['bpm'] == 0)
    resid = cube['flux'] * utils.inverse(cube['sig'])
    # Calculate the statistics
    avg, med = np.mean(resid[ww]), np.median(resid[ww])
    std, mad = np.std(resid[ww]), 1.4826 * np.median(np.abs(np.median(resid[ww]) - resid[ww]))
    # Check the statistics
    assert np.abs(avg) < 0.1, 'residuals (average) is not close to zero for method=subpixel(333)'
    assert np.abs(med) < 0.1, 'residuals (median) is not close to zero for method=subpixel(333)'
    assert np.abs(std-1) < 0.1, 'residuals (std) is not close to 1 for method=subpixel(333)'
    assert np.abs(mad-1) < 0.1, 'residuals (1.4826 * mad) is not close to 1 for method=subpixel(333)'

    ######################################
    # Now check the NGP algorithm
    output_fileNGP = "DOMEFLAT_BH2_NGP.fits"
    parset['reduce']['cube']['output_filename'] = output_fileNGP
    parset['reduce']['cube']['method'] = 'ngp'
    parset['reduce']['cube']['spat_subpixel'] = 1
    parset['reduce']['cube']['spec_subpixel'] = 1
    parset['reduce']['cube']['slice_subpixel'] = 1

    # Extract the options
    ra_offsets = coadd3dfile.options['ra_offset']
    dec_offsets = coadd3dfile.options['dec_offset']
    skysub_frame = coadd3dfile.options['skysub_frame']
    scale_corr = coadd3dfile.options['scale_corr']

    # Instantiate CoAdd3d, and then coadd the frames
    coadd = CoAdd3D.get_instance(coadd3dfile.filenames, parset, skysub_frame=skysub_frame, grating_corr=grating_corr,
                                 scale_corr=scale_corr, sensfile=sensfuncfile,
                                 ra_offsets=ra_offsets, dec_offsets=dec_offsets, spectrograph=spec, overwrite=True)
    coadd.run()

    # Check the files exist
    assert(os.path.exists(output_fileNGP))
    ######################################
    # Check the residuals are OK for method=NGP
    cube = DataCube.from_file(output_fileNGP)
    ww = np.where(cube['bpm'] == 0)
    resid = cube['flux'] * utils.inverse(cube['sig'])
    # Calculate the statistics
    avg, med = np.mean(resid[ww]), np.median(resid[ww])
    std, mad = np.std(resid[ww]), 1.4826 * np.median(np.abs(np.median(resid[ww]) - resid[ww]))
    # Check the statistics
    assert np.abs(avg) < 0.1, 'residuals (average) is not close to zero for method=NGP'
    assert np.abs(med) < 0.1, 'residuals (median) is not close to zero for method=NGP'
    assert np.abs(std-1) < 0.1, 'residuals (std) is not close to 1 for method=NGP'
    assert np.abs(mad-1) < 0.1, 'residuals (1.4826 * mad) is not close to 1 for method=NGP'
    ######################################
    # Remove all of the created files
    os.remove(output_filename)
    os.remove(output_fileNGP)
