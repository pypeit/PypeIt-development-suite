"""
Module to run tests on arcoadd
"""
import os

import pytest

import numpy as np
from astropy.table import Table
from astropy.io import ascii

from pypeit.coadd3d import CoAdd3D, DataCube
from pypeit.scripts.extract_datacube import ExtractDataCube
from pypeit.scripts.sensfunc import SensFunc
from pypeit.spectrographs.util import load_spectrograph
from pypeit import inputfiles, specobjs, utils
from pypeit.core import flux_calib

from IPython import embed

import warnings
warnings.simplefilter("ignore", UserWarning)


def test_coadd_datacube(redux_out):
    """ Test the coaddition of spec2D files into datacubes """
    # Setup the dev path
    dev_path = os.getenv('PYPEIT_DEV')
    # Define the input files
    droot = os.path.join(redux_out,
                         'keck_kcwi', 
                         'small_bh2_4200',
                         'Science')
    files = ['spec2d_KB.20191219.56886-BB1245p4238_KCWI_20191219T154806.538.fits',
             'spec2d_KB.20191219.57662-BB1245p4238_KCWI_20191219T160102.755.fits']
    config = ['[rdx]',
              '  spectrograph = keck_kcwi']
    output_filename = "BB1245p4238_KCWI_20191219.fits"
    # Fake data table
    #tbl = ascii.read([files], header_start=0, data_start=1, delimiter='|', format='basic')
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
    # Test the extraction of a 1D spectrum from the datacube
    # Prepare the output filename
    output1d_filename = output_filename.replace('.fits', '_spec1d.fits')
    pargs = ExtractDataCube.parse_args(["-o", "-s", output1d_filename, output_filename])
    ExtractDataCube.main(pargs)
    # Check the files exist
    assert(os.path.exists(output1d_filename))
    ######################################
    # Using the 1D spectrum, generate a sensitivity function to flux the datacube
    # Prepare the output filename
    outfile_sens = output1d_filename.replace('.fits', '_sens.fits')
    input_senspar = os.path.join(dev_path, 'sensfunc_files', 'keck_kcwi_small_bh2_4200.sens')
    pargs = SensFunc.parse_args(["-o", outfile_sens, "-s", input_senspar, output1d_filename])
    SensFunc.main(pargs)
    # Check the files exist
    assert(os.path.exists(outfile_sens))
    ######################################
    # Now test the fluxing of the datacube
    output_fileflux = "BB1245p4238_KCWI_20191219_fluxing.fits"
    parset['reduce']['cube']['output_filename'] = output_fileflux
    parset['reduce']['cube']['sensfile'] = outfile_sens

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
    assert(os.path.exists(output_fileflux))
    ######################################
    # Finally, test the extraction of a 1D fluxed spectrum from the datacube
    # Prepare the output filename
    output1d_fileflux = output_fileflux.replace('.fits', '_spec1d.fits')
    pargs = ExtractDataCube.parse_args(["-o", "-s", output1d_fileflux, output_fileflux])
    ExtractDataCube.main(pargs)
    # Check the files exist
    assert(os.path.exists(output1d_fileflux))
    ######################################
    # Load in the extracted spec1d file, and compare it to the expected values
    spec1d = specobjs.SpecObjs.from_fitsfile(output1d_fileflux)
    # Generate a spectrum of the standard star that was used to generate the sensitivity function
    # Load in the standard star spectrum
    ra, dec = 191.39844, 42.64016
    std_dict = flux_calib.find_standard_file(ra, dec)
    wave_std, flux_std = std_dict['wave'].value, std_dict['flux'].value
    # Test the optimal extraction
    # Interpolate the standard star spectrum to the same wavelength grid as the spec1d
    flux_std_interp = np.interp(spec1d[0].OPT_WAVE, wave_std, flux_std)
    # Compare the extracted spectrum to the standard star spectrum, and make sure that the residuals are small
    resid = (spec1d[0].OPT_FLAM-flux_std_interp)*utils.inverse(spec1d[0].OPT_FLAM_SIG)
    med, std = np.median(resid), 1.4826*np.median(np.abs(np.median(resid) - resid))
    assert(np.abs(med) < 0.1*std)
    # Test the boxcar extraction
    # Interpolate the standard star spectrum to the same wavelength grid as the spec1d
    flux_std_interp = np.interp(spec1d[0].BOX_WAVE, wave_std, flux_std)
    # Compare the extracted spectrum to the standard star spectrum, and make sure that the residuals are small
    resid = (spec1d[0].BOX_FLAM-flux_std_interp)*utils.inverse(spec1d[0].BOX_FLAM_SIG)
    med, std = np.median(resid), 1.4826*np.median(np.abs(np.median(resid) - resid))
    # The sensitivity function is based on the optimal extraction, so the optimal should be spot on.
    # The boxcar will be worse, so allow a larger tolerance
    assert(np.abs(med) < 0.5*std)
    ######################################
    # Remove all of the created files
    # First remove the non-fluxed files
    os.remove(output_filename)
    os.remove(output1d_filename)
    # Remove the sensitivity function files and the associated QA files
    os.remove(outfile_sens)
    os.remove('sensfunc.par')
    os.remove(outfile_sens.replace('.fits', '_QA.pdf'))
    os.remove(outfile_sens.replace('.fits', '_throughput.pdf'))
    os.remove(outfile_sens.replace('.fits', '_fluxed_std.pdf'))
    # Remove the fluxed files
    os.remove(output_fileflux)
    os.remove(output1d_fileflux)
