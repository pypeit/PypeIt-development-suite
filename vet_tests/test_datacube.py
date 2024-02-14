"""
Module to run tests on arcoadd
"""
import os

import pytest

from astropy.table import Table

from pypeit.coadd3d import CoAdd3D, DataCube
from pypeit.scripts.extract_datacube import ExtractDataCube
from pypeit.scripts.sensfunc import SensFunc
from pypeit.spectrographs.util import load_spectrograph
from pypeit import inputfiles
from astropy.io import ascii

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
    droot = "/Users/rcooke/Work/Research/BBN/Yp/HIIregions/IZw18_KCWI/keck_kcwi_BH2_N1/Science"
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
    parset['reduce']['cube']['combine'] = True
    # Speed up the computation by reducing the number of subpixels
    parset['reduce']['cube']['spat_subpixel'] = 3
    parset['reduce']['cube']['spec_subpixel'] = 1
    parset['reduce']['cube']['slice_subpixel'] = 3

    # Extract the options
    ra_offsets = coadd3dfile.options['ra_offset']
    dec_offsets = coadd3dfile.options['dec_offset']
    skysub_frame = coadd3dfile.options['skysub_frame']
    scale_corr = coadd3dfile.options['scale_corr']

    # Instantiate CoAdd3d, and then coadd the frames
    coadd = CoAdd3D.get_instance(coadd3dfile.filenames, parset, skysub_frame=skysub_frame, scale_corr=scale_corr,
                                 ra_offsets=ra_offsets, dec_offsets=dec_offsets, spectrograph=spec, overwrite=True)
    coadd.run()
    # Check the file exists
    assert(os.path.exists(output_filename))
    ######################################
    # Test the extraction of a 1D spectrum from the datacube
    # Prepare the output filename
    output1d_filename = output_filename.replace('.fits', '_spec1d.fits')
    pargs = ExtractDataCube.parse_args(["-o", output1d_filename, output_filename])
    ExtractDataCube.main(pargs)
    # Check the files exist
    assert(os.path.exists(output1d_filename))
    ######################################
    # Using the 1D spectrum, generate a sensitivity function to flux the datacube
    # Prepare the output filename
    outfile_sens = output1d_filename.replace('.fits', '_sens.fits')
    input_senspar = os.path.join(dev_path, 'sensfunc_files', 'keck_kcwi_small_bh2_4200.sens')
    pargs = SensFunc.parse_args(["--algorithm", "UVIS", "-o", outfile_sens, "-s", input_senspar, output1d_filename])
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
    coadd = CoAdd3D.get_instance(coadd3dfile.filenames, parset, skysub_frame=skysub_frame, scale_corr=scale_corr,
                                 ra_offsets=ra_offsets, dec_offsets=dec_offsets, spectrograph=spec, overwrite=True)
    coadd.run()

    # Check the files exist
    assert(os.path.exists(output_fileflux))

    # Remove all of the created files
    os.remove(output_filename)
    os.remove(output1d_filename)
    os.remove(outfile_sens)
    os.remove(output_fileflux)
