"""
Module to run tests on arcoadd
"""
import os

import pytest

from astropy.table import Table

from pypeit.coadd3d import CoAdd3D
from pypeit.spectrographs.util import load_spectrograph
from pypeit import inputfiles
from astropy.io import ascii

from IPython import embed

import warnings
warnings.simplefilter("ignore", UserWarning)


def test_coadd_datacube(redux_out):
    """ Test the coaddition of spec2D files into datacubes """
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
    parset['reduce']['cube']['combine'] = True    
    # Instantiate CoAdd3d, and then coadd the frames
    coadd = CoAdd3D.get_instance(coadd3dfile.filenames, coadd3dfile.options, spectrograph=spec, par=parset, overwrite=True)
    coadd.coadd()
    # Now test the fluxing - make a shorter set of files to speed it up
    files = ['spec2d_KB.20191219.56886-BB1245p4238_KCWI_20191219T154806.538.fits']
    tbl = Table()
    tbl['filename'] = files
    #tbl = ascii.read([files], header_start=0, data_start=1, delimiter='|', format='basic')

    # Generate a mock coadd3dfile
    coadd3dfile = inputfiles.Coadd3DFile(config=config,
                                         file_paths=[droot],
                                         data_table=tbl,
                                         setup=None)
    output_fileflux = "BB1245p4238_KCWI_20191219_fluxing.fits"
    parset['reduce']['cube']['output_filename'] = output_fileflux
    parset['reduce']['cube']['combine'] = False
    parset['reduce']['cube']['standard_cube'] = output_filename

    # Instantiate CoAdd3d, and then coadd the frames
    coadd = CoAdd3D.get_instance(coadd3dfile.filenames, coadd3dfile.options, spectrograph=spec, par=parset, overwrite=True)
    coadd.coadd()

    # Check the files exist
    assert(os.path.exists(output_filename))
    assert(os.path.exists(output_fileflux))
    # Remove the created files
    os.remove(output_filename)
    os.remove(output_fileflux)


