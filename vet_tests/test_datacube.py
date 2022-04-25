"""
Module to run tests on arcoadd
"""
import os

import pytest

from pypeit.core.datacube import coadd_cube

from IPython import embed

import warnings
warnings.simplefilter("ignore", UserWarning)

def test_coadd_datacube():
    """ Test the coaddition of spec2D files into datacubes """
    droot = os.path.join(os.getenv('PYPEIT_DEV'), 
                             'REDUX_OUT', 
                             'keck_kcwi', 
                             'bh2_4200', 
                             'Science')
    files = [os.path.join(droot,
                          'spec2d_KB.20191219.56886-BB1245p4238_KCWI_20191219T154806.538.fits'),
             os.path.join(droot,
                          'spec2d_KB.20191219.57662-BB1245p4238_KCWI_20191219T160102.755.fits')]
    coadd_cube(files, overwrite=True)
    os.remove('datacube.fits')


