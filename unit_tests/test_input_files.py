""" Tests Reading of PypeIt Input files """

import os
import glob

import pytest

from pypeit import inputfiles


def test_read_fluxing_files():
    """ Test reading fluxing files """
    # Grab em
    fluxing_files = glob.glob(os.path.join(
        os.getenv('PYPEIT_DEV'), 'fluxing_files', '*.flux'))
    # Loop
    for ifile in fluxing_files:
        fluxFile = inputfiles.FluxFile.from_file(ifile)

def test_read_coadd1d_files():
    """ Test reading coadd1d files """
    # Grab em
    coadd1d_files = glob.glob(os.path.join(
        os.getenv('PYPEIT_DEV'), 'coadd1d_files', '*.coadd1d'))
    # Loop
    for ifile in coadd1d_files:
        coadd1dFile = inputfiles.Coadd1DFile.from_file(ifile)

def test_read_coadd2d_files():
    """ Test reading coadd2d files """
    # Grab em
    coadd2d_files = glob.glob(os.path.join(
        os.getenv('PYPEIT_DEV'), 'coadd2d_files', '*.coadd1d'))
    # Loop
    for ifile in coadd2d_files:
        coadd2dFile = inputfiles.Coadd2DFile.from_file(ifile)

def test_read_pypeit_files():
    """ Test reading PypeIt files """
    # Grab em
    pypeit_files = glob.glob(os.path.join(
        os.getenv('PYPEIT_DEV'), 'pypeit_files', '*.pypeit'))
    # Loop
    for ifile in pypeit_files:
        pypeitFile = inputfiles.PypeItFile.from_file(ifile)