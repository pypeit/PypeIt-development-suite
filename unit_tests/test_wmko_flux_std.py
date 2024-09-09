"""
Module to run tests on wmko flux stds
Requires files in PypeIt/pypeit/data and telluric files
"""
import os

import pytest

from IPython import embed

import numpy as np

from pypeit import dataPaths
from pypeit import sensfunc
from pypeit.tests.tstutils import data_output_path
from pypeit.spectrographs.util import load_spectrograph
from pypeit.spectrographs import keck_deimos


def test_wmko_flux_std():

    outfile = data_output_path('tmp_sens.fits')
    if os.path.isfile(outfile):
        os.remove(outfile)

    # Do it
    wmko_file = str(dataPaths.tests.get_file_path('2017may28_d0528_0088.fits'))
    spectrograph = load_spectrograph('keck_deimos')

    # Load + write
    spec1dfile = data_output_path('tmp_spec1d.fits')
    sobjs = keck_deimos.load_wmko_std_spectrum(wmko_file, outfile=spec1dfile)

    # Sensfunc
    #  The following mirrors the main() call of sensfunc.py
    par = spectrograph.default_pypeit_par()
    par['sensfunc']['algorithm'] = "IR"
    par['sensfunc']['multi_spec_det'] = [3,7]
    # Instantiate the relevant class for the requested algorithm
    sensobj = sensfunc.SensFunc.get_instance(spec1dfile, outfile, par['sensfunc'])
    # Generate the sensfunc
    sensobj.run()
    # Write it out to a file
    sensobj.to_file(outfile)

    os.remove(spec1dfile)
    os.remove(outfile)

