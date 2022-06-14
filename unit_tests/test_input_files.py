""" Tests Reading of PypeIt Input files """

import os
import glob

import pytest

from pypeit import inputfiles


def test_fluxing_files():
    # Grab em
    fluxing_files = glob.glob(os.path.join(
        os.getenv('PYPEIT_DEV'), 'fluxing_files', '*.flux'))
    # Loop
    for ifile in fluxing_files:
        fluxFile = inputfiles.FluxFile.from_file(ifile)