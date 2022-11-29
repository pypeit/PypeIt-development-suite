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

def test_shane_kast():
    pass

    # Check standard
    # Check boxcar -- 4.65116279