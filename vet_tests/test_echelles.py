""" Tests on echelle spectrographs """
import os
import sys
import glob
import numpy as np

from pypeit import spec2dobj

# THIS MIGHT NOT WORK
sys.path.append(os.path.abspath("../test_scripts"))
import setups

import pytest

def test_vlt_xshooter_orders(redux_out):
    """ Confirm that all of the orders processed fine for each setup"""

    for setup in setups._all['vlt_xshooter']:
        # Grab a spec2d file
        file_path = os.path.join(redux_out,
            'vlt_xshooter', setup, 'Science',
            'spec2d*fits')
        spec2d_files = glob.glob(file_path)
        assert len(spec2d_files) > 0, f'No spec2d files found for {setup}'
        # Take one
        spec2d_file = spec2d_files[0]
        # Load
        spec2d = spec2dobj.Spec2DObj.from_file(spec2d_file, 'DET01')
        assert np.all(spec2d.slits.mask == 0), f'Bad order(s) for {setup}'
        

