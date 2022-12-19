"""
Module to run tests on scripts
"""
import os
import glob
import numpy as np
import pytest

from pypeit.pypmsgs import PypeItError
from pypeit.inputfiles import PypeItFile
from pypeit import specobjs

def test_bok_bc_manual(redux_out):
    """ Checks that the manual extraction with FWHM is working for Bok BC"""
    instr = 'bok_bc' 
    rdxdir = os.path.join(redux_out, instr, '300')
    scidir = os.path.join(rdxdir, 'Science')

    spec1d_files = glob.glob(os.path.join(scidir, 'spec1d*.fits')) 
    spec1d_files.sort()
    sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_files[1]) # 0 should be the standard

    hand_sobj = sobjs[sobjs.hand_extract_flag]
    # Test
    assert np.isclose(hand_sobj.BOX_RADIUS[0], 4.)  # Value in the pypeit file
