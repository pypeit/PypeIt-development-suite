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

def chk_orders(instr, redux_out, det='DET01', max_bad:int=0):
    for setup in setups._all[instr]:
        # Grab a spec2d file
        file_path = os.path.join(redux_out,
            instr, setup, 'Science',
            'spec2d*fits')
        spec2d_files = glob.glob(file_path)
        assert len(spec2d_files) > 0, f'No spec2d files found for {setup} and instr={instr}'
        # Take one
        spec2d_file = spec2d_files[0]
        # Load
        spec2d = spec2dobj.Spec2DObj.from_file(spec2d_file, det)
        assert np.sum(spec2d.slits.mask != 0) <= max_bad, f'Bad order(s) for {setup}'

def test_vlt_xshooter_orders(redux_out):
    """ Confirm that all of the orders processed fine for each setup"""
    instr = 'vlt_xshooter'
    chk_orders(instr, redux_out)

def test_magellan_mage_orders(redux_out):
    """ Confirm that all of the orders processed fine for each setup"""
    instr = 'magellan_mage'

    chk_orders(instr, redux_out)

#def test_keck_hires_orders(redux_out):
#    """ Confirm that all of the orders processed fine for each setup"""
#    instr = 'keck_hires'
#
#    # Some orders are rightly rejected
#    chk_orders(instr, redux_out, det='MSC01', max_bad=3)

def test_keck_nires_orders(redux_out):
    """ Confirm that all of the orders processed fine for each setup"""
    instr = 'keck_nires'

    chk_orders(instr, redux_out)

def test_gemini_gnirs_orders(redux_out):
    """ Confirm that all of the orders processed fine for each setup"""
    instr = 'gemini_gnirs'

    chk_orders(instr, redux_out)

def test_magellan_fire_orders(redux_out):
    """ Confirm that all of the orders processed fine for each setup"""
    instr = 'magellan_fire'

    chk_orders(instr, redux_out)