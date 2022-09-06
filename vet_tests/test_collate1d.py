import pytest
import os.path
import glob

import numpy as np
from astropy.table import Table

from pypeit.specobjs import SpecObjs

def test_collate_1d(redux_out):

    # Test that coadd files exist
    output_dir = os.path.join(redux_out, 'keck_deimos', '830G_M_8500')
    coadd_pattern = os.path.join(output_dir, 'J*.fits')

    coadd_files = glob.glob(coadd_pattern)

    assert len(coadd_files) > 1

    # Load the collate_report.dat file, and verify that wave rms's are under threshold
    collate_report = Table.read(os.path.join(output_dir, 'collate_report.dat'), format='ipac')

    assert np.sum(collate_report['wave_rms'] >= 0.1) == 0

    # Check that two of the higher s/n  objects are matched and coadded to the same file
    obj_idx = collate_report['maskdef_objname'] == 'ero89'
    assert np.sum(obj_idx) == 2
    assert collate_report[obj_idx][0]['filename'] == collate_report[obj_idx][1]['filename']

    obj_idx = collate_report['maskdef_objname'] == 'ero234'
    assert np.sum(obj_idx) == 2
    assert collate_report[obj_idx][0]['filename'] == collate_report[obj_idx][1]['filename']

    # Also assert the warnings file has at least one warning of rejecting rms
    with open(os.path.join(output_dir, "collate_warnings.txt")) as f:
        collate_warnings = f.read()
        assert "due to wave_rms" in collate_warnings

    # Verify flux and helio correction was applied to spec1d
    spec1d_files = glob.glob(os.path.join(output_dir, "spec1d_DE.20100913.22358*.fits"))

    # We only expect one spec1d to match this
    assert len(spec1d_files) == 1
    sobjs = SpecObjs.from_fitsfile(spec1d_files[0])
    for sobj in sobjs:
        assert sobj['OPT_FLAM'] is not None
        assert sobj['BOX_FLAM'] is not None
        assert sobj['VEL_TYPE'] == 'heliocentric'
        assert sobj['VEL_CORR'] is not None

