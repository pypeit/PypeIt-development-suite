"""
Module to run QL function tests
"""
from pathlib import Path
import os
import shutil

from IPython import embed

import pytest

import numpy as np

from pypeit import pypeitsetup
from pypeit.pypmsgs import PypeItError
from pypeit.scripts import ql

def test_dither_parse():

    dev_root = Path(os.environ['PYPEIT_DEV']).resolve()
    data_root = dev_root / 'RAW_DATA' / 'keck_nires'

    files = sorted(list((data_root / 'ABBA_nostandard').glob('*fits*')))
    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_nires')
    ps.run(setup_only=True)

    # This should be raised because the current directory contains science
    # observations of multiple targets.
    with pytest.raises(PypeItError):
        ql.quicklook_regroup(ps.fitstbl)

    # Remove observations of "eris_28 0620"
    ps.remove_table_rows(ps.fitstbl['target'] == 'eris_28 0620', regroup=True)
    ql.quicklook_regroup(ps.fitstbl) 

    # Test the result
    is_sci = ps.fitstbl.find_frames('science')
    assert np.array_equal(np.unique(ps.fitstbl['comb_id'].data[is_sci]), np.array([0,1])), \
            'There should only be 2 combination IDs.'
    is_std = ps.fitstbl.find_frames('standard')
    assert all(ps.fitstbl['bkg_id'].data[is_std] == -1), \
            'No background frames were taken for the standard.'
    assert all(ps.fitstbl['comb_id'].data[is_std] == 2), \
            'All standards should be combined.'

    # Rebuild
    ps = pypeitsetup.PypeItSetup(files, spectrograph_name='keck_nires')
    ps.run(setup_only=True)
    # Remove observations of "eris_28 0620"
    ps.remove_table_rows(ps.fitstbl['target'] == 'eris_28 0620', regroup=True)
    # Fake putting one ABBA sequence at a different offset
    last_dither = ps.fitstbl.find_frames('science', index=True)[-4:]
    ps.fitstbl['dithoff'].data[last_dither] = [3.0, -3.0, -3.0, 3.0]

    # Regroup
    ql.quicklook_regroup(ps.fitstbl) 

    # Test the result
    is_sci = ps.fitstbl.find_frames('science')
    assert np.array_equal(np.unique(ps.fitstbl['comb_id'].data[is_sci]), np.array([0,1,2,3])), \
            'There should only be 2 combination IDs.'
    is_std = ps.fitstbl.find_frames('standard')
    assert all(ps.fitstbl['bkg_id'].data[is_std] == -1), \
            'No background frames were taken for the standard.'
    assert all(ps.fitstbl['comb_id'].data[is_std] == 4), \
            'All standards should be combined.'

    # TODO:
    #   - Test if all dithers are unique?
    #   - Test for incomplete sequences?


def test_run_ql():
    """
    Test a basic execution of QL that only requires raw files.  This could also
    go in test_scripts.py ...
    """
    rawpath = Path(os.environ['PYPEIT_DEV']).resolve() \
                    / 'RAW_DATA' / 'shane_kast_blue' / '600_4310_d55'
    rdxpath = Path().resolve() / 'test_ql'

    if rdxpath.is_dir():
        shutil.rmtree(rdxpath)

    ql.QL.main(ql.QL.parse_args(['shane_kast_blue', '--skip_display', '--redux_path', str(rdxpath),
                                 '--raw_path', str(rawpath), '--raw_files', 'b1.fits.gz',
                                 'b10.fits.gz', 'b27.fits.gz']))

    shutil.rmtree(rdxpath)



