import os
from pypeit import specobjs

import pytest

def test_lris_slitmask():
    # Check that the LRIS slitmask was read in and used!
    file_path = os.path.join(os.environ['PYPEIT_DEV'],
                             'REDUX_OUT',
                             'keck_lris_blue', 
                             'multi_600_4000_slitmask',
                             'Science', 
                             'spec1d_LB.20200129.48812-frb19071_LRISb_20200129T133332.890.fits')
    # Load                                
    specObjs = specobjs.SpecObjs.from_fitsfile(file_path)

    # Test
    assert len(specObjs.MASKDEF_ID) > 0
    assert 'gal21' in specObjs.MASKDEF_OBJNAME
    assert 'gal49' in specObjs.MASKDEF_OBJNAME # This was "manually" extracted