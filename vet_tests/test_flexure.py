import os
from pypeit import spec2dobj

import pytest

def test_spat_flexure():
    # Check that spatial flexure shift was set!
    file_path = os.path.join(os.environ['PYPEIT_DEV'],
                             'REDUX_OUT',
                             'keck_lris_red', 
                             'multi_600_5000_d560',
                             'Science', 
                             'spec2d_LR.20181206.40617-nR2n25061_LRISr_20181206T111657.418.fits')
    # Load                                
    spec2dObj = spec2dobj.Spec2DObj.from_file(file_path, det=1)
    assert spec2dObj.sci_spat_flexure is not None
    assert spec2dObj.sci_spat_flexure > 0.
