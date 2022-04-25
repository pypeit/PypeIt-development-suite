import os
import numpy as np
from pypeit import spec2dobj

from pypeit.core import flexure
from pypeit.spectrographs.util import load_spectrograph

from pypeit.tests.tstutils import data_path

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
    spec2dObj = spec2dobj.Spec2DObj.from_file(file_path, 'DET01')
    assert spec2dObj.sci_spat_flexure is not None
    assert spec2dObj.sci_spat_flexure > 0.


def test_flex_multi():

    # Set output file
    outfile = data_path('tst_multi_flex.fits')
    if os.path.isfile(outfile):
        # Remove it if it already exists
        os.remove(outfile)

    spec1d_file = os.path.join(os.getenv('PYPEIT_DEV'), 
                             'REDUX_OUT', 
                             'keck_deimos',
                             '830G_M_8500', 
                             'Science', 
                             'spec1d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits')

    msFlex = flexure.MultiSlitFlexure(s1dfile=spec1d_file) 
    # Parameters
    keck_deimos = load_spectrograph('keck_deimos')
    par = keck_deimos.default_pypeit_par()
    # Init                    
    msFlex.init(keck_deimos, par['flexure'])
    # INITIAL SKY LINE STUFF
    msFlex.measure_sky_lines()
    # FIT SURFACES
    msFlex.fit_mask_surfaces()
    # Apply
    msFlex.update_fit()
    # QA
    #mask = header['TARGET'].strip()
    #fnames = header['FILENAME'].split('.')
    #root = mask+'_'+fnames[2]
    #mdFlex.qa_plots('./', root)

    # Write
    msFlex.to_file(outfile, overwrite=True)

    # Read
    msFlex2 = flexure.MultiSlitFlexure.from_file(outfile)
    # Check
    assert np.array_equal(msFlex2.fit_b, msFlex.fit_b), 'Bad read'

    # Try to overwrite
    msFlex2.to_file(outfile, overwrite=True)

    # Clean up
    if os.path.isfile(outfile):
        os.remove(outfile)

