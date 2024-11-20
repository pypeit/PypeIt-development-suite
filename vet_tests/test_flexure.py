import os
import numpy as np
from pypeit import spec2dobj
from pypeit.wavetilts import WaveTilts
from pypeit.images.buildimage import TiltImage
from pypeit.core import flexure
from pypeit.spectrographs.util import load_spectrograph
from pypeit import inputfiles

from pypeit.tests.tstutils import data_output_path
from pathlib import Path
from IPython import embed

import pytest

def test_spat_flexure(redux_out):
    # Check that spatial flexure shift was set!
    file_path = os.path.join(redux_out,
                             'keck_lris_red', 
                             'multi_600_5000_d560',
                             'Science', 
                             'spec2d_LR.20181206.40617-nR2n25061_LRISr_20181206T111657.418.fits')
    # Load                                
    spec2dObj = spec2dobj.Spec2DObj.from_file(file_path, 'DET01')
    assert spec2dObj.sci_spat_flexure is not None
    assert spec2dObj.sci_spat_flexure > 0.

def test_spat_flexure_science(redux_out):
    spec = 'keck_lris_red_mark4'
    setup = 'long_600_10000_d680'

    # Define the path with the raw data
    data_redux = Path(redux_out).resolve() / spec / setup
    assert data_redux.exists(), f'TEST ERROR: REDUX data path does not exist for {spec}/{setup}'

    # find pypeit file
    pypeit_file = data_redux / 'keck_lris_red_mark4_long_600_10000_d680.pypeit'
    assert pypeit_file.exists(), 'Missing test pypeit file'
    pypeitFile = inputfiles.PypeItFile.from_file(pypeit_file)
    _, par, _ = pypeitFile.get_pypeitpar()
    # check flexure parameters
    assert par['scienceframe']['process']['spat_flexure_correct'] == True, 'Spat flexure should be set to True'
    assert par['scienceframe']['process']['spat_flexure_sigdetect'] == 10.0, 'Spat flexure sigdetect should be 10.0'
    assert par['scienceframe']['process']['spat_flexure_maxlag'] == 10, 'Spat flexure maxlag should be 10'

    # read the spec2d file
    spec2d_files = list(data_redux.glob('Science/spec2d*.fits'))
    assert len(spec2d_files) > 0, 'No spec2d files found'
    spat_flex_list = []
    for spec2d_file in spec2d_files:
        spec2dObj = spec2dobj.Spec2DObj.from_file(spec2d_file, 'DET01')
        spat_flex_list.append(spec2dObj.sci_spat_flexure)
    spat_flex_list = np.array(spat_flex_list)
    assert np.all(spat_flex_list != None), 'Spat flexure should not be None'
    assert np.all(spat_flex_list < par['scienceframe']['process']['spat_flexure_maxlag']), \
        'Spat flexure should be less than maxlag'
    assert np.all(np.isclose(spat_flex_list, 6., atol=1.5)), 'Spat flexure should be close to 6 pixels'

def test_spat_flexure_tilts(redux_out):
    spec = 'keck_lris_red_mark4'
    setup = 'long_600_10000_d680'

    # Define the path with the raw data
    data_redux = Path(redux_out).resolve() / spec / setup
    assert data_redux.exists(), f'TEST ERROR: REDUX data path does not exist for {spec}/{setup}'

    # find pypeit file
    pypeit_file = data_redux / 'keck_lris_red_mark4_long_600_10000_d680.pypeit'
    assert pypeit_file.exists(), 'Missing test pypeit file'
    pypeitFile = inputfiles.PypeItFile.from_file(pypeit_file)
    _, par, _ = pypeitFile.get_pypeitpar()
    # check flexure parameters
    assert par['calibrations']['tiltframe']['process']['spat_flexure_correct'] == True, 'Spat flexure should be set to True'
    assert par['calibrations']['tiltframe']['process']['spat_flexure_maxlag'] == 10, 'Spat flexure maxlag should be 10'
    # read the tiltimg file
    tiltimg_file = list(data_redux.glob('Calibrations/Tiltimg_A_0_DET01.fits'))
    assert len(tiltimg_file) == 1, 'No tiltimg files found'
    tiltimg = TiltImage.from_file(tiltimg_file[0])
    assert tiltimg.spat_flexure is not None, 'Spat flexure should not be None'
    assert tiltimg.spat_flexure < par['calibrations']['tiltframe']['process']['spat_flexure_maxlag'], \
        'Spat flexure should be less than maxlag'
    assert np.isclose(tiltimg.spat_flexure, 6., atol=1.5), 'Spat flexure should be close to 6 pixels'
    # read the wavetilts file
    tilts_file = list(data_redux.glob('Calibrations/Tilts_A_0_DET01.fits'))
    assert len(tilts_file) == 1, 'No tilts files found'
    tilts = WaveTilts.from_file(tilts_file[0])
    assert tilts.spat_flexure is not None, 'Spat flexure should not be None'
    assert tilts.spat_flexure < par['calibrations']['tiltframe']['process']['spat_flexure_maxlag'], \
        'Spat flexure should be less than maxlag'
    assert np.isclose(tilts.spat_flexure, 6., atol=1.5), 'Spat flexure should be close to 6 pixels'

def test_flex_multi(redux_out):

    # Set output file
    outfile = data_output_path('tst_multi_flex.fits')
    if os.path.isfile(outfile):
        # Remove it if it already exists
        os.remove(outfile)

    spec1d_file = os.path.join(redux_out,
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

