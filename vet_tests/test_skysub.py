"""
Module to run test user-defined sky-subtraction regions
"""
from pathlib import Path
import os
import shutil

from IPython import embed

import numpy as np

from pypeit.spectrographs.util import load_spectrograph
from pypeit import spec2dobj
from pypeit.calibframe import CalibFrame
from pypeit.images.buildimage import SkyRegions
from pypeit import io
from pypeit import pypeit
from pypeit.core import skysub


def test_skysub(redux_out):

    redux_path = Path(redux_out).resolve() / 'shane_kast_blue' / '600_4310_d55' \
                    / 'shane_kast_blue_A'

    # Load the spec2d file
    spec2d = redux_path / 'Science' / 'spec2d_b24-Feige66_KASTb_20150520T041246.960.fits'
    assert spec2d.exists(), 'spec2D file does not exist'
    spec = load_spectrograph('shane_kast_blue')
    spec2DObj = spec2dobj.Spec2DObj.from_file(spec2d, spec.get_det_name(1), chk_version=True)

    # Set the name for the SkyRegions file
    calib_key, _ = CalibFrame.parse_key_dir(spec2DObj.calibs['EDGES'], from_filename=True)
    regfile = SkyRegions.construct_file_name(calib_key, calib_dir=spec2DObj.head0['CALIBDIR'],
                                             basename=io.remove_suffix(spec2DObj.head0['FILENAME']))
    regfile = Path(regfile).resolve()

    # If it exists, remove it
    if regfile.exists():
        regfile.unlink()

    # Make the mask regions
    region = ':40,60:'
    status, regs = skysub.read_userregions(region, spec2DObj.slits.nslits,
                maxslitlength=np.max(spec2DObj.slits.get_slitlengths(initial=True)))
    assert status == 0, 'Bad user region status'
    # And save them to the SkyRegions file
    left, right, _ = spec2DObj.slits.select_edges(initial=True)
    mask = skysub.generate_mask('MultiSlit', regs, spec2DObj.slits, left, right)
    SkyRegions(image=mask.astype(float), PYP_SPEC='shane_kast_blue').to_file(file_path=regfile)
    assert regfile.exists(), 'Regions not written'

    # Try to re-reduce the standard using the "user-defined" sky-regions

    # Set the pypeit file
    pypeit_file = Path(redux_out).resolve() / 'shane_kast_blue' / '600_4310_d55' \
                    / 'shane_kast_blue_A' / 'shane_kast_blue_A.pypeit'
    assert pypeit_file.exists(), 'PypeIt file missing!'

    # Initialize the main pypeit run
    pypeIt = pypeit.PypeIt(str(pypeit_file), verbosity=2, reuse_calibs=True, overwrite=True,
                           redux_path=str(redux_path))
    assert pypeIt.fitstbl.n_calib_groups == 1, 'Number of calibration groups changed'
    is_standard = pypeIt.fitstbl.find_frames('standard')
    assert np.sum(is_standard) == 1, 'Number of standard frames changed'
    std_frame = np.where(is_standard)[0]

    # Use the SkyRegions file
    pypeIt.par['reduce']['skysub']['user_regions'] = 'user'
    # This should *not* overwrite any existing spec2d or spec1d file
    new_spec2d, new_spec1d = pypeIt.reduce_exposure(std_frame)

    # TODO: Would be nice to have tests that can tell whether or not the code
    # actually used the SkyRegions file...
    assert len(new_spec1d) == 1, 'Should extract 1 spectrum'

    # And try again using a directly defined region
    pypeIt.par['reduce']['skysub']['user_regions'] = region
    # This should *not* overwrite any existing spec2d or spec1d file
    _new_spec2d, _new_spec1d = pypeIt.reduce_exposure(std_frame)
    assert len(_new_spec1d) == 1, 'Should extract 1 spectrum'

    # Result should be identical
    # NOTE: This used to use np.array_equal, but this caused an error in the
    # cloud.
    assert np.allclose(new_spec1d[0].BOX_WAVE, _new_spec1d[0].BOX_WAVE), \
            'wavelength should be the same'
    assert np.allclose(new_spec1d[0].BOX_COUNTS, _new_spec1d[0].BOX_COUNTS), \
            'extracted counts should be the same'



