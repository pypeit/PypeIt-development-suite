"""
Module to run tests on BiasFrame class
Requires files in Development suite and an Environmental variable
"""
from pathlib import Path
import os

from IPython import embed

import pytest
import glob
import numpy as np

from pypeit.pypmsgs import PypeItError
from pypeit.images import buildimage
from pypeit.images import pypeitimage
from pypeit.calibframe import CalibFrame
from pypeit.images.mosaic import Mosaic
from pypeit.tests.tstutils import data_output_path
from pypeit.spectrographs.util import load_spectrograph


# Init a few things
shane_kast_blue = load_spectrograph('shane_kast_blue')
frame_par = shane_kast_blue.default_pypeit_par()['calibrations']['biasframe']
setup = 'A'
calib_id = '1'
det = 1
calib_dir = data_output_path('')


@pytest.fixture
def kast_blue_bias_files():
    root = Path(os.getenv('PYPEIT_DEV')).resolve() / 'RAW_DATA' / 'shane_kast_blue' / '600_4310_d55'
    kast_blue_files = sorted(list(root.glob('b1?.fits*')))
    return [str(f) for f in kast_blue_files[5:]]


def test_process(kast_blue_bias_files):
    # Instantiate
    bias_img = buildimage.buildimage_fromlist(shane_kast_blue, 1, frame_par,
                                              kast_blue_bias_files)
    assert isinstance(bias_img, CalibFrame), 'Class hierarchy error, CalibFrame'
    assert isinstance(bias_img, pypeitimage.PypeItImage), 'Class hierarchy error, PypeItImage'
    assert isinstance(bias_img, pypeitimage.PypeItCalibrationImage), \
            'Class hierarchy error, PypeItCalibrationImage'
    assert isinstance(bias_img, buildimage.BiasImage), 'Class hierarchy error, BiasImage'
    assert isinstance(bias_img.image, np.ndarray)

    with pytest.raises(PypeItError):
        # Will not write it because you haven't defined the attributes necessary
        # to set the output path
        bias_img.to_file()


def test_io(kast_blue_bias_files):
    bias_img = buildimage.buildimage_fromlist(shane_kast_blue, 1, frame_par,
                                              kast_blue_bias_files, calib_dir=calib_dir,
                                              setup=setup, calib_id=calib_id)
    ofile = Path(bias_img.get_path()).resolve()
    bias_img.to_file(overwrite=True)
    assert ofile.exists(), 'Error writing Bias file'

    # Load processed calibration frame
    _bias_img = buildimage.BiasImage.from_file(str(ofile))
    assert np.array_equal(bias_img.image, _bias_img.image), 'Image changed'
    assert np.array_equal(bias_img.ivar, _bias_img.ivar), 'Inverse-variance changed'

    # Clean up
    ofile.unlink()


def test_process_multidet():
    files = glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'gemini_gmos',
                                   'GN_HAM_R400_885', 'N20190205S024*.fits'))
    files.sort()
    spec = load_spectrograph('gemini_gmos_north_ham')
    frame_par = spec.default_pypeit_par()['calibrations']['biasframe']

    det = 1
    bias_img_det1 = buildimage.buildimage_fromlist(spec, det, frame_par, files)

    det = (1,2,3)
    bias_img = buildimage.buildimage_fromlist(spec, det, frame_par, files)

    assert np.array_equal(bias_img_det1.image, bias_img.image[0]) \
            and not np.array_equal(bias_img_det1.image, bias_img.image[1]) \
            and not np.array_equal(bias_img_det1.image, bias_img.image[2]), \
                'Bad multi-detector processing'


def test_mosaic_io():
    files = glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'gemini_gmos',
                                   'GN_HAM_R400_885', 'N20190205S024*.fits'))
    files.sort()
    spec = load_spectrograph('gemini_gmos_north_ham')
    frame_par = spec.default_pypeit_par()['calibrations']['biasframe']

    det = (1,2,3)
    bias = buildimage.buildimage_fromlist(spec, det, frame_par, files, calib_dir=calib_dir,
                                          setup=setup, calib_id=calib_id)
    bias.to_file()
    outfile = bias.get_path()

    _bias = buildimage.BiasImage.from_file(outfile)
    assert isinstance(_bias.detector, Mosaic), 'Bad detector type'
    assert _bias.detector.ndet == bias.detector.ndet, 'Bad number of detectors'
    assert np.array_equal(_bias.detector.tform, bias.detector.tform), 'Bad number of detectors'

    # Clean up
    Path(outfile).unlink()


