"""
Module to run tests on scripts
"""
from pathlib import Path
import os
import numpy as np
import pytest

from IPython import embed
from pypeit import coadd2d
from pypeit.spectrographs.util import load_spectrograph


def test_offsets_and_weights(redux_out):

    # echelle data
    spec_name = 'keck_nires'
    _redux_out = Path(redux_out).resolve() / spec_name / 'ABBA_nostandard'
    sci_dir = _redux_out / 'Science'

    cdir = os.getcwd()
    os.chdir(_redux_out)
    # load spectrograph
    spectrograph = load_spectrograph(spec_name)
    par = spectrograph.default_pypeit_par()
    # grab the spec2d files
    spec2d_files = [str(f) for f in sorted(sci_dir.glob('spec2d*.fits'))]

    # check the number of exposures
    assert len(spec2d_files) == 4, 'Wrong number of exposures'

    # init coadd2d with offsets='auto' and weights='auto'
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par, offsets='auto', weights='auto')

    # check the offsets
    assert coadd.offsets is None, 'Wrong offsets'
    # check the weights
    assert len(coadd.use_weights) == 5, 'The size of the weights should be equal to the number of orders'

    # init coadd2d with offsets=list and weights=list
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par,
                                         offsets=[0.5, 0.5, 0.5, 0.5], weights=[1., 1., 1., 1.])

    # check the offsets
    assert np.all(coadd.offsets == [0.5, 0.5, 0.5, 0.5]), 'Wrong offsets'
    assert isinstance(coadd.offsets, np.ndarray), 'Offsets should be a numpy array'
    # check the weights
    assert np.all(coadd.use_weights == [1.0, 1.0, 1.0, 1.0]), 'Wrong weights'
    assert isinstance(coadd.use_weights, list), 'Weights should be a list'

    # init coadd2d with offsets='header' and weights='uniform'
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par, offsets='header', weights='uniform')

    # check the offsets
    assert np.all(np.round(coadd.offsets) == [0., -67., 0., -67.]), 'Wrong offsets'
    assert isinstance(coadd.offsets, np.ndarray), 'Offsets should be a numpy array'
    # check the weights
    assert np.all(coadd.use_weights == [0.25, 0.25, 0.25, 0.25]), 'Wrong weights'
    assert isinstance(coadd.use_weights, list), 'Weights should be a list'

    # Go back
    os.chdir(cdir)

    # check multislit data
    spec_name = 'keck_mosfire'
    _redux_out = Path(redux_out).resolve() / spec_name / 'long2pos1_H'
    sci_dir = _redux_out / 'Science'

    os.chdir(_redux_out)
    # load spectrograph
    spectrograph = load_spectrograph(spec_name)
    par = spectrograph.default_pypeit_par()
    # grab the spec2d files
    spec2d_files = [str(f) for f in sorted(sci_dir.glob('spec2d*.fits'))]

    # check the number of exposures
    assert len(spec2d_files) == 2, 'Wrong number of exposures'

    # init coadd2d with offsets='maskdef_offsets' and weights='auto'
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par, offsets='maskdef_offsets', weights='auto')

    # check the offsets
    assert np.all(np.round(coadd.offsets) == [0., 78.]), 'Wrong offsets'
    # check the weights
    assert len(coadd.use_weights) == len(spec2d_files), 'Wrong size for the weights'

    # init coadd2d with offsets='auto' and weights=list
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par, offsets='auto', weights=[1., 1.])

    # check the offsets
    assert np.all(np.round(coadd.offsets) == [0., 77.]), 'Wrong offsets'
    # check the weights
    assert np.all(coadd.use_weights == [1.0, 1.0]), 'Wrong weights'
    assert isinstance(coadd.use_weights, list), 'weights should be a list'
    # Go back
    os.chdir(cdir)




