"""
Module to run tests on scripts
"""
from pathlib import Path
import os
import numpy as np
import pytest

from IPython import embed
from pypeit import coadd2d
from pypeit import inputfiles
from pypeit import specobjs
from pypeit.spectrographs.util import load_spectrograph


def test_offsets_and_weights(redux_out):

    # echelle data
    spec_name = 'keck_nires'
    _redux_out = Path(redux_out).absolute() / spec_name / 'ABBA_nostandard'
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
    par['coadd2d']['offsets'] = 'auto'
    par['coadd2d']['weights'] = 'auto'
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par)

    # check the offsets
    assert coadd.offsets is None, 'Wrong offsets'
    # check the weights
    assert len(coadd.use_weights) == 5, 'The size of the weights should be equal to the number of orders'

    # init coadd2d with offsets=list and weights=list
    par['coadd2d']['offsets'] = [0.5, 0.5, 0.5, 0.5]
    par['coadd2d']['weights'] = [1., 1., 1., 1.]
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par)

    # check the offsets
    assert np.all(coadd.offsets == [0.5, 0.5, 0.5, 0.5]), 'Wrong offsets'
    assert isinstance(coadd.offsets, np.ndarray), 'Offsets should be a numpy array'
    # check the weights
    assert np.all(coadd.use_weights == [1.0, 1.0, 1.0, 1.0]), 'Wrong weights'
    assert isinstance(coadd.use_weights, list), 'Weights should be a list'

    # init coadd2d with offsets='header' and weights='uniform'
    par['coadd2d']['offsets'] = 'header'
    par['coadd2d']['weights'] = 'uniform'
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par)

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
    _redux_out = Path(redux_out).absolute() / spec_name / 'long2pos1_H'
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
    par['coadd2d']['offsets'] = 'maskdef_offsets'
    par['coadd2d']['weights'] = 'auto'
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par)

    # check the offsets
    assert np.all(np.round(coadd.offsets) == [0., 78.]), 'Wrong offsets'
    # check the weights
    assert len(coadd.use_weights) == len(spec2d_files), 'Wrong size for the weights'

    # init coadd2d with offsets='auto' and weights=list
    par['coadd2d']['offsets'] = 'auto'
    par['coadd2d']['weights'] = [1., 1.]
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par)

    # check the offsets
    assert np.all(np.round(coadd.offsets) == [0., 77.]), 'Wrong offsets'
    # check the weights
    assert np.all(coadd.use_weights == [1.0, 1.0]), 'Wrong weights'
    assert isinstance(coadd.use_weights, list), 'weights should be a list'
    # Go back
    os.chdir(cdir)

def test_offsets_with_user_obj_ids(redux_out):

    # LDT/DeVeny DV1 data
    spec_name = 'ldt_deveny'
    _redux_out = Path(redux_out).absolute() / spec_name / 'DV1'
    sci_dir = _redux_out / 'Science'

    cdir = os.getcwd()
    os.chdir(_redux_out)
    # Load spectrograph
    spectrograph = load_spectrograph(spec_name)
    par = spectrograph.default_pypeit_par()
    # Grab the spec2d files
    spec2d_files = [str(f) for f in sorted(sci_dir.glob('spec2d*.fits'))]

    # check the number of exposures
    assert len(spec2d_files) == 3, 'Wrong number of exposures'

    # Init coadd2d with offsets='auto', weights='auto', and user_obj_ids
    par['coadd2d']['offsets'] = 'auto'
    par['coadd2d']['weights'] = 'auto'
    par['coadd2d']['user_obj_ids'] = [244, 243, 241]
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par)

    # Check the offsets
    assert np.allclose(coadd.offsets, [0.0, 1.29, 3.91], atol=0.1), 'Wrong offsets'

    # Check that, if no `user_obj_ids` are provided, other bright sources are detected
    par['coadd2d']['user_obj_ids'] = None
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par)
    assert not np.allclose(coadd.offsets, [0.0, 1.29, 3.91], atol=0.02), 'Found dim target'

def test_spat_spec_fract(redux_out):

    # check multislit with slitmask data - gmos data
    _redux_out = Path(redux_out).absolute() / 'gemini_gmos' / 'GS_HAM_B600_MOS'

    # TODO: I think this may be the only time that we had tried to grab the
    # parent of the `redux_out` directory, assuming it is the same as the
    # top-level directory dev-suite directory.
#    coadd2dfname = Path(redux_out).absolute().parent / 'coadd2d_files'/ 'gemini_gmos_gs_ham_b600_mos.coadd2d'
    # I changed it to the below instead, mimicking how the top-level directory
    # is defined in the unit tests.  We may want to define a `redux_in` or
    # `dev_root` fixture for this.
    coadd2dfname = (
        Path(os.getenv('PYPEIT_DEV')).absolute() / 'coadd2d_files'
        / 'gemini_gmos_gs_ham_b600_mos.coadd2d'
    )
    coadd_dir = _redux_out / 'Science_coadd'
    sci_dir = _redux_out / 'Science'

    coadd2dFile = inputfiles.Coadd2DFile.from_file(str(coadd2dfname))
    # update coadd2dFile files paths to point to the correct location
    coadd2dFile.file_paths = [Path(sci_dir).absolute()]
    spectrograph, par, _ = coadd2dFile.get_pypeitpar(pypeit_fits=True)

    # check the spatial and spectral fraction resampling parameters
    assert par['coadd2d']['spat_samp_fact'] == 1.2, 'Wrong spatial sampling factor'
    assert par['coadd2d']['spec_samp_fact'] == 1.2, 'Wrong spectral sampling factor'

    # check the boxcar parameters
    assert par['reduce']['extraction']['boxcar_radius'] == 1.5, 'Wrong boxcar radius in reduce parameters'
    assert par['reduce']['slitmask']['missing_objs_boxcar_rad'] == 1., 'Wrong boxcar radius for missing objects in slitmask reductions'

    # grab the spec1d files and objects
    # non coadded
    spec1d_file = sorted(sci_dir.glob('spec1d*.fits'))[0]
    sobj = specobjs.SpecObjs.from_fitsfile(spec1d_file)
    # coadded
    spec1d_coadd_file = sorted(coadd_dir.glob('spec1d*.fits'))[0]
    sobj_coadd = specobjs.SpecObjs.from_fitsfile(spec1d_coadd_file)

    # check that the boxcar radius (in arcsec) used for the coadd2d is the same as
    # the one in the parameters after applying the spatial sampling factor
    for s in sobj_coadd:
        if s.MASKDEF_EXTRACT:
            assert s.BOX_R_ASEC == par['reduce']['slitmask']['missing_objs_boxcar_rad'], \
                f'Wrong boxcar radius for missing objects in slitmask reductions for slit {s.SLITID}'
        else:
            assert s.BOX_R_ASEC == par['reduce']['extraction']['boxcar_radius'], \
                f'Wrong boxcar radius in reduce parameters for slit {s.SLITID}'

    # check that the platescale of the coadd2d data is equal to
    # the platescale of the non coadded data * the spatial sampling factor
    for s in sobj_coadd:
        assert np.isclose(s.DETECTOR.platescale, sobj[0].DETECTOR.platescale * par['coadd2d']['spat_samp_fact']), \
            f'Wrong platescale for slit {s.SLITID} in coadd2d mosaic container'
        for d in s.DETECTOR.detectors:
            assert d.platescale == s.DETECTOR.platescale, \
                f'Wrong platescale for slit {s.SLITID} in coadd2d single detector container'

    # check the spectral sampling factor by comparing the wavelength arrays of the coadded and non coadded data
    assert np.isclose(sobj.OPT_WAVE.shape[1]/sobj_coadd.OPT_WAVE.shape[1], par['coadd2d']['spec_samp_fact'], atol=0.01), \
        'Wrong spectral sampling factor applied to the coadd2d data'


def test_default_basename_single_spec2d(redux_out):
    # SOAR Goodman Blue M1 produces a single spec2d file.  This exercises
    # CoAdd2D.default_basename's extension stripping for a raw '.fits.fz'
    # file (see spectrographs/soar_goodman.py) in the single-file
    # (first == last) case.
    spec_name = 'soar_goodman_blue'
    _redux_out = Path(redux_out).absolute() / spec_name / 'M1'
    sci_dir = _redux_out / 'Science'

    spec2d_files = sorted(str(f) for f in sci_dir.glob('spec2d*.fits'))
    assert len(spec2d_files) == 1, 'Expected a single spec2d file for soar_goodman_blue/M1'

    basename = coadd2d.CoAdd2D.default_basename(spec2d_files)

    assert basename == ('0141_BB034315_m660205_10-02-2023-'
                        '0141_BB034315_m660205_10-02-2023-BB034315_m660205'), \
        'Wrong basename for soar_goodman_blue/M1'
    # The raw '.fits.fz' extension should be fully stripped, not left as '.fits'
    assert '.fits' not in basename and '.fz' not in basename, \
        'Raw extension not fully stripped from basename'


def test_default_basename_multi_spec2d(redux_out):
    # LDT DeVeny DV1 produces three spec2d files.  This exercises
    # CoAdd2D.default_basename's first/last file handling in the
    # multi-file (first != last) case.
    spec_name = 'ldt_deveny'
    _redux_out = Path(redux_out).absolute() / spec_name / 'DV1'
    sci_dir = _redux_out / 'Science'

    spec2d_files = sorted(str(f) for f in sci_dir.glob('spec2d*.fits'))
    assert len(spec2d_files) == 3, 'Expected three spec2d files for ldt_deveny/DV1'

    basename = coadd2d.CoAdd2D.default_basename(spec2d_files)

    assert basename == '20251027.0210-20251027.0212-2009HC', 'Wrong basename for ldt_deveny/DV1'
