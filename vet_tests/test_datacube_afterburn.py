"""
Vet tests for the pypeit_setup_datacube -> pypeit_coadd_datacube ->
pypeit_extract_datacube afterburn chain (see test_scripts/pypeit_tests.py and
test_scripts/test_setups.py's _setup_datacube/_coadd_datacube/_extract_datacube
registrations).
"""
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units
import numpy as np
import pytest

from pypeit import inputfiles, specobjs
from pypeit.coadd3d import DataCube


_DATASETS = [('keck_kcwi', 'large_bl'), ('keck_kcrm', 'large_rl')]
_TARGETS = ['SDSSJ2222 2745', 'gd50']


def _target_stub(target):
    return target.replace(' ', '')


@pytest.mark.parametrize('instr,setup', _DATASETS)
@pytest.mark.parametrize('target', _TARGETS)
def test_coadd3d_files_reference_correct_spec2d(redux_out, instr, setup, target):
    """The .coadd3d file written by pypeit_setup_datacube for a target must
    list exactly the spec2d files reduced for that target, no more, no less.
    """
    target_stub = _target_stub(target)
    setup_dir = Path(redux_out) / instr / setup
    coadd3d_file = setup_dir / 'sources' / target_stub / f'{target_stub}.coadd3d'
    assert coadd3d_file.is_file(), f'{coadd3d_file} does not exist'

    coadd3d = inputfiles.Coadd3DFile.from_file(str(coadd3d_file))
    listed_files = set(coadd3d.data['filename'])

    expected_files = {f.name for f in (setup_dir / 'Science').glob(f'spec2d_*{target_stub}*.fits')}
    assert listed_files == expected_files, (
        f'{coadd3d_file} lists {listed_files}, but the spec2d files reduced for '
        f'{target} are {expected_files}'
    )

    assert coadd3d.config['reduce']['cube']['output_filename'] == target_stub, \
        f'output_filename does not match the target stub {target_stub}'
    assert coadd3d.config['reduce']['cube']['combine'].lower() == 'true', \
        'combine should be True for a setup_datacube-generated .coadd3d file'
    assert coadd3d.config['reduce']['cube']['save_whitelight'].lower() == 'true', \
        'save_whitelight should be True for a setup_datacube-generated .coadd3d file'


@pytest.mark.parametrize('instr,setup', _DATASETS)
@pytest.mark.parametrize('target', _TARGETS)
def test_datacube_shape(redux_out, instr, setup, target):
    """The combined datacube must be a well-formed, non-degenerate 3D cube."""
    target_stub = _target_stub(target)
    cube_file = Path(redux_out) / instr / setup / 'Science_cube' / f'{target_stub}.fits'
    assert cube_file.is_file(), f'{cube_file} does not exist'

    cube = DataCube.from_file(str(cube_file))
    assert cube.flux.ndim == 3, f'flux is not 3D: {cube.flux.ndim}D'
    assert cube.sig.shape == cube.flux.shape, 'sig shape does not match flux shape'
    assert cube.bpm.shape == cube.flux.shape, 'bpm shape does not match flux shape'

    nwave, ny, nx = cube.flux.shape
    assert nwave > 1 and ny > 1 and nx > 1, f'degenerate cube shape: {cube.flux.shape}'
    assert cube.wave.size == nwave, 'wave array length does not match the spectral axis'

    good_frac = np.sum(cube.bpm == 0) / cube.bpm.size
    assert good_frac > 0.05, f'only {good_frac:.1%} of cube pixels are good'
    assert np.any(np.isfinite(cube.flux) & (cube.flux != 0)), \
        'flux cube is all zero or non-finite'


@pytest.mark.parametrize('instr,setup', _DATASETS)
@pytest.mark.parametrize('target', _TARGETS)
def test_whitelight_matches_cube(redux_out, instr, setup, target):
    """The whitelight image must match the cube's spatial shape and WCS."""
    target_stub = _target_stub(target)
    sci_cube_dir = Path(redux_out) / instr / setup / 'Science_cube'
    cube_file = sci_cube_dir / f'{target_stub}.fits'
    wl_file = sci_cube_dir / f'{target_stub}_whitelight.fits'
    assert cube_file.is_file(), f'{cube_file} does not exist'
    assert wl_file.is_file(), f'{wl_file} does not exist'

    cube_hdr = fits.getheader(cube_file, 'FLUX')
    wl_data, wl_hdr = fits.getdata(wl_file, 0, header=True)

    assert wl_data.ndim == 2, f'whitelight image is not 2D: {wl_data.ndim}D'
    expected_shape = (cube_hdr['NAXIS2'], cube_hdr['NAXIS1'])
    assert wl_data.shape == expected_shape, (
        f'whitelight shape {wl_data.shape} does not match the cube spatial shape '
        f'{expected_shape}'
    )
    assert np.any(np.isfinite(wl_data) & (wl_data != 0)), \
        'whitelight image is all zero or non-finite'

    cube_wcs = WCS(cube_hdr).celestial
    wl_wcs = WCS(wl_hdr)
    assert wl_wcs.has_celestial, 'whitelight header has no valid celestial WCS'
    assert np.allclose(cube_wcs.wcs.crval, wl_wcs.wcs.crval), \
        'whitelight and cube WCS CRVAL do not match'
    assert np.allclose(cube_wcs.wcs.cdelt, wl_wcs.wcs.cdelt), \
        'whitelight and cube WCS CDELT do not match'
    assert list(cube_wcs.wcs.ctype) == list(wl_wcs.wcs.ctype), \
        'whitelight and cube WCS CTYPE do not match'


@pytest.mark.parametrize('instr,setup', _DATASETS)
def test_gd50_extraction_finds_single_source(redux_out, instr, setup):
    """pypeit_extract_datacube on gd50 must find exactly one point source,
    with a non-degenerate extracted spectrum.
    """
    spec1d_file = Path(redux_out) / instr / setup / 'Science_cube' / 'spec1d_gd50_extract.fits'
    assert spec1d_file.is_file(), f'{spec1d_file} does not exist'

    sobjs = specobjs.SpecObjs.from_fitsfile(str(spec1d_file))
    assert len(sobjs) == 1, f'expected exactly 1 extracted source, found {len(sobjs)}'

    sobj = sobjs[0]
    # These test datasets are not flux-calibrated (no sensfunc supplied to
    # pypeit_setup_datacube), so the extraction is in counts, not flam. Pixels
    # outside the covered wavelength range have OPT_COUNTS_SIG == 0 (masked),
    # so a "good" pixel requires a positive, finite error, not just finite.
    good = np.isfinite(sobj.OPT_COUNTS) & np.isfinite(sobj.OPT_COUNTS_SIG) \
        & (sobj.OPT_COUNTS_SIG > 0)
    assert np.mean(good) > 0.5, \
        'fewer than half of the extracted spectrum pixels are good'


def test_gd50_coordinates_match_across_datasets(redux_out):
    """The gd50 point source must be extracted at the same sky position in
    both the KCWI and KCRM reductions of the same night's data.

    This is only a meaningful (non-circular) check because
    `extract_point_source`'s SpecObj.RA/DEC reflect the actual fitted or
    manually-specified extraction position, not the cube's WCS reference
    point (CRVAL).
    """
    coords = {}
    for instr, setup in _DATASETS:
        spec1d_file = Path(redux_out) / instr / setup / 'Science_cube' / 'spec1d_gd50_extract.fits'
        assert spec1d_file.is_file(), f'{spec1d_file} does not exist'
        sobjs = specobjs.SpecObjs.from_fitsfile(str(spec1d_file))
        assert len(sobjs) == 1, f'expected exactly 1 extracted source in {spec1d_file}'
        coords[instr] = SkyCoord(ra=sobjs[0].RA * units.deg, dec=sobjs[0].DEC * units.deg)

    separation = coords['keck_kcwi'].separation(coords['keck_kcrm'])
    assert separation < 2 * units.arcsec, (
        f'gd50 extracted position differs by {separation.arcsec:.2f} arcsec between '
        'keck_kcwi/large_bl and keck_kcrm/large_rl (empirically ~1.6 arcsec on real data)'
    )
