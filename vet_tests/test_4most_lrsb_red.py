"""Check the science products from the seeded 4MOST commissioning reduction."""

import os
from pathlib import Path

from astropy.io import fits
import numpy as np
import pytest

from pypeit.calibframe import CalibFrame


def test_vista_4most_lrsb_red(redux_out):
    if os.environ.get('PYPEIT_4MOST_FIXTURE') != '1':
        pytest.skip('The 4MOST commissioning fixture is opt-in.')

    instrument = 'vista_4most_lrsb_red'
    setup = Path(redux_out) / instrument / 'OB5148227'
    science_dir = setup / 'Science'
    calib_dir = (setup / 'Calibrations').resolve()
    for name in ('Bias_A_0_DET01.fits', 'Arc_A_0_DET01.fits'):
        assert (calib_dir / name).is_file()
    _, slits_calib_dir = CalibFrame.parse_key_dir(str(calib_dir / 'Slits_A_0_DET01.fits.gz'))
    assert Path(slits_calib_dir) == calib_dir

    prefix = 'QMOST.2026-11-27T00:47:46.453-'
    spec1d, = science_dir.glob(f'spec1d_{prefix}*.fits')
    spec2d, = science_dir.glob(f'spec2d_{prefix}*.fits')

    with fits.open(spec1d) as hdus:
        assert hdus[0].header['PYP_SPEC'] == instrument
        assert hdus[0].header['PYPELINE'] == 'Fiber'
        fibers = [hdu for hdu in hdus if hdu.header.get('DMODCLS') == 'SpecObj']
        ids = {hdu.header['MASKDEF_ID'] for hdu in fibers}
        assert len(fibers) == len(ids) == 811
        assert 1052 not in ids  # Inactive commissioning fiber.
        assert sum(hdu.header.get('MASKDEF_OBJNAME') == 'SKY' for hdu in fibers) == 80

        bright, = (hdu for hdu in fibers if hdu.header['MASKDEF_ID'] == 2025)
        spectrum = bright.data
        assert spectrum['OPT_WAVE'].shape == (6144,)
        good = spectrum['OPT_MASK']
        wave = spectrum['OPT_WAVE'][good]
        assert np.all(np.diff(wave) > 0)
        assert np.count_nonzero((wave >= 8440) & (wave < 8860)) > 300
        assert np.count_nonzero((wave >= 9240) & (wave < 9300)) >= 50
        assert np.all(np.isfinite(spectrum['OPT_COUNTS'][good]))
        ivar = spectrum['OPT_COUNTS_IVAR'][good]
        assert np.all(np.isfinite(ivar) & (ivar > 0))

    with fits.open(spec2d) as hdus:
        assert hdus[0].header['PYP_SPEC'] == instrument
        sky = hdus['DET01-SKYMODEL']
        assert sky.shape == (6144, 6160)
        row = sky.section[3072, :]
        assert np.any(np.isfinite(row) & (row != 0))
