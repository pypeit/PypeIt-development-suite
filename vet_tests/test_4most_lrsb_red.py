"""Check the science products from the seeded 4MOST commissioning reduction."""

import os
from pathlib import Path

from astropy.io import fits
import numpy as np
import pytest


def test_vista_4most_lrsb_red(redux_out):
    if os.environ.get('PYPEIT_4MOST_FIXTURE') != '1':
        pytest.skip('The 4MOST commissioning fixture is opt-in.')

    instrument = 'vista_4most_lrsb_red'
    setup = Path(redux_out) / instrument / 'OB5148227'
    science_dir = setup / 'Science'
    for name in ('Bias_A_0_DET01.fits', 'Arc_A_0_DET01.fits'):
        assert (setup / 'Calibrations' / name).is_file()

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
        assert np.all(spectrum['OPT_COUNTS_IVAR'][good] > 0)

    with fits.open(spec2d) as hdus:
        assert hdus[0].header['PYP_SPEC'] == instrument
        sky = hdus['DET01-SKYMODEL']
        assert sky.shape == (6144, 6160)
        assert np.any(np.isfinite(sky.section[3072, :]) & (sky.section[3072, :] != 0))
