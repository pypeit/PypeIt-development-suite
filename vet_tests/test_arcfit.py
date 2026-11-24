"""
Vet tests of the arc-line fitting and the median-width floor of
:func:`pypeit.core.arc.detect_lines` (PR C, speed-up workstream C).

By default, ``detect_lines`` rejects lines whose Gaussian sigma is less than
``[calibrations][wavelengths] min_frac_med_width`` (0.5) times the median sigma
of the lines in the same spectrum, if at least ``min_nlines_med_width`` (5)
lines pass the other cuts.  These tests check, on the dev-suite wavelength
calibrations, that

  - no line used in a final wavelength fit is much narrower than the typical
    line of its spectrum;
  - Keck/DEIMOS 830G_L_8100, whose solution sub-resolution spikes degraded
    when a line fit that accepts them was tried without the floor (RMS
    0.27 px instead of 0.05 px, with one spike given the wavelength of the
    real line 1.7 px away), keeps a low RMS and no duplicate identification.
"""
import os
import re
from pathlib import Path

import numpy as np
import pytest

from pypeit.core import arc
from pypeit.par.pypeitpar import WavelengthSolutionPar
from pypeit.wavecalib import WaveCalib

# Instrument directories reduced for the PR C vetting; others are checked if present
INSTRUMENTS = ['shane_kast_blue', 'shane_kast_red', 'keck_deimos', 'keck_lris_blue']

# Tolerance on the floor: the lines are re-measured with the pipeline's main
# detection parameters, which reproduce its median sigma only approximately
FLOOR_TOL = 0.8


def _wavecalib_files(redux_out):
    files = []
    for inst in INSTRUMENTS:
        root = Path(redux_out) / inst
        if not root.is_dir():
            continue
        files += [f for f in sorted(root.rglob('Calibrations/WaveCalib_*.fits'))
                  if 'step_by_step' not in str(f) and '/ncpu/' not in str(f)]
    return files


def _wavelengths_par(wavecalib_file):
    """
    (fwhm, fwhm_fromlines, sigdetect, min_frac_med_width, min_nlines_med_width)
    of the reduction's [[wavelengths]] parameters.
    """
    dflt = WavelengthSolutionPar()
    pars = sorted(Path(wavecalib_file).parent.parent.glob('*.par'))
    if not pars:
        return 4.0, True, 5.0, dflt['min_frac_med_width'], dflt['min_nlines_med_width']
    txt = pars[-1].read_text()
    m = re.search(r'\[\[wavelengths\]\](.*?)(?=\n\s*\[\[)', txt, re.S)
    block = m.group(1) if m else txt
    get = lambda key: re.search(rf'^\s*{key}\s*=\s*(.+)$', block, re.M)
    fwhm = float(get('fwhm').group(1)) if get('fwhm') else 4.0
    froml = get('fwhm_fromlines').group(1).strip() == 'True' if get('fwhm_fromlines') else True
    sd = get('sigdetect')
    sigdetect = float(re.split(r'[,\s]+', sd.group(1).strip())[0]) if sd else 5.0
    frac = float(get('min_frac_med_width').group(1)) if get('min_frac_med_width') \
                else dflt['min_frac_med_width']
    nmin = int(get('min_nlines_med_width').group(1)) if get('min_nlines_med_width') \
                else dflt['min_nlines_med_width']
    return fwhm, froml, sigdetect, frac, nmin


def _final_lines(wv):
    gpm = wv.pypeitfit.gpm.astype(bool)
    return np.asarray(wv.pixel_fit)[gpm], np.asarray(wv.wave_fit)[gpm]


def test_arcfit_width_floor(redux_out):
    """No line in a final wavelength fit is far below the median-width floor."""
    files = _wavecalib_files(redux_out)
    if len(files) == 0:
        pytest.skip(f'No wavelength calibrations found in {redux_out}')
    bad, nchecked = [], 0
    for f in files:
        fwhm_par, fromlines, sigdetect, frac, nmin = _wavelengths_par(f)
        if frac <= 0:
            continue
        wc = WaveCalib.from_file(str(f), chk_version=False)
        for i, wv in enumerate(wc.wv_fits):
            if wv is None or wv.pypeitfit is None or wv.spec is None:
                continue
            fwhm = float(wv.fwhm) if fromlines and wv.fwhm is not None else fwhm_par
            # Re-measure the lines as the wavelength calibration does
            # (wvutils.arc_lines_from_spec), without the floor
            _, _, cent, wid, _, good, _, _ = arc.detect_lines(
                np.asarray(wv.spec), sigdetect=sigdetect, fwhm=fwhm, max_frac_fwhm=2.0,
                min_frac_med_width=0.)
            if good.size < nmin:
                continue
            floor = frac * np.median(wid[good])
            pix, wave = _final_lines(wv)
            for p, w in zip(pix, wave):
                d = np.absolute(cent[good] - p)
                if d.size == 0 or d.min() > 0.5:
                    continue
                nchecked += 1
                j = good[np.argmin(d)]
                if wid[j] < FLOOR_TOL * floor:
                    bad.append(f'{f.relative_to(redux_out)} slit {i}: pixel {p:.1f} ({w:.2f} A) '
                               f'sigma {wid[j]:.2f} < {FLOOR_TOL} x {floor:.2f}')
    assert nchecked > 1000, f'Only {nchecked} final-fit lines could be checked'
    assert len(bad) == 0, 'Lines narrower than the width floor in final fits:\n' + '\n'.join(bad)


def test_arcfit_deimos_830g_l_8100(redux_out):
    """Sub-resolution spikes must not degrade the 830G_L_8100 solution."""
    f = Path(redux_out) / 'keck_deimos' / '830G_L_8100' / 'Calibrations' \
            / 'WaveCalib_A_0_MSC03.fits'
    if not f.is_file():
        pytest.skip(f'{f} not found')
    wv = WaveCalib.from_file(str(f), chk_version=False).wv_fits[0]
    # 0.051 px in develop; 0.27 px with the (not adopted) vectorized fit and no floor
    assert wv.rms < 0.1, f'RMS of keck_deimos 830G_L_8100 is too high: {wv.rms:.3f}'
    pix, wave = _final_lines(wv)
    assert pix.size >= 65
    u, n = np.unique(np.round(wave, 3), return_counts=True)
    assert np.all(n == 1), f'Duplicate identifications: {u[n > 1]}'
