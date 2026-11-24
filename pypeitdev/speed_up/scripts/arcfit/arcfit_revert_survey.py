"""
"Revert to CurveFit?" item 1 (C.md): can the legacy per-line curve_fit, plus
explicit cuts, reproduce the line lists of option B (the PR default)?

For every arc spectrum stored in the legacy dev-suite WaveCalib files
(REDUX_OUT_C/legacy, 38 files, 370 spectra), run detect_lines with the fwhm
the pipeline used, under
  B   : PR default (vectorized fit, no hand-off, median-width floor 0.5)
  L   : legacy curve_fit loop (USE_VECTORIZED_ARCFIT=False), no floor
  LF  : L + median-width floor 0.5
  LFP : LF + a plateau cut: reject a line if, in its fit window of the
        continuum-subtracted arc, the 3rd-highest pixel exceeds
        PLAT_FRAC x the maximum.  For a sampled Gaussian of the pipeline
        fwhm that ratio is at most exp(-1/(2 sigma^2)) (0.84 for fwhm 4,
        0.95 for fwhm 7.5); plateau windows give ~0.999.
Lines are matched between variants by peak pixel (pixt).  Also records
whether each line is used in the stored final wavelength fit (pixel_fit
within 0.5 px) and its raw peak, and times the fits.

Usage: python arcfit_revert_survey.py [out.npz]
"""
import sys, time, warnings
from pathlib import Path
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.core import arc
from pypeit.wavecalib import WaveCalib
sys.path.insert(0, str(Path(__file__).parent))
from arcfit_pipeline_fwhm import pipeline_fwhm

ROOT = Path('/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C/legacy')
OUT = sys.argv[1] if len(sys.argv) > 1 else \
    '/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C/arcfit_revert_survey.npz'
PLAT_FRAC = 0.98


def run(spec, fwhm, vec, floor):
    calls, orig = [], arc.fit_arcspec
    arc.fit_arcspec = lambda *a: (calls.append(a[2]), orig(*a))[1]
    arc.USE_VECTORIZED_ARCFIT = vec
    t0 = time.perf_counter()
    try:
        out = arc.detect_lines(spec, sigdetect=5., fwhm=fwhm, min_frac_med_width=floor)
    finally:
        arc.fit_arcspec = orig
        arc.USE_VECTORIZED_ARCFIT = True
    dt = time.perf_counter() - t0
    tcent, ww, arcsub = out[2], out[5], out[6]
    pixt = calls[0] if calls else np.array([], int)
    good = np.zeros(pixt.size, bool)
    good[ww] = True
    return pixt, good, tcent, arcsub, dt


def plateau_ratio(arcsub, pixt, half):
    r = np.zeros(pixt.size)
    for j, p in enumerate(pixt):
        w = np.sort(arcsub[max(p - half, 0):p + half + 1])
        r[j] = w[-3] / w[-1] if w.size >= 3 and w[-1] > 0 else 0.
    return r


rows, names, times = [], [], np.zeros(3)
files = sorted(f for f in ROOT.rglob('Calibrations/WaveCalib_*.fits')
               if '/ncpu/' not in str(f) and 'step_by_step' not in str(f))
for fi_, f in enumerate(files):
    setup = str(f.relative_to(ROOT)).split('/Calibrations/')[0] + ':' + f.name.replace('WaveCalib_', '').replace('.fits', '')
    wc = WaveCalib.from_file(str(f), chk_version=False)
    for i, w in enumerate(wc.wv_fits):
        if w is None or w.spec is None or w.fwhm is None:
            continue
        spec, fw = np.asarray(w.spec), pipeline_fwhm(f, w)
        pB, gB, cB, aB, tB = run(spec, fw, True, None)
        pL, gL, cL, aL, tL = run(spec, fw, False, 0.)
        pF, gF, cF, aF, tF = run(spec, fw, False, 0.5)
        times += [tB, tL, tF]
        assert np.array_equal(pB, pL) and np.array_equal(pB, pF)
        nfitpix = int(np.round(1.25 * fw))
        half = (nfitpix + nfitpix % 2) // 2
        ratio = plateau_ratio(aL, pL, half)
        sig = fw / 2.355
        pf = np.asarray(w.pixel_fit)[w.pypeitfit.gpm.astype(bool)] if w.pypeitfit is not None else np.array([])
        for j in range(pL.size):
            if not (gB[j] or gL[j]):
                continue
            c = cL[j] if gL[j] else cB[j]
            infit = pf.size > 0 and np.min(np.abs(pf - c)) < 0.5
            rows.append((fi_, i, pL[j], gB[j], gL[j], gF[j], ratio[j], np.exp(-0.5 / sig**2),
                         spec[pL[j]], infit, fw, cB[j], cL[j]))
    names.append(setup)
    print(f'{fi_+1}/{len(files)} {setup}: {len(rows)} lines; t(B,L,LF) = {np.round(times, 1)}', flush=True)
R = np.array(rows, dtype=float)
np.savez(OUT, rows=R, names=np.array(names), times=times)
print('saved', OUT)
