"""
Task 10 (C.md): what would a lower bound on the fitted sigma do?

For every arc spectrum stored in the legacy dev-suite WaveCalib files
(REDUX_OUT_C/legacy, 38 files), re-run detect_lines with the fwhm the pipeline
used (the slit's measured fwhm, or the parameter when fwhm_fromlines is False), capture the fit_arcspec inputs, fit every line under A (shipped) and B
(ARCFIT_MATCH_CURVEFIT=False), apply detect_lines' cuts, and then evaluate
candidate sigma floors:
  abs085 : sigma > 0.85 px (FWHM > 2 px)
  rel050 : sigma > 0.50 * fwhm/2.355
  rule3  : sigma > ARCFIT_MIN_WIDTH * window half-width (the A hand-off threshold, used as a cut)
For each floor: lines it removes among the A-accepted (= legacy) lines, among
the B-only lines, and among the lines used in the stored final wavelength fit
(pixel_fit matched within 0.5 px).  Writes a per-line table (npz) and prints
per-setup summaries.
Usage: python arcfit_sigma_floor.py [out.npz]
"""
import sys, warnings
from pathlib import Path
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.core import arc
from pypeit.wavecalib import WaveCalib
from arcfit_pipeline_fwhm import pipeline_fwhm

ROOT = Path('/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C/legacy')
OUT = sys.argv[1] if len(sys.argv) > 1 else '/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C/arcfit_sigma_floor.npz'


def fit_slit(spec, fwhm):
    calls, orig = [], arc.fit_arcspec
    arc.fit_arcspec = lambda *a: (calls.append(a), orig(*a))[1]
    try:
        arc.detect_lines(spec, sigdetect=5., fwhm=fwhm)
    finally:
        arc.fit_arcspec = orig
    x, y, pixt, fitp = calls[0]
    res, good = {}, {}
    for v, match in [('A', True), ('B', False)]:
        arc.ARCFIT_MATCH_CURVEFIT = match
        res[v] = orig(x, y, pixt, fitp)
        r = res[v]
        good[v] = np.isfinite(r[2]) & (r[2] > 0) & (r[2] < 3*fwhm/2.35) & (r[1] > 0) & (r[1] < x[-1]) \
                  & (np.abs(r[1] - pixt) < 0.75*fwhm)
    arc.ARCFIT_MATCH_CURVEFIT = True
    return pixt, res, good, (fitp + fitp % 2)//2


rows, names = [], []
files = sorted(f for f in ROOT.rglob('Calibrations/WaveCalib_*.fits') if '/ncpu/' not in str(f) and 'step_by_step' not in str(f))
for fi_, f in enumerate(files):
    setup = str(f.relative_to(ROOT)).split('/Calibrations/')[0] + ':' + f.name.replace('WaveCalib_', '').replace('.fits', '')
    wc = WaveCalib.from_file(str(f), chk_version=False)
    for i, w in enumerate(wc.wv_fits):
        if w is None or w.spec is None or w.fwhm is None:
            continue
        # the fwhm the pipeline used (measured, or the parameter when fwhm_fromlines is False)
        spec, fw = np.asarray(w.spec), pipeline_fwhm(f, w)
        pixt, res, good, half = fit_slit(spec, fw)
        infit = np.zeros(pixt.size, bool)
        if w.pypeitfit is not None:
            pf = np.asarray(w.pixel_fit)[w.pypeitfit.gpm.astype(bool)]
            for j in range(pixt.size):
                c = res['A'][1][j] if good['A'][j] else res['B'][1][j]
                infit[j] = pf.size > 0 and np.min(np.abs(pf - c)) < 0.5
        for j in range(pixt.size):
            if not (good['A'][j] or good['B'][j]):
                continue
            rows.append((fi_, i, pixt[j], res['A'][2][j], res['B'][2][j], good['A'][j], good['B'][j], fw, half,
                         spec[pixt[j]], infit[j]))
    names.append(setup)
    print(f'{fi_+1}/{len(files)} {setup}: {len(rows)} lines so far', flush=True)
R = np.array(rows, dtype=float)
np.savez(OUT, rows=R, names=np.array(names))
print('saved', OUT)
