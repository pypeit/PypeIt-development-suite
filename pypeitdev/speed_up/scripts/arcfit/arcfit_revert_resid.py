"""
"Revert to CurveFit?" item 1: test a goodness-of-fit cut on the legacy
curve_fit + floor (LF) line list: resid = RMS(window - fitted Gaussian) /
fitted amplitude, over the same fit window curve_fit used.  Compares the 46
LF-not-B lines (saturated plateau edges) with the lines used in final fits.
"""
import sys, warnings
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
S = np.load('/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C/arcfit_revert_survey.npz')
R = S['rows']
files = sorted(f for f in ROOT.rglob('Calibrations/WaveCalib_*.fits')
               if '/ncpu/' not in str(f) and 'step_by_step' not in str(f))
out = []
for fi_, f in enumerate(files):
    wc = WaveCalib.from_file(str(f), chk_version=False)
    for i, w in enumerate(wc.wv_fits):
        sel = (R[:, 0] == fi_) & (R[:, 1] == i)
        if w is None or w.spec is None or w.fwhm is None or not sel.any():
            continue
        spec, fw = np.asarray(w.spec), pipeline_fwhm(f, w)
        cap, orig = [], arc.fit_arcspec
        arc.fit_arcspec = lambda *a: (lambda r: (cap.append((a, r)), r)[1])(orig(*a))
        arc.USE_VECTORIZED_ARCFIT = False
        try:
            arc.detect_lines(spec, sigdetect=5., fwhm=fw, min_frac_med_width=0.5)
        finally:
            arc.fit_arcspec = orig
            arc.USE_VECTORIZED_ARCFIT = True
        (x, y, pixt, fitp), (amp, cen, sig, _) = cap[0]
        half = (fitp + fitp % 2) // 2
        idx = {int(p): k for k, p in enumerate(pixt)}
        for r in R[sel]:
            k = idx[int(r[2])]; p = int(pixt[k])
            lo, hi = max(p - half, 0), min(p + half + 1, y.size)
            xx, yy = x[lo:hi], y[lo:hi]
            model = amp[k] * np.exp(-0.5 * ((xx - cen[k]) / sig[k])**2)
            out.append(np.sqrt(np.mean((yy - model)**2)) / abs(amp[k]) if amp[k] != 0 else np.inf)
    print(fi_ + 1, end=' ', flush=True)
res = np.array(out)
np.save('/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C/arcfit_revert_resid.npy', res)
fi, sl, pix, gB, gL, gF, ratio, rmax, raw, infit, fw, cB, cL = R.T
gB, gF, infit = gB > 0, gF > 0, infit > 0
onlyL = gF & ~gB
print('\nLF-not-B (n=%d) resid quantiles' % onlyL.sum(), np.round(np.percentile(res[onlyL], [0, 10, 50, 90, 100]), 3))
fit = gF & infit
print('final-fit lines (n=%d) resid quantiles' % fit.sum(), np.round(np.percentile(res[fit], [50, 99, 99.9, 100]), 3))
for t in [0.05, 0.08, 0.1, 0.15, 0.2]:
    cut = gF & (res > t)
    print(f'cut resid > {t}: removes {cut.sum()} LF lines, {np.sum(cut & onlyL)}/{onlyL.sum()} LF-not-B, '
          f'{np.sum(cut & infit)} in final fits, {np.sum(cut & gB)} accepted by B')
