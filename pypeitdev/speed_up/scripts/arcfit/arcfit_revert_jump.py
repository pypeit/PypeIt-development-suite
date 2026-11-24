"""
"Revert to CurveFit?" item 1: test a step (vertical-edge) cut on the legacy
curve_fit + floor (LF) line list.  For every LF line, jump = max |adjacent-pixel
difference| in the fit window of the continuum-subtracted arc / max of the
window; gjump = the same quantity for a noiseless sampled Gaussian of the
pipeline fwhm (worst case over sub-pixel phase).  Compares the 46 LF-not-B
lines with the lines used in final fits.  Reads arcfit_revert_survey.npz for
the B/LF flags.
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
R, names = S['rows'], list(S['names'])


def gauss_jump(fw):
    s = fw / 2.355
    best = 0.
    for ph in np.linspace(0, 1, 21):
        x = np.arange(-30, 31) + ph
        g = np.exp(-0.5 * (x / s)**2)
        best = max(best, np.max(np.abs(np.diff(g))) / g.max())
    return best


files = sorted(f for f in ROOT.rglob('Calibrations/WaveCalib_*.fits')
               if '/ncpu/' not in str(f) and 'step_by_step' not in str(f))
out = []
for fi_, f in enumerate(files):
    wc = WaveCalib.from_file(str(f), chk_version=False)
    for i, w in enumerate(wc.wv_fits):
        if w is None or w.spec is None or w.fwhm is None:
            continue
        sel = (R[:, 0] == fi_) & (R[:, 1] == i)
        if not sel.any():
            continue
        spec, fw = np.asarray(w.spec), pipeline_fwhm(f, w)
        arc.USE_VECTORIZED_ARCFIT = False
        a = arc.detect_lines(spec, sigdetect=5., fwhm=fw, min_frac_med_width=0.5)[6]
        arc.USE_VECTORIZED_ARCFIT = True
        nfitpix = int(np.round(1.25 * fw)); half = (nfitpix + nfitpix % 2) // 2
        gj = gauss_jump(fw)
        for r in R[sel]:
            p = int(r[2]); win = a[max(p - half, 0):p + half + 1]
            j = np.max(np.abs(np.diff(win))) / win.max() if win.max() > 0 else 0.
            out.append((j, gj))
    print(fi_ + 1, flush=True, end=' ')
J = np.array(out)
fi, sl, pix, gB, gL, gF, ratio, rmax, raw, infit, fw, cB, cL = R.T
gB, gF, infit = gB > 0, gF > 0, infit > 0
np.save('/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C/arcfit_revert_jump.npy', J)
jump, gj = J.T
onlyL = gF & ~gB
print('\nLF-not-B (n=%d): jump quantiles' % onlyL.sum(), np.round(np.percentile(jump[onlyL], [0, 10, 50, 90]), 3),
      ' jump/gj', np.round(np.percentile((jump / gj)[onlyL], [0, 10, 50]), 2))
fit = gF & infit
print('final-fit lines (n=%d): jump max %.3f, 99.9pct %.3f; jump/gj max %.2f 99.9pct %.2f' %
      (fit.sum(), jump[fit].max(), np.percentile(jump[fit], 99.9), (jump / gj)[fit].max(),
       np.percentile((jump / gj)[fit], 99.9)))
for k in [1.5, 2.0, 2.5, 3.0]:
    for rawcut in [0, 3e4]:
        cut = gF & (jump > k * gj) & (raw > rawcut)
        print(f'cut jump > {k} x gauss{" & raw>30k" if rawcut else ""}: removes {cut.sum()} LF lines, '
              f'{np.sum(cut & onlyL)} of the {onlyL.sum()} LF-not-B, {np.sum(cut & infit)} in final fits, '
              f'{np.sum(cut & gB)} accepted by B')
