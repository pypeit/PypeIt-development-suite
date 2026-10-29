"""
For every arc spectrum stored in the WaveCalib files of a dev-suite run,
compare detect_lines with the legacy and the vectorized fit_arcspec, and
report the lines that only one of them accepts, by width relative to the
expected sigma (fwhm/2.355).
Usage: python arcfit_gained_lost.py REDUX_OUT_C/legacy
"""
import sys, warnings
from pathlib import Path
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.wavecalib import WaveCalib
from pypeit.core import arc

root = Path(sys.argv[1])
rat = {'legacy_good': [], 'gained': [], 'lost_new': [], 'lost_old': []}
peak = {'gained': [], 'lost': []}
nspec = 0
for f in sorted(root.rglob('Calibrations/WaveCalib_*.fits')):
    if '/ncpu/' in str(f) or 'step_by_step' in str(f):
        continue
    wc = WaveCalib.from_file(str(f), chk_version=False)
    for w in wc.wv_fits:
        if w is None or w.spec is None or w.fwhm is None:
            continue
        s, fw = np.asarray(w.spec), float(w.fwhm)
        res = {}
        for v in ['legacy', 'new']:
            arc.USE_VECTORIZED_ARCFIT = v == 'new'
            out = arc.detect_lines(s, sigdetect=5., fwhm=fw)
            res[v] = out
        arc.USE_VECTORIZED_ARCFIT = True
        nspec += 1
        sig0 = fw / 2.355
        go, gn = set(res['legacy'][5]), set(res['new'][5])
        rat['legacy_good'] += list(res['legacy'][3][sorted(go)] / sig0)
        for j in sorted(gn - go):
            rat['gained'].append(res['new'][3][j] / sig0); peak['gained'].append(res['new'][0][j])
        for j in sorted(go - gn):
            rat['lost_old'].append(res['legacy'][3][j] / sig0); rat['lost_new'].append(res['new'][3][j] / sig0)
            peak['lost'].append(res['legacy'][0][j])
print(f'{nspec} arc spectra')
lg = np.array(rat['legacy_good']); g = np.array(rat['gained']); lo = np.array(rat['lost_new'])
print(f'legacy-accepted lines: {lg.size}; sigma/sigma_expected percentiles 0.1/1/50: '
      + ' '.join(f'{np.percentile(lg, q):.2f}' for q in [0.1, 1, 50])
      + f'; fraction < 0.5: {np.mean(lg < 0.5):.4f}')
print(f'gained by new fit: {g.size}; width ratio < 0.5: {np.sum(g < 0.5)}; median ratio {np.median(g) if g.size else np.nan:.2f}')
print('  gained width ratios:', np.round(np.sort(g), 2)[:60])
print(f'lost by new fit: {lo.size}; new width ratio > fwhm_max cut (3*fwhm/2.35 -> ratio {3*2.355/2.35:.2f}): {np.sum(lo > 3*2.355/2.35)}; '
      f'median raw peak of lost {np.median(peak["lost"]) if peak["lost"] else np.nan:.0f}')
