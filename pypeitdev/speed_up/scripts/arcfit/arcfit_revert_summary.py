"""
Summarise arcfit_revert_survey.py: line-list differences between option B and
the legacy curve_fit loop with the median-width floor (LF) and with the floor
plus a plateau cut (LFP), for several plateau thresholds.
Usage: python arcfit_revert_summary.py [survey.npz]
"""
import sys
import numpy as np

fn = sys.argv[1] if len(sys.argv) > 1 else \
    '/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C/arcfit_revert_survey.npz'
d = np.load(fn)
R, names, times = d['rows'], d['names'], d['times']
fi, sl, pix, gB, gL, gF, ratio, rmax, raw, infit, fw, cB, cL = R.T
gB, gL, gF, infit = gB > 0, gL > 0, gF > 0, infit > 0
nspec = len({(a, b) for a, b in zip(fi, sl)})
print(f'{nspec} spectra, {len(names)} files; detect_lines time B {times[0]:.1f} s, '
      f'L {times[1]:.1f} s, LF {times[2]:.1f} s')
print(f'accepted: B {gB.sum()}, L {gL.sum()}, LF {gF.sum()}')
both = gB & gF
print(f'shared B & LF: {both.sum()}; centroid |dc| median {np.median(np.abs(cB-cL)[both]):.2e} px, '
      f'99th pct {np.percentile(np.abs(cB-cL)[both], 99):.2e}, max {np.max(np.abs(cB-cL)[both]):.3f}')


def report(label, g):
    onlyL, onlyB = g & ~gB, gB & ~g
    print(f'\n{label}: accepted {g.sum()}; in it but not B {onlyL.sum()} '
          f'(in final fit {np.sum(onlyL & infit)}, raw>30k {np.sum(onlyL & (raw > 3e4))}); '
          f'in B but not it {onlyB.sum()} (in final fit {np.sum(onlyB & infit)})')
    return onlyL, onlyB


onlyL, onlyB = report('LF (curve_fit + floor)', gF)
print('  LF-not-B plateau ratio quantiles:', np.round(np.percentile(ratio[onlyL], [0, 10, 50, 90, 100]), 4)
      if onlyL.any() else '')
print('  final-fit lines (LF) plateau ratio: max', np.round(ratio[gF & infit].max(), 4),
      ' 99.9 pct', np.round(np.percentile(ratio[gF & infit], 99.9), 4),
      ' max ratio/gauss-bound', np.round(np.max(ratio[gF & infit] / rmax[gF & infit]), 4))
for thr in [0.97, 0.98, 0.99, 0.995]:
    g = gF & (ratio <= thr)
    oL, oB = report(f'LFP thr={thr}', g)
    cut = gF & (ratio > thr)
    print(f'  plateau cut removes {cut.sum()} LF lines, {np.sum(cut & infit)} of them in a final fit, '
          f'{np.sum(cut & gB)} of them accepted by B')
    if np.sum(cut & infit):
        for j in np.where(cut & infit)[0][:10]:
            print(f'    {names[int(fi[j])]} slit {int(sl[j])} pix {int(pix[j])} ratio {ratio[j]:.4f} '
                  f'raw {raw[j]:.0f} fwhm {fw[j]:.2f}')
    if thr == 0.98:
        per = {}
        for j in np.where(oL | oB)[0]:
            k = names[int(fi[j])]
            per.setdefault(k, [0, 0])[0 if oL[j] else 1] += 1
        print('  per file (LFP-only, B-only):')
        for k, v in sorted(per.items(), key=lambda x: -sum(x[1])):
            print(f'    {k}: {v}')
