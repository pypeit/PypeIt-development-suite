"""
Summaries for Task 10 from arcfit_sigma_floor.npz (see arcfit_sigma_floor.py).
Usage: python arcfit_sigma_floor_summary.py [npz]
"""
import sys
import numpy as np
npz = sys.argv[1] if len(sys.argv) > 1 else '/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C/arcfit_sigma_floor.npz'
d = np.load(npz, allow_pickle=True)
R, names = d['rows'], d['names']
fi, slit, pixt, sA, sB, gA, gB, fw, half, raw, infit = (R[:, k] for k in range(11))
gA, gB, infit = gA > 0, gB > 0, infit > 0
exp = fw / 2.355
sig = np.where(gB, sB, sA)          # the sigma B would report (A's where only A accepts)
floors = {'abs085 (sigma > 0.85 px)': 0.85 + 0*fw,
          'rel050 (sigma > 0.5 fwhm/2.355)': 0.5*exp,
          'rule3 (sigma > 0.35 half-window)': 0.35*half}
both, onlyA, onlyB = gA & gB, gA & ~gB, gB & ~gA
print(f'{len(names)} files; lines accepted by A (= legacy): {gA.sum()}, by B: {gB.sum()}, both {both.sum()}, only A {onlyA.sum()}, only B {onlyB.sum()}')
print(f'A-accepted lines used in the stored final fits (matched within 0.5 px): {np.sum(gA & infit)}')
inst = np.array([n.split('/')[0] for n in names])
for name, thr in floors.items():
    cutB = gB & (sB <= thr)           # B-accepted lines the floor removes
    cutA = gA & (sA <= thr)           # legacy/A-accepted lines the floor would remove if applied to A's sigma
    print(f'\n== {name}')
    print(f'  removes {cutB.sum()} of {gB.sum()} B-accepted lines: {np.sum(cutB & onlyB)} of the {onlyB.sum()} B-only, '
          f'{np.sum(cutB & both)} accepted by both (legacy too); {np.sum(cutB & both & infit)} of those are in final fits')
    print(f'  applied to A/legacy sigma: removes {cutA.sum()} legacy-accepted lines, {np.sum(cutA & infit)} in final fits')
    # per setup
    print('  per setup (B-accepted removed / of which legacy-accepted / in final fit | legacy lines total):')
    for k, n in enumerate(names):
        m = fi == k
        if np.sum(cutB & m) or np.sum(cutA & m):
            print(f'    {n:55s} {np.sum(cutB & m):4d} / {np.sum(cutB & both & m):4d} / {np.sum(cutB & both & infit & m):3d} | {np.sum(gA & m):5d}   (A-sigma: {np.sum(cutA & m)} removed, {np.sum(cutA & infit & m)} in fit)')
# distribution of sigma/expected for legacy-accepted lines, per instrument
print('\n== sigma_A / expected for legacy-accepted lines, per instrument: p1 p5 p50; fraction < 0.5; fraction sigma_A < 0.85 px')
for ins in np.unique(inst):
    m = gA & np.isin(fi, np.where(inst == ins)[0])
    r = sA[m]/exp[m]
    print(f'  {ins:18s} n {m.sum():5d}  ' + ' '.join(f'{np.percentile(r, q):.2f}' for q in [1, 5, 50]) + f'   {np.mean(r < 0.5):.4f}   {np.mean(sA[m] < 0.85):.4f}   fwhm range {fw[m].min():.2f}-{fw[m].max():.2f}')
# the B-only (spike) population
if onlyB.any():
    print(f'\n== B-only lines: sigma_B/expected p5/p50/p95: ' + ' '.join(f'{np.percentile(sB[onlyB]/exp[onlyB], q):.2f}' for q in [5, 50, 95])
          + f'; max sigma_B {sB[onlyB].max():.2f} px; by instrument: ' + ', '.join(f'{ins} {np.sum(onlyB & np.isin(fi, np.where(inst == ins)[0]))}' for ins in np.unique(inst)))
# narrow lines curve_fit accepts: both-accepted with sigma_B below rule3 threshold
nb = both & (sB < 0.35*half)
print(f'\n== lines accepted by both with sigma_B below the rule-3 threshold (A handed them off; curve_fit accepted): {nb.sum()}; '
      f'in final fits: {np.sum(nb & infit)}; sigma_A/expected p50 {np.percentile(sA[nb]/exp[nb], 50) if nb.any() else np.nan:.2f}; '
      f'|sigma_A - sigma_B| p50 {np.percentile(np.abs(sA[nb]-sB[nb]), 50) if nb.any() else np.nan:.3f}')
