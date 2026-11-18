"""
PR C, "PR" task 1: compare the wavelength calibrations of two full dev-suite
runs (default: REDUX_OUT_PR/develop vs REDUX_OUT_PR/arcfit).

For every WaveCalib file present in both runs (``ncpu/`` and ``step_by_step/``
duplicates are skipped) and every slit/order:
  - the pypeit_chk_wavecalib diagnostics (WaveCalib.wave_diagnostics):
    Nlin, IDs_Wave_cov(%), measured_fwhm, RMS, and the wavelength range;
  - the change of the wavelength solution, in pixels, over the pixel range
    covered by the lines of both fits (max and median |dlambda| / dispersion);
  - the arc-line measurements stored in WaveFit.tcent (every detected line):
    lines matched within 0.5 px, gained, lost, and the centroid shift of the
    matched lines;
  - the lines of the final fits: common identifications (same wavelength
    within 0.01 A), gained, lost.
Writes one row per slit to <out>/pr_metrics.npz and a per-file CSV.

Usage: python arcfit_pr_metrics.py [ROOT_A ROOT_B OUTDIR]
"""
import sys, warnings, csv
from pathlib import Path
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.wavecalib import WaveCalib

DEV = Path('/home/xavier/Projects/PypeIt/PypeIt-development-suite')
A = Path(sys.argv[1]) if len(sys.argv) > 1 else DEV / 'REDUX_OUT_PR' / 'develop'
B = Path(sys.argv[2]) if len(sys.argv) > 2 else DEV / 'REDUX_OUT_PR' / 'arcfit'
OUT = Path(sys.argv[3]) if len(sys.argv) > 3 else DEV / 'REDUX_OUT_PR'


def diag(wc):
    d = wc.wave_diagnostics(print_diag=False)
    return {k: np.asarray(d[k]) for k in d.colnames}


def final(wv):
    if wv is None or wv.pypeitfit is None:
        return None, None
    g = wv.pypeitfit.gpm.astype(bool)
    return np.asarray(wv.pixel_fit)[g], np.asarray(wv.wave_fit)[g]


def match(a, b, tol):
    """Indices of a and b matched within tol (greedy, nearest)."""
    ia, ib = [], []
    if a.size == 0 or b.size == 0:
        return np.array(ia, int), np.array(ib, int)
    used = np.zeros(b.size, bool)
    for i in np.argsort(a):
        d = np.absolute(b - a[i]); d[used] = np.inf
        j = int(np.argmin(d))
        if d[j] < tol:
            ia.append(i); ib.append(j); used[j] = True
    return np.array(ia, int), np.array(ib, int)


SLIT_TOL = 3   # pixels
rows, slitlog = [], []
# Skip duplicated products: --ncpu and step-by-step re-runs, and quick-look
# (QL) re-reductions, which reuse the setup's calibrations
files = sorted(f for f in A.rglob('Calibrations/WaveCalib_*.fits')
               if '/ncpu/' not in str(f) and 'step_by_step' not in str(f)
               and '/QL' not in str(f.relative_to(A)))
missing = []
for f in files:
    rel = f.relative_to(A)
    fb = B / rel
    if not fb.is_file():
        missing.append(str(rel)); continue
    try:
        wa = WaveCalib.from_file(str(f), chk_version=False)
        wb = WaveCalib.from_file(str(fb), chk_version=False)
        da, db = diag(wa), diag(wb)
    except Exception as e:
        print(f'skip {rel}: {e}'); continue
    inst = rel.parts[0]
    setup = str(rel).split('/Calibrations/')[0]
    key = f.name.replace('WaveCalib_', '').replace('.fits', '')
    # Pair slits by spatial position (spat_id, within SLIT_TOL pixels) or by
    # echelle order: the two runs can trace different slits
    ida = np.asarray(da['SpatOrderID'], dtype=int); idb = np.asarray(db['SpatOrderID'], dtype=int)
    tol = 0 if wa.ech_orders is not None else SLIT_TOL
    pairs, usedb = [], set()
    for i in range(ida.size):
        dd = np.absolute(idb - ida[i]).astype(float)
        dd[list(usedb)] = np.inf
        j = int(np.argmin(dd)) if dd.size else -1
        if j >= 0 and dd[j] <= tol:
            pairs.append((i, j)); usedb.add(j)
    slitlog.append(dict(setup=setup, key=key, n_a=ida.size, n_b=idb.size, paired=len(pairs),
                        only_a=[int(ida[i]) for i in range(ida.size) if i not in [p[0] for p in pairs]],
                        only_b=[int(idb[j]) for j in range(idb.size) if j not in usedb]))
    n = len(pairs)
    for i, j in pairs:
        va, vb = wa.wv_fits[i], wb.wv_fits[j]
        pa, la = final(va); pb, lb = final(vb)
        r = dict(inst=inst, setup=setup, key=key, slit=i, sid=int(ida[i]), sid_b=int(idb[j]),
                 echelle=wa.ech_orders is not None)
        for k, src, ii in [('a', da, i), ('b', db, j)]:
            r[f'nlin_{k}'] = int(src['Nlin'][ii]); r[f'cov_{k}'] = float(src['IDs_Wave_cov(%)'][ii])
            r[f'fwhm_{k}'] = float(src['measured_fwhm'][ii]); r[f'rms_{k}'] = float(src['RMS'][ii])
        r['fit_a'], r['fit_b'] = pa is not None, pb is not None
        r['dpix_max'] = r['dpix_med'] = np.nan
        r['id_common'] = r['id_lost'] = r['id_gained'] = 0
        r['tc_a'] = r['tc_b'] = r['tc_match'] = 0
        r['tc_dmed'] = r['tc_dmax'] = r['tc_dp99'] = np.nan
        if pa is not None and pb is not None and pa.size and pb.size:
            disp = np.median(np.absolute(np.diff(va.wave_soln)))
            lo = int(np.floor(max(pa.min(), pb.min()))); hi = int(np.ceil(min(pa.max(), pb.max()))) + 1
            if hi > lo and disp > 0:
                d = np.absolute(va.wave_soln[lo:hi] - vb.wave_soln[lo:hi]) / disp
                r['dpix_max'], r['dpix_med'] = float(d.max()), float(np.median(d))
            ia, ib = match(la, lb, 0.01)
            r['id_common'], r['id_lost'], r['id_gained'] = ia.size, la.size - ia.size, lb.size - ib.size
        if va is not None and vb is not None and va.tcent is not None and vb.tcent is not None:
            ta, tb = np.asarray(va.tcent, float), np.asarray(vb.tcent, float)
            ta, tb = ta[ta > 0], tb[tb > 0]
            ia, ib = match(ta, tb, 0.5)
            r['tc_a'], r['tc_b'], r['tc_match'] = ta.size, tb.size, ia.size
            if ia.size:
                dd = np.absolute(ta[ia] - tb[ib])
                r['tc_dmed'], r['tc_dmax'], r['tc_dp99'] = float(np.median(dd)), float(dd.max()), float(np.percentile(dd, 99))
        rows.append(r)
    print(f'{setup} {key}: {n} slits paired ({ida.size} / {idb.size})', flush=True)

OUT.mkdir(parents=True, exist_ok=True)
keys = list(rows[0].keys())
with open(OUT / 'pr_metrics.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
np.savez(OUT / 'pr_metrics.npz', **{k: np.array([r[k] for r in rows]) for k in keys})
import json
(OUT / 'pr_slitlists.json').write_text(json.dumps(slitlog, indent=1))
nd = [l for l in slitlog if l['only_a'] or l['only_b']]
print(f'files whose slit lists differ: {len(nd)}')
for l in nd:
    print(f"  {l['setup']} {l['key']}: {l['n_a']} vs {l['n_b']} slits; only develop {l['only_a']}; only PR C {l['only_b']}")
print(f'{len(rows)} slits from {len(files) - len(missing)} files; missing in B: {len(missing)}')
for m in missing:
    print('  missing in B:', m)
