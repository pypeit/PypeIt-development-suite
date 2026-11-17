"""
Compare the calibration products of the PR C dev-suite variants (Task 5).

For every WaveCalib and Tilts file found under the reference variant
(default ``legacy``), compare the same file in each other variant:

- wavelength: per-slit RMS (the quantity vet_tests/test_wavelengths.py
  asserts on), number of lines in the final fit, and the change of the
  wavelength solution, expressed in pixels (|dlambda| / dispersion);
- tilts: number of traced lines and of 2D-fit points, and the change of the
  tilt image, in spectral pixels.

Usage: python arcfit_metrics.py REDUX_OUT_C [legacy gn10 gn0] > table.md
"""
import sys, glob, warnings
from pathlib import Path
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.wavecalib import WaveCalib
from pypeit.wavetilts import WaveTilts
from pypeit.slittrace import SlitTraceSet

root = Path(sys.argv[1])
variants = sys.argv[2:] if len(sys.argv) > 2 else ['legacy', 'gn10', 'gn0']
ref = variants[0]

def load(cls, f):
    try:
        return cls.from_file(str(f), chk_version=False)
    except Exception as e:
        return None

def wave_stats(wc):
    out = []
    for w in wc.wv_fits:
        if w is None or w.rms is None or w.pypeitfit is None:
            out.append(None); continue
        n = int(np.sum(w.pypeitfit.gpm)) if w.pypeitfit.gpm is not None else len(w.pixel_fit)
        # Pixel range covered by the fitted lines; the solutions are only
        # compared there (outside it, they are polynomial extrapolations)
        pix = np.asarray(w.pixel_fit)[w.pypeitfit.gpm.astype(bool)] if w.pypeitfit.gpm is not None \
                else np.asarray(w.pixel_fit)
        rng = (int(np.floor(pix.min())), int(np.ceil(pix.max())) + 1) if pix.size else (0, 0)
        out.append((float(w.rms), n, None if w.wave_soln is None else np.asarray(w.wave_soln), rng))
    return out

def short(f):
    p = f.relative_to(root / ref).parts
    sub = '/'.join(p[:-2])
    return f'{sub} {f.name.split("_",1)[1].replace(".fits","")}'

print('## Wavelength solutions\n')
print('| setup / file | slits | ' + ' | '.join(f'{v} RMS med / max' for v in variants)
      + ' | ' + ' | '.join(f'{v} lines' for v in variants) + ' | '
      + ' | '.join(f'{v} max \\|Δλ\\| in common line range (pix)' for v in variants[1:]) + ' | RMS worse by >10% |')
print('|' + '---|' * (2 + 2*len(variants) + len(variants) - 1 + 1))
tot = {v: 0 for v in variants}; worse_tot = 0
for f in sorted((root / ref).rglob('Calibrations/WaveCalib_*.fits')):
    rel = f.relative_to(root / ref)
    st = {v: (lambda wc: wave_stats(wc) if wc is not None else None)(load(WaveCalib, root / v / rel)) for v in variants}
    if st[ref] is None:
        continue
    good = [i for i, s in enumerate(st[ref]) if s is not None]
    cells_rms, cells_n, cells_d = [], [], []
    worse = 0
    for v in variants:
        if st[v] is None:
            cells_rms.append('missing'); cells_n.append('—'); continue
        r = [st[v][i][0] for i in good if st[v][i] is not None]
        n = sum(st[v][i][1] for i in good if st[v][i] is not None)
        tot[v] += n
        cells_rms.append(f'{np.median(r):.4f} / {np.max(r):.4f}' if r else '—')
        cells_n.append(str(n))
    for v in variants[1:]:
        if st[v] is None:
            cells_d.append('—'); continue
        d = []
        for i in good:
            a, b = st[ref][i], st[v][i]
            if b is None:
                worse += 1; continue
            if b[0] > 1.1 * a[0] and b[0] - a[0] > 1e-3:
                worse += 1
            if a[2] is not None and b[2] is not None:
                disp = np.median(np.absolute(np.diff(a[2])))
                # Pixel range covered by the lines of both variants
                lo, hi = max(a[3][0], b[3][0]), min(a[3][1], b[3][1])
                if hi > lo:
                    d.append(np.max(np.absolute(a[2][lo:hi] - b[2][lo:hi])) / disp)
        cells_d.append(f'{max(d):.2e}' if d else '—')
    worse_tot += worse
    print(f'| {short(f)} | {len(good)} | ' + ' | '.join(cells_rms) + ' | ' + ' | '.join(cells_n)
          + ' | ' + ' | '.join(cells_d) + f' | {worse} |')
print('\nTotal lines in the final fits: ' + ', '.join(f'{v} {tot[v]}' for v in variants)
      + f'; slit-variant pairs with RMS worse by >10%: {worse_tot}\n')

print('## Tilts\n')
print('| setup / file | ' + ' | '.join(f'{v} lines / fit pts' for v in variants) + ' | '
      + ' | '.join(f'{v} med / max \\|Δtilt\\| (pix)' for v in variants[1:]) + ' |')
print('|' + '---|' * (1 + len(variants) + len(variants) - 1))
for f in sorted((root / ref).rglob('Calibrations/Tilts_*.fits')):
    rel = f.relative_to(root / ref)
    wts = {v: load(WaveTilts, root / v / rel) for v in variants}
    if wts[ref] is None:
        continue
    cells = []
    for v in variants:
        t = wts[v]
        if t is None or t.tilt_traces is None:
            cells.append('missing'); continue
        nl = sum(len(np.unique(np.asarray(r)[np.isfinite(np.asarray(r, dtype=float))]))
                 for r in t.tilt_traces['good2dfit_lid'])
        npts = int(sum(np.sum(np.isfinite(np.asarray(r, dtype=float))) for r in t.tilt_traces['good2dfit_tilt']))
        cells.append(f'{nl} / {npts}')
    # Tilt images, on the slits of the reference variant
    sfile = sorted(f.parent.glob('Slits_*' + f.name.split('_', 1)[1].replace('.fits', '') + '*'))
    slits = load(SlitTraceSet, sfile[0]) if sfile else None
    dcells = []
    if slits is not None:
        slitimg = slits.slit_img()
        img_ref = wts[ref].fit2tiltimg(slitimg)
        nspec = slitimg.shape[0]
        on = slitimg > 0
        for v in variants[1:]:
            if wts[v] is None:
                dcells.append('—'); continue
            d = np.absolute(wts[v].fit2tiltimg(slitimg) - img_ref)[on] * (nspec - 1)
            dcells.append(f'{np.median(d):.2e} / {np.max(d):.3f}')
    else:
        dcells = ['no slits'] * (len(variants) - 1)
    print(f'| {short(f)} | ' + ' | '.join(cells) + ' | ' + ' | '.join(dcells) + ' |')
