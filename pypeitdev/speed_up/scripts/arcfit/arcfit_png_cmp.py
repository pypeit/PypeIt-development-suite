"""
Task 6 item 2: compare QA PNGs between two dev-suite variants (default legacy vs A):
byte comparison, and for byte-different files the max/mean absolute pixel difference
and the number of differing pixels.
Usage: python arcfit_png_cmp.py REDUX_OUT_C <setup> <png glob> [<setup> <png glob> ...]
"""
import sys, filecmp
from pathlib import Path
import numpy as np
from PIL import Image

root = Path(sys.argv[1])
args = sys.argv[2:]
for setup, pat in zip(args[::2], args[1::2]):
    fa = sorted((root / 'legacy' / setup / 'QA' / 'PNGs').glob(pat))
    nsame = 0
    for f in fa:
        g = root / 'A' / setup / 'QA' / 'PNGs' / f.name
        if not g.exists():
            print(f'  {setup} {f.name}: missing in A'); continue
        if filecmp.cmp(f, g, shallow=False):
            nsame += 1; continue
        ia = np.asarray(Image.open(f).convert('RGB'), dtype=int)
        ib = np.asarray(Image.open(g).convert('RGB'), dtype=int)
        if ia.shape != ib.shape:
            print(f'  {setup} {f.name}: shapes differ {ia.shape} {ib.shape}'); continue
        d = np.abs(ia - ib)
        print(f'  {setup} {f.name}: bytes differ; max|dpix| {d.max()}, '
              f'differing pixels {np.sum(d.max(axis=2) > 0)} / {d.shape[0]*d.shape[1]}')
    print(f'{setup} {pat}: {len(fa)} files, {nsame} byte-identical')
