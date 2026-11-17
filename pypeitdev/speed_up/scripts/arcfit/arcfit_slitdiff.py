"""
Per-slit comparison of one WaveCalib file between two dev-suite variants:
RMS, number of lines, identified wavelengths, and the change of the solution.
Usage: python arcfit_slitdiff.py REDUX_OUT_C <relative WaveCalib path> [legacy gn10]
"""
import sys, warnings
from pathlib import Path
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.wavecalib import WaveCalib
root, rel = Path(sys.argv[1]), sys.argv[2]
a, b = (sys.argv[3:5] if len(sys.argv) > 4 else ['legacy', 'gn10'])
A = WaveCalib.from_file(str(root / a / rel), chk_version=False)
B = WaveCalib.from_file(str(root / b / rel), chk_version=False)
print(f'{rel}: {a} vs {b}')
print(' slit  spat_id | rms_a  rms_b | n_a n_b | common IDs  only_a only_b | max|dlam| pix (line range) | method')
for i, (wa, wb) in enumerate(zip(A.wv_fits, B.wv_fits)):
    sid = A.spat_ids[i] if A.spat_ids is not None else i
    if wa is None or wb is None or wa.pypeitfit is None or wb.pypeitfit is None:
        print(f' {i:4d} {sid:8d} | fit a: {wa is not None and wa.pypeitfit is not None}  fit b: {wb is not None and wb.pypeitfit is not None}')
        continue
    ia = set(np.round(np.asarray(wa.wave_fit)[wa.pypeitfit.gpm.astype(bool)], 2))
    ib = set(np.round(np.asarray(wb.wave_fit)[wb.pypeitfit.gpm.astype(bool)], 2))
    disp = np.median(np.absolute(np.diff(wa.wave_soln)))
    pa = np.asarray(wa.pixel_fit)[wa.pypeitfit.gpm.astype(bool)]
    pb = np.asarray(wb.pixel_fit)[wb.pypeitfit.gpm.astype(bool)]
    lo = int(np.floor(max(pa.min(), pb.min())))
    hi = int(np.ceil(min(pa.max(), pb.max()))) + 1
    d = np.max(np.absolute(wa.wave_soln[lo:hi] - wb.wave_soln[lo:hi])) / disp
    flag = '  <--' if d > 0.05 or wb.rms > 1.1 * wa.rms + 1e-3 else ''
    print(f' {i:4d} {sid:8d} | {wa.rms:.4f} {wb.rms:.4f} | {len(ia):3d} {len(ib):3d} | {len(ia & ib):4d} {len(ia - ib):6d} {len(ib - ia):6d} | {d:10.3e}{flag}')
