"""
Task 6 item 2: compare per-object WAVE_RMS in the spec1d files of two dev-suite
variants (default: REDUX_OUT_C/legacy vs REDUX_OUT_C/A) for the given setups.
Usage: python arcfit_wave_rms.py REDUX_OUT_C setup_glob [setup_glob ...]
 e.g.  python arcfit_wave_rms.py REDUX_OUT_C keck_deimos/600ZD_M_6500 'shane_kast_*/*'
"""
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits

root = Path(sys.argv[1])
a, b = 'legacy', 'A'
for pat in sys.argv[2:]:
    for sdir in sorted((root / a).glob(pat)):
        setup = sdir.relative_to(root / a)
        nobj = ndiff = 0; maxd = 0.; missing = []; diffs = []
        for f in sorted(sdir.glob('Science/spec1d_*.fits')):
            g = root / b / setup / 'Science' / f.name
            if not g.exists():
                missing.append(f.name); continue
            with fits.open(f) as ha, fits.open(g) as hb:
                ea = {h.name: h.header.get('WAVE_RMS') for h in ha[1:] if 'WAVE_RMS' in h.header}
                eb = {h.name: h.header.get('WAVE_RMS') for h in hb[1:] if 'WAVE_RMS' in h.header}
            if set(ea) != set(eb):
                missing.append(f'{f.name}: objects differ ({len(ea)} vs {len(eb)}; '
                               f'only {a}: {sorted(set(ea)-set(eb))[:4]} only {b}: {sorted(set(eb)-set(ea))[:4]})')
            for k in sorted(set(ea) & set(eb)):
                nobj += 1
                d = abs(ea[k] - eb[k])
                maxd = max(maxd, d)
                if d > 1e-4:
                    ndiff += 1; diffs.append((f.name, k, ea[k], eb[k]))
        print(f'{str(setup):35s} objects {nobj:4d}  |dWAVE_RMS|>1e-4: {ndiff:3d}  max|d| {maxd:.2e}')
        for x in diffs:
            print(f'    {x[0]} {x[1]}: {a} {x[2]:.4f} {b} {x[3]:.4f}')
        for m in missing:
            print('    !', m)
