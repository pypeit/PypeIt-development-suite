"""
"Revert to CurveFit?" item 2: which dev-suite setups have a different arc line
list under R (curve_fit + median-width floor) than develop (curve_fit, no
floor)?  For every WaveCalib file of the full develop dev-suite run
(REDUX_OUT_PR/develop), re-run detect_lines on each stored arc spectrum with
the pipeline fwhm and sigdetect, with and without the floor, and count the
lines the floor removes and how many of them are in the stored final fit.
Usage: python arcfit_revert_floor_suite.py [root] [out.txt]
"""
import sys, re, warnings
from pathlib import Path
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.core import arc
from pypeit.wavecalib import WaveCalib

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else
            '/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_PR/develop')
OUT = sys.argv[2] if len(sys.argv) > 2 else \
    '/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_PR/revert_floor_suite.txt'


def wpar(f):
    pars = sorted(Path(f).parent.parent.glob('*.par'))
    if not pars:
        return 4.0, True, 5.0
    txt = pars[-1].read_text()
    m = re.search(r'\[\[wavelengths\]\](.*?)(?=\n\s*\[\[)', txt, re.S)
    block = m.group(1) if m else txt
    get = lambda k: re.search(rf'^\s*{k}\s*=\s*(.+)$', block, re.M)
    fwhm = float(get('fwhm').group(1)) if get('fwhm') else 4.0
    froml = get('fwhm_fromlines').group(1).strip() == 'True' if get('fwhm_fromlines') else True
    sd = get('sigdetect')
    sig = float(re.split(r'[,\s]+', sd.group(1).strip())[0]) if sd else 5.0
    return fwhm, froml, sig


files = sorted(f for f in ROOT.rglob('Calibrations/WaveCalib_*.fits')
               if 'step_by_step' not in str(f) and '/ncpu' not in str(f))
rows = []
for k, f in enumerate(files):
    setup = '/'.join(f.relative_to(ROOT).parts[:2])
    try:
        wc = WaveCalib.from_file(str(f), chk_version=False)
    except Exception as e:
        rows.append((setup, f.name, -1, -1, -1, -1)); continue
    fw_par, froml, sig = wpar(f)
    nall = nrem = nfit = nslit = 0
    for wv in wc.wv_fits:
        if wv is None or wv.spec is None:
            continue
        fw = float(wv.fwhm) if froml and wv.fwhm is not None else fw_par
        spec = np.asarray(wv.spec)
        old = arc.detect_lines(spec, sigdetect=sig, fwhm=fw, min_frac_med_width=0.)
        new = arc.detect_lines(spec, sigdetect=sig, fwhm=fw)
        rem = np.setdiff1d(old[5], new[5])
        nall += old[5].size; nrem += rem.size; nslit += 1
        if rem.size and wv.pypeitfit is not None and wv.pixel_fit is not None:
            pf = np.asarray(wv.pixel_fit)[wv.pypeitfit.gpm.astype(bool)]
            if pf.size:
                nfit += int(np.sum([np.min(np.abs(pf - old[2][j])) < 0.5 for j in rem]))
    rows.append((setup, f.name, nslit, nall, nrem, nfit))
    print(f'{k+1}/{len(files)} {setup} {f.name}: removed {nrem}/{nall}, in final fit {nfit}', flush=True)
with open(OUT, 'w') as fh:
    fh.write('setup file nslit nlines nremoved nremoved_in_final_fit\n')
    for r in rows:
        fh.write(' '.join(map(str, r)) + '\n')
print('saved', OUT)
