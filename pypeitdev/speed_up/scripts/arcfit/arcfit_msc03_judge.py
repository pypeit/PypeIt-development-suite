"""
Task 6 item 1c: which solution of keck_deimos/600ZD_M_6500 MSC03 slit 11 (spat_id
1063) is correct, legacy or A?  Two independent checks:
 1. dispersion dlambda/dpix vs wavelength of slit 11 compared with the neighbouring
    slits (whose solutions are identical in both runs) and the arxiv template;
 2. number of detected lines (tcent) in the blue that fall within 0.5 / 1 px of a
    line-list line when mapped through each solution.
Also makes a figure.  Usage: python arcfit_msc03_judge.py [REDUX_OUT_C] [outdir]
"""
import sys, warnings
from pathlib import Path
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pypeit.wavecalib import WaveCalib
from pypeit.core.wavecal import waveio

root = Path(sys.argv[1] if len(sys.argv) > 1 else
            '/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C')
outdir = Path(sys.argv[2] if len(sys.argv) > 2 else '.')
rel = 'keck_deimos/600ZD_M_6500/Calibrations/WaveCalib_A_0_MSC03.fits'
W = {v: WaveCalib.from_file(str(root / v / rel), chk_version=False) for v in ['legacy', 'A']}
ll, _, _ = waveio.load_line_lists(list(W['A'].par['lamps']), include_unknown=False)
llw = np.asarray(ll['wave'])
tw, ts, tb, *_ = waveio.load_template(W['A'].par['reid_arxiv'], 1)

def disp_vs_wave(ws):
    return 0.5 * (ws[1:] + ws[:-1]), np.diff(ws)

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
wgrid = np.arange(4300., 5600., 100.)
print('dlambda/dpix at wavelength (A/pix):')
print('            ' + ' '.join(f'{w:6.0f}' for w in wgrid))
rows = {}
for lab, ws in [('slit10', W['legacy'].wv_fits[10].wave_soln), ('slit12', W['legacy'].wv_fits[12].wave_soln),
                ('slit9', W['legacy'].wv_fits[9].wave_soln), ('slit13', W['legacy'].wv_fits[13].wave_soln),
                ('template', tw), ('11 legacy', W['legacy'].wv_fits[11].wave_soln),
                ('11 A', W['A'].wv_fits[11].wave_soln)]:
    wm, d = disp_vs_wave(np.asarray(ws))
    ok = (wm > 3000) & np.isfinite(d) & (d > 0)
    rows[lab] = np.interp(wgrid, wm[ok], d[ok], left=np.nan)
    print(f'{lab:11s} ' + ' '.join(f'{x:6.4f}' for x in rows[lab]))
    style = dict(lw=2.5, color='C3' if 'legacy' in lab else 'C0') if lab.startswith('11') else dict(lw=1, color='0.6' if lab != 'template' else 'k', ls='-' if lab != 'template' else ':')
    ax[0].plot(wm[ok], d[ok], label=lab, **style)
ax[0].set_xlim(4100, 6500); ax[0].set_ylim(0.6, 0.68)
ax[0].set_xlabel('wavelength (A)'); ax[0].set_ylabel('dlambda/dpix (A/pix)'); ax[0].legend(fontsize=8)
ax[0].set_title('600ZD MSC03: dispersion, slit 11 vs neighbours')

print('\nDetected lines (tcent, both runs identical) mapped through each solution; blue part (pix < 2000):')
for v in ['legacy', 'A']:
    f = W[v].wv_fits[11]
    tc = np.asarray(f.tcent); tc = tc[(tc > 0) & (tc < 2000)]
    wl = np.interp(tc, np.arange(f.wave_soln.size), f.wave_soln)
    disp = np.interp(tc, np.arange(f.wave_soln.size - 1) + 0.5, np.diff(f.wave_soln))
    dpx = np.array([np.min(np.abs(llw - w)) for w in wl]) / disp
    print(f'  {v:6s}: {tc.size} lines; within 0.5 px of a line-list line: {np.sum(dpx < 0.5)}, within 1 px: {np.sum(dpx < 1)}; '
          f'median offset {np.median(dpx):.2f} px')
    ax[1].plot(tc, dpx, 'o', color='C3' if v == 'legacy' else 'C0', label=v, ms=6 if v == 'legacy' else 4)
ax[1].axhline(0.5, color='k', lw=0.5)
ax[1].set_yscale('log'); ax[1].set_xlabel('detected line centroid (pix)')
ax[1].set_ylabel('|nearest line-list line| (pix)'); ax[1].legend()
ax[1].set_title('slit 11 (spat_id 1063): blue detected lines')
fig.tight_layout()
fig.savefig(outdir / 'msc03_slit11_judge.png', dpi=110)
print('wrote', outdir / 'msc03_slit11_judge.png')
