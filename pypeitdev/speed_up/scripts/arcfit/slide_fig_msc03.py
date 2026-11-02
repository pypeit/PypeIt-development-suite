"""Slide variant of arcfit_msc03_judge.py (DEIMOS 600ZD MSC03 slit 11, the
known limitation), with 20 pt fonts, drawn at the size placed on the slide.
Same inputs and the same two checks: dispersion vs neighbours/template, and
the offset of the blue detected lines from the line list through each solution.

Run (pypeit14b env): python slide_fig_msc03.py
Input : $PYPEIT_DEV/REDUX_OUT_C/{legacy,A}/keck_deimos/600ZD_M_6500/Calibrations/WaveCalib_A_0_MSC03.fits
Output: $PYPEIT_DEV/pypeitdev/speed_up/Slides/figures/slide_msc03_slit11.png
"""
import warnings
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pypeit.wavecalib import WaveCalib
from pypeit.core.wavecal import waveio
from slide_fig_style import setup, save, DEV, INK, C_LEG, C_NEW

C_LEGACY, C_A = '#c0392b', C_LEG
rel = 'keck_deimos/600ZD_M_6500/Calibrations/WaveCalib_A_0_MSC03.fits'
W = {v: WaveCalib.from_file(f'{DEV}/REDUX_OUT_C/{v}/{rel}', chk_version=False) for v in ['legacy', 'A']}
ll, _, _ = waveio.load_line_lists(list(W['A'].par['lamps']), include_unknown=False)
llw = np.asarray(ll['wave'])
tw, *_ = waveio.load_template(W['A'].par['reid_arxiv'], 1)

setup()
fig, ax = plt.subplots(1, 2, figsize=(12.3, 4.3), constrained_layout=True)
first = True
for lab, ws in [('slit10', W['legacy'].wv_fits[10].wave_soln), ('slit12', W['legacy'].wv_fits[12].wave_soln),
                ('slit9', W['legacy'].wv_fits[9].wave_soln), ('slit13', W['legacy'].wv_fits[13].wave_soln),
                ('template', tw), ('legacy', W['legacy'].wv_fits[11].wave_soln),
                ('A', W['A'].wv_fits[11].wave_soln)]:
    ws = np.asarray(ws)
    wm, d = 0.5*(ws[1:] + ws[:-1]), np.diff(ws)
    ok = (wm > 3000) & np.isfinite(d) & (d > 0)
    if lab in ('legacy', 'A'):
        ax[0].plot(wm[ok], d[ok], lw=3.5, color=C_LEGACY if lab == 'legacy' else C_A, label=f'slit 11, {lab}')
    elif lab == 'template':
        ax[0].plot(wm[ok], d[ok], lw=1.5, color='k', ls=':', label='template')
    else:
        ax[0].plot(wm[ok], d[ok], lw=1.5, color='0.65', label='neighbours' if first else None)
        first = False
ax[0].set_xlim(4100, 6500); ax[0].set_ylim(0.6, 0.672)
ax[0].xaxis.set_major_locator(MaxNLocator(3)); ax[0].yaxis.set_major_locator(MaxNLocator(3))
ax[0].set_xlabel('wavelength (Å)'); ax[0].set_ylabel('Å / pixel')
ax[0].set_title('Dispersion vs neighbours', loc='left', color=INK)
ax[0].legend(frameon=False, loc='upper left', ncol=2, handlelength=1.0, columnspacing=0.8,
             borderaxespad=0.1)
for v in ['legacy', 'A']:
    f = W[v].wv_fits[11]
    tc = np.asarray(f.tcent); tc = tc[(tc > 0) & (tc < 2000)]
    wl = np.interp(tc, np.arange(f.wave_soln.size), f.wave_soln)
    disp = np.interp(tc, np.arange(f.wave_soln.size - 1) + 0.5, np.diff(f.wave_soln))
    dpx = np.array([np.min(np.abs(llw - w)) for w in wl]) / disp
    print(v, tc.size, 'within 0.5 px:', np.sum(dpx < 0.5), 'median %.2f' % np.median(dpx))
    ax[1].plot(tc, dpx, 'o', color=C_LEGACY if v == 'legacy' else C_A, label=v, ms=10 if v == 'legacy' else 8)
ax[1].axhline(0.5, color='k', lw=1, ls='--')
ax[1].set_yscale('log'); ax[1].set_ylim(0.02, 60)
ax[1].set_xlabel('blue line centroid (pixel)')
ax[1].set_ylabel('offset to list (px)')
ax[1].xaxis.set_major_locator(MaxNLocator(3))
ax[1].set_title('Blue lines vs line list', loc='left', color=INK)
ax[1].legend(frameon=False, loc='upper right', ncol=2, handlelength=0.8, borderaxespad=0.1)
save(fig, 'slide_msc03_slit11.png')
