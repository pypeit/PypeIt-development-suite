"""
Arc_1dfit-style comparison of one slit's wavelength solution between two
dev-suite runs (rows 1-2: PypeIt's Arc_1dfit QA panels for each run, with
identifications unique to one run in orange; row 3: the change of the
solution in pixels).

Usage: python arcfit_pr_slitfig.py ROOT_A ROOT_B REL_WAVECALIB SLIT OUT.png [LABEL_A LABEL_B]
"""
import sys, warnings
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.wavecalib import WaveCalib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from report_fig_arc1dfit import arc1dfit_row, C_LEG, C_NEW, INK

ra, rb, rel, slit, out = sys.argv[1:6]
la, lb = (sys.argv[6:8] if len(sys.argv) > 7 else ['develop', 'develop + PR C'])
slit = int(slit)
A = WaveCalib.from_file(f'{ra}/{rel}', chk_version=False)
B = WaveCalib.from_file(f'{rb}/{rel}', chk_version=False)
# Pair by spatial position (or echelle order), not by index
ida = A.ech_orders if A.ech_orders is not None else A.spat_ids
idb = B.ech_orders if B.ech_orders is not None else B.spat_ids
j = int(np.argmin(np.absolute(np.asarray(idb) - ida[slit])))
if abs(int(idb[j]) - int(ida[slit])) > 3:
    raise SystemExit(f'slit {slit} (id {ida[slit]}) has no counterpart in {rb}')
wa, wb = A.wv_fits[slit], B.wv_fits[j]
plt.rcdefaults()
plt.rcParams.update({'font.family': 'serif', 'font.size': 9})
fig = plt.figure(figsize=(12, 11.5))
outer = gridspec.GridSpec(3, 1, figure=fig, hspace=0.42, height_ratios=[1, 1, 0.6], top=0.93, bottom=0.05,
                          left=0.07, right=0.98)
for k, (w, o, name, col) in enumerate([(wa, wb, la, C_LEG), (wb, wa, lb, C_NEW)]):
    g = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[k], width_ratios=[1.35, 1], wspace=0.22, hspace=0.08)
    arc1dfit_row(fig, g, w, name, col, o)
ax = fig.add_subplot(outer[2])
d = (wb.wave_soln - wa.wave_soln) / wa.cen_disp
ax.plot(np.arange(d.size), d, color=INK, lw=1.2)
ga, gb = wa.pypeitfit.gpm.astype(bool), wb.pypeitfit.gpm.astype(bool)
ax.plot(wa.pixel_fit[ga], np.interp(wa.pixel_fit[ga], np.arange(d.size), d), '|', color=C_LEG, ms=10, label=f'{la} lines')
ax.plot(wb.pixel_fit[gb], np.interp(wb.pixel_fit[gb], np.arange(d.size), d), '|', color=C_NEW, ms=10, label=f'{lb} lines')
ax.axhline(0, color='k', ls='--', lw=0.7); ax.set_xlim(0, d.size)
ax.set_xlabel('Pixel'); ax.set_ylabel(f'{lb} − {la} (pixels)'); ax.legend(frameon=False, fontsize=8)
ax.set_title('change of the wavelength solution', loc='left', fontsize=9)
sid = A.spat_ids[slit] if A.ech_orders is None else A.ech_orders[slit]
fig.suptitle(f"{rel.split('/Calibrations/')[0]} {rel.split('WaveCalib_')[1].replace('.fits','')}, "
             f"slit {slit} ({'order' if A.ech_orders is not None else 'spat_id'} {sid})", fontsize=12, y=0.995)
fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
print('wrote', out, f'| RMS {wa.rms:.3f} -> {wb.rms:.3f}, lines {ga.sum()} -> {gb.sum()}, max|d| {np.abs(d).max():.2f} px')
