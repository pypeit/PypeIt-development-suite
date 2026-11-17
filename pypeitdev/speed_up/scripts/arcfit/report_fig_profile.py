"""
Figure for Speed_Up_Arcfit_Report.md: serial (--ncpu 1, cProfile) DEIMOS
600ZD_M_6500 profiles of the five PR-C variants.  Numbers are copied from
Reports/speed_up_results.md, "After workstream C" (from the
keck_deimos_600zd_m_6500_C_{legacy,A,B,A_cache,B_cache}.pstats.txt tables).

Run (any env with matplotlib):  python report_fig_profile.py
Output: $PYPEIT_DEV/pypeitdev/speed_up/Reports/figures/arcfit/profile_deimos.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEV = os.environ.get('PYPEIT_DEV', '/mnt/tank/Astronomy/PypeIt/PypeIt-development-suite')
OUT = f'{DEV}/pypeitdev/speed_up/Reports/figures/arcfit/profile_deimos.png'
variants = ['legacy', 'A', 'B', 'A + cache\n(shipped)', 'B + cache']
wall = [10420, 10094, 9554, 9462, 9495]
fit_arcspec = [851, 526, 44, 57, 44]
measure_fwhm = [801, 571, 126, 3, 1]

INK, MUTED, C1, C2 = '#0b0b0b', '#8a8984', '#2a78d6', '#eb6834'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False,
                     'axes.spines.right': False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
y = np.arange(len(variants))[::-1]
a1.barh(y, wall, color=C1, height=0.6)
for yi, w in zip(y, wall):
    a1.text(w + 80, yi, f'{w:,} s', va='center', fontsize=9, color=INK)
a1.set_yticks(y, variants)
a1.set_xlim(0, 12000)
a1.set_xlabel('wall-clock time, serial cProfile run (s)')
a1.set_title('Whole reduction', loc='left', fontsize=10, color=INK)
h = 0.36
a2.barh(y + h/2, fit_arcspec, height=h, color=C1, label='fit_arcspec (all callers)')
a2.barh(y - h/2, measure_fwhm, height=h, color=C2, label='flexure measure_fwhm (archive sky)')
for yi, v1, v2 in zip(y, fit_arcspec, measure_fwhm):
    a2.text(v1 + 10, yi + h/2, f'{v1} s', va='center', fontsize=8.5, color=INK)
    a2.text(v2 + 10, yi - h/2, f'{v2} s', va='center', fontsize=8.5, color=INK)
a2.set_yticks(y, variants)
a2.set_xlim(0, 1000)
a2.set_xlabel('cumulative time (s); the two overlap')
a2.set_title('Arc-line fitting and the archive-sky FWHM', loc='left', fontsize=10, color=INK)
a2.legend(frameon=False, loc='lower right', fontsize=8.5)
fig.savefig(OUT, dpi=170, bbox_inches='tight', facecolor='white')
print('wrote', OUT)
