"""Slide variant of report_fig_profile.py (20 pt fonts, drawn at placed size).
Numbers copied from Reports/speed_up_results.md, "After workstream C".
Run (any env with matplotlib): python slide_fig_profile.py
Output: $PYPEIT_DEV/pypeitdev/speed_up/Slides/figures/slide_profile_deimos.png
"""
import numpy as np
import matplotlib.pyplot as plt
from slide_fig_style import setup, save, INK, C_LEG, C_NEW

variants = ['legacy', 'A', 'B', 'A + cache', 'B + cache']
wall = [10420, 10094, 9554, 9462, 9495]
fit_arcspec = [851, 526, 44, 57, 44]
measure_fwhm = [801, 571, 126, 3, 1]

setup()
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.3, 4.3), constrained_layout=True)
y = np.arange(len(variants))[::-1]
cols = [C_LEG if v != 'A + cache' else '#1b4f91' for v in variants]
a1.barh(y, wall, color=cols, height=0.62)
for yi, w in zip(y, wall):
    a1.text(w - 150, yi, f'{w:,} s', va='center', ha='right', color='white', fontweight='bold')
a1.set_yticks(y, variants)
a1.get_yticklabels()[3].set_fontweight('bold')
a1.set_xlim(0, 11000)
a1.set_xticks([0, 5000, 10000], ['0', '5,000', '10,000'])
a1.set_xlabel('serial wall-clock (s)')
a1.set_title('Whole reduction', loc='left', color=INK)
h = 0.38
a2.barh(y + h/2, fit_arcspec, height=h, color=C_LEG, label='fit_arcspec')
a2.barh(y - h/2, measure_fwhm, height=h, color=C_NEW, label='archive-sky measure_fwhm')
for yi, v1, v2 in zip(y, fit_arcspec, measure_fwhm):
    a2.text(v1 + 12, yi + h/2, f'{v1}', va='center', color=INK)
    a2.text(v2 + 12, yi - h/2, f'{v2}', va='center', color=INK)
a2.set_yticks(y, [''] * len(y))
a2.set_xlim(0, 1000)
a2.set_xticks([0, 250, 500, 750])
a2.set_xlabel('cumulative time (s); overlapping')
a2.set_title('Arc-line fitting vs archive-sky FWHM', loc='left', color=INK)
a2.legend(frameon=False, loc='lower right', handlelength=1.0, bbox_to_anchor=(1.03, -0.02))
save(fig, 'slide_profile_deimos.png')
