"""Slide variant of report_fig_saturated.py: the three DEIMOS 600ZD plateau
features that legacy curve_fit keeps and option B rejects, with 20 pt fonts,
drawn at the size placed on the slide.  Re-uses the data functions of
report_fig_saturated.py (same inputs, same examples).

Run (pypeit14b env): python slide_fig_saturated.py
Output: $PYPEIT_DEV/pypeitdev/speed_up/Slides/figures/slide_sat_examples.png
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import report_fig_saturated as R
from slide_fig_style import setup, save, INK, MUTED, C_LEG, C_NEW, WIN

KIND = {('MSC04', 4): 'edge', ('MSC03', 11): 'edge', ('MSC03', 14): 'flat top'}
setup()
want = {(m, i) for m, i, _ in R.EXAMPLES}
found = {(m, i): (s, fw) for m, i, s, fw in R.slits() if (m, i) in want}
fig, axes = plt.subplots(1, 3, figsize=(12.3, 4.3), constrained_layout=True)
for ax, (msc, i, p0) in zip(axes, R.EXAMPLES):
    s, fw = found[(msc, i)]
    x, y, pixt, fitp, res, good = R.run_slit(s, fw)
    j = int(np.where(pixt == p0)[0][0])
    fi = (fitp + fitp % 2)//2
    lo, hi = p0 - 22, p0 + 12
    xx = np.arange(lo, hi + 1)
    ax.plot(xx, y[lo:hi+1]/1e3, color=MUTED, lw=1.5, drawstyle='steps-mid', zorder=1, label='arc')
    w = np.arange(max(p0-fi, 0), p0+fi+1)
    ax.axvspan(w[0]-0.5, w[-1]+0.5, color=WIN, zorder=0)
    ax.plot(w, y[w]/1e3, 'o', ms=7, color=INK, zorder=4, label='fit window')
    fine = np.linspace(lo, hi, 800)
    L, B = res['legacy'], res['B']
    ax.plot(fine, R.gauss(fine, L[0][j], L[1][j], L[2][j])/1e3, color=C_LEG, lw=3, zorder=3,
            label='legacy (kept)')
    ax.plot(fine, R.gauss(fine, B[0][j], B[1][j], B[2][j])/1e3, color=C_NEW, lw=3, ls='--', zorder=3,
            label='option B (rejected)')
    top = max(y[lo:hi+1].max(), L[0][j], 1)/1e3
    ax.set_ylim(-0.03*top, 1.1*top)
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(MaxNLocator(3, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(3, integer=True))
    ax.set_title(f'{msc} slit {i}: {KIND[(msc, i)]}\n'
                 f'$\\sigma$: legacy {L[2][j]:.1f}, B {B[2][j]:.0f} px', loc='left', color=INK)
    ax.set_xlabel('spectral pixel')
axes[0].set_ylabel('counts ($10^3$)')
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc='outside lower center', ncol=4, frameon=False, handlelength=1.6)
save(fig, 'slide_sat_examples.png')
