"""Slide variant of report_fig_negsigma.py: two Paranal sky lines on which
legacy curve_fit converges to sigma < 0 (rejected by twid > 0), with 20 pt
fonts, drawn at the size placed on the slide.  Re-uses report_fig_negsigma.run().

Run (pypeit14b env): python slide_fig_negsigma.py
Output: $PYPEIT_DEV/pypeitdev/speed_up/Slides/figures/slide_negsigma_example.png
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import report_fig_negsigma as R
from slide_fig_style import setup, save, INK, MUTED, C_LEG, C_NEW, WIN

setup()
flux, x, y, pixt, fitp, res, good, neg, flip = R.run()
fig, axes = plt.subplots(1, 2, figsize=(12.3, 4.3), constrained_layout=True)
fi = (fitp + fitp % 2)//2
for ax, pex in zip(axes, R.EXAMPLE_PIX):
    j = int(np.where(pixt == pex)[0][0])
    p0 = pixt[j]
    lo, hi = p0 - 10, p0 + 10
    xx = np.arange(lo, hi + 1)
    ax.plot(xx, y[lo:hi+1], color=MUTED, lw=1.5, drawstyle='steps-mid', label='sky')
    w = np.arange(p0 - fi, p0 + fi + 1)
    ax.axvspan(w[0]-0.5, w[-1]+0.5, color=WIN, zorder=0)
    ax.plot(w, y[w], 'o', ms=7, color=INK, zorder=4, label='fit window')
    fine = np.linspace(lo, hi, 1200)
    L, B = res['legacy'], res['B']
    g = lambda a, c, s: a*np.exp(-0.5*(fine - c)**2/s**2)
    ax.plot(fine, g(L[0][j], L[1][j], L[2][j]), color=C_LEG, lw=3.5, label='legacy ($\\sigma$<0, rejected)')
    ax.plot(fine, g(B[0][j], B[1][j], B[2][j]), color=C_NEW, lw=3, ls='--', label='option B (kept)')
    ymin, ymax = y[lo:hi+1].min(), y[lo:hi+1].max()
    ax.set_ylim(ymin - 0.05*(ymax-ymin), max(ymax, L[0][j], B[0][j]) + 0.08*(ymax-ymin))
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(MaxNLocator(3, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.set_xlabel('spectral pixel')
    ax.set_title(f'Paranal sky, pixel {p0}\n$\\sigma$: legacy {L[2][j]:+.2f}, B {B[2][j]:+.2f} px', loc='left', color=INK)
axes[0].set_ylabel('flux')
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc='outside lower center', ncol=4, frameon=False, handlelength=1.6)
save(fig, 'slide_negsigma_example.png')
