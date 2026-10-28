"""
Figure for Speed_Up_Arcfit_Report.md: Paranal sky lines on which the legacy
per-line curve_fit converges to a negative sigma (gauss_3deg uses only sigma**2,
so the sign is arbitrary) and detect_lines then rejects the line through
twid > 0.  Compares legacy, the shipped default (A) and option B.

Run (pypeit14b env):
  python report_fig_negsigma.py --list   # statistics and every negative-sigma line
  python report_fig_negsigma.py          # write the figure
Input : pypeit's paranal_sky.fits (as in test_detect_lines, default detect_lines args)
Output: $PYPEIT_DEV/pypeitdev/speed_up/Reports/figures/arcfit/negsigma_example.png
"""
import os, sys, warnings
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.core import arc, skyspec

DEV = os.environ.get('PYPEIT_DEV', '/mnt/tank/Astronomy/PypeIt/PypeIt-development-suite')
OUT = f'{DEV}/pypeitdev/speed_up/Reports/figures/arcfit/negsigma_example.png'
EXAMPLE_PIX = [6479, 16750]   # peak pixels of the examples (chosen from --list)
VARIANTS = {'legacy': (False, True), 'A': (True, True), 'B': (True, False)}
FWHM = 4.0           # detect_lines default


def run():
    flux = skyspec.load_sky_spectrum('paranal_sky.fits').flux
    calls, orig = [], arc.fit_arcspec
    arc.fit_arcspec = lambda *a: (calls.append(a), orig(*a))[1]
    try:
        arc.detect_lines(flux)
    finally:
        arc.fit_arcspec = orig
    x, y, pixt, fitp = calls[0]
    res = {}
    for v, (use, match) in VARIANTS.items():
        arc.USE_VECTORIZED_ARCFIT, arc.ARCFIT_MATCH_CURVEFIT = use, match
        res[v] = orig(x, y, pixt, fitp)
    arc.USE_VECTORIZED_ARCFIT, arc.ARCFIT_MATCH_CURVEFIT = True, True
    good = {v: np.isfinite(r[2]) & (r[2] > 0) & (r[2] < 3*FWHM/2.35) & (r[1] > 0) & (r[1] < x[-1])
               & (np.abs(r[1] - pixt) < 0.75*FWHM) for v, r in res.items()}
    # legacy would pass every cut if its sigma had the other sign
    L = res['legacy']
    neg = (L[2] < 0) & (L[2] != -999.)
    flip = neg & (-L[2] < 3*FWHM/2.35) & (L[1] > 0) & (np.abs(L[1] - pixt) < 0.75*FWHM)
    return flux, x, y, pixt, fitp, res, good, neg, flip


def report():
    flux, x, y, pixt, fitp, res, good, neg, flip = run()
    print(f'{pixt.size} lines; accepted: ' + ', '.join(f'{v} {g.sum()}' for v, g in good.items()))
    print(f'legacy negative sigma: {neg.sum()}; would pass the other cuts with |sigma|: {flip.sum()}')
    for v in ['A', 'B']:
        r = res[v]
        kept = flip & good[v]
        d = np.abs(r[1][kept] - res['legacy'][1][kept])
        ds = np.abs(r[2][kept]/np.abs(res['legacy'][2][kept]) - 1)
        print(f'  {v}: keeps {kept.sum()} of them; |dc| median {np.median(d) if d.size else np.nan:.1e} '
              f'max {d.max() if d.size else np.nan:.1e} px; | |sigma|/|sigma_legacy|-1 | max {ds.max() if ds.size else np.nan:.1e}')
    for j in np.where(flip)[0]:
        print(f'pix {pixt[j]:5d} amp {res["legacy"][0][j]:9.1f} | legacy c {res["legacy"][1][j]:.4f} s {res["legacy"][2][j]:.4f} '
              f'| A c {res["A"][1][j]:.4f} s {res["A"][2][j]:.4f} {good["A"][j]} | B c {res["B"][1][j]:.4f} s {res["B"][2][j]:.4f} {good["B"][j]}')


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    INK, MUTED, C_LEG, C_NEW = '#0b0b0b', '#8a8984', '#2a78d6', '#eb6834'
    plt.rcParams.update({'font.size': 10, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                         'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False,
                         'axes.spines.right': False})
    flux, x, y, pixt, fitp, res, good, neg, flip = run()
    fig, axes = plt.subplots(1, len(EXAMPLE_PIX), figsize=(11, 4.2), constrained_layout=True)
    fi = (fitp + fitp % 2)//2
    for ax, pex in zip(axes, EXAMPLE_PIX):
        j = int(np.where(pixt == pex)[0][0])
        p0 = pixt[j]
        lo, hi = p0 - 12, p0 + 12
        xx = np.arange(lo, hi + 1)
        ax.plot(xx, y[lo:hi+1], color=MUTED, lw=1, drawstyle='steps-mid', label='continuum-subtracted sky')
        w = np.arange(p0 - fi, p0 + fi + 1)
        ax.axvspan(w[0]-0.5, w[-1]+0.5, color='#e9e8e4', zorder=0)
        ax.plot(w, y[w], 'o', ms=6, color=INK, zorder=4, label=f'fit window ({w.size} px)')
        fine = np.linspace(lo, hi, 1200)
        L, B = res['legacy'], res['B']
        g = lambda a, c, s: a*np.exp(-0.5*(fine - c)**2/s**2)
        ax.plot(fine, g(L[0][j], L[1][j], L[2][j]), color=C_LEG, lw=2.5,
                label='legacy curve_fit ($\\sigma$ < 0: rejected)')
        ax.plot(fine, g(B[0][j], B[1][j], B[2][j]), color=C_NEW, lw=2, ls='--',
                label='vectorized, option B (kept)')
        ax.text(0.02, 0.97, f'legacy: c = {L[1][j]:.3f}, $\\sigma$ = {L[2][j]:+.3f} px\n'
                            f'B:      c = {B[1][j]:.3f}, $\\sigma$ = {B[2][j]:+.3f} px\n'
                            f'shipped A: legacy result (rejected)',
                transform=ax.transAxes, va='top', fontsize=8.5, family='monospace', color=INK)
        ymin, ymax = y[lo:hi+1].min(), y[lo:hi+1].max()
        ax.set_ylim(ymin - 0.05*(ymax-ymin), ymax + 0.45*(ymax-ymin))
        ax.set_xlim(lo, hi)
        ax.set_xlabel('spectral pixel')
        ax.set_title(f'Paranal sky, peak pixel {p0}', loc='left', fontsize=10, color=INK)
        print(pex, {v: (round(r[0][j], 3), round(r[1][j], 4), round(r[2][j], 4), bool(good[v][j]))
                    for v, r in res.items()})
    axes[0].set_ylabel('flux (continuum-subtracted)')
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.09))
    fig.savefig(OUT, dpi=170, bbox_inches='tight', facecolor='white')
    print('wrote', OUT)


if __name__ == '__main__':
    report() if '--list' in sys.argv else plot()
