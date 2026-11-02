"""
Figure for Speed_Up_Arcfit_Report.md: fitted Gaussian sigma of every arc line
that detect_lines accepts, against wavelength, for several instruments and
configurations, under the shipped default (A) and option B
(ARCFIT_MATCH_CURVEFIT = False).  Lines that B keeps but A rejects are
highlighted; so are the lines A keeps and B rejects.

As in report_fig_saturated.py, detect_lines is re-run on the arc spectra stored
in the legacy dev-suite WaveCalib files (per-slit spectrum and measured fwhm),
the fit_arcspec inputs are captured, and each line is fit under A and B.  The
wavelength of a line is its B (or A) centroid mapped through the slit's stored
wavelength solution.

Run (pypeit14b env, PYPEIT_DEV set):
  python report_fig_sigma.py          # write the figure
  python report_fig_sigma.py --list   # also print every A/B-only line
Input : $PYPEIT_DEV/REDUX_OUT_C/legacy/<setups below>/Calibrations/WaveCalib_*.fits
Output: $PYPEIT_DEV/pypeitdev/speed_up/Reports/figures/arcfit/sigma_vs_wave.png
"""
import os, sys, warnings
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.core import arc
from pypeit.wavecalib import WaveCalib

DEV = os.environ.get('PYPEIT_DEV', '/mnt/tank/Astronomy/PypeIt/PypeIt-development-suite')
ROOT = f'{DEV}/REDUX_OUT_C/legacy/'
OUT = f'{DEV}/pypeitdev/speed_up/Reports/figures/arcfit/sigma_vs_wave.png'
SETUPS = [
    ('Kast blue 600/4310 d55', 'shane_kast_blue/600_4310_d55/shane_kast_blue_A/Calibrations/WaveCalib_A_0_DET01.fits'),
    ('Kast red 600/7500 d57', 'shane_kast_red/600_7500_d57/Calibrations/WaveCalib_A_0_DET01.fits'),
    ('LRIS blue multi 300/5000 d680, DET01', 'keck_lris_blue/multi_300_5000_d680/Calibrations/WaveCalib_A_0_DET01.fits'),
    ('DEIMOS 600ZD_M_6500 MSC01 (16 slits)', 'keck_deimos/600ZD_M_6500/Calibrations/WaveCalib_A_0_MSC01.fits'),
    ('DEIMOS 830G_L_8100', 'keck_deimos/830G_L_8100/Calibrations/WaveCalib_A_0_MSC03.fits'),
    ('DEIMOS 1200G_M_7750 MSC03', 'keck_deimos/1200G_M_7750/Calibrations/WaveCalib_A_all_MSC03.fits'),
]
VARIANTS = {'A': True, 'B': False}   # ARCFIT_MATCH_CURVEFIT


def fit_slit(spec, fwhm):
    """Capture detect_lines' fit_arcspec call and fit the lines under A and B."""
    calls, orig = [], arc.fit_arcspec
    arc.fit_arcspec = lambda *a: (calls.append(a), orig(*a))[1]
    try:
        arc.detect_lines(spec, sigdetect=5., fwhm=fwhm)
    finally:
        arc.fit_arcspec = orig
    x, y, pixt, fitp = calls[0]
    res, good = {}, {}
    for v, match in VARIANTS.items():
        arc.ARCFIT_MATCH_CURVEFIT = match
        res[v] = orig(x, y, pixt, fitp)
        r = res[v]
        good[v] = np.isfinite(r[2]) & (r[2] > 0) & (r[2] < 3*fwhm/2.35) & (r[1] > 0) & (r[1] < x[-1]) \
                  & (np.abs(r[1] - pixt) < 0.75*fwhm)
    arc.ARCFIT_MATCH_CURVEFIT = True
    fit_interval = (fitp + fitp % 2)//2
    return pixt, res, good, fit_interval


def collect(rel):
    wc = WaveCalib.from_file(ROOT + rel, chk_version=False)
    rows = []   # (wave, sigA, sigB, goodA, goodB, fwhm, fit_interval, slit, pixt)
    for i, w in enumerate(wc.wv_fits):
        if w is None or w.spec is None or w.fwhm is None or w.wave_soln is None:
            continue
        spec, fw = np.asarray(w.spec), float(w.fwhm)
        pixt, res, good, fi = fit_slit(spec, fw)
        pix = np.arange(spec.size)
        for j in range(pixt.size):
            if not (good['A'][j] or good['B'][j]):
                continue
            c = res['B'][1][j] if good['B'][j] else res['A'][1][j]
            rows.append((np.interp(c, pix, w.wave_soln), res['A'][2][j], res['B'][2][j],
                         good['A'][j], good['B'][j], fw, fi, i, pixt[j], spec[pixt[j]]))
    return np.array(rows, dtype=float)


def plot(list_lines=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    INK, MUTED, C_LEG, C_NEW = '#0b0b0b', '#8a8984', '#2a78d6', '#eb6834'
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                         'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False,
                         'axes.spines.right': False})
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.4), constrained_layout=True)
    summary = []
    for ax, (name, rel) in zip(axes.ravel(), SETUPS):
        R = collect(rel)
        wave, sA, sB, gA, gB, fw, fi = R[:, 0], R[:, 1], R[:, 2], R[:, 3] > 0, R[:, 4] > 0, R[:, 5], R[:, 6]
        both, onlyB, onlyA = gA & gB, gB & ~gA, gA & ~gB
        ax.scatter(wave[both], sB[both], s=7, color=MUTED, alpha=0.6, lw=0, label='accepted by A and B')
        ax.scatter(wave[onlyA], sA[onlyA], s=34, facecolor='none', edgecolor=C_LEG, lw=1.2,
                   label='accepted only by A (legacy result; B: width cut)')
        ax.scatter(wave[onlyB], sB[onlyB], s=34, color=C_NEW, zorder=4,
                   label='accepted only by B (A hands off to curve_fit, which fails)')
        fwm, fim = np.median(fw), np.median(fi)
        ax.axhline(fwm/2.355, color=INK, lw=0.8, ls=':')
        ax.axhline(3*fwm/2.35, color=INK, lw=0.8, ls='--')
        ax.axhline(arc.ARCFIT_MIN_WIDTH*fim, color=C_NEW, lw=0.8, ls='-.')
        x1 = wave.max()
        ax.text(x1, fwm/2.355, ' expected $\\sigma$', va='center', ha='left', fontsize=7, color=INK)
        ax.text(x1, 3*fwm/2.35, ' width cut', va='center', ha='left', fontsize=7, color=INK)
        ax.text(x1, arc.ARCFIT_MIN_WIDTH*fim, ' A hand-off', va='center', ha='left', fontsize=7, color=C_NEW)
        ax.set_yscale('log')
        from matplotlib.ticker import ScalarFormatter, NullFormatter
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_yticks([t for t in [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10] if ax.get_ylim()[0] <= t <= ax.get_ylim()[1]])
        ax.set_title(f'{name}\nfwhm {fw.min():.1f}–{fw.max():.1f} px; both {both.sum()}, only B {onlyB.sum()}, only A {onlyA.sum()}',
                     loc='left', fontsize=8.5, color=INK)
        ax.set_xlabel('wavelength (Å)')
        ax.set_ylabel('fitted $\\sigma$ (px)')
        ax.set_xlim(wave.min() - 0.02*(x1 - wave.min()), x1 + 0.18*(x1 - wave.min()))
        summary.append((name, both.sum(), onlyB.sum(), onlyA.sum(), fwm,
                        np.median(sB[onlyB]/(fw[onlyB]/2.355)) if onlyB.any() else np.nan,
                        np.median(sA[onlyA]/(fw[onlyA]/2.355)) if onlyA.any() else np.nan,
                        np.median(sB[onlyA]) if onlyA.any() else np.nan))
        if list_lines:
            for j in np.where(onlyB | onlyA)[0]:
                print(f'{name}: slit {int(R[j, 7])} pix {int(R[j, 8])} wave {wave[j]:.1f} raw {R[j, 9]:.0f} '
                      f'sigA {sA[j]:.3f} sigB {sB[j]:.3f} ' + ('only B' if onlyB[j] else 'only A'))
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.05))
    fig.savefig(OUT, dpi=170, bbox_inches='tight', facecolor='white')
    print('wrote', OUT)
    print('setup | both | only B | only A | median fwhm | median sigB/expected (only B) | median sigA/expected (only A) | median sigB (only A)')
    for s in summary:
        print(' | '.join(str(v) if isinstance(v, (int, np.integer)) else (f'{v:.2f}' if isinstance(v, float) else v) for v in s))


if __name__ == '__main__':
    plot('--list' in sys.argv)
