"""
Figure for Speed_Up_Arcfit_Report.md: saturated, flat-topped DEIMOS 600ZD arc
features that legacy curve_fit keeps and the pure-vectorized fit (option B,
ARCFIT_MATCH_CURVEFIT=False) rejects.

Re-runs detect_lines on the arc spectra stored in the legacy dev-suite
WaveCalib files (wv_fits[i].spec, with the slit's measured fwhm, i.e. what the
pipeline saw), captures the fit_arcspec inputs, and fits each line with
  legacy : USE_VECTORIZED_ARCFIT=False (per-line curve_fit)
  A      : shipped default (ARCFIT_MATCH_CURVEFIT=True)
  B      : ARCFIT_MATCH_CURVEFIT=False
then applies detect_lines' good-line cuts.

Run (pypeit14b env; no git, read-only on the data):
  python report_fig_saturated.py --list     # print every legacy-kept, B-rejected line
  python report_fig_saturated.py            # write the figure
Input : $PYPEIT_DEV/REDUX_OUT_C/legacy/keck_deimos/600ZD_M_6500/Calibrations/WaveCalib_A_0_MSC0*.fits
Output: $PYPEIT_DEV/pypeitdev/speed_up/Reports/figures/arcfit/sat_examples.png
"""
import os, sys, glob, warnings
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.core import arc, fitting
from pypeit.wavecalib import WaveCalib

DEV = os.environ.get('PYPEIT_DEV', '/mnt/tank/Astronomy/PypeIt/PypeIt-development-suite')
ROOT = f'{DEV}/REDUX_OUT_C/legacy/keck_deimos/600ZD_M_6500/Calibrations/'
OUT = f'{DEV}/pypeitdev/speed_up/Reports/figures/arcfit/sat_examples.png'
# (mosaic, slit index, pixt) of the examples shown; chosen from --list
EXAMPLES = [('MSC04', 4, 6573), ('MSC03', 11, 5407), ('MSC03', 14, 5895)]

VARIANTS = {'legacy': dict(USE_VECTORIZED_ARCFIT=False, ARCFIT_MATCH_CURVEFIT=True),
            'A': dict(USE_VECTORIZED_ARCFIT=True, ARCFIT_MATCH_CURVEFIT=True),
            'B': dict(USE_VECTORIZED_ARCFIT=True, ARCFIT_MATCH_CURVEFIT=False)}


def run_slit(spec, fwhm):
    """detect_lines inputs for one spectrum, and fit_arcspec results per variant."""
    calls, orig = [], arc.fit_arcspec
    arc.fit_arcspec = lambda *a: (calls.append(a), orig(*a))[1]
    try:
        arc.detect_lines(spec, sigdetect=5., fwhm=fwhm)
    finally:
        arc.fit_arcspec = orig
    x, y, pixt, fitp = calls[0]
    res = {}
    for v, opts in VARIANTS.items():
        for k, val in opts.items():
            setattr(arc, k, val)
        res[v] = orig(x, y, pixt, fitp)
    arc.USE_VECTORIZED_ARCFIT, arc.ARCFIT_MATCH_CURVEFIT = True, True
    # detect_lines' good-line cuts (max_frac_fwhm=3, 0.75*fwhm centre cut)
    good = {v: np.isfinite(r[2]) & (r[2] > 0) & (r[2] < 3*fwhm/2.35) & (r[1] > 0)
               & (r[1] < x[-1]) & (np.abs(r[1] - pixt) < 0.75*fwhm) for v, r in res.items()}
    return x, y, pixt, fitp, res, good


def slits():
    for f in sorted(glob.glob(ROOT + 'WaveCalib_A_0_MSC0*.fits')):
        msc = f.split('_')[-1][:5]
        wc = WaveCalib.from_file(f, chk_version=False)
        for i, w in enumerate(wc.wv_fits):
            if w is None or w.spec is None or w.fwhm is None:
                continue
            yield msc, i, np.asarray(w.spec), float(w.fwhm)


def list_lines():
    tot = {v: 0 for v in VARIANTS}
    n = 0
    for msc, i, s, fw in slits():
        x, y, pixt, fitp, res, good = run_slit(s, fw)
        n += 1
        for v in tot:
            tot[v] += good[v].sum()
        for j in np.where(good['legacy'] & ~good['B'])[0]:
            p = pixt[j]
            seg = s[max(p-4, 0):p+5]
            print(f'{msc} slit {i:2d} fwhm {fw:.2f} pix {p:4d} raw peak {s[p]:7.0f} | legacy c {res["legacy"][1][j]:.3f} '
                  f's {res["legacy"][2][j]:.3f} | B c {res["B"][1][j]:.3f} s {res["B"][2][j]:.3f} '
                  f'| max step/peak {np.max(np.abs(np.diff(seg)))/s[p]:.2f}')
    print(f'{n} spectra; accepted lines: ' + ', '.join(f'{k} {v}' for k, v in tot.items()))


def gauss(x, a, c, s):
    return a*np.exp(-0.5*(x - c)**2/s**2)


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    INK, MUTED, C_LEG, C_NEW = '#0b0b0b', '#8a8984', '#2a78d6', '#eb6834'
    plt.rcParams.update({'font.size': 10, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                         'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False,
                         'axes.spines.right': False})
    want = {(m, i): p for m, i, p in EXAMPLES}
    found = {}
    for msc, i, s, fw in slits():
        if (msc, i) in want:
            found[(msc, i)] = (s, fw)
    fig, axes = plt.subplots(1, len(EXAMPLES), figsize=(13, 4.4), constrained_layout=True)
    for ax, (msc, i, p0) in zip(axes, EXAMPLES):
        s, fw = found[(msc, i)]
        x, y, pixt, fitp, res, good = run_slit(s, fw)
        j = int(np.where(pixt == p0)[0][0])
        fi = (fitp + fitp % 2)//2
        lo, hi = p0 - 30, p0 + 15
        xx = np.arange(lo, hi + 1)
        ax.plot(xx, y[lo:hi+1]/1e3, color=MUTED, lw=1, drawstyle='steps-mid', zorder=1,
                label='continuum-subtracted arc')
        w = np.arange(max(p0-fi, 0), p0+fi+1)
        ax.axvspan(w[0]-0.5, w[-1]+0.5, color='#e9e8e4', zorder=0)
        ax.plot(w, y[w]/1e3, 'o', ms=6, color=INK, zorder=4, label=f'fit window ({w.size} px)')
        fine = np.linspace(lo, hi, 800)
        L, B = res['legacy'], res['B']
        ax.plot(fine, gauss(fine, L[0][j], L[1][j], L[2][j])/1e3, color=C_LEG, lw=2, zorder=3,
                label='legacy curve_fit (kept)')
        ax.plot(fine, gauss(fine, B[0][j], B[1][j], B[2][j])/1e3, color=C_NEW, lw=2, ls='--', zorder=3,
                label='vectorized, option B (rejected)')
        ax.axvline(L[1][j], color=C_LEG, lw=1, ls=':')
        top = max(y[lo:hi+1].max(), 1)/1e3
        ax.set_ylim(-0.03*top, 1.32*top)
        ax.set_xlim(lo, hi)
        txt = (f'raw peak {s[p0]:,.0f} counts\n'
               f'legacy: c = {L[1][j]:.2f}, $\\sigma$ = {L[2][j]:.2f} px\n'
               f'B:      c = {B[1][j]:.2f}, $\\sigma$ = {B[2][j]:.1f} px\n'
               f'cut: $\\sigma$ < {3*fw/2.35:.2f} px')
        ax.text(0.02, 0.98, txt, transform=ax.transAxes, va='top', ha='left', fontsize=8.5,
                color=INK, family='monospace')
        ax.set_title(f'600ZD {msc} slit {i}, peak pixel {p0}', fontsize=10, color=INK, loc='left')
        ax.set_xlabel('spectral pixel')
    axes[0].set_ylabel('counts (continuum-subtracted, $10^3$)')
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.08))
    fig.savefig(OUT, dpi=170, bbox_inches='tight', facecolor='white')
    print('wrote', OUT)
    for msc, i, p0 in EXAMPLES:
        s, fw = found[(msc, i)]
        x, y, pixt, fitp, res, good = run_slit(s, fw)
        j = int(np.where(pixt == p0)[0][0])
        print(msc, i, p0, 'fwhm %.2f fitp %d raw %.0f' % (fw, fitp, s[p0]),
              {v: (round(r[1][j], 3), round(r[2][j], 3), bool(good[v][j])) for v, r in res.items()})


if __name__ == '__main__':
    list_lines() if '--list' in sys.argv else plot()
