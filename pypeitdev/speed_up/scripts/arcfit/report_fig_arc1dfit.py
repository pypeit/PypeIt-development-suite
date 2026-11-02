"""
Figures for Speed_Up_Arcfit_Report.md: Arc_1dfit-style comparison of one slit's
wavelength solution between the legacy dev-suite reduction and option B
(REDUX_OUT_C/gn10), for the two cases discussed in report section 4:

  830g   : keck_deimos 830G_L_8100, slit 0 (spat_id 762): RMS 0.051 -> 0.27 px
  msc01s2: keck_deimos 600ZD_M_6500 MSC01, slit 2 (spat_id 442): blue-end
           misidentification

Rows 1-2 reproduce PypeIt's Arc_1dfit QA panels (spectrum with the identified
lines, wavelength fit, residuals) for legacy and B, from the stored WaveCalib
products.  Row 3 shows what differs: the lines detect_lines accepts under A and
B (re-run on the stored arc spectrum with the slit's measured fwhm, as in
report_fig_saturated.py), and the change of the solution in pixels.

Run (pypeit14b env, PYPEIT_DEV set):
  python report_fig_arc1dfit.py 830g
  python report_fig_arc1dfit.py msc01s2
Input : $PYPEIT_DEV/REDUX_OUT_C/{legacy,gn10}/keck_deimos/<setup>/Calibrations/WaveCalib_*.fits
Output: $PYPEIT_DEV/pypeitdev/speed_up/Reports/figures/arcfit/arc1dfit_<case>.png
"""
import os, sys, warnings
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.wavecalib import WaveCalib
from report_fig_saturated import run_slit

DEV = os.environ.get('PYPEIT_DEV', '/mnt/tank/Astronomy/PypeIt/PypeIt-development-suite')
ROOT = f'{DEV}/REDUX_OUT_C/'
OUTDIR = f'{DEV}/pypeitdev/speed_up/Reports/figures/arcfit/'
CASES = {
    '830g': dict(rel='keck_deimos/830G_L_8100/Calibrations/WaveCalib_A_0_MSC03.fits', slit=0,
                 label='DEIMOS 830G_L_8100, slit 0 (spat_id 762)', zoom=(1900, 5400)),
    'msc01s2': dict(rel='keck_deimos/600ZD_M_6500/Calibrations/WaveCalib_A_0_MSC01.fits', slit=2,
                    label='DEIMOS 600ZD_M_6500 MSC01, slit 2 (spat_id 442)', zoom=(400, 2500), line=7088),
}
INK, MUTED, C_LEG, C_NEW = '#0b0b0b', '#8a8984', '#2a78d6', '#eb6834'


def load(case):
    c = CASES[case]
    return {v: WaveCalib.from_file(ROOT + v + '/' + c['rel'], chk_version=False).wv_fits[c['slit']]
            for v in ['legacy', 'gn10']}


def in_other(x, lam, other):
    g = other.pypeitfit.gpm.astype(bool)
    return bool(np.any((np.abs(other.pixel_fit[g] - x) < 0.5) & (np.abs(other.wave_fit[g] - lam) < 0.05)))


def arc1dfit_row(fig, gs_row, w, name, color, other, idfont=6):
    """One Arc_1dfit-style row: spectrum + IDs (left), fit (upper right), residuals (lower right)."""
    spec = np.asarray(w.spec)
    gpm = w.pypeitfit.gpm.astype(bool)
    ax = fig.add_subplot(gs_row[:, 0])
    ax.plot(np.arange(spec.size), spec, lw=0.6, color='#1f4e9c')
    ymin, ymax = max(1., spec.min()), spec.max()*4
    yscl = (1.2, 1.5, 1.7)
    for x in w.tcent:
        yl = spec[max(int(x)-2, 0):min(int(x)+2, spec.size-1)].max()
        ax.plot([x, x], [yl*yscl[0], yl*yscl[1]], '-', color='gray', lw=0.6)
    # identified lines absent from the *other* variant (no line within 0.5 px with the same wavelength)
    uniq = gpm & ~np.array([in_other(x, lam, other) for x, lam in zip(w.pixel_fit, w.wave_fit)])
    for x, lam, ion, g, u in zip(w.pixel_fit, w.wave_fit, w.ions, gpm, uniq):
        yl = spec[max(int(x)-2, 0):min(int(x)+2, spec.size-1)].max()
        col = C_NEW if u else 'green'
        ax.plot([x, x], [yl*yscl[0], yl*yscl[1]], '-', color=col, lw=1.2 if u else 0.8)
        ax.text(x, yl*yscl[2], f'{ion} {lam:g}', ha='center', va='bottom', size=idfont, rotation=90.,
                color=col, weight='bold' if u else 'normal')
    ax.set_xlim(0., spec.size)
    ax.set_ylim(ymin, ymax)
    ax.set_yscale('log')
    ax.set_xlabel('Pixel')
    ax.set_ylabel('Flux')
    ax.minorticks_on()
    ax.set_title(f'{name}: {gpm.sum()} lines, RMS = {w.rms:.3f} px'
                 + ('   (orange: identified only in this variant)' if uniq.any() else ''),
                 loc='left', color=color, fontsize=10, weight='bold', y=1.12)
    # fit
    axf = fig.add_subplot(gs_row[0, 1])
    axf.scatter(w.pixel_fit[gpm], w.wave_fit[gpm], marker='x', s=14, color='#1f4e9c')
    if np.any(~gpm):
        axf.scatter(w.pixel_fit[~gpm], w.wave_fit[~gpm], marker='o', edgecolor='gray', facecolor='none', s=18)
    axf.plot(np.arange(spec.size), w.wave_soln, 'r-', lw=0.8)
    axf.set_xlim(0., spec.size)
    axf.set_ylim(w.wave_soln.min()*.95, w.wave_soln.max()*1.05)
    axf.set_ylabel('Wavelength')
    axf.get_xaxis().set_ticks([])
    axf.minorticks_on()
    axf.text(0.05, 0.86, r'$\Delta\lambda$=%.3f$\AA$ (per pix)' % w.cen_disp, size='small', transform=axf.transAxes)
    axf.text(0.05, 0.72, 'RMS=%.3f (pixels)' % w.rms, size='small', transform=axf.transAxes)
    # residuals
    axr = fig.add_subplot(gs_row[1, 1])
    res = (w.wave_fit - w.pypeitfit.eval(w.pixel_fit/w.xnorm)) / w.cen_disp
    axr.scatter(w.pixel_fit[gpm & ~uniq], res[gpm & ~uniq], marker='x', s=14, color='#1f4e9c')
    axr.scatter(w.pixel_fit[uniq], res[uniq], marker='o', s=26, color=C_NEW, zorder=3,
                label='only in this variant')
    axr.plot([0, spec.size], [0, 0], 'k--', lw=0.7)
    axr.set_xlim(0., spec.size)
    axr.set_xlabel('Pixel')
    axr.set_ylabel('Residuals (Pix)')
    axr.minorticks_on()
    return ax, axf, axr


def plot(case):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    c = CASES[case]
    W = load(case)
    L, B = W['legacy'], W['gn10']
    spec, fw = np.asarray(L.spec), float(L.fwhm)
    plt.rcdefaults()
    plt.rcParams.update({'font.family': 'serif', 'font.size': 9})
    fig = plt.figure(figsize=(12, 13.5))
    outer = gridspec.GridSpec(3, 1, figure=fig, hspace=0.38, height_ratios=[1, 1, 0.95], top=0.94, bottom=0.05,
                              left=0.07, right=0.98)
    top = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[0], width_ratios=[1.35, 1], wspace=0.22, hspace=0.08)
    mid = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[1], width_ratios=[1.35, 1], wspace=0.22, hspace=0.08)
    arc1dfit_row(fig, top, L, 'legacy (= shipped A)', C_LEG, B)
    arc1dfit_row(fig, mid, B, 'option B (ARCFIT_MATCH_CURVEFIT = False)', C_NEW, L)

    # Row 3: detect_lines under A and B on the stored spectrum, and the solution change
    bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[2], width_ratios=[1.35, 1], wspace=0.22)
    x, y, pixt, fitp, res, good = run_slit(spec, fw)
    gA, gB = good['A'], good['B']
    cA, cB = res['A'][1], res['B'][1]
    moved = gA & gB & (np.abs(cA - cB) > 0.05)
    ax = fig.add_subplot(bot[0, 0])
    from matplotlib.lines import Line2D
    if 'line' in c:
        # the one line whose centre differs between A and B: its fit window and both Gaussians
        j = int(np.where(pixt == c['line'])[0][0])
        p0 = pixt[j]
        fi = (fitp + fitp % 2)//2
        lo, hi = p0 - 40, p0 + 40
        xx = np.arange(lo, hi + 1)
        ax.plot(xx, y[lo:hi+1]/1e3, color=MUTED, lw=1, drawstyle='steps-mid', label='continuum-subtracted arc')
        wdw = np.arange(p0 - fi, p0 + fi + 1)
        ax.axvspan(wdw[0]-0.5, wdw[-1]+0.5, color='#e9e8e4', zorder=0)
        ax.plot(wdw, y[wdw]/1e3, 'o', ms=5, color=INK, zorder=4, label=f'fit window ({wdw.size} px)')
        fine = np.linspace(lo, hi, 1000)
        g = lambda a, cc, sg: a*np.exp(-0.5*(fine - cc)**2/sg**2)/1e3
        RA, RB = res['A'], res['B']
        ax.plot(fine, g(RA[0][j], RA[1][j], RA[2][j]), color=C_LEG, lw=2, label='A / legacy curve_fit')
        ax.plot(fine, g(RB[0][j], RB[1][j], RB[2][j]), color=C_NEW, lw=2, ls='--', label='option B')
        ax.axvline(RA[1][j], color=C_LEG, lw=1, ls=':')
        ax.axvline(RB[1][j], color=C_NEW, lw=1, ls=':')
        top = y[lo:hi+1].max()/1e3
        ax.set_ylim(-0.05*top, 1.35*top)
        ax.set_xlim(lo, hi)
        ax.text(0.02, 0.97, f'raw peak {spec[p0]:,.0f} counts\n'
                            f'A: c = {RA[1][j]:.2f}, $\\sigma$ = {RA[2][j]:.2f} px\n'
                            f'B: c = {RB[1][j]:.2f}, $\\sigma$ = {RB[2][j]:.2f} px\n'
                            f'$\\Delta$c = {RB[1][j]-RA[1][j]:+.2f} px; both accepted',
                transform=ax.transAxes, va='top', fontsize=8, family='monospace', color=INK)
        ax.set_ylabel('counts (continuum-subtracted, $10^3$)')
        ax.set_xlabel('Pixel')
        ax.legend(loc='upper right', fontsize=7, frameon=False)
        ax.set_title(f'the only detect_lines difference on this slit: the plateau edge at pixel {p0}\n'
                     f'(stored arc, fwhm {fw:.2f} px; A {gA.sum()}, B {gB.sum()} lines; B centre 3.2 px off, dropped by the dev-suite B run)',
                     loc='left', fontsize=8.5)
    else:
        lo, hi = c['zoom']
        ax.plot(np.arange(spec.size), y, lw=0.7, color=MUTED, drawstyle='steps-mid')
        ymax = y[lo:hi].max()
        ax.set_yscale('symlog', linthresh=max(50., 0.002*ymax))
        ax.set_xlim(lo, hi)
        ax.set_ylim(min(-0.02*ymax, y[lo:hi].min()), ymax*6)
        sel = lambda m: (pixt >= lo) & (pixt < hi) & m
        for p in pixt[sel(gA & gB)]:
            ax.plot([p, p], [y[p]*1.3, y[p]*1.8], color='green', lw=0.8)
        for j in np.where(sel(gB & ~gA))[0]:
            p = pixt[j]
            ax.plot([p, p], [y[p]*1.3, y[p]*2.6], color=C_NEW, lw=1.8)
            ax.text(p, y[p]*2.9, f'$\\sigma$={res["B"][2][j]:.2f}', rotation=90, ha='center', va='bottom',
                    fontsize=7, color=C_NEW)
        for j in np.where(sel(gA & ~gB))[0]:
            p = pixt[j]
            ax.plot([p, p], [y[p]*1.3, y[p]*2.6], color=C_LEG, lw=1.8)
        for j in np.where(sel(moved))[0]:
            p = pixt[j]
            ax.plot([p, p], [y[p]*1.3, y[p]*2.6], color='#7a3db8', lw=1.8)
        hs = [Line2D([], [], color='green', lw=0.8, label='accepted by A and B'),
              Line2D([], [], color=C_NEW, lw=1.8, label='accepted only by B (A hands off to curve_fit, which fails)'),
              Line2D([], [], color=C_LEG, lw=1.8, label='accepted only by A'),
              Line2D([], [], color='#7a3db8', lw=1.8, label='centre differs > 0.05 px')]
        ax.legend(handles=hs, loc='upper left', fontsize=7, frameon=False, ncol=2)
        ax.set_xlabel('Pixel')
        ax.set_ylabel('continuum-subtracted flux')
        ax.set_title(f'detect_lines on the stored arc (fwhm {fw:.2f} px, window {fitp + fitp % 2 + 1} px):\n'
                     f'A {gA.sum()}, B {gB.sum()} lines; only A {np.sum(gA & ~gB)}, only B {np.sum(gB & ~gA)}, '
                     f'centre moved {moved.sum()}', loc='left', fontsize=8.5)
    # solution change
    ax2 = fig.add_subplot(bot[0, 1])
    d = (B.wave_soln - L.wave_soln) / L.cen_disp
    ax2.plot(np.arange(spec.size), d, color=INK, lw=1.2)
    gL, gBf = L.pypeitfit.gpm.astype(bool), B.pypeitfit.gpm.astype(bool)
    ax2.plot(L.pixel_fit[gL], np.zeros(gL.sum()) + 0, '|', color=C_LEG, ms=9, label='legacy identified lines')
    ax2.plot(B.pixel_fit[gBf], np.interp(B.pixel_fit[gBf], np.arange(spec.size), d), '|', color=C_NEW, ms=9,
             label='B identified lines')
    ax2.axhline(0, color='k', ls='--', lw=0.7)
    ax2.set_xlim(0, spec.size)
    ax2.set_xlabel('Pixel')
    ax2.set_ylabel('B $-$ legacy solution (pixels)')
    ax2.set_title('change of the wavelength solution\n(B $-$ legacy, in pixels)', loc='left', fontsize=8.5)
    ax2.legend(loc='best', fontsize=7, frameon=False)
    ax2.minorticks_on()
    fig.suptitle(c['label'] + ': legacy vs option B', fontsize=12, y=1.0)
    out = OUTDIR + f'arc1dfit_{case}.png'
    fig.savefig(out, dpi=170, bbox_inches='tight', facecolor='white')
    print('wrote', out)
    # numbers for the report
    print(f'legacy: {gL.sum()} lines rms {L.rms:.4f}; B: {gBf.sum()} lines rms {B.rms:.4f}')
    print('solution change (px) at pixels 0/1000/2000/4000/8000:',
          [round(float(d[p]), 2) for p in [0, 1000, 2000, 4000, 8000]])
    print('max |d| over the range covered by both line sets:',
          round(float(np.max(np.abs(d[int(max(L.pixel_fit[gL].min(), B.pixel_fit[gBf].min())):
                                      int(min(L.pixel_fit[gL].max(), B.pixel_fit[gBf].max()))]))), 3))
    for j in np.where(gB & ~gA)[0]:
        print(f'  only B: pix {pixt[j]} cB {cB[j]:.2f} sigB {res["B"][2][j]:.3f} raw {spec[pixt[j]]:.0f} '
              f'legacy c {res["legacy"][1][j]:.2f} s {res["legacy"][2][j]:.3f} lam {np.interp(cB[j], np.arange(spec.size), B.wave_soln):.2f}')
    for j in np.where(gA & ~gB)[0]:
        print(f'  only A: pix {pixt[j]} cA {cA[j]:.2f} sigA {res["A"][2][j]:.3f} sigB {res["B"][2][j]:.1f} raw {spec[pixt[j]]:.0f}')
    for j in np.where(moved)[0]:
        print(f'  moved: pix {pixt[j]} cA {cA[j]:.3f} cB {cB[j]:.3f} sigA {res["A"][2][j]:.3f} sigB {res["B"][2][j]:.3f} raw {spec[pixt[j]]:.0f}')
    resL = (L.wave_fit - L.pypeitfit.eval(L.pixel_fit/L.xnorm)) / L.cen_disp
    resB = (B.wave_fit - B.pypeitfit.eval(B.pixel_fit/B.xnorm)) / B.cen_disp
    for p, lam, r, g in zip(B.pixel_fit, B.wave_fit, resB, gBf):
        if g and not in_other(p, lam, L):
            print(f'  B-only ID: pix {p:.2f} lam {lam:.2f} resid {r:+.3f} px')
    for p, lam, r, g in zip(L.pixel_fit, L.wave_fit, resL, gL):
        if g and not in_other(p, lam, B):
            print(f'  legacy-only ID: pix {p:.2f} lam {lam:.2f} resid {r:+.3f} px')
    print('  legacy |resid| p50/max', np.round(np.percentile(np.abs(resL[gL]), [50, 100]), 3),
          '; B |resid| p50/max', np.round(np.percentile(np.abs(resB[gBf]), [50, 100]), 3))


if __name__ == '__main__':
    plot(sys.argv[1] if len(sys.argv) > 1 else '830g')
