"""
Figures for the Q17 traces in Speed_Up_Arcfit_DevSuite_Comparison.md (C.md,
"PR" task 2).  Reads the harness outputs in REDUX_OUT_PR/q17 (scripts in
scripts/arcfit/q17) and the develop / arcfitBxa dev-suite outputs.

    q17_maskdesign.png  extra DEIMOS edge in peak_trace (test_add_missing_slits)
    q17_gnirs_fwhm.png  GNIRS 32_SB_SXD order 8: FWHM map, A vs B
    q17_gnirs_tilts.png GNIRS order 8: tilt cliff vs FWHM, spec_order 4 vs 2, spec1d
    q17_mmirs.png       MMIRS nep_as1_MOS slit 18: tripled template, tied xcorr peaks
    q17_ldss3.png       LDSS3 VPH-RED_mos slit 21: box slit, wrong on both branches

Usage: python arcfit_q17_figs.py
"""
import pickle
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits

DEV = Path('/home/xavier/Projects/PypeIt/PypeIt-development-suite')
OUT = DEV / 'REDUX_OUT_PR'
Q = OUT / 'q17'
FIG = DEV / 'pypeitdev/speed_up/Reports/figures/q17'
FIG.mkdir(parents=True, exist_ok=True)
TEMPLATES = Path('/mnt/tank/Astronomy/PypeIt/PypeIt/pypeit/data/arc_lines/reid_arxiv')

INK, MUTED, C_A, C_B, C_FIX = '#0b0b0b', '#8a8984', '#2a78d6', '#eb6834', '#2e9e5b'
plt.rcParams.update({'font.size': 9.5, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 160, 'savefig.bbox': 'tight',
                     'savefig.facecolor': 'white', 'legend.frameon': False,
                     'axes.titlesize': 10, 'axes.titlelocation': 'left'})


def save(fig, name):
    fig.savefig(FIG / name)
    plt.close(fig)
    print('wrote', FIG / name)


def gauss(x, a, c, s):
    return a * np.exp(-0.5 * ((x - c) / s) ** 2)


# --------------------------------------------------------------------------
# 1. maskdesign: B keeps a narrow Sobel spike that curve_fit rejects (sigma<0)
# --------------------------------------------------------------------------
calls = pickle.load(open(Q / 'mdm_calls.pkl', 'rb'))
res = pickle.load(open(Q / 'mdm_res.pkl', 'rb'))
spec = calls[1][0]
tA, tB = res[(1, 'A')], res[(1, 'B')]
jA = np.argmin(np.abs(tA[2] - 298.86))
jB = np.argmin(np.abs(tB[2] - 298.86))
edev = np.load(Q / 'mdm_dev.npy')[0]
eB = np.load(Q / 'mdm_B.npy')[0]

fig, ax = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
x = np.arange(spec.size)
sl = slice(280, 320)
ax[0].plot(x[sl], spec[sl], 'o-', color=INK, ms=3, lw=1, label='rectified, collapsed Sobel image')
xx = np.linspace(292, 306, 300)
ax[0].plot(xx, gauss(xx, tB[0][jB], tB[2][jB], abs(tB[3][jB])), color=C_B, lw=2,
           label=f'B fit: centre {tB[2][jB]:.2f}, σ = +{tB[3][jB]:.2f} px → kept')
ax[0].plot(xx, gauss(xx, tA[0][jA], tA[2][jA], abs(tA[3][jA])), '--', color=C_A, lw=1.5,
           label=f'A / curve_fit: same Gaussian, σ = {tA[3][jA]:.2f} px → twid > 0 fails')
good = tB[5]
medw = np.median(tB[3][good])
ax[0].set_title(f'DEIMOS 1200G_Cooper DET01, peak_trace: a narrow spike next to a deep trough\n'
                f'(σ/median σ = {tB[3][jB] / medw:.2f}; the 0.5 floor is off in peak_trace)')
ax[0].set_xlabel('spatial pixel (rectified)')
ax[0].set_ylabel('Sobel signal')
ax[0].legend(loc='lower left', fontsize=8)

n = 16
ax[1].plot(np.arange(n), edev[:n], 'o', color=C_A, ms=7, label=f'develop: {edev.size} edges')
ax[1].plot(np.arange(n), eB[:n], 's', mfc='none', color=C_B, ms=9, mew=1.5, label=f'PR C (B): {eB.size} edges')
iextra = int(np.argmin(np.abs(eB - 303.4)))
ax[1].annotate('PR-only edge at 303.4 px\n(from the spike)', (iextra, eB[iextra]),
               (iextra + 2.5, 230), arrowprops=dict(arrowstyle='->', color=C_B), color=C_B, fontsize=8.5)
for i, (v, lbl) in enumerate([(7, '[7] 458'), (9, '[9] 496'), (10, '[10] 561')]):
    ax[1].annotate(f'test expects {lbl}', (v, edev[v]), (v + 2.5, 420 - 45 * i), fontsize=8, color=C_A,
                   arrowprops=dict(arrowstyle='->', color=C_A, lw=0.8))
ax[1].set_xlabel('trace index')
ax[1].set_ylabel('edge position at the PCA reference row (px)')
ax[1].set_title('Every later index shifts by one, so the hard-coded\n'
                'indices in test_add_missing_slits pick different edges')
ax[1].legend(loc='upper left')
save(fig, 'q17_maskdesign.png')

# --------------------------------------------------------------------------
# 2. GNIRS order 8: the FWHM map
# --------------------------------------------------------------------------
g7 = np.load(Q / 'gn7.npy', allow_pickle=True).item()
a, b = g7['A'], g7['B']
onlyA = np.array([not np.any((b['samp'] == a['samp'][i]) & (np.abs(b['cen'] - a['cen'][i]) < 0.5))
                  for i in range(a['cen'].size)])
rows = np.arange(a['curve'].size)
fig, ax = plt.subplots(figsize=(8.5, 4.6), constrained_layout=True)
ax.axvspan(0, 512.5, color=MUTED, alpha=0.12, lw=0)
ax.text(10, 5.3, 'not on the detector\nfor order 8', color=MUTED, fontsize=8.5, va='top')
ax.scatter(b['cen'], b['fwhm'], s=14, color=MUTED, lw=0, label='lines kept by both (10 spatial samples)')
ax.scatter(a['cen'][onlyA], a['fwhm'][onlyA], s=40, color=C_A, marker='v',
           label=f'{onlyA.sum()} sub-resolution spikes kept by A only (FWHM 0.6–1.6 px)')
ax.plot(rows, a['curve'], color=C_A, lw=2, label='develop (A) linear fit')
ax.plot(rows, b['curve'], color=C_B, lw=2, label='PR C (B, σ floor) linear fit')
r0 = rows.size // 2
ax.axvline(r0, color=INK, lw=0.8, ls=':')
ax.plot([r0], [a['val']], 'o', color=C_A, ms=8)
ax.plot([r0], [b['val']], 'o', color=C_B, ms=8)
ax.annotate(f'evaluated at row nspec//2 = {r0}:\n develop {a["val"]:.2f} px, PR {b["val"]:.2f} px',
            (r0, (a['val'] + b['val']) / 2), (r0 - 330, 1.0), fontsize=8.5,
            arrowprops=dict(arrowstyle='->', color=INK, lw=0.8))
ax.axhspan(4.02, 4.05, color=C_B, alpha=0.25, lw=0)
ax.text(1015, 4.035, 'tilt cliff\n(Fig. q17_gnirs_tilts)', ha='right', va='center', fontsize=8, color=C_B)
ax.set_ylim(0, 5.5)
ax.set_xlim(0, rows.size)
ax.set_xlabel('spectral row')
ax.set_ylabel('arc-line FWHM (px)')
ax.set_title('GNIRS 32_SB_SXD order 8: map_fwhm (fwhm_spec_order = 1).  The order-8 FWHM is an extrapolation\n'
             'to row 511; develop\'s value is held below 4 px by spikes that the σ floor now removes')
ax.legend(loc='lower right', fontsize=8)
save(fig, 'q17_gnirs_fwhm.png')

# --------------------------------------------------------------------------
# 3. GNIRS order 8: the tilt cliff and its fix
# --------------------------------------------------------------------------
scan = {}
for line in open(Q / 'gn6_scan.txt'):
    m = re.search(r'so=(\d)\].*fwhm8=([0-9.]+).*max\|tilt-row\| ([0-9.]+)', line)
    if m:
        scan.setdefault(int(m.group(1)), []).append((float(m.group(2)), float(m.group(3))))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3), constrained_layout=True)
for so, col, lbl in [(4, C_B, 'spec_order = 4 (GNIRS default)'), (2, C_FIX, 'spec_order = 2 (proposed)')]:
    v = np.array(sorted(scan[so]))
    ax[0].plot(v[:, 0], v[:, 1], 'o-', color=col, label=lbl)
ax[0].axvline(a['val'], color=C_A, ls=':', lw=1)
ax[0].axvline(b['val'], color=C_B, ls=':', lw=1)
ax[0].text(a['val'], 30, ' develop\n 3.73', color=C_A, fontsize=8)
ax[0].text(b['val'], 30, ' PR\n 4.16', color=C_B, fontsize=8)
ax[0].set_yscale('log')
ax[0].set_ylim(0.05, 150)
ax[0].set_xlabel('FWHM passed to the order-8 tilt tracing (px)')
ax[0].set_ylabel('max |tilt − row| along the order centre (px)')
ax[0].set_title('(a) a cliff between 4.02 and 4.05 px')
ax[0].legend(loc='center left', fontsize=8)
for f, so, col, lbl in [('gn6_fwhm3.73_so4.npy', 4, C_A, 'FWHM 3.73, spec_order 4 (develop)'),
                        ('gn6_fwhm4.16_so4.npy', 4, C_B, 'FWHM 4.16, spec_order 4 (PR)'),
                        ('gn6_fwhm4.16_so2.npy', 2, C_FIX, 'FWHM 4.16, spec_order 2 (proposed)')]:
    r, v, m = np.load(Q / f)
    m = m.astype(bool)
    ax[1].plot(r[m], v[m] - r[m], color=col, lw=1.8, label=lbl)
for c, lbl in [(552.1, '552'), (558.9, '559'), (723.2, '723')]:
    ax[1].axvline(c, color=MUTED, lw=0.6, ls='--')
ax[1].text(560, 50, 'bluest tilt line:\n 552 used at 3.73 px;\n 559 found but rejected at 4.16 px,\n next usable line is at 723',
           fontsize=8, color=INK)
ax[1].set_xlabel('spectral row')
ax[1].set_ylabel('tilt − row at the order centre (px)')
ax[1].set_title('(b) 4th-order fit extrapolated over rows 513–723')
ax[1].legend(loc='lower right', fontsize=8)
sp = 'gemini_gnirs_echelle/32_SB_SXD/Science/spec1d_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
for run, col, lbl in [('develop', C_A, 'develop'), ('arcfitBxa', C_B, 'PR C')]:
    with fits.open(OUT / run / sp) as h:
        for hdu in h[1:]:
            if 'ORDER0008' in hdu.name:
                w = hdu.data['OPT_WAVE']
                pix = np.arange(w.size)
                m = w > 0
                ax[2].plot(pix[m], w[m], color=col, lw=1.5,
                           label=f'{lbl}: {w[m].min():.0f}–{w[m].max():.0f} Å, '
                                 f'{np.sum(np.diff(w[m]) <= 0)} non-increasing steps')
ax[2].set_xlim(505, 640)
ax[2].set_ylim(8150, 8650)
ax[2].set_xlabel('spec1d pixel')
ax[2].set_ylabel('OPT_WAVE (Å)')
ax[2].set_title('(c) HIP62745 spec1d, order 8, blue end:\ntest_specutils::test_load_spec1d fails on PR')
ax[2].legend(loc='upper left', fontsize=8)
save(fig, 'q17_gnirs_tilts.png')

# --------------------------------------------------------------------------
# 4. MMIRS: tripled template, tied cross-correlation peaks, fix
# --------------------------------------------------------------------------
t = fits.open(TEMPLATES / 'mmt_mmirs_J_zJ.fits')[1].data
tw, tf = np.asarray(t['wave']).ravel(), np.asarray(t['flux']).ravel()
cdev = pickle.load(open(Q / 'corrcaps_dev.pkl', 'rb'))[0][0]
cpr = pickle.load(open(Q / 'corrcaps_prA.pkl', 'rb'))[0][0]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3), constrained_layout=True)
ax[0].plot(tw, color=INK, lw=1.2)
ax2 = ax[0].twinx()
ax2.plot(tf, color=MUTED, lw=0.5, alpha=0.7)
ax2.set_yticks([])
ax2.spines['right'].set_visible(False)
for s in [2040, 4080]:
    ax[0].axvline(s, color=C_B, lw=1, ls='--')
ax[0].set_xlabel('template pixel')
ax[0].set_ylabel('template wavelength (Å)')
ax[0].set_title('(a) mmt_mmirs_J_zJ.fits = 3 identical copies of a 2040-px\n'
                'spectrum (also HK_zJ, K3000_Kspec)')
lags = np.arange(cdev.size)
for c, col, lbl in [(cdev, C_A, 'develop'), (cpr, C_B, 'PR C')]:
    ax[1].plot(lags, c, color=col, lw=1, label=lbl, alpha=0.9)
for pk in [3832.5, 5872.6, 7912.6]:
    ax[1].axvline(pk, color=MUTED, lw=0.5, ls=':')
ax[1].annotate('develop picks 7912.6\n(continuum-subtracted\nheight 0.137161 vs 0.137154)', (7912.6, 0.253),
               (8400, 0.20), fontsize=8, color=C_A, arrowprops=dict(arrowstyle='->', color=C_A, lw=0.8))
ax[1].annotate('PR picks 5872.6\n(0.136362 vs 0.136351)', (5872.6, 0.253), (4050, 0.03), fontsize=8, color=C_B,
               arrowprops=dict(arrowstyle='->', color=C_B, lw=0.8))
ax[1].text(3050, 0.235, 'the two curves overlap\n(develop is under PR C)', fontsize=7.5, color=MUTED)
ax[1].set_xlim(3000, 9500)
ax[1].set_xlabel('lag index')
ax[1].set_ylabel('normalised cross-correlation')
ax[1].set_title('(b) slit 18: global cross-correlation with the template;\n'
                'two peaks of identical height 0.253480 (max |Δ| 3e-7)')
ax[1].legend(loc='upper left', fontsize=8)
for f, col, ls, lbl in [('ft_0_mmt_mmirs_nep_as1_MOS_18_native_develop_develop.npy', C_A, '-', 'develop'),
                        ('ft_0_mmt_mmirs_nep_as1_MOS_18_native_arcfitBxa_arcfitBxa.npy', C_B, '-', 'PR C'),
                        ('ft_1_mmt_mmirs_nep_as1_MOS_18_native_arcfitBxa_arcfitBxa.npy', C_FIX, '--',
                         'PR C, one-copy template (proposed)'),
                        ('ft_1_mmt_mmirs_nep_as1_MOS_18_native_develop_develop.npy', INK, ':',
                         'develop, one-copy template')]:
    r = np.load(Q / f, allow_pickle=True).item()
    ax[2].plot(r['wave'], color=col, ls=ls, lw=1.8, label=f'{lbl}: RMS {r["rms"]:.2f} px, {r["nlin"]} lines')
ax[2].set_xlabel('spectral pixel')
ax[2].set_ylabel('wavelength (Å)')
ax[2].set_title('(c) nep_as1_MOS slit 18 solution: the blue end of the\n'
                'PR window is filled with the red end of copy 1')
ax[2].legend(loc='lower right', fontsize=8)
save(fig, 'q17_mmirs.png')

# --------------------------------------------------------------------------
# 5. LDSS3 slit 21: a box slit that fails on both branches
# --------------------------------------------------------------------------
import sys
sys.path.insert(0, '/mnt/tank/Astronomy/PypeIt/worktrees/develop')
from pypeit.wavecalib import WaveCalib
from pypeit.slittrace import SlitTraceSet
setup = 'magellan_ldss3/VPH-RED_mos'
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3), constrained_layout=True)
sl = SlitTraceSet.from_file(next((OUT / 'develop' / setup / 'Calibrations').glob('Slits_*')), chk_version=False)
width = np.median(sl.right_init - sl.left_init, axis=0)
box = width < 32
for run, col, mk, lbl in [('develop', C_A, 'o', 'develop'), ('arcfitBxa', C_B, 's', 'PR C')]:
    w = WaveCalib.from_file(OUT / run / setup / 'Calibrations/WaveCalib_A_1_DET01.fits', chk_version=False)
    ns = w.arc_spectra.shape[0]
    xx = np.arange(ns) / (ns - 1)
    w0 = np.array([f.pypeitfit.eval(np.array([0.]))[0] for f in w.wv_fits])
    ax[0].scatter(np.arange(w0.size)[~box], w0[~box], color=col, marker=mk, s=30,
                  facecolors='none' if run != 'develop' else col, label=f'{lbl}, 38-px slits')
    ax[0].scatter(np.arange(w0.size)[box], w0[box], color=col, marker='*', s=110,
                  label=f'{lbl}, 27-px box slits (0, 1, 5, 21)')
    for i, ls in [(20, ':'), (21, '-'), (22, '--')]:
        wv = w.wv_fits[i].pypeitfit.eval(xx)
        if i == 21 or run == 'develop':
            ax[1].plot(np.arange(ns), wv, color=col if i == 21 else MUTED, ls=ls, lw=1.6 if i == 21 else 1,
                       label=f'{lbl} slit {i} ({len(w.wv_fits[i].pixel_fit)} lines)')
ax[0].set_ylim(4000, 9500)
ax[0].set_xlabel('slit index')
ax[0].set_ylabel('wavelength at pixel 0 (Å)')
ax[0].set_title('(a) VPH-RED_mos A_1: slit 21 is a box slit; it starts ~1 100 Å\n'
                '(develop) and ~1 700 Å (PR) redder than its neighbours')
ax[0].legend(loc='upper left', fontsize=7.5, ncol=2)
ax[1].set_ylim(4000, 12500)
ax[1].set_xlabel('spectral pixel')
ax[1].set_ylabel('wavelength (Å)')
ax[1].set_title('(b) slit 21 vs its neighbours: both solutions turn over\n'
                '(develop IDs_Wave_cov = 102%, PR 1.7%)')
ax[1].legend(loc='lower right', fontsize=8)
save(fig, 'q17_ldss3.png')
