"""
Figures for Speed_Up_Arcfit_DevSuite_Comparison.md from pr_metrics.npz
(arcfit_pr_metrics.py).  A = develop, B = develop + PR C.

Usage: python arcfit_pr_figs.py [NPZ FIGDIR]
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEV = Path('/home/xavier/Projects/PypeIt/PypeIt-development-suite')
NPZ = Path(sys.argv[1]) if len(sys.argv) > 1 else DEV / 'REDUX_OUT_PR' / 'pr_metrics.npz'
FIG = Path(sys.argv[2]) if len(sys.argv) > 2 else DEV / 'pypeitdev/speed_up/Reports/figures/pr_compare'
FIG.mkdir(parents=True, exist_ok=True)
INK, MUTED, C_A, C_B, C_HI = '#0b0b0b', '#8a8984', '#2a78d6', '#eb6834', '#7a3db8'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 160, 'savefig.bbox': 'tight',
                     'savefig.facecolor': 'white'})
d = dict(np.load(NPZ, allow_pickle=True))
ok = d['fit_a'] & d['fit_b'] & (d['rms_a'] > 0) & (d['rms_b'] > 0)
inst = d['inst']
same = ok & (d['dpix_max'] < 0.01)
small = ok & (d['dpix_max'] >= 0.01) & (d['dpix_max'] < 0.5)
large = ok & (d['dpix_max'] >= 0.5)
lab = np.array([f"{s.split('/', 1)[1] if '/' in s else s} {k} s{sl}" for s, k, sl in zip(d['setup'], d['key'], d['slit'])])


def save(fig, name):
    fig.savefig(FIG / name); plt.close(fig); print('wrote', FIG / name)


def cat_scatter(ax, x, y):
    ax.scatter(x[same], y[same], s=8, color=MUTED, alpha=0.6, lw=0, label=f'solution unchanged (<0.01 px): {same.sum()}')
    ax.scatter(x[small], y[small], s=16, color=C_A, lw=0, label=f'moved 0.01–0.5 px: {small.sum()}')
    ax.scatter(x[large], y[large], s=26, color=C_B, lw=0, label=f'moved ≥ 0.5 px: {large.sum()}')


# 1. RMS
fig, ax = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True)
cat_scatter(ax[0], d['rms_a'], d['rms_b'])
lo, hi = np.nanmin(d['rms_a'][ok]) * 0.7, np.nanmax(d['rms_a'][ok]) * 1.4
xx = np.array([lo, hi])
ax[0].plot(xx, xx, color=INK, lw=0.8); ax[0].plot(xx, 1.1 * xx, ':', color=INK, lw=0.8); ax[0].plot(xx, xx / 1.1, ':', color=INK, lw=0.8)
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlim(lo, hi); ax[0].set_ylim(lo, hi)
ax[0].set_xlabel('RMS, develop (px)'); ax[0].set_ylabel('RMS, develop + PR C (px)')
ax[0].set_title('per-slit RMS of the wavelength solution (dotted: ±10%)', loc='left')
ax[0].legend(frameon=False, loc='upper left', fontsize=8)
r = d['rms_b'][ok] / d['rms_a'][ok]
ax[1].hist(np.log10(r), bins=np.linspace(-1, 1, 81), color=C_A)
ax[1].set_yscale('log'); ax[1].set_xlabel('log10( RMS PR C / RMS develop )'); ax[1].set_ylabel('slits')
ax[1].set_title(f'{ok.sum()} slits: {np.mean(np.abs(np.log10(r)) < 0.001)*100:.1f}% within 0.2%, '
                f'{np.sum(r > 1.1)} worse > 10%, {np.sum(r < 1/1.1)} better > 10%', loc='left')
save(fig, 'rms.png')

# 2. Nlin and coverage
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
dn = d['nlin_b'][ok] - d['nlin_a'][ok]
b = np.arange(dn.min() - 0.5, dn.max() + 1.5)
ax[0].hist(dn, bins=b, color=C_A); ax[0].set_yscale('log')
ax[0].set_xlabel('Nlin(PR C) − Nlin(develop)'); ax[0].set_ylabel('slits')
ax[0].set_title(f'lines in the final fit: total {d["nlin_a"][ok].sum()} → {d["nlin_b"][ok].sum()};\n'
                f'{np.sum(dn == 0)} of {dn.size} slits unchanged', loc='left')
dc = d['cov_b'][ok] - d['cov_a'][ok]
ax[1].hist(np.clip(dc, -20, 20), bins=np.linspace(-20, 20, 81), color=C_A); ax[1].set_yscale('log')
ax[1].set_xlabel('IDs_Wave_cov(%): PR C − develop (clipped at ±20)'); ax[1].set_ylabel('slits')
ax[1].set_title(f'coverage: {np.sum(np.abs(dc) < 0.05)} unchanged (<0.05%),\n{np.sum(dc <= -1)} lose ≥ 1%, '
                f'{np.sum(dc >= 1)} gain ≥ 1%', loc='left')
save(fig, 'nlin_cov.png')

# 3. measured FWHM
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
m = ok & (d['fwhm_a'] > 0) & (d['fwhm_b'] > 0)
cat_scatter(ax[0], d['fwhm_a'], d['fwhm_b'])
xx = np.array([d['fwhm_a'][m].min() * 0.9, d['fwhm_a'][m].max() * 1.1])
ax[0].plot(xx, xx, color=INK, lw=0.8)
ax[0].set_xscale('log'); ax[0].set_yscale('log')
ax[0].set_xlabel('measured_fwhm, develop (px)'); ax[0].set_ylabel('measured_fwhm, PR C (px)')
ax[0].set_title('per-slit measured arc-line FWHM', loc='left'); ax[0].legend(frameon=False, fontsize=8)
fr = d['fwhm_b'][m] / d['fwhm_a'][m] - 1
ax[1].hist(np.clip(fr * 100, -10, 10), bins=np.linspace(-10, 10, 81), color=C_A); ax[1].set_yscale('log')
ax[1].set_xlabel('FWHM change (%, clipped at ±10)'); ax[1].set_ylabel('slits')
ax[1].set_title(f'median {np.median(fr)*100:+.3f}%; {np.sum(np.abs(fr) > 0.01)} slits change by > 1%', loc='left')
save(fig, 'fwhm.png')

# 4. wavelength-solution change
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
v = np.sort(d['dpix_max'][ok & np.isfinite(d['dpix_max'])])
ax[0].step(np.maximum(v, 1e-9), np.arange(1, v.size + 1) / v.size, where='post', color=C_A)
for t, s in [(0.01, '0.01'), (0.5, '0.5'), (5, '5')]:
    ax[0].axvline(t, color=MUTED, ls=':', lw=0.8); ax[0].text(t, 0.02, f' {s} px', color=MUTED, fontsize=8)
ax[0].set_xscale('log'); ax[0].set_xlim(1e-8, max(v.max() * 2, 10))
ax[0].set_xlabel('max |Δλ| / dispersion over the common line range (px)'); ax[0].set_ylabel('cumulative fraction of slits')
ax[0].set_title('change of the wavelength solution', loc='left')
insts = sorted(set(inst[ok]))
fam = [(i, np.sum(same & (inst == i)), np.sum(small & (inst == i)), np.sum(large & (inst == i))) for i in insts]
fam = [f for f in fam if f[2] + f[3] > 0] or fam[:10]
y = np.arange(len(fam))
ax[1].barh(y, [f[1] for f in fam], color=MUTED, label='< 0.01 px')
ax[1].barh(y, [f[2] for f in fam], left=[f[1] for f in fam], color=C_A, label='0.01–0.5 px')
ax[1].barh(y, [f[3] for f in fam], left=[f[1] + f[2] for f in fam], color=C_B, label='≥ 0.5 px')
ax[1].set_yticks(y); ax[1].set_yticklabels([f[0] for f in fam], fontsize=8); ax[1].invert_yaxis()
ax[1].set_xlabel('slits / orders'); ax[1].legend(frameon=False, fontsize=8, loc='lower right')
ax[1].set_title(f'instruments with any change ({len(fam)} of {len(insts)})', loc='left')
save(fig, 'dlambda.png')

# 5. detected arc lines (tcent)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
t = ok & (d['tc_a'] > 0)
lost = d['tc_a'][t] - d['tc_match'][t]; gained = d['tc_b'][t] - d['tc_match'][t]
ax[0].scatter(lost + np.random.default_rng(0).uniform(-0.2, 0.2, lost.size), gained + np.random.default_rng(1).uniform(-0.2, 0.2, lost.size),
              s=10, color=C_A, alpha=0.6, lw=0)
ax[0].set_xlabel('detected lines lost (develop only)'); ax[0].set_ylabel('detected lines gained (PR C only)')
ax[0].set_title(f'stored tcent: {d["tc_a"][t].sum()} → {d["tc_b"][t].sum()} lines; matched {d["tc_match"][t].sum()} '
                f'(lost {lost.sum()}, gained {gained.sum()})', loc='left', fontsize=9)
c = d['tc_dmed'][t & np.isfinite(d['tc_dmed'])]; c99 = d['tc_dp99'][t & np.isfinite(d['tc_dp99'])]
bins = np.logspace(-10, 0, 61)
ax[1].hist(np.clip(c, 1e-10, None), bins=bins, color=C_A, alpha=0.8, label='median per slit')
ax[1].hist(np.clip(c99, 1e-10, None), bins=bins, histtype='step', color=C_B, lw=1.5, label='99th percentile per slit')
ax[1].set_xscale('log'); ax[1].set_xlabel('|Δ centroid| of matched detected lines (px)'); ax[1].set_ylabel('slits')
ax[1].legend(frameon=False, fontsize=8); ax[1].set_title('centroid change of the lines both runs detect', loc='left')
save(fig, 'tcent.png')

# 6. RMS change against solution change, labelled outliers
fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
mm = ok & np.isfinite(d['dpix_max'])
x = np.maximum(d['dpix_max'][mm], 1e-8); yv = d['rms_b'][mm] / d['rms_a'][mm]
ax.scatter(x, yv, s=12, color=C_A, lw=0, alpha=0.7)
ax.set_xscale('log'); ax.set_yscale('log'); ax.axhline(1, color=INK, lw=0.6); ax.axvline(0.5, color=MUTED, ls=':', lw=0.8)
for j in np.argsort(-(np.abs(np.log(yv)) + 0.3 * np.log10(x + 1e-8)))[:14]:
    ax.annotate(lab[mm][j], (x[j], yv[j]), fontsize=7, color=INK, xytext=(4, 2), textcoords='offset points')
ax.set_xlabel('max |Δλ| / dispersion (px)'); ax.set_ylabel('RMS PR C / RMS develop')
ax.set_title('which slits changed, and did the fit get better or worse?', loc='left')
save(fig, 'rms_vs_dlambda.png')
