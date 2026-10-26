"""
Task 6 item 1: record every detect_lines call made by autoid.full_template for
keck_deimos/600ZD_M_6500 MSC03 slit 11 (spat_id 1063), legacy vs vectorized
fit_arcspec, and report the first call whose outputs differ, the differing lines,
and their fit_arcspec parameters under both paths.
Usage: python arcfit_msc03_trace.py [REDUX_OUT_C]
"""
import sys, warnings
from pathlib import Path
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.wavecalib import WaveCalib
from pypeit.images.buildimage import ArcImage
from pypeit.par.pypeitpar import WavelengthSolutionPar
from pypeit.core import arc
from pypeit.core.wavecal import autoid

root = Path(sys.argv[1] if len(sys.argv) > 1 else
            '/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C')
rel = 'keck_deimos/600ZD_M_6500/Calibrations/WaveCalib_A_0_MSC03.fits'
islit = 11
wc = WaveCalib.from_file(str(root / 'A' / rel), chk_version=False)
msarc = ArcImage.from_file(str(root / 'A' / rel.replace('WaveCalib', 'Arc')), chk_version=False)
nlc = np.max([d.nonlinear_counts() for d in msarc.detector.detectors])
par = WavelengthSolutionPar.from_dict(dict(wc.par))
spec = wc.arc_spectra[:, islit]
fwhm = wc.wv_fits[islit].fwhm

orig_dl, orig_fa = arc.detect_lines, arc.fit_arcspec
rec = {}
def run(vec):
    calls = []
    def dl(*a, **k):
        out = orig_dl(*a, **k)
        calls.append(dict(args=a, kw=k, out=out))
        return out
    def fa(x, y, pixt, fitp):
        out = orig_fa(x, y, pixt, fitp)
        calls[-1:] and None
        fa.last = (x, y, pixt, fitp)
        return out
    arc.detect_lines = dl
    arc.USE_VECTORIZED_ARCFIT = vec
    try:
        out, _ = autoid.full_template(spec, par['lamps'], par, [0], 1, 1, slit_ids=[1063],
                                      measured_fwhms=[fwhm], nonlinear_counts=nlc,
                                      nsnippet=par['nsnippet'], x_percentile=par['cc_percent_ceil'])
    finally:
        arc.detect_lines = orig_dl
        arc.USE_VECTORIZED_ARCFIT = True
    return calls, out['0']

cl, fl = run(False)
cv, fv = run(True)
print(f'{len(cl)} detect_lines calls (legacy), {len(cv)} (vectorized)')
for i, (a, b) in enumerate(zip(cl, cv)):
    sa, sb = a['args'][0], b['args'][0]
    same_in = np.array_equal(sa, sb)
    ga, gb = a['out'][5], b['out'][5]
    pa, pb = a['out'][2], b['out'][2]
    same_out = np.array_equal(ga, gb) and pa.size == pb.size and np.allclose(pa[ga], pb[gb], atol=1e-4)
    print(f' call {i}: nspec {sa.size} kw {sorted(a["kw"].keys())} input same {same_in}; n peaks {pa.size}/{pb.size}; '
          f'n good {ga.size}/{gb.size}; outputs same {same_out}')
    if not same_out and same_in:
        # Re-run fit_arcspec on the very same peaks to characterize the differing lines
        s, kw = sa, a['kw']
        fwhm_i = kw.get('fwhm', 4.0)
        ffw = kw.get('fit_frac_fwhm', 1.25)
        arcsub = a['out'][6]
        x = np.arange(s.size, dtype=float)
        thresh_kw = {}
        pixt = np.round(pa).astype(int)  # placeholder, replaced below
        # recover pixt exactly from detect_peaks as in detect_lines
        tl, tc = a['out'], b['out']
        only_l = sorted(set(ga) - set(gb)); only_v = sorted(set(gb) - set(ga))
        nfitpix = np.round(ffw * fwhm_i).astype(int)
        fwhm_max = kw.get('max_frac_fwhm', 3.0) * fwhm_i
        print(f'   fwhm {fwhm_i:.6f} nfitpix {nfitpix} (half-width {(nfitpix + nfitpix % 2)//2}) '
              f'width cut twid < {fwhm_max/2.35:.4f}; nonlinear {kw.get("nonlinear_counts")}')
        for tag, idxs in [('legacy only', only_l), ('vectorized only', only_v)]:
            for j in idxs:
                print(f'   {tag}: idx {j}: legacy ampl_true {tl[0][j]:.1f} tcent {tl[2][j]:.4f} twid {tl[3][j]:.4f} '
                      f'cerr {tl[4][j]:.3g} | vect tcent {tc[2][j]:.4f} twid {tc[3][j]:.4f} cerr {tc[4][j]:.3g} '
                      f'nsig {tl[7][j]:.1f}')
        rec[i] = dict(arcsub=arcsub, only_l=only_l, only_v=only_v, tl=tl, tc=tc, nfitpix=nfitpix)
np.save('/tmp/claude-1000/-mnt-tank-Astronomy-PypeIt-PypeIt-development-suite/2e03cf89-b763-496b-a135-c28efe22343e/scratchpad/msc03_trace.npy',
        rec, allow_pickle=True)

# ---- All peaks (not just the 'good' ones): get_xcorr_arc and xcorr_shift use
# every peak whose tcent != -999, so compare tcent for every peak.
print('\nPer-call differences over ALL peaks (|dtcent| > 1e-3 px or -999 status differs):')
for i, (a, b) in enumerate(zip(cl, cv)):
    ta, tb = a['out'][2], b['out'][2]
    wa, wb = a['out'][3], b['out'][3]
    if ta.size != tb.size:
        print(f' call {i}: different number of peaks'); continue
    d = np.abs(ta - tb)
    bad = np.where((d > 1e-3) | ((ta == -999.) != (tb == -999.)))[0]
    if a['args'][0].size != b['args'][0].size:
        din = np.nan
    else:
        din = np.max(np.abs(a['args'][0] - b['args'][0]))
    print(f' call {i}: max|d input| {din:.3g}; peaks differing: {bad.size}')
    good_l = set(a['out'][5]); 
    for j in bad[:20]:
        print(f'    peak {j}: legacy tcent {ta[j]:.4f} twid {wa[j]:.4f} ampl_true {a["out"][0][j]:.1f} '
              f'| vect tcent {tb[j]:.4f} twid {wb[j]:.4f} | legacy-good {j in good_l} nsig {a["out"][7][j]:.1f}')

# ---- Is the call-9 difference due to the fit method or the 1e-8 input difference?
# Re-fit the correlation peak of call 9 with both fit_arcspec paths on the SAME input.
from pypeit.core.arc import detect_peaks
print('\nCall 9 correlation-peak fit on identical input (legacy input):')
a = cl[9]
for lab, src in [('legacy input', cl[9]), ('vect input', cv[9])]:
    s = src['args'][0]
    arcsub = src['out'][6]
    x = np.arange(s.size, dtype=float)
    # nfind=1 kept the most significant peak; refit all peaks and pick the same one
    for vec in [False, True]:
        arc.USE_VECTORIZED_ARCFIT = vec
        out = orig_dl(s, **src['kw'])
        arc.USE_VECTORIZED_ARCFIT = True
        print(f'  {lab}, vectorized={vec}: tcent {out[2][0]:.6f} twid {out[3][0]:.6f} ampl {out[1][0]:.5f}')
# What does the vectorized path do with this line: polish converged? width ratio?
s = cl[9]['args'][0]; arcsub = cl[9]['out'][6]
pk = int(np.round(cl[9]['out'][2][0]))
cand = detect_peaks(arcsub, mph=None, mpd=1)
pixt = cand[np.argmin(np.abs(cand - pk))]
x = np.arange(arcsub.size, dtype=float)
nfit = int(np.round(1.5 * 5.0))
for vec in [False, True]:
    arc.USE_VECTORIZED_ARCFIT = vec
    r = orig_fa(x, arcsub, np.array([pixt]), nfit)
    print(f'  fit_arcspec pixt={pixt} nfitpix={nfit} vectorized={vec}: ampl {r[0][0]:.6f} cent {r[1][0]:.6f} sig {r[2][0]:.6f} cerr {r[3][0]:.3g}')
arc.USE_VECTORIZED_ARCFIT = True
half = (nfit + nfit % 2) // 2
print(f'  window half-width {half}; sig/half = {r[2][0]/half:.2f} (>= {arc.ARCFIT_MIN_WIDTH} -> kept by vectorized path)')
# Objective (sum of squares) at both solutions
u = np.arange(pixt - half, pixt + half + 1); yy = arcsub[u]
for vec in [False, True]:
    arc.USE_VECTORIZED_ARCFIT = vec
    r = orig_fa(x, arcsub, np.array([pixt]), nfit)
    m = r[0][0] * np.exp(-0.5 * ((u - r[1][0]) / r[2][0])**2)
    print(f'  vectorized={vec}: chi2 (unweighted SSE) = {np.sum((yy - m)**2):.6e}')
arc.USE_VECTORIZED_ARCFIT = True

# ---- xcorr_shift_stretch outputs in each run
from pypeit.core.wavecal import wvutils
orig_xs = wvutils.xcorr_shift_stretch
for vec in [False, True]:
    res = []
    def xs(*a, **k):
        o = orig_xs(*a, **k); res.append(o[:4]); return o
    wvutils.xcorr_shift_stretch = xs
    run(vec)
    wvutils.xcorr_shift_stretch = orig_xs
    for o in res:
        print(f'  vectorized={vec}: xcorr_shift_stretch success {o[0]} shift {o[1]:.5f} stretch {o[2]!r} ccorr {o[3]!r}')

# ---- Size of the centroid differences feeding the chain, and the full
# xcorr_shift_stretch outputs (incl. final cc and initial shift) for snippet 0.
print('\nmax |dtcent| over all peaks with tcent != -999, per call:')
print('  ' + ', '.join(f'{i}:{np.max(np.abs(a["out"][2]-b["out"][2])[(a["out"][2]!=-999)]) if np.any(a["out"][2]!=-999) else 0:.1e}'
                       for i, (a, b) in enumerate(zip(cl, cv))))
for vec in [False, True]:
    res = []
    def xs(*a, **k):
        o = orig_xs(*a, **k); res.append((o, a)); return o
    wvutils.xcorr_shift_stretch = xs
    run(vec)
    wvutils.xcorr_shift_stretch = orig_xs
    (o, args) = res[0]
    y1, y2 = args[0], args[1]
    print(f'  snippet 0, vectorized={vec}: success {o[0]} shift {o[1]:.4f} stretch {o[2]:.6f} stretch2 {o[3]:.2e} '
          f'cc {o[4]:.5f} shift_init {o[5]:.5f} cc_init {o[6]:.5f}')

# ---- Why does a 1e-8 change of the correlation give a 6.6e-3 px change of its peak?
a, b = cl[9], cv[9]
print('\ncall 9: max|d corr input| %.2e, max|d cont-subtracted corr| %.2e, cont_subtract=%s'
      % (np.max(np.abs(a['args'][0]-b['args'][0])), np.max(np.abs(a['out'][6]-b['out'][6])), a['kw'].get('cont_subtract')))
