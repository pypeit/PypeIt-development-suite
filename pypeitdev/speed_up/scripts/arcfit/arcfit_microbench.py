"""
Micro-benchmark of pypeit.core.arc.fit_arcspec: vectorized vs the legacy
per-line curve_fit, on ~2000 isolated synthetic lines (the same generator as
pypeit/tests/test_arc.py) and on the Paranal sky spectrum.
"""
import time, warnings
warnings.simplefilter('ignore')
import numpy as np
from pypeit.core import arc, skyspec

def synth(nspec, spacing=20., fwhm=4., noise=1., seed=42):
    rng = np.random.default_rng(seed)
    x = np.arange(nspec, dtype=float)
    cen = np.arange(20., nspec - 20., spacing); cen += rng.uniform(-2., 2., cen.size)
    amp = 10.**rng.uniform(1.5, 3.5, cen.size); sig = fwhm/2.355*rng.uniform(0.9, 1.1, cen.size)
    y = np.zeros(nspec)
    for a, c, s in zip(amp, cen, sig):
        lo, hi = int(c) - 30, int(c) + 31
        y[lo:hi] += a*np.exp(-0.5*((x[lo:hi] - c)/s)**2)
    return x, y + rng.normal(scale=noise, size=nspec), np.round(cen).astype(int)

def bench(f, n=5):
    best = 1e9
    for _ in range(n):
        t = time.perf_counter(); f(); best = min(best, time.perf_counter() - t)
    return best

cases = {'synthetic': synth(40960)}
sky = skyspec.load_sky_spectrum('paranal_sky.fits').flux
calls = []; orig = arc.fit_arcspec
arc.fit_arcspec = lambda *a: (calls.append(a), orig(*a))[1]
arc.detect_lines(sky); arc.fit_arcspec = orig
cases['paranal sky'] = calls[0][:3]
for name, (x, y, pixt) in cases.items():
    tv = bench(lambda: orig(x, y, pixt, 5))
    arc.USE_VECTORIZED_ARCFIT = False
    tl = bench(lambda: orig(x, y, pixt, 5), n=2)
    arc.USE_VECTORIZED_ARCFIT = True
    print(f'{name:12s} {pixt.size:5d} lines: legacy {tl*1e3:8.1f} ms  vectorized {tv*1e3:7.1f} ms  '
          f'speedup {tl/tv:6.1f}x  ({tl/pixt.size*1e6:.0f} vs {tv/pixt.size*1e6:.1f} us/line)')
