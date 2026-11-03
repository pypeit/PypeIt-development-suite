import warnings; warnings.simplefilter('ignore')
import numpy as np
from pypeit.core import arc
x = np.arange(64, dtype=float)
def g(c, a=100., s=1.7): return a*np.exp(-0.5*((x-c)/s)**2)
def run(name, y, pixt):
    pixt = np.asarray(pixt, dtype=int)
    try:
        new = arc.fit_arcspec(x, y, pixt, 5)
    except Exception as e:
        new = f'RAISES {type(e).__name__}: {e}'
    arc.USE_VECTORIZED_ARCFIT = False
    try:
        old = arc.fit_arcspec(x, y, pixt, 5)
    except Exception as e:
        old = f'RAISES {type(e).__name__}: {e}'
    arc.USE_VECTORIZED_ARCFIT = True
    fmt = lambda r: r if isinstance(r, str) else 'c=' + np.array2string(r[1], precision=3) + ' w=' + np.array2string(r[2], precision=3) + ' ce=' + np.array2string(r[3], precision=2)
    print(f'{name}\n   new: {fmt(new)}\n   old: {fmt(old)}')
y = np.zeros(64); y[32] = 100.; y[:8] = -5.
run('sketch edge test (pixt 0,2 all-negative; 32 delta; 63 zeros)', y, [0, 2, 32, 63])
run('zero flux window', np.zeros(64), [30])
run('negative wings', g(30.3) - 8., [30])
yn = g(30.3); yn[29] = np.nan
run('NaN pixel in window', yn, [30])
run('peaks at array edges', g(0.4) + g(63.2), [0, 1, 62, 63])
y2 = g(30.2); y2[[27, 28, 32, 33]] = -1.
run('too few positive pixels', y2, [30])
run('upward parabola (a2>=0, dip)', 50 - g(30.0, a=40), [30])
run('singular: constant positive window', np.full(64, 7.), [30])
run('amplitude overflow (huge flux)', g(30.2, a=1e300), [30])
yt = np.minimum(g(30.3, a=1e4, s=3.), 6e3)
run('saturated/flat-topped', yt, [30])
run('centre outside window (pixt off by 4)', g(34.5), [30])
run('sz_p == 0', g(30.), [])
