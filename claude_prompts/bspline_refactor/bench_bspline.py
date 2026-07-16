"""
Execution-time comparison: old :class:`pypeit.bspline.bspline.bspline` vs.
refactored :class:`pypeit.bspline.refactor.BSpline` / :class:`BSpline2D`.

Both classes are run on identical data.  For the 1D case the old class is
initialised with ``npoly=1``; for the 2D case both ``npoly`` and ``funcname``
are set before fitting.

Usage::

    python bench_bspline.py            # default table
    python bench_bspline.py --nrep 20  # more repetitions for smoother timing
"""

import argparse
import time

import numpy as np

from pypeit.bspline.bspline import bspline
from pypeit.bspline.refactor import BSpline, BSpline2D, Knots


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _time(func, nrep):
    """Return median wall-time of *nrep* calls to ``func()`` in seconds."""
    times = []
    for _ in range(nrep):
        t0 = time.perf_counter()
        func()
        times.append(time.perf_counter() - t0)
    return np.median(times)


def _header(title):
    width = 72
    print()
    print('=' * width)
    print(f'  {title}')
    print('=' * width)
    print(f'  {"case":<30}  {"old (ms)":>9}  {"new (ms)":>9}  {"speedup":>9}')
    print('-' * width)


def _row(label, t_old, t_new):
    speedup = t_old / t_new if t_new > 0 else float('inf')
    print(f'  {label:<30}  {t_old*1e3:>9.2f}  {t_new*1e3:>9.2f}  {speedup:>8.2f}x')


# ---------------------------------------------------------------------------
# 1D benchmarks
# ---------------------------------------------------------------------------

def bench_1d(configs, nrep):
    """Benchmark BSpline.fit (1D) for each (N, nbkpts) configuration."""
    _header('1D fit  [npoly=1]')
    rng = np.random.default_rng(0)

    for N, nbkpts in configs:
        x = np.sort(rng.uniform(0, 10, N))
        y = np.sin(x) + 0.1 * rng.standard_normal(N)
        invvar = np.ones(N)

        def run_old():
            spl = bspline(x=x, nord=4, npoly=1, nbkpts=nbkpts)
            spl.fit(x, y, invvar)

        def run_new():
            spl = BSpline(x=x, knots=Knots(count=nbkpts), nord=4)
            spl.fit(x, y, invvar)

        t_old = _time(run_old, nrep)
        t_new = _time(run_new, nrep)
        _row(f'N={N:>6d}  nknots={nbkpts:>4d}', t_old, t_new)

    print()


# ---------------------------------------------------------------------------
# 2D benchmarks
# ---------------------------------------------------------------------------

def bench_2d(configs, nrep, funcname='legendre'):
    """Benchmark BSpline2D.fit for each (N, nbkpts, npoly) configuration."""
    _header(f'2D fit  [funcname={funcname}]')
    rng = np.random.default_rng(1)

    for N, nbkpts, npoly in configs:
        x = np.sort(rng.uniform(0, 10, N))
        x2 = rng.uniform(0, 1, N)
        y = np.sin(x) * (1 + 0.3 * x2)
        invvar = np.ones(N)

        def run_old(npoly=npoly):
            spl = bspline(x=x, nord=4, npoly=npoly, nbkpts=nbkpts,
                          funcname=funcname)
            spl.xmin = 0.0
            spl.xmax = 1.0
            spl.fit(x, y, invvar, x2=x2)

        def run_new(npoly=npoly):
            spl = BSpline2D(x=x, knots=Knots(count=nbkpts), nord=4)
            spl.fit(x, y, invvar, basis=funcname, npoly=npoly, basis_x=x2,
                    xmin=0.0, xmax=1.0)

        t_old = _time(run_old, nrep)
        t_new = _time(run_new, nrep)
        _row(f'N={N:>6d}  nknots={nbkpts:>4d}  npoly={npoly}', t_old, t_new)

    print()


# ---------------------------------------------------------------------------
# Sigma-clipping loop benchmarks (cache benefit)
# ---------------------------------------------------------------------------

def bench_sigma_clip(configs, nrep):
    """Benchmark repeated fit() calls on fixed x (sigma-clipping pattern).

    The new implementation caches the design matrix after the first call;
    the old implementation recomputes ``action`` on every ``workit`` call.
    """
    _header('Sigma-clipping loop  [10 iterations, same x, varying invvar]')
    rng = np.random.default_rng(2)
    niter = 10

    for N, nbkpts in configs:
        x = np.sort(rng.uniform(0, 10, N))
        y = np.sin(x) + 0.1 * rng.standard_normal(N)
        invvar_base = np.ones(N)
        # Vary invvar between iterations (masking random points) to simulate
        # sigma-clipping while keeping x fixed.
        invvar_list = [invvar_base.copy() for _ in range(niter)]
        for iv in invvar_list:
            iv[rng.integers(0, N, N // 20)] = 0.0

        def run_old():
            spl = bspline(x=x, nord=4, npoly=1, nbkpts=nbkpts)
            for iv in invvar_list:
                spl.fit(x, y, iv)

        def run_new():
            spl = BSpline(x=x, knots=Knots(count=nbkpts), nord=4)
            for iv in invvar_list:
                spl.fit(x, y, iv)

        t_old = _time(run_old, nrep)
        t_new = _time(run_new, nrep)
        _row(f'N={N:>6d}  nknots={nbkpts:>4d}', t_old, t_new)

    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nrep', type=int, default=10,
                        help='Number of repetitions per timing (default 10)')
    args = parser.parse_args()
    nrep = args.nrep

    print(f'\nBenchmark: old bspline vs. new BSpline  (nrep={nrep})')
    print('Timing = median wall-time over all repetitions.')

    bench_1d(
        configs=[
            (500,   20),
            (500,  100),
            (2000,  20),
            (2000, 100),
            (10000, 20),
            (10000, 100),
        ],
        nrep=nrep,
    )

    bench_2d(
        configs=[
            (500,   20,  3),
            (500,  100,  3),
            (2000,  20,  3),
            (2000, 100,  3),
            (2000,  20,  6),
            (10000, 20,  3),
            (10000, 100, 3),
        ],
        nrep=nrep,
    )

    bench_sigma_clip(
        configs=[
            (2000,  20),
            (2000, 100),
            (10000, 20),
        ],
        nrep=nrep,
    )


if __name__ == '__main__':
    main()
