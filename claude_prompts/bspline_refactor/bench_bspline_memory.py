"""
Memory-usage comparison: old :class:`pypeit.bspline.bspline.bspline` vs.
refactored :class:`pypeit.bspline.refactor.BSpline` / :class:`BSpline2D`.

Two metrics are reported for each configuration:

* **Peak** – the maximum number of bytes live simultaneously during
  ``__init__`` + ``fit()``, measured by :mod:`tracemalloc`.  This captures
  temporary working arrays (design matrix, normal-equations matrix, etc.)
  as well as the persistent object state.

* **Stored** – the bytes occupied by the numpy arrays that remain on the
  object *after* fitting (``breakpoints``, ``mask``, ``coeff``, ``icoeff``,
  plus the cached design matrix for the new implementation).  Computed
  directly from ``array.nbytes`` so it is not affected by Python-object
  overhead or GC timing.

Usage::

    python bench_bspline_memory.py
    python bench_bspline_memory.py --nrep 5   # average over more runs
"""

import argparse
import gc
import tracemalloc

import numpy as np

from pypeit.bspline.bspline import bspline
from pypeit.bspline.refactor import BSpline, BSpline2D, Knots


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _peak_bytes(func):
    """Return peak heap allocation (bytes) during a single call to ``func``.

    ``func`` should return the fitted spline object so that its lifetime
    spans the measurement window.
    """
    gc.collect()
    tracemalloc.start()
    obj = func()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    del obj
    return peak


def _stored_bytes_old(spl):
    """Sum the nbytes of all numpy arrays stored on an old ``bspline`` object."""
    total = 0
    for attr in ('breakpoints', 'mask', 'coeff', 'icoeff'):
        arr = getattr(spl, attr, None)
        if isinstance(arr, np.ndarray):
            total += arr.nbytes
    return total


def _stored_bytes_new(spl):
    """Sum the nbytes of all numpy arrays stored on a new BSpline/BSpline2D."""
    total = 0
    for attr in ('breakpoints', 'bkpt_gpm', 'coeff', 'icoeff'):
        arr = getattr(spl, attr, None)
        if isinstance(arr, np.ndarray):
            total += arr.nbytes
    # Include the cached design matrix if present
    cached = getattr(spl, '_cached_design', None)
    if cached is not None:
        for item in cached:
            if isinstance(item, np.ndarray):
                total += item.nbytes
    return total


def _fmt(nbytes):
    """Format a byte count as a human-readable string."""
    if nbytes >= 1024 ** 2:
        return f'{nbytes / 1024**2:.2f} MB'
    if nbytes >= 1024:
        return f'{nbytes / 1024:.1f} KB'
    return f'{nbytes} B'


def _header(title):
    width = 88
    print()
    print('=' * width)
    print(f'  {title}')
    print('=' * width)
    print(f'  {"case":<34}  {"old peak":>9}  {"new peak":>9}  '
          f'{"old stored":>10}  {"new stored":>10}')
    print('-' * width)


def _row(label, peak_old, peak_new, stored_old, stored_new):
    print(f'  {label:<34}  {_fmt(peak_old):>9}  {_fmt(peak_new):>9}  '
          f'{_fmt(stored_old):>10}  {_fmt(stored_new):>10}')


# ---------------------------------------------------------------------------
# 1D benchmarks
# ---------------------------------------------------------------------------

def bench_1d_memory(configs, nrep):
    _header('1D fit  [npoly=1]')
    rng = np.random.default_rng(0)

    for N, nbkpts in configs:
        x = np.sort(rng.uniform(0, 10, N))
        y = np.sin(x) + 0.1 * rng.standard_normal(N)
        invvar = np.ones(N)

        def make_old():
            spl = bspline(x=x, nord=4, npoly=1, nbkpts=nbkpts)
            spl.fit(x, y, invvar)
            return spl

        def make_new():
            spl = BSpline(x=x, knots=Knots(count=nbkpts), nord=4)
            spl.fit(x, y, invvar)
            return spl

        peaks_old = [_peak_bytes(make_old) for _ in range(nrep)]
        peaks_new = [_peak_bytes(make_new) for _ in range(nrep)]

        stored_old = _stored_bytes_old(make_old())
        stored_new = _stored_bytes_new(make_new())

        _row(
            f'N={N:>6d}  nknots={nbkpts:>4d}',
            np.median(peaks_old), np.median(peaks_new),
            stored_old, stored_new,
        )

    print()


# ---------------------------------------------------------------------------
# 2D benchmarks
# ---------------------------------------------------------------------------

def bench_2d_memory(configs, nrep, funcname='legendre'):
    _header(f'2D fit  [funcname={funcname}]')
    rng = np.random.default_rng(1)

    for N, nbkpts, npoly in configs:
        x = np.sort(rng.uniform(0, 10, N))
        x2 = rng.uniform(0, 1, N)
        y = np.sin(x) * (1 + 0.3 * x2)
        invvar = np.ones(N)

        def make_old(npoly=npoly):
            spl = bspline(x=x, nord=4, npoly=npoly, nbkpts=nbkpts,
                          funcname=funcname)
            spl.xmin = 0.0
            spl.xmax = 1.0
            spl.fit(x, y, invvar, x2=x2)
            return spl

        def make_new(npoly=npoly):
            spl = BSpline2D(x=x, knots=Knots(count=nbkpts), nord=4)
            spl.fit(x, y, invvar, basis=funcname, npoly=npoly, basis_x=x2,
                    xmin=0.0, xmax=1.0)
            return spl

        peaks_old = [_peak_bytes(make_old) for _ in range(nrep)]
        peaks_new = [_peak_bytes(make_new) for _ in range(nrep)]

        stored_old = _stored_bytes_old(make_old())
        stored_new = _stored_bytes_new(make_new())

        _row(
            f'N={N:>6d}  nknots={nbkpts:>4d}  npoly={npoly}',
            np.median(peaks_old), np.median(peaks_new),
            stored_old, stored_new,
        )

    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nrep', type=int, default=3,
                        help='Repetitions per measurement (default 3)')
    args = parser.parse_args()
    nrep = args.nrep

    print(f'\nMemory benchmark: old bspline vs. new BSpline  (nrep={nrep})')
    print('Peak  = max bytes live simultaneously during __init__ + fit().')
    print('Stored = bytes in numpy arrays persisted on the object after fit().')

    bench_1d_memory(
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

    bench_2d_memory(
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


if __name__ == '__main__':
    main()
