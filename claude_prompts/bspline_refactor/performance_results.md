# BSpline Refactor — Performance Results

**Date:** 2026-06-11
**Commit:** `1f02367ce` (branch `bspline_refactor`)
**Platform:** macOS, Intel Core i5-1038NG7 @ 2.00 GHz (x86\_64)
**Python:** 3.13.13 · **NumPy:** 2.4.4

Both benchmarks compare the legacy `pypeit.bspline.bspline.bspline` class against the
refactored `pypeit.bspline.refactor.BSpline` / `BSpline2D` classes on identical data.

---

## Execution-time benchmark (`bench_bspline.py`, nrep=20)

Timing is the **median wall-time** over 20 repetitions.

### 1D fit (npoly=1)

| N | nknots | old (ms) | new (ms) | speedup |
|------:|-------:|---------:|---------:|--------:|
|    500 |  20 | 1.44 | 0.73 | 1.95× |
|    500 | 100 | 5.16 | 1.94 | 2.66× |
|  2,000 |  20 | 1.96 | 0.82 | 2.39× |
|  2,000 | 100 | 5.94 | 2.07 | 2.88× |
| 10,000 |  20 | 5.24 | 1.86 | 2.81× |
| 10,000 | 100 | 8.78 | 3.30 | 2.66× |

### 2D fit (funcname=legendre)

| N | nknots | npoly | old (ms) | new (ms) | speedup |
|------:|-------:|------:|---------:|---------:|--------:|
|    500 |  20 | 3 |  3.49 | 1.27 | 2.75× |
|    500 | 100 | 3 | 13.06 | 2.87 | 4.55× |
|  2,000 |  20 | 3 |  4.32 | 1.52 | 2.85× |
|  2,000 | 100 | 3 | 13.94 | 3.31 | 4.21× |
|  2,000 |  20 | 6 |  8.84 | 2.31 | 3.82× |
| 10,000 |  20 | 3 |  8.06 | 3.40 | 2.37× |
| 10,000 | 100 | 3 | 17.22 | 5.11 | 3.37× |

### Sigma-clipping loop (10 iterations, same x, varying invvar)

The new class caches the design matrix after the first `fit()` call and reuses it on
subsequent calls with the same `x`.  The old class recomputes it on every call.

| N | nknots | old (ms) | new (ms) | speedup |
|------:|-------:|---------:|---------:|--------:|
|  2,000 |  20 | 19.15 |  5.47 | 3.50× |
|  2,000 | 100 | 55.21 | 18.38 | 3.00× |
| 10,000 |  20 | 52.12 | 10.12 | 5.15× |

---

## Memory benchmark (`bench_bspline_memory.py`, nrep=5)

**Peak** = maximum bytes live simultaneously during `__init__` + `fit()`, measured by
`tracemalloc`.

**Stored** = bytes in numpy arrays persisted on the object *after* fitting:
`breakpoints`/`bkpt_gpm`, `coeff`, `icoeff`, and (new only) the cached design matrix
`_cached_design`.

### 1D fit (npoly=1)

| N | nknots | old peak | new peak | old stored | new stored |
|------:|-------:|---------:|---------:|-----------:|-----------:|
|    500 |  20 |  78.7 KB |  70.5 KB |    586 B |  16.5 KB |
|    500 | 100 |  85.5 KB |  73.0 KB |  2.5 KB |  19.7 KB |
|  2,000 |  20 | 301.2 KB | 269.5 KB |    586 B |  63.4 KB |
|  2,000 | 100 | 305.0 KB | 272.0 KB |  2.5 KB |  66.6 KB |
| 10,000 |  20 |  1.45 MB |  1.30 MB |    586 B | 313.4 KB |
| 10,000 | 100 |  1.46 MB |  1.30 MB |  2.5 KB | 316.6 KB |

### 2D fit (funcname=legendre)

| N | nknots | npoly | old peak | new peak | old stored | new stored |
|------:|-------:|------:|---------:|---------:|-----------:|-----------:|
|    500 |  20 | 3 | 270.2 KB | 178.9 KB |   1.3 KB |  48.4 KB |
|    500 | 100 | 3 | 505.0 KB | 182.0 KB |   5.7 KB |  54.1 KB |
|  2,000 |  20 | 3 | 692.1 KB | 509.1 KB |   1.3 KB | 189.1 KB |
|  2,000 | 100 | 3 | 926.9 KB | 514.7 KB |   5.7 KB | 194.8 KB |
|  2,000 |  20 | 6 |  1.58 MB | 939.7 KB |   2.3 KB | 377.6 KB |
| 10,000 |  20 | 3 |  3.09 MB |  2.24 MB |   1.3 KB | 939.1 KB |
| 10,000 | 100 | 3 |  3.15 MB |  2.27 MB |   5.7 KB | 944.8 KB |

---

## Summary

**Speed:** The new classes are consistently faster across all configurations.

- Single `fit()` call: **2–5× faster**, with larger improvements when there are many
  knots or a large polynomial basis (the refactored solver uses a more efficient linear
  algebra path).
- Sigma-clipping loops (repeated `fit()` on the same `x`): **3–5× faster**, with the
  largest gains at high N, because the design matrix is built once and cached.

**Memory:** The new classes use **10–40% less peak memory** during fitting, owing to
fewer intermediate temporaries.  However, the stored footprint is larger: the legacy
class discards the design matrix after fitting and keeps only a few small arrays (0.6–
5.7 KB), whereas the new class retains `_cached_design` to accelerate future calls
(16 KB – 944 KB depending on N and nknots).  This is an intentional time-vs-space
tradeoff; the cache can be released by calling `reset_coeff()` or by setting
`_cached_design = None` when memory pressure is a concern.
