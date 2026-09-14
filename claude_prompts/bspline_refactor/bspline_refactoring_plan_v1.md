# BSpline Refactoring Plan — v1

**Date**: 2026-06-08  
**Author**: Kyle Westfall (kwestfal@ucsc.edu) / Claude Sonnet 4.6  
**Branch**: `bspline_refactor`

---

## Context

The `bspline` class in `pypeit/bspline/bspline.py` is a direct port of IDL's PYDL bspline
implementation. Its problems: (1) column-major (`order='F'`) memory coercion throughout to
preserve IDL layout, which is foreign to NumPy and harms readability; (2) matrix operations
use `np.dot` and manual loops instead of NumPy/SciPy idioms; (3) the 1D and quasi-2D
fitting paths are entangled in a single class; (4) `fit` and `workit` perform the same
algorithm with slight differences (one assembles the design matrix itself, the other accepts
it pre-built), creating maintainability risk.

---

## Design Goals

1. **C-order memory throughout** — remove all `order='F'` allocations, `flatten('F')`, and
   `reshape(..., order='F')`.
2. **NumPy/SciPy idioms** — `@` for matrix multiplication; `scipy.linalg` banded-SPD
   routines to replace the hand-rolled Cholesky; `scipy.sparse` noted below.
3. **Modular class hierarchy** — `BSpline` handles 1D; `BSpline2D` subclass handles
   second-variable polynomial extension.
4. **One method per purpose** — `fit` absorbs `workit`; normal-equations assembly and
   model evaluation move into single-purpose private methods.

---

## Target File

New code goes in **`pypeit/bspline/refactor.py`**. The existing `bspline.py` and `util.py`
are left completely untouched. This allows direct numerical comparison between old and new
implementations during development and review.

A separate follow-on step (not in scope here) will add `DataContainer`-based wrapper classes
and wire them into the rest of the codebase.

---

## Proposed Class Hierarchy

```
BSpline                         ← 1D weighted least-squares B-spline (no base class)
    BSpline2D(BSpline)          ← extends BSpline with x2/npoly/funcname
```

`BSpline` has **no base class** — it is self-contained and carries no dependency on
`DataContainer` or any other PypeIt infrastructure. FITS serialisation will be added later
via a thin wrapper class (not in this file). `BSpline2D` inherits only from `BSpline`.

---

## Memory Order Changes

Four coupled artefacts to remove, in order of dependency:

| Current code | Proposed replacement |
|---|---|
| `vnikx = np.zeros(..., order='F')` | `np.zeros(...)` (default C) |
| `deltap = vnikx.copy()` / `deltam = vnikx.copy()` | same, no-op (inherits from vnikx) |
| `action = np.zeros(..., order='F')` | `np.zeros(...)` |
| `sol.reshape(npoly, nn, order='F')` → `coeff[jj, k]` | `sol.reshape(nn, npoly)` → `coeff[k, jj]` |
| `coeff.flatten('F')[i*npoly + spot]` | eliminated; see `_evaluate_model` below |

The new `coeff` shape for 2D fits is `(nc, npoly)` (knot index first), not `(npoly, nc)`.
Indexing in `value()` that currently uses `coeff[..., coeffbk]` becomes `coeff[coeffbk, :]`.

The correctness argument: `sol.reshape(nn, npoly)[k, jj] == sol[k*npoly + jj]` (C-order
reshape fills row-major, so `coeff[k, jj] = sol[k*npoly + jj]`).

### Policy on `flatten` / `ravel`

- **`flatten()` is prohibited** — it always copies memory and obscures intent.
- **`ravel()` is discouraged** — use only if no multi-dimensional alternative exists, and
  document why.
- **`reshape()`** is acceptable because it returns a *view* (no memory copy) whenever the
  array is contiguous. It is the preferred tool when a shape change is genuinely needed (e.g.
  for passing to a banded-matrix routine that requires a 2D input).

---

## `BSpline` — 1D Base Class

### Attributes (no DataContainer; plain instance attributes)

| Attribute | Shape | Notes |
|---|---|---|
| `breakpoints` | `(nbkpt,)` | full padded knot vector |
| `mask` | `(nbkpt,)` | active breakpoints |
| `coeff` | `(nc,)` | solution coefficients |
| `icoeff` | `(nc,)` | diagonal of Cholesky factor |
| `nord` | int | B-spline order |

### Private Methods

Each isolates one algorithmic step with a clear docstring.

| Method | Replaces | Purpose |
|---|---|---|
| `_build_breakpoints(x, ...)` (static) | `get_breakpoints` + `_fill_bkpt` | construct + pad knot vector |
| `_find_spans(x, breakpoints, nord)` (static) | `intrv` in util.py | find B-spline interval index for each x |
| `_bspline_basis(self, x, ileft)` | `bsplvn` | de Boor recursion → `(N, nord)` basis matrix, C-order |
| `_build_design_matrix(self, x)` | `action` (1D path) | returns `(A, lower, upper)` |
| `_assemble_normal_equations(self, A, y, w, lower, upper)` | `solution_arrays` in util.py | banded `alpha (bw, nfull)` and `beta (nfull,)` via `A.T @ (W*A)` per span using `@` |
| `_solve_banded(self, alpha, beta)` | `cholesky_band` + `cholesky_solve` | returns `(sol, bad_cols)`; see scipy note below |
| `_update_coefficients(self, sol, goodbk)` | inline in `fit`/`workit` | stores `coeff`, `icoeff` from solution |
| `_evaluate_model(self, A, lower, upper)` | `bspline_model` in util.py | evaluates model per span; see implementation note below |
| `_mask_breakpoints(self, bad_cols)` | `maskpoints` | maps bad Cholesky cols → breakpoint indices; returns -1 (retry) or -2 (too few) |

### `_evaluate_model` — No-flatten Implementation

`A` has shape `(N, bw)` where `bw = nord` (1D) or `bw = nord * npoly` (2D). `coeff` has
shape `(nc,)` (1D) or `(nc, npoly)` (2D).

**1D** — direct matrix–vector multiply, no reshape:
```python
for k in range(nn - nord + 1):
    if lower[k] <= upper[k]:
        yfit[lower[k]:upper[k]+1] = A[lower[k]:upper[k]+1, :] @ coeff[k:k+nord]
```

**2D** — treat `A` as a 3D view per span to avoid flattening `coeff`:
```python
for k in range(nn - nord + 1):
    if lower[k] <= upper[k]:
        sl = slice(lower[k], upper[k] + 1)
        # reshape is a view (no copy); coeff slice is (nord, npoly)
        yfit[sl] = np.einsum('nij,ij->n',
                             A[sl, :].reshape(-1, self.nord, self.npoly),
                             coeff[k:k+nord, :])
```

The `einsum` contracts the `(nord, npoly)` coefficient block against the local action tensor
without ever creating a 1D copy of the coefficients.

### Public Methods

| Method | Notes |
|---|---|
| `fit(self, xdata, ydata, invvar)` | Unified: builds design matrix, assembles normal equations, solves, stores coefficients, returns `(err, yfit)`. Absorbs `workit` entirely. See caller-update note. |
| `value(self, x)` | Calls `_build_design_matrix` then `_evaluate_model` |
| `reinit_coeff(self)` | Reset `coeff` / `icoeff` to zeros |
| `copy(self)` | Deep copy of all attributes |

---

## `BSpline2D(BSpline)` — Quasi-2D Subclass

### Additional Datamodel

| Attribute | Shape | Notes |
|---|---|---|
| `coeff` | `(nc, npoly)` | overrides shape vs 1D case |
| `icoeff` | `(nc, npoly)` | same |
| `npoly` | int | polynomial order in x2 |
| `xmin`, `xmax` | float | normalisation range for x2 |
| `funcname` | str | `'legendre'`, `'chebyshev'`, `'poly'`, `'poly1'` |

### Additional/Overridden Private Methods

| Method | Purpose |
|---|---|
| `_normalize_x2(self, x2)` | maps x2 → `[-1, +1]` via `xmin`/`xmax` |
| `_poly_basis(self, x2norm)` | dispatches on `funcname`; returns `(N, npoly)` polynomial basis matrix |
| `_build_design_matrix(self, x, x2)` | overrides base; builds outer product — see below |

### Additional/Overridden Public Methods

`x2` is a **required positional argument** in `BSpline2D` — not optional. A `BSpline2D`
object always operates on two-variable data.

| Method | Signature |
|---|---|
| `fit` | `fit(self, xdata, ydata, invvar, x2)` |
| `value` | `value(self, x, x2)` |

### `_build_design_matrix` — Vectorised 2D Outer Product

The `ii`/`jj` nested loop is replaced by broadcasting with a single `reshape` (a view):

```python
# B shape (N, nord), P shape (N, npoly)
# Outer product: A_3d[n, ii, jj] = B[n, ii] * P[n, jj]
A_3d = B[:, :, np.newaxis] * P[:, np.newaxis, :]   # (N, nord, npoly), no copy
# Reshape to (N, bw) so that A[:, ii*npoly + jj] = B[:, ii] * P[:, jj]
A = A_3d.reshape(N, self.nord * self.npoly)          # view, no copy
```

This is fully vectorised with C-order arrays and no `flatten`.

---

## `scipy.linalg` Replacement for Banded Cholesky

`_solve_banded` uses `scipy.linalg.cholesky_banded` + `scipy.linalg.cho_solve_banded`
(LAPACK `dpbtrf` / `dpbtrs`):

```python
from scipy.linalg import cholesky_banded, cho_solve_banded

def _solve_banded(self, alpha, beta):
    nfull = alpha.shape[1] - alpha.shape[0]  # actual unknowns
    try:
        chol = cholesky_banded(alpha[:, :nfull], lower=False)
        sol = cho_solve_banded((chol, False), beta[:nfull])
        return sol, np.array([-1])          # success sentinel
    except np.linalg.LinAlgError:
        bad = np.where(alpha[0, :nfull] <= mininf)[0]
        return None, bad
```

The current `alpha` banded storage (`alpha[0, :]` = diagonal, `alpha[k, :]` = k-th
superdiagonal) matches scipy's upper-triangular banded format. **Verify column-alignment
convention** during implementation (off-by-one shifts between scipy's `ab[i, j]` = `a[j-i,
j]` and the current `bo`-indexed accumulation).

**On `scipy.sparse`**: The design matrix `A` is sparse (at most `bw` nonzeros per row), but
for typical problem sizes (N ~ 10⁴, nfull ~ 400) the banded LAPACK path is faster and more
memory-efficient than `scipy.sparse` general sparse solvers. `scipy.sparse` would become
advantageous only for very large `npoly`/`nord`; defer to a future optimisation.

---

## `fit` / `workit` Consolidation

`workit` exists because `core/fitting.py` pre-computes the action matrix once and calls
`workit` in a sigma-clipping loop. In the refactored design, `fit` handles everything
internally. To preserve the optimisation, `fit` will **cache the design matrix** on the
instance:

```python
def fit(self, xdata, ydata, invvar):
    if self._cached_design is None or xdata.shape != self._cached_x_shape:
        self._cached_design = self._build_design_matrix(xdata)
    A, lower, upper = self._cached_design
    alpha, beta = self._assemble_normal_equations(A, ydata, invvar, lower, upper)
    sol, bad_cols = self._solve_banded(alpha, beta)
    ...
```

`_cached_design` is NOT part of the datamodel (not serialised). The cache is invalidated
whenever `x` changes shape (for sigma-clipping the shape stays fixed, so the cache hits
every iteration after the first).

`workit` is **removed**. `action` (the public design-matrix getter) is also removed since
its only non-test caller is `core/fitting.py`, which will be updated.

---

## `util.py` / `bspline.py` Changes

**None.** Both files are left intact. All refactored logic lives in `refactor.py`.

The utility functions `intrv`, `uniq`, `solution_arrays`, `cholesky_band`, `cholesky_solve`,
and `bspline_model` are **re-implemented** inside `refactor.py` as private methods of
`BSpline` (or as module-level helpers within that file). They are not imported from `util.py`.

---

## Caller Updates

**None in this PR.** Because the new classes live in a new file, no existing callers need to
change. Integration with `core/fitting.py`, `flatfield.py`, etc. is a follow-on step once
the refactored classes are validated.

---

## Verification

Existing tests confirm the old code is unchanged:
```bash
pytest pypeit/tests/test_bspline.py pypeit/tests/test_bspline_refactor.py -v
```

New per-method tests live in **`pypeit/tests/test_bspline_refactor_new_classes.py`**.
One test function (or small group) per method, covering both `BSpline` and `BSpline2D`:

| Method under test | What to verify |
|---|---|
| `BSpline._build_breakpoints` | knot count, padding (nord-1 phantoms at each end), all input strategies (bkspace, nbkpts, everyn, bkpt, fullbkpt) |
| `BSpline._find_spans` | correct interval bracketing; clamping at edges; sorted-input assumption |
| `BSpline._bspline_basis` | partition of unity; non-negativity; exact linear case; C-order output |
| `BSpline._build_design_matrix` | shape `(N, nord)`; lower/upper consistency with `_find_spans`; full data coverage |
| `BSpline._assemble_normal_equations` | `alpha` / `beta` match `A.T @ W @ A` / `A.T @ W @ y` for a small known system |
| `BSpline._solve_banded` | known tridiagonal SPD system; success path and failure path (bad breakpoint triggers correct bad-column indices) |
| `BSpline._evaluate_model` | matches direct `A @ c` for a known coefficient vector |
| `BSpline._mask_breakpoints` | neighborhood masking; `-2` return when too few breakpoints remain |
| `BSpline.fit` | cubic polynomial recovery ≤ 1e-10; smooth function residuals ≤ 2e-3 |
| `BSpline.value` | matches `fit` output at training points; out-of-range mask |
| `BSpline2D._normalize_x2` | maps to `[-1, +1]`; edge values |
| `BSpline2D._poly_basis` | correct values for each `funcname`; shape `(N, npoly)` |
| `BSpline2D._build_design_matrix` | outer-product structure `A[:, ii*npoly+jj] = B[:, ii] * P[:, jj]`; shape `(N, bw)` |
| `BSpline2D.fit` | exact recovery of `h0(x) + h1(x)*P1(x2)` (Legendre, Chebyshev, poly); smooth 2D residuals ≤ 2e-3 |
| `BSpline2D.value` | matches `fit` output at training (x, x2) pairs |

**Cross-check tests** (enabled by the parallel-file approach):

For identical input data, `BSpline.fit` and old `bspline.fit` must return identical fitted
values, and `BSpline.coeff` must equal `old.coeff.T` (transposed due to the shape change
from `(npoly, nc)` to `(nc, npoly)`).

---

## Implementation Notes (added 2026-06-08)

### scipy banded Cholesky format: `lower=True`, not `lower=False`

The plan stated that the `alpha` banded storage (`alpha[0, :]` = diagonal, `alpha[k, :]`
= k-th off-diagonal) "matches scipy's upper-triangular banded format."  This was incorrect.

scipy's **upper** banded format (`lower=False`) stores the diagonal in the **last** row
(`alpha[bw-1, :]`) and the highest superdiagonal in row 0.  Our `alpha` stores the
diagonal in **row 0**, which matches scipy's **lower** banded format (`lower=True`).
`_solve_banded` was originally coded with `lower=False`; this caused scipy to interpret
our diagonal (positive) as the off-diagonal and our off-diagonal (variable sign) as the
diagonal, making every Cholesky attempt fail and silently triggering breakpoint masking
in an infinite loop.  The fix is:

```python
chol = cholesky_banded(alpha[:, :nfull], lower=True)
sol  = cho_solve_banded((chol, True), beta[:nfull])
```

After `cholesky_banded(lower=True)`, `chol[0, :]` is the diagonal of the lower Cholesky
factor `L` such that `A = L Lᵀ`, which is what we store in `icoeff`.

### Test file

Per-method tests live in
`pypeit/tests/test_bspline_refactor_new_classes.py` (66 tests, flat function structure,
no test-class wrappers).  The `_assemble_normal_equations` test verifies correctness by
constructing the full `(N, nn)` sparse design matrix and comparing yfit from the banded
solve to `numpy.linalg.lstsq`.  The `_solve_banded` tests use a hand-built tridiagonal
alpha in scipy lower banded format: `alpha[0, :n] = diagonal`,
`alpha[1, :n-1] = off-diagonal` (NOT `alpha[1, 1:n]`, which would be off by one).
