# PypeIt `bspline` Class Summary

## Architecture Overview

The class fits a weighted least-squares B-spline to 1D data, optionally extended to
quasi-2D by multiplying the B-spline basis by a second-variable polynomial (Legendre,
Chebyshev, or power). The core linear algebra problem is:

```
minimize  Σᵢ wᵢ (yᵢ − A[i,:] · c)²
```

where `A` is a sparse **design matrix** (each row has at most `bw = nord * npoly`
non-zero entries), `w = invvar`, and `c` is the coefficient vector. The normal equations
`(AᵀWA) c = AᵀWy` are assembled and solved using a custom banded Cholesky solver.

The class inherits from `DataContainer` solely for FITS serialization; all computation
is in the methods described below.

---

## Datamodel State

| Attribute | Shape | Meaning |
|---|---|---|
| `breakpoints` | `(nbkpt,)` | Full padded knot vector (includes `nord-1` phantom knots at each end) |
| `mask` | `(nbkpt,)` | Which breakpoints are currently active |
| `coeff` | `(nc,)` or `(npoly, nc)` | Solution coefficients; `nc = nbkpt - nord` |
| `icoeff` | same | First row of Cholesky factor `a[0, :]`, i.e., diagonal of `G` |
| `nord` | `int` | B-spline order (4 = cubic) |
| `npoly` | `int` | Polynomial order in the second variable |
| `xmin`, `xmax` | `float` | Normalization range for the second variable `x2` |
| `funcname` | `str` | Basis function type for `x2`: `'legendre'`, `'chebyshev'`, `'poly'`, `'poly1'` |

---

## `get_breakpoints` / `_fill_bkpt` — Knot Vector Setup

Constructs the full B-spline knot vector. The interior breakpoints are determined by one
of four strategies (`fullbkpt`, `bkpt`, `bkspace`, `nbkpts`, or `everyn`). The result is
then **padded** by `_fill_bkpt`:

```
fullbkpt = [t₀ − (k)·Δ, …, t₀ − Δ, t₀, t₁, …, tₘ, tₘ + Δ, …, tₘ + (k)·Δ]
            ←────── nord-1 phantom knots ──────→              ←── nord-1 ──→
```

where `Δ = (t₁ - t₀) * bkspread`. This padding ensures B-spline basis functions are
non-zero at the edges.

---

## `bsplvn` — B-Spline Basis Evaluation (de Boor Recursion)

Given sorted input positions `x` (length `N`) and their breakpoint interval indices
`ileft` (from `intrv`), this implements the **de Boor triangular recursion** to evaluate
the `nord` non-zero B-spline basis functions at each point.

The recurrence relation is:

```
B[i,1](x) = 1        (piecewise constant)

B[i,k](x) =  (x − tᵢ)/(tᵢ₊ₖ₋₁ − tᵢ) · B[i,k−1](x)
           + (tᵢ₊ₖ − x)/(tᵢ₊ₖ − tᵢ₊₁) · B[i+1,k−1](x)
```

The implementation accumulates left- and right-distance arrays:

```
deltap[:, j] = breakpoints[ileft + j + 1] - x     (right distances)
deltam[:, j] = x - breakpoints[ileft - j]          (left distances)
```

The algorithm iterates from order 1 up to `nord`, building each row of `vnikx` in-place.
The result is `vnikx` of shape `(N, nord)`, where `vnikx[i, :]` are the `nord` non-zero
B-spline basis values at `x[i]`.

**Column-major artifact**: `vnikx`, `deltap`, `deltam` are allocated `order='F'`.

---

## `action` — Design Matrix Construction

Returns the **banded design matrix** `A` of shape `(N, bw)` where `bw = nord * npoly`,
along with integer arrays `lower` and `upper`.

**1D case** (`x2=None`): `A = bsplvn(x, indx)` — simply the B-spline basis matrix.

**2D case** (`x2` provided): The design matrix is the row-wise outer product of the
B-spline basis with the second-variable polynomial basis:

```
A[i, ii*npoly + jj] = B[i, ii] * P[i, jj]
```

where `B` is the B-spline basis `(N, nord)` and `P` is the polynomial basis `(N, npoly)`
evaluated on the normalized `x2`:

```
x2norm = 2 * (x2 - xmin) / (xmax - xmin) - 1        ∈ [-1, +1]
P = flegendre(x2norm, npoly)   or fchebyshev(...)
```

The `lower[k]` and `upper[k]` arrays give the first and last row indices of `A` where
data falls in B-spline span `k` (i.e., between `breakpoints[k+nord-1]` and
`breakpoints[k+nord]`). These are used to split the normal equations assembly into
per-span blocks.

**Column-major artifact**: `action` array is allocated `order='F'`.

---

## `fit` / `workit` + `solution_arrays` — Normal Equations Assembly

The weighted normal equations `(AᵀWA) c = AᵀWy` are assembled in **banded storage**.

Because `A` has bandwidth `bw`, the matrix `AᵀWA` is symmetric positive definite with
bandwidth `bw`. It is stored as `alpha` with shape `(bw, nfull + bw)` where
`nfull = nn * npoly`.

The assembly loop iterates over spans `k = 0..nn-nord`:

```python
a2 = A * sqrt(W)[:, None]                  # whitened design matrix
# Per span:
alpha_block += a2[lower[k]:upper[k]+1, :].T  @  a2[lower[k]:upper[k]+1, :]   # shape (bw, bw)
beta_block  += (y * sqrt(w))[lower[k]:upper[k]+1]  @  a2[lower[k]:upper[k]+1, :]
```

The `(bw, bw)` blocks are written into `alpha` via pre-computed index arrays `bi` (flat
indices into the dense block, upper triangle only) and `bo` (flat indices into the banded
storage, offset by `itop * bw`). This accumulation is the most complex and
IDL-heritage-heavy part of the code.

`fit` performs this assembly inline with a Python loop; `workit` delegates to
`solution_arrays` which does the same thing more cleanly.

---

## `cholesky_band` — Banded Cholesky Decomposition

Given the banded symmetric positive definite matrix `l` of shape `(bw, n + bw)`, this
computes the Cholesky factor `G` such that `l = GᵀG` (stored in-place in band format).

For column `j = 0..n-1`:

```
G[0, j]     = sqrt(G[0, j])              ← diagonal element
G[1:, j]   /= G[0, j]                   ← sub-diagonal column
G[i, j+k]  -= G[i, j] * G[k, j]        ← rank-1 update of trailing submatrix
                                            for i, k = 1..bw-1
```

Returns `(-1, G)` on success, or `(bad_indices, l)` if any diagonal element is ≤
`mininf` or non-finite — triggering breakpoint masking.

---

## `cholesky_solve` — Banded Cholesky Solve

Solves `GᵀG x = b` in two passes:

**Forward substitution** (solves `G y = b`):

```
b[j] /= a[0, j]
b[j + 1:j + bw] -= b[j] * a[1:bw, j]
```

**Backward substitution** (solves `Gᵀ x = y`):

```
b[j] = (b[j] - Σₖ a[k, j] * b[j+k]) / a[0, j]
```

Result is stored in `b` in-place; always returns `(-1, b)`.

---

## `bspline_model` — Model Evaluation

After solving, evaluates the fitted model at all `x` values using the pre-computed
`action` matrix:

```
ŷ[lower[k]:upper[k]+1] = A[lower[k]:upper[k]+1, :] · c_flat[k*npoly : k*npoly + bw]
```

where `c_flat = coeff.flatten('F')` — the **Fortran-order** (column-major) flattening of
the coefficient array. This ordering convention is the crux of the IDL heritage: when
`coeff` has shape `(npoly, nn)`, flattening column-major gives the coefficient vector in
the order the banded assembly expects.

---

## `value` — Public Evaluation Entry Point

Sorts `x`, calls `action` to build the design matrix, calls `bspline_model` to compute
`ŷ`, then sets `mask=False` for any data outside the span of good breakpoints (including
gaps created by masked-out breakpoints).

---

## `maskpoints` — Breakpoint Masking

When `cholesky_band` returns a set of bad column indices, this method maps those column
indices back to breakpoint indices and masks a neighborhood of size `nord` centered on
each bad breakpoint. Returns `-1` if masking succeeded (trigger re-fit), `-2` if too few
breakpoints remain.

---

## Column-Major Artifacts (Summary)

The IDL-origin column-major conventions appear in four places:

1. **`vnikx`, `deltap`, `deltam`** in `bsplvn` — allocated `order='F'`
2. **`action`** matrix — allocated `order='F'`
3. **`coeff.flatten('F')`** in `bspline_model` — coefficient lookup uses Fortran order
4. **`sol.reshape(npoly, nn, order='F')`** in `fit`/`workit` — coefficient storage uses Fortran order

Items 3 and 4 are semantically coupled: the solution vector from `cholesky_solve` is
ordered as `[c₀₀, c₁₀, ..., c_{npoly-1,0}, c₀₁, ..., c_{npoly-1,nn-1}]` (polynomial
index fastest) to match how the `bi`/`bo` index mapping accumulates the normal equations.
The `order='F'` allocations in 1 and 2 likely only affect cache performance, not
correctness, but have never been verified to be removable.
