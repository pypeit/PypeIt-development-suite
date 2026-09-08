# `bspline` Test Coverage Assessment

## What is tested

**`test_breakpoints`** — thorough. Exercises `get_breakpoints` and `_fill_bkpt` across all
five breakpoint-specification modes (`bkpt`, `bkspace`, `nbkpts`, `everyn`, `fullbkpt`),
including edge cases (single breakpoint, spacing too large, irregular grids, `bkspread`).

**`test_profile_spec` / `test_profile_spat` / `test_profile_twod`** — integration
regression tests against golden `.npz` files. They drive `fitting.bspline_profile`, which
internally calls `action`, `workit`, `value`, and occasionally `maskpoints`. The 2D test
additionally exercises the `x2` code path in `action`.

**`test_io`** — exercises FITS serialization and a `value()` roundtrip equality check.

---

## What is not tested

### Class methods

| Method | Coverage gap |
|---|---|
| `bsplvn` | Never tested directly. The de Boor recursion and its `order='F'` memory layout are untested in isolation. |
| `action` | Never tested directly. The design matrix structure, the `lower`/`upper` boundary arrays, and the 2D outer-product construction are only verified implicitly by the integration tests. |
| `fit` | Never tested with a controlled synthetic dataset. `fit` is only exercised through `iterfit`; it has a distinct inline normal equations assembly that is separate from `workit` + `solution_arrays`, and neither its correctness nor its equivalence to `workit` is verified. |
| `workit` | Only exercised via `bspline_profile`. Never called with a controlled synthetic dataset where the correct answer is known. |
| `value` | Only checked for roundtrip equality after FITS I/O. Never verified for mathematical correctness, out-of-range masking, gap behavior from masked breakpoints, or the `x2` evaluation path. |
| `maskpoints` | Never tested. The return-code logic (`-1` vs `-2`) and the breakpoint neighborhood masking are completely uncovered. |
| `copy` | Never tested. |
| `reinit_coeff` | Never tested. |

### Utility functions (`util.py`)

| Function | Coverage gap |
|---|---|
| `intrv` | Never tested. Boundary behavior (x at or beyond the first/last breakpoint) is uncovered. |
| `solution_arrays` | Never tested. The `bi`/`bo` index mapping that writes per-span `(bw, bw)` blocks into the banded `alpha` matrix is the most complex IDL-heritage logic in the code and has no direct test. |
| `cholesky_band` | Never tested. Neither the successful factorization path nor the ill-conditioned failure path (returning bad indices) is verified against a matrix with a known answer. |
| `cholesky_solve` | Never tested against a system with a known solution. |
| `bspline_model` | Never tested. The `coeff.flatten('F')` coefficient indexing convention, which couples the model evaluation to the column-major ordering of the solver output, is untested. |

### Module-level

| Item | Coverage gap |
|---|---|
| `uniq` | Has a docstring example but no `pytest` test. |

---

## Critical gaps for refactoring verification

The integration tests will catch gross failures but are insufficient for a refactor, for
these reasons:

1. **No controlled synthetic fit.** No test constructs a function with known coefficients,
   fits it, and verifies that the recovered coefficients and residuals are numerically
   correct. Without this, a refactored solver could introduce a systematic offset or
   scaling error and pass all existing tests.

2. **`fit` and `workit` are never cross-validated.** Both methods solve the same linear
   system — `fit` via an inline assembly loop and `workit` via `solution_arrays` — but no
   test verifies that they produce identical results on the same input. A refactor that
   aligns both paths with standard NumPy conventions needs this equivalence test as a
   baseline.

3. **The `bi`/`bo` index mapping is untested.** This is the most opaque piece of
   IDL-heritage logic, translating dense per-span block contributions into banded storage.
   A refactor that changes this mapping would silently produce wrong coefficients unless
   `solution_arrays` is tested against `cholesky_band`/`cholesky_solve` with a known
   matrix.

4. **The `order='F'` + `flatten('F')` coupling is untested.** `solution_arrays` writes
   coefficients in polynomial-index-fastest order, and `bspline_model` reads them back
   with `coeff.flatten('F')`. A refactor that changes one without the other will silently
   produce zero or garbage output, and no existing test would catch it.

5. **`value` masking is untested.** Any refactor that changes how `value` handles
   out-of-range data or masked-breakpoint gaps would pass all existing tests.
