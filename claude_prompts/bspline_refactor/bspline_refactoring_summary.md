# B-Spline Refactoring Summary

**Date:** 2026-07-15

| Branch | Commit |
|---|---|
| `develop` | `a41980f8ef6aa5ac05812bc17d2b29effc149a44` |
| `bspline_refactor_deprec` | `87e9fa5e2ef31dcc4e3d433b0d92b6302ffb2d40` |

---

## Old architecture (on `develop`)

| Component | Location | Notes |
|---|---|---|
| Core class | `pypeit/core/bspline/bspline.py` — `class bspline(DataContainer)` | ~740 lines; Fortran memory order |
| Low-level helpers | `pypeit/core/bspline/util.py` — `cholesky_band`, `cholesky_solve`, `intrv`, `bspline_model` | ~250 lines |
| Package init | `pypeit/core/bspline/__init__.py` — `from ... import bspline` | |
| Simple iterative fit | `fitting.iterfit()` | Wraps `bspline.bspline`; no polynomial basis |
| Profile-basis fit | `fitting.bspline_profile()` | Main fitting loop; takes `invvar` and `profile_basis` as positional args |

## New architecture (on `bspline_refactor_deprec`)

| Component | Location | Notes |
|---|---|---|
| Core classes | `pypeit/core/bspline.py` — `Knots`, `BSpline`, `BSpline2D` | Single flat module, ~1867 lines |
| FITS containers | `pypeit/containers/bspline.py` — `BSplineContainer`, `BSpline2DContainer` | Inherits DataContainer + BSpline/BSpline2D |
| Package init | `pypeit/containers/__init__.py` | Empty file for new subpackage |
| Unified fitting | `fitting.iterative_bspline_fit()` | Replaces both `iterfit` and `bspline_profile` |

---

## Key design changes in the core classes (`core/bspline.py`)

- **C-order memory layout throughout** — the old code used Fortran column-major order
  (`order='F'`, `flatten('F')`, `reshape(..., order='F')`); all gone.
- **`coeff` shape for 2D fits is `(nc, npoly)`** — the old `bspline` stored it as
  `(npoly, nc)`; callers that slice it (e.g. `coeff[0, :]`) must be updated to
  `coeff[:, 0]` (already done in `spatialprofile.py`).
- **`fit` and `workit` merged** — the old class had separate `fit()`, `action()`, and
  `workit()` methods; `BSpline.fit()` performs the full solve in one call.
- **Design matrix caching** — `BSpline` caches `_cached_design` between `fit()` calls,
  avoiding redundant computation in sigma-clipping loops.
- **scipy banded Cholesky** — `_cholesky_banded()` replaces the hand-rolled
  `cholesky_band` / `cholesky_solve` from `util.py`; it returns the LAPACK `info` flag
  instead of raising to allow degenerate breakpoints to be identified and masked.
- **Knot construction isolated in `Knots`** — formerly scattered across `bspline.__init__`.
  `Knots` maps the old keyword names to new ones:

  | Old (`bspline.bspline`) | New (`Knots`) |
  |---|---|
  | `bkspace=s` | `spacing=s` |
  | `nbkpts=n` | `count=n` |
  | `everyn=n` | `stride=n` |
  | `bkpt=t` | `interior=t` |
  | `bkspread=f` | `spread=f` |
  | `fullbkpt=t` | `full=t` |

- **`bkpt_gpm` replaces `mask`** — the breakpoint active-mask attribute was renamed.
  All callers updated.
- **`breakpoints` is a read-only property** → delegates to `self.knots.breakpoints`;
  FITS containers store it under `bkpt_full` to avoid a naming conflict with the
  property descriptor.
- **Quasi-2D separated into `BSpline2D`** — old `bspline` handled both 1D and 2D via the
  `npoly` parameter; now they are distinct classes. `BSpline2D.fit()` accepts either a
  string basis name (`'legendre'`, `'chebyshev'`, `'poly'`, `'poly1'`) plus `basis_x`,
  or a pre-built array.

---

## Changes to `fitting.py`

`iterfit` (~120 lines) and `bspline_profile` (~200 lines) are deleted entirely and
replaced by a single `iterative_bspline_fit()` (~200 lines). Key differences:

| Aspect | Old (`iterfit` / `bspline_profile`) | New (`iterative_bspline_fit`) |
|---|---|---|
| 1D vs 2D dispatch | Separate functions | Single function; `basis=None` → 1D, else 2D |
| `profile_basis` / `x2` | Positional required args | Optional `basis` and `basis_x` kwargs |
| `invvar` | Positional in `bspline_profile` | `ivar` keyword |
| `ingpm` | `bspline_profile` keyword | `gpm` keyword |
| Knot control | `kwargs_bspline={'bkspace': s}` | `kwargs_knots={'spacing': s}` |
| `quiet` flag | Suppresses all logging | Removed; logging always happens |
| Return on failure | Returns `sset` (old class instance) | Returns `None` for catastrophic failure |
| `bspline_qa` | Referenced `sset.mask` | Updated to `sset.bkpt_gpm` |

---

## Callers updated

All production call sites were migrated to `iterative_bspline_fit`:

- **`flatfield.py`** — spectral 1D fits, 2D spatial flat fits; `spat_bsplines` and
  `illumflat_spat_bsplines` datamodel types updated from `bspline.bspline` to
  `BSplineContainer`.
- **`core/skysub.py`** — 1D and quasi-2D (Legendre) sky subtraction.
- **`core/flux_calib.py`** — 1D flux calibration zeropoint fit.
- **`core/scattlight.py`** — 1D scattered-light fits (formerly `iterfit` with `everyn`).
- **`core/spatialprofile.py`** — 1D and 2D profile fits; coefficient indexing updated
  from `[0, :]` / `[1, :]` (old `(npoly, nc)` layout) to `[:, 0]` / `[:, 1]` (new
  `(nc, npoly)` layout).
- **`core/findobj_skymask.py`** — minor formatting cleanup only (stale import removed).

The stale imports of `from pypeit.core import bspline` were removed from `flatfield.py`,
`skysub.py`, and `flux_calib.py` after those modules were confirmed to use only the new API.

---

## Documentation

`doc/bspline.rst` is a new file on this branch (not present on `develop`). It covers:

- The algorithm: 1D B-spline WLS formulation, quasi-2D polynomial extension,
  iterative sigma-clipping rejection loop.
- Usage examples for `BSpline`, `BSpline2D`, and `iterative_bspline_fit`.
- A migration section mapping the old API to the new one (class instantiation,
  fitting, evaluation, and the profile-fitting function).

**Note:** As of the commit above, `doc/bspline.rst` contains stale module paths
(`pypeit.bspline.refactor.*`, `bspline_profile_refactor`) that reflect an earlier
draft of the module layout. The correct paths are `pypeit.core.bspline.{BSpline,
BSpline2D, Knots}`, `pypeit.containers.bspline.{BSplineContainer, BSpline2DContainer}`,
and `pypeit.core.fitting.iterative_bspline_fit`. These will be corrected in the
upcoming docstring update pass.

---

## Tests

| File | Status | Notes |
|---|---|---|
| `tests/test_bspline.py` | Replaced | ~1915 lines of new per-method unit tests for `BSpline`, `BSpline2D`, `Knots`, `_cholesky_banded`, and `iterative_bspline_fit`; the original 5-test file for `bspline.bspline` is gone |
| `tests/test_bspline_containers.py` | New | ~525 lines; covers `BSplineContainer` and `BSpline2DContainer` construction, round-trip FITS I/O, overwrite, empty serialization, copy, `from_bspline`/`from_bspline2d`, version/datamodel |
| `tests/test_bspline_refactor_new_classes.py` | New (added earlier on branch) | Per-method tests for `Knots` strategies, `BSpline`/`BSpline2D` internals, `reset_coeff`, `_mask_breakpoints`, fit, evaluate |
| `tests/test_flatfield.py` | Updated | Uses `BSplineContainer` instead of `bspline.bspline` |
