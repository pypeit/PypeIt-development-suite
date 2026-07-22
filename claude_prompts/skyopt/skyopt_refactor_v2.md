# Eliminate `.flat` usage in `skysub.py::skyoptimal` (side-by-side comparison pass)

## Context

`pypeit/core/spatialprofile.py::fit_profile` was previously refactored (commit
`4be258b7d`) so that it no longer needs `.flatten()`/`.flat[idx]`/`.reshape()`
scaffolding around calls to `pypeit.core.fitting.iterative_bspline_fit`.
That routine has always been, and remains, strictly 1D-only (`x`, `y`, `ivar`,
`gpm`, `basis` must be flat arrays) — the refactor did not touch
`iterative_bspline_fit` itself. Instead, `fit_profile` moved the
flatten-on-demand logic *inside* the function: it derives a 2D boolean mask
once (`totmask`), immediately extracts 1D "valid pixel" arrays via boolean
fancy-indexing (`image[totmask]`, `ivar[totmask]`, ...), does all its
per-pixel math on those 1D arrays, and only reconstructs a 2D array once at
the very end via `profile_model[totmask] = full_bsp * pb_x`. That refactor
was itself originally developed as a parallel implementation
(`spatialprofile_refactor.py`) alongside the original, and only merged in
(replacing the original, deleting the parallel file) once validated. The
untracked `notes` file in the repo root explicitly lists
`skysub.py::skyoptimal` as the next planned target for the same `.flat`-
elimination treatment.

`skyoptimal` (`pypeit/core/skysub.py:259-379`) is the analogous function for
joint sky/object B-spline fitting. Its docstring requires all of `piximg`,
`data`, `ivar`, `spatial_img` to be pre-flattened 1D arrays and `oprof` to be
pre-reshaped to `(nflat, nobj)` — pushing all the flatten/scatter work onto
its single caller, `local_skysub_extract` (`skysub.py:532`), specifically the
block at `skysub.py:934-1008`. That block builds
`isub = np.where(localmask.flatten())`, uses `.flat[isub]` on `piximg`,
`sciimg`, `modelivar*skymask`, `spatial_img`, `base_var`, `count_scale`,
`sciivar`, reshapes `obj_profiles` (3D) to 2D before indexing it, and
scatters every result back with `.flat[isub] = ...` /
`outmask.flat[isub[igood1]] = ...`. This is the last place in `skysub.py`
using `.flat`/`.flatten()`/`.reshape()` — confirmed by grep, nothing else in
the file uses them.

**This pass follows the same parallel-development approach used for
`fit_profile`, but keeps both versions in the same file rather than a
separate module.** A new function, `skyoptimal_refactor`, is added directly
below the existing `skyoptimal` in `skysub.py`, implementing the
native-shape + explicit-mask idiom. `local_skysub_extract` is **not**
changed in this pass — it keeps calling the existing `skyoptimal` exactly as
today. This lets both versions be exercised side by side in a single test
for direct, in-process comparison, with no need to freeze `.npy` reference
arrays (unlike the `spatialprofile.py` regression test, which needed frozen
references because the old implementation no longer existed once refactored
in place). Switching the caller over, and eventually retiring the old
`skyoptimal` (mirroring the notes file's `bspline_refactor_deprec ->
spatprof_refactor` rename-once-validated pattern), is left as a deliberate
follow-up once this comparison confirms equivalence.

## Changes

### 1. `pypeit/core/skysub.py::skyoptimal` (lines 259-379) — unchanged

Left exactly as-is. It remains the reference implementation for the
comparison test and the one actually used by `local_skysub_extract` in this
pass.

### 2. New `pypeit/core/skysub.py::skyoptimal_refactor` — added directly below `skyoptimal` (i.e. immediately after line 379, before `optimal_bkpts`)

New function, new signature:

```python
def skyoptimal_refactor(piximg, data, ivar, oprof, localmask, sigrej=3.0, npoly=1,
                         spatial_img=None, fullbkpt=None):
```

- `piximg`, `data`, `ivar`, `spatial_img`: native-shape (e.g. `(nspec, nspat)`)
  arrays instead of pre-flattened 1D.
- `oprof`: `(nspec, nspat, nobj)` instead of pre-reshaped `(nflat, nobj)`.
- `localmask`: new required 2D boolean parameter (the sub-region selection
  currently applied by the caller via `isub`).

Body, mirroring `fit_profile`'s `image[totmask]` idiom:
- Extract the 1D "valid pixel" arrays via boolean fancy-indexing:
  ```python
  piximg_x = piximg[localmask]
  data_x = data[localmask]
  ivar_x = ivar[localmask]
  oprof_x = oprof[localmask, :]        # (nspec,nspat,nobj)[mask] -> (n, nobj)
  spatial_x = None if spatial_img is None else spatial_img[localmask]
  ```
- Replace `nc = oprof.shape[0]; nobj = int(oprof.size/nc)` / `nx = data.size`
  with the equivalent computed from the `_x` arrays (`nx = data_x.size`,
  `nobj = oprof_x.shape[1]`), validating shape consistency against
  `piximg.shape`/`oprof.shape[:2]` up front instead of the old `nc != nx`
  check.
- Everything from the current `sortpix = piximg.argsort(...)` through the
  two `iterative_bspline_fit` calls (lines 312-367 in the original) is
  otherwise unchanged, just operating on the `_x`-suffixed arrays instead of
  the caller-flattened ones.
- Final model evaluation currently runs `sset.value(piximg, ...)` over the
  *entire* input (i.e., the whole `isub`/local-mask region, not just the
  `good` ivar>0 subset) — keep that behavior by evaluating over `piximg_x`
  (not the full frame), then scatter into full-shape outputs:
  ```python
  sky_bmodel = np.zeros_like(piximg, dtype=float)
  obj_bmodel = np.zeros_like(piximg, dtype=float)
  gpm = np.zeros(piximg.shape, dtype=bool)
  ...
  sky_bmodel_x, _ = sset.value(piximg_x, basis=profile_basis[:, nobj:], coeff=sset.coeff[:, nobj:])
  obj_bmodel_x, _ = sset.value(piximg_x, basis=profile_basis[:, :nobj], coeff=sset.coeff[:, :nobj])
  sky_bmodel[localmask] = sky_bmodel_x
  obj_bmodel[localmask] = obj_bmodel_x
  gpm_x = np.zeros(nx, dtype=bool)
  gpm_x[good] = gpm_good
  gpm[localmask] = gpm_x
  ```
- The two early-return branches (`ngood == 0`, all pixels rejected after
  first pass) return the same zero/False full-shape arrays.
- Docstring describes `piximg`/`data`/`ivar`/`spatial_img` as native-shape
  arrays, `oprof` as `(nspec, nspat, nobj)`, documents the new `localmask`
  parameter, and describes the returns as full-shape arrays (zero/False
  outside `localmask`).

`local_skysub_extract` is not touched — it keeps its current call to
`skyoptimal` with `.flat[isub]`/`.reshape()` unchanged.

### 3. Comparison test

Add a new test to `pypeit/tests/test_skysub.py` (no existing test covers
`skyoptimal`) that calls **both** functions on the same synthetic data within
a single test — no frozen `.npy` files needed since the old implementation
still exists in the file:

- Build small synthetic 2D `piximg`, `data`, `ivar`, `oprof` (3D), and
  `localmask` arrays (deterministic; fixed seed if randomness is involved),
  exercising both the `npoly==1` and `npoly>1`/`spatial_img` branches, and a
  case where some pixels have `ivar<=0` inside `localmask` (to exercise the
  internal "good" sub-selection).
- Call the current `skyoptimal`, replicating the caller's existing
  convention inline in the test:
  ```python
  isub, = np.where(localmask.flatten())
  sky_old, obj_old, gpm_old = skyoptimal(
      piximg.flat[isub], data.flat[isub], ivar.flat[isub],
      oprof.reshape(nspec * nspat, nobj)[isub, :],
      spatial_img=spatial_img.flat[isub] if spatial_img is not None else None,
      fullbkpt=fullbkpt, sigrej=sigrej, npoly=npoly)
  ```
- Call the new function directly:
  ```python
  sky_new, obj_new, gpm_new = skyoptimal_refactor(
      piximg, data, ivar, oprof, localmask,
      spatial_img=spatial_img, fullbkpt=fullbkpt, sigrej=sigrej, npoly=npoly)
  ```
- Assert equivalence at the masked locations:
  ```python
  assert np.allclose(sky_old, sky_new[localmask])
  assert np.allclose(obj_old, obj_new[localmask])
  assert np.array_equal(gpm_old, gpm_new[localmask])
  ```

## Verification

- Run the new comparison test (and the full `pypeit/tests/test_skysub.py`)
  to confirm `skyoptimal_refactor` matches `skyoptimal` bit-for-bit on
  synthetic data.
- Run the full unit test suite (`pytest pypeit/tests/`) to catch any
  incidental breakage (none expected, since `skyoptimal` and
  `local_skysub_extract` are unmodified in this pass).
- `grep -n "\.flat\b\|\.flatten(\|reshape(" pypeit/core/skysub.py` will still
  show hits — the untouched `skyoptimal`/`local_skysub_extract` pairing and,
  intentionally, the comparison test's inline replication of the caller
  convention. Confirm instead that `skyoptimal_refactor`'s own body contains
  none.

## Deferred follow-up (not part of this pass)

Once the comparison test (and, ideally, a real dev-suite reduction) confirms
`skyoptimal_refactor` is equivalent:
1. Switch `local_skysub_extract` to call `skyoptimal_refactor` with native
   arrays + `localmask`, and replace its post-fit `.flat[isub]` /
   `isub[igood1]` bookkeeping (lines 943-1008) with direct `[localmask]`
   boolean indexing (e.g. `isub[igood1]` becomes the boolean combination
   `localmask & skymask` wherever it indexes the full-shape `outmask`).
2. Remove the old `skyoptimal` and rename `skyoptimal_refactor` ->
   `skyoptimal`, matching the notes file's `bspline_refactor_deprec ->
   spatprof_refactor` rename-once-validated convention.
3. Retire the comparison test (or fold its synthetic cases into a plain
   regression test for the single remaining `skyoptimal`).

## Implementation summary

This pass was implemented exactly as planned above, with no deviations from
the design.

- **`pypeit/core/skysub.py`**: `skyoptimal` (lines 259-379) and
  `local_skysub_extract` were left byte-for-byte unchanged — confirmed via
  `git diff --stat` showing a pure insertion (152 lines added, 0 removed/
  modified). `skyoptimal_refactor` was added directly below `skyoptimal`
  (before `optimal_bkpts`), matching the signature, body structure, and
  docstring described above. `local_skysub_extract` still calls the original
  `skyoptimal` at its one call site.
- **`pypeit/tests/test_skysub.py`**: added a `_make_skyoptimal_inputs()`
  helper (deterministic, `np.random.default_rng(1234)`-seeded) that builds
  synthetic `(50, 16)` `piximg`/`data`/`ivar`, a `(50, 16, 2)` `oprof` (two
  Gaussian object profiles), and a deliberately **irregular** `localmask` — a
  central spatial band with ~10% of its pixels randomly dropped, so the mask
  is not expressible as a rectangular sub-image/`np.ix_` slice, which
  exercises genuine boolean fancy-indexing rather than a simpler bounding-box
  crop. A separate ~10% subset of pixels inside the (trimmed) mask has
  `ivar` set to zero, exercising `skyoptimal`'s internal ivar>0 "good"
  sub-selection.
  - `test_skyoptimal_refactor_matches_skyoptimal`: runs both `npoly=1` and
    `npoly=3`/`spatial_img` cases, calling `skyoptimal` via the caller's
    existing flatten convention (`isub`/`.flat[isub]`/`reshape`) built inline
    in the test, and `skyoptimal_refactor` directly on the native-shape
    arrays + `localmask`. Asserts `sky_new[localmask]`/`obj_new[localmask]`/
    `gpm_new[localmask]` match the old outputs (`np.allclose`/
    `np.array_equal`), and that all three outputs are zero/False everywhere
    outside `localmask`.
  - `test_skyoptimal_refactor_all_masked`: degenerate case with `ivar` zeroed
    everywhere, confirming both functions hit the `ngood == 0` early return
    and produce identical all-zero/all-False results.
  - Both new tests pass.
- **Verification results**:
  - `pytest pypeit/tests/test_skysub.py -v`: 5 passed (the 3 pre-existing
    tests plus the 2 new ones).
  - `pytest pypeit/tests/ --ignore=pypeit/tests/test_runpypeit.py`: 623
    passed, 4 failed. All 4 failures are pre-existing and unrelated to this
    change: `test_fluxspec.py::test_flux_calib` fails with
    `AttributeError: module 'pypeit.scripts' has no attribute 'flux_calib'`
    (an import-ordering issue in the scripts package, untouched by this
    work); the 3 `test_pkgdata.py` failures (`test_fetch_github_files`,
    `test_github_contents`, `test_cache_to_pkg`) all 404 on
    `https://raw.githubusercontent.com/pypeit/PypeIt/skyopt/...` — they
    build the fetch URL from the current local git branch name (`skyopt`),
    which doesn't exist on the `pypeit/PypeIt` remote, an environment
    artifact unrelated to `skysub.py`.
  - `pytest pypeit/tests/test_runpypeit.py`: 1 passed (run separately since
    it's a slow, full end-to-end reduction test; ~5m13s).
  - `grep` confirms `skyoptimal_refactor`'s own code body contains no
    `.flat`/`.flatten()`/`reshape()` — the only match within that function's
    span is a docstring sentence describing the *old* convention for
    comparison purposes.

Net effect: `skyoptimal_refactor` is a drop-in, `.flat`-free alternative to
`skyoptimal`, validated equivalent on synthetic data, sitting unused by any
caller until the deferred follow-up switches `local_skysub_extract` over and
retires the original.
