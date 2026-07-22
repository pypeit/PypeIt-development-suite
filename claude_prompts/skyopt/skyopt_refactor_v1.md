# Eliminate `.flat` usage in `skysub.py::skyoptimal`

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
the very end via `profile_model[totmask] = full_bsp * pb_x`. The untracked
`notes` file in the repo root explicitly lists
`skysub.py::skyoptimal` as the next planned target for the same treatment
("eliminate the need for `.flat`").

`skyoptimal` (`pypeit/core/skysub.py:259-379`) is the analogous function for
joint sky/object B-spline fitting. Its docstring requires all of `piximg`,
`data`, `ivar`, `spatial_img` to be pre-flattened 1D arrays and `oprof` to be
pre-reshaped to `(nflat, nobj)` — pushing all the flatten/scatter work onto
its single caller, `local_skysub_extract` (`skysub.py:532`), specifically the
block at `skysub.py:934-1008`. That block builds `isub = np.where(localmask.flatten())`,
uses `.flat[isub]` on `piximg`, `sciimg`, `modelivar*skymask`, `spatial_img`,
`base_var`, `count_scale`, `sciivar`, reshapes `obj_profiles` (3D) to 2D before
indexing it, and scatters every result back with `.flat[isub] = ...` /
`outmask.flat[isub[igood1]] = ...`. This is a direct mirror of the pre-refactor
`fit_profile` pattern and is the last place in `skysub.py` using `.flat` /
`.flatten()` / `.reshape()` — confirmed by grep, nothing else in the file uses
them.

Goal: apply the same idiom used in `fit_profile` — accept native-shape
(2D/3D) arrays plus an explicit 2D boolean mask, do boolean-index extraction
internally, and scatter results back into freshly allocated full-shape
outputs — so neither `skyoptimal` nor its caller need manual `.flat`/`isub`
bookkeeping. The core two-pass B-spline fit math inside `skyoptimal`
(argsort, `ivar>0` "good" subset, `relative` mask, the two
`iterative_bspline_fit` calls) already operates on flat 1D arrays and does
not need to change — only the entry and exit points do. The
"avoid instantiating new BSpline objects as copies" note is a separate,
unrelated concern and is out of scope here.

## Changes

### 1. `pypeit/core/skysub.py::skyoptimal` (lines 259-379)

New signature:

```python
def skyoptimal(piximg, data, ivar, oprof, localmask, sigrej=3.0, npoly=1,
               spatial_img=None, fullbkpt=None):
```

- `piximg`, `data`, `ivar`, `spatial_img`: now native-shape (e.g. `(nspec, nspat)`)
  arrays instead of pre-flattened 1D.
- `oprof`: now `(nspec, nspat, nobj)` instead of pre-reshaped `(nflat, nobj)`.
- `localmask`: new required 2D boolean parameter (the sub-region selection
  previously applied by the caller via `isub`).

Body changes:
- At the top, extract the 1D "valid pixel" arrays via boolean fancy-indexing,
  mirroring `fit_profile`'s `image[totmask]` idiom:
  ```python
  piximg_x = piximg[localmask]
  data_x = data[localmask]
  ivar_x = ivar[localmask]
  oprof_x = oprof[localmask, :]        # (nspec,nspat,nobj)[mask] -> (n, nobj)
  spatial_x = None if spatial_img is None else spatial_img[localmask]
  ```
- Replace `nc = oprof.shape[0]; nobj = int(oprof.size/nc)` /
  `nx = data.size` with the equivalent computed from the `_x` arrays
  (`nx = data_x.size`, `nobj = oprof_x.shape[1]`), and validate shape
  consistency against `piximg.shape`/`oprof.shape[:2]` up front instead of
  the old `nc != nx` check.
- Everything from the current `sortpix = piximg.argsort(...)` through the
  two `iterative_bspline_fit` calls (lines 312-367) is otherwise unchanged,
  just operating on the `_x`-suffixed arrays instead of the caller-flattened
  ones (`piximg` → `piximg_x`, `data` → `data_x`, `ivar` → `ivar_x`,
  `profile_basis` built from `oprof_x`/`spatial_x`).
- Final model evaluation currently runs `sset.value(piximg, ...)` over the
  *entire* input (i.e., the whole `isub`/local-mask region, not just the
  `good` ivar>0 subset) — keep that behavior by evaluating over `piximg_x`
  (not the full frame), then scatter into full-shape outputs:
  ```python
  sky_bmodel = np.zeros_like(piximg, dtype=float)
  obj_bmodel = np.zeros_like(piximg, dtype=float)
  gpm = np.zeros(piximg.shape, dtype=bool)
  ...
  sky_bmodel[localmask], _ = sset.value(piximg_x, basis=profile_basis[:, nobj:], coeff=sset.coeff[:, nobj:])
  obj_bmodel[localmask], _ = sset.value(piximg_x, basis=profile_basis[:, :nobj], coeff=sset.coeff[:, :nobj])
  gpm_x = np.zeros(nx, dtype=bool)
  gpm_x[good] = gpm_good
  gpm[localmask] = gpm_x
  ```
  (Note: `sset.value(...)` returns a tuple; assigning directly into
  `sky_bmodel[localmask]` from the first tuple element requires an
  intermediate 1D variable — write it as
  `sky_bmodel_x, _ = sset.value(...); sky_bmodel[localmask] = sky_bmodel_x`.)
- The two early-return branches (`ngood == 0`, all pixels rejected after
  first pass) return the same zero/False full-shape arrays instead of
  `np.zeros_like(piximg)` on the (formerly flat) input — this already works
  unchanged since `piximg` is now full-shape.
- Update the docstring to describe `piximg`/`data`/`ivar`/`spatial_img` as
  native-shape arrays, `oprof` as `(nspec, nspat, nobj)`, add the `localmask`
  parameter, and describe the returns as full-shape arrays (zero/False
  outside `localmask`) rather than "same shape as piximg" (which is still
  technically true but should clarify the masking behavior).

### 2. `pypeit/core/skysub.py::local_skysub_extract` (lines 934-1008)

Remove `isub = np.where(localmask.flatten())` and
`obj_profiles_flat = obj_profiles.reshape(nspec*nspat, objwork)`. Call
`skyoptimal` with native-shape arrays and `localmask` directly:

```python
sky_bmodel, obj_bmodel, outmask_opt = skyoptimal(
        piximg, sciimg, modelivar * skymask, obj_profiles, localmask,
        spatial_img=spatial_img, fullbkpt=fullbkpt, sigrej=sigrej_eff, npoly=npoly)
```

Replace every `.flat[isub]` in the post-fit block (lines 960-999) with
`[localmask]` boolean indexing on the same arrays — this is a direct
substitution since `arr.flat[isub]` and `arr[localmask]` produce identical
1D arrays in identical order when `isub = np.where(localmask.flatten())`.
For the compound `isub[igood1]` indexing (selecting, within the
localmask-restricted domain, the subset where `igood1` — itself
`skymask.flat[isub]` — is true), replace with the equivalent 2D boolean
combination `localmask & skymask` wherever it's used to index the full-shape
`outmask` array directly (this is exactly the same pixel set — set algebra:
"positions in `localmask` where `skymask` is also true" — but expressed
without an intermediate integer index array). Keep the intermediate chi2/
rejection arithmetic (`chi2`, `igood`, `chi2_srt`, `sigind`, etc.) working on
1D arrays extracted via `[localmask]`, exactly mirroring today's logic with
`.flat[isub]` → `[localmask]` as the only substitution, to minimize risk of
introducing an ordering/indexing bug.

Also update the failure-fallback branch (currently lines 1006-1008):
```python
skyimage[localmask] = global_sky[localmask]
```

### 3. Regression test

Add a new test to `pypeit/tests/test_skysub.py` (no existing test covers
`skyoptimal` or `local_skysub_extract`) following the pattern of
`pypeit/tests/test_spatialprofile.py::test_refactor_regression`:

- Build small synthetic 2D `piximg`, `data`, `ivar`, `oprof` (3D), and
  `localmask` arrays (deterministic, fixed-seed if any randomness is
  involved) exercising both the `npoly==1` and `npoly>1`/`spatial_img`
  branches.
- Before making the code change, run the **current** `skyoptimal` with the
  caller's existing flatten convention (`piximg.flat[isub]`, etc.) on this
  synthetic data and freeze the `sky_bmodel`, `obj_bmodel`, `gpm` outputs as
  `.npy` reference files (same convention as
  `pypeit/tests/files/spatprof_*_ref.npy`, e.g. via a small
  `freeze_skyoptimal.py` helper script analogous to
  `pypeit/tests/files/freeze_spatprof.py`).
- After refactoring, add `test_skyoptimal_refactor_regression` that calls the
  new `skyoptimal(piximg, data, ivar, oprof, localmask, ...)` signature and
  asserts the full-shape `sky_bmodel[localmask]`/`obj_bmodel[localmask]`/
  `gpm[localmask]` match the frozen 1D reference arrays to `atol=1e-10`.

## Verification

- Run the new `test_skyoptimal_refactor_regression` (and full
  `pypeit/tests/test_skysub.py`) to confirm bit-for-bit equivalence with
  pre-refactor behavior on synthetic data.
- Run the full unit test suite (`pytest pypeit/tests/`) to catch any other
  incidental breakage.
- Run (or ask the user to run) a real `run_pypeit` reduction for a
  spectrograph that exercises local sky subtraction (e.g. an existing
  dev-suite MultiSlit setup) and confirm the reduced 2D/1D spectra are
  unchanged from a pre-refactor run — this is the only way to validate
  `local_skysub_extract`'s full behavior end-to-end since it requires real
  calibrated images that unit tests deliberately avoid.
- `grep -n "\.flat\b\|\.flatten(\|reshape(" pypeit/core/skysub.py` should
  return no results after the change.
