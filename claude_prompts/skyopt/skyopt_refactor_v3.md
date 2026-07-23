# Eliminate `.flat` usage in `skysub.py::skyoptimal` (follow-up: switch-over and retirement)

## Context

`skyopt_refactor_v2.md` added `skyoptimal_refactor` — a native-shape,
explicit-`localmask` version of `skyoptimal` — directly below the original
`skyoptimal` in `pypeit/core/skysub.py`, without touching `skyoptimal` itself
or its one caller, `local_skysub_extract`. That let both versions be
exercised side by side in a single in-process comparison test
(`test_skyoptimal_refactor_matches_skyoptimal`), with no need for frozen
`.npy` reference files.

Before requesting this follow-up, the user validated the v2 state
thoroughly: the full `pypeit/pypeit/tests` suite, the full
`PypeIt-development-suite/unit_tests` suite, and both the `reduce` and
`afterburn` dev-suite tests for the `shane_kast_blue/600_4310_d55` dataset —
all passing cleanly. Separately, the user also fixed an unrelated
pre-existing `develop`-level regression in
`pypeit/tests/test_fluxspec.py::test_flux_calib` (commit `33f729d94`), which
had nothing to do with this refactor (see `project_scripts_flux_calib_attributeerror`
memory note) but was surfaced while running the full suite.

With `skyoptimal_refactor` validated equivalent, this pass executes the
"Deferred follow-up" from `skyopt_refactor_v2.md`: wire `local_skysub_extract`
over to the new implementation, retire the old `skyoptimal`, and fold the
now-redundant comparison test into a plain regression test for the single
remaining `skyoptimal`.

## Changes

### 1. `pypeit/core/skysub.py::local_skysub_extract` — switched over

Removed the `isub = np.where(localmask.flatten())` construction and the
`obj_profiles_flat = obj_profiles.reshape(nspec * nspat, objwork)` line.
The `skyoptimal` call now passes native-shape arrays plus `localmask`
directly:

```python
sky_bmodel, obj_bmodel, outmask_opt = skyoptimal(
        piximg, sciimg, modelivar * skymask, obj_profiles, localmask,
        spatial_img=spatial_img, fullbkpt=fullbkpt, sigrej=sigrej_eff, npoly=npoly)
```

Every `.flat[isub]` in the post-fit block was replaced with `[localmask]`
boolean indexing on the same (now full-shape) arrays returned by
`skyoptimal`:

```python
skyimage[localmask] = sky_bmodel[localmask]
objimage[localmask] = obj_bmodel[localmask]
img_minsky[localmask] = sciimg[localmask] - sky_bmodel[localmask]
igood1 = localmask & skymask
outmask[igood1] = outmask_opt[igood1]
```

The compound `isub[igood1]` indexing (selecting, within the isub domain,
the subset where `igood1 = skymask.flat[isub]` was true) is exactly the set
of pixels where both `localmask` and `skymask` are true, so it collapses to
the single boolean combination `igood1 = localmask & skymask`, reused
everywhere `isub[igood1]` previously indexed the full-shape `outmask`,
`chi2`, and `sciivar` arrays. `chi2` itself became a full-shape array,
computed only within `localmask` (`chi2[localmask] = (img_minsky[localmask]
- obj_bmodel[localmask])**2 * modelivar[localmask]`), and `igood = igood1 &
(chi2 <= chi2_sigrej**2)` replaces the old `(skymask.flat[isub]) & (chi2 <=
chi2_sigrej**2)` — equivalent since `igood1` already *is*
`skymask.flat[isub]` re-expressed as a 2D mask. `base_var`/`count_scale`
became `[localmask]`-indexed, and `procimg.variance_model`'s `counts=`
argument became `(sky_bmodel + obj_bmodel)[localmask]` (previously the raw
flat return values, now sliced from the full-shape return). The
global-sky fallback branch (`skyimage.flat[isub] = global_sky.flat[isub]`)
became `skyimage[localmask] = global_sky[localmask]`.

### 2. `pypeit/core/skysub.py` — old `skyoptimal` removed, `skyoptimal_refactor` renamed

Deleted the original `skyoptimal` (the pre-flattened-1D-input version)
entirely. Renamed `skyoptimal_refactor` -> `skyoptimal`, updating its
signature accordingly:

```python
def skyoptimal(piximg, data, ivar, oprof, localmask, sigrej=3.0, npoly=1,
               spatial_img=None, fullbkpt=None):
```

Its docstring's "This is a refactored version of `skyoptimal`..." framing
paragraph (which only made sense while both versions coexisted) was
removed, since there's no longer an old version to contrast against — it
now reads as a normal, single docstring for the only `skyoptimal`. The two
internal `log.warning(...)` messages that said "in skyoptimal_refactor"
were corrected back to "in skyoptimal".

### 3. `pypeit/tests/test_skysub.py` — comparison test folded into a plain regression test

Replaced `test_skyoptimal_refactor_matches_skyoptimal` /
`test_skyoptimal_refactor_all_masked` (which called both functions and
asserted equivalence) with `test_skyoptimal` / `test_skyoptimal_all_masked`,
which call the single remaining `skyoptimal` directly and check its
properties intrinsically rather than by comparison to a second
implementation:

- output shapes match the input `piximg` shape;
- `sky_bmodel`/`obj_bmodel`/`gpm` are zero/False everywhere outside
  `localmask`;
- `gpm` can only be `True` where `localmask & (ivar > 0)`;
- the recovered `sky_bmodel` tracks the known noise-free synthetic sky
  trend (`4.0 + 0.03 * piximg`) to `atol=0.5` within `localmask` — a
  correctness check that wasn't needed before (equivalence to the old
  implementation was the check), but is meaningful now that there's only
  one implementation to test.
- `test_skyoptimal_all_masked` keeps the degenerate `ivar<=0`-everywhere
  case, now checked directly against the single `skyoptimal`.

The `_make_skyoptimal_inputs()` helper (irregular `localmask`, scattered
zero-`ivar` pixels, both `npoly` branches) is unchanged.

## Verification

- `pytest pypeit/tests/test_skysub.py -v`: 5 passed.
- `pytest pypeit/tests/ --ignore=pypeit/tests/test_runpypeit.py`: 627
  passed, 0 failed (the previously-noted `test_fluxspec.py`/`test_pkgdata.py`
  failures are gone now that the user's separate `test_fluxspec.py` fix,
  commit `33f729d94`, is in place).
- `pytest pypeit/tests/test_runpypeit.py`: 1 passed (~5m13s).
- `grep -n "\.flat\b\|\.flatten(\|reshape(" pypeit/core/skysub.py`: no
  matches anywhere in the file — the elimination is now complete, not just
  in the new function but in `local_skysub_extract` too.
- `grep -rn "skyoptimal_refactor" --include="*.py" .`: no matches — the
  temporary name is gone everywhere.
- `grep -rn "skyoptimal(" --include="*.py" pypeit/`: only the definition,
  the one call site in `local_skysub_extract`, and the two test calls —
  confirms no other caller was missed.

This closes out the `skyopt` refactor: `pypeit/core/skysub.py` no longer
uses `.flat`/`.flatten()`/`.reshape()` anywhere, `skyoptimal` now has the
same native-shape + explicit-mask contract as
`spatialprofile.py::fit_profile`, and there is exactly one implementation
of `skyoptimal` again.
