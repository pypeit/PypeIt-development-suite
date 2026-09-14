# Integrate `construct_basename` into `setup_datacube.py`

## Context

This branch stack is tracked in `kcwi_wcs.md` (same directory):
`develop` <- `kcwi_dec_2024` <- `kcwi_dec_2024_weights` <- `fix-cube-manual-coords-kbw`
(`PypeIt` repo), with the identically-named `PypeIt-development-suite` branches mapped
1:1. `pypeit/scripts/setup_datacube.py` (`pypeit_setup_datacube`) was added on this branch
before the separate `basename` branch existed, and was merged into this stack without ever
being updated to use it.

The `basename` branch (see `basename_refactoring_plan.md`, same directory) extracted
`PypeItMetaData.construct_basename` into a standalone
`outputfiles.construct_basename(filename, target, camera, mjd, allowed_extensions)` (plus
a companion `outputfiles.strip_raw_extension(filename, allowed_extensions)` helper),
specifically fixing basenames for raw files whose extension isn't a literal
`.fits`/`.fits.gz` substring, and routed several existing call sites through it
(`coadd2d.py`, `pypeit_steps.py`, `rawimage.py`, `skysub_regions.py`, `ql.py`).
`setup_datacube.py` was not part of that sweep, since it lives on a different branch.

**Trigger**: running `pypeit_setup_datacube keck_kcwi_D.pypeit "SDSSJ2222 2745"` against
the real OMEARA KCWI test reduction (`/Users/westfall/Work/pypeit/kcwi_wcs/omeara/rdx/blue`)
finds **zero** matching `Science/spec2d` frames, even though four exist.

**Root cause, confirmed**: `find_reduced_spec2d()` (`setup_datacube.py:237`) does
`raw_stem = Path(raw_filename).stem`. The OMEARA raw filenames are compound-extension
(e.g. `KB.20230916.28203.70.fits.gz`); `Path.stem` only strips the *last* suffix, leaving
`raw_stem = "KB.20230916.28203.70.fits"`. The glob pattern `spec2d_{raw_stem}-*.fits` then
becomes `spec2d_KB.20230916.28203.70.fits-*.fits`, which can never match the real product
`spec2d_KB.20230916.28203.70-SDSSJ22222745_KCWI_20230916T075003.725.fits`, whose root was
correctly stripped of the full `.fits.gz` extension by `construct_basename`'s underlying
logic when the file was originally written. This is exactly the class of bug the
`basename` branch fixed everywhere else — this script just wasn't in scope for that sweep.

There's a second, cosmetic instance of the same bug at `setup_datacube.py:288`
(`existing_spec2d_files`), where the same naive `.stem` is used only to name a "missing"
raw frame in a log warning.

## Investigation: can `construct_basename` reconstruct the exact spec2d name?

Rather than just patching the extension-stripping in place (which `strip_raw_extension`
alone would do, preserving the existing glob + per-candidate-header-read design), the
better integration is to ask whether the *full* `outputfiles.construct_basename` call can
reproduce the real spec2d basename **exactly**, using only what's already present in the
`.pypeit` file's data table (`filename`, `target`, `mjd`) plus the spectrograph's
`camera`/`allowed_extensions`. If so, `find_reduced_spec2d` can look for that literal file
directly, instead of a glob followed by opening each candidate's FITS header to verify the
target.

The initial working assumption was that this *wouldn't* be reliable: the `.pypeit` file's
`mjd` column is textually truncated to 6 decimal places (e.g. `60203.326432`), which looked
too coarse to reproduce the sub-second timestamp token embedded in real spec2d filenames
(`..._20230916T075003.725.fits`). **This concern does not hold up, for two independent
reasons**:

1. `astropy`'s `ascii.fixed_width` writer (used by `PypeItMetaData.write_pypeit`/
   `inputfiles.InputFile` to serialize the data block) does not truncate float precision on
   write — verified directly by writing a synthetic 17-significant-figure float through it
   and confirming it comes back in full. So the 6-decimal-place `mjd` values seen in real
   `.pypeit` files are the genuine, exact stored values, not an artifact of ASCII
   serialization.
2. Calling `outputfiles.construct_basename(...)` directly with the OMEARA `.pypeit` file's
   own `filename`/`target`/`mjd` table values, and `'KCWI'` for `camera` (matching
   `KeckKCWISpectrograph.camera` in `spectrographs/keck_kcwi.py:833`), reproduces **all
   four** real spec2d filenames present in `Science/` bit-for-bit:

   | raw filename | `construct_basename(...)` output | matches real file? |
   |---|---|---|
   | `KB.20230916.28203.70.fits.gz` | `spec2d_KB.20230916.28203.70-SDSSJ22222745_KCWI_20230916T075003.725.fits` | yes |
   | `KB.20230916.29633.74.fits.gz` | `spec2d_KB.20230916.29633.74-SDSSJ22222745_KCWI_20230916T081353.731.fits` | yes |
   | `KB.20230916.53503.39.fits.gz` | `spec2d_KB.20230916.53503.39-gd50_KCWI_20230916T145143.373.fits` | yes |
   | `KB.20230916.53597.23.fits.gz` | `spec2d_KB.20230916.53597.23-gd50_KCWI_20230916T145317.203.fits` | yes |

This works because `PypeItMetaData.construct_basename` (used once, during the original
`run_pypeit` reduction that produced these files) and `pypeit_setup_datacube` (reading the
same value back out of the `.pypeit` file afterward) both ultimately read the identical
`mjd` table value — nothing recomputes it with different precision in between, and the
ASCII round-trip doesn't lose any digits.

**Residual risk, and why a fallback is still worth keeping.** This guarantee rests on the
`mjd`/ASCII round-trip being lossless, which was verified for this data and the current
`astropy` version, but is not proven for every spectrograph's `mjd` metadata precision,
`.pypeit` files written by an older PypeIt/`astropy` version, or a row with `mjd = None`
(`metadata.py` explicitly anticipates "possibly None mjds if there were corrupt header
cards"). For that reason, the integration below tries the exact match first and falls back
to the current (bug-fixed) prefix-glob + header-verify approach if it misses, rather than
dropping that logic entirely.

## Design

### 1. `pypeit/scripts/setup_datacube.py`

Add `from pypeit import outputfiles` to the imports. No new spectrograph-loading import is
needed — reuse the existing `inputfiles.PypeItFile.get_spectrograph()` convenience method
already present in the codebase.

- **`find_reduced_spec2d(science_dir, row, spectrograph)`**: change signature to take the
  matched data-table `row` (with `filename`/`target`/`mjd` columns) and the spectrograph
  instance, instead of separate `raw_filename`/`target` strings.
  - *Primary path*: if `row['mjd']` is not `None`, build
    `expected = outputfiles.construct_basename(row['filename'], row['target'],
    spectrograph.camera, row['mjd'], spectrograph.allowed_extensions)` and check whether
    `Path(science_dir) / f'spec2d_{expected}.fits'` exists; if so, return it directly — no
    header read needed, since the row's own `target` was used to build the name, so a hit
    is target-correct by construction.
  - *Fallback path* (`row['mjd'] is None`, or the exact file isn't found): keep today's
    logic, corrected to use `outputfiles.strip_raw_extension(row['filename'],
    spectrograph.allowed_extensions)` for the raw-stem prefix (this is what makes the
    fallback itself correct for compound extensions, same fix as the primary path), glob
    `spec2d_{raw_stem}-*.fits`, and verify each candidate's header target via the existing
    `spec2d_target()`/`target_matches()`, including the existing multiple-candidates
    warning.
  - Update the docstring to describe the new signature and two-tier behavior.

- **`existing_spec2d_files(pypeit_file, target, science_dir, spectrograph)`**: add the
  `spectrograph` parameter; pass `group[0]` (the row, not just its `filename`) and
  `spectrograph` through to `find_reduced_spec2d(...)`; use
  `outputfiles.strip_raw_extension(row['filename'], spectrograph.allowed_extensions)` for
  the `missing` stem (fixes the same cosmetic bug at line 288). Update the docstring.

- **`SetupDataCube.main()`**: after `spectrograph = pypeit_file.config['rdx']['spectrograph']`
  (still needed as the plain string written into the `.coadd3d` file's `[rdx] spectrograph`
  line), add `spec = pypeit_file.get_spectrograph()` and pass `spec` as the new last
  argument to `existing_spec2d_files(...)`.

### 2. `pypeit/tests/test_setup_datacube.py`

- Add a test exercising `find_reduced_spec2d`'s primary (exact-match) path directly: a
  synthetic row with a `.fits.gz` raw filename and a spec2d file on disk named exactly as
  `outputfiles.construct_basename` would produce it, asserting it's found even with no (or
  a mismatched) `TARGET` header — proving the exact-match path doesn't depend on
  `spec2d_target()`.
- Add a test exercising the fallback path: a spec2d file whose name does *not* exactly
  match `construct_basename`'s output (simulating an `mjd`-precision miss or `mjd=None`),
  but whose raw-stem prefix does, with a matching header `TARGET` — asserting it's still
  found via the glob+header route, with a `.fits.gz`-style compound raw extension handled
  correctly there too. This is the regression test for the original bug, covering both code
  paths.
- Extend `_write_pypeit_file`/`test_setup_datacube_write_and_append` (or add a sibling test)
  with at least one `.fits.gz` raw filename row so `SetupDataCube.main()` is exercised
  end-to-end with a compound extension.
- No changes needed to `test_setup_datacube_manual_validation` (unrelated to this fix).

## Verification

- `pytest pypeit/tests/test_setup_datacube.py -v` — confirm existing tests still pass and
  the new exact-match/fallback tests pass.
- `pytest pypeit/tests/test_outputfiles.py pypeit/tests/test_setup_datacube.py -v` — confirm
  no regressions in the `basename` branch's own tests.
- Manually re-run the real-world repro:
  `pypeit_setup_datacube /Users/westfall/Work/pypeit/kcwi_wcs/omeara/rdx/blue/keck_kcwi_D.pypeit "SDSSJ2222 2745"`
  and confirm it now finds and lists the existing
  `Science/spec2d_KB.*-SDSSJ22222745_KCWI_*.fits` files in the generated `.coadd3d` file via
  the primary exact-match path (not just the fallback).
- Cross-reference `kcwi_wcs.md`'s "Addressing inconsistencies" section with a short entry
  once implemented, recording that `setup_datacube.py` was integrated with
  `outputfiles.construct_basename`, the OMEARA real-data bug this fixes, and the corrected
  mjd-precision finding above.

## Open questions / inconsistencies uncovered

None beyond what's already documented above. Both instances of the underlying bug (the
real matching failure at `setup_datacube.py:237`, and the cosmetic log-message instance at
`setup_datacube.py:288`) trace to the same root cause (extension-agnostic `Path.stem` usage
predating the `basename` branch), and the mjd-precision concern that motivated keeping a
fallback path was investigated and found to be a non-issue for the data checked, not a
confirmed limitation — see the Investigation section above for why a fallback is retained
anyway, defensively, rather than because a real failure case was found.

## Implementation

### Item 1 — `pypeit/scripts/setup_datacube.py` (done)

Implemented as designed: added `from pypeit import outputfiles`;
`find_reduced_spec2d()` now takes `(science_dir, row, spectrograph)` and tries the exact
`outputfiles.construct_basename(...)`-derived filename first, falling back to the
(extension-bug-fixed) raw-stem-prefix glob + header-verify path when the exact file isn't
found or `mjd` isn't usable; `existing_spec2d_files()` now takes a `spectrograph` argument
and threads the row (not just its `filename`) through, and uses
`outputfiles.strip_raw_extension()` for the `missing`-list stems as well, fixing the
cosmetic instance of the bug too; `SetupDataCube.main()` now calls
`pypeit_file.get_spectrograph()` once and passes it down.

**Issue surfaced during implementation, not anticipated in the Design section above**:
the existing unit test's synthetic `.pypeit` file (`test_setup_datacube.py`'s
`_write_pypeit_file`) has no `mjd` column at all (`filename | frametype | target |
comb_id`) — the Design section's `row['mjd'] is not None` check assumed the column
exists and is merely possibly-`None`-valued, per `metadata.py`'s "corrupt header cards"
comment, but did not account for it being entirely absent. Against that test file, `row['mjd']`
raised `KeyError` rather than falling back cleanly. This matters beyond the test itself:
a minimal, hand-written `.pypeit` file (as opposed to one written by `run_pypeit` itself,
which always includes the core `mjd` metadata column) could plausibly omit it too. Fixed
with `mjd = row['mjd'] if 'mjd' in row.colnames else None`, so a missing column now
degrades to the fallback path instead of raising.

**Verification performed**:
- `pytest pypeit/tests/test_setup_datacube.py` — both pre-existing tests pass (both
  exercise the fallback path, since the test's synthetic `.pypeit` file has no `mjd`
  column and therefore never takes the primary exact-match path — this is itself a gap;
  see Item 2 below).
- Re-ran the real-world repro —
  `pypeit_setup_datacube keck_kcwi_D.pypeit "SDSSJ2222 2745"` against
  `/Users/westfall/Work/pypeit/kcwi_wcs/omeara/rdx/blue` — and confirmed the generated
  `.coadd3d` file now lists both real `Science/spec2d_KB.*-SDSSJ22222745_KCWI_*.fits`
  files, found via the primary exact-match path (confirmed by inspecting the written
  file's `spec2d read` block directly).

### Item 2 — `pypeit/tests/test_setup_datacube.py` (done)

Implemented as three new tests, added as *siblings* to the existing tests rather than by
modifying the shared `_write_pypeit_file`/`test_setup_datacube_write_and_append` fixture
(that fixture's hand-typed spec2d filenames, e.g. `spec2d_kr260610_00054-J0750+6927_KCRM_test.fits`,
don't match what `construct_basename` actually produces, and its rows have no `mjd`
column at all -- both are exactly why it already exercises the fallback path today, per
Item 1's note above; changing it to hit the primary path would have meant reworking
those hand-typed names and risked disturbing an already-passing append/overwrite/alias
CLI test that isn't otherwise related to this fix):

- `test_find_reduced_spec2d_exact_match_primary_path`: builds a one-row
  `astropy.table.Table` directly (`filename`/`target`/`mjd`) and a real
  `keck_kcrm` spectrograph via `load_spectrograph`, writes the spec2d file under the
  *exact* name `outputfiles.construct_basename(...)` computes for that row, and gives it
  a **deliberately wrong** `TARGET` header. `find_reduced_spec2d` still finds it,
  proving the primary path matches on the reconstructed name alone and never reads the
  header.
- `test_find_reduced_spec2d_fallback_path`: same setup, but the spec2d file is named
  with a timestamp token that does not match the row's `mjd` (simulating an
  mjd-precision miss), with a `.fits.gz` raw extension and a *correct* header `TARGET`.
  `find_reduced_spec2d` finds it only via the raw-stem-prefix glob + header-verify
  fallback. A second row with no `mjd` column at all is checked against the same file to
  confirm that path too. This is the direct regression test for the original bug: it was
  verified to fail against the pre-fix logic (see "Issue surfaced" note below) before
  being fixed.
- `test_setup_datacube_exact_match_end_to_end`: drives the full CLI (`SetupDataCube.main()`)
  against a `.pypeit` file with a `.fits.gz` raw filename and a real `mjd` column, with the
  matching spec2d file's `TARGET` header again deliberately wrong, and asserts the
  generated `.coadd3d` file lists it. Since the fallback path would reject a
  mismatched-header candidate, a passing run here can only be explained by the primary
  exact-match path succeeding -- closing the gap noted in Item 1 (no end-to-end coverage
  of the primary path).

**Issue surfaced during implementation**: to confirm the new fallback-path test
(`test_find_reduced_spec2d_fallback_path`) actually exercises the fix rather than passing
tautologically, `find_reduced_spec2d`'s `raw_stem` line was temporarily reverted to the
original `Path(str(raw_filename).strip()).stem` and the suite re-run. That test failed as
expected (`AssertionError: assert None == ...spec2d_kr260610_00058-...fits`), confirming
it does catch the original bug; the fix was then restored and the full suite re-verified
passing. No other issues surfaced -- the new tests, the pre-existing tests, and
`test_outputfiles.py` (16 tests total across both files) all pass.

**Verification performed**:
- `pytest pypeit/tests/test_setup_datacube.py -v` -- all 5 tests pass (2 pre-existing, 3
  new).
- `pytest pypeit/tests/test_setup_datacube.py pypeit/tests/test_outputfiles.py -v` -- 16
  tests, all pass; no regressions in the `basename` branch's own tests.
- Sanity-checked test validity by reverting the fix and confirming the new fallback test
  fails (see "Issue surfaced" above), then re-confirmed all tests pass with the fix
  restored.
