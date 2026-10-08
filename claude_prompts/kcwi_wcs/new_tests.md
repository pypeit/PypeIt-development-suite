# Datacube afterburn/vet test coverage + RA/DEC correctness fix

## Context

Two new IFU test datasets have been staged in the PypeIt development suite:
`keck_kcwi/large_bl` and `keck_kcrm/large_rl`, each containing two targets
(`SDSSJ2222 2745` and `gd50`, both typed `science` per KCWI/KCRM convention).
Both are already registered in `test_scripts/setups.py`'s `all_setups` dict, so
each already gets a basic `run_pypeit` "reduce" test for free. What's missing
is coverage of the three IFU datacube post-processing scripts
(`pypeit_setup_datacube` -> `pypeit_coadd_datacube` -> `pypeit_extract_datacube`),
which currently have **zero** dev-suite test coverage — confirmed by an
existing, never-implemented placeholder in `test_scripts/test_setups.py`:
`_coadd3d = {}` with a comment `# TODO: Test the pypeit_coadd_datacube setups!`.

While researching how to validate the "extracted gd50 coordinates match between
the two datasets" vet test, we discovered that `pypeit/core/datacube.py`'s
`extract_point_source()` hardcodes the extracted `SpecObj`'s `RA`/`DEC`
attributes to the cube's WCS `CRVAL` (the field-center/pointing reference) —
**never** the actual fitted-or-manual position that was used to build the
extraction aperture. This has been silently wrong for every point-source
extraction from a datacube, for any SlicerIFU instrument, since this feature
existed (confirmed: no existing test or downstream code depends on the current
value, so this is a clean, low-risk bugfix). The user asked to fix this as part
of this work, since it's necessary to make vet test 5 (matching celestial
coordinates across the two datasets) meaningful rather than circular.

This plan has three phases: (1) fix the RA/DEC bug in the main package, with a
permanent regression test; (2) add the three-script afterburn test chain to the
dev suite; (3) add five new vet tests that validate the afterburn chain's
output content. All work happens across two repos: the main package
(`/Users/westfall/Work/packages/pypeit-main/pypeit`) and the dev suite
(`/Users/westfall/Work/packages/pypeit-main/PypeIt-development-suite`).

---

## Phase 1 — Fix `SpecObj.RA`/`DEC` in cube point-source extraction

**File**: `pypeit/core/datacube.py`, function `extract_point_source()`.

Currently (lines ~1176-1180), `sobj.RA`/`DEC` are assigned immediately after
the `SpecObj` is created, using only the WCS reference point:
```python
sobj = specobj.SpecObj(_spectrograph.pypeline, "DET01", SLITID=0)
sobj.RA = wcscube.wcs.crval[0]
sobj.DEC = wcscube.wcs.crval[1]
```
This happens *before* `xobj, yobj` (the actual extraction-aperture center —
either the manual override or the auto 2D-Gaussian-fit peak, determined at
lines ~1220-1223) are even computed.

**Fix**: move the `RA`/`DEC` assignment to *after* `xobj, yobj` are determined,
and compute them via the WCS's celestial sub-transform:
```python
if manual_position is not None:
    xobj, yobj = manual_position
else:
    xobj, yobj = gaussian_position

skycoord = wcscube.celestial.pixel_to_world(xobj, yobj)
sobj.RA = skycoord.ra.deg
sobj.DEC = skycoord.dec.deg
```
- `xobj, yobj` are already 0-indexed, matching astropy's `pixel_to_world`
  convention (confirmed: the whitelight-image pixel grid used for the Gaussian
  fit runs `0` to `numxx-1`/`numyy-1`, and `manual_position` is documented and
  used identically to that same grid).
- No fallback/error-handling is needed for a failed fit: `fitGaussian2D()`
  already raises before returning if the fit doesn't converge or no source is
  found, so by the time this code runs, `xobj, yobj` are guaranteed valid.
- Confirmed no other code depends on the old (CRVAL) value: `grep` for
  `extract_point_source` callers found only `pypeit/coadd3d.py:328`; no
  production code reads `spec1d.RA`/`DEC` from a cube-extracted file to drive
  further automated logic (e.g. no auto standard-star lookup consumes it).

**Regression test**: add a new permanent test to `pypeit/tests/test_datacube.py`
verifying `extract_point_source()` assigns `RA`/`DEC` matching
`wcs.celestial.pixel_to_world(x, y)` for both the auto-fit case and the
manual-position case, using a small synthetic single-point-source cube (follow
the existing synthetic-cube-building helpers already in that file, e.g. the
one used by `test_extract_point_source_manual_position_selects_correct_source`
at line ~393, for a consistent, deterministic construction pattern). Assert
`SkyCoord` separation between the returned `sobj.RA/DEC` and the independently
computed `wcs.celestial.pixel_to_world(xobj, yobj)` is at the sub-mas level
(exact equality modulo floating point).

**Changelog**: record this as a bug fix (`update-changelog` skill) — this
changes the `RA`/`DEC` values written to every future `spec1d` file produced by
`pypeit_extract_datacube` / `DataCube.extract_spec()`, for any SlicerIFU
instrument, so it needs to be visible to users relying on those coordinates.

---

## Phase 2 — Afterburn test chain (dev suite)

**Files**: `test_scripts/pypeit_tests.py`, `test_scripts/test_setups.py`.

### 2.1 Three new `PypeItTest` subclasses (`pypeit_tests.py`)

Add near the other AFTERBURN classes (after `PypeItCoadd2DTest`). Each needs
only `__init__`/`build_command_line()` — the inherited `run()` already does
exit-code-only pass/fail via subprocess, matching every other AFTERBURN class
(e.g. `PypeItFluxSetupTest`, `PypeItFluxTest`).

- **`PypeItSetupDataCubeTest(setup, pargs, target)`**: runs
  `pypeit_setup_datacube <pyp_file> <target> -o`, where `pyp_file` is built via
  the existing `pypeit_file_name(instr, setup_name)` helper. `target` is passed
  as one argv element (no shell involved, so the embedded space in
  `'SDSSJ2222 2745'` needs no quoting). Give each instance a target-specific
  `log_suffix` (e.g. `test_setup_datacube_gd50`) so the two per-dataset
  invocations don't collide.
- **`PypeItCoaddDataCubeTest(setup, pargs, target)`**: runs
  `pypeit_coadd_datacube <coadd3d_file> -o`, where `coadd3d_file =
  self.setup.rdxdir / 'sources' / target_stub / f'{target_stub}.coadd3d'`
  (`target_stub = target.replace(' ', '')`, matching `setup_datacube.py`'s own
  convention). Have `build_command_line()` raise `FileNotFoundError` if that
  file doesn't exist yet — this belongs in `build_command_line()` (caught by
  the base `run()`'s try/except), **not** in `check_for_missing_files()`,
  since that method runs once at construction time, before any test's `run()`
  has executed, and would always see the file as missing.
- **`PypeItExtractDataCubeTest(setup, pargs, target)`**: runs
  `pypeit_extract_datacube <spec3d_file> -e <extract_file> -o`, where
  `spec3d_file = rdxdir/Science_cube/{target_stub}.fits` and `extract_file =
  rdxdir/sources/{target_stub}/{target_stub}.extract`. Same
  `build_command_line()`-level existence-check pattern as above for both
  required inputs.

### 2.2 Registration (`test_setups.py`)

Replace the placeholder
```python
# TODO: Test the pypeit_coadd_datacube setups!
_coadd3d = {}
```
with three dicts (no other file references `_coadd3d`, confirmed by repo-wide
grep, so this is a clean replacement):
```python
_setup_datacube = {
    'keck_kcwi': {'large_bl': [dict(target='SDSSJ2222 2745'), dict(target='gd50')]},
    'keck_kcrm': {'large_rl': [dict(target='SDSSJ2222 2745'), dict(target='gd50')]},
}
_coadd_datacube = {
    'keck_kcwi': {'large_bl': [dict(target='SDSSJ2222 2745'), dict(target='gd50')]},
    'keck_kcrm': {'large_rl': [dict(target='SDSSJ2222 2745'), dict(target='gd50')]},
}
_extract_datacube = {
    'keck_kcwi': {'large_bl': [dict(target='gd50')]},
    'keck_kcrm': {'large_rl': [dict(target='gd50')]},
}
```
The `[dict(...), dict(...)]`-per-setup, multi-invocation pattern is already
used elsewhere (e.g. `_quick_look`'s multi-dict entries), confirmed via
`build_test_setup()` in `test_main.py`, which instantiates one `PypeItTest` per
dict in the list and appends each to one flat `setup.tests` list — no new
harness capability needed.

Then add three new `TestPhase.AFTERBURN` entries to `all_tests`, in this
relative order (right where the old `_coadd3d` entry conceptually belonged):
`PypeItSetupDataCubeTest`/`_setup_datacube`, then
`PypeItCoaddDataCubeTest`/`_coadd_datacube`, then
`PypeItExtractDataCubeTest`/`_extract_datacube`. Because `all_tests`' list
order *is* execution order (per the existing code comment at
`test_setups.py:453-457`), and all three entries are appended after the
existing REDUCE-phase entry for these setups, this guarantees the correct
strict chain order per dataset: both `setup_datacube` invocations, then both
`coadd_datacube` invocations, then the one `extract_datacube` invocation — and
`thread_target()`'s existing fail-fast sequential loop already skips every
later step for a setup the moment one step fails, so no custom
dependency-skipping logic is needed.

Update the module docstring's list of private registration dicts to mention
the three new ones.

---

## Phase 3 — Vet tests (new file)

**File**: new `vet_tests/test_datacube_afterburn.py` (kept separate from the
existing `vet_tests/test_datacube.py`, which drives its own KCWI
`small_bh2_4200` chain end-to-end via the Python API rather than reading
already-produced afterburn output — different enough in style to warrant a
separate file).

Uses the standard `redux_out` fixture (`vet_tests/conftest.py`); no shared
path-builder helper exists in this repo's convention, so construct
`Path(redux_out) / instr / setup / 'Science'` (and `/ 'sources' / target_stub`)
inline per test, matching every other `vet_tests/*.py` file. Use plain
`assert expr, msg` (not `assert(expr, msg)` — the latter is a real,
always-truthy bug present in the existing `test_datacube.py`; do not
replicate it).

Five test functions, covering the five checks the user requested, run for
each of `(keck_kcwi, large_bl)` and `(keck_kcrm, large_rl)` (parametrize with
`@pytest.mark.parametrize` over `(instr, setup)`, matching this repo's plain
pytest-function convention):

1. **`test_coadd3d_files_reference_correct_spec2d`** — for both targets, read
   `sources/<target_stub>/<target_stub>.coadd3d` via
   `inputfiles.Coadd3DFile.from_file(...)` and assert its `data['filename']`
   set exactly equals the set of `Science/spec2d_*<target_stub>*fits` files
   actually on disk. Also sanity-check `output_filename == target_stub`,
   `combine == True`, `save_whitelight == True` in the parsed config.
2. **`test_datacube_shape`** — for both targets, open
   `Science_cube/<target_stub>.fits` via `coadd3d.DataCube.from_file(...)` and
   assert `flux`/`sig`/`bpm` are all 3D with identical shape
   `(nwave, ny, nx)` (per the documented datamodel axis order), `nwave/ny/nx
   > 1`, and a reasonable good-pixel fraction (not degenerate/all-bad).
3. **`test_whitelight_matches_cube`** — for both targets, open
   `Science_cube/<target_stub>_whitelight.fits` and assert its 2D shape equals the
   cube's `(ny, nx)`, and that its header's WCS (`CRVAL`/`CDELT`/`CTYPE` on the
   celestial axes) matches the cube's own celestial WCS, read directly from
   the cube FITS file's `FLUX` extension header via
   `astropy.wcs.WCS(header).celestial` (not the `DataCube` object's private
   `_wcs` attribute).
4. **`test_gd50_extraction_finds_single_source`** — for `gd50` only, load
   `Science_cube/spec1d_gd50_extract.fits` via `specobjs.SpecObjs.from_fitsfile`
   and assert exactly one `SpecObj` is present, with finite, non-degenerate
   `OPT_COUNTS`/`OPT_COUNTS_SIG` over a majority of wavelength pixels (these
   test datasets are not flux-calibrated, so `OPT_FLAM`/`OPT_FLAM_SIG` are
   `None`; `OPT_COUNTS`/`OPT_COUNTS_SIG` are what's actually populated).
5. **`test_gd50_coordinates_match_across_datasets`** — **not** parametrized
   per-setup (needs both datasets at once): read `RA`/`DEC` from
   `spec1d_gd50_extract.fits` in both `keck_kcwi/large_bl` and
   `keck_kcrm/large_rl`, build two `astropy.coordinates.SkyCoord` objects, and
   assert `coord_kcwi.separation(coord_kcrm) < 2 arcsec`. This is only a
   meaningful, non-circular check *because* Phase 1's fix makes these RA/DEC
   values reflect the actual detected source position rather than the WCS
   CRVAL; verified empirically on the real staged data to separate by ~1.6
   arcsec between the two instrument channels.

---

## Verification

The `reduce` step for these two datasets takes ~45 min each, versus ~1 min for
the datacube chain, and the user has already produced and staged the `reduce`
output in `REDUX_OUT/keck_kcwi/large_bl` and `REDUX_OUT/keck_kcrm/large_rl`.
Verification during development must reuse that existing output and must
**not** re-run `reduce` (the user will run the full `reduce afterburn vet`
pipeline themselves once this plan is implemented):

- **Phase 1**: `cd pypeit && python -m pytest -o addopts="" pypeit/tests/test_datacube.py --remote-data=none -q` — confirm the new regression test passes and no existing test regresses.
- **Phase 2**: from the dev suite, with `PYPEIT_DEV` set and `run_pypeit`/`pypeit_setup_datacube`/`pypeit_coadd_datacube`/`pypeit_extract_datacube` on `PATH`, run only the `afterburn` phase against the already-staged `reduce` output:
  `./pypeit_test afterburn -i keck_kcwi keck_kcrm -s large_bl large_rl` — confirm all 5+5=10 new subprocess-based tests pass, and check the produced `sources/`/`Science_cube/*.fits` files match the naming conventions above.
- **Phase 3**: `pytest vet_tests/test_datacube_afterburn.py --redux_out $PYPEIT_DEV/REDUX_OUT -v` — confirm all 5 (parametrized) vet functions pass against the Phase 2 output.
- Do not run the `reduce` phase (or `pypeit_test all`/`reduce afterburn vet` end-to-end) as part of this work's verification — that's the user's own follow-up step.

### Critical files
- `pypeit/pypeit/core/datacube.py` (Phase 1 fix)
- `pypeit/pypeit/tests/test_datacube.py` (Phase 1 regression test)
- `PypeIt-development-suite/test_scripts/pypeit_tests.py` (Phase 2 classes)
- `PypeIt-development-suite/test_scripts/test_setups.py` (Phase 2 registration)
- `PypeIt-development-suite/vet_tests/test_datacube_afterburn.py` (Phase 3, new file)

---

## Implementation status (2026-09-08)

All three phases implemented and verified. Deviations discovered while
implementing, versus what was assumed above:

- **`pypeit_coadd_datacube`/`pypeit_extract_datacube` write to `Science_cube/`,
  not `Science/`** — `CoAdd3D.output_paths()` always appends `_cube` to the
  science directory name (`Science` -> `Science_cube`). This document and the
  actual code were both corrected to use `Science_cube/<target_stub>.fits`,
  `Science_cube/<target_stub>_whitelight.fits`,
  `Science_cube/spec1d_gd50_extract.fits`, etc. The `Science/` directory only
  ever holds the *input* spec2d files.
- **These test datasets are not flux-calibrated** (no `sensfile` given to
  `pypeit_setup_datacube`), so the extracted `SpecObj`'s `OPT_FLAM`/
  `OPT_FLAM_SIG` are `None`; the vet test uses `OPT_COUNTS`/`OPT_COUNTS_SIG`
  instead, and a "good pixel" requires `OPT_COUNTS_SIG > 0` (not just finite),
  since out-of-range wavelength pixels have `OPT_COUNTS_SIG == 0` rather than
  NaN/inf.
- **`Coadd3DFile` config values are strings**, not Python bools (e.g.
  `combine == 'True'`, not `True`) — the vet test compares
  `.lower() == 'true'`.
- **Verified end-to-end for real** on an isolated scratch copy of both
  datasets' `REDUX_OUT` directories (not the user's staged copy): ran the full
  `pypeit_setup_datacube` -> `pypeit_coadd_datacube` -> `pypeit_extract_datacube`
  chain for real for both targets in both datasets, then ran all 5 (15,
  parametrized) vet tests against the real output — all pass. The real,
  independently-measured gd50 sky-position separation between the KCWI and
  KCRM extractions is ~1.6 arcsec, confirming both the RA/DEC fix and the vet
  test's 2 arcsec tolerance are sound.
- The scratch copy used for this verification was deleted after use; the
  user's original staged `REDUX_OUT/keck_kcwi/large_bl` and
  `REDUX_OUT/keck_kcrm/large_rl` directories were never modified.
