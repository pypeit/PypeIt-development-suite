# Refactor `setup_datacube.py` into a thin script

## Context

`pypeit/scripts/setup_datacube.py` currently carries 13 module-level support functions
alongside its `SetupDataCube` script class — everything from target-name matching to
hand-rolled `.coadd3d`/`.extract` text-file writers. This doesn't match the codebase's
preferred convention (seen in `pypeit/scripts/setup_coadd2d.py`) of scripts being thin:
CLI parsing in `get_parser()`, and `main()` gathering args and calling logic that lives in
the package's real modules (`inputfiles`, `outputfiles`, `coadd2d`, etc.).

A full assessment (verified directly against the code, not just by inspection heuristics)
found that most of these 13 functions fall into one of three buckets: (1) they can be
**eliminated entirely** because an existing class/mechanism already does the job better,
(2) they are **genuinely novel** but belong with sibling logic in an existing module
(`inputfiles.py` or `outputfiles.py`), or (3) they are **CLI-argument-parsing helpers**
that legitimately belong in the script itself. No circular-import risk exists for the
moves proposed below (`inputfiles.py` and `outputfiles.py` currently import neither each
other; adding a one-directional `outputfiles → inputfiles` edge is safe — verified by
checking both modules' import lists and their transitive dependencies).

## Bucket 1 — Eliminated: replaced by existing `InputFile`/`ParSet` machinery

- **`write_coadd3d_file`** and **`append_spec2d_files`**: `inputfiles.Coadd3DFile`
  (`pypeit/inputfiles.py:1090`) is exactly the `.coadd3d` file type (`data_block =
  'spec2d'`, requires `[rdx] spectrograph`) and is currently hand-reimplemented as raw
  text lines instead of being used. Confirmed `InputFile.write()` (`inputfiles.py:555`)
  and `InputFile.from_file(..., preserve_comments=True)` round-trip a config's comments
  *and* any user-edited parameter values (via `ConfigObj`), and that the data table is
  regenerated from `self.data` (an `astropy.table.Table`) on every write — so appending a
  new filename is just adding a row to `.data` and calling `.write()` again. This directly
  satisfies `test_setup_datacube_write_and_append`'s requirement that a user's own edit to
  `weight_method` survive an `--append`.
- **`write_extract_file`**: `inputfiles.ExtractFile` (`inputfiles.py:951`) is exactly the
  `.extract` file type (`datablock_required = False`, i.e. parameter-block-only, already
  used this way by `extract_datacube.py`). `InputFile.write()` serializes whatever nested
  `config` dict is passed via `ConfigObj` — the same generic mechanism `Coadd2DFile` uses
  in `setup_coadd2d.py:214-215` — so no hand-formatted f-string lines are needed.
- **`validate_whitelight_range`**: `CubePar`/`CubeExtractionPar.validate()`
  (`pypeit/par/pypeitpar.py:1960-1961`, `:2116-2117`) already enforce
  `len(whitelight_range) == 2`; both classes call `self.validate()` from `__init__`
  (confirmed at `pypeitpar.py:2096`), so constructing `CubeExtractionPar(whitelight_range=
  [...])` raises for free. Only the CLI-string → `[float_or_None, float_or_None]`
  conversion is genuinely script-specific and stays local.
- **`validate_manual`**: `CubeExtractionPar.validate()` (`pypeitpar.py:2133-2142`) already
  rejects a `;` (multi-object) and enforces exactly one `:`-separated `x:y` pair — the
  *identical* restriction `validate_manual` reimplements, with the identical rationale
  documented in that method's own comment block. The one thing it does that `validate()`
  doesn't is confirm `x`/`y` are numeric; that check is cheap to keep as a couple of lines
  in the script (or left to fail later at `ManualCubeExtractionObj.parse()` time — to be
  decided during implementation, not this plan).

## Bucket 2 — Moved to existing modules

> **Superseded.** The `default_science_dir(redux_dir)` consolidation described below
> (and its corresponding `coadd2d.py`/`coadd3d.py`/`state/science_status.py` changes in
> "Files touched" and "Verification") was **not** implemented this way. See
> Implementation → Step 1's "Significant course-correction" for what was actually done
> instead, and why.

### `pypeit/outputfiles.py` (output-product discovery; already hosts `construct_basename`/
`strip_raw_extension`/`spec_output_file`/`science_path`)

- **`spec2d_target(spec2d_file)`** — unchanged.
- **`find_reduced_spec2d(science_dir, row, spectrograph)`** — unchanged; already calls
  `construct_basename`/`strip_raw_extension`, so this just becomes an intra-module call.
- **`existing_spec2d_files(pypeit_file, target, science_dir, spectrograph)`** — unchanged
  logic, but now calls `inputfiles.matching_science_rows`/`inputfiles.group_science_rows`
  (see below), which is the one new `outputfiles → inputfiles` import edge.
- **New `default_science_dir(redux_dir)`** — replaces `science_directory(pypeit_file)`.
  This also **consolidates the same inlined pattern** found at three other sites, all of
  which reduce to `<root>/'Science'` once you look past their differing inputs:
  - `coadd2d.py:520` / `coadd3d.py:790`: `Path(coadd_dir).absolute() / 'Science'` (inside
    an `if coadd_dir is not None:` branch — the sibling `else` branch derives the
    directory a different way, from `spec2d_files[0]`, and is untouched).
  - `state/science_status.py:567`: `Path(redux_dir) / 'Science'`.
  - `setup_datacube.py`'s own `science_directory(pypeit_file)`: `Path(pypeit_file).
    absolute().parent / 'Science'` — the odd one out, since it takes a **file**, not a
    **root directory**.
  The new helper's contract is root-directory-in: `default_science_dir(redux_dir) ->
  Path(redux_dir).absolute() / 'Science'`. `coadd2d.py`/`coadd3d.py`/`science_status.py`
  call it with the directory they already have; `setup_datacube.py`'s `main()` calls it
  with `Path(args.pypeit_file).absolute().parent`.

### `pypeit/inputfiles.py` (input-row selection; already hosts one module-level function,
`grab_rawfiles`, alongside the `InputFile` subclasses)

- **`target_match_key(target)`**, **`target_matches(value, target)`** — unchanged; no
  overlap found anywhere else in the package (checked `coadd2d.py`, `setup_coadd2d.py`,
  `ql.py` — all do cruder, inlined, less-permissive string comparisons, not a reusable
  function). Kept as two small functions rather than merged — normalize vs. compare are
  distinct, independently-testable responsibilities, and `target_match_key` may be useful
  on its own to future callers now that it lives in a shared module.
- **`matching_science_rows(pypeit_file, target)`**, **`group_science_rows(rows)`** —
  unchanged; genuinely novel (nothing else filters a bare `PypeItFile.data` table by
  frametype+target, or groups plain `astropy.table.Row` objects by `comb_id` with a
  filename fallback — the closest neighbors, `PypeItMetaData.get_frames_from_combid` and
  `ql.py`'s `quicklook_regroup`, operate on a full `PypeItMetaData`/fitstbl and
  assign/require `comb_id` rather than discover groups from a plain row list). Kept as two
  functions (select vs. group) since grouping is useful independent of target-filtering.

No changes are needed to `inputfiles.py`'s imports — `pypeitpar`, `spectrographs.util`,
and `utils` (its own dependencies) were all checked and none import `outputfiles.py`, so
the new `outputfiles → inputfiles` edge does not create a cycle.

## Bucket 3 — Stay in `setup_datacube.py` (CLI-parsing only)

Trimmed to genuinely CLI-specific string coercion once the semantic checks above move to
`CubePar`/`CubeExtractionPar`:
- A small `whitelight_range` CLI-string → `[float_or_None, float_or_None]` parser (the
  comma-split step `CubePar`/`CubeExtractionPar` don't do themselves).
- A small numeric sanity check on `--manual`'s `x`/`y` fields (the one thing
  `CubeExtractionPar.validate()` doesn't already check).

## `setup_datacube.py`'s new shape

Following `setup_coadd2d.py`'s style exactly (no dedicated `write_*`/`validate_*` helper
functions — `main()` builds config dicts inline and calls the `InputFile` subclasses
directly):

- `get_parser()`: unchanged (already just argparse setup).
- `main()`: parses/validates CLI strings inline, calls
  `inputfiles.matching_science_rows`/`inputfiles.group_science_rows` (via
  `outputfiles.existing_spec2d_files`), `outputfiles.default_science_dir`, builds a nested
  `cfg` dict for the cube-construction block (using `utils.add_sub_dict`, the same helper
  `setup_coadd2d.py:183` uses) and calls `inputfiles.Coadd3DFile(config=cfg,
  file_paths=[...], data_table=tbl).write(...)`; for `--append`, loads the existing file
  with `Coadd3DFile.from_file(coadd3d_file, preserve_comments=True)`, appends new rows to
  `.data`, and re-`.write()`s it; builds a second `cfg` dict for the extraction block and
  calls `inputfiles.ExtractFile(config=cfg).write(...)`.

## Files touched

> **The `coadd2d.py`/`coadd3d.py`/`state/science_status.py` bullet below is superseded**
> -- see the note atop Bucket 2 and Implementation → Step 1.

- `pypeit/outputfiles.py`: add `spec2d_target`, `find_reduced_spec2d`,
  `existing_spec2d_files`, `default_science_dir`; add `from pypeit import inputfiles` and
  `from astropy.io import fits`.
- `pypeit/inputfiles.py`: add `target_match_key`, `target_matches`,
  `matching_science_rows`, `group_science_rows` as module-level functions (alongside
  `grab_rawfiles`).
- `pypeit/coadd2d.py`, `pypeit/coadd3d.py`, `pypeit/state/science_status.py`: replace their
  inlined `Path(coadd_dir).absolute() / 'Science'` / `Path(redux_dir) / 'Science'` with
  `outputfiles.default_science_dir(...)`.
- `pypeit/scripts/setup_datacube.py`: reduced to `SetupDataCube` (`get_parser`/`main`) plus
  the two small CLI-parsing bits from Bucket 3.
- `pypeit/tests/test_setup_datacube.py`: keeps only script-level (`SetupDataCube.main()`
  end-to-end) tests. Tests for moved functions migrate to `pypeit/tests/test_outputfiles.py`
  (`spec2d_target`, `find_reduced_spec2d`, `existing_spec2d_files`, `default_science_dir`)
  and `pypeit/tests/test_inputfiles.py` (`target_match_key`, `target_matches`,
  `matching_science_rows`, `group_science_rows`).

## Verification

- `pytest pypeit/tests/test_setup_datacube.py pypeit/tests/test_inputfiles.py
  pypeit/tests/test_outputfiles.py -v` — all pass, including the migrated tests in their
  new locations.
- `pytest pypeit/tests/test_coadd2d.py pypeit/tests/test_coadd3d.py
  pypeit/tests/test_dashboard.py -v` (or whatever covers `science_status.py`) — confirm the
  `default_science_dir` swap in the other 3 sites is behavior-preserving.
- Re-run the real-world repro against the OMEARA KCWI reduction (as in prior sessions) to
  confirm `pypeit_setup_datacube` still finds and writes the expected `.coadd3d`/`.extract`
  files end-to-end after the refactor.
- Manually diff a freshly-generated `.coadd3d`/`.extract` file (via the new
  `Coadd3DFile`/`ExtractFile`-based writer) against one generated by the current
  hand-rolled writer, to confirm the `ConfigObj`-based serialization produces an
  equivalent, still-valid file (formatting may differ slightly; content must not).

## Implementation

### Step 1 — `outputfiles.py`/`inputfiles.py` additions (done)

Implemented `spec2d_target`, `find_reduced_spec2d`, `existing_spec2d_files` in
`pypeit/outputfiles.py` and `target_match_key`, `target_matches`,
`matching_science_rows`, `group_science_rows` in `pypeit/inputfiles.py`, all as designed
in Bucket 2 above (unchanged logic, moved as-is). `setup_datacube.py` itself was left
untouched in this step -- it still carries its own copies of these functions, to be
removed once the script itself is refactored in a later step. One deviation from the
plan text: `group_science_rows`'s internal `group_key` nested function was pulled out to
a module-level `science_row_group_key(row)` function instead, per a no-nested-functions
preference.

**Significant course-correction: the `default_science_dir(redux_dir)` consolidation
described in Bucket 2/Files touched above was not implemented as planned.** Investigating
it further turned up two things the original assessment missed:

1. `pypeit.spectrographs.spectrograph.Spectrograph.get_meta_value()` (`spectrograph.py:1701-1708`)
   already reads a config-specific value directly from a `.pypeit` file's own table row
   first, and only opens the raw file as a fallback if a needed column is missing. This
   means `PypeItFile.get_pypeitpar()` -- which the original setup_datacube.py deliberately
   avoided, assuming it required raw data to be present -- usually works fine without raw
   data at all. Given that, the better fix for `setup_datacube.py`'s own Science-directory
   lookup (per direction) is to call `get_pypeitpar()` and the already-existing
   `outputfiles.science_path(par)` (which correctly respects a customized
   `par['rdx']['scidir']`, unlike a hardcoded `'Science'` literal), rather than adding a
   new, cruder, `par`-independent helper.
2. To keep this robust for the rarer case where the raw file genuinely is required and
   missing, `PypeItFile.get_pypeitpar()` (and the `InputFile` base class version) gained a
   new `require_rawfile:bool=True` parameter. When `False`, a `PypeItError` from reading
   the config-specific raw file is caught and the configuration-specific parameters fall
   back to `spec.default_pypeit_par()`, so a fully populated `PypeItPar` is still
   returned either way. Default behavior (`True`) is unchanged for every existing caller.
3. The three other sites originally assumed to share one duplicated pattern turned out not
   to: `coadd2d.py`'s `output_paths` already had `par` in scope, so its hardcoded
   `'Science'` literal was simply changed to `par['rdx']['scidir']` (no new function
   needed). `coadd3d.py`'s `output_paths` already receives a `science_dir` parameter (the
   scidir *name*) but was hardcoding `'Science'` anyway, ignoring its own parameter -- an
   independent, pre-existing bug, fixed by using the parameter it already has.
   `state/science_status.py`'s `derive_science_from_disk` can be called with
   `fitstbl=None`, in which case no `par` is available at all; per direction, it was left
   untouched, since its constraint is genuinely different from the other three sites, not
   a shared duplicate.

Net effect: no `default_science_dir` function exists anywhere, and none is planned.
`setup_datacube.py`'s own Science-directory logic (not yet implemented -- that's part of
the still-pending script rewrite) will use `pypeit_file.get_pypeitpar(require_rawfile=False)`
followed by `outputfiles.science_path(par)` when that step happens.

**Verification performed**:
- Confirmed directly (not just via static analysis) that `outputfiles.py` and
  `inputfiles.py` import cleanly together with no circular import, after adding the new
  `outputfiles → inputfiles` edge.
- `pytest pypeit/tests/test_setup_datacube.py pypeit/tests/test_outputfiles.py
  pypeit/tests/test_inputfiles.py -v` -- 59 tests, all pass, including 7 new tests for the
  three `outputfiles.py` additions and 6 new tests for the four `inputfiles.py` additions
  plus the `require_rawfile` parameter (one of which revives the previously-disabled
  `test_get_pypeitpar_selects_science_file` scenario in `test_inputfiles.py`, using it to
  demonstrate the new fallback: `require_rawfile=True`, the default, still raises on a
  broken/empty raw file exactly as before; `require_rawfile=False` now succeeds).
- Confirmed `pypeit.scripts.setup_datacube`, `pypeit.coadd2d`, and `pypeit.coadd3d` all
  still import together cleanly after the `coadd2d.py`/`coadd3d.py` directory-lookup
  fixes.
- The `coadd2d.py`/`coadd3d.py` fixes were not exercised against a dev-suite reduction in
  this step (no `pypeit/tests/test_coadd2d.py`/`test_coadd3d.py` exists in the main repo);
  that verification is deferred, as originally planned, to the dev suite.

### Step 2 — `setup_datacube.py` rewrite (done)

Rewrote `pypeit/scripts/setup_datacube.py` down to the `SetupDataCube` script class
(`get_parser`/`main`) plus one small CLI-parsing helper, `_parse_whitelight_range(value)`
-- used directly as the `--whitelight_range` argparse `type=` callable, so
`args.whitelight_range` is already a validated `[min, max]` list (numeric entries or
`None`) by the time `main()` sees it, rather than a raw string parsed later. `main()`
follows `setup_coadd2d.py`'s style exactly: no dedicated `write_*`/`validate_*` helpers --
it builds plain nested `cfg` dicts inline and hands them to
`inputfiles.Coadd3DFile`/`inputfiles.ExtractFile` directly. All 13 of the original
module-level functions are gone from this file: the 8 moved in Step 1 are called via
`outputfiles.existing_spec2d_files`/`inputfiles.matching_science_rows` etc. as designed;
`write_coadd3d_file`/`append_spec2d_files`/`write_extract_file` are replaced by
`Coadd3DFile`/`ExtractFile` construction; `validate_whitelight_range` is now
`_parse_whitelight_range` plus `CubeExtractionPar`'s own `len==2` check;
`validate_manual`'s structural check (reject `;`, require one `x:y` pair) is now
`CubeExtractionPar`'s own `validate()`, with just the comma-rejection and numeric-cast
checks (the one thing `validate()` doesn't cover) kept as a few inline lines in `main()`.

**Two real bugs surfaced by prototyping the write/append logic in a scratch script before
touching `setup_datacube.py` itself** (rather than assuming the design from Step 1 would
just work):

1. **`whitelight_range` must be passed as an actual Python `list`, not a joined string.**
   An early prototype passed `'None,None'` (a string) into the `cfg` dict; `ConfigObj`
   quoted it on write (`whitelight_range = "None,None"`), which reads back as a single
   literal string rather than the two-element list `CubePar`/`CubeExtractionPar` require --
   silently breaking downstream parsing with no error at write time. Passing an actual
   list (`[None, None]` or `[9400.0, 10000.0]`) writes correctly as
   `whitelight_range = None, None` (unquoted, comma-separated), which reads back as
   `['None', 'None']` and coerces normally. This is exactly why `_parse_whitelight_range`
   returns a list rather than a re-joined string.
2. **`outputfiles.science_path(par)` alone is not enough -- it silently pointed at the
   current working directory, not the `.pypeit` file's directory**, when
   `par['rdx']['redux_path']` isn't explicitly set in the `.pypeit` file (the common
   case). This is the exact same problem `setup_coadd2d.py` already has a documented,
   two-step fix for (`pypeit/scripts/setup_coadd2d.py:117-127`): try
   `outputfiles.science_path(par)` first, then fall back to
   `pypeit_path.parent / par['rdx']['scidir']` if that directory doesn't exist. Reused
   that exact pattern rather than inventing a new one -- caught by running the existing
   `test_setup_datacube_write_and_append` test against the new code (it failed with
   `Expected Science directory does not exist: <cwd>/Science`), not by static review.

**Comment preservation**: the commented-out `# weights_init_obj_pos = x:y` example block
(asserted on by the pre-existing test) has no natural home in a plain dict-to-`ConfigObj`
construction, since comments only normally come from parsing an existing file. Confirmed
`configobj.Section.comments[key]` can be set programmatically on a freshly-built
`Coadd3DFile.config` object before `.write()`, and used it to attach the same four-line
comment block immediately before `weight_method` in the `[[cube]]` section, reproducing
the original file's helpful hint.

**Test file rewrite** (`pypeit/tests/test_setup_datacube.py`): trimmed to script-level
tests only, per the plan -- the `find_reduced_spec2d`/etc. tests added in Step 1 already
cover those functions in their new homes. `test_setup_datacube_write_and_append` and
`test_setup_datacube_manual_validation` were updated for the new (cosmetically different,
functionally equivalent) `ConfigObj`-based output format (e.g. `whitelight_range = None,
None` with a space, vs. the old hand-written `None,None`) and to exercise
`SetupDataCube.main()` directly instead of calling the now-removed `write_extract_file`.
Added a small `test_parse_whitelight_range` unit test, and extended the manual-validation
test to also cover the too-many-fields and non-numeric cases (previously untested).

**Verification performed**:
- `pytest pypeit/tests/test_setup_datacube.py pypeit/tests/test_outputfiles.py
  pypeit/tests/test_inputfiles.py -v` -- 57 tests, all pass.
- Repo-wide grep for the 13 removed function names (as `setup_datacube.<name>`
  references) and for `setup_datacube` in `doc/` -- no stale references found; the API
  doc page (`doc/api/pypeit.scripts.setup_datacube.rst`) is auto-generated via
  `automodule`/`:members:`, so it needs no manual update.
- Re-ran the real-world repro against the OMEARA KCWI reduction one more time, end to
  end, after the full rewrite: both real spec2d files are still found and written
  correctly to `.coadd3d`, and the `.extract` file content was inspected directly and
  confirmed correct.
- The `--append` and `-o`/`--overwrite` paths, including a simulated user edit to
  `weight_method` surviving `--append`, were validated both in a standalone prototype
  script (before editing `setup_datacube.py`) and via the updated
  `test_setup_datacube_write_and_append` test.

### Post-implementation review (done)

A line-by-line comparison of this document's original plan sections (Bucket 1-3, Files
touched, Verification) against the actual code found two things:

1. **Documentation-only**: Bucket 2, Files touched, and the original Verification section
   still describe the abandoned `default_science_dir` consolidation as the real plan.
   Added short pointer notes atop Bucket 2 and Files touched directing readers to
   Implementation → Step 1's course-correction, rather than rewriting the original plan
   text (kept as the historical record of what was proposed).
2. **A real, previously-unnoticed regression**: `--append` against a target with no
   existing `.coadd3d` file used to raise `PypeItError('Cannot append to missing .coadd3d
   file: {path}')`; after the Step 2 rewrite it instead raised a generic
   `FileNotFoundError('Input file {path} does not exist!')` from deep inside
   `inputfiles.InputFile.from_file`, since `main()`'s `--append` branch had no error
   handling around the missing-file case. Confirmed no test, before or after the rewrite,
   ever covered this scenario -- it slipped through both the original review and the test
   suite. **Fixed**: `main()` now checks `coadd3d_file.is_file()` before calling
   `Coadd3DFile.from_file(...)` and raises the original `PypeItError` message directly if
   it's missing. Added `test_setup_datacube_append_missing_coadd3d_file_raises` (new --
   this scenario had no test before this fix either).

Two additional, lower-priority items were noted but left as-is: the "manually diff a
freshly-generated file against the old writer" verification step was done via direct
inspection/assertions rather than a literal diff (functionally equivalent); and
`--spatial_delta`'s help text claims a `0.678924` default for `keck_kcrm` that no code
(old or new) has ever implemented -- a pre-existing inaccuracy, unrelated to this
refactor, not addressed here.

**Verification performed**: `pytest pypeit/tests/test_setup_datacube.py
pypeit/tests/test_outputfiles.py pypeit/tests/test_inputfiles.py -v` -- 58 tests, all
pass (57 from Steps 1-2, plus the one new regression test above).
