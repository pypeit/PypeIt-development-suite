# Extract `construct_basename` into `pypeit/outputfiles.py`

## Context

`PypeItMetaData.construct_basename` (`pypeit/metadata.py:484-505`) builds the
root name used for most PypeIt output files. It currently requires a full
`PypeItMetaData` instance (row-indexed table access) even though the actual
string-building logic only needs four primitive values: the raw filename, the
target name, the camera name, and the observation MJD. This blocks reuse of
the basename logic anywhere a full metadata table isn't available.

Separately, the current implementation strips the raw file extension with a
hard-coded `self['filename'][row].split('.fits')[0]`. This only works for
extensions containing the literal substring `.fits` (`.fits`, `.fits.gz`,
`.fits.bz2`). It silently fails to strip the extension for spectrographs whose
`Spectrograph.allowed_extensions` are different, e.g. `aat_uhrf.py` (`.FTS`),
`soar_goodman.py` (`.fz`), `mdm_modspec.py`/`wht_isis.py` (`.fit`,
`.fit.gz`). This work fixes that bug as part of the extraction.

The goal is to move the string-building logic into a standalone function in
`pypeit/outputfiles.py` (same function name, `construct_basename`), have
`PypeItMetaData.construct_basename` become a thin wrapper around it, and have
the new function accept the spectrograph's actual `allowed_extensions` list so
extension-stripping is correct for every supported instrument.

## Design

### 1. New function in `pypeit/outputfiles.py`

Add (near the other basename/path helpers, e.g. above `spec_output_file`):

```python
def construct_basename(filename, target, camera, mjd, allowed_extensions):
    """
    Construct the root name primarily for PypeIt file output.

    Args:
        filename (:obj:`str`):
            The name of the raw file (may include a compound extension,
            e.g. ``.fits.gz``).
        target (:obj:`str`):
            The name of the target/object observed.
        camera (:obj:`str`):
            The camera name (e.g. :attr:`~pypeit.spectrographs.spectrograph.Spectrograph.camera`).
        mjd (:obj:`float`):
            The MJD of the observation.
        allowed_extensions (:obj:`list`):
            List of recognized raw-file extensions for the relevant
            spectrograph (e.g.
            :attr:`~pypeit.spectrographs.spectrograph.Spectrograph.allowed_extensions`),
            used to correctly strip the extension from ``filename``.

    Returns:
        :obj:`str`: The root name for file output.
    """
    root = filename
    for ext in sorted(allowed_extensions, key=len, reverse=True):
        if filename.endswith(ext):
            root = filename[:-len(ext)]
            break
    tiso = time.Time(mjd, format='mjd', out_subfmt='date_hms')
    dtime = datetime.datetime.strptime(tiso.isot, '%Y-%m-%dT%H:%M:%S.%f')
    return '{0}-{1}_{2}_{3}{4}'.format(root, target.replace(' ', ''), camera,
                                        datetime.datetime.strftime(dtime, '%Y%m%dT'),
                                        tiso.isot.split('T')[1].replace(':', ''))
```

Requires two new imports at the top of `outputfiles.py`: `import datetime` and
`from astropy import time`. `allowed_extensions` is a required argument (no
default) so the base `Spectrograph.allowed_extensions` default
(`['.fits', '.fits.gz']`) isn't duplicated/drifted here — every real caller
already has a spectrograph object to pull it from.

### 2. `PypeItMetaData.construct_basename` becomes a thin wrapper

In `pypeit/metadata.py`, replace the body (keep the same signature/docstring
so all 5 existing call sites — `pypeit_steps.py:104`, `outputfiles.py:120`,
`pypeit.py:398`, `science_status.py:377,416` — need no changes):

```python
def construct_basename(self, row, obstime=None):
    """
    Construct the root name primarily for PypeIt file output.
    ...
    """
    _obstime = self.construct_obstime(row) if obstime is None else obstime
    return outputfiles.construct_basename(self['filename'][row], self['target'][row],
                                           self.spectrograph.camera, _obstime.mjd,
                                           self.spectrograph.allowed_extensions)
```

Add `from pypeit import outputfiles` to `metadata.py`'s import block
(alongside the existing `from pypeit import inputfiles`, `pypeit/metadata.py:22`).

### 3. Dependency-tree / circular-import check

This adds one new edge: `metadata.py -> outputfiles.py`. Verified safe:

- `outputfiles.py`'s only intra-package imports are `from pypeit import log`
  and `from pypeit import PypeItError` — both names are fully defined by
  `pypeit/__init__.py` (lines 20 and 29) using only `pypeit/pkg/logger.py`,
  `pypeit/pkg/pypeitdata.py`, and `pypeit/pkg/exceptions.py`, none of which
  import `metadata.py` or `outputfiles.py`.
- `outputfiles.py` does not import `metadata.py`, `spectrographs/*`, or
  anything that transitively imports `metadata.py`.
- `pypeit/__init__.py` itself does not import `metadata.py` or
  `outputfiles.py` at package-init time.

So the new edge is one-directional (`metadata.py` depends on
`outputfiles.py`, never the reverse) — no cycle is introduced.

### 4. No other call sites need changes

All 5 real call sites and the `test_dashboard.py` mock stub call
`fitstbl.construct_basename(...)` on a `PypeItMetaData` instance; the wrapper
keeps that interface identical. `outputfiles.spec_output_file` (which itself
calls `fitstbl.construct_basename(frame)` at line 120) is left as-is — it
already has a `fitstbl` object handy and doesn't need to call the new
standalone function directly.

## Verification

- Add a unit test in `pypeit/tests/test_outputfiles.py` (new file, since none
  exists yet) exercising `outputfiles.construct_basename` directly with
  plain arguments — no `PypeItMetaData` needed — covering: a plain `.fits`
  filename, a `.fits.gz` filename, and an extension that doesn't contain
  the substring `.fits` (e.g. `.fz` or `.fit`) to confirm the
  extension-stripping fix. Use a fixed MJD per repo testing conventions.
- Run `pytest pypeit/tests/test_outputfiles.py pypeit/tests/test_dashboard.py -v`
  to confirm the new test passes and the existing `_FakeFitstbl` mock-based
  dashboard tests are unaffected.
- The user will separately run the broader unit test suite
  (`pytest pypeit/tests/`) to catch any regression in code that indirectly
  depends on `construct_basename` (`pypeit_steps.py`, `pypeit.py`,
  `state/science_status.py`) — no action needed here.
- Sanity-check in a Python shell that
  `outputfiles.construct_basename('file.fits.gz', 'M31', 'CAM', 58000.5, ['.fits', '.fits.gz'])`
  and the equivalent `.fz`/`.fit` cases produce the expected root name with
  the extension fully stripped.

## Summary of Completed Work

### Main repo (`pypeit`)

- `pypeit/outputfiles.py`: added the standalone `construct_basename()`
  function as designed above. It strips the extension using
  `allowed_extensions` (longest-suffix-first), accepts `str` or `Path`
  input via `Path(filename).name`, and logs a warning (without raising) if
  no extension matches, rather than silently leaving the input unmodified.
- `pypeit/metadata.py`: `PypeItMetaData.construct_basename` is now a thin
  wrapper delegating to `outputfiles.construct_basename`, with the
  now-unused `datetime` import removed. Verified this produces identical
  output to the pre-refactor implementation.
- `pypeit/tests/test_outputfiles.py` (new): unit tests exercising the new
  function directly with plain arguments, covering `.fits`, `.fits.gz`, a
  non-`.fits`-containing extension, `Path` input, and the
  unrecognized-extension warning path.

**Bug found and fixed along the way:** `pypeit/spectrographs/soar_goodman.py`
declared `allowed_extensions = ['.fz']`, but real SOAR Goodman raw files are
always named `*.fits.fz`. The old `construct_basename` stripped `.fits.fz`
correctly only by the accident of its naive `.split('.fits')[0]` logic; once
extension-stripping was switched to be driven by `allowed_extensions`, this
mismatch surfaced as a real regression (basenames retained a trailing
`.fits`). Fixed by changing `allowed_extensions` to
`['.fits', '.fits.fz']` (the `.fits` entry also guards against users who
have already decompressed their files before running PypeIt). Updated the
one existing test that asserted the old value
(`pypeit/tests/test_spectrographs.py`).

All touched main-repo test modules pass: `test_outputfiles.py`,
`test_metadata.py`, `test_spectrographs.py`, and `test_dashboard.py` (92
tests total).

### Dev suite (`PypeIt-development-suite`)

- `unit_tests/test_construct_basename.py` (new): follows the pattern in
  `unit_tests/test_load_images.py` (module-level `test_*` functions, a
  small `grab_*` helper, real raw files under `RAW_DATA`). Each test
  instantiates the real spectrograph, reads the file header via
  `spectrograph.get_headarr`/`get_meta_value`, and calls
  `outputfiles.construct_basename` directly — no `PypeItMetaData` table
  needed. Files were chosen from each setup's `.pypeit` file to be actual
  uncommented `science`-frame rows, and together cover every distinct raw
  extension pattern used across the spectrographs in the codebase:

  | Spectrograph            | Extension    | Raw file (under `RAW_DATA/`)                              |
  |--------------------------|-------------|------------------------------------------------------------|
  | `aat_uhrf`               | `.FTS`      | `aat_uhrf/3875/RUN0034.FTS`                                 |
  | `wht_isis_blue`          | `.fit.gz`   | `wht_isis_blue/long_R300B_d5300/r2324653.fit.gz`             |
  | `gemini_gnirs_echelle`   | `.fits`     | `gemini_gnirs_echelle/10_LB_SXD/N20160717S0114.fits`         |
  | `gemini_gmos_south_ham`  | `.fits.gz`  | `gemini_gmos/GS_HAM_R400_700/S20181005S0079.fits.gz`         |
  | `gemini_gmos_south_ham`  | `.fits.bz2` | `gemini_gmos/GS_HAM_R400_795_SENS/S20221207S0183.fits.bz2`   |
  | `soar_goodman_red`       | `.fits.fz`  | `soar_goodman_red/M2/0320_FRB210320_host_05-04-2021.fits.fz` |
  | `mdm_modspec`            | `.fit`      | `mdm_modspec/Echelle/MDM_Science_Test.fit`                   |

  All 7 tests pass and assert the exact expected basename string, including
  confirming the raw extension is fully stripped.

## Follow-up: `coadd2d.default_basename` and other basename call sites

### Issue

`pypeit.coadd2d.CoAdd2D.default_basename` (`pypeit/coadd2d.py:445-473`) has
the same class of bug as the original `construct_basename`: it strips the
raw-file extension with a hard-coded

```python
frsthdr['FILENAME'].split('.fits')[0]
```

`frsthdr`/`lasthdr` are headers read directly off `spec2d` output files, and
the `FILENAME` card in them is populated from `row_fitstbl['filename']` when
the `spec2d` file is written (`Spectrograph.initialize_header`,
`pypeit/spectrographs/spectrograph.py:475-476` sets header key `filename`,
which is the same case-insensitive FITS card as `FILENAME`). So this is the
raw file's original name, extension and all, and `.split('.fits')` fails the
same way as before for `.FTS`, `.fz`, `.fit`-style raw extensions.

### Broader search

Searching the main repo for `basename`, `.split('.fits')`, and
`io.remove_suffix` turned up a second, related root cause:
`pypeit.io.remove_suffix` (`pypeit/io.py:205-233`) only special-cases `.gz`
as a "strip two suffix levels" case; for any other compound extension
(`.fits.bz2`, `.fits.fz`) it only strips the last suffix, leaving a stray
`.fits` behind. Four call sites use it on **raw** filenames (as opposed to
PypeIt's own always-`.fits` output products) and are affected the same way:

| Call site | Context |
|---|---|
| `pypeit/pypeit_steps.py:685` (`load_skyregions`) | `io.remove_suffix(scifile)`, where `scifile = fitstbl.frame_paths(frames[0])` (a raw file path) |
| `pypeit/images/rawimage.py:798` (`RawImage`, spatial flexure QA) | `io.remove_suffix(self.filename)`, the raw exposure file |
| `pypeit/scripts/skysub_regions.py:95` | `io.remove_suffix(file_base)`, derived from the `FILENAME` card in a `spec2d` header (same raw-name provenance as the coadd2d bug above) |
| `pypeit/scripts/ql.py:104-121` (`folder_name_from_scifiles`) | `io.remove_suffix(...)` over a list of raw science files |

For all four, the failure only shows up for GMOS `.fits.bz2` and SOAR
Goodman `.fits.fz` raw files (single-suffix extensions like `.FTS`/`.fit`
and the `.gz`-suffixed ones are already handled correctly by
`Path.stem`/the existing `.gz` special case).

Four more call sites use the equally fragile `.split('.fits')` pattern, but
operate on PypeIt's own output products (always literally `.fits`, never
compressed), so they are not exposed to the spectrograph-extension bug.
They're still worth cleaning up (see Design, item 3, below), since the same
brittle string-splitting pattern is being retired everywhere else in this
effort: `pypeit/coadd1d.py:183`, `pypeit/scripts/chk_noise_2dspec.py:227`,
`pypeit/scripts/chk_noise_1dspec.py:284`, `pypeit/scripts/sensfunc.py:215`.

Also out of scope: `pypeit/state/science_status.py`'s several `basename`
uses already route through the now-fixed `fitstbl.construct_basename`.

**Related but different bug, not fixed here:**
`pypeit/scripts/multislit_flexure.py:87` does
`header['FILENAME'].split('.')` and indexes `fnames[2]` to pull out a
frame-ID token, assuming a fixed dot-delimited raw-naming convention. This
is a different failure mode (a naming-convention assumption, not a
mis-stripped extension) and is left alone.

### Design

Three parts: (1) harden the generic `io.remove_suffix` utility itself as
defense-in-depth, (2) route every genuinely spectrograph-dependent call
site through spectrograph-aware stripping (consistent with the
`outputfiles.construct_basename` fix above), and (3) retire the remaining
`.split('.fits')` occurrences that touch PypeIt's own (always-`.fits`)
output products.

**1. Generalize `pypeit/io.py:remove_suffix`.** Currently it only
special-cases `.gz` as a double-suffix. Widen the special-cased suffix set
so it also handles `.bz2` and `.fz`:

```python
def remove_suffix(file):
    _compression_suffixes = ('.gz', '.bz2', '.fz')
    _file = Path(file)
    while _file.suffix in _compression_suffixes:
        _file = _file.with_suffix('')
    return _file.stem
```

This is behavior-preserving for every existing case (verified by tracing
`.txt`, `.fits`, `.fits.gz`, `dot.separated.file.name.txt`) and additionally
fixes `.fits.bz2` -> `...` and `.fits.fz` -> `...` generically, without
needing a `Spectrograph` instance. It remains a generic utility, though, so
it still doesn't know about spectrograph-specific single-token extensions
like `.FTS` in a case where those genuinely need spectrograph-level
handling -- hence part 2 below for the call sites that need that.

**2. Route the spectrograph-dependent raw-filename call sites through
`allowed_extensions`.** Since every affected call site already has, or can
cheaply obtain, the relevant `Spectrograph` instance: **Done** -- all six
sub-items below are implemented. Verified: all five touched modules
(`coadd2d.py`, `pypeit_steps.py`, `images/rawimage.py`,
`scripts/skysub_regions.py`, `scripts/ql.py`) import cleanly (no circular
imports); `coadd2d.CoAdd2D.default_basename` smoke-tested against real
`spec2d` output (`ldt_deveny/DV1`); `folder_name_from_scifiles`
smoke-tested with synthetic GMOS `.fits.bz2` and SOAR Goodman `.fits.fz`
filenames, both correctly stripped with no stray extension remnants;
`pytest pypeit/tests/test_ql.py pypeit/tests/test_outputfiles.py
pypeit/tests/test_io.py pypeit/tests/test_metadata.py
pypeit/tests/test_dashboard.py pypeit/tests/test_spectrographs.py` --
100 tests, all passing (no dedicated unit tests exist yet for `coadd2d.py`,
`pypeit_steps.py`, `rawimage.py`, or `skysub_regions.py` in the main repo).
The full `pytest pypeit/tests/` run is left to the user.

   a. **`pypeit/outputfiles.py`**: extract the extension-stripping loop
      currently inlined in `construct_basename` into its own function,
      `strip_raw_extension(filename, allowed_extensions)`, and have
      `construct_basename` call it. Same behavior (longest-suffix-first
      match, `Path(filename).name` handling, warn-and-return-unmodified
      fallback), just reusable on its own.

   b. **`pypeit/coadd2d.py`**: in `default_basename`, after the existing
      `FILENAME`/`TARGET` presence checks, also require `PYP_SPEC` (raise
      `PypeItError` if missing, matching the existing style), then:
      ```python
      spec = load_spectrograph(frsthdr['PYP_SPEC'])
      first = outputfiles.strip_raw_extension(frsthdr['FILENAME'], spec.allowed_extensions)
      last = outputfiles.strip_raw_extension(lasthdr['FILENAME'], spec.allowed_extensions)
      return f"{first}-{last}-{frsthdr['TARGET'].replace(' ', '')}"
      ```
      Add `from pypeit import outputfiles` and
      `from pypeit.spectrographs.util import load_spectrograph` to the
      module's imports (verified no circular import: neither module
      imports `coadd2d.py`, directly or transitively).

   c. **`pypeit/pypeit_steps.py`**: in `load_skyregions`, replace
      `io.remove_suffix(scifile)` with
      `outputfiles.strip_raw_extension(scifile, spectrograph.allowed_extensions)`
      (the function already takes `spectrograph` as its first argument).
      Add `from pypeit import outputfiles`; remove `from pypeit import io`,
      which becomes unused (it's only used at this one call site in the
      file).

   d. **`pypeit/images/rawimage.py`**: replace
      `io.remove_suffix(self.filename)` with
      `outputfiles.strip_raw_extension(self.filename, self.spectrograph.allowed_extensions)`
      (`RawImage` already carries `self.spectrograph`). Add
      `from pypeit import outputfiles`; remove `from pypeit import io`
      (likewise its only use in the file).

   e. **`pypeit/scripts/skysub_regions.py`**: `specname = hdr['PYP_SPEC']`
      is already read (line 69). Load the spectrograph
      (`load_spectrograph(specname)`) and replace
      `io.remove_suffix(file_base)` with
      `outputfiles.strip_raw_extension(file_base, spec.allowed_extensions)`.
      Add `from pypeit import outputfiles` and
      `from pypeit.spectrographs.util import load_spectrograph` to the
      existing lazy-import block in `main()`; remove the now-unused
      `from pypeit import io`.

   f. **`pypeit/scripts/ql.py`**: `folder_name_from_scifiles` currently
      only takes `sci_files`; add an `allowed_extensions` parameter (plain
      list, consistent with the primitives-first design used elsewhere in
      this effort) and use `outputfiles.strip_raw_extension` in place of
      `io.remove_suffix`. Update its one call site (line 274) to pass
      `ps_sci.fitstbl.spectrograph.allowed_extensions`. `outputfiles` is
      already imported in this file; remove the now-unused
      `from pypeit import io`.

**3. Retire the remaining `.split('.fits')` sites.** These touch PypeIt's
own output products, so no `Spectrograph` lookup is needed -- but a plain
`io.remove_suffix`/`Path.stem` swap is only safe where the original code
didn't rely on the path's directory being preserved. **Done** -- both
sub-items below are implemented as designed. Verified: all four modules
(`sensfunc.py`, `coadd1d.py`, `chk_noise_2dspec.py`, `chk_noise_1dspec.py`)
import cleanly; a grep of the main repo for `split('.fits')` now returns no
matches outside `pypeit/tests/`; manually confirmed the directory-preserving
behavior is unchanged for the `coadd1d`/`chk_noise_*` cases and the
basename-only behavior is unchanged for the `sensfunc` case. No dedicated
main-repo unit tests exist for these four scripts (their coverage lives in
the dev suite's `vet_tests/`), so nothing needed updating there. Checked
each individually:

   - **`pypeit/scripts/sensfunc.py:215`**: `_names` are already
     directory-stripped at line 212 (`Path(f).name`), so
     `_names[0].split('.fits')[0]` can become `Path(_names[0]).stem`
     directly -- no behavior change, and `Path.stem` is simpler than
     `io.remove_suffix` here since these are always plain `spec1d_*.fits`
     files (never compressed).

   - **`pypeit/coadd1d.py:183`**, **`pypeit/scripts/chk_noise_2dspec.py:227`**,
     **`pypeit/scripts/chk_noise_1dspec.py:284`**: in all three, the
     string being split (`coaddfile`, `file`, `file`) is a full/relative
     path that may include a directory (a user-supplied output path in
     `coadd1d.py`; a CLI-supplied file argument in the `chk_noise_*`
     scripts), and the result is used to build a sibling output
     path/folder next to the original file. **Neither `Path.stem` nor
     `io.remove_suffix` is appropriate here** -- both discard the parent
     directory (`Path.stem`/`remove_suffix` only ever return the final
     path component), which would silently relocate the output file/folder
     to the current working directory instead of alongside the input.
     This was the key finding of exploring this cleanup: the correct,
     directory-preserving replacement is `Path(file).with_suffix('')`
     (as a string, e.g. `str(Path(coaddfile).with_suffix(''))`), not
     `io.remove_suffix`.

### Verification

- [x] Extend `pypeit/tests/test_io.py::test_remove_suffix` with `.fits.bz2`
  and `.fits.fz` cases to lock in the generalized compression-suffix
  handling. **Done**: `io.remove_suffix` generalized to loop over
  `('.gz', '.bz2', '.fz')`; `test_remove_suffix` extended with both cases;
  docstring examples updated and manually verified.
- [x] Extend `pypeit/tests/test_outputfiles.py` with tests for the new
  `strip_raw_extension` function directly (already exercised indirectly
  through the existing `construct_basename` tests, but direct cases make
  the contract explicit). **Done**: `strip_raw_extension` extracted out of
  `construct_basename` in `pypeit/outputfiles.py` (which now just calls
  it), with 6 direct unit tests added, including a `.fits.bz2` case.
  `pytest pypeit/tests/test_io.py pypeit/tests/test_outputfiles.py
  pypeit/tests/test_metadata.py pypeit/tests/test_dashboard.py` -- 87
  tests, all passing. (Design items 2 and 3 above were completed in
  subsequent work, including `coadd2d.py`.)
- [x] Add a regression test to the dev suite's `vet_tests/` directory (in
  `vet_tests/test_coadd2d.py`, alongside the existing coadd2d tests),
  exercising `coadd2d.CoAdd2D.default_basename` against real `spec2d`
  output. **Done**, once both datasets were reduced: two tests were added,
  `test_default_basename_single_spec2d` (`soar_goodman_blue/M1`, a single
  `spec2d` file, raw extension `.fits.fz` -- the same pattern that exposed
  the `soar_goodman.py` bug above, and exercises the `first == last`
  branch) and `test_default_basename_multi_spec2d` (`ldt_deveny/DV1`,
  three `spec2d` files, exercising the `first != last` branch). Both were
  run against the real reduced output
  (`pytest vet_tests/test_coadd2d.py::test_default_basename_single_spec2d
  vet_tests/test_coadd2d.py::test_default_basename_multi_spec2d
  --redux_out <REDUX_OUT>`) and pass, confirming the expected basenames
  with no stray `.fits`/`.fz` remnants.
- Run `pytest pypeit/tests/` to confirm no regressions from the removed
  `io` imports, the generalized `remove_suffix`, and the
  `folder_name_from_scifiles` signature change (check for any existing
  callers/tests of that function beyond `ql.py:274`).
- The user will separately run relevant parts of the development suite
  that exercise `coadd2d`, `load_skyregions`, `RawImage`, `coadd1d`,
  `chk_noise_1dspec`/`chk_noise_2dspec`, and quick-look (`ql.py`)
  end-to-end.

### Additional cleanup: drop `os` from the `chk_noise_*` scripts

While touching `pypeit/scripts/chk_noise_1dspec.py` and
`chk_noise_2dspec.py` in item 3 above, `import os` was left in place for
`os.path.exists`/`os.makedirs`. Since `Path` was already newly imported in
both files, this was replaced with the `pathlib` equivalent --
`folder = Path(...); folder.mkdir(parents=True, exist_ok=True)` -- removing
the `os` import entirely from both scripts. Behavior-preserving (parent
directories still get created, and re-running with an existing folder still
doesn't raise); the `folder` variable remains safe to use downstream in the
existing `'{}/...'.format(folder, ...)` path-joining, since `Path` stringifies
correctly there. Both modules import cleanly with no `os` references left.
