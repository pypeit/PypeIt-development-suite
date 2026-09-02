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
