---
name: add-devsuite-setup
description: Add a new instrument/setup or a new test to the PypeIt development suite — stage raw data and input files, then register it in test_scripts/test_setups.py (and pypeit_tests.py for a new test type). Use when adding dev-suite coverage for a new spectrograph, grating mode, or post-reduction tool.
---

# Add a setup or test to the development suite

This skill operates on **this** repo (`PypeIt-development-suite/`).

## Reference material

- `test_scripts/test_setups.py` — **single source of truth**; its header
  docstring documents exactly how to add instruments, setups, tests, and test
  types. Read it first.
- `test_scripts/pypeit_tests.py` — the `PypeItTest` subclasses (one per test
  type) that run each test as a child process.
- `README.rst` — overall layout and the test-type table.

## Add a new instrument/setup (reduce coverage)

1. **Stage the data and input files** (none of this is committed except the
   small input files):
   - Raw frames: `$PYPEIT_DEV/RAW_DATA/<instrument>/<setup>/` (and on the
     Google Drive).
   - Reduction input: `pypeit_files/<instrument>_<setup>.pypeit`. **Note the
     filename is lower-cased**: the harness builds it as
     `f'{instr}_{setup.lower()}.pypeit'` (so setup `TSPEC` →
     `<instrument>_tspec.pypeit`). The `path` line inside the file is rewritten
     to the real RAW_DATA dir at run time, so its literal value does not matter;
     the data-block filenames must match what you staged.
   - Afterburn inputs as needed, same naming convention:
     `coadd1d_files/*.coadd1d`, `coadd2d_files/*.coadd2d`,
     `fluxing_files/*.flux`, `sensfunc_files/*.sens`, `tellfit_files/*`,
     `flexure_files/*`.

2. **Register the setup:**
   - Add the setup to the `all_setups` dict under its instrument. This dict
     lives in **`test_scripts/setups.py`** (a simple `instrument -> [setups]`
     map); `test_scripts/test_setups.py` imports it and auto-derives the default
     `reduce` coverage (`_reduce_setups`), so a basic reduce test runs
     automatically with no further edit. (There is no longer a
     `supported_instruments` list in code, despite older docstrings — instruments
     are just the `all_setups` keys.)
   - For extra tests, add `'<instrument>/<setup>'` to the relevant per-test-type
     lists/dicts in `test_setups.py` (e.g. `telluric_tests`, `quick_look_tests`,
     coadd/flux/sensfunc tests). When a test needs arguments, use the dict form
     (e.g. `{'coadd_file': 'pisco_coadd.fits', 'redshift': 7.52, ...}`).
   - If the raw-data directory does not follow the `RAW_DATA/<instrument>/<setup>`
     convention, add an entry to `_raw_data_dirs` in `test_setups.py`.

3. **Load-image unit test** (the convention for "supported" instruments): add an
   entry to `unit_tests/test_load_images.py`. Not all instruments have one (e.g.
   `p200_tspec` does not), so mirror whatever the analog instrument you are
   copying does.

## Add a new test type

1. Add a `PypeItTest` subclass in `pypeit_tests.py` that launches the test as a
   child process.
2. Create a setup list with at least one `instrument/setup` that runs it.
3. Register it in the `all_tests` list with its `factory`, `type`
   (`TestPhase.PREP|REDUCE|AFTERBURN|QL`), and `setups` — placed in the correct
   run order (later tests can depend on earlier outputs).

## Verify

- `./pypeit_test list` shows the new instrument/setup.
- Run it narrowly (see the `run-dev-suite` skill):
  ```console
  ./pypeit_test reduce -s <instrument>/<setup>
  ./pypeit_test all -s <instrument>/<setup>      # incl. afterburn/vet if added
  ```
- Confirm the small input files (not raw/FITS data) are the only things staged
  for commit (`git status`); raw data, FITS, logs are git-ignored.
