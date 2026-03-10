# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

The PypeIt-development-suite is the testing and development framework for [PypeIt](https://github.com/pypeit/PypeIt), a Python spectroscopic data reduction pipeline. It provides raw test data references, configuration files, and tests for 60+ astronomical instruments.

**Required environment variable:** `export PYPEIT_DEV=/path/to/PypeIt-development-suite`

Raw data (`RAW_DATA/`) and calibrations (`CALIBS/`) are hosted on Google Drive and symlinked into the repo — they are not committed to git.

## Running Tests

### Full test suite (via custom runner)
```bash
$PYPEIT_DEV/pypeit_test all                          # All test phases
$PYPEIT_DEV/pypeit_test reduce vet                   # Reductions + verification
$PYPEIT_DEV/pypeit_test all -i shane_kast_blue       # Single instrument
$PYPEIT_DEV/pypeit_test all -t 4                     # 4 parallel threads
$PYPEIT_DEV/pypeit_test all -r report.txt            # With report
$PYPEIT_DEV/pypeit_test list                         # List all instruments/setups
```

### Test phases (in execution order)
1. **pypeit_tests** — PypeIt's own pytest suite (CI-compatible)
2. **unit** — Dev-suite unit tests (require RAW_DATA)
3. **reduce** — Full reduction tests via `run_pypeit`
4. **afterburn** — Post-reduction tools (flux cal, coadding, telluric)
5. **ql** — Quick-look tests
6. **vet** — Verification pytest tests on reduction outputs

### Running pytest directly
```bash
cd $PYPEIT_DEV
pytest unit_tests                              # All unit tests
pytest vet_tests                               # All vet tests
pytest unit_tests/test_scripts.py              # Single test file
pytest vet_tests --redux_out /custom/path      # Custom output dir
```

The `--redux_out` fixture (defined in `vet_tests/conftest.py`) defaults to `$PYPEIT_DEV/REDUX_OUT`.

## Architecture

### Key directories
- **`test_scripts/`** — Test orchestration engine
  - `setups.py` — Master dict of all instruments → setup names
  - `test_setups.py` — Defines which setups run which test types, test phase ordering, and how to add new tests (read its module docstring)
  - `pypeit_tests.py` — `PypeItTest` subclasses that execute each test type as child processes
  - `test_main.py` — Test execution engine with parallel scheduling
- **`pypeit_files/`** — `.pypeit` configuration files (one per instrument/setup)
- **`unit_tests/`** — pytest tests requiring RAW_DATA
- **`vet_tests/`** — pytest tests that validate reduction outputs in REDUX_OUT
- **`pypeitdev/`** — Development modules: wavelength templates, filter curves, flux calibration, sensitivity functions
- **`nautilus/`** — Kubernetes job generation for distributed testing on Nautilus cluster

### Adding a new instrument/setup
1. Add raw data to Google Drive under `RAW_DATA/instrument/setup`
2. Add the instrument and setup to `test_scripts/setups.py` (`all_setups` dict)
3. Add a `.pypeit` file to `pypeit_files/`
4. Register test types in `test_scripts/test_setups.py` (see its docstring for details)
5. Test: `$PYPEIT_DEV/pypeit_test reduce -i <instrument> -s <setup>`

### Test priority
`test_priority_list` orders setups slowest-to-fastest for optimal parallel scheduling. It is auto-regenerated on full suite passes and should be periodically committed.

## Headless testing
Set `QT_QPA_PLATFORM=offscreen` or `source $PYPEIT_DEV/source_headless_test.sh` for environments without a display.
