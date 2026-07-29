# PypeIt-development-suite

This repository (`PypeIt-development-suite/`) provides the data, input files, and
test harness used to validate [PypeIt](https://github.com/pypeit/PypeIt), the
pure-Python spectroscopic data-reduction pipeline.  It is **not** a Python
package users install for processing; it is a testing/integration repository that
exercises a working PypeIt installation end-to-end against real raw data from
many supported spectrographs.

## Relationship to the main PypeIt repository

- The main code repository is `PypeIt/` (cloned as a sibling directory here at
  `/home/xavier/Projects/PypeIt/PypeIt`).  It has its own `CLAUDE.md` describing
  the architecture (MultiSlit / Echelle / SlicerIFU paths, `Spectrograph`
  classes, `PypeItPar`, `DataContainer`, etc.).  Consult it when a test failure
  points at pipeline behavior rather than test setup.

- This dev suite **assumes PypeIt is already installed** (e.g. `pip install -e
  ".[dev]"` in the `PypeIt/` repo) so that `run_pypeit` and the other
  `pypeit/scripts` entry points are on `PATH`.

- The two repos are versioned independently and are developed in tandem: a
  change to PypeIt that affects reduction results often requires a matching
  change to a `.pypeit` file, expected output, or test here.  Branches are
  frequently coordinated by name (see Nautilus `-p`/`-d` options below).

- A subset of the tests run here are the **same pytest unit tests** that live in
  `PypeIt/pypeit/tests/` (invoked via the `pypeit_tests` test type); the rest
  are integration/regression tests unique to this repo.

## Prerequisites and environment

- **`PYPEIT_DEV`** environment variable must point to the top level of this repo.
  Both `pypeit_test` and PypeIt's own data-dependent tests read it.

- **`RAW_DATA/`** and **`CALIBS/`** are large and are **not** committed to git
  (see `.gitignore`).  They are hosted on a
  [Google Drive](https://drive.google.com/drive/folders/1oh19siB1-F0jjmY-F_jr73eA-TQYEiFW?usp=sharing)
  and are typically symlinked into this directory.  Here, `RAW_DATA` is a symlink
  to an external drive (`/media/xavier/SamsungT7/RAW_DATA`).  Never add raw data,
  FITS, or reduction outputs to the repo.

- `REDUX_OUT/` is the default output directory for reductions and is also
  git-ignored.

- For headless/GUI tests, `source source_headless_test.sh` to set
  `QT_QPA_PLATFORM=offscreen` (QT5/QT6 must be installed).

## Running the tests

The primary entry point is the `pypeit_test` script (a thin wrapper around
`test_scripts/test_main.py`).  It requires `run_pypeit` on `PATH` and `PYPEIT_DEV`
set.

```console
# Run everything (in order: pypeit_tests, unit, reduce, afterburn, ql, vet)
$PYPEIT_DEV/pypeit_test all

# Run only specific test types
$PYPEIT_DEV/pypeit_test reduce afterburn ql

# Restrict to instruments and/or setups
./pypeit_test reduce -i shane_kast_blue shane_kast_red
./pypeit_test reduce -s shane_kast_blue/600_4310_d55 shane_kast_red/600_7500_d57

# List every supported instrument/setup
./pypeit_test list

# Parallelism, reports, performance, coverage
./pypeit_test -t 2 all
./pypeit_test all -r test_report.txt
./pypeit_test all --csv performance.csv
./pypeit_test all --coverage coverage_report_file.txt
```

### Test types

| Test type      | What it does |
|----------------|--------------|
| `pypeit_tests` | The pytest unit tests bundled with PypeIt (`PypeIt/pypeit/tests/`). Self-contained; CI-safe. |
| `unit`         | Pytest tests in `unit_tests/` that **require** `RAW_DATA`. |
| `reduce`       | Full reductions that call `run_pypeit` directly on the `.pypeit` files. |
| `afterburn`    | Post-reduction tools: flux calibration, 1D/2D coadd, telluric, sensfunc, collate, etc. |
| `ql`           | Quick-Look reduction tests. |
| `vet`          | Pytest tests in `vet_tests/` that verify the products created by `reduce`/`afterburn`. |
| `all`          | Runs all of the above, in the listed order (later phases depend on earlier ones). |
| `list`         | Lists supported instruments/setups; runs nothing. |

The `unit_tests` and `vet_tests` directories can also be driven directly with
pytest; `vet_tests` expects a `--redux_out` option (defaults to
`$PYPEIT_DEV/REDUX_OUT`, see `vet_tests/conftest.py`).

```console
cd $PYPEIT_DEV
pytest unit_tests
pytest vet_tests --redux_out /path/to/REDUX_OUT
```

## Repository layout

- **`pypeit_test`** — top-level test runner (wrapper for `test_scripts.test_main`).
- **`test_scripts/`** — the test framework:
  - `test_setups.py` — **the single source of truth** for what the dev suite
    contains: the `all_setups` dict (instrument → setups), `supported_instruments`,
    and the per-test-type setup lists (`telluric_tests`, `quick_look_tests`, etc.)
    plus `all_tests`, which defines test ordering and phases (PREP, REDUCE,
    AFTERBURN, QL). **Edit this file to add/modify instruments, setups, or tests.**
  - `pypeit_tests.py` — the `PypeItTest` subclasses that launch each kind of
    test as a child process; add a subclass here to define a new test type.
  - `test_main.py` — argument parsing, orchestration, threading, reporting.
  - `setups.py`, `test_pypeit_test.py`, `test_setups.py` — supporting/self-test code.
- **`RAW_DATA/`** (symlink, not committed) — raw frames organized as
  `RAW_DATA/<instrument>/<setup>/`, e.g. `RAW_DATA/shane_kast_blue/600_4310_d55/`.
- **`pypeit_files/`** — `.pypeit` reduction input files used by `reduce` tests, one
  per setup, named `<instrument>_<setup>.pypeit`.
- **Afterburn input files**, each named `<instrument>_<setup>.<ext>`:
  - `coadd1d_files/` (`.coadd1d`), `coadd2d_files/` (`.coadd2d`),
    `fluxing_files/` (`.flux`), `sensfunc_files/` (`.sens`),
    `tellfit_files/`, `flexure_files/`.
- **`unit_tests/`** — pytest tests needing `RAW_DATA` (calibrations, metadata,
  frame typing, mosaics, telluric, setup GUI, scripts, …).
- **`vet_tests/`** — pytest tests that vet reduction products (coadd2d, datacube,
  echelle, edgetrace, extraction, flexure, flux, sensfunc, skysub, wavelengths,
  wavetilts, …); `vet_tests/files/` holds reference/fixture data.
- **`sensfunc_archive/`** — archived sensitivity functions.
- **`nautilus/`** — Kubernetes/Nautilus cluster job generation (`gen_kube_devsuite`),
  Docker definitions, and S3 report retrieval config for running the full suite
  remotely.
- **`pypeitdev/`** — loose development scratch/analysis notebooks and scripts
  organized by topic (echelle, flats, flexure, fluxing, wavelengths, soar_tspec,
  …). Not part of the automated test suite.
- **Helper scripts**: `build_ql_calibs` (build Quick-Look calibrations),
  `pypeit_syncraw` (sync raw data), `test_priority_list` (slowest→fastest setup
  ordering for parallel scheduling; regenerated on a full passing run), and
  `fix_xshooter_names.py`.

## Adding to the dev suite

To add a new instrument/setup or test (per the header docs in
`test_scripts/test_setups.py`):

1. Place the raw data under `$PYPEIT_DEV/RAW_DATA/<instrument>/<setup>/` (and on
   the Google Drive), plus any needed `.pypeit`/afterburn input files in the
   corresponding directory above.
2. Add the instrument/setup to `all_setups` in `test_scripts/test_setups.py`, and
   add new instruments to `supported_instruments`.
3. Add `instrument/setup` entries to the relevant per-test-type lists for any
   additional tests desired.
4. To add a wholly new *type* of test: add a `PypeItTest` subclass in
   `pypeit_tests.py`, create a setup list for it, and register it in `all_tests`
   with its factory, phase, and setups.

## Conventions and constraints

- **Deterministic tests**: random-number generators must use fixed seeds (same
  rule as the main PypeIt repo).
- **No large/binary data in git**: raw data, FITS, PNG, JSON, YAML, logs, and
  notebooks are git-ignored (see `.gitignore`); a few wavelength template files
  under `pypeitdev/wavelengths/template_files/` are explicit exceptions.
- **Naming**: setup-specific input files follow `<instrument>_<setup>.<ext>`; this
  mapping is relied on by the test framework.
- The full suite takes 10+ hours and tens of GiB of RAM serially; use `-t` for
  parallelism (memory scales roughly per the table in `README.rst`). The pytest
  phases cannot run in parallel.
- Python 3.11–3.12 (see `setup.cfg`).

## Documentation

- This repo's user-facing instructions live in `README.rst`.
- The development-suite workflow is also documented in the main repo at
  `PypeIt/doc/dev/development.rst` and on
  [ReadTheDocs](https://pypeit.readthedocs.io/).
