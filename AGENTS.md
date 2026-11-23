# PypeIt-development-suite

This repository (`PypeIt-development-suite/`) provides the data, input files, and
test harness used to validate [PypeIt](https://github.com/pypeit/PypeIt), the
pure-Python spectroscopic data-reduction pipeline.  It is **not** a Python
package users install for processing; it is a testing/integration repository that
exercises a working PypeIt installation end-to-end against real raw data from
many supported spectrographs.

## Relationship to the main PypeIt repository

- The main code repository (https://github.com/pypeit/PypeIt) is typically
  cloned as a sibling directory of this one; the clone's directory name varies
  between developers.  It has its own `AGENTS.md` describing its architecture,
  coding conventions, and development workflow.  Consult it when a test failure
  points at pipeline behavior rather than test setup.

- This dev suite **assumes PypeIt is already installed** (e.g. `pip install -e
  ".[dev]"` in the main repo) so that `run_pypeit` and the other
  `pypeit/scripts` entry points are on `PATH`.  The supported Python versions
  are those of the PypeIt installation; see `pyproject.toml` in the main repo.

- The two repos are versioned independently and are developed in tandem: a
  change to PypeIt that affects reduction results often requires a matching
  change to a `.pypeit` file, expected output, or test here.  Branches are
  frequently coordinated by name (see the Nautilus branch options under
  "Repository layout").

- A subset of the tests run here are the **same pytest unit tests** that live in
  the main repo's `pypeit/tests/` (invoked via the `pypeit_tests` test type);
  the rest are integration/regression tests unique to this repo.

## Prerequisites and environment

- **`PYPEIT_DEV`** environment variable must point to the top level of this repo.
  Both `pypeit_test` and PypeIt's own data-dependent tests (and its full
  documentation build) read it.

- **`RAW_DATA/`** is large and is **not** committed to git (see `.gitignore`).
  It is hosted on a
  [Google Drive](https://drive.google.com/drive/folders/1oh19siB1-F0jjmY-F_jr73eA-TQYEiFW?usp=sharing)
  and is placed (or symlinked) at the top level of this repo.  Never add raw
  data, FITS files, or reduction outputs to the repo.

- Git-ignored working directories: `REDUX_OUT/` (the default output directory
  for reductions) and `QL_CALIB/` (the default output directory of
  `build_ql_calibs`).

- For headless/GUI tests, `source source_headless_test.sh` to set
  `QT_QPA_PLATFORM=offscreen` (a Qt backend must be installed; PypeIt installs
  `pyqt6` by default).

## Running the tests

**`README.rst` is the authority on how to run the tests**: the `pypeit_test`
script (`test_scripts/test_main.py`), its test types and options, running the
`unit_tests` and `vet_tests` directly with pytest, reports, coverage,
performance statistics, parallel and headless testing, and running on
Nautilus.  Read it before running or changing the tests.

**Keep `README.rst` consistent with the code.**  Any change to the test
framework (`pypeit_test`, `test_scripts/`, `vet_tests/conftest.py`,
`nautilus/gen_kube_devsuite`, `source_headless_test.sh`, etc.) must be
verified against `README.rst`, and the README updated in the same change if it
no longer matches; e.g., the `pypeit_test -h` listing must match the current
output of `./pypeit_test -h`.  Note that `README.rst` is also included in the
PypeIt documentation on ReadTheDocs.

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
  - `setups.py`, `test_pypeit_test.py` — supporting/self-test code.
- **`RAW_DATA/`** (not committed) — raw frames organized as
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
- **`wavecalib/`** — wavelength-calibration tests, with reference data in
  `wavecalib/TEST_DATA/`.
- **`sensfunc_archive/`** — archived sensitivity functions.
- **`nautilus/`** — Kubernetes/Nautilus cluster job generation
  (`gen_kube_devsuite`), Docker definitions, and S3 report retrieval config for
  running the full suite remotely.  `gen_kube_devsuite` selects the branches to
  test with `-p/--pypeit_branch` and `-d/--dev_branch`; note that `-p` means
  something different (`--prep_only`) in `pypeit_test`.
- **`pypeitdev/`** — loose development scratch/analysis notebooks and scripts
  organized by topic (echelle, flats, flexure, fluxing, wavelengths, …). Not
  part of the automated test suite.
- **`.claude/skills/`** — shared agent skills for PypeIt development (e.g.
  `add-devsuite-setup`, `run-dev-suite`, `new-spectrograph`, `add-parameter`),
  each in `<name>/SKILL.md`.  Skills live here regardless of which repository
  they operate on.
- **`claude_prompts/`** — a common space for developers to share agent prompts
  and planning documents, kept as a historical record of agent-aided
  development activity.  Documents are grouped in subdirectories named for the
  topic or the git branch the work is on (e.g. `claude_prompts/bspline_refactor/`).
- **Helper scripts**: `build_ql_calibs` (build Quick-Look calibrations; writes to the
  directory given by the `QL_CALIB` environment variable, or `./QL_CALIB` if it
  is not set),
  `test_priority_list` (slowest→fastest setup ordering for parallel scheduling;
  regenerated on a full passing run), `fix_xshooter_names.py`, and
  `pypeit_syncraw` (syncs `RAW_DATA` via rclone; its remotes are hard-coded per
  developer, so it only works for the developers it lists).

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

## When changing PypeIt, also update the dev suite

- **A change that alters reduction results, parameters, or input-file
  formats**: update the affected `.pypeit`/afterburn input files and the
  `vet_tests` that check the changed products.
- **A new spectrograph, or a new mode of an existing one**: add a dev-suite
  setup (see "Adding to the dev suite" above).
- **A renamed or removed PypeIt script, option, or API**: update the callers in
  `test_scripts/pypeit_tests.py`, `unit_tests/`, and `vet_tests/`.

## Coding Conventions

The code in this repo falls into two tiers:

- **The test framework and tests** (`test_scripts/`, `unit_tests/`,
  `vet_tests/`, `wavecalib/`) follow the conventions below.
- **Scratch areas** (`pypeitdev/`, `claude_prompts/`, and the notebooks in
  `nautilus/`) have no style requirements.

As in the main repo, **these conventions apply to new code.**  Much of the
existing code does not yet follow them; it is migrated gradually, and only when
a developer asks for it.  Do not reformat or restyle code that you are not
otherwise changing.

### Inherited from the main repo

Follow the "Coding Conventions" section of `AGENTS.md` in the main PypeIt repo
(formatting, logic and numerics, structure and imports, and American English),
with these exceptions:

- **Docstrings:** only the NumPy style and the "opening `"""` alone on its
  line" rules apply, and only in `test_scripts/`.  The Sphinx-specific rules
  (79-character wrapping, `r"""`, intersphinx roles, `ln`/`log`) do not apply,
  because this repo has no Sphinx documentation.  Tests do not require
  docstrings.
- **Not required:** the rules on module-level globals, `import ... as`, and
  attribution notes for adapted code.

### Specific to the dev suite

- **Tests are plain module-level functions**, not `TestX` classes.
- **Tests are deterministic**: random-number generators use fixed seeds.
- **Every `assert` carries a message describing the expected behavior.**
- **Test data:** raw data goes in `RAW_DATA` (and on the Google Drive); small
  reference files live next to the tests that use them (e.g.
  `vet_tests/files/`, `wavecalib/TEST_DATA/`).  This replaces the main repo's
  rule that test data never goes under a tests directory.
- **Sort directory listings** (`sorted(Path.glob(...))`) even more strictly
  than in the main repo: tests that walk `REDUX_OUT` or `RAW_DATA` must behave
  the same on macOS, Linux, and Nautilus.
- **When test code uses PypeIt's API:**
  - select frames with the `PypeItMetaData` query methods (`find_frames`,
    `find_frame_files`, ...) rather than operating on `fitstbl.table`
    directly;
  - use `from pypeit import log` and `PypeItError`, not the removed `msgs`
    interface;
  - do not call PypeIt's private (`_`-prefixed) functions.  The two repos are
    versioned independently, so such calls break silently when PypeIt is
    refactored.

## Constraints

- **No large/binary data in git**: raw data, FITS, PNG, JSON, YAML, TXT, logs,
  and notebooks are git-ignored (see `.gitignore`); a few wavelength template
  files under `pypeitdev/wavelengths/template_files/` are explicit exceptions.
  Because these patterns are ignored, a new file of one of these types will not
  appear in `git status`; a small file that genuinely must be tracked has to be
  added with `git add -f`.
- **Naming**: setup-specific input files follow `<instrument>_<setup>.<ext>`; this
  mapping is relied on by the test framework.
- **Runtime and memory**: the full suite is long-running and memory-intensive;
  see `README.rst` before running it.

## Documentation

- This repo's user-facing instructions live in `README.rst`.
- The development-suite workflow is also documented in the main repo at
  `doc/dev/development.rst` and on
  [ReadTheDocs](https://pypeit.readthedocs.io/).
