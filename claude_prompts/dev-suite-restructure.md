# Restructure PypeIt-development-suite: separate testing from development material

## Context

The dev-suite repo mixes two things at its top level: (1) the regression/integration
test harness and its inputs (`test_scripts/`, `unit_tests/`, `vet_tests/`, `wavecalib/`,
seven `*_files/` input directories, `test_priority_list`, `build_ql_calibs`, `nautilus/`),
and (2) development-only material (`pypeitdev/` — 1554 tracked files of notebooks and
analysis scripts, `sensfunc_archive/`, `claude_prompts/`, `reports/`, `deprecated/`, ...).
The goal is to put all testing code and inputs under one top-level directory, put the
development material under another, and keep `./pypeit_test` at the top level.

Decisions already made:
- Testing directory is **`tests/`**; development material goes under a sibling **`dev/`**.
- **`RAW_DATA` stays at the top level** (along with the other git-ignored data/working
  dirs `Telluric`, `REDUX_OUT`, `QL_CALIB`).  It is shared by the tests, the main repo's
  docs build (`doc/scripts/make_example_files.py`), and the dev notebooks, and every
  developer's symlink and the Nautilus rclone step point there.  `PYPEIT_DEV` keeps
  meaning "top of the dev-suite repo".

## Proposed layout

```
PypeIt-development-suite/
├── pypeit_test                  # unchanged location; imports tests.framework.runner
├── pypeit_syncraw               # unchanged location; syncs RAW_DATA (also top level)
├── README.rst  AGENTS.md  CLAUDE.md  CHANGES.rst  .gitignore
├── .claude/skills/              # unchanged (shared agent skills)
├── RAW_DATA -> ...              # git-ignored, unchanged
├── Telluric -> ...  REDUX_OUT/  QL_CALIB/      # git-ignored, unchanged
│
├── tests/                       # EVERYTHING the test suite runs or reads
│   ├── __init__.py
│   ├── conftest.py              # shared fixtures: dev_root, raw_data_dir, input_file_dir
│   ├── framework/               # was test_scripts/
│   │   ├── __init__.py
│   │   ├── paths.py             # NEW: the one place dev-suite paths are defined
│   │   ├── runner.py            # was test_main.py
│   │   ├── setups.py            # unchanged: all_setups (instrument → RAW_DATA setup dirs)
│   │   ├── registry.py          # was test_setups.py (imports setups; all_tests, ...)
│   │   ├── pypeit_tests.py      # unchanged name (PypeItTest subclasses)
│   │   ├── test_runner.py       # was test_pypeit_test.py (self-tests of the harness)
│   │   └── test_priority_list   # was top-level; data file read/written by the runner
│   ├── input_files/             # was the seven top-level *_files/ dirs
│   │   ├── pypeit/              # was pypeit_files/
│   │   ├── coadd1d/  coadd2d/  flux/  sensfunc/  tellfit/  flexure/
│   ├── unit/                    # was unit_tests/  (+ __init__.py)
│   ├── vet/                     # was vet_tests/   (+ __init__.py; conftest.py, files/ kept)
│   ├── wavecalib/               # was wavecalib/   (+ __init__.py; TEST_DATA/ kept)
│   ├── tools/
│   │   ├── build_ql_calibs      # invoked by the QL tests
│   │   └── source_headless_test.sh
│   └── nautilus/                # was nautilus/ (remote execution of the suite)
│
└── dev/                         # development-only material, not run by any test
    ├── pypeitdev/               # was pypeitdev/ (topic notebooks/scripts, unchanged inside)
    ├── sensfunc_archive/        # scripts that build archived sensfuncs (no test uses them)
    ├── claude_prompts/          # agent prompts and plans
    ├── reports/
    ├── deprecated/              # build_cooked
    ├── fix_xshooter_names.py
    └── pyproject.toml  setup.cfg  setup.py   # only package pypeitdev; move with it
```

Rationale for the less obvious choices:
- **Rename `test_main.py`/`test_setups.py`**: under `tests/`, pytest would try to collect
  every `test_*.py`.  Only real test modules keep the `test_` prefix.
- **`__init__.py` in `tests/`, `framework/`, `unit/`, `vet/`, `wavecalib/`**: `unit_tests`
  and `vet_tests` both contain `test_scripts.py` and `test_slitmask.py`, which collide if
  `pytest tests` is run as one session; packages give them distinct module names.  A
  regular `tests` package also cannot be shadowed by a stray `tests` package installed in
  site-packages (a namespace package can be).  With pytest's default `prepend` import
  mode, the repo root goes on `sys.path`, so tests can `from tests.framework import paths`.
- **`input_files/<kind>/`**: drops the redundant `_files` suffix and groups all inputs;
  file names (`<instrument>_<setup>.<ext>`) are unchanged.
- **`test_priority_list` into `framework/`**: it is harness state, not a user-facing tool.

## Needed decisions

1. **Make the `wavecalib` tests independent of the working directory.**
   `wavecalib/test_longslit.py:15` and `wavecalib/test_ThAr.py:17` build the path to
   their reference data from the current directory
   (`str(Path().absolute()) + '/TEST_DATA/'`), so they only find `TEST_DATA/` when pytest
   is launched from inside `wavecalib/`.  After the move that directory becomes
   `tests/wavecalib/`, and `pytest tests/wavecalib` run from the top level would fail.
   - *Proposed:* resolve it relative to the test file,
     `Path(__file__).resolve().parent / 'TEST_DATA'`, which works from any directory.
   - *Alternative:* leave the code alone and document that these tests must be run from
     `tests/wavecalib/`.

2. **Remove the tracked `REDUX_OUT_C/*.md` files.**  `REDUX_OUT_C/arcfit_metrics.md` and
   `REDUX_OUT_C/arcfit_metrics_A.md` were committed in `9e1b6d3b` (profxj, 2026-10-26,
   "mo") together with `pypeitdev/speed_up/` profiling changes.  They appear to be
   scratch output from that arc-fitting/speed-up work, and nothing in the repo
   references them.  `REDUX_OUT/` is git-ignored, but `REDUX_OUT_C/` is not.
   - *Proposed:* ask profxj whether they are still needed.  If so, move them to
     `dev/pypeitdev/speed_up/` next to the scripts that produced them; if not, delete
     them.  Either way, add `REDUX_OUT*/` to `.gitignore`.
   - *Alternative:* leave them in place; they are outside both `tests/` and `dev/`, so the
     top level would still contain non-test, non-data material.

## Implementation steps

1. **Moves with `git mv`** (preserves history) for every row of the layout above.
2. **Add `tests/framework/paths.py`** with small functions (no module-level constants
   beyond what is needed): `dev_root()` (from `PYPEIT_DEV`, raising a clear error if
   unset), `raw_data_dir()`, `tests_dir()`, `input_file_dir(kind)`, `priority_list_file()`,
   `tool_path(name)`.
3. **Replace hard-coded paths** with calls to `paths`:
   - `tests/framework/runner.py`: `RAW_DATA` (was test_main.py:654), `unit_tests`
     and `vet_tests` (was :910, :1007), `test_priority_list` (was :918).
   - `tests/framework/pypeit_tests.py`: `sensfunc_files`, `fluxing_files`,
     `flexure_files`, `coadd1d_files`, `tellfit_files` (was :494–747), `pypeit_files`
     and `coadd2d_files` (was :902, :910), `build_ql_calibs` (was :870).
   - `tests/unit/*.py` and `tests/vet/*.py`: ~150 `os.path.join(os.getenv('PYPEIT_DEV'),
     'RAW_DATA', ...)` / `... 'pypeit_files' ...` patterns (test_load_images, test_setups,
     test_frametype, test_metadata, test_scienceimage, test_runcalibonly, test_setup_gui,
     vet test_wavelengths, test_slitmask).  Convert them to the `paths` helpers or to the
     `tests/conftest.py` fixtures.  Do it as a mechanical substitution, in pathlib form
     (per AGENTS.md), and touch only the path expressions.
   - `tests/vet/test_ncpu.py:94` (`parent.parent / 'test_scripts'`).
   - `tests/wavecalib/test_*.py`: `TEST_DATA` path, per "Needed decisions" item 1.
   - Imports: `from test_scripts...` → `from tests.framework...`
     (`pypeit_test`, runner, registry, the self-tests).
4. **`pypeit_test`**: change only the import to `from tests.framework.runner import main`.
5. **Nautilus** (`tests/nautilus/gen_kube_devsuite`): update the default YAML path
   (was `$PYPEIT_DEV/nautilus/...`, :105).  `./pypeit_test` and the `RAW_DATA/` rclone
   target are unchanged.  Update `nautilus/rclone.conf` references.
6. **`.gitignore`**: re-point `pypeitdev/...`, `nautilus/...`, `unit_tests/*.sorted`
   patterns to their new locations.
7. **`REDUX_OUT_C/`**: per "Needed decisions" item 2.
8. **Documentation in this repo**: `README.rst` (layout, `pytest tests/unit`,
   `pytest tests/vet --redux_out ...`, `source tests/tools/source_headless_test.sh`,
   `cd $PYPEIT_DEV/tests/nautilus`; regenerate the `pypeit_test -h` listing),
   `AGENTS.md` "Repository layout"/"Adding to the dev suite", and the skills that name old
   paths (`.claude/skills/add-devsuite-setup`, `run-dev-suite`, `new-spectrograph`, ...;
   8 references).
9. **Coordinated main-repo PR (PypeIt)** — `RAW_DATA` stays put, so only these change:
   - `doc/scripts/make_example_files.py`: `DEV_ROOT / 'pypeit_files'` →
     `DEV_ROOT / 'tests' / 'input_files' / 'pypeit'`.
   - `doc/dev/*config.rst`, `*frames.rst`, `slitmask_ids.rst`: `${PYPEIT_DEV}/unit_tests`
     → `${PYPEIT_DEV}/tests/unit`.
   - `doc/calibrations/construct_template.rst`: stale `dev_algorithms/wavelengths/...`
     → `dev/pypeitdev/wavelengths/template_files/`.
   - `doc/include/dev_suite_readme.rst` (copy of README.rst) and a release-notes entry.

Execute it as one dev-suite branch (moves in one commit, path fixes in a second, so
`git log --follow` stays clean) plus a matching main-repo branch of the same name, so
Nautilus `-p/-d` branch pairing works.

## Verification

- `grep -rnE "test_scripts|unit_tests|vet_tests|_files['/]|PYPEIT_DEV'\), *'(pypeit|coadd|flux|sens|tell|flex)"`
  over the dev suite returns nothing outside history/plans.
- `./pypeit_test -h` and `./pypeit_test list` run from the top level; the `-h` output
  matches README.rst.
- `pytest tests/framework` (harness self-tests).
- `./pypeit_test unit` and `pytest tests/unit` (needs RAW_DATA).
- `./pypeit_test reduce afterburn ql -i shane_kast_blue -s 600_4310_d55` (and one
  setup that exercises each afterburn input dir: flux/sensfunc/coadd1d/coadd2d/tellfit/
  flexure).  The vet tests are not run here: they need output from every spectrograph,
  so at most confirm `pytest tests/vet --collect-only` imports cleanly.
- `cd tests/wavecalib && pytest` **and** `pytest tests/wavecalib` from the root (cwd
  independence).
- In the main repo: `cd doc ; make html` with `PYPEIT_DEV` set (example files are built
  from the moved `.pypeit` files).
- `gen_kube_devsuite --help` / dry-generate a YAML and inspect the command string.
