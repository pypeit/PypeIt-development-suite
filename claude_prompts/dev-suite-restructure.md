# Separating and refactoring the PypeIt test suite

This is the plan for moving the PypeIt test suite out of `PypeIt-development-suite` into
a new repository, `PypeIt-test`, and refactoring it so that each instrument/setup is
tested on its own.  It supersedes the earlier version of this file, which covered only
the restructure, and the separate `dev-suite-test-refactor.md`.  Decisions made so far,
with the alternatives considered, are logged in section 8.

## 1. Context

Two problems motivate the work:

1. **The dev suite mixes two things at its top level.**
   - The regression/integration test harness and its inputs: `test_scripts/`,
     `unit_tests/`, `vet_tests/`, seven `*_files/` input directories,
     `test_priority_list`, `build_ql_calibs` and `nautilus/`.
   - Development-only material: `pypeitdev/` (1554 tracked files of notebooks and
     analysis scripts), `sensfunc_archive/`, `claude_prompts/`, `reports/`,
     `deprecated/`, and so on.

   The testing code moves to its own repository.
2. **The test suite can't run one instrument on its own.**
   - **The vet tests can't be run for one instrument.**  Each vet test reaches for
     specific products of specific instruments through hard-coded paths, and nothing
     declares that mapping.  So the only reliable way to run the vet tests is to run the
     whole dev suite first, which takes many hours.
   - **Cloud runs download everything.**  The Nautilus job copies all of `RAW_DATA`
     (90 GB) and `CALIBS` before every run, even for a single setup.  Running
     `shane_kast_blue` alone needs 40 MB.

The goal is to make every instrument/setup ("dataset") a self-contained testing path:
raw data → reduction → selected afterburn scripts → QL → that dataset's own unit and vet
tests, with the cloud downloading only what the requested datasets need.
`./pypeit_test reduce` (and `all`, etc.) with no instrument still runs every dataset,
but as a set of independent dataset paths.

## 2. How each instrument/setup is tested today

### The registry
- `setups.py:all_setups` lists 69 instruments and 247 setups as plain strings.
- `test_setups.py` holds 19 per-test-type dicts of the form
  `{instr: {setup: [kwargs, ...]}}`.  Each kwargs dict becomes one test instance, e.g.
  `_sensfunc['gemini_gmos']['GS_HAM_R400_860'] = [dict(std_file=..., sens_file=...)]`.
- `all_tests` orders them into phases (PREP → REDUCE → AFTERBURN → QL) and pairs each
  dict with a `PypeItTest` factory.
- `_raw_data_dirs` remaps the four `p200_ngps_{u,g,r,i}` instruments onto one shared
  `RAW_DATA/p200_ngps`.

### The runner (`test_main.py`)
- It builds one `TestSetup` per dataset containing its ordered `PypeItTest` list.
- `-t N` workers each pull a whole dataset from a priority queue (ordered by
  `test_priority_list`) and run its tests serially.  The rest of a dataset is skipped
  after its first failure.
- `-i`/`-s` filter datasets with odd rules: a bare setup name that matches nothing in
  an instrument runs *all* of that instrument's setups, and an unknown `instr/setup` is
  ignored silently.
- `PypeItTest` subclasses (`pypeit_tests.py`) wrap one script each.  The scripts are
  `pypeit_setup`, `run_pypeit` (plus its calib-only, std, ncpu and step-by-step
  variants), `pypeit_sensfunc`, `flux_setup`, `flux_calib`, `multislit_flexure`,
  `collate_1d`, `coadd_1dspec`, `coadd_2dspec`, the datacube setup/coadd/extract
  scripts, `tellfit` and `pypeit_ql`.  Each runs as a subprocess in
  `REDUX_OUT/<instr>/<setup>/` and passes if the return code is 0.

### What is already per-dataset
No reduce, afterburn or QL step reads another dataset's output.  The three couplings are:
- the NGPS shared raw directory;
- the global `REDUX_OUT/QL_CALIB` built by `build_ql_calibs` from
  `keck_nires/ABpat_wstandard` and `keck_mosfire/Y_long`;
- the unused `vlt_xshooter_nir-vis_*.coadd1d` files, which point at `../` of other
  setups.

### What is global
`pypeit_tests`, `unit` and `vet` each run as a single pytest call over a whole directory.
None of them is filtered by `-i`/`-s`, their results are scraped from the console
output, and their failures don't change the exit code.

### How the vet and unit tests depend on data

| Category | Examples |
|---|---|
| Single dataset (most vet tests) | collate1d, edgetrace, findobj, fluxspec, quicklook, scattlight, step_by_step, wavetilts; `test_ncpu` and `test_datacube_afterburn` are already parametrized per dataset |
| Several setups of one instrument | `test_echelles` (loops over `all_setups[instr]`), most of `test_wavelengths`, hamspec extraction, binospec joint fit, xshooter/hires flux+coadd1d |
| Several instruments | coadd2d offsets (nires+mosfire), spat_flexure_science (5 instruments), assign_maskinfo (deimos+mosfire), `test_scripts` collate_1d and setup_coadd2d, the specutils loaders (kast_blue+gnirs), kcwi+kcrm coordinates; unit: maskdesign_matching, processimage, chk_calibs, setup_gui |
| Generic sweeps over many instruments | unit `test_load_images`, `test_setups`, `test_spectrographs` (42 dataset/file-count pairs), `test_construct_basename` |
| No raw data or reduction products | unit `test_mosaic::transform`, `test_opticalmodel` (2 tests) and `test_telluric`, which need no dev-suite files at all; the `test_input_files` readers, which parse every dev-suite input file; vet `ncpu_registration`, which checks the registry against `test_ncpu.py` |

Other findings:
- **Unregistered raw data.**  The unit tests read about 15 raw setups that are not
  registered datasets (e.g. `keck_deimos/830G_L_8400`, `1200G_Cooper`, `wht_isis_blue`).
- **Mini-pipeline tests.**  About 10 vet tests re-run PypeIt scripts and write into
  `REDUX_OUT` or the working directory: datacube, flux_coadd1d, sensfunc_from_onespec,
  skysub, redoslits (temporarily corrupts a Slits file), run_to_calibstep, and others.
- **Fragile references.**  Tests pick files by index into a sorted glob (`spec1d[1]`),
  use full file names containing timestamps, and parse log text.  The specutils loaders
  depend on whether the flux step already ran.

### Bugs found along the way
These are live, and are fixed in the dev suite's `develop` independently of this work
(see section 6, Phase 0):
- **Duplicate dict keys** in `_coadd1d['gemini_gmos']` and `_coadd2d['keck_nires']` mean
  the GS_HAM_R400_860 coadd1d and the ABBA_wstandard coadd2d never run.
- **`test_datacube` asserts can never fail.**  They are written as
  `assert(cond, msg)`, which asserts a non-empty tuple.
- **`test_frametype::test_hires` typo.**  It has `'lJ0306…'`, so it checks only one
  setup.
- **`PypeItSensFuncTest.run` raises a NameError** (it uses an undefined
  `search_pattern`) when the standard file is missing.
- **pytest failures don't set the exit code.**
- **`build_ql_calibs` has broken argument handling** in its fallback `run_devsuite`
  path.
- **`wavecalib/` is not runnable**: `test_longslit` takes arguments that pytest can't
  supply, and `test_ThAr` contains no test function.  It is not carried into
  `PypeIt-test` (decision 8), so this needs fixing only if the dev suite wants to keep
  running it.

## 3. Where and how the work is done

The work happens in a new repository, `PypeIt-test`, seeded from the dev suite with its
history intact (decision 6).  It is kept locally, plus a temporary public copy for cloud
runs, until it is ready to publish.

### The three branches involved
- **`PypeIt-development-suite:develop`** is the connection point (decision 13).  The
  team keeps developing there throughout, and `PypeIt-test` regularly merges from it.
- **`PypeIt-development-suite:dev-suite-restructure`** holds only this plan
  (decision 12).
  - It is edited as needed.  Nothing else depends on its commits, so it can be updated
    from `develop` in any way.
  - Once `agent_file`, which it branches from, is merged into `develop` (Phase 0), it
    only adds files under `claude_prompts/`, so it doesn't conflict with `develop`.
  - It is merged into `develop` at cut-over (Phase 10), just before the dev suite's
    testing code is removed.  Merging it early by mistake would only add this document.
- **`PypeIt-test:develop`** is where the restructure and refactor are built.

The dev suite and `PypeIt-test` share one VS Code workspace.  So this plan and the agent
skills, both of which stay in the dev suite (decisions 11 and 12), are reached from a
single place.

### Seeding `PypeIt-test`

Seed it by cloning the dev suite's `develop`, not by rewriting it:
```
git clone --branch develop https://github.com/pypeit/PypeIt-development-suite.git PypeIt-test
cd PypeIt-test
git remote rename origin devsuite
```
`PypeIt-test` then shares its commits with the dev suite, which is what makes the sync
below possible.  Rewriting the history (e.g. with `git filter-repo`) would give every
commit a new ID and break that connection.  The restructure (section 4) is then done
directly in `PypeIt-test`.

**The temporary public repository** (decision 10):
- **Why.**  Nautilus jobs clone from GitHub, so cloud runs of `PypeIt-test` need it
  hosted somewhere they can reach.  `PypeIt-test` is therefore pushed to a public
  repository in a personal GitHub account, e.g. `<user>/PypeIt-test`.
  `gen_kube_devsuite` already takes the repository URL (`--dev_fork`) and branch
  (`-d/--dev_branch`).
- **How it is used.**  It accepts no pull requests; work is merged locally and pushed.
- **What happens to it.**  It is deleted at cut-over, after the cleaned history is
  pushed to the `pypeit` organization.

Its history is the dev suite's history, which is already public, so nothing new is
exposed.  Forks of it would survive its deletion, which is one more reason not to invite
contributions to it.

**History.**
- Input files keep their history through both moves: into `tests/input_files/` and
  then into per-dataset directories.
- Most test files are split across many datasets, so their history stays with the
  original topic files rather than following the new per-dataset tests.
- The shared history makes `PypeIt-test` as large as the dev suite (1.7 GB packed),
  including the history of the deleted development material.  It is cleaned once, at
  cut-over, which brings it down to about 2 MB (see "Cleaning the history").  Until
  then, the Nautilus jobs use shallow clones, so the size does not affect them.

### Keeping up with development in the dev suite

Other developers keep adding datasets, input files and tests to the dev suite's
`develop` while the work is under way.  Because the two repositories share history,
those changes are **merged** into `PypeIt-test` rather than copied by hand:

```
git fetch devsuite
git merge devsuite/develop
```

Git follows renames when merging, so a change made in the dev suite to a file at its
old path is applied to the file at its new path.  What each kind of upstream change
turns into:

| Upstream change in the dev suite | Result of the merge in `PypeIt-test` |
|---|---|
| Edit to an existing `.pypeit` or afterburn input file | Applied automatically to the moved file, including in a per-dataset directory |
| Edit to a harness file, or to a unit or vet test not yet migrated | Applied automatically at its new path under `tests/` |
| New entry in `setups.py` or `test_setups.py` | Applied automatically to `tests/framework/setups.py`/`registry.py`.  Re-running the converter regenerates the affected `dataset.toml` files, and the equivalence check flags anything it misses. |
| New input file (e.g. a new `pypeit_files/<instr>_<setup>.pypeit`) | Lands at its old path, which no longer exists in `PypeIt-test`.  Move it with the converter. |
| Edit to a unit or vet test file that has already been split into per-dataset tests | Merge conflict (modified upstream, deleted here).  Port the change by hand into the same-named topic file in the relevant dataset's `unit/` or `vet/` directory, then resolve by keeping the deletion. |
| Edit to development-only material, `wavecalib/` or the agent skills (`pypeitdev/`, `.claude/`, …) | Merge conflict (modified upstream, deleted here).  Resolve by keeping the deletion. |
| New file in development-only material | Lands in a directory `PypeIt-test` deleted.  Delete it again. |
| Edit to `README.rst`, `AGENTS.md` or `.gitignore` | Content conflict, because `PypeIt-test` rewrote these files.  Keep `PypeIt-test`'s version, and port anything about testing by hand. |

So most changes come across automatically.  The rest show up as merge conflicts or as
files at old paths, which is a reliable list of what needs porting by hand, rather than
something a person has to notice in the dev suite's log.  The plans branch never enters
this flow, so edits to this document never cause conflicts.

To keep this working:
- **Merge often**, e.g. weekly and before starting each phase, so each merge stays small.
- **Move files and edit them in separate commits.**  Git matches renamed files by
  content, so a file that is moved and heavily edited in the same commit may not be
  recognized as moved.
- **Keep the old registry and the converter until cut-over.**  They are what let
  upstream registry changes be merged and then converted, rather than re-typed.
- **Never rewrite `PypeIt-test` history** (no `filter-repo`, no rebase of `develop`)
  while syncing.  Squash-merging `PypeIt-test`'s own feature branches is fine.
- **A small `sync_devsuite` script** in `tests/tools/` does the routine part of each
  merge:
  - run the merge;
  - resolve every conflict on development-only paths by keeping the deletion, and delete
    new files there;
  - run the converter;
  - list the remaining conflicts and any files left at old paths.

### Cleaning the history

**Why bother.**  `PypeIt-test` inherits the dev suite's whole history.  Almost all of it
is material the test repository doesn't need:
- raw data committed in 2017–2018 under `RAW_DATA/` and `RAW_DATA_old/`, about 1.2 GB;
- the development material in `dev_algorithms/`, later renamed `pypeitdev/`, about
  0.5 GB;
- old pixel-flat calibrations under `CALIBS/` and `PYPIT_CALIBS/`.

**Measured result.**  The cleanup was tested on a scratch mirror of the dev suite on
2026-11-19, using the keep-list below.
- The packed repository went from **1.71 GiB to 1.88 MiB**.
- 2,167 of 3,712 commits remain.  Commits that touched only removed paths are dropped.
- The test files keep their full history.  For example, `test_scripts/test_main.py`
  still has 172 commits back to 2018-08-15, and
  `pypeit_files/keck_deimos_830g_m_8500.pypeit` still has 17.
- An earlier run that also kept `wavecalib/` gave 23.8 MiB.  Almost all of the
  difference is `wavecalib/TEST_DATA/`, which is not carried over (decision 8).

**When.**  At cut-over, after the final merge from the dev suite and before the first
push to the `pypeit` organization.
- Rewriting the history gives every commit a new ID, which ends any further merging from
  the dev suite.
- Done before publishing, nobody else is working in `PypeIt-test`, so nothing breaks.
- Done after the team has started working in `PypeIt-test`, it would force everyone to
  re-clone and rebuild their open branches.  It might not even shrink the hosted repo,
  because GitHub keeps old commits reachable through pull-request refs.

**How.**
1. **Run it in a fresh clone**, as `filter-repo` requires.  Install the tool with
   `pip install git-filter-repo`.
2. **Measure first** with `git filter-repo --analyze`.  It only writes reports, to
   `.git/filter-repo/analysis/`, and changes nothing.
3. **Filter with a keep-list**, not a drop-list:
   `git filter-repo --paths-from-file keep-paths.txt`.  `filter-repo` matches paths as
   they were in each commit, so the list must name *every* path the test material has
   had.  Otherwise its history before a move is lost.  The keep-list used in the test:
   ```
   tests/
   pypeit_files/        pypit_files/
   coadd1d_files/       coadd2d_files/      fluxing_files/     sensfunc_files/
   tellfit_files/       flexure_files/
   test_scripts/        unit_tests/         vet_tests/
   nautilus/
   pypeit_test          pypit_test          pypeit_syncraw     source_headless_test.sh
   build_ql_calibs      build_ql_masters    build_nires_masters
   test_priority_list   README.rst          README.md          AGENTS.md
   CLAUDE.md            CHANGES.rst         .gitignore         pyproject.toml
   ```
   (one path per line in the actual file).
   - `pypit_files/`, `pypit_test`, `build_ql_masters`, `build_nires_masters` and
     `README.md` are earlier names of current files.
   - **`pyproject.toml` is kept for its new history.**  `PypeIt-test` adds a new
     top-level `pyproject.toml` for the pytest configuration (section 5.1).  Keeping the
     path also keeps the old `pypeitdev` packaging file's history at that path, which
     is small and harmless.
   - **Extend the list** with any other top-level paths added before cut-over.
4. **Repack and check.**  Run `git reflog expire --expire=now --all` and
   `git gc --prune=now --aggressive`.  Then compare `git log --follow` counts for a few
   files against the unfiltered repo, and run `./pypeit_test` from the cleaned clone.
5. **Push to a new, empty repository** in the `pypeit` organization (Phase 10).

**Side effects to expect:**
- `filter-repo` removes the `origin` remote.
- It updates commit IDs mentioned in commit messages, and writes the old-to-new mapping
  to `.git/filter-repo/commit-map`.
- Signed commits lose their signatures.
- Commit IDs quoted in issues, pull requests or documents elsewhere no longer resolve in
  `PypeIt-test`; they still resolve in the dev suite, which keeps the full history.

**Without rewriting.**  If the team decides against cleaning:
- day-to-day clones can stay small with `git clone --filter=blob:none`, which downloads
  file contents only as they are needed;
- the Nautilus jobs already use shallow clones;
- the hosted repository keeps the full 1.7 GB.

## 4. The restructure (done in Phase 2)

Immediately after seeding, `PypeIt-test` is restructured so that all testing code and
inputs sit under one top-level `tests/` directory (decision 14).  `./pypeit_test` stays
at the top level.  The development-only material is deleted rather than moved
(decision 7).

### Layout after the restructure

```
PypeIt-test/
├── pypeit_test                  # unchanged location; imports tests.framework.runner
├── pypeit_syncraw               # unchanged location; syncs RAW_DATA
├── README.rst  AGENTS.md  CLAUDE.md  CHANGES.rst  .gitignore   # rewritten for PypeIt-test
├── RAW_DATA -> ...              # git-ignored, unchanged (decision 15)
├── Telluric -> ...  REDUX_OUT/  QL_CALIB/      # git-ignored, unchanged
│
└── tests/                       # EVERYTHING the test suite runs or reads
    ├── __init__.py
    ├── conftest.py              # shared fixtures: test_root, raw_data_dir, input_file_dir
    ├── framework/               # was test_scripts/
    │   ├── __init__.py
    │   ├── paths.py             # NEW: the one place test-suite paths are defined
    │   ├── runner.py            # was test_main.py
    │   ├── setups.py            # unchanged: all_setups (instrument → RAW_DATA setup dirs)
    │   ├── registry.py          # was test_setups.py (imports setups; all_tests, ...)
    │   ├── pypeit_tests.py      # unchanged name (PypeItTest subclasses)
    │   ├── test_runner.py       # was test_pypeit_test.py (self-tests of the harness)
    │   └── test_priority_list   # was top-level; data file read/written by the runner
    ├── input_files/             # was the seven top-level *_files/ dirs
    │   ├── pypeit/              # was pypeit_files/
    │   ├── coadd1d/  coadd2d/  flux/  sensfunc/  tellfit/  flexure/
    ├── unit/                    # was unit_tests/  (+ __init__.py)
    ├── vet/                     # was vet_tests/   (+ __init__.py; conftest.py, files/ kept)
    ├── tools/
    │   ├── build_ql_calibs      # invoked by the QL tests
    │   ├── source_headless_test.sh
    │   └── sync_devsuite        # NEW (section 3)
    └── nautilus/                # was nautilus/ (remote execution of the suite)
```

Deleted, not carried over: `pypeitdev/`, `sensfunc_archive/`, `claude_prompts/`,
`reports/`, `deprecated/`, `wavecalib/` (decision 8), `fix_xshooter_names.py`,
`REDUX_OUT_C/`, the `pyproject.toml`, `setup.cfg` and `setup.py` that package
`pypeitdev`, and `.claude/` (decision 11).

Rationale for the less obvious choices:
- **Rename `test_main.py`/`test_setups.py`**: under `tests/`, pytest would try to collect
  every `test_*.py`.  Only real test modules keep the `test_` prefix.
- **`__init__.py` in `tests/`, `framework/`, `unit/`, `vet/`**:
  - `unit_tests` and `vet_tests` both contain `test_scripts.py` and `test_slitmask.py`,
    which collide if `pytest tests` is run as one session.  Packages give them distinct
    module names.
  - A regular `tests` package also cannot be shadowed by a stray `tests` package
    installed in site-packages; a namespace package can be.
  - With pytest's default `prepend` import mode, the repo root goes on `sys.path`, so
    tests can `from tests.framework import paths`.
  - Phase 5 later switches to pytest's `importlib` import mode (section 5.1).  These
    files are harmless under that mode.
- **`input_files/<kind>/`**: drops the redundant `_files` suffix and groups all inputs.
  File names (`<instrument>_<setup>.<ext>`) are unchanged.  Phase 3 later moves these
  files into per-dataset directories.
- **`test_priority_list` into `framework/`**: it is harness state, not a user-facing tool.

### Commit sequence

Each numbered item is its own commit.  The moves and the edits stay separate, so that
`git log --follow` and the merge from the dev suite both recognize the moved files.

1. **Delete the development-only material** listed above, and remove its patterns
   (`pypeitdev/...`, `dev_algorithms/...`) from `.gitignore`.
2. **Move the testing files with `git mv`** to the layout above.
3. **Fix the paths:**
   - Add `tests/framework/paths.py` with small functions (no module-level constants
     beyond what is needed):
     - `test_root()`, from the environment variable, raising a clear error if it is
       unset;
     - `raw_data_dir()`, `tests_dir()`, `input_file_dir(kind)`, `priority_list_file()`,
       `tool_path(name)`.
   - Replace the hard-coded paths with calls to `paths`:
     - `tests/framework/runner.py`: `RAW_DATA` (was test_main.py:654), `unit_tests`
       and `vet_tests` (was :910, :1007), `test_priority_list` (was :918).
     - `tests/framework/pypeit_tests.py`: `sensfunc_files`, `fluxing_files`,
       `flexure_files`, `coadd1d_files`, `tellfit_files` (was :494–747), `pypeit_files`
       and `coadd2d_files` (was :902, :910), `build_ql_calibs` (was :870).
     - `tests/unit/*.py` and `tests/vet/*.py`: about 150 `os.path.join(...)` patterns
       such as `os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', ...)` and
       `... 'pypeit_files' ...`.  They are in test_load_images, test_setups,
       test_frametype, test_metadata, test_scienceimage, test_runcalibonly,
       test_setup_gui, and vet test_wavelengths and test_slitmask.  Convert them to the
       `paths` helpers or the `tests/conftest.py` fixtures, as a mechanical substitution
       in pathlib form (per AGENTS.md), touching only the path expressions.
     - `tests/vet/test_ncpu.py:94` (`parent.parent / 'test_scripts'`).
     - Imports: `from test_scripts...` → `from tests.framework...`
       (`pypeit_test`, runner, registry, the self-tests).
   - `pypeit_test`: change the import to `from tests.framework.runner import main`.
   - `tests/nautilus/gen_kube_devsuite`: update the default YAML path (was
     `$PYPEIT_DEV/nautilus/...`, :105), the `nautilus/rclone.conf` reference, and the
     `source source_headless_test.sh` in the job's command string (:165), which becomes
     `tests/tools/source_headless_test.sh`.
   - `.gitignore`: re-point the test-related patterns (`nautilus/...`,
     `unit_tests/*.sorted`) to their new locations.
4. **Switch from `PYPEIT_DEV` to `PYPEIT_TST`** (decision 4) and rewrite the
   documentation.
   - `PYPEIT_TST` points to the top level of `PypeIt-test`, where `RAW_DATA` is placed or
     symlinked.  With `paths.py` in place, the switch touches:
     - `paths.test_root()` and the environment check in `pypeit_test`;
     - in `tests/nautilus/`: the `env` block of `kube_dev_suite.yaml`, and
       `gen_kube_devsuite`'s own use of the variable (:105).  Also rename the directory
       the job clones into, which is hard-coded as `PypeIt-development-suite`
       (:163–164), to `PypeIt-test`, so that it matches the `PYPEIT_TST` value in the
       YAML.
   - Rewrite `README.rst` for `PypeIt-test`, and check that the `pypeit_test -h` listing
     in it matches the current output.
   - Rewrite `AGENTS.md` for `PypeIt-test`.  It reaches `develop` with `agent_file` in
     Phase 0, so the seeded clone has it.  `CLAUDE.md` refers to it.
   - Add a note to the top of `README.rst`: the repository is pre-release, the dev suite
     remains authoritative until cut-over, and pull requests are not accepted.

### Verification

- **No old paths remain.**
  `grep -rnE "test_scripts|unit_tests|vet_tests|(pypeit|coadd1d|coadd2d|fluxing|sensfunc|tellfit|flexure)_files|PYPEIT_DEV|PypeIt-development-suite" tests pypeit_test`
  returns nothing.  (The pattern names the seven old directories explicitly, because a
  bare `_files/` would also match the new `tests/input_files/`.)
- **The runner works from the top level**, with only `PYPEIT_TST` set:
  `./pypeit_test -h` and `./pypeit_test list` both run.
- **Harness self-tests pass**: `pytest tests/framework`.
- **Unit tests run** (they need RAW_DATA): `./pypeit_test unit` and `pytest tests/unit`.
- **Reductions and afterburn steps run.**
  - Run `./pypeit_test reduce afterburn ql -i shane_kast_blue -s 600_4310_d55`.
  - Also run one setup that exercises each afterburn input directory: flux, sensfunc,
    coadd1d, coadd2d, tellfit and flexure.
- **Vet tests import cleanly.**  They are not run here, because they need output from
  every spectrograph.  Confirm `pytest tests/vet --collect-only` succeeds.
- **The Nautilus job runs.**  Dry-generate a YAML with `gen_kube_devsuite`, inspect the
  command string, and run one setup with `--dev_fork` pointing at the temporary public
  repository.

## 5. Target design

### 5.1 The testing path for one dataset

**Core idea: everything about a dataset lives in one directory.**  That covers its
configuration, its input files, and its own unit and vet tests.  The dataset expands into
a `TestPlan`, an ordered list of `Step`s that ends by running the tests in its
directory.  Every test uses data from that one dataset only.

```
Dataset            instrument, setup, raw_dir, workdir (REDUX_OUT/<instr>/<setup>), config
 └─ TestPlan       the Steps for the requested stages, in dependency order
     ├─ UNIT       pytest <dataset dir>/unit/            (raw data only)
     ├─ PREP       setup                         (pypeit_setup, or copy + fix .pypeit)
     ├─ REDUCE     reduce, reduce_std, calib_only, ncpu, step_by_step, calib_step_by_step
     ├─ AFTERBURN  sensfunc, flux_setup, flux, flexure, collate1d, coadd1d, coadd2d,
     │             datacube_setup, datacube_coadd, datacube_extract, tellfit
     ├─ QL         ql (one Step per configured QL run)
     └─ VET        pytest <dataset dir>/vet/             (this dataset's products)
```

**`Step`** (`tests/framework/steps.py`): one class per script.  The current
`PypeItTest` subclasses port almost 1:1.  Each Step declares:
- `kind` and `stage`;
- `requires`: the step kinds it consumes.  For example, flux requires sensfunc and
  tellfit requires coadd1d.  This turns today's implicit order into explicit
  dependencies.
- `inputs()`: files it needs.  Used both for the preflight check and for the data
  manifest in section 5.4.
- `command()` and `outputs()`: the expected products, checked after the run, so "exit
  code 0 but no spec1d" fails at the step instead of later in a vet test.

After running, a Step records a `StepResult` (status, timing, memory, log path) in
`<workdir>/test_state.json`.  That file is how later stages, and later invocations, know
what has been done.

#### Per-dataset tests and shared helpers

Each dataset directory holds its own tests, in two subdirectories, and the tests call
shared helper functions (decision 5).

- **`unit/`** holds tests that need only the dataset's raw data.
- **`vet/`** holds tests that check the dataset's reduction products.
- **Each subdirectory holds several test files, one per topic.**  The file names come
  from a fixed list that mirrors today's topic files, e.g.:
  - `unit/`: `test_load_images.py`, `test_frametype.py`, `test_setups.py`,
    `test_metadata.py`;
  - `vet/`: `test_wavelengths.py`, `test_flexure.py`, `test_coadd2d.py`,
    `test_quicklook.py`.

  A dataset has only the files for the topics it covers.  The fixed names keep the
  datasets consistent and searchable.  They also allow one topic to be run across every
  dataset by path, e.g. `pytest tests/datasets/*/*/vet/test_wavelengths.py`, provided
  those datasets' products exist.
- **The helpers live in `tests/helpers/`**, organized by topic: `images.py`,
  `wavecalib.py`, `flexure.py`, `coadd.py`, `datacube.py`, `slitmask.py`, and so on.
  They hold the logic that today is duplicated or generalized inside topic-organized test
  files.
- **Dataset-specific details stay in the per-dataset tests**: the file names, detectors,
  expected values and regression thresholds (decision 2).

For example, today's `unit_tests/test_load_images.py` sweeps about 20 instruments.  It
becomes a `grab_img` helper plus a short test in each dataset directory (the file name
and shape below are illustrative):

```python
# tests/helpers/images.py
def grab_img(ds, filename, det=1):
    """
    Read and process one raw frame of the dataset ``ds``.
    """
    ...

# tests/datasets/keck_deimos/830G_L_8400/unit/test_load_images.py
from tests.helpers.images import grab_img

def test_load_image(dataset):
    img = grab_img(dataset, 'DE.20100913.57161.fits', det=3)
    assert img.shape == (4096, 2128), 'DEIMOS detector 3 should have the full-frame shape'
```

The vet tests follow the same pattern.  For example, the RMS table in
`vet_tests/test_wavelengths.py` becomes a `check_wavecalib_rms(dataset, rms_max)` helper.
Each dataset's `vet/test_wavelengths.py` then calls it with its own threshold.

**Vet files declare the steps they depend on.**  Each vet file names the steps whose
products it checks, with a module-level marker:

```python
# tests/datasets/gemini_gmos/GS_HAM_R400_860/vet/test_fluxing.py
pytestmark = pytest.mark.requires_steps('sensfunc', 'flux')
```

- When a step a file depends on did not succeed in this run (according to
  `test_state.json`), the plugin skips that file and reports the reason.  The failure is
  then reported once, at the step, rather than again by every downstream check.
- With `--steps` (section 5.3), only the vet files that depend on the selected steps
  run.  For example, `--steps flux` runs `vet/test_fluxing.py` but not
  `vet/test_wavelengths.py`.
- A vet file without the marker depends on the reduction only.

This marker records which steps a file depends on.  It does not bind tests to datasets;
the directory does that.

**Each test uses data from one dataset only** (decision 3).
- **Cross-dataset tests are split.**  Today's tests that span several setups or
  instruments become one test per dataset.  For example, coadd2d offsets for nires+mosfire
  become a nires test and a mosfire test, and `test_echelles` becomes one test in each
  echelle dataset.
- **Comparisons between datasets are dropped or replaced** by checks against stored
  reference values.  The one example is `gd50_coordinates_match_across_datasets`
  (kcwi vs kcrm), which becomes a coordinate check in each dataset.
- **Tests that need no dev-suite files at all** (`test_mosaic::transform`,
  `test_opticalmodel`'s two data-free tests, `test_telluric`) move to the main repo's
  `pypeit/tests`.
- **The `test_input_files` readers become part of config validation.**  Loading a
  `dataset.toml` parses each of its input files with PypeIt's own input-file classes,
  so every input file is checked whenever the registry is loaded.
- **`ncpu_registration` is dropped.**  It checks that the registry and `test_ncpu.py`
  agree.  That link is now made by the `[[ncpu]]` table in `dataset.toml` and the
  `requires_steps('ncpu')` marker in the dataset's vet file.

**The `dataset` fixture** (`tests/conftest.py`) identifies the dataset by searching
upward from the requesting test file for the nearest `dataset.toml`.  No per-dataset
`conftest.py` is needed.
- **What it returns**: a `DatasetContext` with `workdir`, `raw_dir`, the parsed
  `config`, and a `scratch` directory.
- **`products_dir`**, where the reduction products are.  This is usually `workdir`.  For
  datasets whose PREP step runs `pypeit_setup` (e.g. `shane_kast_blue/600_4310_d55`), it
  is `workdir/<instrument>_A/`, where today's `PypeItSetupTest` sends the outputs.  The
  QL outputs stay directly under `workdir`.  Tests use `products_dir` rather than
  knowing this rule.
- **File helpers.**  `find_one(pattern)` and `find_all(pattern)` raise a clear error when
  the wrong number of files match, in place of `[0]`/`[1]` on a glob.
- **Paths stay inside the dataset.**  The helpers refuse paths outside the dataset's
  `raw_dir` and `workdir`, which enforces the one-dataset rule.
- **Mini-pipeline tests write to `scratch`**, never into the workdir.
- **`--redux_out`** remains a plugin option, and `workdir` is resolved from it.

**pytest configuration.**  Because of the fixed topic names, many datasets hold a test
file with the same name, e.g. `vet/test_wavelengths.py`.  Setup names such as `1.5_2x3`
and `J0100+2802_…` are not valid Python package names, so the test directories cannot
be made into packages.  Instead, a new top-level `pyproject.toml` sets:
- `--import-mode=importlib`, so test modules with the same basename can coexist;
- `pythonpath = ["."]`, so `tests.helpers` and `tests.framework` are importable.

`tests/unit/` and `tests/vet/` go away once their tests are migrated.

### 5.2 Configuring a dataset: one TOML file per dataset

Each dataset is configured by one `dataset.toml` (decision 1), kept in the dataset's
directory together with its input files and tests:

```
tests/datasets/<instrument>/<setup>/
    dataset.toml
    <instr>_<setup>.pypeit  [<instr>_<setup>_std.pypeit]
    <instr>_<setup>.sens  .flux  .coadd1d  .coadd2d  .tell  .flex     (as needed)
    unit/                                                              (if any)
        test_load_images.py  test_frametype.py  ...
    vet/                                                               (if any)
        test_wavelengths.py  test_fluxing.py  ...
```

This replaces the per-kind `tests/input_files/<kind>/` directories from the restructure,
and the Python registry dicts.

- **Templates for new datasets.**  `tests/framework/templates/dataset.toml` is a
  commented template listing every step kind and option.  The same directory has template
  `unit/` and `vet/` directories with a starting test file for each common topic.  To add a
  dataset, a developer copies `dataset.toml` and the test files for the topics they need
  into a new `tests/datasets/<instrument>/<setup>/`, and edits them.  The list of allowed
  topic file names is kept in the templates directory's README.
- **TOML is read with `tomllib`** (standard library for Python ≥3.11).  Duplicate keys
  are a parse error, so mistakes like the duplicate dict keys in section 2 can't recur.
- **Validation.**  A dataclass schema (`tests/framework/config.py`) rejects unknown keys
  and unknown step kinds, and checks that every referenced file exists.

Example, for a dataset with most afterburn steps:

```toml
# tests/datasets/gemini_gnirs_echelle/32_SB_SXD/dataset.toml
instrument = "gemini_gnirs_echelle"
setup      = "32_SB_SXD"
raw        = "gemini_gnirs_echelle/32_SB_SXD"   # relative to RAW_DATA; this is the default
tags       = ["echelle", "nir"]                 # optional, for -k style selection

[[reduce]]                                       # one table per reduce run
pypeit_file = "gemini_gnirs_echelle_32_sb_sxd.pypeit"

[[sensfunc]]
std_file  = "Science/spec1d_*HIP62745*.fits"
sens_file = "gemini_gnirs_echelle_32_sb_sxd.sens"

[[coadd1d]]
file = "gemini_gnirs_echelle_32_sb_sxd.coadd1d"

[[coadd2d]]
file = "gemini_gnirs_echelle_32_sb_sxd.coadd2d"

[[tellfit]]
coadd_file = "pisco_coadd.fits"
file       = "gemini_gnirs_echelle_32_sb_sxd.tell"

[[ql]]
name    = "cooked"
files   = ["N20160127S0001.fits", "N20160127S0002.fits"]
options = { "--setup_calib_dir" = "{calib_dir}" }   # {calib_dir} etc. replace the sentinels
```

Conventions:
- Each `[[kind]]` array of tables is one Step instance.  Several datasets have more than
  one: KCWI has two datacube targets and the QL tests have several runs.
- Stage order is fixed by the step kind.  An optional `after = [...]` adds
  dependencies within the dataset.
- **A dataset with no `[[reduce]]` is data-only.**  It registers raw data used only by
  its `unit/` tests (e.g. `keck_deimos/830G_L_8400`), so every raw directory any test
  touches is a registered dataset.
- **Shared raw data**: the four NGPS datasets each set
  `raw = "p200_ngps/1.5_2x3"`.
- **QL calibrations**: `build_ql_calibs` writes into the dataset's own workdir instead of
  the global `REDUX_OUT/QL_CALIB`.  The QL runs that use it are configured on the same
  datasets that produce the calibrations (`keck_nires/ABpat_wstandard`,
  `keck_mosfire/Y_long`), so no dataset depends on another.
- **Migration**: a converter generates every `dataset.toml` from `setups.py`,
  `test_setups.py` and the input-file directories, and moves the input files.  An
  equivalence check compares the step list from the old dicts with the step list from
  the new files, so nothing is dropped.
- **Which is authoritative until cut-over.**  While `PypeIt-test` still merges changes
  from the dev suite (section 3), the old dicts remain authoritative for every dataset
  they define.
  - The converter regenerates those datasets' `dataset.toml` files after each merge, so
    those files are not edited by hand.
  - A dataset added only in `PypeIt-test` has no dict entry, so its `dataset.toml` is
    written by hand and the converter leaves it alone.
  - At cut-over the dicts are deleted, and every `dataset.toml` becomes authoritative.

### 5.3 Choosing which steps to run for a dataset

These run on the data selected by `-i`/`-s`:
- The test-type positionals map to stages: `unit`, `reduce` (PREP+REDUCE),
  `afterburn`, `ql`, `vet`, and `all`.
- **`pypeit_tests` stays as the one global stage** (decision 9).  It runs the main PypeIt
  repo's own pytest suite (`pypeit/tests`) once, before the dataset jobs, when
  `pypeit_tests` or `all` is requested.  It belongs to no dataset, and `-i`/`-s` don't
  affect it.
- **New `--steps sensfunc,flux,coadd1d`** narrows the run to those step kinds and
  automatically adds what they require.  The vet stage then runs only the vet files
  whose `requires_steps` marker names one of the selected steps (section 5.1).
- **New `--reuse`** skips a required step whose `test_state.json` records a success
  with the same PypeIt version and unchanged inputs (compared by a hash of the input
  files).  This makes "re-run only flux and its vet tests for one setup" a minutes-long
  task.
- **New `--no-deps`** runs exactly the named steps on whatever output already exists.
- **Selection rules are made strict:**
  - `-s` with a bare setup name matches setups of that name across instruments.
  - A filter that matches nothing is an error.
  - `instr/setup` must exist.
- **New `./pypeit_test plan ... [--json]`** prints the expanded plans and the
  `RAW_DATA` subdirectories they need, without running anything.  Section 5.4 builds on
  this.

### 5.4 Running all requested tests

`./pypeit_test all` (or any subset of stages, instruments and setups):

1. **Load and validate** every `tests/datasets/*/*/dataset.toml`, in sorted order, and
   apply the `-i`/`-s` filters.
2. **Expand** each selected dataset into a TestPlan for the requested stages and steps.
   If `pypeit_tests` was requested, run it now, before the dataset jobs.
3. **Fetch and preflight.**  The data manifest is the union of the `raw` directories of
   the selected datasets.  Because every test uses only its own dataset, nothing else is
   needed.
   - Locally, check that they exist; missing data fails just that dataset, not the
     whole run.
   - In the cloud, download only those directories (see "Cloud" below).
4. **Schedule** each dataset's plan as one job on a pool of `-t N` workers, the slowest
   first.
   - Timings live in a `timings.json` updated after every successful dataset, not just
     after full passing runs.  It replaces `test_priority_list`.
   - Within a job, the steps run in dependency order, and the dataset's `unit/` and
     `vet/` tests run as their stages come up.
   - pytest runs as
     `pytest <dataset dir>/vet --redux_out … -p no:cacheprovider --junitxml=<workdir>/vet.xml`
     (and likewise for `unit/`).  Results are read from the JUnit XML rather than
     scraped from the console, so several datasets can run their tests in parallel.
5. **Report.**  Aggregate the per-dataset `test_state.json` and JUnit files into the
   text report and CSV.  The exit code counts step *and* pytest failures.

**Cloud** (`tests/nautilus/gen_kube_devsuite`):
- The job runs `./pypeit_test plan <args> --json` inside the container and copies only
  the listed directories with `rclone copy --files-from`.
- It drops the `CALIBS` copy (nothing uses it).
- It lets the PypeIt cache fetch telluric grids on demand instead of copying all seven.
- Because jobs are per dataset, a follow-up can shard a full run across several pods,
  each with its own subset of datasets and downloads, and merge the reports from S3.

## 6. Implementation plan

Each phase ends with a check that must pass before the next starts.  From Phase 2 on,
each phase is one or more branches in `PypeIt-test`, merged locally into its `develop`
and pushed to the temporary public repository.  Each keeps `./pypeit_test all` working.

### Phase 0: Commit the plan and merge `agent_file`
1. Commit this document to `dev-suite-restructure` under `claude_prompts/`, and delete
   `claude_prompts/dev-suite-test-refactor.md`.  Share it with the team.
2. **Merge `agent_file` into the dev suite's `develop`**, by the normal PR process.
   - `dev-suite-restructure` branches from `agent_file`.  That branch adds `AGENTS.md`
     and updates `CLAUDE.md`, `README.rst`, `.gitignore`, `source_headless_test.sh` and
     `test_scripts/pypeit_tests.py`.
   - **This step gates Phase 1.**  `PypeIt-test` is seeded from `develop`, so these
     changes, `AGENTS.md` in particular, must be there first.
   - Afterward, merge `develop` into `dev-suite-restructure`.  The branch then differs
     from `develop` only by this plan (section 3).  This works whether the PR is merged
     or squash-merged, because the content is the same.
3. **Check:** `git diff --stat origin/develop...dev-suite-restructure` lists only
   `claude_prompts/` files.

**Work that runs alongside, and does not gate the phases below:**
- **Live bug fixes.**  The bugs from section 2 are fixed on a separate branch from the
  dev suite's `develop`, and merged into `develop` by the normal PR process.  The next
  sync then brings them into `PypeIt-test`.  Re-enabling the two afterburn tests may
  expose failures; fix them, or record them as known.  The `wavecalib/` fix is optional,
  since `wavecalib/` is not carried into `PypeIt-test` (decision 8).
- **The `construct_template.rst` path** in the main PypeIt repo.  It points at
  `$PYPEIT_DEV/dev_algorithms/wavelengths/template_files/`, which no longer exists; in
  the current layout it should be `$PYPEIT_DEV/pypeitdev/wavelengths/template_files/`.
- **Regular dev-suite runs.**  The team keeps running the existing suite on the dev
  suite's `develop` throughout.  Those runs provide the reference results that the new
  engine is compared against (Phases 4 and 9).

### Phase 1: Seed `PypeIt-test`
1. Clone the dev suite's `develop` as `PypeIt-test` and rename the remote to `devsuite`
   (section 3).
2. Add `PypeIt-test` to the VS Code workspace, alongside the dev suite.
3. Create the temporary public repository in a personal GitHub account.  Add it as
   `PypeIt-test`'s `origin`, but push only after Phase 2.
4. **Check:** `git merge-base HEAD devsuite/develop` is `devsuite/develop` itself, i.e.
   the two share history.

### Phase 2: Restructure (section 4)
1. Make commits 1–4 of section 4's commit sequence.
2. Write `tests/tools/sync_devsuite` and do a first sync from `devsuite/develop`.
3. Push to the temporary public repository.
4. **Check:** section 4's verification list, including one Nautilus job against the
   temporary repository.

### Phase 3: Dataset configuration
1. Add `tests/framework/config.py`: the dataclass schema, the TOML loader and the
   validation.  Validation includes parsing every input file with PypeIt's input-file
   classes, which replaces the `test_input_files` readers.
2. Add the templates in `tests/framework/templates/`: `dataset.toml`, `unit/` and `vet/`
   starter files, and the README with the allowed topic file names.
3. Write the converter and the equivalence check.  Generate every `dataset.toml`.
4. Move the input files into `tests/datasets/<instr>/<setup>/` in their own commit.
   Then point the old runner's file lookups at the new locations.
5. Add `sync_devsuite`'s converter step.
6. **Check:**
   - the equivalence check passes for all 247 setups;
   - the old runner's `reduce afterburn ql` passes on a handful of setups that cover
     every input-file kind;
   - a test sync with a dummy upstream change to a `.pypeit` file applies it in the
     dataset directory.

### Phase 4: New engine
1. Add `tests/framework/steps.py`, porting each `PypeItTest` subclass to a `Step` with
   `requires`, `inputs()`, `command()` and `outputs()`.
2. Add `TestPlan`, the worker-pool scheduler, `test_state.json`, `timings.json`
   (seeded from `test_priority_list`), and JUnit-based reporting.  The exit code counts
   step and pytest failures.
3. Add the CLI behind `--engine new`:
   - the strict `-i`/`-s` rules;
   - `--steps`, `--reuse` and `--no-deps`;
   - `plan [--json]`;
   - the global `pypeit_tests` stage.
4. Move the QL calibrations into each dataset's workdir (`build_ql_calibs`).
5. Rewrite the harness self-tests in `tests/framework/` for the new engine.
6. **Check:**
   - `pytest tests/framework` passes;
   - on a full local run, `--engine new` reproduces the `reduce afterburn ql` results of
     `--engine old`, or of a recent dev-suite `develop` run with the same PypeIt
     version: the same steps pass, with products in the same places;
   - `plan --json` lists the right raw directories for a few selections.

### Phase 5: pytest infrastructure and a pilot dataset
1. Add a top-level `pyproject.toml` with the pytest configuration:
   `--import-mode=importlib` and `pythonpath = ["."]`.
2. Add `tests/conftest.py` with:
   - the `dataset` fixture and `DatasetContext` (`workdir`, `products_dir`, `raw_dir`,
     `config`, `scratch`, `find_one`/`find_all`, path checks);
   - `--redux_out`;
   - the `requires_steps` marker handling.
3. Create `tests/helpers/` with the first helper modules.
4. Pilot with `shane_kast_blue/600_4310_d55`, the dataset with the most vet coverage:
   - migrate all of its unit and vet tests into its `unit/` and `vet/` directories;
   - split the tests it shares with gnirs and deimos;
   - wire the UNIT and VET stages into the new engine.
5. Start a migration ledger, a table in `tests/MIGRATION.md`.  It maps every old test
   function to its new location(s), marked as dropped or moved where applicable.
6. **Check:**
   - `./pypeit_test all --engine new -i shane_kast_blue -s 600_4310_d55` runs this
     dataset's unit and vet tests, in minutes;
   - `--steps flux --reuse` runs only the flux step and its vet files;
   - the old global vet run still passes for the rest.

### Phase 6: Migrate the vet tests
1. One instrument, or a few, per branch.  For each test function:
   - move its logic into `tests/helpers/` if more than one dataset uses it;
   - write the per-dataset tests in the topic files, with `requires_steps` markers;
   - split cross-dataset tests;
   - replace fragile references with `find_one` and reference values;
   - move mini-pipeline tests to `scratch`;
   - remove the old test from `tests/vet/`;
   - update the ledger.
2. Sync from the dev suite before each branch.  Port any conflicting upstream vet-test
   changes into the new files.
3. **Check, per branch:** the migrated datasets pass under `--engine new`, and the
   ledger has no unmapped entries for them.  **At the end:** `tests/vet/` holds no
   tests.

### Phase 7: Migrate the unit tests
1. Same procedure, into `unit/`.  Add data-only datasets for the roughly 15 raw
   directories used only by unit tests.  Split the generic sweeps (`test_load_images`,
   `test_setups`, `test_spectrographs`, `test_construct_basename`) into per-dataset
   files.
2. Move the data-free tests (`test_mosaic::transform`, the two `test_opticalmodel`
   tests, `test_telluric`) to the main PypeIt repo, in a main-repo PR.
3. **Check:** `tests/unit/` holds no tests, the ledger is complete, and
   `./pypeit_test unit --engine new` passes.

### Phase 8: Cloud
1. Change `gen_kube_devsuite` to run `pypeit_test plan --json` in the container and copy
   only the listed raw directories (`rclone copy --files-from`).
   - Drop the `CALIBS` copy.
   - Let the PypeIt cache fetch telluric grids.
2. Optionally, add sharding: split the selected datasets across several pods and merge
   the reports from S3.
3. **Check, on Nautilus against the temporary repository:**
   - a `shane_kast_blue` job downloads about 40 MB and passes;
   - a full `all` job passes and matches the local new-engine results.

### Phase 9: New engine by default
1. Make `--engine new` the default, keeping `--engine old` for the remaining
   comparisons.
2. **Check:** full local and cloud `all` runs pass.  The ledger shows every old check
   accounted for.

### Phase 10: Cut-over
1. **Freeze test development in the dev suite.**  Announce a date after which new
   datasets and tests go to `PypeIt-test`.
2. **Do the final sync** from `devsuite/develop` and port what it flags.
3. **Remove the old harness from `PypeIt-test`**: the dicts, `setups.py`, `registry.py`,
   the old runner, the converter, `test_priority_list`, `sync_devsuite`, `--engine`, and
   the empty `tests/unit/` and `tests/vet/`.  Rewrite the README and AGENTS for the final
   layout.
4. **Clean the history** in a fresh clone (section 3, "Cleaning the history"), with the
   keep-list extended for any new top-level paths.  Check the `git log --follow` counts
   and run `./pypeit_test` on one dataset from the cleaned clone.  This permanently ends
   merging from the dev suite.
5. **Publish** by pushing the cleaned history to a new, empty `pypeit/PypeIt-test`.
   - Change the default `--dev_fork` in `gen_kube_devsuite` to the new repository.
   - Remove the pre-release note from the README.
   - Delete the temporary personal repository.
6. **Main PypeIt repo.**  Switch every reference to dev-suite data, input files or test
   paths to `PYPEIT_TST` and the per-dataset layout, and add a release-notes entry.  The
   references:
   - `doc/scripts/make_example_files.py`, whose `.pypeit` files are now at
     `tests/datasets/<instr>/<setup>/`;
   - the `doc/dev/*config.rst`, `*frames.rst` and `slitmask_ids.rst` pages, which point
     at `${PYPEIT_DEV}/unit_tests/...` and now need the matching dataset `unit/`
     directories;
   - `doc/include/dev_suite_readme.rst`, which becomes `PypeIt-test`'s README;
   - `pypeit/data/tests/copy_cooked_files.src`, which copies from
     `$PYPEIT_DEV/REDUX_OUT`.
7. **Dev suite.**
   1. Merge `dev-suite-restructure` into `develop`, which brings this plan into
      `claude_prompts/` as part of the historical record (decision 12).  Then delete the
      branch.
   2. Remove the testing code and inputs; nothing else in the dev suite is moved or
      reorganized.  Rewrite the testing sections of `README.rst` and `AGENTS.md` to point
      at `PypeIt-test`, and drop the test-related patterns (`nautilus/...`,
      `unit_tests/*.sorted`) from `.gitignore`.  The dev suite keeps `PYPEIT_DEV` for the
      development-only material.
      - **Removed**, because they now live in `PypeIt-test`: `test_scripts/`,
        `unit_tests/`, `vet_tests/`, the seven `*_files/` directories, `nautilus/`,
        `pypeit_test`, `pypeit_syncraw`, `build_ql_calibs`, `source_headless_test.sh`
        and `test_priority_list`.
      - **Top level afterward**, as tracked on `develop` (git-ignored local data such as
        `RAW_DATA` is not shown):
        ```
        PypeIt-development-suite/
        ├── .claude/                  # agent skills (decision 11)
        ├── .gitignore
        ├── AGENTS.md
        ├── CHANGES.rst
        ├── CLAUDE.md
        ├── README.rst                # points at PypeIt-test for testing
        ├── REDUX_OUT_C/
        ├── claude_prompts/           # including this plan, from step 7.1
        ├── deprecated/
        ├── fix_xshooter_names.py
        ├── pypeitdev/
        ├── pyproject.toml            # these three package pypeitdev
        ├── setup.cfg
        ├── setup.py
        ├── reports/
        ├── sensfunc_archive/
        └── wavecalib/                # stays (decision 8)
        ```
   3. Rewrite the skills that operate on the test suite (`add-devsuite-setup`,
      `run-dev-suite`, and the dev-suite parts of `new-spectrograph`) for `PypeIt-test`
      and the new engine.
8. **Check:**
   - a fresh clone of `pypeit/PypeIt-test` passes `./pypeit_test all` locally and on
     Nautilus;
   - the main repo's `cd doc ; make html` builds with `PYPEIT_TST` set;
   - the acceptance criteria in section 7 are all met.

## 7. Acceptance criteria

- **Coverage parity.**
  - The new engine runs the same step set as the old dicts, checked by the equivalence
    test.
  - Every check in the old `unit_tests/` and `vet_tests/` has a counterpart in some
    dataset's `unit/` or `vet/` directory, or has moved to the main PypeIt repo.
  - A full `./pypeit_test all` passes on the new engine.
- **Single-dataset turnaround.**
  `./pypeit_test reduce afterburn vet -i shane_kast_blue -s 600_4310_d55` finishes in
  minutes and runs exactly that dataset's tests.
- **Narrow re-runs.**  `./pypeit_test afterburn vet --steps flux --reuse -i gemini_gmos`
  re-runs only the flux step and the vet files that depend on it.
- **Small cloud downloads.**  A cloud job for `shane_kast_blue` downloads about 40 MB,
  not 90 GB.
- **Self-containment.**
  - Every unit and vet test lives in a dataset directory.  The harness's own self-tests
    in `tests/framework/` are the only other tests.
  - None of them reads outside its own dataset; the `DatasetContext` path checks
    enforce this.

## 8. Decision points

### 1. Configuration format
- **Question:** How should each dataset's test configuration be written?
- **Options considered:**
  - one TOML file per dataset;
  - Python modules (e.g. one `tests/datasets/<instr>.py` per instrument);
  - the current Python dicts in `test_setups.py`.
- **Decision:** one TOML file per dataset, plus a template that developers copy and edit
  to add a new dataset.
- **Rationale:** a per-dataset TOML file is easy to read and edit, and the template makes
  adding a dataset a copy-and-edit task.

### 2. Regression thresholds
- **Question:** Should the expected values in the vet tests (wavelength RMS tables,
  flexure offsets, slit counts) move into `dataset.toml` as data, or stay in the test
  code?
- **Options considered:**
  - move them into `dataset.toml`, so that one generic vet test could serve many
    datasets;
  - keep them in the vet tests.
- **Decision:** keep them in the vet tests for now, and revisit as development proceeds.
- **Rationale:** this follows from decision 5.  Each dataset has its own vet test files,
  which are the natural place for its expected values, and the shared logic lives in the
  helpers.

### 3. Cross-dataset tests
- **Question:** Should tests that use data from more than one instrument/setup be kept?
- **Options considered:**
  - keep them, in a cross-dataset stage that runs after all datasets finish and only when
    every dataset involved succeeded;
  - split them so that each test needs data from only one instrument/setup.
- **Decision:** split all tests so that each test needs data from only one
  instrument/setup.
- **Rationale:** with every test bound to a single dataset, each dataset is a
  self-contained testing path.  Its tests can run as soon as its reduction finishes, and a
  run (or cloud download) for one dataset never needs another.

### 4. Data environment variable
- **Question:** Once the tests leave the dev suite, should `PYPEIT_DEV` still mean "the
  directory that contains `RAW_DATA`", or should `PypeIt-test` get its own variable?
- **Options considered:**
  - keep `PYPEIT_DEV` for both repositories;
  - give `PypeIt-test` its own variable.
- **Decision:** migrate now; `PypeIt-test` uses `PYPEIT_TST`.
- **Rationale:** the two repositories will coexist, so the two locations have to be told
  apart in any case.

### 5. Where the unit and vet tests live
- **Question:** How are unit and vet tests bound to a dataset?
- **Options considered:**
  - keep the tests organized by topic, and bind each to its datasets with a pytest
    marker (e.g. `@pytest.mark.dataset('instr/setup')`) that a collection hook uses to
    select them;
  - place the tests inside each dataset directory, with shared helper functions for
    common functionality.
- **Decision:** for both unit and vet tests, place them inside each dataset directory,
  with a set of helper functions in `tests/helpers/`.  For example, `grab_img` becomes a
  central helper, and each dataset has its own `test_load_image`.
- **Refinement:** each dataset directory has a `unit/` and a `vet/` subdirectory, each
  holding several test files named by topic from a fixed list (e.g.
  `vet/test_wavelengths.py`).  This replaces the single `test_unit.py` and `test_vet.py`
  files proposed at first.
- **Rationale:** a dataset directory then holds everything about that dataset, its
  configuration, input files and tests.  The helpers keep the shared logic in one place.
  A single test file per stage would grow too large for datasets that are tested in many
  ways.  Topic files keep each dataset's tests organized, and the fixed names keep them
  consistent across datasets.

### 6. Where the development happens
- **Question:** Should the refactor be done as a series of PRs in
  `PypeIt-development-suite`, or in the new `PypeIt-test` repository?
- **Options considered:**
  - a series of PRs in the dev suite, splitting out the testing code afterward;
  - a new `PypeIt-test` started from scratch;
  - a new `PypeIt-test` seeded from the dev suite with its history intact.
- **Decision:** develop in `PypeIt-test`, seeded from the dev suite with its history
  intact, kept locally at first, and synced with the dev suite by merging until cut-over.
  The team has agreed.  How it is seeded is decision 13.
- **Rationale:**
  - Refactoring in the dev suite would move roughly 300 input files and 50 test files
    into per-dataset directories, then move them a second time when the repo is split.
  - Starting from scratch would lose the input files' history.  It would also require
    every change made in the dev suite during the refactor to be copied by hand.
  - Seeding from the dev suite keeps that history and lets its ongoing development be
    merged in.

### 7. Scope of the restructure
- **Question:** `PypeIt-test` deletes the development-only material.  Should the
  restructure still move that material into a `dev/` directory and update the dev
  suite's documentation?
- **Options considered:**
  - the full restructure: `tests/` and `dev/`, plus the dev-suite README, AGENTS and
    skills updates;
  - only the `tests/` half.
- **Decision:** restructure only the `tests/` half, and delete the development-only
  material from `PypeIt-test`.
- **Rationale:** the `dev/` moves, their `.gitignore` patterns, the `REDUX_OUT_C/`
  question and the dev-suite documentation updates would all be thrown away.  At
  cut-over, the dev suite's testing code is simply removed; its development-only material
  stays where it is (Phase 10, step 7.2).

### 8. `wavecalib/`
- **Question:** What happens to `wavecalib/` (wavelength-calibration tests with 23 MB of
  reference data, no RAW_DATA)?  Neither of its tests is runnable, and its data is not a
  RAW_DATA dataset.
- **Options considered:**
  - drop it from `PypeIt-test`;
  - fix it and keep it as data-only datasets;
  - treat it as development material.
- **Decision:** drop it from `PypeIt-test`.  It stays in the dev suite, including after
  cut-over: it is not among the testing files removed from the dev suite in Phase 10,
  step 7.2.
- **Rationale:** it hasn't run in a long time, and the per-dataset vet wavelength tests
  cover wavelength solutions.  Dropping it also takes the cleaned history from 23.8 MiB
  to 1.88 MiB.

### 9. The `pypeit_tests` test type
- **Question:** The current `pypeit_tests` type, also run by `all`, runs the main PypeIt
  repo's own pytest suite.  It belongs to no dataset.  Keep it?
- **Options considered:**
  - keep it as a global stage;
  - drop it and leave the PypeIt suite to developers and CI.
- **Decision:** keep it as the one global stage, run once before the dataset jobs.
- **Rationale:** `./pypeit_test all` stays a single command that checks everything.

### 10. Cloud testing during development
- **Question:** Nautilus jobs clone from GitHub, but `PypeIt-test` stays unpublished
  until cut-over, and its history is cleaned before it is published.  How is the cloud
  path tested in the meantime?
- **Options considered:**
  - a temporary private GitHub repository used only by Nautilus;
  - a temporary public repository in a personal GitHub account, accepting no pull
    requests;
  - local runs only, testing the cloud path at cut-over.
- **Decision:** a temporary public repository in a personal GitHub account, accepting no
  pull requests.  It is deleted at cut-over.
- **Rationale:** it isn't clear how Nautilus would be given access to a private
  repository.  A public repository needs no credentials.  Its history is the dev suite's,
  which is already public.

### 11. Where the agent skills live
- **Question:** `.claude/skills/` in the dev suite holds skills for all PypeIt
  development (`add-devsuite-setup`, `run-dev-suite`, `new-spectrograph`,
  `add-parameter`, …).  `PypeIt-test` would inherit a copy, leaving two copies to
  diverge.  Where do they live after the split?
- **Options considered:**
  - keep them in one repository, with the others pointing at it;
  - move them to the main PypeIt repo;
  - split them by the repository each one operates on.
- **Decision:** they stay in the dev suite and are not copied into `PypeIt-test`.  The
  dev suite stays in the VS Code workspace, and `PypeIt-test` is added to it once seeded.
- **Rationale:** a single copy, available to both repositories through the shared
  workspace.  The skills that operate on the test suite (`add-devsuite-setup`,
  `run-dev-suite`) are rewritten in place for `PypeIt-test` at cut-over.

### 12. Where this plan lives
- **Question:** Where is this plan kept, and should `PypeIt-test` carry a copy?
- **Options considered:**
  - commit it to the dev suite's `develop`;
  - commit it to the `dev-suite-restructure` branch;
  - also carry a copy in `PypeIt-test`, e.g. under `doc/`.
- **Decision:** commit it to the `dev-suite-restructure` branch, with no copy in
  `PypeIt-test`.  The branch holds only the plan, is kept until cut-over, and is merged
  into `develop` just before the dev suite's testing code is removed.  It keeps its
  name, since it is already pushed.
- **Rationale:** with both repositories in one workspace, the plan is reached from a
  single source.  Merging the branch at cut-over keeps the plan's edit history in
  `claude_prompts/`, which the dev suite's AGENTS.md describes as the historical record
  of agent-aided development.

### 13. How `PypeIt-test` is seeded, and what connects the two repositories
- **Question:** Should `PypeIt-test` be seeded from a branch on which the restructure has
  already been done, or seeded first and restructured itself?  And which dev-suite branch
  does `PypeIt-test` merge from?
- **Options considered:**
  - restructure on `dev-suite-restructure`, seed `PypeIt-test` from it, and merge from
    `devsuite/develop`;
  - seed first, restructure in `PypeIt-test`, and use `dev-suite-restructure` as a relay
    branch, kept up to date with `develop` and merged into `PypeIt-test`;
  - seed from `develop`, restructure in `PypeIt-test`, merge directly from
    `devsuite/develop`, and keep `dev-suite-restructure` for the plan only.
- **Decision:** seed from `develop`, restructure in `PypeIt-test`, merge directly from
  `devsuite/develop`, and keep `dev-suite-restructure` for the plan only.
- **Rationale:**
  - Upstream changes reach `PypeIt-test` as soon as they land on `develop`, with no relay
    branch to keep current.
  - The plan's commits never enter `PypeIt-test`'s history, so editing it never causes
    sync conflicts.
  - Merging `dev-suite-restructure` by mistake would only add this document, not
    restructure the dev suite.
  - Git follows the moved files equally well whether the moves happen on a dev-suite
    branch or in `PypeIt-test`.

### 14. The testing directory name
- **Question:** What should the single top-level directory for testing code and inputs be
  called?
- **Options considered:** `tests/`, `devsuite/`, `testing/`.
- **Decision:** `tests/`.
- **Rationale:** it is conventional and matches the main repo's `pypeit/tests`.

### 15. Where `RAW_DATA` lives
- **Question:** Should `RAW_DATA`, the git-ignored raw-frame directory or symlink, stay
  at the top level or move into the testing directory?
- **Options considered:**
  - keep it at the top level, with the other git-ignored data and working directories
    (`Telluric`, `REDUX_OUT`, `QL_CALIB`);
  - move it into the testing directory.
- **Decision:** keep it at the top level.
- **Rationale:** every developer's symlink and the Nautilus rclone step point there.
  Moving it would break both and gain nothing.
