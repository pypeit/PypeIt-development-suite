# Dev-suite bug fixes (November 2026)

This document seeds a branch that fixes the live bugs found while planning the test-suite
restructure (`claude_prompts/dev-suite-restructure.md`, on the `dev-suite-restructure`
branch).  The fixes are independent of that plan and do not gate any of its phases.

## Where and how

- **Branch** from the dev suite's `develop`, and merge back into `develop` through the
  normal PR process.  The restructure's next sync then brings the fixes into
  `PypeIt-test`.
- **Conventions.**  Follow the dev suite's `AGENTS.md`, "Coding Conventions", for the
  lines each fix touches.  Do not restyle surrounding code.
- **Keep `README.rst` consistent with the code** (dev suite `AGENTS.md`, "Running the
  tests").  Fix 5 changes `pypeit_test`'s exit status, so check whether the README
  describes it.
- **Release notes.**  The dev suite has no in-development release notes; describe each
  fix in the PR.

Line numbers below refer to `develop` at `c2568392` (the merge of `agent_file`, PR #443).
All seven bugs were confirmed there on 2026-11-23.

## The fixes

### 1. Duplicate dict keys disable two afterburn tests
- **Where:** `test_scripts/test_setups.py`, `_coadd1d` (L288–292) and `_coadd2d`
  (L308–311).
- **Problem:** `'gemini_gmos'` appears twice in `_coadd1d`, and `'keck_nires'` twice in
  `_coadd2d`.  The second entry silently replaces the first, so the `gemini_gmos`
  `GS_HAM_R400_860` coadd1d test and the `keck_nires` `ABBA_wstandard` coadd2d test never
  run.
- **Fix:** merge each pair into one key holding both setups, e.g.
  ```python
  'gemini_gmos': {
      'GS_HAM_R400_860': [{}], 'GS_HAM_R400_700': [{}]},
  ```
- **Check:**
  - Scan the module for duplicate keys with `ast`; the scan should find none.
  - Run `./pypeit_test reduce afterburn -s gemini_gmos/GS_HAM_R400_860 keck_nires/ABBA_wstandard`.
- **Expect:** these two tests have not run in a long time.  The input files
  (`coadd1d_files/gemini_gmos_gs_ham_r400_860.coadd1d`,
  `coadd2d_files/keck_nires_abba_wstandard.coadd2d`) may be stale, and the tests may fail.
  Fix them, or record them as known failures in the PR.

### 2. `test_datacube` asserts can never fail
- **Where:** `vet_tests/test_datacube.py`, L78, 86, 87, 96, 116, 124, 125, 140 and 149.
- **Problem:** each is written as `assert(cond, msg)`, which asserts a non-empty tuple and
  therefore always passes.
- **Fix:** `assert cond, msg`.  Touch only the assert statements.
- **Check:** `pytest vet_tests/test_datacube.py` after a KCWI reduction.  Also confirm that
  forcing one condition false makes the test fail.
- **Expect:** the residual checks (L140, L149) may now fail on current data.  Fix them, or
  record them as known.

### 3. `test_frametype::test_hires` checks only one setup
- **Where:** `unit_tests/test_frametype.py`, L476.
- **Problem:** the first setup name is `'lJ0306+1853_U074_RED_C2_ECH_0.72_XD_1.42_1x3'`,
  with a stray leading `l`.  The `glob` on it matches nothing, silently, so only the
  `J1218+2951…` setup is checked.
- **Fix:**
  - Correct the name to `'J0306+1853_U074_RED_C2_ECH_0.72_XD_1.42_1x3'`, as registered in
    `test_scripts/setups.py`.
  - Assert that each name matches exactly one directory, so a typo can't hide again.
- **Check:** `pytest unit_tests/test_frametype.py::test_hires`, and confirm both setups are
  compared.

### 4. `PypeItSensFuncTest.run` raises a `NameError`
- **Where:** `test_scripts/pypeit_tests.py`, L502.
- **Problem:** when no standard file matches, the error message uses `search_pattern`,
  which is undefined.  The test crashes with a `NameError` instead of reporting the
  missing file.
- **Fix:** report the path and pattern actually searched, e.g.
  `f"Could not find std file matching {search_path / self.std_file}"`.  `run` should also
  return the failed status, as the success branch does.
- **Check:** run `PypeItSensFuncTest.run` with a `std_file` pattern that matches nothing.
  It should report the failure, not raise.  Consider adding this case to
  `test_scripts/test_pypeit_test.py`.

### 5. pytest failures don't set the exit code
- **Where:** `test_scripts/test_main.py`: `run_pytest` (L575–622) and `main` (returns
  `test_report.num_failed`, L1036).
- **Problem:**
  - `pypeit_tests`, `unit` and `vet` run as pytest subprocesses.  Their results are only
    scraped from the console output into `TestReport.pytest_results`.
  - `num_failed` counts only `PypeItTest` failures (L315).  So `pypeit_test` exits 0 when
    pytest tests fail, and CI or Nautilus jobs report success.
- **Fix:**
  - Record each pytest run's return code (`p.returncode`) in `TestReport`.
  - Make `main` return non-zero if any pytest run failed, as well as for failed
    `PypeItTest`s.
  - Treat pytest's "no tests collected" (return code 5) explicitly, e.g. as a failure
    when tests were expected.
- **Also check:** the read loop (`while p.poll() is None: readline()`) stops when the
  process exits.  Any lines still buffered at that point, which can include pytest's
  final summary line, are never read.  Drain the remaining output after the loop.
- **Check:**
  - Add a deliberately failing test to a scratch directory and run it through
    `run_pytest`.  `pypeit_test` should exit non-zero.
  - A passing run still exits 0.

### 6. `build_ql_calibs` fallback path is broken
- **Where:** `build_ql_calibs`, `run_devsuite` (L225–255).  It runs when the QL
  calibrations can't be copied from a previous dev-suite run.
- **Problem:** it fakes `pypeit_test`'s parsed arguments with a `namedtuple`, but two of
  the values are wrong:
  - `outputdir` is a `str` (`os.path.abspath`).  `build_test_setup` computes
    `pargs.outputdir / instr / setup_name` (`test_main.py`, L1080), which raises a
    `TypeError`.  `main` converts it to a `Path` first (L785); `run_devsuite` doesn't.
  - `coverage=False`.  `PypeItTest` tests `pargs.coverage is not None`
    (`pypeit_tests.py`, L55), so every test would run under `coverage run`.
- **Fix:**
  - Pass `outputdir=Path(redux_dir).absolute()` and `coverage=None`.
  - Check that the namedtuple has every attribute that `build_test_setup` and the
    `PypeItTest` subclasses read.
- **Check:** with no existing `REDUX_OUT` products, run
  `./build_ql_calibs keck_nires --force_build` (or the equivalent) and confirm it builds
  the calibrations.

### 7. `wavecalib/` is not runnable (optional)
- **Where:** `wavecalib/test_longslit.py` and `wavecalib/test_ThAr.py`.
- **Problem:**
  - `test_longslit.py::test_wavecalib(name, spec_file, lines, wv_cen, disp, score, fidx, test=...)`
    takes arguments that pytest can't supply, so collecting it errors.
  - `test_ThAr.py` has no test function.
- **Fix, only if the dev suite wants to keep running it:** parametrize `test_wavecalib`
  over its cases, and turn `test_ThAr.py`'s checks into test functions or delete it.
  `wavecalib/` is not carried into `PypeIt-test` (decision 8 of the restructure plan),
  and nothing else depends on it.
- **Check:** `pytest wavecalib` collects and runs.

## Known failures (record in the PR)

- **`test_scripts/test_pypeit_test.py`:** all 10 tests that predate fix 4 fail on `develop`
  (checked 2026-11-23); only fix 4's new `test_sensfunc_missing_std_file` passes.  The module
  has not been part of the development workflow and has drifted:
  - Hard-coded `main()` return values no longer match, e.g. `test_main_with_test_failure`
    expects 3 and gets 6, because `shane_kast_blue` has gained setups.
  - The mocked `subprocess.run` output reaches the space-usage code, which fails to parse
    it as an integer.

  Out of scope: `pypeit_test` is due to be refactored.  Fix 5 changes `main`'s return value
  to 1 on any failure, so these return-value assertions will be out of date in a new way.

## Related, out of scope

- **`construct_template.rst`, in the main PypeIt repo**, points at
  `$PYPEIT_DEV/dev_algorithms/wavelengths/template_files/`.  That should be
  `$PYPEIT_DEV/pypeitdev/wavelengths/template_files/`.  It is a main-repo doc fix,
  listed alongside these fixes in the restructure plan's Phase 0, and needs its own PR.

## Status

| # | Fix | Status |
|---|---|---|
| 1 | Duplicate dict keys | fixed; reductions pending |
| 2 | `test_datacube` asserts | fixed; KCWI reduction pending |
| 3 | `test_hires` typo | fixed |
| 4 | `PypeItSensFuncTest` `NameError` | fixed |
| 5 | pytest exit code | open |
| 6 | `build_ql_calibs` fallback | open |
| 7 | `wavecalib/` (optional) | open |
