# Speed up PypeIt -- PR B (Detector parallelism v1)

## Goals

Implement **workstream B** of the locked speed-up plan — the headline win. Add a
small helper that maps a per-detector function over the detector list with a
`concurrent.futures.ProcessPoolExecutor` when `ncpu>1` (and the *literal* current
serial loop when `ncpu<=1`), then convert the three clean flat per-detector
stages of `pypeit/exposure.py:reduce_exposure` to use it: stage 1 (calibration
build), stage 2 (`process_exposure`), stage 4 (`extract_exposure`). Stage 3
(`findobj_on_exposure`) stays serial because of its cross-detector slitmask
barrier. Together the converted stages are ~70–77% of the Keck/DEIMOS wall-clock.
Also fixes the latent `detectors.remove(det)`-during-iteration bug in the
calibration loop.

## Claude

### Skills

Consider using the skills in `PypeIt-development-suite/.claude/skills/` — in
particular `run-dev-suite`, `build-docs`, `update-changelog`, and
`add-devsuite-setup` (for the 2-detector identical-output regression).

## Context

Authoritative documents (read the coding doc section first; it is the
implementation contract):

- **`PypeIt-development-suite/pypeitdev/speed_up/Reports/speed_up_coding.md` §2**
  — "PR B — Detector parallelism v1". Sub-anchors §B.1 (`pypeit/parallel.py`,
  with a complete code sketch), §B.2 (stage 1 + the failed-detector bug), §B.3
  (stage 2), §B.4 (stage 4), §B.5 (why stage 3 stays serial), §B.6 (tests),
  §B.7 (validation/re-profiling), §B.8 (docs). §0.2 of that doc is a grounded
  code map (file → symbol → line numbers on `speed_up` at `e9ed85c1a`); §7 is
  the risks/revert table.
- **`Reports/speed_up_design.md` v0.5 (locked)** — §0 locked decisions, §4 Tier 1
  (the 4-stage structural map and the stage-3 barrier), §8 roadmap.
- **Profile reports** — `Reports/keck_deimos_600zd_m_6500_profile_report.md`
  (12 390 s total, 4 mosaics `[(1,5),(2,6),(3,7),(4,8)]`; calibration build
  ~5 800 s ≈ 47%, local sky + extract 2 952 s ≈ 24%, image proc/mosaic
  ~700–900 s, object find + global sky ~2 400 s **serial**) and
  `Reports/shane_kast_blue_A_profile_report.md` (119 s, single detector — gains
  nothing from this PR by design).
- **Profiling scripts** — `pypeitdev/speed_up/scripts/profile_kast_blue.py`,
  `profile_deimos.py`, `analyze_profile.py`.
- **PR A results** — `Reports/speed_up_results.md` (baseline / after-A numbers),
  and `claude_prompts/speed_up/A.md` Logs.

Repository and branch:

- Code repo: `/mnt/tank/Astronomy/PypeIt/PypeIt` (same as
  `/home/xavier/Projects/PypeIt/PypeIt`).
- **Branch stack:** A branched from `speed_up` as `speed_up_qa`; **B branches
  from `speed_up_qa` as `speed_up_detpar`**; C will branch from `speed_up_detpar`
  as `speed_up_arcfit`. **PR B targets `speed_up_qa`** (retarget to `speed_up`
  once A merges).
- PR A must be complete on `speed_up_qa` before starting: this PR *reuses* the
  `par['rdx']['ncpu']` parameter, the `run_pypeit --ncpu` flag, and
  `qa.init_qa_pool`.

Locked decisions relevant to PR B (a fresh session needs no other context):

- **`ncpu=1` must run the literal current serial code.** The pool is only
  entered when `ncpu>1`. Inside `map_over_detectors`, the `ncpu<=1` branch is a
  plain `for` loop in the calling process — no pool, no pickling, no log capture.
- **Identical output is a hard requirement**: `ncpu=1` and `ncpu>1` products must
  match **to machine precision** (use exact array comparison in the tests, not
  `allclose`). This is why workers pin BLAS/OpenMP to 1 thread — thread count can
  change floating-point reduction order.
- **Worker thread pinning** via `OMP_NUM_THREADS` etc. *and*
  `threadpoolctl.threadpool_limits(1)` (needed for runtimes already initialised
  at fork time). `threadpoolctl` is already present transitively via
  `scikit-learn`; add it explicitly to `pyproject.toml` dependencies. Per the
  locked design: **be prepared to revert the pinning if it backfires.**
- **Linux `fork` start method** for v1. `spawn`/cross-platform is v2. Guard:
  if `'fork' not in multiprocessing.get_all_start_methods()`, warn and fall back
  to the serial loop.
- **Stage 1 workers return only `(det, success, failed_step)`** — the
  `Calibrations` object never crosses the process boundary; calibrations persist
  to `Calibrations/*.fits` and downstream stages reload them via `reuse_calibs`.
  The mapper supplies the `det` key, so the worker payload is
  `(success, failed_step)`.
- **Stage 3 (`findobj_on_exposure`) stays serial** — fan-out → `adjust_for_slitmask`
  cross-detector barrier → fan-out. Deferred to v2. Add the explanatory comment.
- **Per-worker buffered logging**: each worker captures its records in memory and
  ships them back; the parent replays them **grouped in detector order** so each
  detector's block stays contiguous (the dashboard parses this log). The parent
  also emits coarse live "started/finished detector X" lines.
- **Result assembly is always in `detectors` order**, never completion order,
  even though futures are reaped with `as_completed`.
- **Q4 (accepted):** stage-4 workers return the `Spec2DObj` as today. **Measure**
  the pickling cost in this PR (a `Spec2DObj` carries ~8 full-frame float arrays,
  ~250 MB per DEIMOS mosaic); if marshalling exceeds **~5% of the stage**, switch
  the worker to `spec2DObj.to_file(Intermediate/…)` and return the path.
- **Q5 (accepted):** QA is written **serially inside detector worker processes**
  — `parallel._worker_init` calls `qa.init_qa_pool(1)`, since the parent's
  threads do not survive the fork. QA threading remains only in the
  single-process serial stages.
- **Q6 (accepted):** the real multi-detector identical-output regression lives in
  the dev suite on a **cheap 2-detector setup** (e.g.
  `keck_lris_blue/multi_600_4000_d560`) as a `reduce` variant plus a `vet_tests`
  comparison; a DEIMOS `--ncpu 4` run is done **manually** before merging to
  `develop`, not in the suite.
- **Out of scope:** the B-spline kernel (`pypeit/core/bspline/`) — another branch.
  Do not touch it. Do not change any algorithm; the loop bodies (`calib_one`,
  `process_one_det`, `extract_det`) are untouched.

## Running

If you need to run Python or PypeIt, use the **`pypeit` conda environment**.
The shell's default environment shadows it, so use the absolute env binaries:

```bash
/home/xavier/miniconda3/envs/pypeit/bin/python
/home/xavier/miniconda3/envs/pypeit/bin/run_pypeit
/home/xavier/miniconda3/envs/pypeit/bin/pytest
```

## Prompts

1. Read this doc.  Perform the 1st task under Tasks.
2. Read this doc.  Perform the 2nd task under Tasks.
3. Read this doc.  Perform the 3rd task under Tasks.
4. Read this doc.  Perform the 4th task under Tasks.
5. Read this doc.  Perform the 5th task under Tasks.
6. Read this doc.  Perform the 6th task under Tasks.
7. Read this doc.  Perform the 7th task under Tasks.

## Tasks

1. **Prepare.** Create the branch and ground yourself; write no production code
   yet.
    - Confirm PR A is complete on `speed_up_qa` (the `ncpu` parameter, the
      `--ncpu` flag and `qa.init_qa_pool` all exist). Then create and check out
      **`speed_up_detpar`** from **`speed_up_qa`**; report the HEAD commit.
    - Read `Reports/speed_up_coding.md` **§0**, **§2** and **§7** (risks/revert)
      in full.
    - Read the code you are about to change: all of `pypeit/exposure.py`
      (especially `reduce_exposure`'s stage-1 loop, `process_exposure`,
      `findobj_on_exposure` including the `adjust_for_slitmask` barrier, and
      `extract_exposure`), the relevant functions in `pypeit/pypeit_steps.py`
      (`calib_one`, `process_one_det`, `extract_det`,
      `load_calibrations_for_frame`), `pypeit/pkg/logger.py` (the logger is a
      stock `logging.Logger` subclass — buffered worker logging is a handler
      swap), and `pypeit/pypeit.py:reduce_calibID`.
    - Verify the line numbers quoted in coding-doc §0.2/§2 still hold after PR A;
      note any drift.
    - Verify `threadpoolctl` imports in the `pypeit` env and record its version.
    - Ask questions in the Q&A section below. Log your work in the Logs section.

2. **Create `pypeit/parallel.py` and its unit tests.** Coding-doc **§B.1** and
   the `test_parallel.py` part of **§B.6**.
    - Implement the module exactly as sketched: `THREAD_ENV_VARS`,
      `_FORK_PAYLOAD`, `nworkers`, `_RecordBuffer`, `_worker_init` (BLAS env
      pinning + detach the log handlers inherited from the parent + detach the
      `py.warnings` handlers + `qa.init_qa_pool(1)`), `_run_one` (buffer attach,
      `threadpool_limits(1)`, `qa.flush_qa()`, traceback capture),
      `map_over_detectors`, and `surviving_detectors`.
    - Note the two subtleties: the pool is entered whenever `ncpu>1` **even for a
      single detector** (so the parallel path gets CI coverage on the
      single-detector shane_kast_blue reduction), and a **fresh executor is
      created per call** so `_FORK_PAYLOAD` is never stale.
    - Add the `'fork' not in mp.get_all_start_methods()` guard that falls back to
      the serial loop with a warning.
    - Add `threadpoolctl` to `[project] dependencies` in `pyproject.toml`.
    - Write `pypeit/tests/test_parallel.py` with all the tests in §B.6:
      `test_nworkers`, `test_serial_and_parallel_are_bit_identical`,
      `test_parallel_is_deterministic`, `test_payload_is_per_detector`,
      `test_worker_threads_are_pinned`,
      `test_worker_logs_replayed_in_detector_order`,
      `test_worker_exception_names_the_detector`, and
      `test_surviving_detectors_rebuild`. Run them; they must be green.
    - Ask questions in Q&A. Log your work in the Logs section.

3. **Add the detector-first adapters in `pypeit_steps.py`.** Coding-doc **§B.2**,
   **§B.3**, **§B.4**.
    - `calib_status_one(det, ...)` — wraps `calib_one`, logs
      `f'Calibrating detector {det}'`, returns `(success, failed_step)` only.
    - `process_one_det_bydet(det, ...)` — wraps `process_one_det`, logs
      `f'Reducing detector {det}'`, returns `(sciImg, bkg_redux_sciimg)`.
    - `extract_det_bydet(det, ...)` — wraps `extract_det`, taking the
      per-detector objects (`sciImg`, `final_sky`, `sobjs_obj`, `calib_slits`,
      `bkg_redux_final_sky`) as keyword arguments supplied through the mapper's
      `payload`.
    - All three must be **module-level** functions with full Numpy-style
      docstrings (a closure happens to work under `fork`, but must not be relied
      on — `spawn` support is a v2 goal).
    - Do not convert any caller yet.
    - Ask questions in Q&A. Log your work in the Logs section.

4. **Convert stage 1 (the calibration loop) and fix the failed-detector bug.**
   Coding-doc **§B.2**.
    - In `pypeit/exposure.py:reduce_exposure`, replace the calibration loop with
      `parallel.map_over_detectors(pypeit_steps.calib_status_one, …,
      ncpu=par['rdx']['ncpu'], …)`, add `from pypeit import parallel`, then
      **rebuild the surviving-detector list after the loop** with
      `parallel.surviving_detectors`, preserving the existing warning message
      verbatim.
    - This removes the latent bug: `detectors.remove(det)` was called *while
      iterating* `detectors`, so with the DEIMOS mosaics `[(1,5),(2,6),(3,7),(4,8)]`
      a failure on `(2,6)` silently skipped `(3,7)` — which then stayed in the
      list and was reduced with calibrations that were never built.
    - Add the early return when every detector fails (an empty
      `AllSpec2DObj`/`SpecObjs`; `reduce_calibID` already guards on
      `len(this_spec2d.detectors) > 0`).
    - Add the explanatory comment to `findobj_on_exposure` recording why stage 3
      is deliberately **not** parallelized (coding-doc §B.5).
    - Smoke-test: run `shane_kast_blue` end-to-end at `--ncpu 1` and `--ncpu 2`
      and confirm both complete.
    - Ask questions in Q&A. Log your work in the Logs section.

5. **Convert stages 2 and 4.** Coding-doc **§B.3** and **§B.4**.
    - `process_exposure`: replace the detector loop with the mapper and assemble
      `sciImg_dict` / `bkg_redux_sciimg_dict` in the fixed `detectors` order.
    - `extract_exposure`: slice the per-detector inputs **in the parent** into a
      `payload` dict (`sciImg`, `final_sky`, the `all_specobjs_objfind` selection
      for that `detname`, `calib_slits[i]`, `bkg_redux_final_sky`), call the
      mapper, then assemble `all_spec2d` / `all_specobjs_extract` in `detectors`
      order. Preserve the existing `all_specobjs_extract.calibs` assignment
      verbatim, including the fact that it overwrites each iteration — changing
      that is out of scope.
    - **Q4 measurement:** instrument the stage-4 return path (e.g. time
      `pickle.dumps` on one `Spec2DObj`, or compare summed worker time against
      wall time for the stage) on a DEIMOS mosaic. Record the fraction of the
      stage spent marshalling. **If it exceeds ~5%,** switch
      `extract_det_bydet` to write the per-detector `spec2d` to `Intermediate/`
      and return the path, with the parent reloading it. Report the measurement
      either way.
    - Smoke-test `shane_kast_blue` again at `--ncpu 1` and `--ncpu 2`.
    - Ask questions in Q&A. Log your work in the Logs section.

6. **Tests, dev-suite regression, docs and changelog.** Coding-doc **§B.6** and
   **§B.8**.
    - Add `test_run_pypeit_ncpu` to `pypeit/tests/test_runpypeit.py`: reduce
      shane_kast_blue twice, serially and with `--ncpu 2`, into **separate**
      output directories, and require the `spec1d` arrays to match with
      `np.testing.assert_array_equal` (exact — any difference means the thread
      pinning is not working). Guard it with `PYPEIT_SKIP_SLOW_TESTS`.
    - Per **Q6**, add the dev-suite multi-detector regression: a `reduce`
      variant on a cheap **2-detector** setup (e.g.
      `keck_lris_blue/multi_600_4000_d560`) run with `--ncpu 2` alongside the
      serial run, registered in `test_scripts/test_setups.py`, plus a `vet_tests`
      test that compares the serial and parallel `spec1d`/`spec2d` arrays
      element-by-element. Use the `add-devsuite-setup` skill.
    - Run `pytest pypeit/tests` in the `pypeit` env; must be green.
    - Docs: extend the "Running on multiple CPUs" section of `doc/running.rst`
      (what is and is not parallelized, the RAM-per-detector tradeoff, the
      Linux-only `fork` caveat); regenerate `doc/api` (`./update_docs`, then
      `git add doc/api`); check `doc/sphinx_warnings.out` is empty.
    - Changelog: add the §B.8 bullets to `doc/releases/2.1.0dev.rst` under
      **Functionality/Performance Improvements and Additions**, **Bug Fixes**
      (the `detectors.remove()` bug), **Under-the-hood Improvements**
      (`pypeit.parallel`), **Dependency Changes** (`threadpoolctl`) and
      **Testing**.
    - Ask questions in Q&A. Log your work in the Logs section.

7. **Validate, re-profile and summarize.** Coding-doc **§B.7** and the
   cross-cutting checklist in **§4** of the coding doc.
    - **Measurement caveat:** the profiling scripts run cProfile in-process, so
      child processes are **not** profiled and the `--ncpu>1` `.prof` is blind to
      most of the work. Use `wall_clock_s` from `Reports/*.runmeta.json` as the
      primary metric; keep cProfile for the `--ncpu 1` baseline only.
    - Profile DEIMOS at `--ncpu 1` (baseline, wall + cProfile) and `--ncpu 4`
      (wall only). Profile Kast at `--ncpu 1` and `--ncpu 2` to confirm no
      regression beyond a second or two of fork overhead.
    - **Expected (§B.7):** parallel fraction p ≈ 0.77, ideal 4-worker speedup
      2.4×, realistic **1.9–2.3×** after load imbalance (per-mosaic calibration
      spans 1 338–1 584 s) and marshalling — i.e. **DEIMOS ~11 700 s (post-A) →
      ~5 200–6 200 s**.
    - Record **peak RSS** at `--ncpu 1` and `--ncpu 4` on DEIMOS
      (`/usr/bin/time -v` or a `psutil` sampler) and document the RAM-vs-`ncpu`
      tradeoff in `doc/running.rst`.
    - **Identical-output check:** diff the `--ncpu 1` and `--ncpu 4` DEIMOS
      `spec1d`/`spec2d` arrays — they must be bit-identical. **Determinism
      check:** run `--ncpu 4` twice and diff — also bit-identical. If either
      fails, suspect the BLAS pinning first (see the revert plan in coding-doc
      §7).
    - Run the relevant dev-suite tests: `./pypeit_test reduce -i shane_kast_blue
      keck_deimos keck_lris_blue`, then
      `pytest vet_tests --redux_out $PYPEIT_DEV/REDUX_OUT`.
    - Append the numbers to `Reports/speed_up_results.md` (new "after-B"
      section) and summarize in the Logs section below. Ask any remaining
      questions in Q&A.

## Q&A

Claude poses questions here (as `**Qn — title.** body` followed by a `>A:` line);
the user answers inline beneath each `>A:`.

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs
