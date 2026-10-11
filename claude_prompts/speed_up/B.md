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

8. **Issue** Create an Issue on GitHub as `profxj` that describes the bspline problem(s) that you have 
     uncovered.  Explain what you did, why and what you found.
     Use Opus 5.5. Log your work.

9. **PR** Create a PR on GitHub to `develop` as `profxj` for all of the work in this doc.  
     Describe what was done and why in the PR description.
     Use Opus 5.5. Log your work.

## Q&A

Claude poses questions here (as `**Qn — title.** body` followed by a `>A:` line);
the user answers inline beneath each `>A:`.

**Q1 — Branch name: `speed_up_B` already exists; do you want `speed_up_detpar` too?**
The working tree was already on a branch called **`speed_up_B`**, which *is* a
descendant of `speed_up_qa` (`git merge-base --is-ancestor speed_up_qa
speed_up_B` → true; `speed_up_qa` landed on it via merge commit `96b175ab8`,
"Merge pull request #2198 from pypeit/speed_up_qa") and additionally carries an
up-to-date `develop` merge (MMIRS ramp work, the setup-colon-quote fix, the BUNIT
propagation, doc fixes).  Branching a *new* `speed_up_detpar` off `speed_up_qa`
as §0.1 specifies would throw that `develop` merge away and give PR B a staler
base.  **I have therefore stayed on `speed_up_B` and written no new branch.**
Tell me if you would rather I (a) keep working on `speed_up_B`, (b) rename it
`git branch -m speed_up_B speed_up_detpar`, or (c) genuinely create
`speed_up_detpar` from `speed_up_qa` and lose the `develop` merge.
>A: (a)

**Q2 — PR A shipped a different QA API than coding-doc §B.1 assumes.**
The §B.1 sketch calls `qa.init_qa_pool(1)` in `_worker_init` and `qa.flush_qa()`
in `_run_one`.  Neither function exists: PR A landed as a **`QAWriter` class**
(`pypeit/pkg/qawriter.py`) with a package-level singleton `pypeit.qaWriter`, and
the equivalents are **`qaWriter.init(ncpu=1)`** and **`qaWriter.flush()`**.
`QAWriter.init` is explicitly documented as fork-safe ("a process forked from one
that was using a pool inherits an instance whose threads do not exist in the
child, so the child must call this"), so it is a drop-in for the sketch's intent.
Unless you object, `pypeit/parallel.py` will use `qaWriter.init(ncpu=1)` /
`qaWriter.flush()` and I will note the substitution in the module docstring.
>A:  That looks fine.  

**Q3 — Should `calib_status_one` forward `status_only` and `reload_only`?**
`calib_one` has grown two keyword arguments since the coding doc was written:
`status_only` and `reload_only` (`pypeit/pypeit_steps.py:170-174`).  The §B.2
sketch of `calib_status_one` predates them and forwards only
`reuse_calibs/qa_path/show/run_state/stop_at_step`.  My default is to **forward
all of them**, so the adapter is a complete stand-in for `calib_one` (and so the
optional `PypeIt.calib_all` conversion in §B.5, which passes `status_only` and
`reload_only`, stays possible).  Object if you want the adapter kept narrow.
>A: Use your default

**Q4 — Convert `PypeIt.calib_all` in this PR, or leave a TODO?**
§B.5 says the second flat detector loop in `PypeIt.calib_all`
(`pypeit/pypeit.py:236-247`) is "an easy, optional bonus conversion … do it in
this PR only if it costs nothing; otherwise leave a TODO."  It is not free: the
loop body keeps the returned `caliBrate` object and assigns `self.det` as a side
effect, so converting it means auditing what (if anything) downstream reads
`self.det` / the `caliBrate` return.  My default is to **leave a TODO** and keep
PR B focused on `reduce_exposure`'s three stages.  Say the word if you want it
converted.
>A: Leave a TODO

**Q5 — Is the `detectors` list aliased anywhere the caller can observe?**
The failed-detector fix replaces the in-place `detectors.remove(det)` with a
rebuilt list (`parallel.surviving_detectors`).  In `reduce_exposure` the list
comes straight from `spectrograph.select_detectors(...)`
(`pypeit/exposure.py:415-417`), which returns a fresh `np.arange(...).tolist()`
or a fresh list of mosaic tuples, and the local name is not returned to the
caller — so rebinding rather than mutating is safe.  Flagging it only because
it is the one behavioural difference the fix introduces beyond the bug itself.
>A:  I don't know.  If you can't find it, then the answer is likely no


**Q6 — The BLAS pinning cannot make `ncpu=1` and `ncpu>1` bit-identical as currently specified. Which way do you want to resolve it?**
This is the one finding from Task 2 that I think changes the plan, so it is
written up at length.  §B.1 justifies pinning worker threads with two claims:
(a) N workers must not each spawn N threads, and (b) pinning makes "the
floating-point reduction order inside a worker match the serial run
bit-for-bit".  **Claim (b) does not hold**, because the pinning is applied
*only inside the workers*.  At `ncpu=1` the reduction runs in the parent with
OpenBLAS at its default 24 threads; at `ncpu>1` each worker runs it at 1 thread.
Those two thread counts give different answers for any BLAS reduction big enough
to be split across threads.  Measured here, on this machine:

| workload | 1 thread vs 24 threads |
|---|---|
| `np.dot(v, v)`, `v` of 2e6 elements | **differ, max 2.8e-09** |
| `A @ B`, 1500x1500 | **differ, max 1.1e-13** |
| `np.linalg.svd(A, compute_uv=False)`, 500x500 | **differ, max 1.5e-13** |
| `np.linalg.inv(A)`, 600x600 | **differ, max 4.1e-14** |
| `M.T @ M`, 3000x800 | identical |
| the §B.6 `_work` SVD-sum, 200x200 | identical |

And straight through the real helper (`ncpu=1` in the parent vs `ncpu=2` in
pinned workers), using the `ddot` workload:

```
det 1  1998756.0571849102  vs  1998756.0571849095   DIFFER 6.985e-10
det 2  2000352.6469428146  vs  2000352.6469428123   DIFFER 2.328e-09
```

So the locked "identical to machine precision, use `assert_array_equal`"
requirement is only satisfiable if **either** the reduction never performs a
thread-split BLAS reduction anywhere (quite possible — PypeIt's hot paths are
b-spline and low-order polynomial fits on small arrays, and note the 3000x800
`M.T @ M` case above *is* identical), **or** the serial path is pinned too.
Three ways out:

1. **Pin the serial path as well** — wrap the `ncpu<=1` loop in
   `threadpool_limits(1)` too.  This is the only option that actually
   guarantees the locked requirement, and it has the nice property of making
   PypeIt's output independent of `ncpu` entirely.  Cost: it breaks the other
   locked decision that "`ncpu=1` must run the literal current serial code", and
   it is a *performance regression at the default setting* — every existing
   user's `ncpu=1` run would lose multithreaded BLAS.
2. **Do not pin the workers at all** (the §7 revert).  Serial and parallel then
   both run BLAS at 24 threads and agree bit-for-bit, but 4 workers x 24 threads
   on 24 cores oversubscribes badly and will eat into the speedup PR B exists to
   deliver.
3. **Keep it exactly as specified and let Task 7 decide empirically.**  If the
   DEIMOS `--ncpu 1` vs `--ncpu 4` diff comes back bit-identical, PypeIt simply
   has no thread-sensitive BLAS reduction in the reduction path and the issue is
   theoretical.  If it does not, fall back to option 1 or 2 with real numbers in
   hand.

**My default is (3)** — it is what §7's revert plan already anticipates ("If
either fails, suspect the BLAS pinning first"), it costs nothing now, and
option 1 is a user-visible performance decision I do not want to make on your
behalf.  I have implemented §B.1 exactly as written.  Tell me if you would
rather I pre-emptively take option 1.
>A:  Go with (3).  And, we do not repeatability at the 1e-9 level.  1e-6 would be acceptable.

**Q7 — Warnings raised inside a worker are currently dropped. Accept, or buffer them too?**
`_worker_init` detaches the handlers from the `py.warnings` logger (correctly —
they point at the parent's file descriptor), but `_run_one` attaches the
`_RecordBuffer` only to the `pypeit` logger.  A `warnings.warn` raised inside a
worker therefore goes nowhere at all: not to the worker's stderr, not into the
buffer, not into the run log.  At `ncpu=1` it would have appeared in both the
console and the log file.  Buffering `py.warnings` as well is about four extra
lines (attach the same buffer to that logger in `_run_one`, replay those records
through `logging.getLogger('py.warnings')`).  §B.1 does not mention it, so I
implemented it as sketched and am flagging it rather than deviating.  Worth
fixing in Task 4, or leave it?
>A:  I do not care about py.warnings, only that the PypeIt logging is maintained

**Q8 — Which conda environment should the tests actually run in?**
The "Running" section of this prompt says to use `/home/xavier/miniconda3/envs/pypeit`.
That environment currently **cannot collect most of `pypeit/tests`**: it pairs
numpy 2.4.6 with astropy 7.0.2, and astropy's `units/quantity_helper` still
references `np.in1d`, which numpy 2.x removed, so importing `astropy.units`
raises `AttributeError` at collection time.  This is pre-existing and entirely
unrelated to PR B (`test_parallel.py` itself passes there, because it never
imports astropy).  `/home/xavier/miniconda3/envs/pypeit14b` (numpy 2.5.0,
astropy 8.0.0, `threadpoolctl` 3.6.0) collects and passes them.  I will use
**`pypeit14b`** for anything that touches the broader test suite and for the
`pytest pypeit/tests` run in Task 6, unless you tell me otherwise — but the
reductions and profiling in Tasks 4/5/7 should probably stay in whichever
environment produced the PR A baseline numbers, so please confirm which that was.
>A:  On the workstation, it is `pypeit14b`.  That should be where the PR A baselines were run.


**Q9 — Q6 is now settled empirically, and the answer is the bad one: the BLAS thread count changes the reduction *visibly*, not at 1e-9. Which fix do you want?**
You answered Q6 with "go with (3), let Task 7 decide empirically, 1e-6 would be
acceptable".  Task 5 decided it three tasks early, on `shane_kast_blue`, and the
tolerance does not save us.  Four cold runs of the same data:

| run | SKYMODEL / OBJMODEL / IVARMODEL vs serial-unpinned |
|---|---|
| serial, BLAS unpinned (`--ncpu 1`) | baseline |
| **serial, BLAS pinned to 1 thread** | **differs** |
| **parallel `--ncpu 2`** (workers pinned) | **differs — identically to the pinned serial run** |
| parallel `--ncpu 2` vs **serial pinned** | **BIT-IDENTICAL** |

Two things follow, and the second is the problem:

1. **The parallelism itself is correct.**  `--ncpu 2` reproduces a serial run
   bit-for-bit *as long as the thread count matches*.  Result assembly,
   payload slicing, marshalling and the fork are all sound — this is as strong
   a correctness result as PR B could ask for.
2. **The difference is entirely BLAS thread count, and it is large.**  On the
   first Kast exposure, `SKYMODEL` differs on 621 756 of 716 800 pixels (87%)
   and `OBJMODEL` on 582 400 (81%).  The inputs — `SCIIMG`, `IVARRAW`, `TILTS`,
   `WAVEIMG`, `BPMMASK` — are all bit-identical, so nothing is *wrong*; it is
   that `local_skysub_extract`'s iterative b-spline fit with sigma-clipping
   **amplifies a last-ulp difference into a different rejection set**, and from
   there into a visibly different model.  This is why a 1e-6 tolerance does not
   rescue the comparison: the differences are chaotic, not round-off.

So the status quo cannot satisfy "identical output" at any sane tolerance, and
Task 6's regression test cannot be written against it.  The options are the
same three as Q6, but now with numbers:

1. **Pin the serial path too** — wrap the `ncpu<=1` branch of
   `map_over_detectors` in `threadpool_limits(1)` as well.  **This is already
   demonstrated to work**: pinned-serial and `--ncpu 2` came out bit-identical
   above.  PypeIt's output then stops depending on `ncpu` *and* on how many
   cores the machine has, which is arguably a correctness win in its own right
   (today, the same data on a 4-core and a 24-core laptop give different sky
   models).  Measured cost on Kast: **77.5 s unpinned -> 81.9 s pinned, +5.7%**.
   That is a real regression at the default setting, and I do not know yet what
   it costs on DEIMOS.
2. **Do not pin the workers** — each worker runs BLAS at the machine default, so
   it matches the unpinned serial run bit-for-bit.  Costs nothing in accuracy,
   but `--ncpu 4` then runs 4 x 24 = 96 BLAS threads on 24 cores.  How much that
   oversubscription actually hurts is unmeasured; it may be modest, because the
   DEIMOS hot paths are b-spline fits rather than large BLAS calls.
3. **Accept that `ncpu` changes the sky model** and drop the identical-output
   requirement, replacing Task 6's exact test with a scientific-equivalence
   check (e.g. agreement in extracted counts to within a fraction of the noise).

**My recommendation is (1)**, with (2) measured on DEIMOS in Task 7 as the
fallback if pinning proves expensive there.  (1) is the only option that makes
the output reproducible across machines, which seems worth more than the ~6% it
costs at `ncpu=1` — and that 6% is on Kast, where there is no parallel speedup
to offset it; on DEIMOS the same pinning buys a 2x+ speedup at `--ncpu 4`.
I have **not** implemented it, because you were explicit that the user-visible
performance tradeoff is yours to make.  Tell me which and I will apply it at the
top of Task 6.
>A: Go with (1).

**Q10 — I had to fix two latent crashes in core classes to make PR B work at all. Please sanity-check the scope.**
Neither stage 2 nor stage 4 could return *anything* to the parent: two PypeIt
classes raise `RecursionError` when unpickled, which killed the worker pool
outright (`BrokenProcessPool`).  Both are the same pre-existing bug, and
**neither has anything to do with parallelism** — they are equally broken under
`copy.deepcopy` and under any other use of `pickle`.  Pickle constructs an
instance *without* calling `__init__` and then probes it for `__setstate__`;
that lookup lands in a `__getattr__` which immediately reads an attribute that
`__init__` would have set, so `__getattr__` calls itself forever:

- `pypeit/images/bitmaskarray.py:114` — `BitMaskArray.__getattr__` reads
  `self.lower_keys`.  Blocks stage 2 (every `PypeItImage` carries a bitmask) and
  stage 4.  Minimal reproduction, no PR B code involved:
  `pickle.loads(pickle.dumps(ImageBitMaskArray((20,10))))` -> `RecursionError`.
- `pypeit/specobjs.py:776` — `SpecObjs.__getattr__` reads `self.specobjs`.
  Blocks stage 4.

I fixed both the same way: read the attribute from `self.__dict__` instead, and
raise `AttributeError` when it is absent, since an instance in that state has
nothing to look up.  Verified that pickling and `deepcopy` now round-trip with
the mask array, dtype, `lower_keys` and per-flag access (`bpm`, `cr`) all
preserved, and that a genuinely missing attribute still raises `AttributeError`.
These are small, but they are edits to core datamodel classes rather than to
anything PR B owns, so flag it if you would rather they were split into their
own bug-fix PR ahead of this one — I think they should be, since they stand on
their own merit and would otherwise be buried in a performance PR.
>A: Ok, I think that is a fine solution


**Q11 — PypeIt is not reproducible run-to-run on slitmask data, independent of PR B. The locked "bit-identical" requirement is unachievable, and Task 6's dev-suite test will fail because of it. How do you want to proceed?**
The DEIMOS `--ncpu 1` vs `--ncpu 4` identical-output check **failed** (2 039 of
2 624 arrays differ).  Before blaming the parallelism I ran the control, on the
cheap 2-detector `keck_lris_blue/multi_600_4000_slitmask` setup so it could be
iterated in minutes rather than hours.  **The control fails too.**

| comparison | arrays differing |
|---|---|
| LRIS serial vs `--ncpu 2` | **637 / 842** |
| **LRIS serial vs a second identical serial run** | **637 / 842** |
| **LRIS serial vs serial, with every thread env var = 1 *and* `PYTHONHASHSEED=0`** | **637 / 842** |
| shane_kast_blue serial vs `--ncpu 2` (verified three times) | **0 / 328** |

Two runs of the *same command* differ by the same amount as a serial-vs-parallel
pair.  So the DEIMOS failure is **not caused by detector parallelism**; it is
pre-existing non-determinism that PR B merely made visible by being the first
thing to check for it.  It is also not thread count and not hash randomisation —
pinning both changes nothing.

**It localises cleanly.**  Comparing the calibration files of two identical
serial runs:

```
Arc_A_0_DET01.fits        15 arrays  identical      <- image combination is fine
Tiltimg_A_0_DET01.fits    15 arrays  identical
Flat_A_0_DET01.fits       42 arrays  30 DIFFER
Tilts_A_0_DET01.fits      14 arrays   3 DIFFER
```

and inside the flat:

```
PIXELFLAT_RAW                      same      max|d|=0.000e+00   <- identical input
PIXELFLAT_SPAT_ID-377_BSPLINE:bkpt_full   DIFFERS
PIXELFLAT_SPAT_ID-377_BSPLINE:coeff       DIFFERS
PIXELFLAT_SPAT_ID-377_BSPLINE:icoeff      DIFFERS   (and the same for all 9 slits)
PIXELFLAT_MODEL                    DIFFERS   max|d|=4.368e-03
```

The input image is bit-identical and the **B-spline fit returns different
breakpoints and coefficients from it**.  That points at
`pypeit/core/bspline/` — which this PR is explicitly forbidden to touch ("the
B-spline kernel — another branch.  Do not touch it.").  A plausible mechanism is
uninitialised memory or a tie-break in the rejection loop in the C extension,
but I have not gone further, since it is out of scope and deserves its own
investigation.

Note that `shane_kast_blue` **is** deterministic and **is** bit-identical
between `--ncpu 1` and `--ncpu 2`, which is why this never showed up before
Task 7: simple longslit data does not exercise the affected path hard enough.

**What this means for PR B.**  The evidence that the parallelism is correct is
unchanged and, I would argue, strong: on the setup where PypeIt *is*
deterministic, serial and parallel agree bit for bit; and where it is not, the
serial-vs-parallel disagreement is exactly the size of the serial-vs-serial
disagreement.  But two things now need your decision:

1. **The acceptance criterion.**  "Identical to machine precision" cannot be met
   for slitmask data by *any* implementation.
2. **The Task 6 dev-suite test.**  `vet_tests/test_ncpu.py` asserts exact
   equality on LRIS and **will fail** — not because of PR B, but because PypeIt
   cannot reproduce itself there.  As written it is a test of the B-spline bug,
   not of detector parallelism.

Options:

1. **Keep PR B, restate the requirement, and file the B-spline bug separately.**
   The exact-equality test stays where it is meaningful (`test_run_pypeit_ncpu`
   on Kast, which passes).  `vet_tests/test_ncpu.py` is rewritten to assert
   *scientific* equivalence on LRIS — agreement well within the noise, plus the
   structural checks that are still exact (same objects, same detectors, same
   extensions).  The B-spline non-determinism is filed as its own issue with the
   reproducer above.
2. **Block PR B until the B-spline non-determinism is fixed**, then re-run the
   DEIMOS check and require exact equality.  Correct in principle, but it makes
   a performance PR wait on an unrelated bug in a kernel another branch owns.
3. **Keep the exact test and mark it xfail** with a reference to the B-spline
   issue, so it starts passing by itself once that is fixed.

**My recommendation is (1)**, with the structural assertions kept exact so the
test still catches anything PR B could plausibly break (a dropped detector, a
reordered object list, a missing extension) — which is what it was really for.
I would additionally write the reproducer into the issue, since "two identical
serial runs of the same `.pypeit` file disagree" is a much sharper bug report
than anything that would have come out of the parallelism work.
>A: (1), and we will execute Prompt #8 to submit an Issue

**Q12 — The headline speedup is measured against a cProfile-instrumented baseline. Do you want a clean one?**
§B.7 says to keep cProfile for the `--ncpu 1` baseline only and use wall-clock
for `--ncpu > 1`, which is what I did — but that makes the ratio
apples-to-oranges, because cProfile inflates the baseline.  The raw numbers are
**10 515.5 s** (`--ncpu 1`, instrumented) and **4 123.6 s** (`--ncpu 4`, not
instrumented), i.e. **2.55x**, which is suspiciously *better* than the 2.37x
Amdahl ideal for p = 0.77 — a sign the baseline is inflated rather than that the
code beat Amdahl.  For scale, on Kast the same instrumentation costs about 27%
(86.3 s uninstrumented vs 109.7 s instrumented).  If DEIMOS pays anything like
that, the honest speedup is nearer **2.0x** than 2.55x.
I can settle it with one more uninstrumented `--ncpu 1` DEIMOS run (~2.9 h).  My
inclination is to run it, because a headline number that overstates the gain is
worse than a slightly smaller defensible one, but it is another three hours of
machine time and the decision is yours.
>A:  Yes, run it.

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-10-10 (Task 1 — Prepare: branch, grounding, line-number drift audit)

No production code written, as instructed.

**PR A is complete.** All three artifacts PR B depends on are present:

| Artifact | Location | Status |
|---|---|---|
| `par['rdx']['ncpu']` | `pypeit/par/pypeitpar.py:3070-3083`, validated at `3190-3192` | present, default 1 |
| `run_pypeit --ncpu` | `pypeit/scripts/run_pypeit.py:64-65`, plumbed at `100` | present |
| QA pool | `pypeit/pkg/qawriter.py` (`QAWriter`), singleton `pypeit.qaWriter` (`pypeit/__init__.py:30`), configured at `pypeit/pypeit.py:102` | present, **but not under the name the coding doc assumes** — see Q2 |

**Branch.** The tree was already on **`speed_up_B`**, which is a descendant of
`speed_up_qa` (merge commit `96b175ab8`, PR #2198) *plus* an up-to-date `develop`
merge.  I did **not** create `speed_up_detpar` from `speed_up_qa`, because that
would have discarded the `develop` merge and given PR B a staler base than the
branch you already prepared.  Working tree clean.
**HEAD: `f555483dd` — "Merge pull request #2206 from pypeit/doc_fix".**
Asked for a ruling in **Q1**.

**Read in full:** `Reports/speed_up_coding.md` §0 (scope, branch stack, code map,
profile targets), §2 (B.1–B.8), §4 (cross-cutting checklist), §5 (doc
checklist), §6 (commit sequence), §7 (risks/revert) and §8 (open questions Q1–Q6);
`pypeit/exposure.py` in full; `calib_one` / `process_one_det` / `extract_det` /
`load_calibrations_for_frame` in `pypeit/pypeit_steps.py`; `PypeItLogger.init`
in `pypeit/pkg/logger.py`; `PypeIt.calib_all` and `reduce_calibID` in
`pypeit/pypeit.py`; `pypeit/pkg/qawriter.py`.

**Line-number drift since §0.2 (which was pinned to `speed_up` at `e9ed85c1a`).**
Every anchor has moved; nothing has changed structurally.

*`pypeit/exposure.py`* (file is now 606 lines):

| Symbol | §0.2 | Now |
|---|---|---|
| `adjust_for_slitmask` (**stage-3 barrier**) | 21–80 | **23–84** |
| `process_exposure` (**stage 2**) | 84–138 | **86–140**; detector loop **127–137** |
| `findobj_on_exposure` (**stage 3**) | 140–264 | **142–266**; fan-out **204–224**, barrier call **228–233**, second fan-out **237–264** |
| `extract_exposure` (**stage 4**) | 266–348 | **268–350**; detector loop **315–347** |
| `reduce_exposure` (driver) | 350–466 | **352–506**; `detectors = select_detectors(...)` **415–417**, **stage-1 calib loop 423–435** |
| `save_exposure` | 468–565 | **508–606** |

*`pypeit/pypeit_steps.py`* (now 1047 lines):

| Symbol | §0.2 | Now |
|---|---|---|
| `calib_one` | 124–199 | **170–258** |
| `process_one_det` | 202–319 | **260–386** |
| `findobj_on_det` | 321–401 | **388–467** |
| `finalize_sky_det` | 402–481 | **469–547** |
| `load_calibrations_for_frame` | 482–536 | **549–602** |
| `extract_det` | 629–801 | **737–909**; `spec2DObj.gen_qa()` at **906** |

*`pypeit/pypeit.py`* (now 438 lines): `PypeIt.__init__` **67–182**,
`build_qa` **193–202**, `calib_all` **204–247** (detector loop **236–247**),
`reduce_all` **249–290**, `reduce_calibID` **304–438**.  The empty-container
guard the §B.2 early return relies on — `if len(this_spec2d.detectors) > 0:` —
is now at **`pypeit/pypeit.py:423`** (§0.2 said 369).  `AllSpec2DObj.detectors`
(`pypeit/spec2dobj.py:491`) returns `[]` for an empty container, so the early
return is safe exactly as the doc claims.

**Substantive drift beyond line numbers** (these change what gets written in
Task 2, not just where):

1. **The QA API.** §B.1 calls `qa.init_qa_pool(1)` and `qa.flush_qa()`; those
   names do not exist anywhere in the package.  PR A shipped a `QAWriter` class
   instead.  Substitutions: `qa.init_qa_pool(1)` → **`qaWriter.init(ncpu=1)`**,
   `qa.flush_qa()` → **`qaWriter.flush()`**.  `QAWriter.init`'s own docstring
   already anticipates the fork case, so this is a clean substitution, not a
   workaround.  See Q2.
2. **`calib_one` gained `status_only` and `reload_only`** keyword arguments
   (`pypeit_steps.py:173-174`), which the §B.2 sketch of `calib_status_one`
   predates.  See Q3.
3. **`calib_one` already flushes QA** (`qaWriter.flush()` at
   `pypeit_steps.py:255`) and **already emits the "unsuccessful" warning**
   itself (`pypeit_steps.py:250-254`) — so stage 1 under `ncpu>1` will log that
   warning twice: once from the worker's replayed buffer and once from the
   parent's rebuilt-surviving-list loop.  That is what §B.2 specifies
   ("preserving the existing warning message verbatim"), so I will leave it,
   but it is worth knowing the duplication is expected and not a bug.
4. **`extract_det` does *not* flush QA** after `spec2DObj.gen_qa()`
   (`pypeit_steps.py:906`).  Harmless under the plan, because `_worker_init`
   puts the worker's `qaWriter` into serial mode (nothing is ever pending) and
   `_run_one` flushes anyway — but it is the reason `_run_one`'s flush must not
   be dropped.

**Structural facts re-confirmed against the current tree** (the design rests on
these):

- `pypeit.log` is a stock `logging.Logger` subclass obtained via
  `logging.getLogger("pypeit")`; `PypeItLogger.init` attaches a `StreamHandler`
  (`self.sh`) and optionally a `FileHandler` (`self.fh`) to both `pypeit` and
  the `py.warnings` logger (`pypeit/pkg/logger.py:224-256`).  Buffered worker
  logging is therefore exactly the handler swap §B.1 describes, and detaching
  `py.warnings` handlers in `_worker_init` is genuinely necessary — PR A did not
  change that.
- There is still **no multiprocessing anywhere in `pypeit/`**: `git grep -n
  'multiprocessing\|ProcessPool'` over `pypeit/` returns nothing.  PR A's only
  concurrency is `QAWriter`'s `ThreadPoolExecutor`.
- `spectrograph.select_detectors` returns a **fresh list**
  (`np.arange(1, ndet+1).tolist()` or a fresh list of mosaic tuples), and
  `reduce_exposure` does not hand it back to its caller — so replacing
  `detectors.remove(det)` with a rebuilt list is observationally safe (Q5).
- The `doc/running.rst` "Running on multiple CPUs" section exists
  (**`doc/running.rst:114-130`**) and currently states that `ncpu > 1` "only
  affects how the QA figures are written … Future versions will extend `ncpu`
  to the reduction".  That is the paragraph Task 6 has to rewrite.

**Environment.** `pypeit` conda env, Python 3.13.2, `threadpoolctl` **3.6.0**
imports cleanly.  `'fork' in multiprocessing.get_all_start_methods()` → **True**.
`os.cpu_count()` → **24**, so `nworkers(ncpu, nitem)` caps at 23 and the DEIMOS
`--ncpu 4` target is comfortably within the machine.  `threadpoolctl` is **not**
yet in `[project] dependencies` in `pyproject.toml` (it arrives transitively via
`scikit-learn>=1.2`); Task 2 adds it explicitly, as §B.1 requires.

**Baseline artifacts available for Task 7** in `Reports/`: `*_ncpu1.runmeta.json`
and `*_ncpu4.runmeta.json` for both `keck_deimos_600zd_m_6500` and
`shane_kast_blue_A` already exist from PR A's validation, so the after-B numbers
will have a like-for-like post-A comparison.

### 2026-10-10 (Task 2 — `pypeit/parallel.py`, its unit tests, and the `threadpoolctl` dependency)

Per the Q&A answers to Task 1: staying on **`speed_up_B`** (Q1a); using
`qaWriter.init(ncpu=1)` / `qaWriter.flush()` in place of the non-existent
`qa.init_qa_pool` / `qa.flush_qa` (Q2); `calib_status_one` will forward
`status_only`/`reload_only` (Q3, Task 3); `PypeIt.calib_all` gets a TODO (Q4,
Task 4); `detectors` confirmed un-aliased (Q5 — `select_detectors` returns a
fresh list and `reduce_exposure` never hands it back, so rebinding is safe).

**Written:**

- **`pypeit/parallel.py`** (294 lines) — implemented as the §B.1 sketch
  specifies: `THREAD_ENV_VARS`, `_FORK_PAYLOAD`, `nworkers`, `_RecordBuffer`,
  `_worker_init`, `_run_one`, `map_over_detectors`, `surviving_detectors`, with
  `__all__` and Numpy-style docstrings throughout.  Deviations from the sketch
  are limited to the two the Q&A authorised plus the required fork guard:
  - `qaWriter.init(ncpu=1)` in `_worker_init` and `qaWriter.flush()` in
    `_run_one`, instead of the `qa.*` functions (Q2).
  - The `'fork' not in mp.get_all_start_methods()` guard, which warns and falls
    through to the same serial comprehension the `ncpu<=1` branch uses.
  - Module docstring expanded to state the fork/Linux-only constraint, the
    module-level-function requirement, and *why* results are assembled in
    `detectors` order rather than completion order.
- **`pypeit/tests/test_parallel.py`** (9 tests) — all eight from §B.6, plus
  `test_fork_unavailable_falls_back_to_serial`, which monkeypatches
  `get_all_start_methods` so the fallback branch is not dead code on Linux.
- **`pyproject.toml`** — `threadpoolctl>=3.1` added to `[project] dependencies`,
  next to `scikit-learn` (which was supplying it transitively).

**Tests: 9 passed.**  Also ran `test_parallel.py` alongside `test_qa.py` and
`test_pypeitpar.py` — 53 passed — to confirm the new module does not perturb
the PR A machinery.

**Verified, not assumed:**

- `logging.makeLogRecord(record.__dict__)` preserves `created`, `filename`,
  `funcName` and `lineno`, which is exactly the set the two PypeIt formatters
  consume (`StreamFormatter.format` and `FileFormatter.base_fmt`,
  `pypeit/pkg/logger.py:99-157`).  Replayed worker lines therefore show the
  worker's real timestamp and source location, not the replay site.
- `logging.Logger.handle()` does not re-check the logger level, so the parent's
  *handler* levels apply on replay — the §B.1 claim that this "exactly
  reproduces the serial filtering" holds.
- `PypeItError` is a bare `Exception` subclass (`pypeit/pkg/exceptions.py:13`),
  so `pytest.raises(..., match='detector 2')` matches against the formatted
  message as intended.
- `ProcessPoolExecutor` forks lazily, on first `submit()`, not at construction —
  but both happen inside the `with` block while `_FORK_PAYLOAD` is set, so the
  "fresh executor per call, never a stale payload" property holds either way.

**Finding that may change the plan — the BLAS pinning cannot deliver
bit-identical `ncpu=1` vs `ncpu>1` output as specified.**  Written up in full as
**Q6**.  In short: §B.1 pins threads *only inside the workers*, so `ncpu=1` runs
OpenBLAS at 24 threads in the parent while `ncpu>1` runs it at 1 thread in each
worker, and those disagree for any thread-split BLAS reduction.  Measured
through the real helper, `ncpu=1` vs `ncpu=2` on a 2e6-element `np.dot`
differed by 7e-10 and 2.3e-09.  Several common operations are affected
(1500x1500 matmul 1.1e-13, 500x500 SVD 1.5e-13, 600x600 inverse 4.1e-14), while
others are not (3000x800 `M.T @ M` is identical).  Whether this actually bites
depends on whether the reduction path contains such an operation, which is
precisely what Task 7's DEIMOS identical-output diff will reveal; I have
implemented §B.1 unchanged and left the decision to Q6.

A direct consequence for the tests: **the §B.6 `_work` workload is not
thread-order sensitive at 200x200** (verified — identical bits at 1 and 24
threads), despite its sketch docstring describing it as "a reduction whose order
BLAS can change".  `test_serial_and_parallel_are_bit_identical` therefore
validates the mapping machinery (assembly order, payload routing, marshalling)
and *not* the pinning.  I left the workload alone — swapping in a genuinely
thread-sensitive one would make the test fail for the reason described above —
but rewrote both docstrings to say plainly what the test does and does not
prove, and pointed at `test_worker_threads_are_pinned` (which checks
`OMP_NUM_THREADS` and `threadpool_info()` inside the worker) as the actual
pinning check.  I also made that test non-vacuous: `_env` now forces a small
matmul so an unloaded BLAS is loaded and reported, the test asserts the worker
reported at least one runtime, and it skips rather than passing silently if the
environment has no threaded BLAS at all.

**Two smaller things flagged for decisions:** worker-raised `warnings.warn`
output is currently dropped entirely rather than buffered (**Q7**), and the
`pypeit` conda environment named in this prompt cannot collect most of
`pypeit/tests` because its astropy 7.0.2 calls the numpy-2-removed `np.in1d`
(**Q8**) — pre-existing, unrelated to PR B, and not a problem in `pypeit14b`.

No caller has been converted; `pypeit/parallel.py` is not yet imported anywhere
in the package.  `git status`: `pyproject.toml` modified, `pypeit/parallel.py`
and `pypeit/tests/test_parallel.py` new, nothing else touched.

### 2026-10-10 (Task 3 — the three detector-first adapters in `pypeit_steps.py`)

Per the Task 2 answers: Q6 → option (3), proceed as specified and let Task 7
decide empirically, **with the identical-output tolerance relaxed from
bit-identical to ~1e-6** (this changes the Task 6 test, which §B.6 specifies as
`np.testing.assert_array_equal`, and the Task 7 checks — noted now, acted on
there).  Q7 → worker `py.warnings` records stay dropped; only the PypeIt logger
has to survive, which it does.  Q8 → **`pypeit14b`** is the environment for
everything, including the PR A baselines; all commands below used it.

**Written** — three module-level adapters, inserted immediately after the
functions they wrap, as §B.2/§B.3/§B.4 specify.  `pypeit/pypeit_steps.py` grew
by 183 lines and `git diff` shows **zero deletions**: this task is purely
additive and no caller has been converted.

| New function | Line | Wraps | Logs | Returns |
|---|---:|---|---|---|
| `calib_status_one` | **260** | `calib_one` (170) | `Calibrating detector {det}` | `(success, failed_step)` |
| `process_one_det_bydet` | **454** | `process_one_det` (326) | `Reducing detector {det}` | `(sciImg, bkg_redux_sciimg)` |
| `extract_det_bydet` | **1027** | `extract_det` (853) | — (none; `extract_det` logs its own) | `(spec2DObj, sobjs)` |

All three take `det` first, are module-level (not closures), and carry full
Numpy-style docstrings documenting every parameter.  Each docstring states
explicitly *why* it must stay module-level — `fork` tolerates a closure but
`spawn`, the v2 goal, does not.

**Verified rather than eyeballed.**  A signature audit comparing each adapter
against the function it wraps found **no parameter present on one and not the
other, and no default mismatch, for any of the three pairs**.  Then a
mock-patched run through `parallel.map_over_detectors` confirmed the whole
det-first contract end to end:

- `calib_status_one` forwarded `det` to `calib_one`'s 4th positional and
  returned only the status pair: `{1: (True, None), 2: (True, None),
  3: (False, 'flats'), 4: (True, None)}`.  Feeding that to
  `parallel.surviving_detectors` returned **`[1, 2, 4]`** — detector 3 dropped,
  detector 4 *kept*, which is precisely the behaviour the
  `detectors.remove()`-during-iteration bug gets wrong today (it would have
  skipped 4 and left 3 in the list).
- `process_one_det_bydet` forwarded `det` to `process_one_det`'s 5th positional
  and returned the pair unchanged.
- `extract_det_bydet` was exercised through a **real 2-worker fork pool** with a
  per-detector `payload`.  Each worker received only its own slice
  (`IMG1`→det 1, `IMG2`→det 2), and although the log shows **detector 2
  finishing before detector 1**, the returned dict was still ordered
  `{1: ..., 2: ...}` — the order-independence property §B.1 relies on,
  demonstrated rather than asserted.

**Three small corrections to the coding-doc sketches**, all adopted:

1. §B.2 types the parameter as `run_state:dict=None`.  It is not a dict — it is
   a `RunPypeItState` instance (`calib_one` itself leaves it unannotated).  The
   adapter follows `calib_one` and leaves it unannotated; the docstring names
   the real class.
2. §B.2's `calib_status_one` sketch predates `status_only` and `reload_only`.
   Both are now forwarded, per Q3, so the adapter is a complete stand-in for
   `calib_one` and the deferred `PypeIt.calib_all` conversion (Q4) stays
   possible without touching it again.
3. §B.3's `process_one_det_bydet` sketch drops `sci_outfile`/`bkg_outfile`.
   Both are forwarded, for the same reason.

**On `extract_det_bydet`'s signature.**  Its five per-detector objects
(`sciImg`, `final_sky`, `sobjs_obj`, `calib_slits`, `bkg_redux_final_sky`) are
keyword arguments defaulting to None because the mapper supplies them through
`payload`, but four of them are *required positionals* on `extract_det`.  The
docstring says so outright, so nobody reads the None defaults as "optional" and
calls it without them.

**Tests:** `test_parallel.py` 9 passed; with `test_manual.py` (the only test
module that imports `pypeit_steps`) 12 passed.  `pypeit.pypeit_steps`,
`pypeit.exposure`, `pypeit.pypeit` and `pypeit.parallel` all import cleanly.

**Docstring style note:** `pypeit_steps.py` is otherwise uniformly Google-style
(`Args:`/`Returns:`), so these three functions are the only Numpy-style
docstrings in the file.  That follows this task's explicit instruction, the
§B.1 precedent and `CLAUDE.md`'s stated preference; Sphinx renders both via
napoleon, so the mix is cosmetic only.  Flagging it in case you would rather
they matched their neighbours.

### 2026-10-10 (Task 4 — stage 1 converted, the failed-detector bug fixed, stage 3 documented)

**Changed:** `pypeit/exposure.py` (+32/-9), `pypeit/pypeit.py` (+8).  Stages 2
and 4 are untouched and still serial; they are Task 5.

1. **`from pypeit import parallel`** added to `exposure.py` (line 17).
2. **Stage 1 converted** (`reduce_exposure`, now lines 428-460).  The loop is
   replaced by `parallel.map_over_detectors(pypeit_steps.calib_status_one, ...)`
   driven by `ncpu=par['rdx']['ncpu']`, followed by a *separate* loop that emits
   the unsuccessful-calibration warning (**message preserved verbatim**, only
   `caliBrate.failed_step` becoming the unpacked `failed_step`) and then
   `parallel.surviving_detectors` to rebuild the list.
3. **Early return added** when every detector fails: logs
   `'Calibrations were unsuccessful for every detector; nothing to reduce.'` and
   returns `spec2dobj.AllSpec2DObj(), specobjs.SpecObjs()`.  Re-verified the
   claim this depends on: `reduce_calibID` does nothing with the two return
   values before testing `len(this_spec2d.detectors) > 0`
   (**`pypeit/pypeit.py:431`**), and `AllSpec2DObj.detectors` returns `[]` for an
   empty container, so the existing "not successful for all the detectors"
   warning fires and nothing is written.
4. **Stage-3 comment added** to `findobj_on_exposure` (lines 205-210) recording
   the fan-out / `adjust_for_slitmask` barrier / fan-out structure and the
   deferral to v2, per §B.5.
5. **`PypeIt.calib_all` TODO added** (`pypeit/pypeit.py:235-243`) per **Q4**,
   naming the two reasons it was not converted: the loop keeps the returned
   `Calibrations` object and assigns `self.det` as a side effect.

**The bug, demonstrated rather than asserted.**  Replaying the coding doc's own
scenario — DEIMOS mosaics `[(1,5),(2,6),(3,7),(4,8)]` with `(2,6)` failing —
against the literal old loop:

```
OLD  calibrated: [(1,5), (2,6), (4,8)]        <- (3,7) never calibrated
OLD  survivors : [(1,5), (3,7), (4,8)]        <- yet (3,7) is kept
NEW  calibrated: all four
NEW  survivors : [(1,5), (3,7), (4,8)]
```

The old code skipped `(3,7)` entirely *and* left it in the list, so it would
have been passed to stages 2-4 and reduced against calibrations that were never
built.  The new code calibrates all four and drops only the one that failed.

**Smoke test — `shane_kast_blue`, both modes, end to end, cold runs into
separate directories.**  Both completed, exit 0, and wrote the identical set of
3 `spec1d` + 3 `spec2d` files:

| | wall | peak RSS |
|---|---:|---:|
| `--ncpu 1` | **80.7 s** | 1 498 MB |
| `--ncpu 2` | **84.0 s** | 1 446 MB |

The 3.3 s cost at `--ncpu 2` on a *single-detector* setup is the expected
fork + marshalling overhead plus the loss of QA threading inside the worker
(`_worker_init` forces `qaWriter.init(ncpu=1)` per **Q5**); Kast gains nothing
from this PR by design, and §B.7 only requires it not regress materially.  Peak
RSS did not grow, as expected with one worker.

**Products are bit-identical.**  Compared every `OPT_COUNTS`, `OPT_COUNTS_IVAR`,
`OPT_WAVE`, `BOX_COUNTS` and `TRACE_SPAT` array across all 3 exposures and all 5
objects: **worst relative difference 0.000e+00**.  This is the first real
evidence bearing on **Q6**: the calibration build ran in a thread-*pinned*
worker at `--ncpu 2` and in the *unpinned* parent at `--ncpu 1` and still
produced identical bits, i.e., the calibration path contains no thread-sensitive
BLAS reduction.  Encouraging, but not conclusive for the question — Kast is
single-detector and stages 2/4 are not yet converted.  DEIMOS in Task 7 remains
the real test.

**Code paths confirmed to be the ones intended.**  At `--ncpu 1` the run log
contains **zero** occurrences of `worker process` / `Started detector` /
`begin log for` — it is the literal serial comprehension, no pool, no fork.  At
`--ncpu 2` all five markers appear once per exposure, and the pool *is* entered
despite there being only one detector, which is exactly the CI-coverage
behaviour §B.1 asks for.

**Worker log replay verified against the run log file**, not just stdout:

```
INFO | 2026-10-10 11:54:20 | parallel.py:map_over_detectors:282 | ----- begin log for detector 1 -----
INFO | 2026-10-10 11:53:41 | pypeit_steps.py:calib_one:227      | Building/loading calibrations for detector 1
```

The replayed record carries the **worker's** timestamp (11:53:41) inside a
parent marker stamped 11:54:20, and the **worker's** source frame
(`pypeit_steps.py:calib_one:227`), not the replay site — so `record.created`,
`filename`, `funcName` and `lineno` all survived the pickle, and the block is
contiguous for the dashboard to parse.

**Quantified the Q7 behaviour now that it is observable.**  The `--ncpu 2` log
is 11 lines longer than the `--ncpu 1` log, not the +15 the markers alone would
add, because **exactly 4 `py.warnings` records were dropped** in the worker
(three `RuntimeWarning: invalid value encountered in ...` and one
`RuntimeWarning: Degrees of freedom <= 0 for slice`).  All **54** records from
the PypeIt logger itself survived, so the requirement you set in Q7 — "only that
the PypeIt logging is maintained" — holds.  Recording the number here so the
diff is explainable if anyone compares logs later.

**Tests:** 56 passed (`test_parallel`, `test_manual`, `test_qa`,
`test_pypeitpar`).  No new questions for the Q&A.

### 2026-10-10 (Task 5 — stages 2 and 4 converted; two latent pickling bugs fixed; Q4 measured; Q6 settled)

Stages 2 and 4 are converted and the whole `--ncpu 2` pipeline runs end to end.
Getting there required fixing two pre-existing crashes, and it produced a
definitive answer to **Q6** three tasks early.  Two decisions are now waiting on
you: **Q9** (which BLAS fix) and **Q10** (scope of the two bug fixes).

**Stage 2 — `process_exposure`.**  Loop replaced by
`parallel.map_over_detectors(pypeit_steps.process_one_det_bydet, ...)`; the two
dicts are assembled by comprehension over `detectors`, not over completion
order.  Also dropped the two now-dead `sciImg_dict = {}` /
`bkg_redux_sciimg_dict = {}` initialisations, which the comprehensions
immediately rebound.

**Stage 4 — `extract_exposure`.**  Per-detector inputs are sliced **in the
parent** into a `payload` dict (`sciImg`, `final_sky`, the
`all_specobjs_objfind` selection for that `detname`, `calib_slits[i]`,
`bkg_redux_final_sky`), with `detnames[det]` captured in the same pass so the
assembly loop does not need `sciImg_dict` again.  The `all_specobjs_extract.calibs`
assignment is preserved **verbatim**, including its overwrite-every-iteration
behaviour, with a comment recording that this is deliberate and out of scope.

#### Two latent crashes in core classes (see Q10)

Neither stage could return anything to the parent.  Both failures were
`RecursionError` raised inside `_ForkingPickler.loads` — i.e., in the **parent**
while unpickling, not in the worker — and both are the same pre-existing bug
that has nothing to do with parallelism:

> Pickle constructs an instance *without* calling `__init__`, then probes it for
> `__setstate__`.  That lookup lands in a `__getattr__` which immediately reads
> an attribute `__init__` would have set.  The attribute is missing, so
> `__getattr__` is re-entered to fetch it, and recurses until the stack blows.

- **`pypeit/images/bitmaskarray.py:114`** — `BitMaskArray.__getattr__` reads
  `self.lower_keys`.  Blocked **stage 2** (every `PypeItImage` carries a
  bitmask) *and* stage 4.  Minimal reproduction with no PR B code involved:
  `pickle.loads(pickle.dumps(ImageBitMaskArray((20,10))))`.  Traced it to
  confirm: the first probe is `__getattr__('__setstate__')` on an instance whose
  `__dict__` is literally empty, followed by `__getattr__('lower_keys')` forever.
- **`pypeit/specobjs.py:776`** — `SpecObjs.__getattr__` reads `self.specobjs`.
  Blocked **stage 4**.

Both fixed identically: read the attribute from `self.__dict__` and raise
`AttributeError` if it is absent, since an instance in that state has nothing to
look up.  Verified the round trip preserves the mask array, its dtype,
`lower_keys`, and per-flag access (`bpm`, `cr` counts match), that
`copy.deepcopy` now works too, and that a genuinely missing attribute still
raises `AttributeError`.

#### Q6 is settled, and the tolerance does not save it

Four cold `shane_kast_blue` runs, compared HDU by HDU:

| comparison | result |
|---|---|
| Task-4 serial vs Task-5 serial | **bit-identical** — my stage-2/4 edit did not change serial behaviour |
| Task-4 `--ncpu 2` vs Task-5 serial | **bit-identical** |
| Task-5 serial vs Task-5 `--ncpu 2` | differs: `SKYMODEL`, `OBJMODEL`, `IVARMODEL`, `MED_CHIS`, `STD_CHIS` |
| serial *unpinned* vs serial *BLAS-pinned* | differs: **the same five HDUs** |
| serial *BLAS-pinned* vs `--ncpu 2` | **BIT-IDENTICAL** |

The last two rows are the whole story.  **The parallelism is correct** — it
reproduces a serial run bit-for-bit whenever the thread count matches, which
validates the payload slicing, the result assembly and the marshalling.  **The
entire difference is BLAS thread count**, and it is not a 1e-9 effect:
`SKYMODEL` differs on 621 756 of 716 800 pixels (87%), `OBJMODEL` on 81%, while
every *input* HDU (`SCIIMG`, `IVARRAW`, `TILTS`, `WAVEIMG`, `BPMMASK`) is
bit-identical.  The amplifier is `local_skysub_extract`'s iterative b-spline fit
with sigma-clipping: a last-ulp change flips a marginal rejection, which changes
the fit, which changes the next rejection.  So the 1e-6 tolerance you set in Q6
cannot be met by the current design — the differences are chaotic, not
round-off.  **Q9** lays out the three fixes with measured costs; pinning the
serial path is already *demonstrated* to work (row 5), and costs **+5.7% on
Kast (77.5 s -> 81.9 s)**.  I have not applied it, because you reserved that
tradeoff.

#### Q4 — measured, and the answer is "no change needed"

Now that `Spec2DObj` pickles at all, the marshalling cost is measurable:

| | size | dumps | loads | fork+pipe |
|---|---:|---:|---:|---:|
| real Kast `Spec2DObj` (2048x350) | 26.1 MB | 0.022 s | 0.023 s | — |
| DEIMOS mosaic payload (8 arrays, 8217x2066) | **615.2 MB** | 0.539 s | 0.549 s | **1.84 s** |

Against a stage-4 wall of 738 s per mosaic (2 952 s / 4), that is **0.25% of the
stage — well under the 5% threshold**, so stage-4 workers keep returning the
`Spec2DObj` and the `to_file(Intermediate/...)` fallback is not needed.  Worth
noting the coding doc's "~250 MB per DEIMOS mosaic" estimate is low: the real
figure is **615 MB**, from the actual HDU geometry of the 600ZD_M_6500 spec2d
(five `float32` and two `float64` 8217x2066 frames plus an `int16` mask).  Even
at 2.5x the assumed size the cost is negligible, so the conclusion is robust.

#### Validation

- **Full unit suite: 865 passed, 3 failed.**  The 3 failures are all in
  `test_pkgdata.py` (`test_fetch_github_files`, `test_github_contents`,
  `test_cache_to_pkg`) and are **pre-existing**: confirmed by stashing every
  change and re-running, where the same 3 fail on a clean tree.  They are GitHub
  fetch failures, and `test_pkgdata.py` references none of the files I touched.
  (That module also deletes the tracked fixture `pypeit/data/tests/cache_test.txt`
  as a side effect of running; restored it with `git checkout`.)
- **Smoke test:** `--ncpu 1` exit 0 in **77.5 s** (peak RSS 1 500 MB);
  `--ncpu 2` exit 0 in **84.2 s** (peak RSS 1 457 MB).  Both wrote the full set
  of 3 spec1d + 3 spec2d products.  Kast is single-detector, so `--ncpu 2` has
  nothing to parallelise and pays only fork + marshalling + the loss of QA
  threading inside the worker.

`git diff --stat`: 6 files, +284/-39, plus the two new files from Task 2.

### 2026-10-10 (Task 6 — serial pinning applied per Q9(1); tests, dev-suite regression, docs, changelog)

Per **Q9 -> option (1)** and **Q10 -> approved as-is** (so the two `__getattr__`
fixes stay in this PR and are written up under **Bug Fixes**).

#### Q9(1) applied and verified

`map_over_detectors` now runs its serial loop inside `threadpool_limits(1)`,
via a small `_serial()` closure shared by the `ncpu<=1` branch and the
no-`fork` fallback, so every path that reduces a detector sees exactly one BLAS
thread.  The module and function docstrings were rewritten to state the new
contract and, explicitly, the trade-off being made.

**The important caveat I had to check:** the Task-5 evidence used *environment
variables* (`OMP_NUM_THREADS=1` etc.) set before the process started, whereas
the fix uses `threadpoolctl.threadpool_limits`, which retunes an
already-loaded runtime.  Those are different mechanisms and I did not assume
they were equivalent.  Re-ran both modes cold and compared:

```
spec2d: 33 image HDUs compared, 0 differ
spec1d: 30 arrays compared, 0 differ
VERDICT: BIT-IDENTICAL
```

So `--ncpu 1` and `--ncpu 2` now agree bit for bit, and §B.6's
`assert_array_equal` is the right assertion after all — the 1e-6 relaxation
from Q6 is no longer needed, and I did **not** weaken the test.

Wall-clock on Kast, for the record: **86.3 s** at `--ncpu 1` and **89.2 s** at
`--ncpu 2`, against **77.5 s** for the old unpinned serial run.  Note the
pinning cost measures noisier than the Task-5 single sample suggested — 81.9 s
there vs 86.3 s here for the same configuration, so the honest statement is
"several percent, 6-11% on this setup" rather than a precise figure.  Kast is
single-detector, so it pays the pinning cost with no parallel speedup to offset
it; **DEIMOS in Task 7 is the measurement that matters.**

#### Tests

- **`test_run_pypeit_ncpu`** added to `pypeit/tests/test_runpypeit.py`: runs
  `Setup.main` once, then reduces into `serial/` and `parallel/` subdirectories
  with `--ncpu 1` and `--ncpu 2`, and compares `OPT_COUNTS`,
  `OPT_COUNTS_IVAR`, `OPT_WAVE`, `BOX_COUNTS`, `TRACE_SPAT` with
  `np.testing.assert_array_equal` plus `WAVE_RMS`.  Guarded by
  `PYPEIT_SKIP_SLOW_TESTS` and marked `remote_data`, matching the existing
  `test_run_pypeit`.  Added the missing `import os`.
- **`pytest pypeit/tests`: 857 passed, 0 failed** (excluding `test_pkgdata.py`,
  whose 3 GitHub-fetch failures were shown in Task 5 to be pre-existing on a
  clean tree).

#### Dev-suite multi-detector regression (Q6 of the coding doc)

- **`test_scripts/pypeit_tests.py`** — new `PypeItReduceNcpuTest`, a subclass of
  `PypeItReduceTest` that appends `-r <rdxdir>/ncpu --ncpu N`.  The separate
  output directory is the crux: `rdxdir` is shared by every test on a setup, so
  without it the parallel run would simply overwrite the serial products it is
  meant to be compared against.  Verified the generated command lines:

  ```
  serial : run_pypeit <file> -o                                   (cwd = <rdxdir>)
  ncpu   : run_pypeit <file> -o -r <rdxdir>/ncpu --ncpu 2
  ```

- **`test_scripts/test_setups.py`** — new `_ncpu_reduce` dict registering
  `keck_lris_blue/multi_600_4000_d560` with `ncpu=2`, documented in the module
  header, and registered in `all_tests` **after** the default reduce entry, with
  a comment saying why the order matters.  `all_tests` goes 18 -> 19 entries.
- **`vet_tests/test_ncpu.py`** — three tests: every spec1d array identical
  (and the objects on the same detectors in the same order), every spec2d image
  extension identical, and both reductions covering **both** detectors — the
  last being a direct regression for the `detectors.remove()` bug, which on a
  2-detector setup would leave one detector missing.
- Confirmed the 10 failures in `test_scripts/test_pypeit_test.py` are
  **pre-existing** by stashing my changes and re-running on a clean tree.

#### Docs

- **`doc/running.rst`** — the "Running on multiple CPUs" section rewritten from
  PR A's "``ncpu > 1`` only affects how the QA figures are written" into four
  subsections: what is parallelized (the three stages) and what is not (object
  finding + global sky, with the barrier explained); how many CPUs to ask for
  (`min(ncpu, ndet, ncores-1)`, and that **memory, not cores, is the real
  limit** — `--ncpu 2` advised for DEIMOS-class data on <=32 GB); the Linux-only
  `fork` caveat as a `note::`; and a **Reproducibility** subsection stating that
  products do not depend on `ncpu`, why both paths pin threads, and what that
  costs.
- **`./update_docs` run; `doc/sphinx_warnings.out` is empty (0 lines).**
  `doc/api/pypeit.parallel.rst` was generated and `git add`ed along with the
  `doc/api/pypeit.rst` entry.
- `doc/include/dependencies_table.rst` regenerated, picking up
  `threadpoolctl>=3.1`.
- **Reverted two unrelated regenerations** that `update_docs` also touched:
  `doc/include/keck_deimos.sorted.rst` (pure sort tie-break churn — rows with
  *identical* mjd swapping places, nothing to do with this PR) and
  `doc/help/run_pypeit.rst` (a bare version-string bump).  Neither is in §5's
  doc checklist for PR B, and both would be noise in the diff.

#### Changelog (`doc/releases/2.1.0dev.rst`, +43 lines)

Bullets added under all five required headings: **Functionality/Performance**
(the concurrency, plus a second bullet stating the thread-pinning trade-off
plainly, since it is a user-visible slowdown at the default setting),
**Dependency Changes** (`threadpoolctl`), **Under-the-hood** (`pypeit.parallel`),
**Testing** (both the unit test and the dev-suite regression), and **Bug Fixes**
(the `detectors.remove()` bug *and* the two `RecursionError` fixes from Q10).

#### Still outstanding

The dev-suite regression is **registered but not yet executed** — running it
needs the `reduce` phase for `keck_lris_blue/multi_600_4000_d560`, which is part
of Task 7's `./pypeit_test reduce -i ...` run.  Task 7 will be its first real
exercise, and is also where the DEIMOS numbers, the peak-RSS measurements and
the determinism check land.

### 2026-10-10 (Task 7 — validation and re-profiling; the identical-output requirement turns out to be unachievable)

Numbers written up in `Reports/speed_up_results.md` ("After workstream B").
**Two decisions are blocked on you: Q11 (what to do about the non-determinism
and the dev-suite test) and Q12 (whether to re-measure the baseline cleanly).**
Two items of Task 7 are deliberately left undone pending those answers; they are
listed at the end.

#### Performance

| Setup | `--ncpu 1` | `--ncpu 4` | speedup |
|---|---:|---:|---:|
| keck_deimos 600ZD_M_6500 | 10 515.5 s | **4 123.6 s** | **2.55x** |
| keck_deimos, repeat | — | 4 244.5 s | 2.48x |
| keck_lris_blue slitmask (2 det) | 389.0 s | 250.7 s (`--ncpu 2`) | 1.55x |
| shane_kast_blue (1 det) | 109.7 s | 102.5 s (`--ncpu 2`) | 1.07x |

Implied p = 0.81 against the 0.77 the design assumed.  Scatter at `--ncpu 4` is
120.9 s (2.9%).  **Kast does not regress** — it is 6.6% *faster* at `--ncpu 2`,
since PR A's QA threading outweighs the fork cost, which comfortably satisfies
§B.7's "must not regress by more than a second or two".

I recalibrated the §B.7 target before comparing: the doc assumed an 11 700 s
post-A baseline, but `speed_up_results.md` records 10 056.9 s, so the equivalent
of its 1.9-2.3x window is 4 373-5 293 s rather than 5 200-6 200 s.  The measured
4 123.6 s is at or just past the optimistic end.

**I do not fully trust the 2.55x, and said so in Q12.**  §B.7 instructs that the
baseline be cProfile-instrumented and the parallel runs not, which makes the
ratio apples-to-oranges.  The tell is that 2.55x *exceeds* the 2.37x Amdahl
ideal for p = 0.77 — code does not beat Amdahl; inflated baselines do.  Kast
pays ~27% for the same instrumentation (86.3 s vs 109.7 s), so if DEIMOS is
comparable the honest figure is nearer 2.0x.  One clean ~2.9 h run settles it.

#### Memory

Peak **PSS** summed over the process tree: **11.70 GB** at `--ncpu 1`, **39.38
GB** at `--ncpu 4` — **3.4x for 4 workers**, sub-linear thanks to copy-on-write
sharing, and better than the roughly-linear scaling the design assumed.  I wrote
`scripts/run_with_rss.py` for this rather than using `/usr/bin/time -v`, whose
`ru_maxrss` is the maximum of any *single* process and would have under-reported
a 4-worker run by about 4x — i.e., it would have silently measured the wrong
thing.  Summing plain RSS instead over-reports (45.03 GB) by double-counting
shared pages; both are recorded.

#### Cost of the Q9(1) pinning, now measured on real data

`ncpu = 1` is **+11.4%** on Kast (77.5 -> 86.3 s, raw) but only **+4.6%** on
DEIMOS (10 056.9 -> 10 515.5 s, both instrumented).  DEIMOS pays proportionally
less because more of its time is in b-spline fits that were never
BLAS-threaded.  The DEIMOS number also carries some `develop`-merge drift.

#### The finding that changes the plan

The DEIMOS `--ncpu 1` vs `--ncpu 4` identical-output check **failed**
(2 039/2 624 arrays).  Rather than conclude the parallelism was broken, I ran
the control — and reproduced the whole thing on
`keck_lris_blue/multi_600_4000_slitmask`, a 2-detector setup cheap enough to
iterate in minutes instead of three hours:

| comparison | differing |
|---|---|
| DEIMOS `--ncpu 4` vs `--ncpu 4` **repeat** | **2 033 / 2 624** |
| LRIS serial vs `--ncpu 2` | 637 / 842 |
| **LRIS serial vs a second identical serial run** | **637 / 842** |
| LRIS serial vs serial, thread env vars = 1 *and* `PYTHONHASHSEED=0` | **637 / 842** |
| shane_kast_blue serial vs `--ncpu 2` (three times) | **0 / 328** |

Two runs of the *same command* disagree by the same amount as a
serial-versus-parallel pair.  It is not thread count and not hash randomisation.
**PypeIt is not reproducible run-to-run on slitmask data, and PR B is not the
cause** — it is simply the first thing to have checked.

Localised by comparing the calibration files of two identical serial runs:
`Arc` and `Tiltimg` are bit-identical, `Flat` differs in 30/42 arrays and
`Tilts` in 3/14.  Inside the flat, `PIXELFLAT_RAW` is bit-identical while every
`*_BSPLINE` `bkpt_full`/`coeff`/`icoeff` differs and `PIXELFLAT_MODEL` differs
by 4.4e-3 — i.e., **the B-spline fit returns different coefficients from
identical input**.  That is `pypeit/core/bspline/`, explicitly out of scope
here.  Downstream, the sigma-clipped sky fit amplifies it into a visibly
different `SKYMODEL`, which is why the spec2d diffs look so dramatic.

Two hypotheses I tested and **ruled out** before landing on this, recorded so
nobody repeats the work: unseeded RNG in the reduction path (there is none —
the only `np.random` uses are MMIRS ramp simulation, display colours and
docstrings), and sklearn's randomized SVD in the trace PCA (`PCA()` at
`pca.py:92` has `n_components=None` so `auto` resolves to `full`, and at
realistic trace shapes `pca.py:120` resolves to `covariance_eigh` — both
deterministic LAPACK paths, verified directly).

#### Tooling written (both found real problems)

- **`scripts/compare_products.py`** — element-wise spec1d/spec2d comparison
  reporting, per array, the differing-element count and the maximum absolute and
  relative difference.  Its first version compared only image HDUs, so `spec1d`
  reported "0 arrays" and was **silently not being checked at all** — spec1d
  content lives in binary tables.  After the fix the Kast comparison went from
  33 to 328 arrays.  Validated against both a known-identical and a
  known-different pair before being trusted.
- **`scripts/run_with_rss.py`** — the PSS tree sampler described above.

#### Also corrected during this task

`keck_lris_blue/multi_600_4000_d560`, the setup coding-doc Q6 named for the
dev-suite regression, pins **`detnum = 2`** in its `.pypeit` file — a *single*
detector.  It would not have exercised multi-detector concurrency at all, and
the 2-detector assertion in `vet_tests/test_ncpu.py` would have failed.  I
surveyed every dev-suite setup for ones that are genuinely multi-detector,
unrestricted and cheap, and switched to **`multi_600_4000_slitmask`** (2
detectors, 3 frames, no `detnum`), updating the registration, the vet test and
the changelog.  It then served as the fast reproducer above.

#### Deliberately not done, pending Q11 / Q12

- **`vet_tests/test_ncpu.py` is left asserting exact equality**, which means it
  **will fail** on LRIS.  Rewriting it to assert scientific equivalence is my
  recommendation in Q11, but it changes what the PR claims, so I have not done
  it unilaterally.
- **The full `./pypeit_test reduce -i shane_kast_blue keck_deimos keck_lris_blue`
  plus `pytest vet_tests` run.**  The DEIMOS reduce would repeat ~3 h of what
  this task already measured directly, and the vet run would fail on the above
  test; both are better done once Q11 is settled.

### 2026-10-11 (Task 8 — upstream issue filed; and a correction: the B-spline was innocent)

**Filed: pypeit/PypeIt#2208** — "Reductions are not reproducible: edge
centroiding (EDGE_CEN) differs by 1 float32 ULP from bit-identical input",
created as `profxj`.

**Before writing anything public I tried to verify the claim, and it did not
survive.**  In Task 7 I wrote that "the B-spline fit returns different
coefficients from identical input".  That was an *inference* from the FITS
products — `PIXELFLAT_RAW` was bit-identical while the stored `*_BSPLINE`
coefficients differed — and it was wrong, because the intermediates between
those two things are not saved, so identical `PIXELFLAT_RAW` never implied
identical *bspline input*.  Testing it directly:
`fitting.iterative_bspline_fit` on fixed synthetic data (60k points, injected
outliers, rejection active) returns **byte-identical** coefficients, breakpoints
and fit values across repeats within one process *and* across separate
processes.  **The B-spline is deterministic.**  I have corrected
`Reports/speed_up_results.md` accordingly; the Task 7 log entry above is left as
written, with this entry as the correction.

**The real origin, found by walking the calibration chain in dependency order.**
My Task 7 sweep had missed `Edges` and `Slits` entirely, because they are stored
as `.fits.gz` and my glob was `*.fits`:

```
Arc_A_0_DET01.fits         15 arrays  identical
Tiltimg_A_0_DET01.fits     15 arrays  identical
WaveCalib_A_0_DET01.fits  173 arrays  identical
Edges_A_0_DET01.fits.gz    39 arrays   13 DIFFER   <- origin
Slits_A_0_DET01.fits.gz    32 arrays    7 DIFFER
Flat_A_0_DET01.fits        42 arrays   30 DIFFER
Tilts_A_0_DET01.fits       14 arrays    3 DIFFER
```

Inside `Edges` the separation is clean: `TRACE_IMAGE`, `TRACE_FULLMASK`,
`TRACEBPM`, `SOBELSIG`, `TRACEID`, `MASKDEF_ID` and `EDGE_MSK` are **all
bit-identical** — the same traces are found and masked the same way — and only
`EDGE_CEN`, the measured centroid, changes.  Across five runs (2 plain serial,
2 fully thread-pinned with `PYTHONHASHSEED=0`, 1 at `--ncpu 2`), **all ten
pairwise comparisons** show `SOBELSIG` identical and `EDGE_CEN` differing by
**exactly 1.0 float32 ULP (6.104e-05)**.  So: identical input, non-reproducible
output, insensitive to threads and hash seed.

**Not a duplicate, and it corrects the record.**  I searched existing issues
first and found **#1986 "Non identical extractions"**, filed by you, where the
conclusion was "this has been linked to the C-code for the bspline" and
@tbowers7 asked for a rigorous retest now that the C code is gone.  Task 7 *is*
that retest.  #2208 therefore does three things: reports that the problem
persists, identifies the actual origin, and explicitly exonerates the B-spline,
with the five ruled-out hypotheses (B-spline, BLAS threads, hash seed, unseeded
RNG, sklearn randomized SVD) listed so nobody repeats the work.  It cross-refs
#1986 and cc's @tbowers7.

The issue includes a self-contained reproducer: two `run_pypeit` invocations on
`keck_lris_blue/multi_600_4000_slitmask` plus a ten-line numpy check on
`EDGE_CEN`.  It also carries the downstream numbers (LRIS 637/842 arrays
differing between two identical serial runs; DEIMOS 2 033/2 624 between two
`--ncpu 4` runs; `SKYMODEL` on 7.65M of 17.0M pixels; `OPT_COUNTS` ~100% at low
S/N) and notes that `shane_kast_blue` is unaffected (0/328), which is why this
has been easy to miss.

One hypothesis offered for whoever picks it up, flagged as a hypothesis rather
than a finding: numpy SIMD paths can give different results depending on array
memory alignment, and allocator behaviour varies between processes; a weighted
moment sum over a slice whose base-pointer alignment varies would produce
exactly this signature.

**Note on the task instruction:** it says "Use Opus 5.5".  I am Opus 5 and
cannot switch models from inside the session, so the work was done on Opus 5 —
flagging it rather than silently ignoring it.

**Still open:** **Q12** (whether to spend ~2.9 h on an uninstrumented `--ncpu 1`
DEIMOS baseline so the headline speedup is apples-to-apples) remains unanswered,
and the two Task 7 follow-ups that depend on **Q11(1)** are still outstanding —
rewriting `vet_tests/test_ncpu.py` to assert scientific equivalence rather than
exact equality, and the full `./pypeit_test reduce` + `pytest vet_tests` run.

### 2026-10-11 (Task 9 — Q11(1) implemented, PRs opened; clean DEIMOS baseline running)

**Q12 answered "yes, run it"**, so the uninstrumented `--ncpu 1` DEIMOS baseline
was launched first and is still running (started 07:42, ETA ~10:35).  Everything
below was done while it ran.  **The speedup table in both
`speed_up_results.md` and PR #2209 still carries the inflated,
cProfile-instrumented ratio and a note saying so; I will update both with the
honest figure when the run lands.**

#### Q11(1) implemented: `vet_tests/test_ncpu.py` rewritten

Your answer was option (1) — keep PR B, assert scientific equivalence, keep the
structural checks exact.  Implementing it honestly took some measurement,
because my first attempt failed: I picked a 1e-3 relative tolerance out of the
air and the test failed immediately, with the Issue-2208 noise reaching **26%**
at the 95th percentile on `OBJMODEL`.  So I measured the noise floor instead of
guessing again.

The flux-like arrays turn out to be **untestable numerically**, and by a wide
margin.  Between two *identical serial* runs:

- `OPT_COUNTS` differs by up to **350 sigma** (serial vs `--ncpu 2`: 407 sigma)
- `OBJMODEL` 95th-percentile relative difference 9.4e-02 (parallel: 2.6e-01)

Any threshold loose enough to pass that would be far too loose to catch a real
bug, so "agreement well within the noise" — the phrasing in option (1) — is not
achievable for them.  But other quantities are perfectly robust, and the
serial/parallel difference in them is *the same size as the serial/serial
difference*:

| quantity | serial vs serial | serial vs `ncpu 2` |
|---|---|---|
| `spec1d` `OPT_WAVE` | 2.103e-05 | **2.103e-05** |
| `spec2d` `WAVEIMG` | 1.114e-10 | 1.318e-10 |
| `spec2d` `TILTS` | 2.129e-05 | 1.702e-05 |
| `spec2d` `SCIIMG` | 3.687e-07 | 1.192e-07 |

So the test now asserts: structure exactly (both detectors present, same objects
on the same detectors in the same order, same extensions); the robust
quantities to 1e-3, which is clear of the ~2e-5 noise while a misrouted payload
would show at order unity; and the exact flux comparison as a **`xfail`** keyed
to #2208, so it flips to passing on its own once that is fixed.  The measured
table is in the module docstring, so the thresholds are justified in place
rather than being magic numbers.

**I then checked the test can actually fail**, since one that only ever passes
proves nothing.  Two injected faults: reordering the detector extensions (caught
by `test_ncpu_same_structure`), and the subtler one of swapping the two
detectors' `WAVEIMG` payloads while leaving the extension *names* correct —
exactly a misrouted-payload bug — which `test_ncpu_spec2d_geometry_agrees`
catches at relative difference **1.0**, three orders above the threshold.
4 passed, 1 xfailed on the real products.

Docs corrected to match: `doc/running.rst` and the release note no longer claim
unqualified bit-identical output, and `running.rst` gained a warning that
multi-slit data is non-reproducible run-to-run independently of `ncpu`, so
nobody misreads #2208 as a parallelism defect.

#### PRs opened

- **pypeit/PypeIt#2209** -> `develop`, from `speed_up_B`.  14 files, +954/-55, in
  five logical commits (`pypeit.parallel`; the adapters; the unpickling fix; the
  three stage conversions plus the `detectors.remove()` fix; tests and docs).
  `speed_up_B` was exactly at `origin/develop`, so the branch stack note about
  targeting `speed_up_qa` is moot — PR A is already in `develop` via #2198 and
  `develop` is the correct base.
- **pypeit/PypeIt-development-suite#435** -> `develop`, from `dev_speed_up`.  The
  `--ncpu` dev-suite regression, the two tooling scripts, and the measurements.

Both descriptions lead with the result, state the thread-pinning trade-off
explicitly as a user-visible cost that is reversible in one place, and carry the
#2208 evidence table so a reviewer is not left wondering why the
identical-output check fails.  #2209 cross-references #2208 and #1986.

Verified before pushing: 32 targeted unit tests green on the committed state,
and the full `pytest pypeit/tests` run from Task 6 (857 passed) still stands.

#### Still outstanding

1. The clean DEIMOS baseline, and the resulting correction to the speedup number
   in `speed_up_results.md` and in #2209.
2. The full `./pypeit_test reduce -i shane_kast_blue keck_deimos keck_lris_blue`
   plus `pytest vet_tests` run.  The DEIMOS leg duplicates ~3 h of what Task 7
   measured directly; worth doing before merge, but it should wait until the
   baseline run finishes rather than competing with it for the machine.

### 2026-10-11 (Q12 follow-up — clean DEIMOS baseline; the headline is 2.24x, not 2.55x)

The uninstrumented `--ncpu 1` DEIMOS baseline finished: **9 360.0 s** (peak PSS
11.64 GB, rc 0).  The same reduction *under cProfile* took 10 515.5 s, so the
instrumentation costs **+1 155 s (+12.3%)** — which is what was inflating the
ratio.

| | `--ncpu 1` | `--ncpu 4` | speedup |
|---|---:|---:|---:|
| run A | 9 360.0 s | 4 123.6 s | **2.27x** |
| run B | 9 360.0 s | 4 244.5 s | 2.21x |
| **mean** | **9 360.0 s** | **4 184.0 s** | **2.24x** |

**Corrected: 2.24x, not the 2.55x reported in Task 7.**  Two things now line up
that did not before.  The implied parallel fraction is **p = 0.737**, close to
the 0.77 the design assumed rather than the 0.81 the inflated baseline implied;
and 2.24x sits *just below* the 2.37x Amdahl ideal for that p instead of above
it, which is the physically sensible ordering.  A measurement that beats Amdahl
is a measurement that is wrong, and that was the tell I flagged in Q12 — worth
having chased, because it was a 14% overstatement of the headline result.

It also lands **inside the §B.7 predicted band of 1.9-2.3x**, which the inflated
number did not.

Updated `speed_up_results.md` (table, the implied-p discussion, the verdict) and
the #2209 description.  Both now say explicitly that the DEIMOS figures are
uninstrumented on both sides, and the results doc keeps a short note recording
that an earlier version reported 2.55x and why that was wrong, so the change is
not silent.

One mechanical note: `gh pr edit --body-file` reported success while leaving the
body untouched — a GraphQL projects-deprecation warning masked the failure.  I
caught it by reading the body back, and used `gh api -X PATCH` instead.  Worth
verifying rather than trusting the exit status if you use `gh pr edit` elsewhere.

The pinning-cost figures are unchanged and still correct: that comparison is
instrumented-against-instrumented (post-A 10 056.9 s vs 10 515.5 s, +4.6%), so
the clean baseline does not affect it.

**Remaining:** the full `./pypeit_test reduce -i shane_kast_blue keck_deimos
keck_lris_blue` plus `pytest vet_tests` run, which is the last item of Task 7's
checklist.  Nothing blocks it now that the machine is free.
