# Issue #2208 -- Edge trace irreproducibility

## Goals

We wish to address Issue #2208 for PypeIt

## Context

See the issue description in the PypeIt repository: `https://github.com/pypeit/PypeIt/issues/2208`

Also see some discussion on this in `https://github.com/pypeit/PypeIt/pull/2209`

## Running

If you need to run Python, use the `pypeit14` environment on my laptop
or `pypeit14b` on my workstation

## Prompts

### Plan

1. Read the URLs described in the Context.  Then discuss with me a plan
to examine this Issue in depth on my workstation using the DevSuite.
Put your questions in the Q&A/Plan section below.
Use Opus 5.5. Log your work.

2. I answered your questions in Q&A/Plan below.  Read them and discuss further with me whether we can move ahead.    If you have no further questions, 
generate a series of prompts in the Testing section below to perform the work
and I will have you execute them.
Use Opus 5.5. Log your work.

### Testing

These rules apply to every prompt below:

- **Env and repos:** env `pypeit14b`, called by absolute path
  (`/home/xavier/miniconda3/envs/pypeit14b/bin/...`).  PypeIt is at
  `/mnt/tank/Astronomy/PypeIt/PypeIt` (branch `edge_2208`); the dev suite is
  this repo (branch `dev_edge_2208`).
- **No git:** run no git commands.  When a step reaches a natural commit point,
  say so in the log and list the files to commit.
- **Where things go:** scripts go in `pypeitdev/issue_2208/` (committed).
  Reduction outputs go in `/mnt/tank/Astronomy/PypeIt/REDUX/issue_2208/<name>/`,
  never in the repo.
- **Thread pinning:** unless a step says otherwise, export
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` for every
  reduction and harness run, so that thread count is never a confounder (#2208
  showed it is not the cause).
- **Long runs:** anything over ~10 min is launched detached
  (`setsid nohup ... &` or a named `tmux` session), with its log on disk, so
  the user can leave the Claude session and log out.  Report the command,
  session name and log path.
- **Logging:** log each prompt in Logs as `### <date> — Testing #N (summary)`,
  with exact commands, timings and numbers.  Put questions in Q&A/Testing.

1. **Baseline and reproduce.**
    - Record `np.show_config()`, `threadpoolctl.threadpool_info()`, the CPU
      model and flags (`lscpu`), and the glibc version.
    - Write `pypeitdev/issue_2208/compare_calibs.py`.  Given two output
      directories, it walks every FITS file in `Calibrations/` (and, when
      present, `Science/`), and for every extension/array reports: identical
      or not, number of differing elements, max |diff|, and max float32/float64
      ULP distance.  It also lists files in dependency order (Arc, Tiltimg,
      WaveCalib, Edges, Slits, Flat, Tilts, ...).  It works on `.fits.gz`.
    - Run `keck_lris_blue/multi_600_4000_slitmask` **calibrations only**
      (`run_pypeit -c -o`) cold, twice (`base_a`, `base_b`), and time each run.
      Compare them.  For `EDGE_CEN` on both detectors, report how many elements
      differ, in which traces (columns), and at which rows, and whether those
      traces are left or right edges, inserted, or slitmask-matched.
    - Cheap discriminator: two more cold `-c` runs with
      `OPENBLAS_CORETYPE=Haswell` (`hsw_a`, `hsw_b`; first check in the log
      that the core type actually took effect).  Compare `hsw_a` with `hsw_b`,
      and `hsw_a` with `base_a`.  Record the result; do not draw conclusions
      beyond it yet.

2. **Stepwise harness; find the first diverging step.**
    - Write `pypeitdev/issue_2208/edge_harness.py`.  It rebuilds the
      `EdgeTraceSet` for one detector of the LRIS setup from **saved inputs**
      (the `TraceImage`/BPM saved in an Edges file from Testing #1, or rebuilt
      from the calibration frames if that is not faithful).  It needs the same
      spectrograph, par and slitmask information as `run_pypeit`.  It then runs
      the `auto_trace` sequence step by step (wrap each method `auto_trace`
      calls, including the nested ones such as `fit_refine` inside `sync`).
      After every call it writes a SHA-256 of `edge_cen`, `edge_err`,
      `edge_msk`, `edge_fit`, and the PCA arrays to a JSON trace.
    - First check that the harness reproduces `run_pypeit`: its final
      `EDGE_CEN` must be bit-identical to that of one of the Testing #1 runs, or
      it must show the same kind of run-to-run difference.  If it does not,
      stop and report.
    - Run it 5× within one process (fresh `EdgeTraceSet` each time), and once
      in each of 5 fresh processes.  Report the first step whose hash
      differs, and whether the within-process runs differ at all.  A
      difference across processes but not within one points to allocation or
      alignment.

3. **Pinpoint the operation and the mechanism; then stop for review.**
    - Inside the first diverging step, record the inputs to each numeric call
      and recursively narrow down to the **single function call** whose output
      differs for bit-identical inputs (candidates: `moment.moment1d`, the
      fitting inside `trace.fit_trace` / `fitting.robust_fit` / `legfit` /
      `np.linalg.lstsq`, the PCA, a `np.sum`/`np.dot`).  Save its inputs with
      `np.save`.  Write `pypeitdev/issue_2208/replay.py`, which calls it alone.
    - Alignment test: replay with copies of the inputs placed at byte offsets
      0, 8, 16, ..., 56 (and 64) within an over-allocated buffer.  Report the
      distinct outputs per offset.
    - BLAS test: replay under `OPENBLAS_CORETYPE=Haswell`, `=Prescott` and the
      default.  If needed, create a throwaway env `pypeit14b_mkl` (a clone with
      MKL-backed numpy; Q4 approved) and replay with and without
      `MKL_CBWR=COMPATIBLE`.  Never modify `pypeit14b`.
    - Check whether the same call also explains the Binospec tilt-tracing
      entry point (`wavetilts`/`core/tracewave.py`).  A code read is enough;
      run something only if it is cheap.
    - **Stop.**  Write up the mechanism, with evidence, and 1–3 candidate fixes
      with their expected cost, numerical impact and scope (only the call
      site, or a shared helper).  Wait for my choice in Q&A/Testing.  A
      rounding/quantization fix needs my explicit approval.

4. **Implement the fix and the CI unit test (Q5a).**
    - Implement the chosen fix on `edge_2208`, matching the surrounding code
      style.  Keep the change minimal and do not change the algorithm's
      results beyond rounding level.
    - Add a test to the relevant `pypeit/tests/test_*.py`.  On a small,
      deterministic synthetic input (fixed seed) that exercises the same
      path, call the function on inputs at several byte offsets and assert
      bit-identical outputs.  Show that the test **fails** on the unfixed code
      (temporarily revert the fix by hand, without git) and passes with the fix.
    - Rerun the harness from Testing #2 (5 processes) and `replay.py`: they
      must now be identical.
    - Run the PypeIt unit tests that touch the changed modules, plus the
      edge/trace/tilt tests (`pytest pypeit/tests -k "trace or edge or tilt or
      fit or moment or pca"`), and then the full `pypeit/tests` suite.
    - Add a changelog entry with the `update-changelog` skill.  Log the
      commit point.

5. **Validate end to end (detached).**
    - Write `pypeitdev/issue_2208/run_pairs.sh`.  It runs, sequentially and
      cold, **two full serial reductions** each of
      `keck_lris_blue/multi_600_4000_slitmask`, `keck_deimos/830G_M_8500`
      and `mmt_binospec/Longslit_G1000`, then runs `compare_calibs.py` on each
      pair and writes a summary file.  Launch it detached and report the
      session and log paths.
    - When it finishes (a later invocation of this prompt may just collect
      the results), report for each setup the first differing array, if any.
      Goal: LRIS and DEIMOS are fully bit-identical.  If Binospec still
      differs in `Tilts`, report whether it is the same mechanism.  Per Q1,
      fix it here only if it is; otherwise draft a new issue in the log and do
      not post it.

6. **DevSuite `reduce` + `vet` on an edge-tracing subset (detached).**
    - From `./pypeit_test list` and `test_setups.py`, choose about 15–20
      setups that exercise edge tracing hard: multislit (LRIS blue/red,
      DEIMOS, MOSFIRE, GMOS, Binospec, ...), echelle (e.g. X-Shooter, NIRES,
      HIRES, MagE) and one or two long-slit controls.  Check the memory
      table in `README.rst` before choosing `-t 4`.  Record the list and why in
      the log.
    - Write `pypeitdev/issue_2208/run_devsuite_2208.sh`.  It runs
      `pypeit_test reduce afterburn -t 4 -s <setups>` into a dedicated
      `--outdir`, then `pytest vet_tests --redux_out <that dir>`, each with a
      report file (`-r`).  Launch it detached (named `tmux` session plus
      `nohup`), so the user can log out.  Report the commands, session and log
      paths, and the expected wall time.

7. **Collect, summarize, and draft the GitHub text.**
    - Collect the results of Testing #5 and #6.  Triage any vet failure
      (`diagnose-reduction` skill): is it caused by the fix, or pre-existing
      on `develop`?
    - Write `pypeitdev/issue_2208/Reports/issue_2208_results.md`: mechanism,
      evidence, fix, test, and validation numbers.
    - Draft (in the log, **not posted**) a comment for #2208 and a PR
      description for `edge_2208` → `develop`.  I will review them, and do the
      git work, before anything is posted with `gh`.

## Q&A

### Plan

#### What I read, and what it changes

- **#2208** (yours): from a bit-identical `SOBELSIG`, `EDGE_CEN` differs by
  exactly 1 float32 ULP (6.1e-5 px, so values of 512-1024 px) between two
  cold serial runs of `keck_lris_blue/multi_600_4000_slitmask`.  Pinning
  BLAS/OpenMP threads and `PYTHONHASHSEED=0` do not help.  The B-spline,
  unseeded RNGs and the sklearn PCA solver were ruled out.  The issue points
  at `trace.masked_centroid`; the stated hypothesis is SIMD results that
  depend on alignment.
- **#2209**: the same signature shows up on macOS for LRIS.  Binospec has a
  *second* entry point, in tilt tracing (`Tilts`, ~1e-6 px), with
  `EDGE_CEN` identical there.  Kast is bit-identical.  The PR pins BLAS to 1
  thread on both its serial and parallel paths.

One point from the code: **`masked_centroid` is not necessarily where the
divergence starts.**  `EDGE_CEN` is overwritten several times in
`EdgeTraceSet.auto_trace` (`edgetrace.py:663`):

| step | writes `edge_cen` via | code that could vary |
|---|---|---|
| `initial_trace` | `edge_img` (ints) | none |
| `centroid_refine` (`:1882`) | `follow_centroid` -> `masked_centroid` -> `moment.moment1d` | weighted sums |
| `fit_refine` (`:3084`) | `trace.fit_trace`: an iterative centroid -> Legendre fit -> recentroid loop | **least squares (LAPACK/BLAS)** |
| `pca_refine` / `peak_refine` (`:3630`) | `peak_trace` (centroids seeded by the PCA model) | PCA, fits |
| `sync` / insertions (`:4353`, `:4363`) | inserted traces from the PCA/offsets | PCA |

`SOBELSIG` is fixed, so what differs is the *starting position* passed into
each recentroid (`edge_fit` from a fit), or a reduction inside one of these
steps.  The saved `EDGE_CEN` is only the last of these writes.

My leading hypothesis is a more specific form of yours: the numpy wheel in
`pypeit14b` (numpy 2.5.0, Python 3.14.6) uses `scipy_openblas64` 0.3.33 with
the **SkylakeX (AVX-512) kernels**.  OpenBLAS's AVX-512 `dgemv`/`ddot`/`dgemm`
kernels split off a scalar prologue and epilogue depending on the pointer
alignment, so a **single-threaded** call can round differently when the
allocator hands back a differently aligned buffer.  That would ignore thread
pinning and the hash seed, and it would appear in the `lstsq`/`dot`
calls in `fit_trace` and the PCA.  numpy's own pairwise-sum SIMD loops are a
second candidate (`moment1d`).

#### Proposed plan (workstation `profx`, env `pypeit14b`)

The PypeIt branch `edge_2208` is currently at `origin/develop` (`0fd778f07`),
and the dev-suite branch is `dev_edge_2208`.

**Phase 0: baseline and environment (≈30 min)**
1. Record the numpy/scipy/BLAS configuration (`np.show_config()`,
   `threadpoolctl`) in the log.
2. Reproduce #2208 exactly: two cold `run_pypeit -c` (calibrations only)
   runs of `keck_lris_blue/multi_600_4000_slitmask` into scratch dirs, then
   compare `Edges_*_DET0{1,2}`.  Also record *how many* `EDGE_CEN` elements
   differ and which traces they belong to; the issue gives only the max.

**Phase 1: a fast standalone harness (no `run_pypeit`)**
3. A script (kept in `pypeitdev/issue_2208/`) that loads `TRACE_IMAGE`/BPM
   from the saved Edges file (or rebuilds the `TraceImage` from the
   calibration frames), builds an `EdgeTraceSet` with the same par and
   `auto_trace`s it.  After each step it records a hash of `edge_cen`,
   `edge_fit`, `edge_msk`, and the PCA arrays.
4. Run it N times **within one process** and in **N fresh processes**.  This
   separates within-process from across-process variation, a distinction
   that points to the allocator and alignment.

**Phase 2: bisect to the first diverging step, then to the line**
5. Find the first `auto_trace` step whose output differs.
6. Inside that step, dump its inputs (pickle) and replay them alone in fresh
   processes.  Narrow down to the function (`moment1d`, `fit_trace`'s fitting
   call, `np.linalg.lstsq`/`polyfit`/`legfit`, PCA, ...).
7. Test the cause directly:
   - **Alignment:** replay with deliberately misaligned copies of the
     inputs (views into an over-allocated buffer offset by 8, 16, ... 56
     bytes).  If the output changes with the offset, that confirms it.
   - **BLAS kernel:** rerun with `OPENBLAS_CORETYPE=Haswell` (and `=Prescott`)
     to swap out the AVX-512 kernels.  If the result becomes reproducible,
     it is OpenBLAS.
   - If useful, an MKL-backed numpy in a *cloned* env with `MKL_CBWR=COMPATIBLE`
     as the control (this needs your OK, see Q4).

**Phase 3: fix**
8. Prefer a targeted algorithmic fix in PypeIt over an environment knob, for
   example (depending on what Phase 2 finds):
   - replace an alignment-sensitive BLAS call on small arrays with an
     order-fixed computation (e.g. the normal equations or a `solve` on
     contiguous, freshly copied arrays, or `np.einsum` with
     `optimize=False`);
   - accumulate the moment sums in a fixed order (e.g. `math.fsum` or an
     explicit float64 loop over the window; the windows are only ~6 px).
   If the fix is only environmental (e.g. `OPENBLAS_CORETYPE`), I will report
   back rather than adopt it.
9. Add a regression test.  If the culprit is a pure function, add a
   PypeIt unit test (`pypeit/tests/`) that calls it on misaligned copies of
   the same input and asserts bit-identical output.  It needs no RAW_DATA,
   so it runs in CI.

**Phase 4: validate with the DevSuite**
10. Two cold serial LRIS slitmask runs.  The goal: `Edges`/`Slits`/`Flat`/
    `Tilts` and then spec1d/spec2d bit-identical (if the Binospec-style tilts
    entry point also affects LRIS, it shows up here as the *next* diff).
11. The same check on `keck_deimos/600ZD_M_6500` (or the cheaper
    `830G_M_8500`) and `mmt_binospec/Longslit_G1000`.
12. `./pypeit_test reduce` (with `-t`) on a set of multislit and echelle
    setups whose edge tracing changes, to show the fix moves results only at
    the rounding level.  Then the matching `vet` tests.

Everything gets logged below, with timings and the exact commands.

#### Questions

**Q1 — Scope.** Is this issue only the `EDGE_CEN` entry point, or should I
also chase the Binospec tilt-tracing entry point (and any others found on
the way), with full spec1d/spec2d run-to-run reproducibility as the goal?
I propose: fix `EDGE_CEN` here; if the cause turns out to be the same
mechanism (e.g. the same BLAS call), fix the tilts too, otherwise file it
separately.
>A: Use your proposal

**Q2 — Base branch.** `edge_2208` is off `develop`, without #2209's thread
pinning.  The issue says pinning makes no difference to `EDGE_CEN`, so I
plan to stay on `develop` and export the `*_NUM_THREADS=1` variables
myself in the harness, which keeps the fix independent of #2209.  Is that OK,
or do you want it stacked on `speed_up_B`?
>A: Yes, stay on `edge_2208` which is off `develop` without #2209's thread pinning.

**Q3 — Acceptable result changes.** A determinism fix will change edge
positions at the ~1e-12 px level relative to `develop`, and that can still
cascade into visible spec1d differences in the dev suite (the same
amplification, just once).  Is that acceptable as long as the vet tests
pass?  Also, is a quantization-style fix (e.g. rounding centroids to
1e-6 px at a checkpoint) off the table?  I would only use it as a last
resort, since it only makes divergence less likely.
>A: 1e-12 px is acceptable as long as the vet tests pass. 

**Q4 — Diagnostic environments.** Am I OK to create a throwaway clone of
`pypeit14b` (e.g. `pypeit14b_mkl`) to try an MKL- or reference-BLAS-backed
numpy as a control?  `pypeit14b` itself would not be touched.
>A: Yes, that is fine.

**Q5 — Where the regression test lives.** Options: (a) a PypeIt CI unit
test on the culprit function with misaligned inputs (cheap, my
preference); (b) a dev-suite `unit_tests/` test that runs LRIS edge tracing
twice in separate subprocesses (~minutes, needs RAW_DATA); (c) both.
>A: (a) 

**Q6 — DevSuite budget.** For Phase 4, is running `reduce`+`vet` on a
multislit/echelle subset at `-t 4` overnight fine, or do you want the full
`all` run?  And should I run anything else on the workstation at the same
time?
>A: Yes, that is fine.  But set it up so I can leave my Claude session and log out 
from my laption into my workstation and have it keep running.  The workstation will remain on.

### Testing

(Questions and answers for the Testing prompts go here.)

## Logs

### 2026-10-18 — Prompt Plan #1

- Read #2208 (body; no comments yet) and #2209 (body, the review by
  @tepickering, and both replies) with `gh ... --json`.  Plain `gh issue view`
  fails on the Projects-classic GraphQL deprecation.
- Read `trace.follow_centroid`/`masked_centroid` (`pypeit/core/trace.py:527,697`)
  and `EdgeTraceSet.centroid_refine` and `auto_trace` (`pypeit/edgetrace.py`).
  Found that `edge_cen` is rewritten by `fit_refine` (`fit_trace`),
  `peak_refine` and the trace insertions, not only by `centroid_refine`.
- Environment on `profx`: PypeIt branch `edge_2208` = `origin/develop`
  (`0fd778f07`); `pypeit14b` is Python 3.14.6, numpy 2.5.0, scipy 1.18.0,
  BLAS `scipy_openblas64` 0.3.33 using the SkylakeX kernels, 24 threads by
  default.  `RAW_DATA/keck_lris_blue/multi_600_4000_slitmask` is present.
  `$PYPEIT_DEV` (`/home/xavier/Projects/...`) resolves to this repo.
- Wrote the plan and Q1–Q6 above.  Nothing run yet.

### 2026-10-18 — Prompt Plan #2 (answers reviewed; Testing prompts written)

- Read the answers to Q1–Q6.  Decisions: fix the `EDGE_CEN` entry point, and
  the Binospec tilts entry point only if it is the same mechanism (Q1).  Stay
  on `edge_2208` off `develop`, with no #2209 pinning; the harness exports
  the thread variables itself (Q2).  ~1e-12 px changes are fine if the vet
  tests pass (Q3).  A throwaway MKL clone env is OK (Q4).  The regression test
  is a PypeIt CI unit test (Q5a).  `reduce`+`vet` on a subset at `-t 4`,
  detached so the user can log out (Q6).
- One gap: Q3's question on rounding/quantization was not answered
  explicitly, so Testing #3 requires explicit approval before any
  quantization fix.  No other open questions.
- Wrote 7 prompts under Prompts/Testing: (1) baseline/reproduce and the
  `compare_calibs.py` tool; (2) stepwise `auto_trace` harness;
  (3) pinpoint the call and mechanism, then **stop for review**; (4) fix and
  CI unit test; (5) end-to-end pairs, detached; (6) DevSuite subset,
  detached; (7) results report and draft GitHub text (not posted).
- Added an empty Q&A/Testing section.  Nothing run yet.

