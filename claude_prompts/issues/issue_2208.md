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

**Q7 — Which fix?** (Testing #3; details in its log entry.)  The cause is
the unseeded randomized SVD in `pypeit/core/pca.py:120`.  Options:

- **(A) `svd_solver='full'` in both `PCA(...)` calls of
  `pca_decomposition` — my recommendation.**
  - An exact LAPACK SVD.  It is bit-reproducible in and across processes,
    at every memory offset, under all three OpenBLAS kernel types, and even
    with unpinned (24) threads.
  - Same cost (~4 ms per call at edge shapes; the tilt shapes are smaller).
  - It removes the dependence on sklearn's `'auto'` policy, which has
    changed between versions (`covariance_eigh` was added in 1.5).
  - Two-line change, in one function.  It covers edge tracing, tilt
    tracing (the Binospec entry point) and the echelle object PCA
    (`findobj_skymask`).
- **(B) `random_state=0`, keeping `'auto'`.**  Also deterministic in all of
  the same tests.  But it keeps an approximate solver, and the result still
  depends on sklearn's solver policy and on how the randomized algorithm is
  implemented, both of which can change across sklearn versions.  The same
  one-function scope.
- **(C) Replace sklearn with `np.linalg.svd` plus sklearn's sign
  convention (`svd_flip`).**  Deterministic, and it drops sklearn from
  PCA.  More code, and it has to match sklearn's sign/centering exactly.
  No benefit over (A) for this issue.
>A: (A)

**Q8 — Size of the result change (revisits Q3).**  Any of these fixes
moves the edges by up to **~4e-6 px** (median ~3e-8 px) relative to
current `develop`, not ~1e-12 px.  That is the **same size as the
difference between two unpatched runs today** (max 3.3e-6, median 3e-8
px), because the current output is one random draw from that scatter.
Downstream products (Flat, Tilts, spec1d/spec2d) will change by about as
much as they already do between two identical runs.  Is that acceptable,
with the vet tests as the gate?
>A: Yes, that is acceptable.

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


### 2026-10-18 — Testing #1 (baseline: #2208 reproduced; Haswell kernels do not help)

**Environment (profx).**
- CPU: Intel Xeon w5-2455X (Sapphire Rapids), 12 cores / 24 threads, full
  AVX-512 (incl. fp16, bf16, vnni).  755 GB RAM.  glibc 2.39 (Ubuntu).
- `pypeit14b`: Python 3.14.6, numpy 2.5.0 (SIMD baseline X86_V2; found X86_V3,
  X86_V4, AVX512_ICL, AVX512_SPR), scipy 1.18.0.
- **Two OpenBLAS libraries are loaded**, both DYNAMIC_ARCH and both
  selecting the **SkylakeX** kernels:
  - numpy's `libscipy_openblas64_` 0.3.33.112.0
  - scipy's `libscipy_openblas` 0.3.31.dev

  sklearn also brings in `libgomp` (OpenMP).
- `OPENBLAS_CORETYPE=Haswell` switches **both** libraries to Haswell
  (checked with `threadpoolctl`).
- `pypeit.__file__` is `/mnt/tank/Astronomy/PypeIt/PypeIt/pypeit` (the
  `edge_2208` tree).  The reported version, `2.0.2.dev1216+gf3a1f1d27`, is
  stale install metadata from `speed_up_B`, not the code that runs.

**Scripts written** (`pypeitdev/issue_2208/`):
- `compare_calibs.py <dir_a> <dir_b> [--all] [--files ...]`: compares every
  image HDU and table column in `Calibrations/` and `Science/`, in dependency
  order.  It reports the count of differing elements, max |diff|, and max ULP
  in the native dtype; floats are compared bitwise, NaN==NaN.  It takes ~7 s
  for the LRIS calibrations.
- `edge_cen_diff.py <Edges_a> <Edges_b>`: for each trace with differences,
  lists side, MASKDEF_ID, insertion flags, differing rows and max ULP, and
  how many EDGE_FIT rows differ.
- `run_calib.sh <name> [VAR=val ...]`: a cold `run_pypeit -c -o` of LRIS
  slitmask into `REDUX/issue_2208/<name>`, with threads pinned to 1.  It
  writes `blas_info.txt` (the threadpoolctl view under that environment),
  `run.log` and `timing.txt`.

**Smoke test on existing outputs** (`REDUX_OUT/.../multi_600_4000_slitmask`
serial vs its `ncpu/` run, both older than this branch).  It reproduces
#2208:
- Arc, Tiltimg and WaveCalib are identical.  In Edges DET01, 13 arrays
  differ.  `SOBELSIG`, `TRACE_IMAGE`, `TRACEID`, `MASKDEF_ID` and
  `EDGE_MSK` are identical.
- `EDGE_CEN` differs in **93 / 45 056** elements (DET01) and 1 041 / 57 344
  (DET02, max 9 ULP there).  In DET01 the differences are spread over 14 of
  the 22 traces, both left and right, all slitmask-matched, none inserted,
  at scattered rows.
- `EDGE_ERR` differs in many more elements (4 125), and so does `EDGE_FIT`
  (217).
- Reading: the in-memory float64 values probably differ almost everywhere
  by a tiny amount, and only occasionally cross a float32 rounding
  boundary on output.  This is not yet tested.

**Runs launched**, detached, in tmux session `i2208_t1`:
```
cd pypeitdev/issue_2208
tmux new-session -d -s i2208_t1 "bash -c './run_calib.sh base_a; ./run_calib.sh base_b; \
  ./run_calib.sh hsw_a OPENBLAS_CORETYPE=Haswell; ./run_calib.sh hsw_b OPENBLAS_CORETYPE=Haswell; \
  echo ALLDONE' > /mnt/tank/Astronomy/PypeIt/REDUX/issue_2208/testing1.log 2>&1"
```
- Progress log: `/mnt/tank/Astronomy/PypeIt/REDUX/issue_2208/testing1.log`
  (one `status=/wall=` line per run, then `ALLDONE`).
- `base_a`: status 0, **220.7 s**.  `blas_info.txt` confirms 1 thread and
  SkylakeX for both OpenBLAS libraries.  The other three runs were still
  going when the session ended (~4 min each, ~15 min total).

**Run results** (all status 0; every `blas_info.txt` shows 1 thread, and
`hsw_*` show **Haswell** for both OpenBLAS libraries):

| run | env | wall |
|---|---|---:|
| base_a | default (SkylakeX) | 220.7 s |
| base_b | default (SkylakeX) | 217.2 s |
| hsw_a | `OPENBLAS_CORETYPE=Haswell` | 217.2 s |
| hsw_b | `OPENBLAS_CORETYPE=Haswell` | 213.1 s |

Comparison outputs are saved in `/mnt/tank/Astronomy/PypeIt/REDUX/issue_2208/`:
`cmp_<a>_<b>.txt` (compare_calibs.py) and `edgecen_<a>_<b>_DET0N.txt`
(edge_cen_diff.py).

**File-level summary (number of arrays that differ):**

| file | base_a vs base_b | hsw_a vs hsw_b | hsw_a vs base_a |
|---|---|---|---|
| Arc DET01/02 | identical | identical | 1 / 1 (`ARC_IMAGE`, ~all px, ≤9e-8) |
| Tiltimg DET01/02 | identical | identical | 1 / 1 |
| WaveCalib DET01/02 | identical | identical | 82 / 118 |
| Edges DET01/02 | **13 / 16** | **13 / 16** | 15 / 18 |
| Slits DET01/02 | 7 / 7 | 7 / 7 | 7 / 7 |
| Flat DET01/02 | 30 / 42 | 30 / 42 | 31 / 43 |
| Tilts DET01/02 | 3 / 3 | 3 / 3 | 5 / 5 |

**Edges, base_a vs base_b:**

| | DET01 | DET02 |
|---|---|---|
| SOBELSIG, TRACE_IMAGE, TRACEID, MASKDEF_ID, EDGE_MSK | identical | identical |
| `EDGE_CEN` differing / total | **66 / 45 056** (17 of 22 traces), max 1 ULP | **1 393 / 57 344** (20 of 28 traces), max **14 ULP** |
| `EDGE_FIT` | 236, max 1 ULP | 2 406, max 128 ULP |
| `EDGE_ERR` | 4 103, max 18 ULP | 5 787, max 24 ULP |
| `PCA:*`, `PCA_MODEL_*` | all differ (trace_coo up to 1.3e-6 px) | all differ (trace_coo up to 2.9e-6 px) |

**Edges, hsw_a vs hsw_b:** the same picture.  DET01 `EDGE_CEN` 102 / 45 056
(14 traces, max 1 ULP); DET02 1 788 / 57 344 (22 traces, max 20 ULP).
SOBELSIG and the other inputs are identical.

**Per-trace `EDGE_CEN` pattern (both pairs):**
- No differing trace is inserted (no `*INSERT` flags).  All of them are
  slitmask-matched (`MASKDEF_ID` set), and left and right edges are affected
  alike.
- Almost every trace shows a handful of **scattered single rows**, each
  exactly 1 ULP.  If the float64 values in memory differ by a small,
  roughly uniform δ, the chance that a value crosses a float32 rounding
  boundary is about δ/ULP.  With DET01's 66/45 056 ≈ 1.5e-3 and 1 ULP =
  6.1e-5 px, that gives **δ of order 1e-7 px** (an estimate, not a
  measurement).
- **Exception: DET02 column 0** (TRACEID −1, a left edge, MASKDEF_ID 7)
  differs on **1 342 / 2 048 rows** (base pair) and 1 680 rows (hsw pair).
  Its max is 14–20 ULP (≈1e-3 px), and its `EDGE_FIT` differs on 2 034–2 048
  rows.  This trace has a real, trace-wide difference, so it is a good
  first target for Testing #2/#3.
- The same rows recur across pairs (e.g. DET01 trace −2 rows 7, 41, 43, 87,
  147, ...).  The values there sit close to a float32 rounding boundary.

**Haswell discriminator (recorded, no conclusion beyond this):**
- `hsw_a` vs `hsw_b`: still **not reproducible**, with the same structure and
  magnitude as the SkylakeX pair.  Arc/Tiltimg/WaveCalib are identical; the
  first difference is again in Edges, from an identical SOBELSIG.  So
  swapping OpenBLAS's AVX-512 kernels for AVX2 ones does **not** remove the
  run-to-run difference.
- `hsw_a` vs `base_a`: the kernel choice changes the processed images
  themselves (`ARC_IMAGE`, `TRACE_IMAGE` at ≤9e-8, so SOBELSIG differs too),
  and everything downstream.  This comparison therefore cannot be used to
  study Edges.  It does confirm that image processing is
  **deterministic for a fixed kernel** but kernel-dependent.

**Tool note.** In `compare_calibs.py`, a ULP of `9223372036854775807`
(INT64_MAX) means the two values have opposite signs (SOBELSIG near 0); it is
not a real distance.

**Status:** Testing #1 is complete.  Next is Testing #2 (stepwise harness),
with DET02 trace 0 as the primary target.

Commit point (for the user): `pypeitdev/issue_2208/{compare_calibs.py,edge_cen_diff.py,run_calib.sh}`.

### 2026-10-19 — Testing #2 (harness: first divergence is `build_pca` in `pca_refine`, not centroiding)

**Harness:** `pypeitdev/issue_2208/edge_harness.py`
- `run`: rebuilds the spectrograph and `EdgeTracePar` from the `.pypeit`
  file with `PypeItFile.get_pypeitpar()`, as `run_pypeit` does.  It loads the
  `TraceImage` from a saved Edges file (fresh for every repetition) and sets
  `traceimg.files` to the raw trace frame (needed by `maskdesign_matching`).
  It then runs `EdgeTraceSet(..., auto=True)`.
  - Every method defined on `EdgeTraceSet` is wrapped.  After each call
    returns, a SHA-256 is taken of the **in-memory** `edge_img/cen/err/msk/fit`,
    `traceid`, `maskdef_id`, `sobelsig` and every PCA array (including the
    coefficient-model `xval/yval/fitc/gpm`).  Data-pointer alignment mod 64 is
    also recorded.
  - Writes `<tag>_repN.json` (the call sequence and hashes), `_steps.npz`
    (the full state after each top-level `auto_trace` step) and `_final.npz`.
  - It refuses to run unless the three `*_NUM_THREADS=1` variables are
    exported.
- `compare`: reports the first differing call and the status of each
  top-level step.  `check`: compares a final state with Edges files.
- `copy.deepcopy(TraceImage)` hits the `BitMaskArray.__getattr__` recursion
  on `develop` (fixed in #2209), so the image is re-read from FITS instead.
- About 10 s per repetition for DET02 and 6 s for DET01; 96 / 110 wrapped
  calls.

**Faithfulness check** (float32 cast of the harness `edge_cen` vs the
`run_pypeit` files; SOBELSIG identical in every case):

| det | harness vs base_a | harness vs base_b | base_a vs base_b (Testing #1) |
|---|---:|---:|---:|
| DET01 | 88 | 78 | 66 |
| DET02 | 1 965 | 1 809 | 1 393 |

The per-trace pattern on DET02 is the same as in Testing #1 (trace 0 on most
rows, up to 19–24 ULP; scattered 1-ULP rows elsewhere).  The harness
therefore reproduces the same kind of run-to-run difference.

**Runs** (DET02, `REDUX/issue_2208/harness/det02/`): 5 repetitions in one
process (`inproc`, 86 s total), then 5 fresh processes (`proc0..4`).  All
10 final `edge_cen` hashes are **different**, including the 5 repetitions
within one process.  DET01 (`harness/det01/`, 2 in-process plus 1 fresh) is
the same.

**First diverging call** (identical for all 9 comparisons on DET02 and both
on DET01): **call #25, `build_pca` (depth 2, inside `pca_refine`)**.
- Differing: `pca_coeffs`, `pca_components`, and the coefficient models'
  `yval`/`fitc`.
- **Identical: its inputs `trace_coo` and `pca_mean`, and all of `edge_cen`
  and `edge_fit`.**

Top-level steps, inproc rep0 vs rep1 (DET02), with float64 max |diff| from
`_steps.npz`:

| step | status | size |
|---|---|---|
| `initial_trace` | identical | |
| `centroid_refine` | identical | |
| `fit_refine` | identical | |
| `merge_traces` | identical | |
| `pca_refine` | **differs** | `edge_fit` 1.1e-13, `pca_coeffs` 1.9e-13, `pca_components` 7.8e-16 |
| `peak_refine` | differs (amplified) | `edge_cen` 8.8e-7 px, `edge_fit` 4.3e-6, `pca_coeffs` 1.1e-5 |
| `maskdesign_matching`, `sync` | differ (inherited) | same as `peak_refine` |

- **`masked_centroid`/`centroid_refine` and `fit_trace` are deterministic.**
  The issue's "Where to look" is off by two steps.  The saved `EDGE_CEN`
  differs only because `peak_refine` re-measures the edges from PCA-
  predicted positions.
- `peak_refine` turns a 1e-13 difference into 1e-6 px.  After it, **87 % of
  the float64 `edge_cen` elements differ** (median |d| ≈ 1.8e-8 px, max
  8.3e-7 px).  The float32 file shows only the elements that cross a
  rounding boundary, consistent with the Testing #1 estimate.
- **The DET02 trace 0 "outlier" is not a separate mechanism.**  That left
  edge lies at x ≈ 0.48–11.1 px, where the float32 ULP is 6e-8 to 1e-6 px,
  so a ~1e-7 px difference shows up as many ULPs on every row.

**Lead for Testing #3 (not yet confirmed on the real arrays):**
- Variation *within* one process rules out the "fresh process / allocator
  alignment" explanation; it points to a random number generator.
- `pypeit/core/pca.py:120` calls `PCA(n_components=npca)` with no `svd_solver` and no
  `random_state`.  On a synthetic array of PypeIt's shape (28 × 2048,
  sklearn 1.9.0), `svd_solver='auto'` resolves to **`full`** for
  `n_components=None` (line 92) but to **`randomized`** for
  `n_components=3` (line 120).  Two unseeded fits differ (max 3.1e-5), and
  two fits with `random_state=0` are bit-identical.
- If this holds on the real arrays, #2208's ruled-out item 5 ("sklearn's
  randomized SVD … resolves to `full` or `covariance_eigh`") is wrong for
  the line-120 call.  It would also explain why pinning threads and
  `PYTHONHASHSEED` changed nothing, and why `OPENBLAS_CORETYPE` did not help.
- Testing #3 will confirm this with a replay of the real `build_pca` inputs,
  then stop for review as planned.

Commit point (for the user): `pypeitdev/issue_2208/edge_harness.py`.

### 2026-10-19 — Testing #3 (mechanism: unseeded randomized SVD in `pca_decomposition`; STOP for review)

**Summary.**  The run-to-run irreproducibility in multi-slit edge tracing
comes from an unseeded random number generator.  It is not floating-point
or SIMD/BLAS nondeterminism.

`pypeit.core.pca.pca_decomposition` calls `sklearn.decomposition.PCA` twice:
- `PCA()` (line 92) to choose the number of components.  sklearn's
  `'auto'` resolves this to the deterministic `full` solver.
- `PCA(n_components=npca)` (line 120) for the actual decomposition.  At
  PypeIt's edge-tracing shapes, `'auto'` resolves this to **`randomized`**,
  and PypeIt passes no `random_state`, so sklearn seeds it from global
  entropy on every call.  This is also against the dev suite's fixed-seed
  rule.

sklearn 1.9.0's rule (`_pca.py:524-536`) picks `randomized` when the array
is not tall-and-skinny, `max(shape) > 500`, and `1 ≤ n_components <
0.8·min(shape)`.

**New tools** (`pypeitdev/issue_2208/`):
- `pca_patch.py` patches PypeIt at runtime, without editing it.  It replaces
  `pypeit.core.pca.PCA` with a factory (`fix=None|'full'|'seed'`) and wraps
  `pca_decomposition` to record its inputs, the solver sklearn actually used,
  and the caller.  It can also save each call to `pca_call_<pid>_<n>.npz`.
- `site/sitecustomize.py` applies the patch to any process started with
  `PYTHONPATH=site` and `I2208_PCA_FIX`/`I2208_PCA_CAPTURE` set.  That
  covers `run_pypeit`.
- `edge_harness.py run` gained `--pca_fix {full,seed}` and `--capture <dir>`.
- `replay.py` replays captured calls `--nrep` times at each of `--offsets`
  byte offsets from 64-byte alignment.  It reports the distinct outputs,
  max |d|, an output hash comparable across processes, and ms per call.

**1. The real calls** (capture from the unpatched harness, `harness/capture_det0N/`):

| det | call | shape (vectors) | npca → used | solvers (`PCA()`, `PCA(n)`) | caller |
|---|---|---|---|---|---|
| 01 | 0 | 21 × 2048 | None → 2 | full, **randomized** | `pca_refine` → `build_pca` |
| 01 | 1, 2 | 21 × 2048 | None → 2 | full, **randomized** | `_reset_pca` → `build_pca` |
| 02 | 0 | 28 × 2048 | None → 2 | full, **randomized** | `pca_refine` → `build_pca` |
| 02 | 1 | 28 × 2048 | None → 3 | full, **randomized** | `_reset_pca` → `build_pca` |

`pca_decomposition` is the **only** sklearn use in PypeIt (`grep sklearn`).
It is reached from `EdgeTraceSet.build_pca`, from
`tracepca.pca_trace_object` in tilt tracing (`core/tracewave.py:642`, which
always runs it, with `npca=2`), and from the echelle object PCA
(`core/findobj_skymask.py:980`).

**2. Replay of the real inputs** (SkylakeX, 1 thread; 3 repetitions × 8
offsets {0,8,...,56} = 24 calls each, run in 2 separate processes):

| fix | distinct outputs | per offset | across the 2 processes | max |d| | ms/call |
|---|---|---|---|---|---|
| none | **24 / 24** (every call) | 3 / 3 | different | 1.9e-13 – 2.0e-12 | 3.3–4.3 |
| `svd_solver='full'` | **1 / 24** | 1 | **same hash** | 0 | 3.1–4.3 |
| `random_state=0` | **1 / 24** | 1 | **same hash** | 0 | 3.5–4.3 |

- **Memory alignment:** it plays no part.  Unpatched, every call differs
  at every offset; fixed, every offset gives the same output.
- **BLAS kernel** (`OPENBLAS_CORETYPE=Haswell` and `=Prescott`; the latter
  reports `Katmai`): the same pattern under every kernel.  Unpatched
  varies, fixed is deterministic.  Only the fixed *values* change with the
  kernel, which explains the `hsw_a` vs `base_a` image differences in
  Testing #1 but is not a reproducibility problem on one machine.
- **Threads:** with `full` and **unpinned** BLAS (24 threads) the hash
  equals the pinned one, so at these shapes it is also independent of the
  thread count.
- An MKL clone env (Q4) was **not needed**: the mechanism is an RNG, not
  BLAS, and every OpenBLAS kernel type behaves the same way.

**3. Edge tracing with each fix** (harness; 3 in-process repetitions + 2
fresh processes per detector; `harness/fix_{full,seed}_det0N/`):
- With either fix, on both detectors, **"all calls identical"**: every one
  of the 110 (DET01) / 96 (DET02) wrapped `EdgeTraceSet` calls has an
  identical hash in all 5 runs.  So there is **no second nondeterministic
  source** in edge tracing.

**4. Full LRIS calibrations with the `full` patch** (cold `run_pypeit -c`
pair `full_a`/`full_b` through `site/sitecustomize.py`; 225.1 s and 226.9 s;
tmux `i2208_t3`; `cmp_full_a_full_b.txt`):
- **Every file is bit-identical** on both detectors: Arc, Tiltimg,
  WaveCalib, **Edges, Slits, Flat, Tilts**.  The unpatched pair differed in
  every file from Edges on.
- Patch-active check: the patched harness `edge_cen` equals `full_a`'s
  `EDGE_CEN` exactly (0 differing elements on both detectors), against
  127 / 1 050 differences from unpatched `base_a`.
- Cost: none measurable (225–227 s against 213–221 s unpatched, within
  run-to-run timing noise).

**5. Size of the change** (float64 harness finals, max / median |d| in px):

| | DET01 `edge_cen` | DET02 `edge_cen` |
|---|---|---|
| unpatched rep0 vs rep1 (today's scatter) | 3.3e-6 / 3.0e-8 | 8.8e-7 / 1.5e-8 |
| `full` vs unpatched | 4.3e-6 / 2.7e-8 | 1.2e-6 / 2.8e-8 |
| `seed` vs unpatched | 2.5e-6 / 3.0e-8 | 3.4e-6 / 2.2e-8 |
| `full` vs `seed` | 4.3e-6 / 2.1e-8 | 2.4e-6 / 2.0e-8 |

`edge_fit` shows the same picture (≤5e-6 px).  Any fix changes the
results by the size of the existing scatter, **not** ~1e-12 px; see Q8.

**6. Binospec tilt-tracing entry point** (#2209; code read plus shapes, not
run):
- `trace_tilts` always calls `pca_trace_object(..., npca=2)`, which leads
  to `pca_decomposition(vectors=(n_lines, nsub))`.
- Binospec Longslit_G1000 has one slit per detector spanning ~2 460
  spatial px (DET01 832–3292, DET02 820–3282, from the existing
  `REDUX_OUT` Tilts files).  So `max(shape) > 500`, and with ≥3 lines
  `2 < 0.8·n_lines`: the solver is **randomized**, the same mechanism.
- LRIS slits are tens of pixels wide, so its tilts PCA stays at
  `max(shape) ≤ 500` and uses `full`.  That is why LRIS Tilts only
  inherited the Edges differences, and why Binospec's first difference is
  in Tilts with `EDGE_CEN` identical (Binospec long-slit edges cannot use
  the PCA).
- Per Q1, the same fix covers it.  Testing #5 checks it end to end.

**Corrections to #2208's text** (for the eventual comment):
- "What it is not" item 5 is wrong for line 120.  `svd_solver='auto'`
  resolves there to `randomized`, not `full`/`covariance_eigh`, for every
  PCA in this data set.
- "Where to look" (`masked_centroid`) is two steps too early.
  `centroid_refine` and `fit_refine` are deterministic; `EDGE_CEN`
  differs because `peak_refine` re-measures from PCA-predicted positions.
- The alignment hypothesis is ruled out: outputs differ even between
  repetitions inside one process, and are independent of byte offset.

**Candidate fixes:** see **Q7** and **Q8** in Q&A/Testing.  My
recommendation is (A), `svd_solver='full'`.  The test for Testing #4:
in `pypeit/tests/test_pca.py`, call `pca_decomposition` repeatedly, at
several byte offsets, on a fixed-seed synthetic (28, 2048) array that
`'auto'` sends to `randomized`, and assert bit-identical outputs (it fails on
the current code).

**STOPPED here as Testing #3 requires**; waiting for the answers to Q7/Q8.

Commit point (for the user): `pypeitdev/issue_2208/{edge_harness.py,pca_patch.py,replay.py,site/sitecustomize.py}`.

