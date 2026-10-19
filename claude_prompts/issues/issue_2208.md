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
>A:

**Q2 — Base branch.** `edge_2208` is off `develop`, without #2209's thread
pinning.  The issue says pinning makes no difference to `EDGE_CEN`, so I
plan to stay on `develop` and export the `*_NUM_THREADS=1` variables
myself in the harness, which keeps the fix independent of #2209.  Is that OK,
or do you want it stacked on `speed_up_B`?
>A:

**Q3 — Acceptable result changes.** A determinism fix will change edge
positions at the ~1e-12 px level relative to `develop`, and that can still
cascade into visible spec1d differences in the dev suite (the same
amplification, just once).  Is that acceptable as long as the vet tests
pass?  Also, is a quantization-style fix (e.g. rounding centroids to
1e-6 px at a checkpoint) off the table?  I would only use it as a last
resort, since it only makes divergence less likely.
>A:

**Q4 — Diagnostic environments.** Am I OK to create a throwaway clone of
`pypeit14b` (e.g. `pypeit14b_mkl`) to try an MKL- or reference-BLAS-backed
numpy as a control?  `pypeit14b` itself would not be touched.
>A:

**Q5 — Where the regression test lives.** Options: (a) a PypeIt CI unit
test on the culprit function with misaligned inputs (cheap, my
preference); (b) a dev-suite `unit_tests/` test that runs LRIS edge tracing
twice in separate subprocesses (~minutes, needs RAW_DATA); (c) both.
>A:

**Q6 — DevSuite budget.** For Phase 4, is running `reduce`+`vet` on a
multislit/echelle subset at `-t 4` overnight fine, or do you want the full
`all` run?  And should I run anything else on the workstation at the same
time?
>A:

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