# Speeding up PypeIt — Design / Plan

**Version:** 0.5 (locked — ready to implement)
**Date:** 2026-06-30
**Author:** JXP and Claude

**Changelog**
- 0.1 (2026-06-29): Initial draft from the two profiles; ranked plan + questions.
- 0.2 (2026-06-30): Incorporated round-1 Q&A answers (decisions locked below);
  added the concrete 4-stage per-detector structural map of `reduce_exposure`,
  the slitmask cross-detector barrier, a strengthened identical-output design,
  and split the QA cheap wins into their own early PR.
- 0.3 (2026-06-30): Incorporated round-2 Q&A (CLI `--ncpu` flag; pin worker
  threads; v1 = clean stages 1/2/4 only, stage 3 deferred to v2; `fork`). Added
  the disk-backed lightweight-return design for workers and the failed-detector
  handling fix.
- 0.4 (2026-06-30): Incorporated round-3 Q&A (workstream order A→B→C; parallel
  logging design chosen; lightweight returns confirmed). Added §8 Implementation
  roadmap. Design considered converged.
- 0.5 (2026-06-30): Incorporated round-4 Q&A (QA PR scope finalized; PRs stack on
  the `speed_up` branch). All questions answered — design locked, ready to
  implement.

---

## 0. Decisions locked (round-1 → round-4 Q&A)

- **Mechanism:** stdlib `concurrent.futures.ProcessPoolExecutor`; new parameter
  **`par['rdx']['ncpu']`, default `1`** (opt-in; default run is byte-for-byte the
  current serial path). Smart cap `min(n_detectors, os.cpu_count()-1)`.
- **Also a CLI flag:** `run_pypeit --ncpu N`, which overrides the par value.
- **Granularity:** **detector/mosaic-level first**; coarser levels deferred.
- **v1 scope:** parallelize only the **clean flat stages 1 (calib), 2 (process),
  4 (extract)** — ≈70% of DEIMOS. **Stage 3 (objfind+sky) deferred to v2**
  because of its cross-detector slitmask barrier.
- **Worker threading:** pin each worker to single-threaded BLAS/OpenMP
  (`OMP_NUM_THREADS=1` etc.) for the identical-output guarantee + to avoid
  oversubscription. *Be prepared to revert if it backfires.*
- **Serial path untouched:** `ncpu=1` runs the literal current serial code; the
  pool is only entered when `ncpu>1`.
- **Start method:** Linux `fork` for v1 (`spawn`/cross-platform later).
- **Vectorized fits:** noise-level changes are acceptable, **vetted by dev-suite
  RMS checks**.
- **Parallel vs serial output:** must be **identical within machine precision**
  (a test will enforce this).
- **Validation:** Kast + DEIMOS + relevant `vet_tests` now; full dev suite later.
- **QA cheap wins:** ship as a **separate early PR**, ahead of the parallel work.
- **Workstream order:** **A (QA cheap-wins PR) → B (detector parallelism v1) →
  C (arc-fit vectorization).**
- **Parallel logging (my call, per Q12):** each worker captures its log records
  into an in-memory buffer (per-process handler) instead of writing the shared
  run log directly; on return it ships the buffer back with its result and the
  main process replays the lines **grouped in detector order**. This keeps each
  detector's "Calibrating detector X → … → done" block contiguous so the run log
  stays readable and the dashboard's per-mosaic timeline parsing keeps working.
  Tradeoff: a detector's detailed log appears when it finishes rather than
  streaming live; the main process still emits a coarse "det X started/finished"
  line for live progress. Revisit only if it proves awkward.
- **Worker returns (per Q13):** no large `DataContainer`s cross the process
  boundary — stage 1 returns `(det, success, failed_step)` (calibrations persist
  to disk, reload via `reuse_calibs`); stages 2/4 return as today, keyed by det.
- **QA cheap-wins scope (per Q14):** force matplotlib `Agg` unconditionally for
  headless reductions (no new param); parallelize the QA PNG writes under the same
  `ncpu` control; **no** skip-QA switch.
- **Branching (per Q15):** the workstreams **stack on the current `speed_up`
  branch** (A → B → C, each built on the prior), rather than branching
  independently off `develop`.

## 1. Purpose

A working plan to make `run_pypeit` faster under a deliberately conservative
("KISS") philosophy. It is grounded in the two cold-run profiles we have already
collected — `shane_kast_blue` (simplest: single detector, single slit) and
`keck_deimos 600ZD_M_6500` (heavy: 4 mosaics, slitmask, ~3.5 h) — whose reports
live alongside this file in `Reports/`.

This is a living document to iterate on; it intentionally proposes more than we
will necessarily do, ranked so we can pick.

## 2. Constraints (the brief)

- **KISS — Keep It Simple.** Avoid any significant refactoring; favor localized,
  low-risk changes over architectural rework.
- **Use multiple CPUs.** The pipeline is fully serial today.
- **Push more work into `numpy` array math** (replace Python-level loops and
  per-element `scipy` calls with vectorized array operations).
- **Ignore the B-spline kernel.** `pypeit/core/bspline/` is being optimized on
  another branch, so it is out of scope here — even though it is the single
  largest hotspot in both profiles. *Everything below excludes B-splines.*

## 3. What the profiles tell us (B-spline excluded)

Both profiles point at the same three classes of cost. Numbers below are from
the DEIMOS run (the realistic, multi-detector case; total wall **12 390 s**)
unless noted; the full tables are in `keck_deimos_600zd_m_6500.pstats.txt` and
`shane_kast_blue_A.pstats.txt`.

| Cost class | Evidence (non-B-spline) | KISS lever |
|---|---|---|
| **Serial per-detector work** | 4 mosaics reduced one after another; calibration build ≈ even ~1 300–1 600 s/mosaic; the whole science loop is serial | **Multi-CPU** |
| **Per-element Python/scipy loops** | arc-line fit: **1.11 M** `curve_fit` calls, `gauss_3deg` 332 s over **93 M** calls; `moment1d` 212 s over **1.05 M** calls | **Vectorize** |
| **Avoidable overhead** | QA: matplotlib text-metrics cum 894 s + PNG encode 276 s (971 PNGs); interactive `qtagg` backend used for headless PNGs; NumPy churn (`flatten` 351 s/80 M calls, `reduce` 312 s/233 M calls) | **Cheap wins** |
| **Mosaic resampling** | scipy.ndimage `geometric_transform` 462 s, `correlate` 577 s; `build_image_mosaic` cum 662 s — per detector | **Multi-CPU** (parallelizes for free) |

Key structural fact that makes multi-CPU *simple* here: the reduction is already
factored (on this `core_refactor`/`speed_up` branch) into clean module-level
functions that loop over detectors and return per-detector dicts —
`exposure.process_exposure`, `objfind_exposure`, `extract_exposure`
(`pypeit/exposure.py`), each calling `pypeit_steps.{process_one_det, objfind_one,
extract_det}`. The detector loop bodies are independent and side-effect-light.

## 4. The plan (ranked)

Ranked by (impact × how KISS). Effort/impact are rough.

### Tier 1 — Multi-CPU over detectors  ★ biggest win, lowest algorithmic risk

**The structure that makes this KISS.** `pypeit/exposure.py:reduce_exposure`
drives the whole per-exposure reduction as **four sequential per-detector
stages**, each looping `for det in detectors`:

| # | Stage | Code | Shape | DEIMOS cost |
|---|-------|------|-------|-------------|
| 1 | **Calibration build/load** | `reduce_exposure` lines ~421–433 → `pypeit_steps.calib_one` | flat per-det loop, writes independent FITS per det | flat+edges+tilts+wave ≈ **45%** |
| 2 | **Image process / mosaic** | `process_exposure` → `process_one_det` | flat per-det loop | mosaic resampling, img proc |
| 3 | **Object find + sky** | `findobj_on_exposure` → `findobj_on_det`, then **`adjust_for_slitmask`** (cross-detector), then `finalize_sky_det` | **fan-out → barrier → fan-out** | objfind + global sky |
| 4 | **Local sky + extract** | `extract_exposure` → `extract_det` | flat per-det loop | local-sky+extract ≈ **24%** |

Stages 1, 2, 4 are **clean flat loops** over independent detectors that assemble
results into dicts keyed by detector — swap the loop for a pool `map`, keyed the
same way, loop body untouched. **Together they are ~70% of the DEIMOS runtime
(calibration ~45% + extraction ~24% + image/mosaic).**

**Stage 3 is the one wrinkle:** `findobj_on_exposure` is *not* a flat loop. It
fans out `findobj_on_det` per detector, then hits a **cross-detector barrier** —
`adjust_for_slitmask` (`exposure.py:228`) aggregates objects/slits across *all*
detectors to compute the slitmask offset (DEIMOS-class instruments) — then fans
out `finalize_sky_det` per detector. Parallelizing it means fan-out → serial
barrier → fan-out, which is more than a one-line `map` swap. (See Q&A: candidate
to defer to a v2.)

**KISS v1 scope (locked).** Parallelize **stages 1, 2, and 4** (the clean flat
loops) — already ~70% of DEIMOS wall-clock — and leave stage 3 serial. Stage 3
(its slitmask barrier) is a **v2** item. This is the minimal-risk change that
captures most of the win.

**Workers return lightweight status, not big objects (disk-backed).** The code
already persists per-detector products to disk, which keeps the parallel design
simple and pickling cheap:
- *Stage 1:* `reduce_exposure`'s calib loop calls `calib_one`, which builds the
  calibrations (written to `Calibrations/*.fits`, one set per detector/mosaic) and
  returns a `caliBrate` object that the loop **discards** — it only checks
  `.success`. Downstream stages reload calibrations from disk via
  `reuse_calibs=True`. So a stage-1 worker need only return
  `(det, success, failed_step)`; the heavy object never crosses the process
  boundary.
- *Stage 2:* `process_one_det` already does `sciImg.to_file(...)`. Workers can
  return the (small) per-detector `sciImg` (or a path), consistent with today.
- *Stage 4:* returns the `spec2d`/`specobjs` to be assembled and saved by
  `save_exposure`, as today — keyed by detector.

**Failed-detector handling (fix a latent bug).** The current serial calib loop
does `detectors.remove(det)` *while iterating* `detectors` — which skips the next
element. The parallel rewrite collects `(det, success)` results and rebuilds the
surviving-detector list afterward (no mutation during iteration), which is both
correct and necessary for the pool.

**Why it fits the brief.** No algorithm changes; the loop bodies
(`calib_one`, `process_one_det`, `extract_det`) are untouched. Pure
parallelization of work already isolated into per-detector functions.

**Expected payoff.** DEIMOS has 4 independent mosaics → approaching ~4× on the
parallelized stages (~70% of the run). Single-detector instruments (Kast) gain
nothing from *this* lever — expected and fine; `ncpu=1` leaves them on today's
exact path.

**Risks / things to handle.**
- **Identical-output guarantee (locked requirement).** `ncpu=1` must be the
  literal current serial code path (no behavior change by default). For `ncpu>1`,
  results recombine by detector key (order-independent) and any cross-detector
  aggregation (the slitmask barrier) must iterate detectors in a fixed order. A
  regression test will assert `ncpu=1` vs `ncpu>1` outputs match to machine
  precision. **BLAS/OpenMP thread count can change floating-point reduction order
  inside a worker** → pin workers to single-threaded BLAS (also prevents N×threads
  oversubscription). (See Q&A.)
- **Memory.** DEIMOS already uses tens of GB; N detectors in flight ≈ N× peak
  RAM. The `min(n_detectors, os.cpu_count()-1)` cap (per Q1) bounds this; document
  the RAM-vs-`ncpu` tradeoff.
- **Picklability / start method.** Per-detector inputs and the returned
  `DataContainer`s must pickle; on Linux `fork` is the simplest start method and
  lets workers inherit the large read-only inputs cheaply. Per-process logging
  must funnel back to the main run log.

### Tier 1 — Vectorize arc-line Gaussian fitting  ★ biggest single vectorization win

**Idea.** `fit_arcspec` (`pypeit/core/arc.py:1138`) fits every detected arc line
with `scipy.optimize.curve_fit` (via `fitting.fit_gauss`) in a Python loop — the
source of the 1.11 M `curve_fit` calls and 93 M `gauss_3deg` evals. A 3-parameter
Gaussian has a closed-form fit: take `ln(flux)` over the line window and do a
weighted parabola (degree-2 polynomial) least-squares; amplitude/center/sigma fall
out of the 3 coefficients. This vectorizes across *all* lines at once with NumPy
(`np.polyfit`-style batched lstsq) — no per-line scipy call.

**Why it fits the brief.** Localized to one helper; replaces a slow per-element
loop with array math. Benefits wavelength calibration *and* tilt tracing (both go
through arc-line fitting).

**Expected payoff.** Targets ~840 s cum on DEIMOS; the log-parabola fit is
typically 10–100× faster per line and removes the 93 M-call `gauss_3deg` overhead.

**Risk.** Results change at the noise level (analytic vs iterative fit). Needs a
tolerance decision (Q&A) and a vet against existing wavelength/tilt RMS.

### Tier 2 — Cheap, near-zero-risk wins

- **Force a headless matplotlib backend (`Agg`) for QA.** The run uses the
  interactive `qtagg` backend to write non-interactive PNGs. Setting `Agg` for
  the reduction removes GUI/text-metric overhead (~part of the ~900 s QA cost).
  One-line, low-risk.
- **Parallelize / defer QA PNG writes.** 971 PNGs on DEIMOS. QA generation is
  independent per figure → a thread/process pool, or an option to defer/skip QA,
  recovers a chunk of the ~1.2 ks QA budget. (Pairs naturally with Tier 1.)
- **Reduce obvious NumPy churn** in hot helpers (`flatten`/`copy`/masked-array
  calls in fitting/rejection inner loops): operate in place, avoid masked-array
  where a boolean index suffices. Small, opportunistic, per-hotspot.

### Tier 2 — Vectorize / de-loop `moment1d` centroiding

`moment1d` (`pypeit/core/moment.py`) is called 1.05 M times via trace centroiding
(`trace.masked_centroid` / `follow_centroid`) and boxcar extraction. The function
is already NumPy-internally; the cost is call volume from row-by-row trace
following. Vectorizing across spectral rows is possible but more invasive than the
arc-fit change — flagged as a candidate, not a commitment, pending a closer look.

### Tier 3 — Coarser parallelism (optional, later)

Calibration groups, comb-id groups, and the standard/science frame loops in
`pypeit.py:reduce_all` / `reduce_calibID` are also independent and could be
parallelized. Lower priority: most single-config runs have one calib group, and
detector-level parallelism (Tier 1) already captures the within-exposure win with
finer memory granularity.

## 5. Validation strategy

- **Regression on both profiled setups.** Re-run Kast (single detector — guards
  the vectorization + cheap wins) and DEIMOS (multi-detector — guards the
  parallelism) cold, and compare against the existing reductions.
- **Numerical vetting via the dev suite.** Use the existing `vet_tests`
  (wavelengths, wavetilts, extraction, skysub, flux) to confirm products stay
  within tolerance after the arc-fit vectorization.
- **Re-profile** after each change with the existing
  `scripts/profile_*.py` + `analyze_profile.py` to measure the actual speedup and
  catch regressions.
- **Determinism check** for the parallel path (same outputs at `ncpu=1` and
  `ncpu>1`).

## 6. Explicitly out of scope

- The B-spline kernel (`pypeit/core/bspline/`) — other branch.
- Any large architectural refactor of the reduction driver, datamodel, or IO.
- GPU offload, Cython/C extensions (beyond what the b-spline branch may do).

## 7. Open questions

**None — all questions (rounds 1–4, Q1–Q15) are answered and folded into §0
"Decisions locked."** The design is locked; see §8 for the implementation
roadmap.

## 8. Implementation roadmap

Three workstreams in the locked order **A → B → C**, each its own PR, **stacked on
the `speed_up` branch** (each branches from the prior), each gated by the
validation in §5.

### A. QA cheap-wins PR (smallest, immediate)
- Force matplotlib `Agg` **unconditionally** for headless reductions (the run
  currently uses the interactive `qtagg` backend to write PNGs); no new param.
- Parallelize the QA PNG writes under the `ncpu` control (independent per figure).
  No skip-QA switch. Since A is first in the stack, the **`ncpu` parameter +
  `--ncpu` flag plumbing lands here** (a thread pool suffices for the
  render/write-bound PNGs); workstream B then reuses `ncpu` for its process pool.
- Re-profile Kast + DEIMOS; expect to recover a chunk of the ~0.9–1.2 ks QA cost.

### B. Detector parallelism v1 (the headline win)
1. Reuse the `ncpu` parameter + `--ncpu` flag introduced in A (default 1).
2. Add a small helper to map a per-detector function over `detectors` with a
   `ProcessPoolExecutor` when `ncpu>1`, else the literal serial loop (`ncpu=1`
   path untouched). Pin worker BLAS/OpenMP threads to 1; `fork` start method.
3. Convert the **three clean flat loops** to use it:
   - Stage 1 — calib loop in `reduce_exposure` (workers return
     `(det, success, failed_step)`; rebuild surviving-detector list, fixing the
     `detectors.remove()`-during-iteration bug).
   - Stage 2 — `process_exposure`.
   - Stage 4 — `extract_exposure`.
   Leave stage 3 (`findobj_on_exposure`) serial (v2).
4. Per-worker buffered logging, replayed in detector order (per §0).
5. Add the **identical-output regression test** (`ncpu=1` vs `ncpu>1`, machine
   precision) + the determinism check. Re-profile DEIMOS for the speedup.

### C. Arc-fit vectorization (self-contained)
- Replace the per-line `scipy.curve_fit` in `fit_arcspec` (`core/arc.py:1138`,
  via `fitting.fit_gauss`) with a vectorized analytic weighted log-parabola fit
  across all lines.
- Vet wavelength/tilt RMS via `vet_tests`; confirm within tolerance. Re-profile.

### v2 / later
- Stage 3 parallelism (fan-out → slitmask barrier → fan-out).
- `moment1d` de-looping; coarser (calib-group/frame) parallelism; `spawn` support.
