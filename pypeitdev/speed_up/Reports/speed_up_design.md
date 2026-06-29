# Speeding up PypeIt — Design / Plan

**Version:** 0.1 (draft — for iteration)
**Date:** 2026-06-29
**Author:** JXP and Claude

---

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

**Idea.** Parallelize the existing per-detector loops with a process pool
(`concurrent.futures.ProcessPoolExecutor` / `multiprocessing.Pool`). The three
loops in `pypeit/exposure.py` (`process_exposure`, `objfind_exposure`,
`extract_exposure`) already iterate `for det in detectors` and assemble results
into dicts keyed by detector — swap the loop for a pool `map`, keyed the same way.

**Why it fits the brief.** No algorithm changes; the loop bodies
(`process_one_det`, `objfind_one`, `extract_det`) are untouched. Pure
parallelization of work that is already isolated.

**Expected payoff.** DEIMOS has 4 independent mosaics → up to ~4× on the
per-detector portion, which is the bulk of the run (calibration build + objfind +
extraction + mosaic resampling). Single-detector instruments (Kast) gain nothing
from *this* lever — that's expected and fine.

**Risks / things to handle.**
- **Memory.** DEIMOS already uses tens of GB; N detectors in flight ≈ N× peak
  RAM. Needs a user-controllable worker count and a sane default.
- **Determinism.** The dev suite requires deterministic output; process pools
  must not perturb RNG seeding or result ordering (results are recombined by
  detector key, so ordering is fine; seeds are per-call).
- **Picklability / logging.** Worker results and the per-detector inputs must
  pickle; the per-process logger needs to funnel back to the main log.
- **A new parameter.** Likely `par['rdx']['ncpu']` (or `nproc`), default `1` to
  preserve today's behavior exactly, opt-in to >1. (See Q&A.)

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

Tracked in the prompt doc's **Q&A / Planning** section; summarized here:

1. Parallelism backend, default worker count, and the new parameter
   (`ncpu`/`nproc`, default 1?).
2. Parallelism granularity to target first (detector-level recommended).
3. Tolerance for *non-identical* numerical output from vectorized fits.
4. Determinism requirements under multiprocessing.
5. Which instruments/setups to validate against first.
6. Appetite for the cheap QA wins (force `Agg`; option to skip/defer/parallelize
   QA) as an early, independent PR.
