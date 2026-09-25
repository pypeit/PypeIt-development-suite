# PypeIt speed-up — measured results

Running record of the before/after wall-clock numbers for the speed-up
workstreams (design: `speed_up_design.md` v0.5; implementation:
`speed_up_coding.md`).  One section per landed workstream.

All times are the cProfile-instrumented wall-clock of a **cold** run
(`run_pypeit <file> -o`, single run, artifacts in this directory).  Cold-run
methodology as in the two baseline profile reports.

## Baseline (June 2026 profiles)

Environment: `pypeit` conda env (numpy 1.x-era stack), `speed_up` branch at
`e9ed85c1a`, interactive `qtagg` matplotlib backend (pre-`Agg` fix).

| Setup | Wall (s) | QA PNGs |
|---|---:|---:|
| shane_kast_blue 600_4310_d55 (1 detector) | 121.2 | 20 |
| keck_deimos 600ZD_M_6500 (4 mosaics, slitmask) | 12 393.8 | 971 |

Full analyses: `shane_kast_blue_A_profile_report.md`,
`keck_deimos_600zd_m_6500_profile_report.md`.

## After workstream A (QA cheap wins + `ncpu` plumbing) — 2026-09-21

Branch `speed_up_qa` at `c39818589`; env **`pypeit14b`** (numpy 2.5.0,
astropy 8.0.0, matplotlib 3.10.3, python 3.14) — the June `pypeit` env can no
longer run (numpy/astropy incompatibility), so the `--ncpu 1` runs below are
the **fresh baselines**; the June numbers are not directly comparable (they
differ by env, branch drift, and the pre-`Agg` backend).

Artifacts: `*_ncpu{1,4}.{prof,pstats.txt,timeline.txt,runmeta.json,run.log}`,
QA spot-check samples in `task6_qa_samples/`, PNG counts in
`task6_png_counts.txt`.

| Setup | `--ncpu 1` (s) | `--ncpu 4` (s) | Δ | QA PNGs (both) |
|---|---:|---:|---:|---:|
| shane_kast_blue | 88.5 | 82.6 | **−5.9 s (−6.7%)** | 20 = 20 ✓ |
| keck_deimos | 10 056.9 | 9 834.0 | **−222.9 s (−2.2%)** | 971 = 971 ✓ |

Both runs completed with return code 0 and healthy outputs; sample QA PNGs
(first/middle/last, per run) decode correctly and threaded output is
pixel-identical to serial in the unit/integration tests.

### Attribution

- **`Agg` backend + env/branch refresh** (already in the fresh `ncpu=1`
  baseline): Kast 121.2 → 88.5 s, DEIMOS 12 393.8 → 10 056.9 s.  This bundles
  the `Agg` fix (`bbfe2675a`, expected 1–3%) with the numpy-2.5/python-3.14
  env change and ~3 months of branch drift — the individual contributions are
  not separable without re-running the old env, which no longer works.
- **Threaded QA writes (`--ncpu 4` vs `--ncpu 1`, same code, same env)** — the
  clean PR-A measurement:
  - DEIMOS −222.9 s.  The profile confirms the mechanism: PNG encode
    (`ImagingEncoder.encode`) is 275.6 s of main-thread self-time at `ncpu=1`
    and is moved to worker threads at `ncpu=4` (the small residue is the
    unconverted/PdfPages writers, which stay synchronous by design).
  - Kast −5.9 s (encode budget ~4.9 s at `ncpu=1` plus write overlap).

### Design deviation found during validation (recorded in A.md Q&A)

The original design handed the full `Figure.savefig` (render + encode) to
worker threads.  The DEIMOS `--ncpu 4` run **crashed** within minutes:
matplotlib's mathtext parser is process-global and not thread-safe, and
worker-thread rendering races the main thread's figure layout
(`ParseFatalException: Unknown symbol: \mathdefault` on log-axis tick
labels).  Fix (commit `c39818589`): `qa.save_figure` renders on the main
thread and defers only the **PIL PNG encode** (thread-safe, releases the
GIL).  This lands exactly the "encode half" that the accepted planning answer
Q2 anticipated as the recoverable fraction; per Q2, the machinery is kept and
QA threading is not escalated further in this PR.  **The design change was
accepted in A.md Q&A (Q2, answered 2026-09-22: "Let's follow your
recommendation").**

### Post-review refactor (2026-09-24) — performance-neutral

Review of PR #2198 asked that the deferred-write machinery follow the
package's `log`/`dataPaths` convention rather than living in module-level
globals, so `qa.{init_qa_pool,save_figure,flush_qa}` became
`pypeit.qaWriter.{init,save_figure,flush}` (a `QAWriter` instance created on
import; commit `b8e5f9321`).  The write path itself is unchanged, and the
Kast pair was re-run to confirm it:

| Setup | `--ncpu 1` (s) | `--ncpu 4` (s) | Δ |
|---|---:|---:|---:|
| shane_kast_blue, before refactor | 88.5 | 82.6 | −5.9 (−6.7%) |
| shane_kast_blue, after refactor | 89.8 | 85.2 | −4.6 (−5.2%) |

The ~1 s differences are within the run-to-run scatter of a single cold run;
20 QA PNGs in every case.  **The DEIMOS numbers above were not re-measured**
— they predate this pure refactor, which touches no per-figure work beyond an
attribute lookup.

### Verdict vs the PR-A target

The coding-doc §A.7 combined target (DEIMOS −4–7%, Kast −4–9% from `Agg` +
threading) is **met or exceeded relative to the June baselines** (DEIMOS
−20.7% to the fresh `ncpu=4` number, Kast −31.8%, albeit with env/branch
confounds), while the *isolated* threading gain (−2.2% DEIMOS, −6.7% Kast) is
at the lower edge of the §A.7 threading range, consistent with the Q2
"encode-half" expectation.  Per Q2: measurement recorded, machinery kept; the
render-side QA cost (matplotlib text metrics) remains main-thread and is a
candidate to revisit only if QA is still a top-5 hotspot after PR B.

## After workstream B (detector parallelism v1)

*(pending)*

## After workstream C (arc-fit vectorization)

*(pending)*
