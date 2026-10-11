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

## After workstream B (detector parallelism v1) — 2026-10-10

Branch `speed_up_B`; env `pypeit14b` (numpy 2.5.0, astropy 8.0.0, python 3.14),
24 cores / 755 GB.  Artifacts: `keck_deimos_600zd_m_6500_ncpu{1,4,4b}.{runmeta,rss}.json`,
`shane_kast_blue_A_ncpu{1,2}.runmeta.json`.

Stages 1 (calibration build), 2 (`process_exposure`) and 4 (`extract_exposure`)
are distributed one worker process per detector; stage 3 (object finding +
global sky) stays serial behind its cross-detector slitmask barrier.

### Wall-clock

All DEIMOS numbers below are **uninstrumented** wall-clock, so the ratio is
like-for-like.  Kast is cProfile-instrumented on both sides (also like-for-like).

| Setup | `--ncpu 1` | `--ncpu 4` | speedup |
|---|---:|---:|---:|
| keck_deimos 600ZD_M_6500 (4 mosaics) | **9 360.0 s** | **4 123.6 s** | **2.27x** |
| keck_deimos, repeat at `--ncpu 4` | — | 4 244.5 s | 2.21x |
| keck_deimos, mean of the two `--ncpu 4` runs | 9 360.0 s | 4 184.0 s | **2.24x** |
| keck_lris_blue multi_600_4000_slitmask (2 det) | 389.0 s | 250.7 s (`--ncpu 2`) | **1.55x** |
| shane_kast_blue (1 det) | 109.7 s | 102.5 s (`--ncpu 2`) | 1.07x |

For reference, the same DEIMOS `--ncpu 1` reduction **under cProfile** takes
10 515.5 s, i.e. the instrumentation costs **+1 155 s (+12.3%)**.  An earlier
version of this section compared that instrumented baseline against the
uninstrumented parallel runs and reported **2.55x**; that was inflated, and the
tell was that it exceeded the Amdahl ideal.  The corrected figure is **2.24x**.

Run-to-run scatter at `--ncpu 4` on DEIMOS is 120.9 s (2.9%).  The implied
parallel fraction is **p = 0.737**, close to the 0.77 assumed in coding-doc
§B.7, and the measured 2.24x sits just below the 2.37x Amdahl ideal for that p
— which is the physically sensible ordering, and a useful check that the
measurement is now self-consistent.  **The §B.7 expectation was a realistic
1.9-2.3x, so this lands inside the predicted band.**

`shane_kast_blue` is single-detector and gains nothing by design; §B.7 only
required that it not regress, and it does not — it is 6.6% *faster* at
`--ncpu 2`, because PR A's QA threading more than pays for the fork.

### Peak memory

Measured as the peak **PSS** summed over the whole process tree
(`scripts/run_with_rss.py`).  PSS rather than RSS because the workers are forked
and share a large amount of copy-on-write memory, which summing RSS would
double-count; and a tree sum rather than `/usr/bin/time -v`, whose `ru_maxrss`
is the maximum of any *single* process and so under-reports a 4-worker run by
roughly 4x.

| DEIMOS | peak PSS | peak RSS sum | processes |
|---|---:|---:|---:|
| `--ncpu 1` | 11.70 GB | 11.69 GB | 1 |
| `--ncpu 4` | **39.38 GB** | 45.03 GB | 5 |

Memory grows **3.4x for 4 workers** — sub-linear, because the forked children
share the parent's read-only pages.  The practical implication documented in
`doc/running.rst` holds: on DEIMOS-class data, memory is the binding constraint
on `ncpu`, and a 32 GB machine should use `--ncpu 2`.

### Cost of the thread pinning

To make the products independent of `ncpu`, both the serial and the parallel
paths pin BLAS/OpenMP to one thread (see the `doc/running.rst`
"Reproducibility" section).  Measured cost at the default `ncpu = 1`:

| Setup | unpinned | pinned | cost |
|---|---:|---:|---:|
| shane_kast_blue (raw, no cProfile) | 77.5 s | 86.3 s | +11.4% |
| keck_deimos (vs the post-A baseline, both instrumented) | 10 056.9 s | 10 515.5 s | **+4.6%** |

DEIMOS pays less than Kast proportionally, because more of its time is in
b-spline fits that were never BLAS-threaded.  The DEIMOS figure also carries
some branch drift (`develop` was merged since the post-A measurement).

### Identical-output and determinism: PypeIt is not reproducible run-to-run

Both checks **fail on DEIMOS**, and the control shows this is **not caused by
the parallelism**:

| comparison | arrays differing |
|---|---|
| DEIMOS `--ncpu 1` vs `--ncpu 4` | 2 039 / 2 624 |
| DEIMOS `--ncpu 4` vs `--ncpu 4` (repeat) | **2 033 / 2 624** |
| LRIS serial vs `--ncpu 2` | 637 / 842 |
| **LRIS serial vs a second identical serial run** | **637 / 842** |
| LRIS serial vs serial, all thread env vars = 1 *and* `PYTHONHASHSEED=0` | **637 / 842** |
| shane_kast_blue serial vs `--ncpu 2` (three times) | **0 / 328** |

Two runs of the *same command* disagree by the same amount as a
serial-versus-parallel pair, and pinning threads and fixing the hash seed change
nothing.  `shane_kast_blue`, where PypeIt *is* deterministic, is bit-identical
between `--ncpu 1` and `--ncpu 2`.

The source localises to **edge centroiding**, not — as first suspected — to the
B-spline.  Walking the calibration chain of two identical serial LRIS runs in
dependency order:

```
Arc_A_0_DET01.fits         15 arrays  identical
Tiltimg_A_0_DET01.fits     15 arrays  identical
WaveCalib_A_0_DET01.fits  173 arrays  identical
Edges_A_0_DET01.fits.gz    39 arrays   13 DIFFER   <- origin
Slits_A_0_DET01.fits.gz    32 arrays    7 DIFFER
Flat_A_0_DET01.fits        42 arrays   30 DIFFER
Tilts_A_0_DET01.fits       14 arrays    3 DIFFER
```

Inside `Edges`, the input is bit-identical and only the measurement changes:
`TRACE_IMAGE`, `TRACE_FULLMASK`, `TRACEBPM`, `SOBELSIG`, `TRACEID`,
`MASKDEF_ID` and `EDGE_MSK` all match exactly, while **`EDGE_CEN` differs by
exactly one float32 ULP (6.104e-05)** — in all ten pairwise comparisons of five
runs, including pinned-versus-pinned.  From there it propagates to `EDGE_FIT`,
the trace PCA, `Slits`, the flat and the tilts.

**The B-spline was a red herring** and is explicitly exonerated:
`iterative_bspline_fit` on fixed synthetic data with rejection active returns
byte-identical coefficients, breakpoints and fit values across repeats both
within one process and across separate processes.  The B-spline fits stored in
`Flat` differ only because their inputs already differ.  Also ruled out: BLAS
thread count, `PYTHONHASHSEED`, unseeded RNG (there is none in the reduction
path), and sklearn's randomized SVD in the trace PCA (`auto` resolves to `full`
or `covariance_eigh` at PypeIt's shapes — both deterministic).

Filed upstream as **pypeit/PypeIt#2208**, which is also the rigorous retest of
**#1986** that was requested there after the C B-spline code was removed.

**Consequence:** the locked "identical to machine precision" requirement is
unachievable for slitmask data by any implementation of PR B, and it is not a
property this PR can deliver or restore.  See **Q11** in
`claude_prompts/speed_up/B.md` for the options; the B-spline non-determinism
should be filed as its own issue, with "two identical serial runs of the same
`.pypeit` file disagree" as the reproducer.

### Verdict vs the PR-B target

Wall-clock: **met** — 2.24x on DEIMOS at `--ncpu 4`, inside the 1.9-2.3x band
§B.7 predicted, measured uninstrumented on both sides.
Memory: measured and documented, and better than the linear scaling the design
assumed.  Identical output: **not met, and not attainable** — but the
parallelism is exonerated by the serial-versus-serial control, and is exact on
the setup where PypeIt itself is exact.

## After workstream C (arc-fit vectorization)

*(pending)*
