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

### Dev-suite validation (2026-10-12)

`./pypeit_test reduce -i shane_kast_blue keck_deimos keck_lris_blue -t 6`:
**27 passed, 2 failed of 29**, 2h53m.

The `--ncpu 2` regression on `keck_lris_blue/multi_600_4000_slitmask` ran for
the first time and passed, and its log confirms it genuinely used workers rather
than falling back — `Processing 2 detector(s) using 2 worker process(es)` three
times, once per parallelised stage, against zero such lines in the serial run.
`vet_tests/test_ncpu.py` gives 4 passed, 1 xfailed (the exact flux comparison
keyed to Issue #2208) against those products.

The 2 failures, `shane_kast_blue/452_3306_d57` and `830_3460_d46`, are stale
dev-suite state and not regressions: their `Science/` directories held
June-2025 products carrying superseded datamodel versions (`DMODVER` 1.0.1 /
1.1.0 / 1.1.5), written before the `SlitTraceSet` bump to 1.1.6 that landed on
`develop` in `5d5e5d75d`.  Neither `-o` nor `-m` clears them, because
`save_exposure` opens the existing spec2d to merge detectors into it.  With
those directories moved aside the two setups pass **2/2**.  `pypeit/slittrace.py`
is not in PR B's diff.

`pytest vet_tests` over the whole `REDUX_OUT` reports 71 failures, which is an
artefact of scope rather than a result: that suite expects `./pypeit_test all`,
whereas only the **reduce** phase for **3 of the 8** instruments present was
run.  Scoped to the modules vetting those products, the result is 10 passed,
1 xfailed and 3 failed, and all three failures need instruments outside the run
(`bok_bc`, an echelle setup, and `keck_lris_red`).

### Verdict vs the PR-B target

Wall-clock: **met** — 2.24x on DEIMOS at `--ncpu 4`, inside the 1.9-2.3x band
§B.7 predicted, measured uninstrumented on both sides.
Memory: measured and documented, and better than the linear scaling the design
assumed.  Identical output: **not met, and not attainable** — but the
parallelism is exonerated by the serial-versus-serial control, and is exact on
the setup where PypeIt itself is exact.

## After the B port (`forkserver`/`spawn`, per Tim Pickering's review) — 2026-10-18

Branch `speed_up_B` with the port to `initargs` delivery (`forkserver` on
Linux, `spawn` elsewhere), on top of #2211's memory reductions, which were
merged into the branch after the `fork` measurements above.  Env `pypeit14b`
(Python 3.14.6), same 24-core / 755 GB workstation.  Every run is cold, into
its own directory, inside `systemd-run --scope -p MemoryMax=80G`.  Artifacts
are in `REDUX_OUT/_port_validation/` (one directory per run, with `rss.json`,
`run.log`, `stderr.txt`).  The runs used `instr_run_pypeit.py`, a `run_pypeit`
wrapper that logs each stage's wall time and the size and pickle cost of each
stage-4 payload; it adds about 1 s per detector.

### Wall-clock and peak memory

Peak PSS is summed over the whole process tree (`scripts/run_with_rss.py`).

| Setup | run | wall | speedup | peak PSS |
|---|---|---:|---:|---:|
| keck_lris_blue multi_600_4000_slitmask (2 det) | `--ncpu 1` | 377.8 s | — | 1.25 GB |
| | `--ncpu 1`, repeat | 375.7 s | — | — |
| | `--ncpu 2` forkserver | 240.0 s | 1.57x | 3.08 GB |
| | `--ncpu 2` spawn | 247.6 s | 1.53x | 3.07 GB |
| mmt_binospec Longslit_G1000 (2 det) | `--ncpu 1` | 1 162.3 s | — | 10.05 GB |
| | `--ncpu 2` forkserver | 647.1 s | 1.80x | 20.12 GB |
| | `--ncpu 2` spawn | 656.3 s | 1.77x | 19.43 GB |
| keck_deimos 600ZD_M_6500 (4 mosaics) | `--ncpu 4` forkserver | **3 391.2 s** | 2.76x\* | **36.48 GB** |

\*Against the 9 360.0 s serial baseline above, which predates #2211.  The
`fork` mean on the same data was 4 184.0 s (39.38 GB), so the port is 19%
faster and 7% lighter.  **That gain should not be credited to the start
method**: #2211 (memory and B-spline/fit changes) landed in between, and no
serial DEIMOS run was repeated on the current code.  The like-for-like
comparison of start methods is Binospec and LRIS, where `forkserver` and
`spawn` agree to within 1.5-3% in wall-clock and memory, `spawn` being slightly
slower because each worker imports PypeIt afresh.

Peak memory did **not** rise as feared when copy-on-write sharing was given up.
On the two-detector setups the parallel peak is about 2x the serial one, i.e.
linear in workers; on DEIMOS it is 3.1x serial for 4 workers.

DEIMOS stage times at `--ncpu 4`: calibration build 1 483.2 s, image
processing 78.1 s, extraction 654.6 s; the remaining ~1 175 s is the serial
object finding / global sky stage and I/O.

### Cost of sending the payload (Q4 revisited)

Under `fork` the stage-4 payload was shared copy-on-write; now it is pickled on
submit.  Measured in the parent:

| Setup | payload per detector | `dumps` | `loads` | stage-4 wall |
|---|---:|---:|---:|---:|
| Binospec (2 det) | 977 MB | 0.43-0.46 s | 0.44-0.48 s | 88.6-92.4 s |
| DEIMOS (4 mosaics) | 1 122-1 130 MB | 0.53-0.55 s | 0.31 s | 654.6 s |

That is about 1% of the stage on Binospec and under 0.15% on DEIMOS, far below
the 5% threshold, so there is still no case for writing payloads to disk.

### Workers ran, and nothing was logged twice

Every parallel run logged `Processing N detector(s) using N worker process(es)
(start method: …)` exactly three times (one per parallel stage) with the
expected start method.  For each per-detector line a worker emits
(`Calibrating detector`, `Reducing detector`, `begin log for detector`), the
count on the console equals the count in the log file and equals the serial
run's count: no worker wrote to the inherited stderr.

### Products

`compare_products.py` on the `Science` directories:

| Setup | comparison | arrays differing | worst 99.9th-pct rel. diff (spec1d) |
|---|---|---:|---:|
| LRIS | serial vs serial (control) | **638 / 842** | **0.61** |
| | serial vs `--ncpu 2` forkserver | 638 / 842 | 0.61 |
| | serial vs `--ncpu 2` spawn | 638 / 842 | 0.61 |
| | forkserver vs spawn | 637 / 842 | 0.58 |
| Binospec | serial vs `--ncpu 2` forkserver | 235 / 392 | 3.9e-6 |
| | serial vs `--ncpu 2` spawn | 236 / 392 | 3.0e-7 |
| | forkserver vs spawn | 234 / 392 | 3.8e-6 |
| shane_kast_blue | serial vs `--ncpu 2`, both start methods (`test_run_pypeit_ncpu`) | **0** (bit-identical) | 0 |

The LRIS control is decisive: two identical serial runs differ by exactly as
much as serial vs parallel, so the differences are #2208 (edge centroiding), not
the parallelism.  The large LRIS figure comes from a few pixels of a handful of
low-S/N objects (median S/N 2-3) and is the downstream amplification already
documented in #2208.  `WAVE_RMS` is identical in every comparison.  On the Mac
(PR prompt 6 in B.md) the same LRIS comparisons give only 3e-7, though
`EDGE_CEN` still differs by 1 ULP there; how far #2208 propagates evidently
depends on the platform.

### Dev suite

`./pypeit_test reduce -t 1 -m` on the three `--ncpu` setups (`keck_deimos/830G_M_8500`,
`keck_lris_blue/multi_600_4000_slitmask`, `mmt_binospec/Longslit_G1000`): **6/6
reductions passed**, all under the default start method (`forkserver`).
`vet_tests/test_ncpu.py`: **13 passed, 3 xfailed** (the xfails are the exact-flux
comparisons tied to #2208).  The first vet run had two failures on DEIMOS 830G,
caused by stale state, not the port: the *serial* spec1d had been written into a
`Science/` left over from Task 12, so every object appeared twice.  After that
directory was cleared and the setup re-reduced, both passed.

### Verdict

The port changes no result that the `fork` version did not, works under both
start methods with the same speed-up (to within 3%), keeps memory linear in
workers, and costs about 1% or less to ship payloads.  Like-for-like speed-ups at
`--ncpu 2`: LRIS 1.57x, Binospec 1.80x.

## After workstream C (arc-fit vectorization)

*(pending)*
