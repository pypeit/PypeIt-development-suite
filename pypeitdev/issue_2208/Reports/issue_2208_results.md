# PypeIt Issue #2208: irreproducible reductions — results

Investigation of [pypeit/PypeIt#2208](https://github.com/pypeit/PypeIt/issues/2208)
on the `profx` workstation, 2026-10-18/19.  PypeIt branch `edge_2208` (off
`develop` at `0fd778f07`), dev-suite branch `dev_edge_2208`.  The full
working log, with every command, is in
`claude_prompts/issues/issue_2208.md`.

## Summary

- **Cause.** `pypeit.core.pca.pca_decomposition` called
  `sklearn.decomposition.PCA(n_components=npca)` with the default
  `svd_solver='auto'` and no `random_state`.  At the array shapes PypeIt
  uses, sklearn resolves `'auto'` to its **randomized** SVD, which is
  seeded from fresh entropy on every call.  So every PCA decomposition
  differed slightly (≈1e-13 to 1e-12) from run to run.
- **Amplification.** In slit-edge tracing, `peak_refine` re-measures the
  edges from PCA-predicted positions.  That turns the ~1e-13 difference into
  ~1e-6 px in `EDGE_CEN`, then through Slits → Flat → Tilts → sky
  subtraction into large differences in low-S/N extracted pixels.
- **Fix.** Use `svd_solver='full'` (an exact LAPACK SVD) in both `PCA` calls
  of `pca_decomposition`.  It is a two-line change with no measurable cost.
- **Result.** Two cold, full, serial reductions are now **bit-identical in
  every calibration and science file** for LRIS slitmask, DEIMOS 830G and
  Binospec long slit (the latter was the second, tilt-tracing, entry point).
  A DevSuite subset of 16 setups passes all 29 reduce/afterburn tests and 18
  vet tests.

## What #2208 got right, and what it did not

| #2208 statement | finding |
|---|---|
| `SOBELSIG` is bit-identical; `EDGE_CEN` differs by 1 float32 ULP | **Confirmed.**  The in-memory float64 values actually differ in ~87 % of elements, by ~1e-8 px (median) to ~1e-6 px (max).  The float32 file shows only the elements that cross a rounding boundary. |
| Not threads, not `PYTHONHASHSEED`, not the B-spline | **Confirmed.** |
| Item 5: sklearn's PCA solver resolves to `full`/`covariance_eigh` | **Wrong for the second call.**  `PCA()` (line 92) resolves to `full`, but `PCA(n_components=npca)` (line 120) resolves to `randomized` for every edge PCA in this data (21×2048, 28×2048). |
| "Where to look": `trace.masked_centroid` | **Two steps too early.**  `centroid_refine` and `fit_refine` are deterministic; the first divergence is `build_pca` inside `pca_refine`. |
| Hypothesis: SIMD results depending on memory alignment | **Ruled out.**  Outputs differ between repetitions in one process, and are independent of byte offset and of the OpenBLAS kernel type. |

## Evidence

### 1. Baseline (Testing #1)

Two cold `run_pypeit -c` runs of `keck_lris_blue/multi_600_4000_slitmask`,
with BLAS/OpenMP pinned to one thread:
- Arc, Tiltimg and WaveCalib are identical.  Every file from Edges onward
  differs.
- `EDGE_CEN`: 66 / 45 056 elements differ on DET01 (max 1 ULP), and
  1 393 / 57 344 on DET02 (max 14 ULP).  The DET02 outlier, trace 0, sits at
  x ≈ 0.5–11 px, where the float32 ULP is tiny; it is not a separate
  mechanism.
- `OPENBLAS_CORETYPE=Haswell` (AVX2 instead of AVX-512) gives the same
  irreproducibility.

### 2. Stepwise harness (Testing #2)

`edge_harness.py` re-runs `EdgeTraceSet.auto_trace` from the saved
`TraceImage`, hashing the in-memory state after each of its 96–110 method
calls.  It reproduces `run_pypeit`'s run-to-run behaviour.

- 10 runs (5 in one process, 5 in fresh processes) give 10 different results.
- The first difference is **call #25, `build_pca` in `pca_refine`**.  Its
  outputs (`pca_coeffs`, `pca_components`) differ while its inputs
  (`trace_coo`, `pca_mean`) are identical.

### 3. Mechanism (Testing #3)

- sklearn 1.9.0 (`_pca.py:524-536`) picks `randomized` when the array is not
  tall-and-skinny, `max(shape) > 500`, and `1 ≤ n_components < 0.8·min(shape)`.
- Captured calls: all 5 `pca_decomposition` calls on the two LRIS detectors
  run `['full', 'randomized']`.
- Replay of the captured inputs (`replay.py`; 3 calls × 8 byte offsets, in 2
  processes):

  | solver | distinct outputs | across processes | ms/call |
  |---|---|---|---|
  | unpatched (`auto` → randomized) | 24 / 24 | different | 3.3–4.3 |
  | `svd_solver='full'` | 1 / 24 | identical | 3.1–4.3 |
  | `random_state=0` | 1 / 24 | identical | 3.5–4.3 |

  The pattern is the same under the SkylakeX, Haswell and Prescott kernels.
  With `full`, the result is also identical with unpinned (24) BLAS threads.
- `pca_decomposition` is PypeIt's only use of sklearn.  It is reached from
  edge tracing (`EdgeTraceSet.build_pca`), from tilt tracing
  (`core/tracewave.py`, always run, with `npca=2`), and from the echelle
  object PCA (`core/findobj_skymask.py`).  Tilt tracing hits `randomized`
  when a slit is > 500 px wide, as for Binospec long slit (~2 460 px);
  that is why Binospec first diverged in `Tilts` (#2209).

### 4. Size of the change

The fix moves the edges by about as much as two unpatched runs already
differ from each other:

| `edge_cen` (float64), max / median \|d\| (px) | DET01 | DET02 |
|---|---|---|
| unpatched run vs unpatched run | 3.3e-6 / 3.0e-8 | 8.8e-7 / 1.5e-8 |
| fixed vs unpatched | 4.3e-6 / 2.7e-8 | 1.2e-6 / 2.8e-8 |

## Fix and test (Testing #4)

- `pypeit/core/pca.py`: `PCA(svd_solver='full')` and
  `PCA(n_components=npca, svd_solver='full')`, with a code comment and a
  docstring note.
- `pypeit/tests/test_pca.py::test_decomposition_reproducible`:
  - It uses a fixed-seed (28, 2048) edge-like input and first asserts that
    `'auto'` would select `randomized` for it, so a future change in
    sklearn's policy cannot make the test pass vacuously.
  - It then asserts bit-identical output across 12 calls at 4 byte offsets.
  - It **fails** on the old code and passes with the fix.
- Full `pypeit/tests`: **877 passed**, 0 failed, 0 skipped.
- `doc/releases/2.1.0dev.rst`: a Bug Fixes entry and a Testing entry.

## Validation

### End-to-end pairs (Testing #5)

Two cold, full, serial `run_pypeit` reductions per setup:

| setup | files differing | science arrays compared | wall (each) |
|---|---|---|---|
| keck_lris_blue/multi_600_4000_slitmask | **0 / 16** | 852 | 6.2–6.4 min |
| keck_deimos/830G_M_8500 | **0 / 16** | 1 528 | 34.1–34.4 min |
| mmt_binospec/Longslit_G1000 | **0 / 16** | 398 | 24.1–25.0 min |

Before the fix: LRIS had 637 of 842 science arrays differing (#2208), and
Binospec had 235 of 392 (#2209).

### DevSuite subset (Testing #6)

16 setups at `-t 4` (161 min):
- multi-slit/slit-mask: LRIS blue/red/red-mark4, DEIMOS ×3 including
  `600ZD_M_6500`, MOSFIRE, GMOS, Binospec;
- Binospec long slit;
- echelle: ESI, MagE, X-Shooter VIS, NIRES;
- Kast blue as a long-slit control.

Results:
- `pypeit_test reduce afterburn`: **29 / 29 passed**.
- Vet: **18 / 18 passed**: `test_slitmask` (6), `test_edgetrace` (2),
  `test_findobj` (2), `test_wavetilts` (3), `test_skysub`, `test_fluxspec`
  (2), `test_collate1d`, `test_echelles::test_magellan_mage_orders`.
- Not run: `test_edgetrace::test_chk_*_slits` (needs `$PYPEIT_DEV/Cooked`),
  and vet tests whose inputs were outside the subset.

## Not covered

- **macOS/Accelerate and Windows** were not tested.  The mechanism (an
  unseeded RNG) is platform-independent, so the fix should apply there too,
  but it is unverified.
- **`--ncpu` (#2209):** `edge_2208` does not contain the detector
  parallelism, so serial-vs-parallel equality was not re-measured with the
  fix.
- **The full DevSuite** (`all`) was not run; only the 16-setup subset.
- **Bit-identity across machines or BLAS kernels is not expected**: e.g.
  Haswell vs SkylakeX kernels already change the processed images.  The
  guarantee is run-to-run on one machine and software stack.

## Tools (`pypeitdev/issue_2208/`)

| file | purpose |
|---|---|
| `compare_calibs.py` | array-by-array comparison of two reductions (counts, max \|d\|, ULPs) |
| `edge_cen_diff.py` | per-trace `EDGE_CEN` differences between two Edges files |
| `run_calib.sh` | cold `run_pypeit -c` of LRIS slitmask with threads pinned |
| `edge_harness.py` | stepwise `auto_trace` harness with state hashes (`--pca_fix`, `--capture`) |
| `pca_patch.py`, `site/sitecustomize.py` | runtime patch/capture of `pca_decomposition` (no PypeIt edits) |
| `replay.py` | replay of captured PCA calls at chosen byte offsets |
| `run_pairs.sh` | paired cold full reductions plus comparison |
| `run_devsuite_2208.sh` | the DevSuite subset run and vet selection |
