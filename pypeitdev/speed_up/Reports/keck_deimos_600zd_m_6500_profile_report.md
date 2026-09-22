# Profiling report — `keck_deimos_600zd_m_6500` (cold run)

**Date:** 2026-06-28
**Author:** JXP and Claude
**PypeIt version:** 1.17.5.dev1018+g42ed28072 (branch `speed_up`)
**Dataset:** `RAW_DATA/keck_deimos/600ZD_M_6500` (MultiSlit with a slitmask —
many slits/detector; `ndet = 8` paired into the default **4 mosaics**
`[(1,5),(2,6),(3,7),(4,8)]`). A much heavier target than shane_kast_blue.
**Command profiled:** `run_pypeit keck_deimos_600zd_m_6500.pypeit -o`

---

## 1. Method

Same methodology as the Kast report (per Q3), extended with a per-mosaic
breakdown.

- **Cold run** (per Q2). The redux dir was clean apart from the `.pypeit` and a
  stale `.test.log`, so the cold-clean was effectively a no-op; the full
  calibration build is included in the timing.
- **cProfile**, in-process, single run (per Q4). The keck_deimos MultiSlit path
  uses **no multiprocessing** and has no `nproc` parameter (the Q7 finding is
  global), so the profile captures the entire workload — no child-process blind
  spots.
- A **per-phase wall-clock timeline** was reconstructed from the timestamped run
  log, plus a **per-mosaic** breakdown from the "Calibrating detector (a, b)"
  markers.
- One science frame (`d1010_0056`, 1800 s; `0057` commented out), the full
  4-mosaic default (`detnum` commented out in the `.pypeit`).
- Code: `pypeitdev/speed_up/scripts/profile_deimos.py` (clean + run) and
  `analyze_profile.py --stem keck_deimos_600zd_m_6500` (pstats tables + timeline
  + per-mosaic). Raw artifacts in this folder:
  `keck_deimos_600zd_m_6500.prof`, `.pstats.txt`, `.timeline.txt`, `.run.log`,
  `.runmeta.json`, `.driver.out`.

**Reference numbers.** PypeIt wall-clock (runmeta): **12 393.8 s ≈ 3 h 26 m**;
log span 12 390 s; cProfile-instrumented total 12 389.7 s; **4.10 billion**
function calls. This is ~104× the Kast runtime (119 s). As with Kast, `cumtime`
values overlap when a shared kernel (B-splines) is called from several phases —
they show *where* time goes, not a partition that sums to the wall-clock.

This run was **healthy**: 1 `spec1d` + 1 `spec2d` written for the science frame,
calibrations built for all 4 mosaics (Arc/Bias/Edges/Flat/Tilts/WaveCalib ×4),
971 QA PNGs. The 8 "error" matches in the log are benign INFO lines
(`Max centroid error: None`), not failures. No exceptions; rc 0.

---

## 2. Headline result

> **The same pure-Python/NumPy B-spline kernel in `pypeit/core/bspline/` is again
> the dominant engine — ≈2 750 s of self-time, ~22 % of the entire 12 390 s
> runtime** — reached through sky subtraction and extraction. The Kast finding
> (≈32 % there) holds at scale and is the highest-leverage single change, helping
> *every* spectrograph because sky/extract are common to all paths.

B-spline self-time breakdown (cProfile `tottime`):

| Function (`pypeit/core/bspline/`)        | self-time | ncalls |
|------------------------------------------|----------:|-------:|
| `util.py:cholesky_band`                  |  856.3 s  |  6,818 |
| `util.py:solution_arrays`                |  654.4 s  |  5,463 |
| `util.py:cholesky_solve`                 |  437.2 s  |  5,276 |
| `util.py:intrv`                          |  312.3 s  |  6,663 |
| `bspline.py:bsplvn`                      |  252.2 s  |  6,663 |
| `util.py:bspline_model`                  |  237.9 s  |  9,170 |
| **B-spline core total (self)**           | **≈2 750 s** |     |

Driven from the top by `core/fitting.py:bspline_profile` (cum **3 437 s**) →
`extraction.local_skysub_extract` (cum **2 952 s**) and
`find_objects.global_skysub` (cum **1 677 s**). These hand-written
linear-algebra/interpolation loops are pure Python + small NumPy ops with high
call counts — the canonical profile for Cython/Numba/vectorization.

---

## 3. Time by workflow phase

Aggregate wall-clock per phase from the log timeline (mapped onto
`pypeit_workflow.md` §3). This is a wall-clock partition (sums to ~12 390 s),
distinct from the overlapping cProfile `cumtime`.

| Phase (workflow §)                         | Wall (s) | % of run |
|--------------------------------------------|---------:|---------:|
| **Sky subtraction** (local + global, §3.3) | **4 284** | **35 %** |
| **Flat field** (§3.5)                      | 2 173    | 18 % |
| **Slit edges** (§3.5)                      | 1 896    | 15 % |
| **Extraction** (§3.3)                      | 1 097    | 9 % |
| **Tilts** (§3.5)                           | 779      | 6 % |
| **Wavelength** (§3.5)                      | 771      | 6 % |
| **Object finding / tracing** (§3.3)        | 728      | 6 % |
| **Image combine** (§3.5)                   | 299      | 2 % |
| Calib build/load + finalize + misc         | ~355     | 3 % |
| Initialization (§3.2)                      | ~7       | <1 % |

**Sky subtraction + extraction ≈ 5 380 s (~44 % of wall-clock)** — the
B-spline-heavy half — while calibration building (flat + edges + tilts + wave +
combine ≈ 5 920 s) is the other ~48 %. The slitmask (many slits/detector) and
4 mosaics inflate both halves relative to single-slit Kast.

---

## 4. Secondary hotspots — new at DEIMOS scale

These are negligible on single-detector Kast but substantial here:

- **scipy.ndimage / mosaic resampling.** `_nd_image.correlate` 577 s,
  `geometric_transform` (affine) 462 s, `rank_filter` (median, CR) 434 s,
  `correlate1d` 383 s, `spline_filter1d` 122 s. `core/mosaic.build_image_mosaic`
  cum 662 s and `rawimage.build_mosaic` cum 678 s — resampling 8 chips into 4
  mosaics is a DEIMOS-specific tax absent on Kast.
- **Arc / wavelength fitting churn.** `core/fitting.fit_gauss` cum 840 s over
  **1.11 M** `scipy curve_fit` calls; `gauss_3deg` self 332 s over **93 M**
  calls; `curve_fit` cum 648 s; `leastsq` cum 586 s. Arc-line Gaussian fitting
  via per-line `curve_fit` is an enormous call-count hotspot
  (`detect_lines`/`fit_arcspec` underneath wavecal + tilt tracing).
- **`core/moment.moment1d`** self 212 s / cum 799 s over **1.05 M** calls — the
  object-trace centroiding loop (`trace.masked_centroid`/`follow_centroid`).
- **Flexure.** `spec_flexure_slit` cum 825 s (`measure_fwhm` →
  `get_archive_spectrum`) — a sizeable new cost not seen on Kast.
- **QA / plotting.** matplotlib `_get_text_metrics_with_cache_impl` cum 894 s;
  PNG `ImagingEncoder.encode` 276 s (971 PNGs); `zlib` compress 73 s. Reinforces
  the Kast "force `Agg`, make QA optional/deferred/parallel" win — larger here.
- **Generic NumPy churn.** `ndarray.flatten` 351 s (80.7 M calls),
  `ufunc.reduce` 312 s (233 M calls), `ndarray.copy` 152 s, masked-array
  `ma.core.__call__` 94 s (4.07 M calls), `numpy.array` 111 s (27.3 M calls) —
  death by a billion small array ops in the fitting/rejection inner loops.
- **gzip.** `zlib decompress` 154 s reading the large `.fits.gz` raw frames
  (42–124 MB each).

---

## 5. Per-mosaic breakdown (with caveat)

| mosaic | start | wall (s) |
|--------|-------|---------:|
| (1, 5) | 05:22:48 | 1 338 |
| (2, 6) | 05:45:06 | 1 584 |
| (3, 7) | 06:11:30 | 1 472 |
| (4, 8) | 06:36:02 | 7 989 |

**Caveat:** the breakpoints come from the "Calibrating detector" markers, which
cleanly bound the **calibration-building** phase only. Calibration building is
fairly even per mosaic (~1 300–1 600 s each). The `(4,8)` row reads 7 989 s
because the entire post-calibration **science reduction loop** (object find + sky
+ extract, run after all four mosaics are calibrated) falls after the last
marker and is lumped into it. So the science half (~6 200 s) is shared across
mosaics, **not** a (4,8)-specific cost. A cleaner per-mosaic science split would
need finer instrumentation (e.g. a marker at each detector's reduce entry).

---

## 6. Speed-up opportunities (ranked, for the design doc)

1. **Optimize the B-spline kernel (`pypeit/core/bspline/`)** — *highest leverage,
   ~2 750 s / 22 %.* Cythonize / Numba-JIT / vectorize `cholesky_band`,
   `solution_arrays`, `cholesky_solve`, `intrv`, `bsplvn`, `bspline_model`.
   Accelerates global sky, local sky, and extraction at once, on every
   spectrograph. As noted in the Kast report, these run as pure Python today
   (status quo, not a refactor regression) and a stale `_bspline.*.so` shows they
   were historically compiled — strong precedent for re-introducing a build.
2. **Mosaic resampling (DEIMOS/multi-detector-specific)** — ~1 100 s+ in
   scipy.ndimage `geometric_transform`/`correlate`/`spline_filter1d` via
   `build_image_mosaic`. Cache spline filters, reduce redundant resampling, or
   use a cheaper interpolation order where accuracy permits.
3. **Arc-line fitting (`fit_gauss`/`curve_fit`)** — 1.11 M `curve_fit` calls,
   ~840 s cum. Replace per-line scipy `curve_fit` with a vectorized/analytic
   Gaussian fit (3-param moment or batched least-squares); huge call-count win,
   benefits wavecal + tilt tracing.
4. **Headless / optional / parallel QA** — force matplotlib `Agg`; defer or
   parallelize the 971 PNGs (~900 s of text-metrics + encode).
5. **`moment1d` centroiding** — 1.05 M calls, ~212 s self; vectorize or Cythonize
   the centroid loop.
6. **Flexure (`spec_flexure_slit`)** — ~825 s; cache/limit archive-spectrum FWHM
   measurement.
7. **Reduce NumPy/masked-array allocation churn** in fitting/rejection inner
   loops (preallocate, in-place, avoid masked-array overhead).
8. **Detector/mosaic- and frame-level parallelism** — the pipeline is fully
   serial. Unlike single-detector Kast, DEIMOS has 4 independent mosaics whose
   calibration *and* reduction could run concurrently — the obvious structural
   win for DEIMOS/Echelle/IFU, and it should be designed into the dashboard/
   refactor.

Items 1–4 target roughly half the runtime; the B-spline win (item 1) generalizes
across all spectrographs, while items 2–3 are where DEIMOS-class, multi-detector
slitmask setups pay extra over the simple Kast case.

---

## 7. Reproduce

```bash
# pypeit env python; from this scripts dir
cd PypeIt-development-suite/pypeitdev/speed_up/scripts
/home/xavier/miniconda3/envs/pypeit/bin/python profile_deimos.py            # cold run + cProfile (~3.5 h)
/home/xavier/miniconda3/envs/pypeit/bin/python analyze_profile.py \
    --stem keck_deimos_600zd_m_6500                                          # pstats tables + timeline + per-mosaic
```

Artifacts written to `pypeitdev/speed_up/Reports/`:
`keck_deimos_600zd_m_6500.prof`, `.pstats.txt`, `.timeline.txt`, `.run.log`,
`.runmeta.json`, `.driver.out`, and this report.
