# Profiling report — `shane_kast_blue_A` (cold run)

**Date:** 2026-06-27
**Author:** JXP and Claude
**PypeIt version:** 1.17.5.dev1018+g42ed28072 (branch `speed_up`)
**Dataset:** `RAW_DATA/shane_kast_blue/600_4310_d55` (one of the simplest
spectrographs in PypeIt — single detector, single slit, MultiSlit path)
**Command profiled:** `run_pypeit shane_kast_blue_A.pypeit -o`

---

## 1. Method

- **Cold run.** Before profiling, the redux dir was wiped of all products
  (`Calibrations/`, `Science/`, `QA/`, `Intermediate/`), the stale
  `*_state.json`, `*.calib`, prior run log, and old `*_UTC_*.par` dumps, so the
  full calibration build is included in the timing.
- **cProfile**, in-process, single run (per Q4). The shane_kast_blue MultiSlit
  path uses **no multiprocessing** and has no `nproc` parameter, so the profile
  captures the entire workload — no child-process blind spots.
- A **per-phase wall-clock timeline** was also reconstructed from the timestamped
  run log (this branch's `FileFormatter` includes `%(asctime)s`).
- Code: `pypeitdev/speed_up/scripts/profile_kast_blue.py` (clean + run) and
  `analyze_profile.py` (pstats tables + log timeline). Raw artifacts in this
  folder: `shane_kast_blue_A.prof`, `.pstats.txt`, `.timeline.txt`, `.run.log`,
  `.runmeta.json`.

**Caveat on numbers.** cProfile adds ~3% instrumentation overhead. Reference
wall-clock from PypeIt's own clock: **Execution time = 1m 58.72 s (≈119 s)**;
cProfile-instrumented wall = 121.2 s; log span = 120 s. Function `cumtime`
values overlap when a shared kernel (B-splines) is called from several phases —
they attribute *where* time is spent, not a partition that sums to 119 s.

This run was **healthy** (Q5 sanity check): 3 `spec1d` + 3 `spec2d` (1 standard
+ 2 science) written, 20 QA PNGs spanning arc/tilt/wave/flat/object-trace/
profile/flexure. No errors.

---

## 2. Headline result

> **A single numerical kernel — the pure-Python/NumPy B-spline fitting code in
> `pypeit/core/bspline/` — accounts for ≈39 s of self-time, ~32 % of the entire
> 119 s runtime.** It is the shared engine behind both sky subtraction and
> spectral extraction. Speeding up B-splines is the highest-leverage single
> change available.

B-spline self-time breakdown (cProfile `tottime`):

| Function (`pypeit/core/bspline/`)        | self-time | ncalls |
|------------------------------------------|----------:|-------:|
| `util.py:intrv`                          |   10.7 s  |  1,244 |
| `util.py:cholesky_band`                  |    8.7 s  |    618 |
| `bspline.py:bsplvn`                      |    6.8 s  |  1,244 |
| `util.py:solution_arrays`                |    4.6 s  |    289 |
| `util.py:cholesky_solve`                 |    3.9 s  |    586 |
| `util.py:bspline_model`                  |    2.3 s  |  1,348 |
| **B-spline core total (self)**           | **≈39 s** |        |

These are hand-written linear-algebra/interpolation loops; they are pure
Python + small NumPy ops with high call counts — exactly the profile of code
that benefits from Cython/Numba/C or vectorization.

---

## 3. Time by workflow phase

Mapped onto the phases in `pypeit_workflow.md` §3. cumtime = total time spent in
that phase's entry point (inclusive of the shared B-spline kernel).

| Phase (workflow §)                         | Entry point                         | cumtime | % of 119 s |
|--------------------------------------------|-------------------------------------|--------:|-----------:|
| **Extraction** (§3.3, innermost loop)      | `extraction.local_skysub_extract` (×3 frames) | **40.5 s** | **34 %** |
| **Global sky subtraction** (§3.3)          | `find_objects.global_skysub` (×8)   | 17.2 s | 14 % |
| **Build calibration images** (§3.5)        | `buildimage.buildimage_fromlist` (×8: bias/arc/tiltimg/flat raw combine) | 13.6 s | 11 % |
| **Slit edges** (§3.5)                       | `calibrations.get_slits` → `edgetrace` | 8.1 s | 7 % |
| **Per-frame image processing** (§3.3)      | `rawimage.process` (×37)            | 7.2 s | 6 % |
| ↳ cosmic-ray rejection                     | `procimg.lacosmic` / `build_crmask` (×7) | 5.9 s | 5 % |
| **Object finding / tracing**               | `core/moment.moment1d` (×10,414)    | 8.2 s | 7 % |
| **Tilts** (§3.5)                            | `wavetilts.trace_tilts`             | 4.4 s | 4 % |
| **Wavelength** (§3.5)                       | `autoid.full_template`              | 2.2 s | 2 % |
| **Flat field** (§3.5)                       | `flatfield.spatial_fit`             | 1.4 s | 1 % |
| **QA figure generation**                   | matplotlib + PIL + GUI (see §4)     | ≈12 s | 10 % |
| **Initialization** (§3.2)                  | load `.pypeit` + metadata           | ≈1 s | <1 % |

Cross-checking the log timeline: the longest single uninterrupted stretch is the
local-sky + extraction block for the science frames (the 31 s / 5 s / 5 s gaps in
`shane_kast_blue_A.timeline.txt`), consistent with extraction being the dominant
phase. (The timeline's "sky subtraction" aggregate folds the local-sky+extract
block into one label; the cProfile table above is the authoritative split.)

**Reduction (calibrations + objfind + sky + extract) ≈ 90 % of runtime;
calibration *building* ≈ 25–30 %; the standard + 2 science extractions ≈ 35 %.**

---

## 4. Secondary hotspots

- **QA / plotting ≈ 12 s self-time, spread across steps.** matplotlib self-time
  5.9 s, PNG encode (PIL `ImagingEncoder`) 4.8 s (24 images), GUI toolkit
  (`_tkinter` etc.) 1.5 s. **The default matplotlib backend is `qtagg`** — an
  *interactive GUI* backend — even though all QA output is non-interactive PNGs.
  Forcing a headless `Agg` backend would remove GUI-init/toolkit overhead.
- **Generic NumPy churn:** `numpy.ufunc.reduce` 3.2 s (1.88 M calls),
  `bottleneck.nanmedian` 2.7 s, `numpy.linalg.lstsq` 2.9 s, `legendre.legval`
  2.4 s, `numpy.array` 1.7 s (263 k calls), `ndarray.copy` 1.3 s — death by a
  thousand small array allocations in the fitting/rejection loops.
- **Cosmic-ray rejection:** `scipy.ndimage.rank_filter` (median filter) 5.5 s in
  `lacosmic`.
- **GC:** 9 explicit `gc.collect()` calls cost 1.2 s.
- **gzip:** `zlib` decompress 1.1 s reading the `.fits.gz` raw frames.

---

## 5. Speed-up opportunities (ranked, for the design doc)

1. **Optimize the B-spline kernel (`pypeit/core/bspline/`)** — *highest leverage,
   ~39 s / 32 %.* Cythonize / Numba-JIT / vectorize `intrv`, `bsplvn`,
   `cholesky_band`, `cholesky_solve`, `solution_arrays`. This simultaneously
   accelerates global sky, local sky, and extraction. **These functions run as
   pure Python today** — confirmed on both `release` (`pypeit/bspline/util.py`)
   and this branch (`pypeit/core/bspline/util.py`); neither imports a compiled
   path, so this is the status quo, *not* a refactor regression. A stale leftover
   `pypeit/bspline/_bspline.cpython-*.so` (no longer wired in, exposes nothing)
   shows these kernels were **historically Cython-compiled** — strong precedent
   that re-introducing a compiled build is feasible and likely the original
   design intent.
2. **Cut redundant B-spline work** — fewer reject iterations, reuse of basis/
   breakpoints across the 3 frames, caching of `action`/`intrv` results.
3. **Headless QA plotting** — force matplotlib `Agg`; consider making QA PNG
   generation optional, deferred, or parallel (≈12 s, 10 %).
4. **Image processing / CR rejection** — `lacosmic` + median filter ≈ 6 s; allow
   skipping or a faster CR algorithm where appropriate.
5. **Reduce NumPy allocation churn** in the fitting/rejection inner loops
   (preallocate, operate in-place, avoid masked-array overhead).
6. **Parallelism** — the pipeline is fully serial. Kast is single-detector so it
   gains little, but detector/mosaic- and frame-level parallelism is the obvious
   structural win for DEIMOS/Echelle/IFU and should be designed in.

Items 1–3 alone target ~55 s (~45 %) of the runtime on this simple setup; the
B-spline win generalizes to *every* spectrograph because sky/extract are common
to all paths.

---

## 6. Reproduce

```bash
# pypeit env python; from this scripts dir
cd PypeIt-development-suite/pypeitdev/speed_up/scripts
/home/xavier/miniconda3/envs/pypeit/bin/python profile_kast_blue.py   # cold run + cProfile (~2 min)
/home/xavier/miniconda3/envs/pypeit/bin/python analyze_profile.py     # pstats tables + timeline
```

Artifacts written to `pypeitdev/speed_up/Reports/`:
`shane_kast_blue_A.prof`, `.pstats.txt`, `.timeline.txt`, `.run.log`,
`.runmeta.json`, and this report.
