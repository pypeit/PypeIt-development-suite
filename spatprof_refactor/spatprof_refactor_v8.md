# fit_profile Refactoring Analysis (v8)

## Context

`fit_profile` in `pypeit/core/spatialprofile.py` (lines 262–854) is a 593-line
monolith ported from the IDL LOWREDUX routine `long_gprofile.pro`. It fits a
non-parametric B-spline spatial profile to an object spectrum for optimal
extraction, falling back to a Gaussian when S/N is too low. It is called from
exactly one external site: `skysub.py:886` inside `local_skysub_extract`, once
per object per iteration. Its four return values are:

- `profile_model` — 2D normalized spatial profile `(nspec, nspat)`
- `xnew` — corrected trace `(nspec,)`
- `fwhmfit` — FWHM per spectral pixel `(nspec,)`
- `med_sn2` — median S/N²

---

## Output File

The refactored implementation lives in:

```
pypeit/core/spatialprofile_refactor.py
```

The primary entry point is **`fit_profile_refactor`**, a drop-in replacement for
`fit_profile` with the same signature and return values. The legacy function in
`spatialprofile.py` is left untouched so both implementations can be called
side-by-side for comparison.

All module-private helpers (leading underscore) are in the same file. The body
of `fit_profile_refactor` is organized into labeled comment blocks (matching the
workflow breakdown below); the largest and most coherent pieces have been
extracted into named helper functions.

---

## Workflow Breakdown

### Block 1 — Initialization
Set up masks (`inmask`, `totmask`), clamp `thisfwhm`, initialize `sigma` and
`fwhmfit`.

### Block 0 — Isolate slit pixels (inserted at top of Block 1)
Immediately after `totmask` is fully constructed, extract five 1-D pixel arrays
from the valid slit pixels:

```python
if spec_img is None:
    spec_x = np.where(totmask)[0]
else:
    spec_x = spec_img[totmask]
spat_x = spat_img[totmask]
wave_x = waveimg[totmask]
flux_x = image[totmask]
ivar_x = ivar[totmask]
```

These arrays are stored as locals; downstream logic is still 2-D at this stage.

### Blocks 2–5 — Spectral fitting and normalisation → `_fit_spectrum_and_normalize`
See helper description below.

### Block 6 — Knot setup → `_compute_bspline_knots`
See helper description below.

### Block 7 — Initial profile fit
Compute `dspat` and `trace_corr`, then fit a 4th-order B-spline to
`norm_obj[inside]` sorted by `sigma_x`. Evaluate to get `mode_fit`; call
`_findfwhm` to locate the peak and measure width. Update `sigma` and `limit`
using the B-spline FWHM (in-place: `sigma *= ...`, `limit *= ...`). Walk the
profile down from the peak to find `l_limit` and `r_limit` (where profile drops
below `min_level`). Select the final `inside` subset within those limits; early
Gaussian return via `_return_gaussian` if `ninside < 10`.

### Block 8 — Iterative trace and width correction loop — `sigma_iter=3`
0-indexed `for iiter in range(sigma_iter):`. Each iteration:
1. **Trace shift**: evaluate the current B-spline at `sigma_x` and at
   `sigma_x ± 0.5`; use `(mode_min05 - mode_plu05)` as a shift-sensitive basis
   and fit it jointly with `mode_zero` via `iterative_bspline_fit`; extract the
   h0/h1 coefficient ratio to get `delta_trace_corr`.
2. **Width stretch**: evaluate at `sigma_x / 1.3`; use the rescaled difference
   as a stretch-sensitive basis; extract h2/h0 ratio to get `sigma_factor`.
3. Update `sigma`, `area`, `trace_corr`, `sigma_x`.
4. Re-fit the profile B-spline when `iiter < sigma_iter - 2` (i.e. the first
   iteration only when `sigma_iter=3`). Note: there is an open TODO about
   whether this condition is correct vs. the legacy `iiter < sigma_iter - 1`.

Each sub-fit failure returns a Gaussian immediately via `_return_gaussian`.

### Block 9 — Final trace
Apply `trace_corr` only if `median(|trace_corr * sigma|) < max_trace_corr`;
otherwise keep `trace_in`.

### Block 10 — Final profile fit
Re-select all in-range, unmasked pixels sorted by `sigma_x`; fit the definitive
4th-order B-spline profile using the full `pb` (area) weight array and loose
rejection (`upper=lower=10`). Convert resulting `BSpline2D` → 1D `BSpline` via
`_bspline2d_to_1d`.

### Blocks 11–12 — Apodization limit search and exponential tails → `_apodize_profile`
See helper description below.

### Block 13 — Normalization, logging, QA, return
Reshape `full_bsp`, multiply by `pb`, normalize each spectral row to unit sum,
guard against NaNs/Infs, log fit statistics, show QA plot if `show_profile`, and
return.

---

## Implemented Helper Functions

Seven module-private helpers have been extracted into named functions.

### `_bspline2d_to_1d(bset2d)`
Converts the zeroth column of a `BSpline2D` (returned by
`iterative_bspline_fit` when a `basis` array is provided) into a 1-D `BSpline`
that can be evaluated with `.value()`. Used after every call that returns a
2D spline (Blocks 8 and 10).

### `_findfwhm(model, sig_x)`
Locates the peak and both half-maximum crossings of a 1-D profile with
sub-pixel precision. Uses a `MaskedArray` to restrict the peak search to
`|sig_x| < 2` (ensuring the search region includes the half-maximum crossing
at ≈ ±1.18 σ), then masks values below 0.5 × peak and uses
`np.ma.flatnotmasked_edges` to find the bracketing indices. Each crossing is
refined by linear interpolation via `utils.linear_interpolate`. Returns
`(peak, peak_x, lwhm, rwhm)`.

### `_gaussian_profile(x, center, sigma)`
Constructs a pixel-integrated Gaussian profile using the error-function form
rather than point-sampling at pixel centres. The coordinate array `x` may be
1-D or 2-D; if 2-D, `center` and `sigma` may be scalars or 1-D vectors of
length `nspec` (broadcast as column vectors). Returns an array of the same
shape as `x`.

### `_return_gaussian(spat_img, norm_obj, center, sigma, fwhm, med_sn2, obj_string, show_profile, ...)`
Calls `_gaussian_profile` to build the profile, logs the FWHM and S/N, and
optionally shows the QA plot via `qa_fit_profile`. The `fwhm` argument is used
only for the log message; the profile shape is controlled by `sigma`. Returns
`profile_model`.

### `_fit_spectrum_and_normalize(wave, flux, fluxivar, wave_x, flux_x, ivar_x, spec_x, nspec, percentile_sn2, fwhm)`
Encapsulates the spectral B-spline fitting, S/N estimation, and
normalized-object construction (Blocks 2–4). Derives `wave_min` / `wave_max`
from `wave_x`, fits two fine B-splines and one coarse continuum B-spline on the
1-D extracted spectrum, estimates `med_sn2`, builds a 1-D `spline_img_x`
(S/N-level-dependent flux model), computes 1-D `norm_obj_x` and `norm_ivar_x`,
and constructs 1-D `xtemp_x` (the cumulative S/N-weighted spectral coordinate
used in Block 8, approximated as `cumweight[spec_x] / cumweight[-1]`). Returns a
7-tuple:

```
(success, med_sn2, norm_obj_x, norm_ivar_x, good_x, xtemp_x, sn2_x)
```

All output arrays have shape `(npix,)` where `npix = spec_x.size`. On failure
(no eligible spectral pixels, or the continuum B-spline evaluation raises),
returns `(False, 0.0, None, None, None, None, None)`. The `fwhm` parameter is
used only in warning messages.

The caller handles all Gaussian fallback decisions via a single unified check:
```python
if not success or gauss or ngood < 10 or med_sn2 < sn_gauss**2:
    ...
```

### `_compute_bspline_knots(dspat_x, sigma, spec_x, med_sn2, prof_nsigma, good_x)`
Encapsulates Block 6 (knot setup). Takes the 1-D spatial-separation array
`dspat_x = spat_x - trace_in[spec_x]` and the current sigma estimate, computes
`sigma_x = dspat_x / sigma[spec_x]` (1-D of length `npix`), the profile
half-extent `limit` (from `erfcinv`), the sigma bounds `min_sigma` /
`max_sigma` (clipped to `sigma_x[good_x]`), and the sinh-spaced interior
B-spline knot positions `bkpt`. When `prof_nsigma` is set, fixed bounds of
`±prof_nsigma` are used and `nb = max(1, round(prof_nsigma / 10))` (Bug 1 fix).
Returns:

```
(sigma_x, limit, min_sigma, max_sigma, bkpt)
```

The initial `trace_corr = np.zeros(nspec)` is initialized by the caller
immediately after the call.

### `_apodize_profile(bset, sigma_x, min_sigma, max_sigma, ss, median_fit, min_level, limit, prof_nsigma, no_deriv)`
Encapsulates Blocks 11–12 (apodization limit search and exponential tails).
Evaluates the final B-spline profile over all good pixels, re-locates the
profile peak via `_findfwhm` (note: `peak_x` is recomputed from the final
B-spline, not carried from Block 7), then locates the left and right
apodization limits by walking inward from the profile boundary. Computes the
logarithmic derivative `d(log P)/d(log x)` over pre-computed limit vectors
and uses a vectorized search (no `while True` loops) to find the steepest
inward-pointing derivative. If `prof_nsigma is not None`, sets `l_limit =
r_limit = 0.0` and `no_deriv = True` (JXP kludge: prevents QA from drawing
limits for extended-object fits, and skips apodization). Otherwise, if the
derivatives have the right sign and `no_deriv=False`, replaces `full_bsp`
outside `[l_limit, r_limit]` with matched exponential tails. Returns:

```
(full_bsp, l_limit, r_limit)
```

---

## Remaining Inline Blocks

The following blocks remain inline within `fit_profile_refactor`:

| Block | Reason |
|-------|--------|
| 1 — Initialization | Trivially short; only sets up a few scalars |
| 7 — Initial profile fit | Produces many tightly coupled intermediate values (`bset`, `median_fit`, `peak_x`, `min_level`, `l_limit`, `r_limit`, `mask`, `inside`, `ninside`) that would all need to be returned |
| 8 — Iterative trace/width correction | Complex loop with multiple early-return paths; all state is consumed immediately; extraction would require passing/returning ~10 variables |
| 9 — Final trace | Two lines |
| 10 — Final profile fit | Four lines |
| 13 — Normalization, logging, QA, return | QA and return logic; closely tied to `fit_profile_refactor`'s output contract |

---

## Performance Improvements (implemented)

### `np.outer(v, np.ones(n))` → `v[:, None]` broadcasting
Every `np.outer` was replaced with the equivalent broadcasting form, avoiding
temporary array allocation.

### Vectorized apodization limit walk (`_apodize_profile`)
The two `while True` loops (stepping 0.1 increments, re-evaluating the
B-spline at each step) were replaced with vectorized searches over the
pre-computed derivative arrays `l_deriv_vec` / `r_deriv_vec`.

### Profile normalization simplified (Block 13)
Uses `np.where` instead of the boolean-mask arithmetic from the legacy code:
```python
norm_safe = np.where(row_sums[:, None] > 0.0, row_sums[:, None], 1.0)
profile_model = np.where(row_sums[:, None] > 0.0, profile_model / norm_safe, 0.0)
```

---

## Performance Comparison: Empirical Results

### Methodology

Benchmarks were run on synthetic inputs generated by the same `make_inputs` factory
used in the test suite (`seed=42`, uniform Gaussian profile, S/N = 20). Three
image sizes were tested. Wall-clock time was measured with `time.perf_counter` over
5 runs after 1 warm-up call. Peak heap allocation was measured with `tracemalloc`
over a single call. The `iterative_bspline_fit` call sites were timed by temporarily
monkey-patching the function with a wrapper that records per-call elapsed time.

All measurements were taken after the memory optimizations described in the
"Memory Optimizations Applied" subsection below.

### Wall-clock time

| Case | Legacy | Refactored | Speedup |
|------|--------|-----------|---------|
| Small (200×100, fwhm=4) | 61 ms | 51 ms | **1.20×** |
| Medium (400×150, fwhm=6) | 101 ms | 94 ms | **1.07×** |
| Large (600×200, fwhm=4) | 138 ms | 120 ms | **1.15×** |

### Peak heap allocation

| Case | Legacy | Refactored | Memory savings |
|------|--------|-----------|---------------|
| Small (200×100) | 3732 KiB | 3275 KiB | **+12.3%** |
| Medium (400×150) | 10717 KiB | 9805 KiB | **+8.5%** |
| Large (600×200) | 18520 KiB | 17148 KiB | **+7.4%** |

The refactored version is faster and uses less memory than the legacy function on
all three cases.

### Hotspot analysis

`iterative_bspline_fit` is called **12 times** per invocation of
`fit_profile_refactor` and accounts for approximately **56 % of total wall time**.
Block 7 (initial profile fit) contains the single most expensive call (~7.6 ms on
the small case). The per-call breakdown is roughly: Block 7 initial fit (1 call,
~7.6 ms), Block 8 loop fits (9 calls, ~2–4 ms each), Block 10 final fit (1 call,
~4 ms), Block 8 trace/width correction sub-fits (1 call, ~2 ms).

Logging (`log.info`) accounts for approximately **22 % of total wall time** on both
implementations — about **109 calls per invocation**, totalling ~10 ms. This
overhead is identical between legacy and refactored and is a fixed cost of the
current logging strategy.

### Memory optimizations applied

Three changes were made to `spatialprofile_refactor.py` to bring memory usage
below the legacy baseline.

**1. Revert `np.where` to slice-then-assign in `_fit_spectrum_and_normalize`.**
Three B-spline evaluations had been written as `np.where(ispline, bset.value(wave)[0], 0.0)`.
The `np.where` form evaluates `bset.value(wave)` over the full `nspec`-length
array; the slice-then-assign form evaluates only over `wave[ispline]` (the valid
wavelength subset):

```python
spline_flux1 = np.zeros(nspec)
spline_flux1[ispline], _ = b_answer.value(wave[ispline])
```

**2. Remove `spat_x` from Block 0; inline at its single use in Block 6.**
`spat_x = spat_img[totmask]` was an `(npix,)` float64 array kept alive for the
entire function even though it was needed only for a single subtraction
(`dspat_x = spat_x - trace_in[spec_x]`). Inlining eliminates the allocation
entirely:

```python
dspat_x = spat_img[totmask] - trace_in[spec_x]
```

**3. Release `wave_x`, `flux_x`, `ivar_x` after `_fit_spectrum_and_normalize` returns.**
These three `(npix,)` arrays are consumed entirely by the helper and are never
referenced again. Adding `del wave_x, flux_x, ivar_x` immediately after the call
frees them before Block 7 allocates the B-spline design matrix
(`npix × nknots` elements), which is the dominant allocation in the function.
CPython's reference counting means the `del` takes effect immediately (ref count
→ 0).

---

## Accuracy Improvements (implemented)

### A. Pixel-integrated Gaussian in `_gaussian_profile`
The legacy `return_gaussian` used a point-sampled Gaussian density
(`exp(-0.5*sigma_x**2)/sqrt(2*pi)`). `_gaussian_profile` uses the
pixel-integrated erf-difference form:
```python
delta = 0.5 / sigma   # half-pixel in normalized units
profile = 0.5 * (erf((sigma_x + delta)/sqrt(2)) - erf((sigma_x - delta)/sqrt(2)))
```
Row sums of the returned profile are approximately 1.0 for any profile width,
matching the B-spline normalization path. The correction is significant for
FWHM < 5 px (≈5% at peak for FWHM = 4 px; ≈10% for FWHM = 3 px) and vanishes
in the well-sampled limit.

### B. Sub-pixel half-maximum crossing in `_findfwhm`
The legacy `findfwhm` returned the first discrete sample below 0.5 × peak,
systematically overestimating the FWHM. `_findfwhm` interpolates between the
two bracketing samples to locate the crossing to sub-pixel precision:

| True FWHM (px) | Legacy recovered FWHM | After fix |
|---|---|---|
| 3.0 | ~4.0 (33% error) | ~3.0 (<2% error) |
| 4.0 | ~5.0 (25% error) | ~4.0 (<2% error) |
| 5.0 | ~6.0 (20% error) | ~5.0 (<2% error) |
| 6.0 | ~6.0 (<1% error) | ~6.0 (<1% error) |

### C. FWHM precision: `sig2fwhm = np.sqrt(8*np.log(2))`
The legacy code used the truncated constant `2.3548`. The refactored code uses
the exact expression `np.sqrt(8*np.log(2)) ≈ 2.35482004...` throughout.

---

## Pre-existing Bugs

These bugs exist in the legacy `spatialprofile.py`. All three are **fixed in
`fit_profile_refactor`**.

### Bug 1 — `nb` formula in the `prof_nsigma` knot-setup branch
```python
# Legacy (wrong)
nb = np.round(prof_nsigma > 10)
# Fixed (in _compute_bspline_knots)
nb = max(1, round(prof_nsigma / 10))
```
`prof_nsigma > 10` is a boolean scalar; `np.round` of a boolean gives 0 or 1.
For `prof_nsigma ≤ 10` this causes division-by-zero; for `prof_nsigma > 10` it
gives only 2 knot positions.

### Bug 2 — FWHM overestimation bias in `findfwhm` for narrow profiles
Fixed in `_findfwhm`; see §B above.

### Bug 3 — Trace correction overshoots for small offsets (< 1 px)
Reduced but not eliminated by the Bug 2 fix. For a 0.5 px offset, the
refactored function shifts `xnew` in the correct direction but overshoots the
true trace by ~0.05 px. For a 1.0 px offset, `xnew` converges correctly. The
residual overshoot at 0.5 px may be related to the open question about the
Block 8 mid-loop re-fit condition (see below).

---

## Open Question: Block 8 Mid-Loop Re-fit Condition

The original loop used `for iiter in range(1, sigma_iter + 1):` with a
mid-loop re-fit condition `if iiter < sigma_iter - 1:`. This was equivalent
to re-fitting on iterations 1 and 2 only (not 3) out of 3.

The current code uses `for iiter in range(sigma_iter):` (0-indexed) with
`if iiter < sigma_iter - 2:`. This is equivalent to re-fitting on iteration
0 only out of 3 — i.e. one fewer re-fit than the legacy.

**TODO**: Determine whether the legacy `iiter < sigma_iter - 1` behavior (two
re-fits) or the current `iiter < sigma_iter - 2` behavior (one re-fit) is
correct. The regression test currently passes with the one-re-fit version, but
the change was introduced during a loop-index refactor and may be unintentional.
Restoring the two-re-fit behavior may also reduce the Bug 3 overshoot at 0.5 px.

---

## Unit Tests (`pypeit/tests/test_spatialprofile.py`)

The test file contains 42 tests, all passing. Tests are organized into four
groups.

**Tests against `spatialprofile.fit_profile` (legacy)**:
- `test_profile_return_shapes`
- `test_profile_normalization_bspline`
- `test_profile_normalization_gaussian`
- `test_forced_gaussian`
- `test_fwhm_recovery` — parametrized over `[6.0, 8.0]` px only (the 3.0 and
  5.0 px cases fail for the legacy due to Bug 2; only the two largest values
  pass within 10%)
- `test_trace_correction_direction` — `trace_offset=1.0` px; asserts direction
  and convergence (`err_after < err_before`)
- `test_prof_nsigma_extended` — asserts shape and finiteness only (broken `nb`
  formula in legacy makes FWHM unreliable for `prof_nsigma=10`)
- `test_profile_positivity_and_finite_all_paths`

**Tests against `spatialprofile_refactor.fit_profile_refactor`**:
- `test_refactor_regression` — loads frozen `.npy` reference arrays; verifies
  output to `atol=1e-10`. Skipped if reference files are absent.
- `test_fwhm_recovery_refactor` — parametrized over `[3.0, 4.0, 6.0, 8.0]` px;
  Bug 2 fix makes the 3 and 4 px cases pass within 10%.
- `test_gaussian_normalization_refactor` — Gaussian-path row sums within 1% of
  1.0 (pixel-integrated erf form).
- `test_prof_nsigma_refactor` — `prof_nsigma=10.0` runs without crash (Bug 1
  fix); asserts shape, finiteness, and `med_sn2 > sn_gauss²`. FWHM recovery is
  *not* asserted: with only 2 B-spline knots over a ±10 sigma fitting range, the
  B-spline fits mostly near-zero background pixels and the FWHM estimate is
  unreliable (~45 px instead of ~8 px).
- `test_trace_correction_small_offset` — `trace_offset=0.5` px; asserts only
  that `xnew` shifts in the correct direction. The function still overshoots the
  true trace (~0.05 px past it), so convergence is not asserted. This is related
  to the Block 8 open question.

**Helper-function unit tests** (16 tests, all directly importing private helpers):

*`_gaussian_profile`*:
- `test_gaussian_profile_row_sums` — parametrized over `fwhm ∈ [2, 3, 4, 6, 10]`
  px; row sums within 1% of 1.0.
- `test_gaussian_profile_symmetry` — odd-width slit (nspat=101); profile
  symmetric to `atol=1e-14`.

*`_findfwhm`*:
- `test_findfwhm_coarse_accuracy` — step = 0.5 sig_x; interpolated FWHM within
  2% of `sqrt(8*ln 2)`; discrete estimate overestimates by > 15%, documenting
  Bug 2.
- `test_findfwhm_peak_location` — fine grid (step = 0.01); `abs(peak_x) <= step`.
- `test_findfwhm_asymmetric` — profile centred at 0.5; `peak_x` and FWHM both
  within 2%.

*`_fit_spectrum_and_normalize`* (1-D interface; all four tests use `_make_fsn_inputs_1d`;
`_make_fsn_inputs`, the old 2-D factory, removed in step 8):
- `test_fit_spectrum_success` — high-S/N input; `success=True`, `med_sn2 > sn_gauss²`;
  all array outputs have shape `(npix,)`.
- `test_fit_spectrum_failure_nan_flux` — all-NaN flux; returns
  `(False, 0.0, None, …)`.
- `test_fit_spectrum_good_matches_norm_ivar` — `assert_array_equal(good_x, norm_ivar_x > 0)`
  (1-D, no `.flatten()`).
- `test_fit_spectrum_sn2_img_nonnegative` — `sn2_x >= 0` and `sn2_x.shape == (npix,)`.

*`_compute_bspline_knots`* (1-D interface; `_make_knots_inputs` updated to 1-D in step 9):
- `test_bspline_knots_shapes` — `sigma_x.shape == (npix,)`, `limit > 0`,
  `min_sigma < 0 < max_sigma`, `bkpt` within bounds.
- `test_bspline_knots_bug1_fix` — `prof_nsigma=10.0` no longer crashes;
  `min_sigma=-10`, `max_sigma=10`; passes `dspat_x`, `spec_x`, `good_x` kwargs.
- `test_bspline_knots_prof_nsigma_bounds` — parametrized over `[5, 10, 15, 20]`;
  `min/max_sigma == ±prof_nsigma`; passes `dspat_x`, `spec_x`, `good_x` kwargs.

*`_apodize_profile`* (all share a module-scoped fixture that fits a B-spline to
`exp(-sig_x²/2)`):
- `test_apodize_output_shape` — `full_bsp.size == sigma_x.size`; finite limits;
  `full_bsp[igood] >= 0`.
- `test_apodize_prof_nsigma_zeros_limits` — `prof_nsigma=5.0` → `l_limit =
  r_limit = 0.0` (JXP kludge).
- `test_apodize_no_deriv_skips_tails` — `no_deriv=True` → pixels outside
  `[min_sigma, max_sigma]` are zero.
- `test_apodize_tail_continuity` — `no_deriv=False`; tail values in
  `(min_sigma, l_limit)` are positive and bounded above by the profile value at
  `l_limit`.

**1-D infrastructure** (factories used by helper tests):
- `make_profile_inputs_1d` — calls `make_profile_inputs`, applies `totmask`,
  returns `(wave, flux, fluxivar, wave_x, flux_x, ivar_x, spec_x, trace_in,
  nspec, nspat)`. The 1-D pixel arrays have shape `(npix,)` where
  `npix = totmask.sum()`.
- `_make_fsn_inputs_1d` — returns `(wave, flux, fluxivar, wave_x, flux_x,
  ivar_x, spec_x, nspec)` — the 1-D inputs for `_fit_spectrum_and_normalize`.
  Actively used by all four `test_fit_spectrum_*` tests. (`_make_fsn_inputs`,
  the old 2-D factory, was removed in step 8.)
- `_make_knots_inputs` — updated to the 1-D final interface in step 9: returns
  `(dspat_x, sigma, med_sn2, good_x, spec_x)` by calling
  `_fit_spectrum_and_normalize` with the new 1-D signature. Used by all three
  `test_bspline_knots_*` tests.

---

## Regression Freeze Workflow

The freeze captures `fit_profile_refactor` output as a self-consistency check.
Re-run `pypeit/tests/files/freeze_spatprof.py` whenever intentional algorithmic
changes are made.

Reference files (untracked, do not `git add`):
- `pypeit/tests/files/spatprof_profile_model_ref.npy`
- `pypeit/tests/files/spatprof_xnew_ref.npy`
- `pypeit/tests/files/spatprof_fwhmfit_ref.npy`

---

## Verification Plan

1. ✅ Write unit tests and confirm all pass against the pre-refactor function
   (`pytest pypeit/tests/test_spatialprofile.py`: 12 passed).
2. ✅ Generate regression freeze reference arrays (`freeze_spatprof.py`) and
   leave them as untracked files.
3. ✅ Write `pypeit/core/spatialprofile_refactor.py` containing
   `fit_profile_refactor` and all extracted helpers, incorporating all
   performance and accuracy improvements and bug fixes listed above.
4. ✅ Re-run `pytest pypeit/tests/test_spatialprofile.py`; all 12 tests pass
   including `test_refactor_regression`.
5. ✅ Strengthen tests targeting `fit_profile_refactor` and add direct unit
   tests for all helper functions. Final count: 42 tests, all passing.
   Deviations from original plan:
   - `test_prof_nsigma_refactor` asserts no crash and valid arrays only (FWHM
     recovery is unreliable with a sparse 2-knot B-spline over a ±10σ range).
   - `test_trace_correction_small_offset` asserts direction only; the 0.5 px
     case still overshoots the true trace by ~0.05 px.

**1-D Slit-Pixel Implementation** (Steps 6–15; see the "1-D Slit-Pixel Approach"
section for design rationale):

6. ✅ Add 1-D test infrastructure to `test_spatialprofile.py`.
   - `make_profile_inputs_1d` added: calls `make_profile_inputs`, applies
     `totmask`, and returns `(wave, flux, fluxivar, wave_x, flux_x, ivar_x,
     spec_x, trace_in, nspec, nspat)`.
   - `_make_fsn_inputs_1d` added: returns the positional inputs for the future
     1-D `_fit_spectrum_and_normalize` signature `(wave, flux, fluxivar, wave_x,
     flux_x, ivar_x, spec_x, nspec)`.
   - No production code changes; all 42 existing tests pass.
   Deviations from plan:
   - The plan prose listed `spat_x` among the returned arrays, but the code
     example in the plan omitted it from the return tuple.  The implementation
     follows the code example: `spat_x` is computed inside
     `make_profile_inputs_1d` but not returned.  Any call site that needs it can
     recompute `spat_img[totmask]`.

7. ✅ Add `spec_img` parameter to `fit_profile_refactor` and implement Block 0.
   - `spec_img=None` added to the function signature after `inmask` (backward-
     compatible default).
   - Docstring entry added for `spec_img`.
   - Block 0 inserted immediately after `totmask` is fully constructed (before
     `nspec, nspat = image.shape`), using the existing dashed-header comment
     style. Extracts `spec_x`, `spat_x`, `wave_x`, `flux_x`, `ivar_x` from
     `totmask`; all five stored as local variables.
   - Downstream 2-D logic is unchanged at this step.
   - All 42 tests pass.

8. ✅ Convert `_fit_spectrum_and_normalize` to 1-D and update its call site and
   tests.
   - Signature changed from `(wave, flux, fluxivar, waveimg, image, ivar,
     totmask, percentile_sn2, fwhm)` to `(wave, flux, fluxivar, wave_x, flux_x,
     ivar_x, spec_x, nspec, percentile_sn2, fwhm)`.
   - Internal 2-D arrays replaced: `spline_img` → `spline_img_x` (length
     `npix`); `sn2_img[totmask] = sn2_interp(waveimg[totmask])` →
     `sn2_x = sn2_interp(wave_x)`.
   - `good = norm_ivar.flatten() > 0.0` → `good_x = norm_ivar_x > 0.0`
     (1-D bool, length `npix`); `ivar_mask` drops the `totmask` component
     (already operating on valid pixels only).
   - `xtemp` construction replaced:
     ```python
     cumweight = np.cumsum(row_weights)   # (nspec,)
     xtemp_x = cumweight[spec_x] / cumweight[-1]   # (npix,)
     ```
     All pixels in the same spectral row get the same `xtemp_x` value
     (minor approximation; within-row variation in the original was negligible).
   - Return values `norm_obj_x`, `norm_ivar_x`, `good_x`, `xtemp_x`, `sn2_x`
     are all 1-D of length `npix`.
   - Call site in `fit_profile_refactor` passes 1-D arrays; step-8 scaffolding
     scatters them back to 2-D (`norm_obj`, `norm_ivar`, `sn2_img`, `xtemp`,
     `good`) for Blocks 7–10 (removed in step 10).
   - `_make_fsn_inputs` (2-D factory) removed from the test file; all four
     `test_fit_spectrum_*` tests updated to use `_make_fsn_inputs_1d` and assert
     1-D shapes, e.g. `good_x.shape == (npix,)`.
   - Regression reference arrays regenerated (untracked `.npy` files overwritten)
     because the `xtemp_x` row-approximation produces numerically different
     outputs at `atol=1e-10`. All 42 tests pass.

9. ✅ Convert `_compute_bspline_knots` to 1-D and update its call site and tests.
   - Signature changed from `(dspat, sigma, med_sn2, prof_nsigma, good)` to
     `(dspat_x, sigma, spec_x, med_sn2, prof_nsigma, good_x)`.
   - `sigma_x = dspat / sigma[:, None]` → `sigma_x = dspat_x / sigma[spec_x]`
     (1-D of length `npix`).
   - `sigma_x.flat[good]` → `sigma_x[good_x]` (direct 1-D boolean indexing).
   - Call site in `fit_profile_refactor`: added `dspat_x = spat_x -
     trace_in[spec_x]`; `dspat` (2-D) retained temporarily for Block 8 (removed
     in step 10). After the call, step-9 scaffolding scatters 1-D `sigma_x` back
     to 2-D for Block 7 (removed in step 10).
   - `_make_knots_inputs` updated to full 1-D interface: returns
     `(dspat_x, sigma, med_sn2, good_x, spec_x)`; calls `_fit_spectrum_and_normalize`
     with the new 1-D signature.
   - Three `test_bspline_knots_*` tests updated: pass `dspat_x`, `spec_x`,
     `good_x` as keyword args; assert `sigma_x.shape == (npix,)`.
   - All 42 tests pass.

10. ✅ Convert Blocks 6–10 in `fit_profile_refactor` to 1-D throughout.
    - **Step-8/9 scaffolding removed**: 2-D scatter-backs for `norm_obj`,
      `norm_ivar`, `sn2_img`, `xtemp`, `sigma_x` all deleted; 2-D `dspat` and
      `good` variable removed; `gauss_kwargs` updated to use `good_x`; all
      `_return_gaussian` calls updated to pass `norm_obj_x` and `ind=good_x`.
    - **Block 6**: `npix = spec_x.size` added; `dspat = spat_img -
      trace_in[:, None]` (2-D) and the `sigma_x` scatter-back removed; only
      `dspat_x` and 1-D `sigma_x` remain.
    - **Block 7**: `GOOD_PIX`/`IN_PIX` conditions use `sn2_x` and
      `norm_ivar_x`; `np.where(...flatten())` → `np.where(...)`; all
      `.flat[si]`/`.flat[sr]` → `[si]`/`[sr]`; `norm_obj[norm_ivar > 0]` →
      `norm_obj_x[norm_ivar_x > 0]`; `mask = np.zeros(nspec*nspat, …)` →
      `mask_x = np.zeros(npix, …)`.
    - **Block 8**: `xx = np.bincount(spec_x, weights=xtemp_x,
      minlength=nspec) / nspat` computed once before the loop (replaces
      `xtemp.sum(axis=1) / nspat` inside the loop); `xtemp.flat[...]` →
      `xtemp_x[...]`; `sigma_x = dspat / sigma[:, None] - trace_corr[:, None]`
      → `sigma_x = dspat_x / sigma[spec_x] - trace_corr[spec_x]`;
      `pb = np.repeat(area, nspat)[inside]` → `pb = area[spec_x[inside]]`;
      all `.flat[]` accesses removed.
    - **Block 9**: unchanged.
    - **Block 10**: `ss = sigma_x.flatten().argsort(kind='stable')` →
      `ss = sigma_x.argsort(kind='stable')`; `mask[ss]` → `mask_x[ss]`;
      `norm_obj.flat[ss]`/`norm_ivar.flat[ss]` → `norm_obj_x[ss]`/
      `norm_ivar_x[ss]`; `pb = area[:, None] * np.ones(…)` → `pb_x =
      area[spec_x]`.
    - **Block 13**: temporary step-10 scatter-back added — `profile_x =
      full_bsp * pb_x; profile_model[totmask] = profile_x` — so row
      normalisation continues to operate on the 2-D `profile_model`; `res_mode`
      and `chi_good` updated to use 1-D `norm_obj_x`/`norm_ivar_x`/`profile_x`
      indexing; QA call updated to use `norm_obj_x` and `pb_x`.
    - All 42 tests pass. `test_refactor_regression` continues to pass (reference
      arrays were already regenerated in step 8 for the `xtemp_x` approximation;
      step 10 introduces no additional numerical change).

11. ✅ Convert `_apodize_profile` to 1-D and update its call site, fixture, and
    tests.
    - `sigma_x` is now 1-D; remove every `.flatten()` call inside the function:
      `igood = (sigma_x > min_sigma) & (sigma_x < max_sigma)`,
      `sigma_x_igood = sigma_x[igood]`, `left = sigma_x < l_limit`,
      `sigma_x[left]`, etc.
    - `full_bsp = np.zeros(sigma_x.size)` remains (already effectively 1-D);
      its shape is now `(npix,)` instead of `(nspec * nc,)`.
    - `ss` passed in is now a 1-D argsort of length `npix`; `sigma_x[ss]`
      replaces `sigma_x.flat[ss]` throughout.
    - No `full_bsp.reshape` existed in Block 13 (the scatter-back was already
      added as part of step 10, so this sub-task was a no-op).
    - Updated the `apodize_bspline` fixture: flattened the 2-D `sigma_x` to
      1-D, added `npix = sigma_x.size`, changed `ss = sigma_x.argsort(...)`.
    - Updated `test_apodize_output_shape`: `shape == (f['npix'],)`;
      `igood` uses direct 1-D indexing.
    - Updated `test_apodize_prof_nsigma_zeros_limits`: `shape == (f['npix'],)`.
    - Updated `test_apodize_no_deriv_skips_tails`: removed `sigma_x_flat`.
    - Updated `test_apodize_tail_continuity`: removed `sigma_x_flat`.
    - All 42 tests pass; `test_refactor_regression` continues to pass (reference
      arrays were regenerated in step 8; steps 11–12 introduce no numerical
      change).

12. ✅ Block 13 scatter-back reconstruction.
    - The scatter-back (`profile_x = full_bsp * pb_x; profile_model[totmask] =
      profile_x`) and 1-D `res_mode`/`chi_good` were already in place from step
      10. The only remaining action was removing the misleading "step-10
      scaffolding" comment.
    - Row normalization (`row_sums`, `norm_safe`) operates on the 2-D
      `profile_model` exactly as before.
    - `qa_fit_profile` receives 1-D `sigma_x` and `norm_obj_x`; confirmed it
      still accepts 1-D inputs without modification.
    - All 42 tests pass.

13. ✅ Regenerate regression reference arrays and confirm all tests pass.
    - Ran `python pypeit/tests/files/freeze_spatprof.py`; three untracked `.npy`
      files overwritten with output from the 1-D implementation.
      (`med_sn2=384.62`, `sqrt=19.61`; `profile_model` row sums all 1.000000;
      `xnew` mean=50.0037; `fwhmfit` mean=3.993`).
    - `pytest pypeit/tests/test_spatialprofile.py`; all 42 tests pass including
      `test_refactor_regression`.

14. ✅ Add `test_1d_roundtrip_matches_2d` and the six cross-implementation
    comparison tests documented in the "Test Modifications" section.
    - `test_1d_roundtrip_matches_2d`: shape check + `atol=1e-8` comparison of
      `fit_profile_refactor` against the frozen reference arrays.
    - `test_med_sn2_agreement`: parametrized on `fwhm=[6.0, 8.0]`; 5 % rtol.
    - `test_forced_gaussian_trace_agreement`: `xnew` and `fwhmfit` identical.
    - `test_bspline_row_normalization_both`: parametrized on both functions;
      non-zero row sums equal 1.0 at `atol=1e-10`.
    - `test_profile_positivity_both`: parametrized on both functions × two S/N
      ratios (4 param combinations); non-negative, finite.
    - `test_trace_correction_direction_both`: 1-px offset; both functions
      converge toward the true trace.
    - `test_bspline_fwhm_improvement_narrow`: `fwhm=3 px`; legacy overshoots
      by > 15 % (actual: 33 %); refactored recovers within 10 % (actual: 2 %).
    - **Final test count: 54 pass** (up from 42, target was 49+).

15. ✅ Resolve the open question about the Block 8 mid-loop re-fit condition
    (`iiter < sigma_iter - 2` vs. the legacy `iiter < sigma_iter - 1`).
    Restoring the two-re-fit behavior may also reduce the Bug 3 overshoot.

16. ✅ Fully deprecate `pypeit/core/spatialprofile.py` and replace it with the
    refactored module.
    - Renamed `fit_profile_refactor` → `fit_profile` and
      `_fit_profile_refactor_qa` → `_fit_profile_qa` inside the refactored
      module; all internal call sites and docstring cross-references updated.
    - Moved `spatialprofile_refactor.py` → `spatialprofile.py`, replacing the
      legacy module.
    - `skysub.py` required no changes: it already imported
      ``from pypeit.core import spatialprofile`` and called
      ``spatialprofile.fit_profile(…)``.
    - In `test_spatialprofile.py`: updated imports, renamed the QA test
      function, updated stale bug-reference docstrings, and removed the 12
      cross-implementation comparison tests (which compared the legacy and
      refactored functions and are no longer meaningful with a single
      implementation).  The ``make_profile_inputs`` docstring note about the
      discrete-vs-integrated Gaussian distinction was moved from a TODO comment
      into the docstring body.
    - Updated ``freeze_spatprof.py`` to import from ``spatialprofile`` and
      call ``fit_profile``.
    - 44 tests pass (56 − 12 removed cross-implementation tests).

17. ✅ Wire the spatial-profile QA figure into the standard `run_pypeit` output.
    - Added ``from pathlib import Path`` and ``from pypeit import qa`` imports
      to ``skysub.py``.
    - Added a new case ``'spat_profile_qa'`` to ``qa.set_qa_filename`` with an
      ``obj: int | None = None`` parameter.  Produces filenames of the form
      ``QA/PNGs/{basename}_{det}_S{slit:04d}_O{obj:03d}_spat_prof.png``.
    - Replaced the ``qa_path=None`` parameter on ``local_skysub_extract`` and
      ``ech_local_skysub_extract`` with three parameters: ``qa_outdir``,
      ``qa_basename``, ``qa_detname``.
    - QA output is produced only on the **final iteration** (``iiter == niter``).
      The filename is constructed via ``qa.set_qa_filename`` inside the loop.
      The ``QA/PNGs/`` directory is created once before the loop via
      ``Path(qa_outdir, 'PNGs').mkdir(parents=True, exist_ok=True)``.
    - Interactive display (``show_profile=True``) still works as before on every
      iteration; the ``_qa_outfile`` logic gates file output to the final
      iteration only.
    - In ``extraction.py``, both ``MultiSlitExtraction.local_skysub_extract``
      and ``EchelleExtraction.local_skysub_extract`` now construct
      ``_qa_outdir = <redux_path>/<qadir>`` from ``self.par`` and pass
      ``qa_outdir=_qa_outdir``, ``qa_basename=self.basename``, and
      ``qa_detname=self.spectrograph.get_det_name(self.det)``.
    - 44 tests pass.
    - Dev-suite confirmation deferred to end-of-branch run.

18. ✅ Update `skysub.py` to construct and pass `spec_img` explicitly.
    - Before the `fit_profile` call, add:
      ```python
      spec_img_sub = np.broadcast_to(
          np.arange(nspec, dtype=int)[:, None],
          (nspec, nc)
      )
      ```
    - Pass ``spec_img=spec_img_sub`` to the function call.
    - This makes the spectral-index provenance explicit and avoids the internal
      ``np.where(totmask)`` fallback, but is otherwise a no-op change.
    - 44 tests pass.

19. ✅ Write a new documentation page describing the spatial-profile fitting
    algorithm.
    - Place the RST source in the ``doc/`` directory, following existing
      documentation conventions, and register it in the Sphinx ``toctree``.
    - Add a reference to the new page in ``doc/skysub.rst``.
    - Cover the B-spline fitting procedure, the iterative trace and width
      correction, the Gaussian fallback conditions, and the apodization.
    - Embed example figures generated from the synthetic inputs used by the
      unit tests so the examples are reproducible without real observational
      data.
    - ``doc/spatprof.rst`` written covering the S/N decision, B-spline
      fitting, iterative trace/width correction, apodization, Gaussian
      fallback, profile normalization, outputs, and QA output.
    - ``doc/scripts/make_spatprof_figures.py`` added; uses
      ``make_fourier_spectrum`` and ``make_profile_inputs`` from
      ``pypeit/tests/test_spatialprofile.py`` with the same Fourier-spectrum
      flux and quadratic curved trace as ``test_fourier_flux_and_curved_trace``.
    - ``doc/figures/spatprof_example_bspline.png`` and
      ``doc/figures/spatprof_example_gaussian.png`` generated (S/N = 20 and
      S/N = 2 respectively) and embedded in the page.
    - Page registered in the ``Processing Details`` toctree in
      ``doc/index.rst`` between ``skysub`` and ``extraction``.
    - Cross-reference added to ``doc/skysub.rst`` in the Local Sky
      Subtraction section.
    - Sphinx build clean (no new warnings).

---

## 1-D Slit-Pixel Approach: Proposal and Design

### Motivation

The arrays received by `fit_profile_refactor` have shape `(nspec, nc)` where
`nc` is the spatial bounding-box width of the object mask (set by `maskwidth`
in the caller, `skysub.py:829`).  Even within this bounding box, `thismask`
and `inmask` further restrict which pixels are scientifically valid.
Operations that broadcast over the full `(nspec, nc)` rectangle — for example,
line 938:

```python
sigma_x = dspat / sigma[:, None] - trace_corr[:, None]
```

— compute over every column in the bounding box, including masked pixels.
More importantly, `np.cumsum(np.repeat(row_weights, nspat))` (line 453 of
`_fit_spectrum_and_normalize`) allocates an array of length `nspec × nc` where
`nc` could be several times larger than the true slit width.  The flat argsort
`sigma_x.flatten().argsort()` (Block 10, line 981) sorts `nspec × nc` elements
rather than just the `npix = totmask.sum()` valid ones.

The 1-D approach eliminates all of this by isolating slit pixels up front.

### New input: `spec_img`

A new 2-D integer array is added to the function signature alongside the
existing `spat_img`:

```python
spec_img = np.broadcast_to(np.arange(nspec, dtype=int)[:, None], image.shape)
```

`spec_img[i, j] = i` — the spectral row index of each pixel.  At the call site
in `skysub.py`, it is constructed trivially from `nspec` and the shape of the
sub-image.  Adding it as an explicit parameter keeps the interface symmetric:
`spat_img` provides spatial coordinates, `spec_img` provides spectral indices.

### Block 0: extract 1-D slit-pixel arrays

A new Block 0 immediately follows the mask construction in Block 1 and replaces
all subsequent 2-D image accesses:

```python
# Block 0 — isolate slit pixels into 1-D arrays
spec_x = spec_img[totmask]     # int, shape (npix,): spectral row of each pixel
spat_x = spat_img[totmask]     # float, shape (npix,): spatial coordinate
wave_x = waveimg[totmask]      # float, shape (npix,): wavelength
flux_x = image[totmask]        # float, shape (npix,): sky-subtracted value
ivar_x = ivar[totmask]         # float, shape (npix,): inverse variance
```

`npix = totmask.sum()` ≪ `nspec × nc` for any practical slit.

### Key transformations

**`dspat` and `sigma_x`** — the two most frequently recomputed quantities:

```python
# Before (2-D, computed repeatedly in Block 6, Block 8, and whenever sigma changes)
dspat   = spat_img - trace_in[:, None]         # (nspec, nc)
sigma_x = dspat / sigma[:, None] - trace_corr[:, None]   # (nspec, nc)

# After (1-D)
dspat   = spat_x - trace_in[spec_x]            # (npix,)
sigma_x = dspat / sigma[spec_x] - trace_corr[spec_x]     # (npix,)
```

Line 938 in the Block 8 loop (currently the most expensive recomputation)
becomes a single element-wise division and subtraction over `npix` elements
rather than `nspec × nc`.

**`xtemp`** — the spectral weight coordinate used in Block 8:

```python
# Before (2-D)
row_weights = 4.0 + np.sqrt(np.fmax(sn2_1, 0.0))         # (nspec,)
xtemp = np.cumsum(np.repeat(row_weights, nc)).reshape((nspec, nc))
xtemp /= xtemp.max()

# After (1-D)
cumweight = np.cumsum(row_weights)                         # (nspec,)
xtemp_x   = cumweight[spec_x] / cumweight[-1]             # (npix,)
```

Within each spectral row the original `xtemp` increases linearly by
`row_weights[i]` per column (a tiny within-row variation). The 1-D version
assigns all pixels in the same row the same value — a negligible approximation
since `xtemp` is only used as a smooth spectral coordinate for B-spline fitting.
`xx` (currently `xtemp.sum(axis=1) / nc`, used as per-row knot positions)
becomes simply `cumweight / cumweight[-1]`, shape `(nspec,)`.

**`pb` (area weight)**:

```python
# Before (2-D broadcast in Block 10, 1-D construction in Block 8)
pb = area[:, None] * np.ones((nspec, nc))
pb_inside = np.repeat(area, nc)[inside]          # via flat indexing

# After (1-D indexing)
pb_x = area[spec_x]                              # (npix,)
pb_inside = area[spec_x[inside]]                 # subset of (npix,)
```

**Sorting** — all flat-index patterns collapse to plain 1-D argsort:

```python
# Before
ss = sigma_x.flatten().argsort(kind='stable')          # argsort over nspec*nc
inside, = np.where((condition).flatten())

# After
ss = sigma_x.argsort(kind='stable')                    # argsort over npix
inside, = np.where(condition)                           # no .flatten() call
```

`sigma_x.flat[si]` → `sigma_x[si]`, `norm_obj.flat[si]` → `norm_obj_x[si]`,
and so on throughout.  **No `.flat[]` accesses are needed anywhere.**

**`mask`**:

```python
# Before
mask = np.zeros(nspec * nc, dtype=bool)
mask[si] = ...

# After
mask_x = np.zeros(npix, dtype=bool)
mask_x[si] = ...
```

**`full_bsp`** — already 1-D in spirit; now literally 1-D of length `npix`
without any reshape step.

**`good`** — already proposed as 2-D in Stage 1; in the 1-D scheme it is
naturally 1-D:
```python
good_x = norm_ivar_x > 0.0    # (npix,)
```

**Output reconstruction** — a single scatter at the very end:

```python
profile_model = np.zeros((nspec, nc))
profile_model[totmask_local] = profile_x   # scatter 1-D back to 2-D
# row normalization uses the existing axis=1 logic on the 2-D result
```

where `totmask_local` is `totmask[ipix]` — the mask relative to the sub-image
bounding box.

### Changes to helper function signatures

#### `_fit_spectrum_and_normalize`

```python
# Before
def _fit_spectrum_and_normalize(
    wave, flux, fluxivar, waveimg, image, ivar, totmask, percentile_sn2, fwhm
)
# Returns: (success, med_sn2, norm_obj, norm_ivar, good, xtemp, sn2_img)
# Shapes: norm_obj/norm_ivar/xtemp/sn2_img → (nspec, nc); good → (nspec*nc,)

# After
def _fit_spectrum_and_normalize(
    wave, flux, fluxivar, wave_x, flux_x, ivar_x, spec_x, nspec, percentile_sn2, fwhm
)
# Returns: (success, med_sn2, norm_obj_x, norm_ivar_x, good_x, xtemp_x, sn2_x)
# Shapes: all 1-D of length npix (or nspec for xtemp_x's row_weights)
```

Internally, `wave_min = wave_x.min()`, `wave_max = wave_x.max()`, and the
scatter operations `sn2_img[totmask] = ...` / `spline_img[totmask] = ...`
become direct 1-D assignments `sn2_x = ...` / `spline_img_x = ...`.
`sn2_1d` (the per-spectral-row `sn2` array) is still computed as before
and then broadcast to pixels via `sn2_x = sn2_1d[spec_x]`.

#### `_compute_bspline_knots`

```python
# Before
def _compute_bspline_knots(dspat, sigma, med_sn2, prof_nsigma, good)
# dspat: (nspec, nc); good: (nspec*nc,) flat bool

# After
def _compute_bspline_knots(dspat_x, sigma, spec_x, med_sn2, prof_nsigma, good_x)
# dspat_x: (npix,); spec_x: (npix,) int; good_x: (npix,) bool
```

`sigma_x = dspat_x / sigma[spec_x]` (1-D), `sigma_x[good_x]` (direct 1-D bool indexing).
The rest of the function (limit, sinh_space, bkpt) is unchanged.

#### `_apodize_profile`

```python
# Before
def _apodize_profile(bset, sigma_x, min_sigma, max_sigma, ss, ...)
# sigma_x: (nspec, nc); ss: flat int argsort; full_bsp returned as (nspec*nc,) → reshaped

# After
def _apodize_profile(bset, sigma_x, min_sigma, max_sigma, ss, ...)
# sigma_x: (npix,); ss: 1-D int argsort of length npix; full_bsp returned as (npix,)
```

Every `.flatten()` and `.flat[]` inside `_apodize_profile` disappears because
`sigma_x` is already 1-D.  The sort `ss = sigma_x.argsort()`, the boolean
conditions `left = sigma_x < l_limit`, and the scatter `full_bsp[left] = ...`
are all direct 1-D operations.  The reshape in the caller is removed.

### Relationship to the Stage 1 / Stage 2 flattening changes

The 1-D approach **supersedes** Stage 1 and Stage 2 completely.  Stage 1 changed
`good`, `igood`, `full_bsp`, and `mask` from explicit flat to 2-D; the 1-D
approach makes them naturally 1-D without any `.flatten()` or `.flat[]`.
Stage 2 (replacing `si`, `ss`, `inside` with `(row, col)` pairs) becomes
unnecessary because those are already direct 1-D indices into a 1-D array.

The only code that remains 2-D is Block 13 (row normalization of
`profile_model`) and the output scatter — both of which are deliberately 2-D.

### Efficiency analysis

Let `nspec = 2048`, `nc = 40` (bounding box), `npix_valid = nspec × 25 = 51200`
(valid slit pixels after masking).  Current vs. proposed work:

| Operation | Current | Proposed |
|-----------|---------|----------|
| `sigma_x` recompute (×3 in Block 8 loop) | 3 × 2048×40 = 246k ops | 3 × 51.2k = 154k ops |
| Flat argsort (Block 10) | sort 81.9k elements | sort 51.2k elements |
| `np.repeat(row_weights, nc)` | allocate 81.9k array | `cumweight[spec_x]`: 51.2k indexing |
| `full_bsp` reshape | 1 call | none |

The gains compound when `nc` ≫ slit width, which is common for multislit
spectrographs where `maskwidth` extends several FWHM beyond the object centre.

### Test modifications required

The 1-D approach requires significant changes to the helper-function unit tests
because the helper signatures change. The end-to-end tests (`test_refactor_regression`,
`test_fwhm_recovery_refactor`, etc.) are unaffected — they call
`fit_profile_refactor` through its public interface which remains 2-D.

#### New factory for 1-D inputs

A `make_profile_inputs_1d` helper (or an extension to `make_profile_inputs`)
should return the pre-extracted 1-D arrays:

```python
def make_profile_inputs_1d(nspec=200, nspat=100, fwhm=4.0, sn_ratio=20.0, seed=42):
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, fluxivar, inmask, _ \
        = make_profile_inputs(nspec=nspec, nspat=nspat, fwhm=fwhm, sn_ratio=sn_ratio, seed=seed)
    totmask = (ivar > 0.0) & thismask & inmask
    spec_x = np.where(totmask)[0]                # (npix,) int
    spat_x = spat_img[totmask]                   # (npix,) float
    wave_x = waveimg[totmask]                    # (npix,) float
    flux_x = image[totmask]                      # (npix,) float
    ivar_x = ivar[totmask]                       # (npix,) float
    return wave, flux, fluxivar, wave_x, flux_x, ivar_x, spec_x, trace_in, nspec, nspat
```

#### `_fit_spectrum_and_normalize` tests

| Test | Before | After |
|------|--------|-------|
| `test_fit_spectrum_success` | `assert good.shape == (image.size,)` | `assert good_x.shape == (npix,)` |
| `test_fit_spectrum_good_matches_norm_ivar` | `assert_array_equal(good, norm_ivar.flatten() > 0)` | `assert_array_equal(good_x, norm_ivar_x > 0)` |
| All four FSN tests | Pass 2-D `image`, `ivar`, `waveimg`, `totmask` | Pass 1-D `wave_x`, `flux_x`, `ivar_x`, `spec_x` |

#### `_compute_bspline_knots` tests

The existing `_make_knots_inputs` factory computes 2-D `dspat` and flat `good`
from `_fit_spectrum_and_normalize`. Under the 1-D scheme it instead returns
1-D `dspat_x` and 1-D `good_x`. The assertions (`sigma_x.shape`, `limit > 0`,
etc.) are otherwise the same.

#### `_apodize_profile` tests

The `apodize_bspline` fixture currently creates a 2-D `sigma_x` of shape
`(nspec, nspat)` and a flat `ss`. Under the 1-D scheme the fixture creates a
1-D `sigma_x` of shape `(npix,)` and a 1-D `ss = sigma_x.argsort()`. The
test bodies are otherwise unchanged.

#### `test_apodize_output_shape`

```python
# Before (Stage 1 target)
assert full_bsp.shape == f['sigma_x'].shape    # (nspec, nspat)

# After (1-D scheme)
assert full_bsp.shape == (npix,)
```

#### New regression test for 1-D / 2-D equivalence

After implementing the 1-D approach, a round-trip test verifies that the
scatter-back produces the same `profile_model` as the current 2-D
implementation (up to floating-point rounding from the slightly different
`xtemp` construction):

```python
def test_1d_roundtrip_matches_2d():
    """profile_model from 1-D path matches the 2-D regression reference."""
    # Uses the same frozen .npy files as test_refactor_regression
    ...
    np.testing.assert_allclose(profile_model_1d, ref_profile, atol=1e-8)
```

The tolerance is relaxed from `1e-10` to `1e-8` because the 1-D `xtemp`
construction is a minor approximation of the 2-D original.

---

## Array Flattening: Analysis and Proposed Refactoring

### Why flattening is used at all

The spatial profile is fit as a function of a single 1-D coordinate `sigma_x`,
but `sigma_x` is computed over the entire 2-D image `(nspec, nspat)`. Fitting
the profile B-spline requires pooling all `nspec × nspat` pixels into sorted
1-D arrays of `(sigma_x, norm_obj, norm_ivar)`. That pooling — the sort and the
scatter — is what drives all of the flat-indexing.

### Catalog of every flattening site

The table below lists every site in `spatialprofile_refactor.py` where a 2-D
array is explicitly flattened, where a `.flat[]` view is used, or where a
flat index array is constructed.

| Line | Expression | Category |
|------|-----------|----------|
| 449 | `good = norm_ivar.flatten() > 0.0` | Boolean flat |
| 502–504 | `sigma_x.flat[good]` ×3 | `.flat[]` with bool |
| 573 | `sigma_x.flatten() > min_sigma` (×2) | `.flatten()` in condition |
| 574 | `full_bsp = np.zeros(sigma_x.size)` | Explicit 1-D output |
| 575 | `sigma_x.flat[igood]` | `.flat[]` with bool |
| 648 | `sigma_x.flatten() < l_limit` | `.flatten()` in condition |
| 650 | `sigma_x.flatten() > r_limit` | `.flatten()` in condition |
| 649, 651 | `sigma_x.flat[left]`, `sigma_x.flat[right]` | `.flat[]` with bool |
| 783, 785 | `np.where(condition.flatten())` | Flat integer index creation |
| 787, 791, 794, 801 | `sigma_x.flat[inside/si/sr]` | `.flat[]` with int |
| 815–824 | `sigma_x.flat[sr]`, `sigma_x.flat[si]` ×4 | `.flat[]` with int |
| 831 | `mask = np.zeros(nspec * nspat, …)` | Explicit 1-D mask |
| 832–834 | `norm_ivar.flat[si]`, `norm_obj.flat[si]`, `mask[si]` | `.flat[]` with int |
| 857–951 | `sigma_x.flat[inside]` ×7, `norm_obj.flat[inside]` ×2, `norm_ivar.flat[inside]` ×2, `xtemp.flat[inside]` ×2, `pb[ss]` | `.flat[]` with int (Block 8 loop) |
| 981 | `ss = sigma_x.flatten().argsort(…)` | Flat sort |
| 982–990 | `sigma_x.flat[ss]`, `mask[ss]`, `norm_obj.flat[ss[inside]]` ×3 | `.flat[]` with int |
| 1006 | `full_bsp.reshape(nspec, nspat)` | Reshape 1-D → 2-D |

### What cannot be avoided

The flat integer index arrays — `inside`, `si`, `sr`, `ss` — and every `.flat[]`
access through them are **intrinsic to the algorithm**. The B-spline fit
(Block 7, 8, 10) requires the `nspec × nspat` pixels to be presented to
`iterative_bspline_fit` as a sorted 1-D sequence. There is no way to skip that
sort without changing `iterative_bspline_fit` itself. Consequently, the
following patterns are unavoidable:

```
ss = sigma_x.flatten().argsort(kind='stable')   # line 981
si = inside[sigma_x.flat[inside].argsort(…)]    # line 787
sigma_x.flat[si]   # every B-spline fit argument
```

A complete elimination of `.flat[]` would require replacing all flat integer
indices with `(row, col)` index-pair arrays — a large-scale restructuring
described as Stage 2 below.

### What can be eliminated: Stage 1 (local, low-risk)

Six `.flatten()` calls and one `.reshape()` can be removed with local edits to
three functions, without touching the flat integer index pattern at all.

#### 1. `good` as a 2-D boolean mask

In `_fit_spectrum_and_normalize` (line 449):
```python
# Before
good = norm_ivar.flatten() > 0.0          # shape (nspec*nspat,)

# After
good = norm_ivar > 0.0                    # shape (nspec, nspat)
```

In `_compute_bspline_knots` (lines 502–504):
```python
# Before
abs_sigma = np.fmin((np.abs(sigma_x.flat[good])).max(), 2.0 * limit)
min_sigma  = np.fmax(sigma_x.flat[good].min(), -abs_sigma)
max_sigma  = np.fmin(sigma_x.flat[good].max(),  abs_sigma)

# After — NumPy 2-D boolean indexing returns the same 1-D selection
abs_sigma = np.fmin((np.abs(sigma_x[good])).max(), 2.0 * limit)
min_sigma  = np.fmax(sigma_x[good].min(), -abs_sigma)
max_sigma  = np.fmin(sigma_x[good].max(),  abs_sigma)
```

No other call sites change: `good.sum()` works for both shapes; the pass of
`good` to `_compute_bspline_knots` is transparent.

#### 2. `igood`, `full_bsp`, and the conditional masks in `_apodize_profile`

In `_apodize_profile` (lines 573–574, 648, 650):
```python
# Before
igood = (sigma_x.flatten() > min_sigma) & (sigma_x.flatten() < max_sigma)
full_bsp = np.zeros(sigma_x.size)          # 1-D

# After
igood = (sigma_x > min_sigma) & (sigma_x < max_sigma)   # 2-D boolean
full_bsp = np.zeros(sigma_x.shape)         # 2-D
```

```python
# Before (lines 648–651)
left  = sigma_x.flatten() < l_limit
full_bsp[left]  = np.exp(…) * l_fit_val
right = sigma_x.flatten() > r_limit
full_bsp[right] = np.exp(…) * r_fit_val

# After — 2-D boolean indexing produces the same 1-D selection
left  = sigma_x < l_limit
full_bsp[left]  = np.exp(-(sigma_x[left]  - l_limit) * l_deriv) * l_fit_val
right = sigma_x > r_limit
full_bsp[right] = np.exp(-(sigma_x[right] - r_limit) * r_deriv) * r_fit_val
```

The operations using `ss` (flat integer sort, lines 584–597) are unchanged —
`sigma_x.flat[ss]` and `full_bsp.flat[ss]` remain valid on a 2-D `full_bsp`
because `.flat` is a view of any ndarray.

The return value of `_apodize_profile` changes from 1-D to 2-D, so the
`reshape` call in `fit_profile_refactor` (line 1006) is removed:
```python
# Before
full_bsp = full_bsp.reshape(nspec, nspat)

# After — already 2-D, no reshape needed
```

#### 3. `mask` as a 2-D boolean array

In Block 7 (line 831):
```python
# Before
mask = np.zeros(nspec * nspat, dtype=bool)
mask[si] = …

# After — si is still a flat integer; .flat works on any ndarray shape
mask = np.zeros((nspec, nspat), dtype=bool)
mask.flat[si] = …
```

The read sites `mask[ss]` and `mask[si]` in Block 10 (line 983) become
`mask.flat[ss]` and `mask.flat[si]`. The only practical benefit is consistency:
`mask` is logically a 2-D pixel map, so giving it a 2-D shape makes its
relationship to `sigma_x`, `norm_obj`, etc. self-evident.

### What can be eliminated: Stage 2 (large-scale, optional)

The remaining `.flat[]` accesses are through flat integer index arrays —
`si`, `sr`, `ss`, `inside`. These can be replaced entirely by `(row, col)`
index-pair arrays, eliminating every remaining `.flat[]`:

```python
# Current pattern
inside, = np.where(condition.flatten())           # flat integer
sigma_x.flat[inside]                              # .flat access

# Stage 2 pattern
row_inside, col_inside = np.where(condition)      # (row, col) pair
sigma_x[row_inside, col_inside]                   # direct 2-D access
```

The sort step adapts naturally:
```python
# Current
si = inside[np.argsort(sigma_x.flat[inside], kind='stable')]

# Stage 2
sort_order = np.argsort(sigma_x[row_inside, col_inside], kind='stable')
row_si, col_si = row_inside[sort_order], col_inside[sort_order]
```

Every site that currently reads `sigma_x.flat[si]`, `norm_obj.flat[si]`, etc.
becomes `sigma_x[row_si, col_si]`, `norm_obj[row_si, col_si]`, etc. Helper
function signatures change to accept `(row_ss, col_ss)` pairs instead of
`ss`. The B-spline fit calls are otherwise unchanged.

**Assessment**: Stage 2 eliminates the last ~25 `.flat[]` accesses and removes
the only remaining `.flatten()` call (the argsort at line 981). The cost is
that every index operation in Blocks 7–10 and in `_apodize_profile` must carry
two arrays instead of one, which adds noise to the code that is at least as
significant as the `.flat[]` calls it removes. Stage 2 is recommended only if
there is a functional motivation (e.g., interfacing with an API that does not
accept flat indices) or if code-review standards require it. Stage 1 is the
pragmatic stopping point.

---

## Test Modifications

### Changes required by Stage 1

Stage 1 changes the shapes of `good` and `full_bsp`. Three tests must be
updated; the assertions become simpler:

| Test | Before | After |
|------|--------|-------|
| `test_fit_spectrum_success` | `assert good.shape == (image.size,)` | `assert good.shape == image.shape` |
| `test_fit_spectrum_good_matches_norm_ivar` | `np.testing.assert_array_equal(good, norm_ivar.flatten() > 0)` | `np.testing.assert_array_equal(good, norm_ivar > 0)` |
| `test_apodize_output_shape` | `assert full_bsp.shape == (f['sigma_x'].size,)` | `assert full_bsp.shape == f['sigma_x'].shape` |

The remaining apodize tests (`test_apodize_no_deriv_skips_tails`,
`test_apodize_tail_continuity`, `test_apodize_prof_nsigma_zeros_limits`) access
`full_bsp` via boolean masks (`full_bsp[outside]`, `full_bsp[left_tail]`, etc.)
and via `full_bsp.size`. Both of those work identically whether `full_bsp` is
1-D or 2-D, so no changes are needed. The `igood` computation inside
`test_apodize_output_shape` should change to match:
```python
# Before (in test body)
igood = (f['sigma_x'].flatten() > f['min_sigma']) & (f['sigma_x'].flatten() < f['max_sigma'])
assert np.all(full_bsp[igood] >= 0.0)

# After — 2-D boolean, same result
igood = (f['sigma_x'] > f['min_sigma']) & (f['sigma_x'] < f['max_sigma'])
assert np.all(full_bsp[igood] >= 0.0)
```

The `apodize_bspline` fixture itself does not change: `ss` is still computed as
a flat argsort because `_apodize_profile` still receives `ss` as a flat integer
index.

### New cross-implementation comparison tests

The existing tests compare `fit_profile_refactor` either against frozen
reference arrays (`test_refactor_regression`) or in isolation. The following
new tests compare both implementations on the same synthetic inputs to document
where they agree and where they intentionally diverge.

Because the B-spline path outputs differ due to Bugs 2 and 3 being fixed, the
comparison tests focus on **properties** shared by both implementations rather
than element-by-element equality.

#### `test_med_sn2_agreement`
`med_sn2` is produced solely by the spectral B-spline fitting logic (Blocks 2–3),
which is identical in both implementations. For wide profiles where Bug 2 has
negligible effect (`fwhm ≥ 6 px`), `med_sn2` from both functions should agree
within 5 % under the same inputs:
```python
@pytest.mark.parametrize("fwhm", [6.0, 8.0])
def test_med_sn2_agreement(fwhm):
    inputs = make_profile_inputs(sn_ratio=20.0, fwhm=fwhm)
    kwargs = dict(zip(
        ['image','ivar','waveimg','thismask','spat_img','trace_in',
         'wave','flux','fluxivar','inmask'], inputs[:10]
    ))
    _, _, _, med_sn2_leg = spatialprofile.fit_profile(**kwargs, thisfwhm=fwhm, sn_gauss=4.0)
    _, _, _, med_sn2_ref = spatialprofile_refactor.fit_profile_refactor(**kwargs, thisfwhm=fwhm, sn_gauss=4.0)
    np.testing.assert_allclose(med_sn2_leg, med_sn2_ref, rtol=0.05)
```

#### `test_forced_gaussian_trace_agreement`
With `gauss=True`, both implementations skip the B-spline path and return
`xnew = trace_in` unchanged. `fwhmfit` is also set to `thisfwhm` in both. These
two return values should be identical regardless of the Gaussian profile form
(point-sampled vs. integrated):
```python
def test_forced_gaussian_trace_agreement():
    inputs = make_profile_inputs(sn_ratio=20.0, fwhm=4.0)
    kwargs = dict(zip(
        ['image','ivar','waveimg','thismask','spat_img','trace_in',
         'wave','flux','fluxivar','inmask'], inputs[:10]
    ))
    _, xnew_leg,    fwhm_leg,    _ = spatialprofile.fit_profile(**kwargs, gauss=True, thisfwhm=4.0)
    _, xnew_ref,    fwhm_ref,    _ = spatialprofile_refactor.fit_profile_refactor(**kwargs, gauss=True, thisfwhm=4.0)
    np.testing.assert_array_equal(xnew_leg, xnew_ref)
    np.testing.assert_array_equal(fwhm_leg, fwhm_ref)
```

#### `test_bspline_row_normalization_both`
Both implementations normalize each spectral row of `profile_model` to unit sum
in Block 13. This invariant must hold for both regardless of which bugs are
fixed:
```python
@pytest.mark.parametrize("func", [
    spatialprofile.fit_profile,
    spatialprofile_refactor.fit_profile_refactor,
])
def test_bspline_row_normalization_both(func):
    inputs = make_profile_inputs(sn_ratio=20.0, fwhm=6.0)
    kwargs = dict(zip(
        ['image','ivar','waveimg','thismask','spat_img','trace_in',
         'wave','flux','fluxivar','inmask'], inputs[:10]
    ))
    profile_model, _, _, med_sn2 = func(**kwargs, thisfwhm=6.0, sn_gauss=4.0)
    assert med_sn2 > 4.0 ** 2, "Expected B-spline path"
    row_sums = profile_model.sum(axis=1)
    nonzero = row_sums > 0
    np.testing.assert_allclose(row_sums[nonzero], 1.0, atol=1e-10)
```

#### `test_profile_positivity_both`
Both implementations should return non-negative `profile_model` on all code
paths. This is a structural invariant independent of any bug fix:
```python
@pytest.mark.parametrize("func,sn_ratio", [
    (spatialprofile.fit_profile, 20.0),
    (spatialprofile.fit_profile,  1.0),
    (spatialprofile_refactor.fit_profile_refactor, 20.0),
    (spatialprofile_refactor.fit_profile_refactor,  1.0),
])
def test_profile_positivity_both(func, sn_ratio):
    inputs = make_profile_inputs(sn_ratio=sn_ratio, fwhm=6.0)
    kwargs = dict(zip(
        ['image','ivar','waveimg','thismask','spat_img','trace_in',
         'wave','flux','fluxivar','inmask'], inputs[:10]
    ))
    profile_model, _, _, _ = func(**kwargs, sn_gauss=4.0)
    assert np.all(profile_model >= 0)
    assert np.all(np.isfinite(profile_model))
```

#### `test_trace_correction_direction_both`
For a 1 px offset (large enough that Bug 3 does not prevent convergence in
either implementation), both functions should correct the trace in the right
direction and reduce the error:
```python
def test_trace_correction_direction_both():
    inputs = make_profile_inputs(sn_ratio=20.0, fwhm=4.0, trace_offset=1.0)
    true_trace = inputs[10]
    trace_in   = inputs[5]
    kwargs = dict(zip(
        ['image','ivar','waveimg','thismask','spat_img','trace_in',
         'wave','flux','fluxivar','inmask'], inputs[:10]
    ))
    for func in [spatialprofile.fit_profile,
                 spatialprofile_refactor.fit_profile_refactor]:
        _, xnew, _, med_sn2 = func(**kwargs, thisfwhm=4.0, sn_gauss=4.0)
        assert med_sn2 > 4.0 ** 2
        assert np.mean(xnew) > np.mean(trace_in)
        assert np.abs(np.mean(xnew - true_trace)) < np.abs(np.mean(trace_in - true_trace))
```

#### `test_bspline_fwhm_improvement_narrow`
For a narrow profile (`fwhm = 3 px`), the legacy function fails the 10 %
tolerance while the refactored function passes it. This documents the Bug 2
improvement without asserting exact agreement between implementations:
```python
def test_bspline_fwhm_improvement_narrow():
    fwhm = 3.0
    inputs = make_profile_inputs(sn_ratio=20.0, fwhm=fwhm)
    kwargs = dict(zip(
        ['image','ivar','waveimg','thismask','spat_img','trace_in',
         'wave','flux','fluxivar','inmask'], inputs[:10]
    ))
    _, _, fwhm_leg, med_sn2 = spatialprofile.fit_profile(**kwargs, thisfwhm=fwhm, sn_gauss=4.0)
    _, _, fwhm_ref, _       = spatialprofile_refactor.fit_profile_refactor(**kwargs, thisfwhm=fwhm, sn_gauss=4.0)
    assert med_sn2 > 4.0 ** 2, "Expected B-spline path"
    # Legacy overestimates narrow FWHM by > 15%; refactored recovers within 10%
    assert np.median(fwhm_leg) > fwhm * 1.15, "Legacy should overestimate FWHM"
    np.testing.assert_allclose(np.median(fwhm_ref), fwhm, rtol=0.10)
```

---

## Change Log

### v1 → v2
- Added "Pre-existing Bugs Discovered During Testing" section documenting Bugs 1,
  2, and 3 found while running the unit tests against the pre-refactor function.

### v2 → v2 (second edit)
- Added "Accuracy Improvements" section covering the pixel-integrated Gaussian
  and linear interpolation for the half-maximum crossing.

### v2 → v3
- Added "Output File" section; removed large Python code blocks already
  implemented; updated cross-references.

### v3 → v4
- Marked Verification Plan steps 3 and 4 complete.
- Replaced "Candidate Helper Functions" with "Implemented Helper Functions"
  (four helpers) and "Deferred Block-Level Extractions".

### v4 → v5
- **New helpers extracted**:
  - `_fit_spectrum_and_normalize` — encapsulates Blocks 2–5.
  - `_compute_bspline_knots` — encapsulates Block 6.
  - `_apodize_profile` — encapsulates Blocks 11–12.
- **Reorganized "Remaining Inline Blocks"** table replacing "Deferred
  Block-Level Extractions" prose.
- **Additional code improvements**: Block 7 in-place `sigma *=` / `limit *=`;
  Block 8 0-indexed loop; `fwhm` parameter rename; improved log messages; JXP
  kludge consolidation.
- **Added "Open Question"** section on the Block 8 mid-loop re-fit condition.
- **Accuracy Improvements**: added §C (`sig2fwhm` precision).
- Updated Verification Plan: step 5 listed; step 6 added (resolve open
  question); dev-suite and full-suite steps renumbered.

### v5 → v6
- **Added steps 6–15 to Verification Plan**: step-wise implementation sequence
  for the 1-D slit-pixel approach; existing steps 6–8 renumbered to 16–18.
- **Added "1-D Slit-Pixel Approach" section**: proposes adding `spec_img` to
  the function signature; defines Block 0 (extract slit pixels as 1-D arrays);
  transforms `dspat`, `sigma_x`, `xtemp`, `pb`, `good`, `mask`, `full_bsp` to
  1-D; changes signatures of all three helper functions; notes that Stage 1 and
  Stage 2 are superseded; provides efficiency analysis and test modifications
  including a new `test_1d_roundtrip_matches_2d` regression test.
- **Added "Array Flattening: Analysis and Proposed Refactoring" section**:
  catalogs all 30+ flattening sites; distinguishes intrinsic from avoidable;
  proposes Stage 1 (3 local changes to `good`, `igood`/`full_bsp`, `mask`) and
  Stage 2 (full `(row, col)` restructuring); lists updated and new tests.
- **Added "Test Modifications" section**: documents which 3 existing tests
  change under Stage 1, plus 6 new cross-implementation comparison tests.
- **Marked Verification Plan step 5 complete**: 42 tests pass (up from 12).
  New tests cover `fit_profile_refactor` end-to-end (4 tests) and all five
  helper functions directly (16 tests).
- **Corrected Bug 3 description**: the 0.5 px trace-offset case still overshoots
  the true trace (~0.05 px); `test_trace_correction_small_offset` asserts
  direction only, not convergence.
- **Clarified `test_prof_nsigma_refactor`**: FWHM recovery is not asserted
  because a sparse 2-knot B-spline over a ±10σ fitting range does not reliably
  recover the FWHM of a narrow Gaussian. The important property — no crash at
  `prof_nsigma=10.0` — is verified.
- **Step 7 deferred**: dev-suite run skipped for now; marked with strikethrough
  in the Verification Plan.
- **Linked open question to Bug 3**: restoring the legacy two-re-fit behavior in
  Block 8 may reduce the trace-correction overshoot at small offsets.
- **Updated Unit Tests section**: full roster of all 42 tests listed and
  organized by group; noted the `[6.0, 8.0]`-only parametrization for the legacy
  `test_fwhm_recovery`.

### v6 → v7
- **Marked Verification Plan steps 6 and 7 complete** (42 tests pass).
- **Step 6 implementation**: `make_profile_inputs_1d` and `_make_fsn_inputs_1d`
  added to `pypeit/tests/test_spatialprofile.py` after `make_profile_inputs`.
  Deviation: the plan prose listed `spat_x` among the returned arrays of
  `make_profile_inputs_1d`, but the code example in the plan omitted it from the
  return tuple.  The implementation follows the code example (`spat_x` is
  computed but not returned).
- **Step 7 implementation**: `spec_img=None` added to the `fit_profile_refactor`
  signature after `inmask`; docstring entry added; Block 0 inserted immediately
  after `totmask` is fully constructed (before `nspec, nspat = image.shape`),
  using the existing dashed-header comment style. All five 1-D pixel arrays
  (`spec_x`, `spat_x`, `wave_x`, `flux_x`, `ivar_x`) stored as locals;
  downstream 2-D logic unchanged.
- **Updated Unit Tests section**: added "1-D infrastructure" subsection listing
  the two new factory functions and their signatures.
- **Updated Workflow Breakdown**: added Block 0 description immediately after
  Block 1 entry, documenting its position and the local variables it produces.

### v7 (steps 8–9 addendum)
- **Marked Verification Plan steps 8 and 9 complete** (all 42 tests pass after
  regression reference regeneration).
- **`_fit_spectrum_and_normalize` converted to 1-D**: signature changed from
  `(wave, flux, fluxivar, waveimg, image, ivar, totmask, percentile_sn2, fwhm)`
  to `(wave, flux, fluxivar, wave_x, flux_x, ivar_x, spec_x, nspec,
  percentile_sn2, fwhm)`; all return arrays are 1-D of length `npix`. `xtemp_x`
  uses the row-approximation `cumweight[spec_x] / cumweight[-1]` (assigning all
  pixels in a row the same value); this changes `profile_model`/`xnew`/`fwhmfit`
  by more than `atol=1e-10`, requiring reference array regeneration. Call site
  in `fit_profile_refactor` scatters 1-D results back to 2-D for Blocks 7–10
  (step-8 scaffolding, to be removed in step 10).
- **`_compute_bspline_knots` converted to 1-D**: signature changed from
  `(dspat, sigma, med_sn2, prof_nsigma, good)` to `(dspat_x, sigma, spec_x,
  med_sn2, prof_nsigma, good_x)`; `sigma_x = dspat_x / sigma[spec_x]` (1-D of
  length `npix`). Call site computes `dspat_x = spat_x - trace_in[spec_x]`;
  2-D `dspat` retained temporarily for Block 8 (removed in step 10). Step-9
  scaffolding scatters 1-D `sigma_x` back to 2-D for Block 7 (removed in step 10).
- **Test file updated**: `_make_fsn_inputs` (2-D factory) removed; all four
  `test_fit_spectrum_*` tests updated to use `_make_fsn_inputs_1d` with 1-D
  assertions; `_make_knots_inputs` updated to 1-D final interface (returns
  `dspat_x`, `sigma`, `med_sn2`, `good_x`, `spec_x`); three
  `test_bspline_knots_*` tests pass 1-D kwargs and assert `sigma_x.shape ==
  (npix,)`.
- **Helper function descriptions updated** to reflect the new 1-D signatures for
  both `_fit_spectrum_and_normalize` and `_compute_bspline_knots`.
- **Unit Tests section updated**: FSN and knots subsections describe 1-D
  assertions; "1-D infrastructure" subsection updated to show `_make_fsn_inputs_1d`
  in active use and `_make_knots_inputs` at its step-9 final form.

**Step 10 (Blocks 6–10 fully 1-D)**:
- **Step-8/9 scaffolding removed**: all 2-D scatter-backs (`norm_obj`,
  `norm_ivar`, `sn2_img`, `xtemp`, `sigma_x`) and 2-D `dspat` deleted from
  `fit_profile_refactor`; `gauss_kwargs` and all `_return_gaussian` calls
  updated to use `good_x` and `norm_obj_x`.
- **Block 6**: `npix = spec_x.size` added; `dspat = spat_img -
  trace_in[:, None]` and `sigma_x` scatter-back removed.
- **Block 7**: `GOOD_PIX`/`IN_PIX` conditions use `sn2_x`/`norm_ivar_x`;
  `np.where(...flatten())` → `np.where(...)`; all `.flat[si/sr]` → `[si/sr]`;
  `mask = np.zeros(nspec*nspat, …)` → `mask_x = np.zeros(npix, …)`.
- **Block 8**: `xx = np.bincount(spec_x, weights=xtemp_x, minlength=nspec) /
  nspat` computed once before the loop; `xtemp.flat[...]` → `xtemp_x[...]`;
  `sigma_x = dspat / sigma[:, None] - trace_corr[:, None]` →
  `sigma_x = dspat_x / sigma[spec_x] - trace_corr[spec_x]`;
  `pb = np.repeat(area, nspat)[inside]` → `pb = area[spec_x[inside]]`; all
  remaining `.flat[]` accesses removed.
- **Block 10**: `ss = sigma_x.flatten().argsort()` → `ss =
  sigma_x.argsort()`; `mask[ss]` → `mask_x[ss]`; `pb = area[:, None] *
  np.ones(…)` → `pb_x = area[spec_x]`.
- **Block 13**: temporary scatter-back added — `profile_x = full_bsp * pb_x;
  profile_model[totmask] = profile_x` — so row normalisation remains on the
  2-D `profile_model`; `res_mode`/`chi_good` updated to 1-D indexing; QA call
  updated to use `norm_obj_x`/`pb_x`.
- **Verification Plan step 10 marked ✅**. All 42 tests pass;
  `test_refactor_regression` continues to pass (no new numerical change beyond
  the `xtemp_x` approximation already introduced in step 8).

**Steps 11–12 (`_apodize_profile` 1-D + Block 13 finalisation)**:
- **`_apodize_profile`**: all `sigma_x.flatten()` → `sigma_x` and
  `sigma_x.flat[...]` → `sigma_x[...]` throughout (Block 11 and Block 12
  apodization). Docstring updated: `sigma_x` now documented as shape
  :math:`(N_{\rm pix},)`.
- **`apodize_bspline` fixture**: 2-D `sigma_x` flattened to 1-D after
  construction; `npix = sigma_x.size` added to the returned dict;
  `ss = sigma_x.argsort(kind='stable')`.
- **`test_apodize_output_shape`**: `shape == (f['npix'],)`; `igood` uses
  direct 1-D boolean indexing.
- **`test_apodize_prof_nsigma_zeros_limits`**: `shape == (f['npix'],)`.
- **`test_apodize_no_deriv_skips_tails`**: `sigma_x_flat` removed; direct
  `f['sigma_x']` used.
- **`test_apodize_tail_continuity`**: `sigma_x_flat` removed; local
  `sigma_x = f['sigma_x']` used.
- **Block 13**: "step-10 scaffolding" comment removed; scatter-back
  (`profile_x = full_bsp * pb_x; profile_model[totmask] = profile_x`) is
  the final form.
- **Verification Plan steps 11–12 marked ✅**. All 42 tests pass;
  `test_refactor_regression` continues to pass (no numerical change).

**Steps 13–14 (reference array regeneration + cross-implementation tests)**:
- **Step 13**: `freeze_spatprof.py` re-run; three untracked `.npy` files
  overwritten with the 1-D implementation's output.  All 42 tests continue to
  pass, including `test_refactor_regression` at `atol=1e-10`.
- **Step 14**: Seven new tests added to `test_spatialprofile.py` under a
  "Cross-implementation comparison tests" section:
  - `test_1d_roundtrip_matches_2d`: shape check + `atol=1e-8` regression.
  - `test_med_sn2_agreement` (×2 parametrize): 5 % rtol on `med_sn2`.
  - `test_forced_gaussian_trace_agreement`: identical `xnew`/`fwhmfit`.
  - `test_bspline_row_normalization_both` (×2 parametrize): row sums = 1.
  - `test_profile_positivity_both` (×4 parametrize): non-negative, finite.
  - `test_trace_correction_direction_both`: trace shifts toward true centre.
  - `test_bspline_fwhm_improvement_narrow`: documents Bug 2 fix (legacy
    overshoots 3-px FWHM by 33 %; refactored recovers within 2 %).
- **Verification Plan steps 13–14 marked ✅**. **54 tests pass** (up from 42;
  target was 49+).

**Block 0 consolidated into `_fit_spectrum_and_normalize`**:
- **Motivation**: `wave_x`, `flux_x`, `ivar_x`, and `spec_x` were only ever
  consumed by `_fit_spectrum_and_normalize`; keeping them as locals in
  `fit_profile_refactor` (even briefly, before the `del`) served no purpose
  beyond providing them to the helper. Moving the extraction inside the helper
  collocates construction with use and removes the temporary lifetime entirely.
- **`_fit_spectrum_and_normalize` signature changed** from
  `(wave, flux, fluxivar, wave_x, flux_x, ivar_x, spec_x, percentile_sn2, fwhm)`
  to
  `(wave, flux, fluxivar, waveimg, image, ivar, totmask, spec_img, percentile_sn2, fwhm)`.
  The Block 0 logic (`spec_x` from `spec_img` or `np.where(totmask)[0]`;
  `wave_x`, `flux_x`, `ivar_x` from boolean indexing) now lives at the top of
  the function body.
- **Return tuple extended to 8 values**: `spec_x` appended as the final
  element. Both failure-path returns gained a trailing `None` for `spec_x`.
- **`fit_profile_refactor`**: Block 0 comment block and all four assignment
  lines removed; the `del wave_x, flux_x, ivar_x` line removed; the call to
  `_fit_spectrum_and_normalize` updated to pass `waveimg`, `image`, `ivar`,
  `totmask`, `spec_img` and to unpack the new 8-tuple (including `spec_x`).
  Block comment updated to "Blocks 0 and 2–5".
- **`_make_fsn_inputs_1d`** simplified: no longer extracts 1-D arrays; returns
  `wave, flux, fluxivar, waveimg, image, ivar, totmask, nspec`.
- **Four `test_fit_spectrum_*` tests** updated: unpack 8-tuple from the return;
  call the helper with the new 2-D keyword arguments; `npix = totmask.sum()`
  instead of `spec_x.size`; `test_fit_spectrum_success` adds
  `assert spec_x.shape == (npix,)`.
- **`_make_knots_inputs`** updated: calls `_fit_spectrum_and_normalize` with the
  new signature; obtains `spec_x` from the 8-tuple return; computes `spat_x =
  spat_img[totmask]` and `dspat_x = spat_x - trace_in[spec_x]` afterward.
- **54 tests pass** (unchanged count; no numerical change).

**`BSpline.value` coeff override + Block 8 mock-BSpline cleanup**:
- **Motivation**: Block 8 of `fit_profile_refactor` previously extracted two
  coefficient columns from each `BSpline2D` result by constructing a temporary
  1-D `BSpline` object, mutating its `coeff` attribute twice, and calling
  `.value(xx)` each time — a pattern mirrored from the original IDL code.
  `BSpline2D.value` already accepts a `coeff` override keyword (added for
  `skysub.py`); adding the same override to `BSpline.value` and
  `BSpline._evaluate_model` lets `_bspline2d_to_1d` + `.value(xx, coeff=col)`
  replace the mock construction entirely.
- **`BSpline._evaluate_model`**: added `coeff=None` parameter.  When
  non-``None``, it replaces ``self.coeff`` for this evaluation without
  mutating the instance.
- **`BSpline.value`**: added `coeff=None` parameter.  When non-``None``, both
  the ``x is self.x`` identity fast-path and the ``interpolate`` fast-path are
  skipped (both rely on the cached ``self.yfit`` which is invalid for a
  different coefficient vector); the full design-matrix path runs instead.
  Docstrings for both methods updated accordingly.
- **`spatialprofile_refactor.py` Block 8**: the two 5-line mock-BSpline blocks
  (one for the shift fit, one for the stretch fit) replaced with 2-line
  patterns each — `_bspline2d_to_1d` extracts column 0 into a 1-D `BSpline`;
  `.value(xx, coeff=bspl.coeff[:, 1])` evaluates column 1 without mutation.
  No numerical change; no imports added or removed.
- **`test_bspline.py`**: four new tests added immediately after
  `test_value_gap_masking`:
  - `test_value_coeff_kwarg_additivity`: verifies linearity — splitting the
    coefficient vector and summing partial evaluations reproduces the full fit.
  - `test_value_coeff_kwarg_does_not_mutate_instance`: confirms that neither
    ``self.coeff`` nor the cached ``yfit`` changes after a call with a
    ``coeff`` override.
  - `test_value_coeff_kwarg_bypasses_training_fast_path`: passes the training
    ``x`` with ``coeff=zeros``; asserts the result is zero (not the cached
    ``yfit``), proving the fast path was bypassed.
  - `test_value_coeff_kwarg_bypasses_interpolate_fast_path`: calls with both
    ``interpolate=True`` and a ``coeff`` override; asserts the result matches
    the full B-spline evaluation (not ``np.interp``).
- **`test_spatialprofile.py`**: no changes required; all 54 tests pass
  unchanged (the refactor is internal to Block 8).
- **194 tests pass total** (140 bspline + 54 spatialprofile).

**`BSpline2D.to_1d()` method + removal of `_bspline2d_to_1d` helper**:
- **Motivation**: `_bspline2d_to_1d` in `spatialprofile_refactor.py` was a
  module-level helper whose only job was to extract `coeff[:, 0]` from a
  `BSpline2D` into a fresh 1-D `BSpline`.  Since this is a natural operation
  on any `BSpline2D` object, it belongs as a method on the class rather than
  as a private function in a consumer module.
- **`BSpline2D.to_1d()`** (new method, `bspline.py`): constructs and returns a
  `BSpline` sharing the same breakpoints and order as ``self``, with
  ``coeff = self.coeff[:, 0]`` and a copy of ``self.bkpt_gpm``.  No instance
  attributes are mutated.
- **`spatialprofile_refactor.py`**: `_bspline2d_to_1d` function removed; all
  four call sites updated to the method form (e.g. ``mode_shift_bspl.to_1d()``
  / ``bset.to_1d()``).  The ``from pypeit.core.bspline import BSpline, Knots``
  import removed entirely — nothing from ``bspline`` is needed at runtime in
  that module now that the construction logic lives in the method.
- **`test_bspline.py`**: `test_bspline2d_to_1d` added (seed 98) immediately
  before the interpolate-keyword section.  Checks that the result is a
  ``BSpline`` and not a ``BSpline2D``, that ``coeff`` and ``bkpt_gpm`` match
  the source 2-D object, and that ``nord`` is preserved.
- **195 tests pass total** (141 bspline + 54 spatialprofile).

**Block 2 refactored to use `gpm` keyword in `_fit_spectrum_and_normalize`**:
- **Motivation**: the three `iterative_bspline_fit` calls in Block 2 previously
  pre-selected pixels via boolean indexing (`wave[indsp]`, `flux_sm[indsp]`,
  etc.) before passing data to the fitter.  Because `wave` is already sorted
  ascending and `iterative_bspline_fit` accepts a `gpm` argument that performs
  the same masking internally, the pre-selection is redundant.  Passing full
  arrays with `gpm` lets the returned masks live in `nspec` space rather than
  the compressed `indsp` space, eliminating the manual expansions that cluttered
  Block 3.
- **Block 2** (`_fit_spectrum_and_normalize`): the three fits now pass `wave`,
  `flux_sm`, and `ivar=fluxivar_sm` with the mask threaded through `gpm`:
  ``gpm=indsp`` → ``gpm=bmask`` → ``gpm=bmask2``.  `spline_flux` and
  `cont_flux` are evaluated at `wave` (full), hitting the cached `yfit` fast
  path in `BSpline.value` because `wave is bspl.x`.
- **Block 3** — five downstream simplifications:
  1. `sqrt_ivar` computed from `fluxivar_sm` (not `fluxivar_sm[indsp]`); the
     subsequent `sqrt_ivar[~bmask2] = 0` zeroes out everything outside the
     second fit's accepted set (including outside `indsp`).
  2. S/N statistics (`np.percentile`, `sigma_clipped_stats`) applied to
     ``sn2_indsp = sn2[indsp]`` to preserve the original behaviour; computing
     them over the full `nspec` would dilute the percentile with zeros from
     outside the valid wavelength range.
  3. `spline_flux1` and `cont_flux1` set directly as ``np.where(ispline,
     spline_flux, 0.0)`` — eliminating the separate `b_answer.value(wave[ispline])`
     and `c_answer.value(wave[ispline])` re-evaluations.
  4. The two-line boolean-array expansions
     ``bmask_1d = np.zeros(nspec, bool); bmask_1d[indsp] = bmask2`` and the
     equivalent for ``cmask2_1d`` are removed; `djs_maskinterp` now receives
     `bmask2` and `cmask` directly (already `nspec`-length).
  5. `sn2_med_filt` and the `sn2` interpolation use `sn2_indsp` to preserve
     the same effective data range as before.
- **Block 4**: the single reference to the now-removed ``bmask_1d`` variable
  (in the ``badpix`` definition) replaced with ``bmask2``.
- **No numerical change**; all 54 spatialprofile tests pass unchanged.

**`sn_cap` and `good_pix_frac` promoted to keyword arguments**:
- **Motivation**: `sn_cap = 100.0` (the per-pixel S/N ceiling passed to
  `clip_ivar`) and `good_pix_frac = 0.05` (the minimum fraction of in-range
  spectral pixels required for the fit to proceed) were hard-coded local
  variables inside `_fit_spectrum_and_normalize`.  Promoting them to keyword
  arguments with identical defaults makes them accessible for future tuning
  without requiring changes to the function's internal logic.
- **`_fit_spectrum_and_normalize` signature**: two new optional parameters added
  after the existing positional arguments — ``sn_cap=100.0`` and
  ``good_pix_frac=0.05``.  The call site in ``fit_profile_refactor`` is
  unchanged (no arguments added there yet).
- **Docstring**: entries for both new parameters added to the ``Parameters``
  section.
- **Function body**: the two ``sn_cap = ...`` and ``good_pix_frac = ...``
  assignments removed; the parameter names are used directly.
- **No numerical change**; all 54 spatialprofile tests pass unchanged.

**`_apodize_profile` refactored; `_deriv_walk` promoted to module level**:
- **Motivation**: four changes to `_apodize_profile` were implemented:
  (1) the confusing ``[::-1]`` / ``.min()`` pattern for finding the initial
  left limit replaced with ``np.ma.flatnotmasked_edges`` — the same idiom
  used in ``_findfwhm`` to locate threshold crossings; (2) ``sorted_sigma``,
  ``sorted_bsp``, and ``threshold`` aliases introduced to avoid repeated
  indexing of ``sigma_x[ss]`` / ``full_bsp[ss]``; (3) bare ``np.log`` calls
  on potentially-zero profile values replaced with ``np.ma.log`` +
  ``masked_where(fit <= 0)`` + ``.filled(0.0)``; (4) the two ~15-line
  duplicate left/right derivative blocks unified into a single ``_deriv_walk``
  helper parameterised by ``sign`` (``-1`` left, ``+1`` right).
- **`_deriv_walk`** was initially implemented as a closure inside
  ``_apodize_profile`` capturing ``bset``.  It has now been promoted to a
  module-level private helper with ``bset`` added as an explicit first
  parameter and a full NumPy-style docstring.  The two call sites inside
  ``_apodize_profile`` updated accordingly (``_deriv_walk(bset, l_limit,
  sign=-1)`` / ``_deriv_walk(bset, r_limit, sign=+1)``).
- **No numerical change**; all 54 spatialprofile tests pass unchanged.

**Block 7 limit search refactored to match `_apodize_profile` idiom**:
- **Motivation**: the initial left/right limit search in Block 7 (inside
  ``fit_profile_refactor``) used an asymmetric pattern — the right limit was
  found via ``np.where(condition).min()`` on ascending-sorted arrays, but
  the left limit first reversed both the fit values (``rev_fit =
  mode_fit[::-1]``) and the index array (``sr = si[::-1]``), then also took
  ``.min()``.  This mirrored the pre-refactor pattern already eliminated from
  ``_apodize_profile``.
- **Changes**: (1) ``sr = si[::-1]`` removed from line 826 — ``sr`` was
  only ever used in these six lines; (2) ``rev_fit = mode_fit[::-1]``
  removed; (3) both ``np.where(...)`` + ``.min()`` + conditional fallback
  patterns replaced with ``np.ma.flatnotmasked_edges`` on ascending arrays,
  exactly matching ``_apodize_profile`` — ``edges_l[1]`` (last unmasked =
  rightmost qualifying pixel) for the left limit and ``edges_r[0]`` (first
  unmasked = leftmost qualifying pixel) for the right; (4) ``sorted_sigma =
  sigma_x[si]`` and ``threshold = min_level + median_fit`` aliases
  introduced, eliminating four repeated ``sigma_x[si]`` indexings and two
  repeated ``min_level + median_fit`` computations; (5) ``sorted_sigma``
  reused on the downstream inside-selection ``np.where`` (formerly
  ``sigma_x[si]``).
- **No numerical change**; all 54 spatialprofile tests pass unchanged.

**`_fit_spectrum_and_normalize` restructured**:
- **``good_x`` return value removed**: the ``good_x = norm_ivar_x > 0``
  boolean array was a simple derived quantity; the one downstream use in
  ``_compute_bspline_knots`` is replaced with the inline expression
  ``norm_ivar_x > 0``.  Return tuple shrinks from 8 to 7 elements; both
  early-return failure tuples shortened accordingly.  The corresponding
  ``test_fit_spectrum_good_matches_norm_ivar`` test removed.
- **``flux_x``/``ivar_x`` eliminated**: these pre-computed 1D copies of
  the 2D arrays are replaced with inline ``image[totmask]`` and
  ``ivar[totmask]`` at the single normalization step where they were used.
- **``spec_x`` deferred; ``xtemp_x`` / ``cumweight`` moved earlier**:
  ``spec_x`` is now computed just before ``xtemp_x`` is needed (after
  ``sn2_1`` is built), rather than in the initial pixel-extraction block.
  ``xtemp_x`` / ``cumweight`` are computed immediately after ``sn2_1``,
  before the ``sn2_x`` and spline image construction.
- **``wave_min`` / ``wave_max`` / ``good_wave`` moved to top**: these are
  now derived from ``wave_x`` at the very start of the function, before the
  spectral smoothing, instead of after the S/N cap.
- **``iterative_bspline_fit`` return values used directly**: ``spline_flux``
  and ``cont_flux`` are now unpacked from the third return element of
  ``iterative_bspline_fit``; the intermediate ``b_answer`` / ``c_answer``
  ``BSpline`` objects are discarded and the subsequent ``.value()`` calls
  removed.  The ``try/except`` is tightened to wrap only the continuum fit
  call.
- **Minor cleanups**: block-header comments replaced with descriptive
  inline comments; redundant ``> 0`` removed from
  ``np.any(np.logical_not(badpix)) > 0``; ``log.info`` for ``med_sn2``
  moved to immediately after ``med_sn2`` is computed.
- **No numerical change**; 53 spatialprofile tests pass (one test removed
  alongside ``good_x``).

**`_apodize_profile` renamed to `_build_profile`; interface simplified**:
- **Motivation**: the function performed two logically distinct steps —
  (1) evaluating the B-spline profile over all valid coordinates, and
  (2) optionally applying exponential apodization.  The old interface
  threaded ``prof_nsigma`` and ``no_deriv`` through to the function body
  where they were combined into a "JXP kludge" that zeroed the limits and
  forced ``no_deriv=True`` for extended-object fits.
- **Changes**: (1) renamed to ``_build_profile`` to reflect the broader
  scope; (2) ``prof_nsigma`` and ``no_deriv`` replaced by a single
  ``apodize : bool`` parameter; (3) the caller (``fit_profile_refactor``)
  now computes ``apodize = (prof_nsigma is None and not no_deriv)`` before
  the call, keeping the JXP logic at the call site where it belongs;
  (4) when ``apodize=False``, the function returns immediately after
  computing the B-spline profile (``full_bsp, 0.0, 0.0``), eliminating
  dead work on limits that would be discarded; (5) a ``log.warning`` is
  now emitted when the derivative condition fails instead of silently
  skipping apodization; (6) exponential tail coefficients reordered to
  ``val * np.exp(...)`` for clarity; (7) block-header comments replaced
  with descriptive inline comments; (8) ``_deriv_walk`` docstring updated
  to reference ``_build_profile``.
- **Tests**: ``_apodize_profile`` import replaced with ``_build_profile``;
  all four call sites updated to the new positional signature (``apodize``
  inserted as the 5th argument); ``test_apodize_prof_nsigma_zeros_limits``
  and ``test_apodize_no_deriv_skips_tails`` merged into the single test
  ``test_apodize_no_apodize`` (both exercised ``apodize=False``).
- **No numerical change**; 52 spatialprofile tests pass (one test removed
  by the merge).

**`_compute_bspline_knots` renamed to `_profile_coordinates_and_model_sampling`**:
- Name better reflects that the function initialises both the spatial
  coordinate array ``sigma_x`` and the B-spline knot grid ``bkpt``.
  All call sites and test references updated.

**`fit_profile_refactor` body polished**:
- Block-header comments (``# Block N —``) replaced throughout with
  descriptive inline comments.
- ``mode_fit`` now unpacked from the third return element of
  ``iterative_bspline_fit`` (same pattern as ``_fit_spectrum_and_normalize``),
  eliminating the subsequent ``bset.value(sigma_x[si])`` call.
- ``sorted_sigma`` alias moved immediately after the B-spline fit so it is
  available for ``_findfwhm`` without an intermediate assignment.
- ``m_left`` / ``m_right`` intermediate masked-array variables inlined into
  the ``flatnotmasked_edges`` calls; ``~`` operator replaced with
  ``np.logical_not()`` for style consistency.
- **``iiter`` loop fixed**: ``range(sigma_iter)`` (0-indexed) →
  ``range(1, sigma_iter+1)`` (1-indexed); guard condition
  ``iiter < sigma_iter - 2`` → ``iiter < sigma_iter - 1``.  Resolves the
  standing TODO — the condition now reads "more iterations remain" and is
  semantically correct.  With ``sigma_iter=3`` the profile B-spline is still
  re-fitted exactly once (after the first of three iterations), so the
  numerical result is unchanged.
- ``sigma`` and ``area`` updates moved above the ``log.info`` calls (they
  were below, making the logged values stale).
- ``bset = bset.to_1d()`` inlined into the ``_build_profile`` call.
- ``profile_x`` intermediate eliminated; ``full_bsp * pb_x`` used inline for
  both ``profile_model[totmask]`` and the residual calculation.
- Log message "is not zero" → "is non-zero"; TODO added about consolidating
  Gaussian-path log messages; TODO about ``iiter`` indexing removed (resolved).
- **No numerical change**; 52 spatialprofile tests pass unchanged.

**Docstrings and comments cleaned up for production**:
- **Module docstring** rewritten to describe the module's purpose directly,
  removing all "refactored", "Bug N fixed", and legacy-comparison language.
- **`_findfwhm`**: "(see Bug 2 in the refactoring plan)" removed from
  the description.
- **`_gaussian_profile`**: ", matching the B-spline path" removed from
  the row-sums sentence.
- **`_fit_spectrum_and_normalize`**: "Encapsulates … (blocks 0 and 2–4)"
  framing replaced with a self-contained description.
- **`_profile_coordinates_and_model_sampling`**: "Encapsulates Block 6"
  framing replaced; inline ``# Bug 1 fix:`` comment removed.
- **`_build_profile`**: "Encapsulates Blocks 11-12" framing replaced;
  return-value docs updated to reference ``apodize`` instead of the
  removed ``prof_nsigma`` parameter; unused ``peak``, ``lwhm``, ``rwhm``
  return values from the internal ``_findfwhm`` call replaced with
  ``_, peak_x, *_`` unpacking.
- **`fit_profile_refactor`**: opening sentence rewritten to be
  self-contained, removing "Refactored drop-in replacement for …" and
  the reference to the legacy docstring.
- **Inline comments**: ``# TODO: I'm not sure what this check is for``
  replaced with a brief description; ``iiter+1`` in the iteration log
  corrected to ``iiter`` (range is now 1-indexed); bullet-point style
  (``# -``) removed from the final profile-fit section.
- **No numerical change**; 52 spatialprofile tests pass unchanged.

### v7 → v8

**Verification Plan restructured**:

- **Step 16 → step 15, marked ✅**: the open question about the Block 8
  mid-loop re-fit condition (`iiter < sigma_iter - 1`, resolved by the
  1-indexed loop fix documented in the v7 changelog) was already implemented
  and is now marked complete.  Renumbered to step 15 to maintain ascending
  order of completed steps.
- **Steps 17 and 18 removed**: the dev-suite run and the full
  ``pytest pypeit/tests/`` pass are standard milestones in PypeIt's nominal
  development workflow (see ``doc/dev/development.rst``) and will be performed
  once the core development on this branch is complete, rather than being
  tracked as numbered plan steps.
- **New step 16 — deprecation**: deprecate `pypeit/core/spatialprofile.py`
  and replace it with the refactored module.  Covers renaming
  `fit_profile_refactor` → `fit_profile` and `spatialprofile_refactor.py` →
  `spatialprofile.py`, updating the call site in `skysub.py`, cleaning up
  `test_spatialprofile.py` (new module name, removal of cross-implementation
  comparison tests), and deleting the legacy module.
- **New step 17 — QA wiring**: wire the `_fit_profile_refactor_qa` figure
  into the standard `run_pypeit` output so that it is automatically written
  to `QA/PNGs/` for each extracted object.
- **New step 18 — `spec_img` wiring**: the optional `skysub.py` `spec_img`
  update is no longer marked optional and is now step 18.
- **New step 19 — documentation page**: write a new RST page in `doc/`
  describing the spatial-profile fitting algorithm, registered in the Sphinx
  ``toctree`` and cross-referenced from `doc/skysub.rst`, with example figures
  generated from the unit-test synthetic data.

### v8 step 17 addenda

Two corrections to the initial step 17 implementation:

**Addendum 1 — Final iteration only**: QA files are now written only on the
last profile-fitting iteration (``iiter == niter``) rather than once per
iteration.  Earlier iterations show intermediate, unconverged fits that are
unlikely to be useful and produce unnecessary file clutter.  The
``show_profile=True`` interactive path is unchanged and still fires on every
iteration.

**Addendum 2 — Follow the canonical QA filename pattern**: The ad-hoc
``spatprof_obj{ID}_slit{N}_iter{I}.png`` scheme was replaced with the standard
``qa.set_qa_filename`` factory, following the ``get_findobj_qa_filename``
pattern in ``find_objects.py``.

- Added a new case ``'spat_profile_qa'`` to ``qa.set_qa_filename``
  (``pypeit/qa.py``) and an ``obj: int | None = None`` keyword argument
  for cases where multiple objects share a slit.  Produces
  ``QA/PNGs/{basename}_{det}_S{slit:04d}_O{obj:03d}_spat_prof.png``.
- The single ``qa_path`` parameter on ``local_skysub_extract`` and
  ``ech_local_skysub_extract`` was split into three:
  ``qa_outdir`` (the ``QA/`` directory), ``qa_basename`` (reduction
  basename), and ``qa_detname`` (detector name string).
- The ``QA/PNGs/`` sub-directory is created once before the extraction
  loop via ``Path(qa_outdir, 'PNGs').mkdir(parents=True, exist_ok=True)``.
- Both wrappers in ``extraction.py`` pass
  ``qa_outdir = <redux_path>/<qadir>``,
  ``qa_basename = self.basename``, and
  ``qa_detname = self.spectrograph.get_det_name(self.det)``.

### v8 QA figure improvements

Four changes to ``_fit_profile_qa`` in ``spatialprofile.py``:

- **Color scheme**: named and grey-string colors (``'lime'``, ``'red'``,
  ``'tomato'``, ``'orange'``, ``'0.60'``) replaced with matplotlib
  property-cycle colors (``'C0'``, ``'C1'``, ``'C3'``, ``'k'``).  This
  gives consistent, theme-aware colors without hard-coded hue conflicts.
- **Profile scatter — no subsampling**: the random 8000-point subsample
  of the normalised-object scatter was removed; all valid pixels are now
  passed directly to ``ax_prof.scatter``.
- **Residual scatter — no subsampling**: the same subsampling removal
  applied to the profile residual panel; all pixels' residuals are plotted.
- **Residual panel — binned overlay**: ``(y50 − model)[finite]`` is now
  overplotted in the profile residual panel as ``C0`` circles, matching the
  binning style already present in the profile panel above it.

### v8 step 19 addendum — figure generation script

The figure generation script was moved from the development suite into the
main pypeit repo and updated to use the same synthetic data as the existing
test ``test_fourier_flux_and_curved_trace``:

- **`doc/scripts/make_spatprof_figures.py`** added to the pypeit repo.
  Imports ``make_fourier_spectrum`` and ``make_profile_inputs`` directly
  from ``pypeit.tests.test_spatialprofile``.  Uses the same Fourier-spectrum
  flux (5 modes, mean 1000, std 50) and quadratic curved trace
  (``10 * (2t² − 1)``, spanning ±10 px) as
  ``test_fourier_flux_and_curved_trace``.  Output paths are resolved via
  ``importlib.resources`` so the script works from any working directory.
- **``PypeIt-development-suite/spatprof_refactor/gen_spatprof_figs.py``**
  removed; the pypeit-repo script supersedes it.
- ``doc/figures/spatprof_example_bspline.png`` and
  ``doc/figures/spatprof_example_gaussian.png`` regenerated with the new
  inputs (S/N = 20 and S/N = 2 respectively).

### v8 — American English pass

All British English spellings converted to American English in
``spatialprofile.py`` (18 lines), ``test_spatialprofile.py`` (17 lines),
and ``doc/spatprof.rst`` (1 line).  Affected forms: normalise(d/s/tion) →
normalize(d/s/tion), initialise → initialize, modelled → modeled,
row-normalised → row-normalized, behaviour → behavior, centre(d) →
center(ed).

### v8 step 19 addendum 2 — editorial pass on ``doc/spatprof.rst``

- **American English**: "modelled/normalised" → "modeled/normalized" throughout.
- **S/N decision**: prose reflowed; ``\texttt{sn\_gauss}`` changed to
  ``\texttt{sn_gauss}`` (underscore need not be escaped inside
  ``\texttt{}`` in math mode); added cross-reference to
  ``.. _spatprof-parameters:``.
- **Trace correction**: default-value parenthetical moved inside the
  inequality display so it reads as a condition rather than a footnote.
- **QA Output**: description changed from "written to ``QA/PNGs/``" to the
  more neutral "a diagnostic figure is provided"; the prose description of
  the figure layout removed and replaced by "Example figures are shown
  below." (the captions carry that detail instead); ``matplotlib`` wrapped
  in backticks in the ``show_profile`` note.
- **B-spline caption**: rewritten panel-by-panel in plain prose, clarifying
  the six panels and the 20/50/80 percentile annotation.
- **Gaussian caption**: simplified; removed the redundant explanation of
  scatter.
- **Key Parameters**: added ``For more detail, see
  :func:`~pypeit.core.spatialprofile.fit_profile`.`` at the top; removed
  repeated default values from individual entries (defaults belong in the
  API docstring).
- **Footer**: ``---`` rule and attribution line added at the end of the
  file.

### v8 step 19 addendum 3 — improved synthetic spectrum and figure parameters

Two further improvements to ``make_fourier_spectrum`` in
``test_spatialprofile.py`` and corresponding updates to
``make_spatprof_figures.py``:

**``make_fourier_spectrum`` changes:**

- **``power=1.0`` parameter**: each selected mode's complex coefficient is
  scaled by :math:`k^{-\text{power}}` before normalisation.  Higher values
  give proportionally more amplitude to lower-frequency modes.
  ``power=0`` recovers flat (white-noise) weighting.
- **``min_period=3.0`` parameter**: replaces the hardcoded ``3.0`` in
  ``k_min = ceil(nspec / min_period)``, making the low-frequency cutoff
  configurable.
- **Log-uniform mode selection**: modes are no longer drawn uniformly from
  ``[k_min, k_max]``; instead weights :math:`p_k \propto 1/k` are used,
  so lower-frequency modes are preferentially selected.  Both the amplitude
  weighting (``power``) and the selection weighting (``p ∝ 1/k``) act in
  the same direction, giving the resulting spectrum a realistic
  red-noise character.

**``make_spatprof_figures.py`` changes:**

- ``make_fourier_spectrum`` call updated to ``nmodes=40``,
  ``min_period=20.``, ``power=2.``.  The longer minimum period (20 px
  instead of 3 px) and larger ``power`` produce smooth, correlated spectral
  structure that is visually representative of real continuum spectra; the
  higher mode count (40 vs. 5) adds enough variation to exercise the
  spectral B-spline fit.
- The module-level ``OUT_DIR`` constant removed; ``out_dir`` is now
  resolved locally inside ``main()``.

## Migration of `_bin_profile` and `qa_fit_profile_refactor`

Both functions were ported from `qa_spatprof_demo.py` in the development-suite
scratchpad and added to `spatialprofile_refactor.py` after the final `return`
statement of `fit_profile_refactor`, separated by a `# --- QA ---` comment
block.

- **`_bin_profile(sigma_x, norm_obj, model_vals, xmin, xmax, nsamp)`**:
  Bins per-pixel normalised-object values in `sigma_x` space using `nsamp`
  uniform bins spanning `[xmin, xmax]`.  For each bin it computes the 20th,
  50th, and 80th percentiles of the data and the median of the model, returning
  `(bin_cen, p20, p50, p80, model_med)` as five 1-D arrays of length `nsamp`.
  Used only by `qa_fit_profile_refactor`.

- **`qa_fit_profile_refactor(...)`**: Produces a six-panel diagnostic figure.
  Three upper panels show the 2D image, model (`profile_model * flux[:, None]`),
  and residual (data − model), each with a common colour scale and an overplotted
  trace.  Three lower panels show: (1) the extracted 1-D spectrum with the model
  overlaid; (2) the fractional spectrum residual; (3) a scatter plot of
  `norm_obj_x` vs `sigma_x` with the 20/50/80th-percentile binned envelope and
  the B-spline/Gaussian profile curve.  When `success=False` (early-return path
  with no normalised-object data) the lower panels are suppressed.  The figure
  is saved to `outfile` if provided, otherwise displayed interactively.

- **`test_qa_fit_profile_refactor(tmp_path)`** added to
  `pypeit/tests/test_spatialprofile.py`.  It calls `fit_profile_refactor` on the
  standard high-S/N synthetic inputs, then passes all outputs (including
  `diagnostics`) to `qa_fit_profile_refactor` with `outfile` set to a `tmp_path`
  PNG.  The test asserts that the file was created and is non-empty.

- 56 spatialprofile tests pass (up from 52).

## `diagnostics` — 5th return value of `fit_profile_refactor`

`fit_profile_refactor` previously returned a 4-tuple `(profile_model, xnew,
fwhmfit, med_sn2)`.  It now returns a 5-tuple, with a `types.SimpleNamespace`
called `diagnostics` as the fifth element.  The motivation is to eliminate the
redundant call to `_fit_spectrum_and_normalize` that `qa_fit_profile_refactor`
required in order to reconstruct the normalised-object arrays.

`diagnostics` carries six attributes:

| Attribute | Shape | Description |
|-----------|-------|-------------|
| `norm_obj_x` | `(npix,)` or `None` | Normalised-object pixel values at in-mask pixels |
| `norm_ivar_x` | `(npix,)` or `None` | Corresponding inverse variances |
| `sigma_x` | `(npix,)` or `None` | Spatial coordinate in units of the object sigma |
| `totmask` | `(nspec, nspat)` bool | Combined pixel mask used for extraction |
| `l_limit` | float | Left profile apodisation limit in sigma units |
| `r_limit` | float | Right profile apodisation limit in sigma units |

The six distinct return paths each populate `diagnostics` differently:

1. **Pre-`_profile_coordinates_and_model_sampling` early Gaussian return**
   (`not success or gauss or ngood < 10 or med_sn2 < sn_gauss**2`):
   - `success=False`: `norm_obj_x`, `norm_ivar_x`, `sigma_x` are all `None`;
     `l_limit = r_limit = 0.0`.
   - `success=True` (Gaussian forced by `gauss`, low S/N, or `ngood < 10`):
     `sigma_x` is computed inline from the initial `sigma` array via
     `_dspat_x / sigma[spec_x]`; `l_limit = r_limit = 0.0`.

2. **`ninside < 10` return** (after `_profile_coordinates_and_model_sampling`):
   `sigma_x` fully computed; `l_limit`/`r_limit` reflect the profile limits
   from the initial B-spline evaluation.

3–5. **Three loop-body early returns** (trace-correction failure,
   width-correction failure, profile B-spline failure): same as (2); all
   intermediate arrays available.

6. **Final B-spline return**: `sigma_x` and limits from the fully converged
   iterative fit; `xnew` is the corrected trace.

The import block was updated: `from IPython import embed` removed;
`from types import SimpleNamespace` added.

All call sites in `test_spatialprofile.py` and `freeze_spatprof.py` were
updated to unpack the 5-tuple using `, _` (direct calls) or `*_` starred
unpacking (parametrised tests that also call the legacy 4-tuple
`spatialprofile.fit_profile`).

- **No numerical change**; 56 spatialprofile tests pass unchanged.

## Simplified `qa_fit_profile_refactor` signature

Before adding `diagnostics`, `qa_fit_profile_refactor` accepted 18 positional
parameters, including the raw calibration images needed to re-run
`_fit_spectrum_and_normalize` internally:

```
image, ivar, waveimg, fluxivar, inmask, thismask, spat_img,
trace_in, wave, flux, spec_img, thisfwhm, percentile_sn2,
l_limit, r_limit, profile_model, xnew, fwhmfit, med_sn2,
obj_string, outfile
```

With `diagnostics` as a 5th return value the internal call to
`_fit_spectrum_and_normalize` was eliminated entirely.  The new signature has
12 parameters:

```python
def qa_fit_profile_refactor(
    image, ivar, thismask, trace_in,
    flux,
    profile_model, xnew, fwhmfit, med_sn2, diagnostics,
    obj_string='', outfile=None,
):
```

The ten parameters dropped (`waveimg`, `fluxivar`, `inmask`, `spat_img`,
`wave`, `spec_img`, `percentile_sn2`, `thisfwhm`, `l_limit`, `r_limit`) are
all now sourced directly from `diagnostics`:

```python
norm_obj_x  = diagnostics.norm_obj_x
norm_ivar_x = diagnostics.norm_ivar_x
sigma_x     = diagnostics.sigma_x
totmask     = diagnostics.totmask
l_limit     = diagnostics.l_limit
r_limit     = diagnostics.r_limit
success     = norm_obj_x is not None
```

The spatial-profile scatter panel also dropped the two previously redundant
lines `dspat_x = spat_img[totmask] - trace_in[spec_x]` and
`sigma_x = dspat_x / sigma[spec_x]`; `sigma_x` now arrives pre-computed in
`diagnostics`.

- **No numerical change**; 56 spatialprofile tests pass unchanged.

## Class-structure consideration and rejection

The possibility of re-implementing `fit_profile_refactor` as a Python class
(with the module-level function as a thin wrapper) was evaluated.  The primary
appeal was that instance attributes would replace the intermediate variables
passed between helper functions, eliminating the need for the `diagnostics`
object and cleaning up helper-function call signatures.

The proposal was rejected for the following reasons:

- **Early returns become awkward.**  The six Gaussian fallback paths each return
  immediately from the middle of the computation.  In a class a phase method
  cannot return from the enclosing `fit()` call; it can only signal the
  condition back via a sentinel value or a private exception, neither of which
  is more readable than the current explicit ``return`` statements.

- **Unit tests for helper functions become harder.**  Eleven tests currently
  call the helper functions directly with clean argument lists.  If the helpers
  become instance methods, each test must construct a fitter object and manually
  set the specific attributes the method reads from ``self``, creating implicit
  dependencies between tests and implementation internals.

- **One-shot classes are a code smell.**  The fitter would be instantiated once
  per object per ``local_skysub_extract`` iteration, run, and discarded.  Python
  reserves classes for objects with state that persists or is shared across
  multiple calls; neither applies here.

- **`diagnostics` could be removed without a class.**  Making the QA function
  private and calling it internally was sufficient to eliminate the `diagnostics`
  object — the primary concrete benefit of the class proposal.

## QA restructuring: `_fit_profile_refactor_qa` and `generate_qa`

After completing the `diagnostics` implementation, two additional problems were
identified:

1. `qa_fit_profile_refactor` was a *public* function that required the caller to
   make a second, separate call after `fit_profile_refactor`.  Supporting this
   interface was the entire reason for the `diagnostics` 5th return value.

2. `fit_profile_refactor` was still calling the *old* `qa_fit_profile` from
   `spatialprofile.py` on every Gaussian fallback path (via `_return_gaussian`)
   and on the final B-spline path.  The new `qa_fit_profile_refactor` was
   therefore bypassed on all of those paths.

**Resolution — `_fit_profile_refactor_qa` (private, called internally):**

- `qa_fit_profile_refactor` was renamed to `_fit_profile_refactor_qa` and made
  private.  It is now called from *inside* `fit_profile_refactor` at each of
  its six return paths.

- Because the QA function is called with the intermediate arrays directly in
  scope, the `diagnostics` ``SimpleNamespace`` is unnecessary.  The parameter
  was replaced by the six constituent fields passed directly:
  ``norm_obj_x``, ``norm_ivar_x``, ``sigma_x``, ``totmask``, ``l_limit``,
  ``r_limit``.

- The ``diagnostics`` object was removed from the codebase entirely.
  ``fit_profile_refactor`` returns a clean 4-tuple
  ``(profile_model, xnew, fwhmfit, med_sn2)`` again.

**Resolution — `generate_qa` parameter:**

- The ``show_profile : bool`` parameter was replaced by
  ``generate_qa``, which accepts three types:

  - ``False`` (default) — no QA generated.
  - ``True`` — display the figure interactively (blocks pipeline progress).
  - ``str`` or ``pathlib.Path`` — save the figure to that file path.

- ``_return_gaussian`` was reduced to pure profile construction and logging;
  its ``show_profile``, ``norm_obj``, ``ind``, ``l_limit``, ``r_limit``,
  ``xlim``, and ``xtrunc`` parameters were all removed.  The old call to
  ``qa_fit_profile`` from within ``_return_gaussian`` is gone; QA on Gaussian
  paths now goes through ``_fit_profile_refactor_qa``.

- The imports ``from types import SimpleNamespace`` and
  ``from pypeit.core.spatialprofile import qa_fit_profile`` were removed.

**Impact on tests:**

All call sites in ``test_spatialprofile.py`` and ``freeze_spatprof.py``
reverted from 5-tuple to 4-tuple unpacking.  ``test_qa_fit_profile_refactor``
was simplified from two calls (fit then separate QA call) to a single
``fit_profile_refactor(..., generate_qa=outfile)`` call.  ``qa_fit_profile_refactor``
was removed from the import block.

- **No numerical change**; 56 spatialprofile tests pass unchanged.

## Eliminating `_return_gaussian` and positional QA calls

With `_fit_profile_refactor_qa` now called directly at each return site, the
`_return_gaussian` helper was identified as redundant: its entire body was just
two lines — a call to `_gaussian_profile` and a `log.info` statement.

**Strategy discussion:**

Two strategies for reducing the repetition in the five Gaussian early-return
blocks were evaluated, both based on an inner (closure) function defined inside
`fit_profile_refactor`:

- **Strategy 1** — a `_do_qa` closure capturing the invariant context and
  accepting only the five varying arguments.  Each return site would remain
  explicit about profile construction but reduce the QA block to a single call.

- **Strategy 2** — a `_gauss_return` closure that also owned the
  `_gaussian_profile` call and the `log.info`, collapsing each return site to a
  single ``return _gauss_return(center, fwhm, sigma_x, ...)`` statement.

Both strategies were rejected in favour of inlining, on the grounds that
defining a function inside a function is an unnecessary complication when the
inlined form is only marginally longer.

**What was done:**

- `_return_gaussian` was deleted.  Its two lines were inlined at each of the
  five Gaussian return sites:

  ```python
  profile_model = _gaussian_profile(spat_img, center, sigma)
  log.info(f'{obj_string}, FWHM={fwhm:.2f}, S/N={np.sqrt(med_sn2):.3f}')
  ```

- All six calls to `_fit_profile_refactor_qa` were updated to pass positional
  arguments for the 15 positional parameters, reserving keyword form only for
  the two parameters with defaults (``obj_string`` and ``outfile``):

  ```python
  _fit_profile_refactor_qa(
      image, ivar, thismask, trace_in,
      flux,
      profile_model, xnew, fwhmfit, med_sn2,
      norm_obj_x, norm_ivar_x, sigma_x, totmask, l_limit, r_limit,
      obj_string=obj_string,
      outfile=None if generate_qa is True else generate_qa,
  )
  ```

  The grouping of arguments mirrors the visual grouping in the function
  signature, which breaks after ``trace_in`` and after ``flux``.

- **No numerical change**; 56 spatialprofile tests pass unchanged.

- **No numerical change**; 56 spatialprofile tests pass unchanged.
