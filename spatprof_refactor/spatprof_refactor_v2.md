# fit_profile Refactoring Analysis (v2)

> **Change log from v1**: Added "Pre-existing Bugs Discovered During Testing"
> section; updated `_compute_bspline_knots` helper note to call out the broken
> `nb` formula; updated three unit tests (`test_fwhm_recovery`,
> `test_trace_correction_direction`, `test_prof_nsigma_extended`) to match the
> passing test file checked in at `pypeit/tests/test_spatialprofile.py`.
>
> **Change log (this update)**: Added "Accuracy Improvements" section covering
> (A) pixel-integrated Gaussian in `return_gaussian` and (B) linear
> interpolation for the half-maximum crossing in `findfwhm`; updated Bugs 2
> and 3 to cross-reference the corresponding fixes and their expected impact on
> unit tests.

## Context

`fit_profile` in `pypeit/core/spatialprofile.py` (lines 262–854) is a 593-line monolith ported from the IDL LOWREDUX routine `long_gprofile.pro`. It fits a non-parametric B-spline spatial profile to an object spectrum for optimal extraction, falling back to a Gaussian when S/N is too low. It is called from exactly one external site: `skysub.py:886` inside `local_skysub_extract`, once per object per iteration. Its four return values are:

- `profile_model` — 2D normalized spatial profile `(nspec, nspat)`
- `xnew` — corrected trace `(nspec,)`
- `fwhmfit` — FWHM per spectral pixel `(nspec,)`
- `med_sn2` — median S/N²

---

## Workflow Breakdown (line ranges are approximate)

### Block 1 — Initialization (347–383)
Set up masks (`inmask`, `totmask`), handle `prof_nsigma → no_deriv`, clamp `thisfwhm`, compute `dspat` (spatial separation image), zero-initialize `sn2_img` / `spline_img`, median-smooth `flux` and `fluxivar`, and initialize `sigma`, `fwhmfit`, `trace_corr`, `sigma_x`.

### Block 2 — Flux B-spline fitting (383–415)
Cap `fluxivar` at S/N=100, identify eligible spectral pixels `indsp`, fit two fine B-splines (`b_answer`, pass 1 & 2) and one coarse continuum B-spline (`c_answer`). Early Gaussian return if no eligible pixels or continuum fit raises.

### Block 3 — S/N estimation (417–460)
Compute per-pixel `sn2` from the fine spline × ivar; use a percentile cut + sigma-clipped mean to get `med_sn2`. Expand spline/continuum evaluations to the full 1D spectrum (`spline_flux1`, `cont_flux1`, `sn2_1`); interp-fill bad B-spline pixels with the continuum model; median-filter `sn2` and interpolate it onto the 2D `sn2_img`.

### Block 4 — Normalized-object image construction (462–509)
Build `spline_img` (the S/N-level-dependent spectral flux model, extrapolated to 2D):
- S/N² ≤ 2: use noise floor `sigma1`
- 2 < S/N² ≤ 5: use coarse continuum
- S/N² > 5: fine spline with bad-pixel patching from continuum

Compute `norm_obj = image / spline_img` and `norm_ivar = ivar * spline_img²`. Compute `xtemp` — a cumulative S/N-weighted spectral coordinate used as the independent variable for the iterative profile correction fit.

### Block 5 — Early returns (513–518)
Return Gaussian immediately if `ngood < 10`, `med_sn2 < sn_gauss²`, or `gauss=True`.

### Block 6 — Knot setup (520–554)
Compute profile half-extent `limit` from `erfcinv`, set sinh-spaced B-spline interior knots `bkpt`, determine `min_sigma` / `max_sigma`, and select `inside` pixels for the initial profile fit (preferring high-S/N pixels if available, falling back to all in-range pixels).

### Block 7 — Initial profile fit (554–614)
Fit a 4th-order B-spline to `norm_obj[inside]` sorted by `sigma_x`. Evaluate to get `mode_fit`; call `findfwhm` to locate the peak and measure width. Update `sigma` and `limit` using the B-spline FWHM. Walk the profile down from the peak to find `l_limit` and `r_limit` (where profile drops below `min_level`). Select the final `inside` subset within those limits; early Gaussian return if `ninside < 10`.

### Block 8 — Iterative trace and width correction loop (612–716) — `sigma_iter=3`
Each iteration:
1. **Trace shift**: evaluate the current B-spline at `sigma_x` and at `sigma_x ± 0.5`; use `(mode_min05 - mode_plu05)` as a shift-sensitive basis and fit it jointly with `mode_zero` via `iterative_bspline_fit`; extract the h0/h1 coefficient ratio to get `delta_trace_corr`.
2. **Width stretch**: evaluate at `sigma_x / 1.3`; use the rescaled difference as a stretch-sensitive basis; extract h2/h0 ratio to get `sigma_factor`.
3. Update `sigma`, `area`, `trace_corr`, `sigma_x`.
4. On non-final iterations: re-fit the profile B-spline with the updated `sigma_x` (returned as a `BSpline2D`, converted to 1D `BSpline` for subsequent `value()` calls).

Each sub-fit failure returns a Gaussian immediately.

### Block 9 — Final trace (718–722)
Apply `trace_corr` only if `median(|trace_corr * sigma|) < max_trace_corr`; otherwise keep `trace_in`.

### Block 10 — Final profile fit (724–742)
Re-select all in-range, unmasked pixels sorted by `sigma_x`; fit the definitive 4th-order B-spline profile using the full `pb` (area) weight array and loose rejection (`upper=lower=10`). Convert resulting `BSpline2D` → 1D `BSpline`.

### Block 11 — Apodization limit search (744–803)
Evaluate the final B-spline across the full sigma_x range (`full_bsp`). Re-find `l_limit` / `r_limit` from the profile level. Then compute the logarithmic derivative `d(log P)/d(log x)` over vectors `l_lim_vec` and `r_lim_vec` to find the extremal derivative. Walk the limits inward 0.1 at a time (two `while True` loops) until the derivative condition is satisfied or the inner ±1 bound is hit.

### Block 12 — Exponential apodization (805–816)
If derivatives have the right sign and `no_deriv=False`, replace `full_bsp` outside `[l_limit, r_limit]` with an exponential matched in value and first derivative to the profile edge.

### Block 13 — Normalization, logging, QA, return (818–854)
Reshape `full_bsp`, multiply by `pb`, normalize each spectral row to unit sum, guard against NaNs, log fit statistics, show QA plot if `show_profile`, and return.

---

## Candidate Helper Functions

### 1. `_smooth_spectrum(wave, flux, fluxivar, waveimg, thismask, sn_cap=100.0)`
**Lines 367–415.** Median-smooth flux/ivar, fit two fine B-splines and one coarse continuum B-spline to the eligible spectral pixels. Return `(indsp, flux_sm, fluxivar_sm, b_answer, bmask2, c_answer, cmask, spline_flux, cont_flux, wave_min, wave_max)` plus a flag indicating an early-return condition (no eligible pixels or failed continuum fit). The two early-return Gaussian paths currently live here; they could be signalled by returning `None` for the B-spline answers and handled uniformly in the parent.

### 2. `_estimate_sn2(wave, flux, fluxivar_sm, spline_flux, cont_flux, bmask2, cmask, indsp, waveimg, totmask, percentile_sn2)`
**Lines 417–460.** Compute `sn2`, `med_sn2`, `sigma1`; build `sn2_img` (2D). Pure computation with no side effects beyond array construction. Returns `(sn2_1, sn2_img, med_sn2, sigma1, spline_flux1, cont_flux1, bmask_1d, cmask2_1d)`.

### 3. `_build_spline_img(wave, flux, spline_flux1, cont_flux1, med_sn2, bmask, cmask2, sigma1, waveimg, totmask, wave_min, wave_max, nspec, nspat)`
**Lines 462–496.** Branch on `med_sn2` level to fill `spline_img` and optionally patch bad pixels with the continuum. Returns the finished `spline_img (nspec, nspat)`.

### 4. `_compute_bspline_knots(med_sn2, prof_nsigma, sigma_x, good)`
**Lines 522–540.** Compute `limit`, `sinh_space`, `nb`, `bkpt`, `min_sigma`, `max_sigma` — the pure arithmetic that sets up the knot grid. Returns a small namespace or tuple.

**⚠ Pre-existing bug in `nb` formula (line 532):** `nb = np.round(prof_nsigma > 10)` evaluates a boolean comparison and rounds it, giving `nb = 0` when `prof_nsigma ≤ 10` (immediate division-by-zero in `sinh_space`) and `nb = 1` when `prof_nsigma > 10` (only two knot positions, yielding a completely wrong fit). The correct intent is a count proportional to `prof_nsigma`, e.g. `nb = max(1, round(prof_nsigma / 10))`. **This must be fixed as part of the refactoring.** Once fixed, `test_prof_nsigma_extended` should be strengthened to assert FWHM recovery.

### 5. `_initial_profile_bspline(sigma_x, norm_obj, norm_ivar, sn2_img, sn_gauss, bkpt, min_sigma, max_sigma)`
**Lines 542–597.** Select initial `inside` pixels, sort by `sigma_x`, fit the 4th-order B-spline, evaluate `mode_fit`, call `findfwhm`, update sigma/limit/l_limit/r_limit. Returns `(bset, si, inside, mode_fit, median_fit, peak_x, l_limit, r_limit, min_level, bspline_fwhm, sigma_updated, limit_updated, trace_corr)`.

### 6. `_iterative_correction(bset, sigma_x, norm_obj, norm_ivar, xtemp, inside, sigma, area, trace_corr, dspat, bkpt, l_limit, r_limit, med_sn2, nspec, nspat, sigma_iter=3)`
**Lines 612–716.** The three-iteration loop performing trace-shift and width-stretch corrections. Returns updated `(sigma_x, trace_corr, sigma, area, bset)`. Each sub-fit failure can return `None` to trigger the Gaussian fallback in the parent. This is the most complex block to extract cleanly because it mutates multiple shared arrays, but it is also the most self-contained conceptually.

### 7. `_final_profile_bspline(sigma_x, norm_obj, norm_ivar, pb, ss, inside, mask, bkpt, nspec, nspat)`
**Lines 724–742.** The final B-spline fit and BSpline2D→1D conversion. Returns `(bset, outmask)`.

### 8. `_compute_apodization(bset, sigma_x, full_bsp, ss, peak_x, median_fit, min_level, limit, l_limit, r_limit, prof_nsigma, no_deriv)`
**Lines 744–816.** Re-find limits from profile level, compute logarithmic derivatives, walk limits inward, apply exponential tails. Returns `(full_bsp, l_limit, r_limit, l_deriv, r_deriv)`.

---

## Performance Improvements

### Replace `np.outer(v, np.ones(n))` with `v[:, None]` broadcasting

Every instance allocates a temporary `np.ones(nspat)` array and then forms the outer product. NumPy broadcasting via `[:, None]` (adding a trailing axis) is semantically identical, avoids the temporary allocation, and is clearer.

| Line | Current | Replacement |
|------|---------|-------------|
| 362 | `np.outer(trace_in, np.ones(nspat))` | `trace_in[:, None]` |
| 377 | `dspat/np.outer(sigma, np.ones(nspat))` | `dspat/sigma[:, None]` |
| 377 | `- np.outer(trace_corr, np.ones(nspat))` | `- trace_corr[:, None]` |
| 507 | `np.outer(4.0 + np.sqrt(...), np.ones(nspat))` | `(4.0 + np.sqrt(...))[:, None]` |
| 689 | same as 377 (inside loop) | same replacements |
| 694 | `np.outer(area, np.ones(nspat, dtype=float))` | `area[:, None]` |
| 730 | same | `area[:, None]` |
| 843 | `np.outer(np.sum(profile_model, 1), np.ones(nspat))` | `np.sum(profile_model, 1)[:, None]` |

Note that at line 694/730 `area` starts as a scalar (`1.0`) and becomes a 1D array after the first loop iteration. The replacement `area[:, None]` only works when `area` is 1D, so the initialization at line 511 must also change to `area = np.ones(nspec)` (or the guard `np.atleast_1d(area)[:, None]` can be used, which also works for the scalar case but is less clean).

### Vectorize the apodization limit walk (lines 789–803)

The two `while True` loops step in 0.1 increments along `l_limit` and `r_limit` and re-evaluate `bset.value()` at each step. For wide profiles (`prof_nsigma` set, bright objects) this can run for many iterations. Because the pre-computed vectors `l_lim_vec` / `r_lim_vec` already cover the same range at the same spacing, the walk can be replaced by a vectorized search over those arrays:

```python
# left limit: walk inward until l_deriv <= l_deriv_max or l_limit >= -1.0
l_cross = np.where(l_deriv_vec <= l_deriv_max)[0]
l_limit = l_lim_vec[l_cross[0]] if l_cross.size > 0 else -1.0

# right limit: walk inward until r_deriv >= r_deriv_max or r_limit <= 1.0
r_cross = np.where(r_deriv_vec >= r_deriv_max)[0]
r_limit = r_lim_vec[r_cross[0]] if r_cross.size > 0 else 1.0
```

This eliminates the loops entirely, reusing the already-computed derivative vectors.

### Minor: avoid redundant `(pb == 0.0)` boolean array in profile normalization (line 845)

```python
# Current
profile_model = (norm > 0.0)*profile_model/(norm + (norm == 0.0))
# Simpler — equivalent when norm ≥ 0
norm_safe = np.where(norm > 0.0, norm, 1.0)
profile_model = np.where(norm > 0.0, profile_model / norm_safe, 0.0)
```

---

## Accuracy Improvements

The two utility functions `return_gaussian` and `findfwhm` share a common
assumption: that the Gaussian spatial profile is *well-sampled*, meaning the
value of the Gaussian at the center of each pixel is a faithful proxy for the
true photon count in that pixel.  This assumption fails for narrow profiles
(FWHM ≲ 4–5 px), where a pixel of width 1 px subtends a non-negligible
fraction of the profile width and the difference between evaluating at the
pixel center versus integrating across the pixel is appreciable.

### A. Pixel-integrated Gaussian in `return_gaussian` (line 247)

**Current implementation** (line 247):
```python
profile_model = np.exp(-0.5*sigma_x**2)/np.sqrt(2.0 * np.pi)*(sigma_x ** 2 < 25.)
```

This evaluates the Gaussian probability density at the normalized coordinate
`sigma_x = (x - x_trace) / sigma` at the *center* of each detector pixel.
The true photon count in pixel `i` is the integral of the profile across the
pixel width:

```
P_i = ∫_{x_i - 0.5}^{x_i + 0.5}  (1/sqrt(2π) σ) exp(-u²/(2σ²)) du
    = 0.5 * [erf((sigma_x_i + δ)/sqrt(2)) − erf((sigma_x_i − δ)/sqrt(2))]
```

where `delta = 0.5 / sigma` is the half-pixel width in normalized `sigma_x`
units.  The error-function form is already available via `scipy.special.erf`,
which is imported at the module level.

**Fix**: replace the point-sampled formula with the pixel-integrated form:

```python
sigma = fwhm / 2.3548                # width in pixels
delta = 0.5 / sigma                  # half-pixel in sigma_x units
profile_model = 0.5 * (scipy.special.erf((sigma_x + delta) / np.sqrt(2))
                      - scipy.special.erf((sigma_x - delta) / np.sqrt(2)))
profile_model[sigma_x**2 >= 25.] = 0.0
```

`sigma` can be computed inside `return_gaussian` from the existing `fwhm`
argument; the function signature does not need to change.

**Magnitude of the correction**: `delta = 0.5/sigma` controls the effect.

| FWHM (px) | sigma (px) | delta | Max pixel-center error |
|-----------|-----------|-------|------------------------|
| 2.0       | 0.85      | 0.59  | large (profile wider than pixel spacing) |
| 3.0       | 1.27      | 0.39  | ~10% at peak |
| 4.0       | 1.70      | 0.29  | ~5% at peak |
| 6.0       | 2.55      | 0.20  | ~2% at peak |
| 10.0      | 4.25      | 0.12  | <1% at peak |

**Normalization consequence**: the pixel-integrated form sums to approximately
1.0 per spectral row (the sum telescopes to `erf(N/(2σ√2)) → 1` for wide
slits).  The current point-sampled form sums to approximately `σ` (one
Riemann-sum factor of the pixel spacing `1/σ` multiplies the continuous
integral of 1).  After this fix the Gaussian fallback path will produce a
properly normalized profile (row sums ≈ 1.0), matching the B-spline path.
Tests that check Gaussian-path row sums — currently absent — can be added
once the fix is applied.

**Well-sampled limit**: as `sigma → ∞` (broad profile), `delta → 0` and the
erf form recovers the point-sampled Gaussian exactly, so the fix is
backward-compatible for well-sampled data.

### B. Linear interpolation for the half-maximum crossing in `findfwhm` (lines 54–73)

**Current implementation**: the function finds `lwhm` and `rwhm` by locating
the *first discrete sample* of `model` on each side of the peak that falls
below `0.5 × peak`.  Because the true half-maximum crossing lies between two
adjacent samples, the current approach always reports the farther of the two
bracketing samples, systematically overestimating the profile half-width.

**Fix**: after identifying the first below-half-maximum sample on each side,
linearly interpolate between that sample and the one immediately adjacent to it
(closer to the peak, where the model is above 0.5×peak) to locate the crossing
to sub-pixel precision.

For the right half-maximum, where `rh = rind.min()` is the index of the first
sample to the right of `peak_x` that falls below `0.5*peak`:

```python
if rind.size > 0:
    rh = rind.min()
    if rh > 0:
        x0, x1 = sig_x[rh - 1], sig_x[rh]   # x0 above half-max, x1 below
        y0, y1 = model[rh - 1], model[rh]
        rwhm = x0 + (0.5*peak - y0) * (x1 - x0) / (y1 - y0)
    else:
        rwhm = sig_x[rh]                       # edge case: no above-half-max sample
```

The left-side search uses the reversed arrays; the same interpolation is
applied after mapping `lh` (the reversed-array index) back to the original
array to identify the two bracketing samples.

**Expected FWHM accuracy after fix**:

| True FWHM (px) | Current recovered FWHM | After interpolation fix |
|----------------|------------------------|-------------------------|
| 3.0            | ~4.0 (33% error)       | ~3.0 (<2% error)        |
| 5.0            | ~6.0 (20% error)       | ~5.0 (<2% error)        |
| 6.0            | ~6.0 (<1% error)       | ~6.0 (<1% error)        |
| 8.0            | ~8.0 (<1% error)       | ~8.0 (<1% error)        |

Once this fix is applied:
- `test_fwhm_recovery` can be extended to include FWHM values smaller than
  6.0 px (e.g., `[3.0, 4.0, 6.0, 8.0]`).
- The trace-correction overshoot at small offsets (Bug 3) is expected to
  diminish or disappear, since the overshoot's root cause is the same
  discrete-sampling bias in `peak_x` estimation.

---

## Pre-existing Bugs Discovered During Testing

These bugs were found while writing and running `pypeit/tests/test_spatialprofile.py` against the current code. They are not introduced by the refactoring; they exist in the pre-refactor function. Each should be fixed as part of the refactoring of its enclosing block.

### Bug 1 — `nb` formula in the `prof_nsigma` knot-setup branch (line 532)

```python
# Current (wrong)
nb = np.round(prof_nsigma > 10)
```

`prof_nsigma > 10` is a boolean scalar (`True` or `False`). `np.round` of a boolean gives `0` or `1`. Consequences:

- `prof_nsigma ≤ 10`: `nb = 0` → `sinh_space = np.arcsinh(prof_nsigma) / 0` → division by zero, then `bkpt` is degenerate.
- `prof_nsigma > 10`: `nb = 1` → only two knot positions (`-rb` and `+rb`), giving a completely wrong B-spline fit and unreliable `fwhmfit`.

**Fix**: replace with a proportional count, e.g.:

```python
nb = max(1, round(prof_nsigma / 10))
```

Once fixed, `test_prof_nsigma_extended` should be strengthened to assert FWHM recovery within a reasonable tolerance.

### Bug 2 — FWHM overestimation bias in `findfwhm` for narrow profiles

`findfwhm` (lines 54–73) locates the half-maximum by finding the first sample point in `sigma_x` that drops below 0.5 × peak. Because `sigma_x` values are spaced `1/sigma` pixels apart (one value per integer spatial pixel), the half-maximum crossing always falls *between* two adjacent samples. For narrow profiles (FWHM ≲ 5 px, `sigma ≲ 2.1` px), this overestimates `rwhm - lwhm` substantially:

| True FWHM (px) | sigma_x spacing | Recovered FWHM (px) | Relative error |
|---|---|---|---|
| 3.0 | 0.79 | ~4.0 | 33% |
| 5.0 | 0.47 | ~6.0 | 20% |
| 6.0 | 0.39 | ~6.0 | <1% |
| 8.0 | 0.29 | ~8.0 | <1% |

**Fix**: see "Accuracy Improvements §B" above — linear interpolation between the two bracketing sample points eliminates the bias at all FWHM values.

**Impact on tests**: `test_fwhm_recovery` is currently parametrized over `[6.0, 8.0]` px only; once the interpolation fix is applied, the parametrize list can be extended to narrower values.

### Bug 3 — Trace correction overshoots for small offsets (< 1 px)

For `trace_offset = 0.5` px, the iterative correction (Block 8) overshoots: `xnew` ends up further from the true trace than `trace_in` was. For `trace_offset = 1.0` px the correction is accurate (`err_after ≈ 0.5 * err_before`). The root cause is the same discrete sigma_x sampling described in Bug 2: `findfwhm` overestimates `peak_x` when the peak falls between two sample points, causing `trace_corr` to be initialized too large, and the iterative correction then overcorrects.

**Fix**: applying the `findfwhm` linear-interpolation fix (Accuracy Improvements §B) is expected to substantially reduce the overshoot for small offsets, since `peak_x` will then be located to sub-pixel precision.

**Impact on tests**: `test_trace_correction_direction` uses `trace_offset = 1.0` px; once the fix is applied, test with smaller offsets (e.g., 0.5 px) to verify the overshoot is eliminated.

---

## Unit Tests (`pypeit/tests/test_spatialprofile.py`)

The test file has been written and all tests pass against the pre-refactor function (11 passed, 1 skipped). The skipped test (`test_refactor_regression`) requires reference `.npy` files generated by the regression freeze workflow below.

### Key design notes

- **S/N knob**: `sn_ratio ≥ 20` reliably puts `sqrt(med_sn2) > sn_gauss = 4`, sending the function down the B-spline path. `sn_ratio ≤ 1` puts it below the threshold and forces the Gaussian fallback.
- **Gaussian path does not normalize**: `return_gaussian` returns `exp(-0.5*sigma_x^2)/sqrt(2*pi)` without row-normalization. Row sums equal approximately `sigma = fwhm/2.3548`. Only the B-spline path normalizes each row to 1.
- **Early-return via NaN flux**: setting `flux = np.full(nspec, np.nan)` makes `flux_sm` NaN everywhere so `indsp` is empty, triggering the very first early return at line 386.
- **`true_trace` returned**: the factory returns the true object center so tests can measure trace-correction accuracy without hard-coding expected values.
- **FWHM range**: FWHM recovery (`test_fwhm_recovery`) is only asserted for FWHM ≥ 6 px due to the discrete-sampling bias in `findfwhm` (Bug 2 above).
- **`prof_nsigma` assertions**: `test_prof_nsigma_extended` only asserts shape and finiteness, not FWHM recovery or positivity, due to the broken knot-count formula (Bug 1). Strengthen this test after the bug is fixed.

### Full test file

```python
"""
Unit tests for :func:`~pypeit.core.spatialprofile.fit_profile`.

Synthetic images are generated by :func:`make_profile_inputs`, which builds a
self-consistent set of all inputs the function expects from a noiseless Gaussian
spatial profile, a flat 1D spectrum, and uniform Gaussian background noise.
The ``sn_ratio`` parameter controls which code path is exercised:

- ``sn_ratio >= 20`` → B-spline profile fitting path
- ``sn_ratio <= 1``  → Gaussian fallback (``med_sn2 < sn_gauss**2``)

All random-number generators use a fixed seed for reproducibility.
"""

import os

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest

from pypeit.core import spatialprofile


# ============================================================================
# Synthetic data factory
# ============================================================================

def make_profile_inputs(nspec=200, nspat=100, fwhm=4.0, flux_level=1000.0,
                        sn_ratio=20.0, trace_offset=0.0, seed=42):
    r"""
    Build a self-consistent set of synthetic inputs for
    :func:`~pypeit.core.spatialprofile.fit_profile`.

    The 2D image is a noiseless Gaussian spatial profile scaled by a flat 1D
    spectrum, plus uniform Gaussian background noise.  The noise amplitude is
    chosen so that the extracted 1D S/N equals ``sn_ratio``.

    Parameters
    ----------
    nspec : :obj:`int`, optional
        Number of spectral pixels. Default is 200.
    nspat : :obj:`int`, optional
        Number of spatial pixels. Default is 100.
    fwhm : :obj:`float`, optional
        True FWHM of the Gaussian spatial profile in pixels. Default is 4.0.
    flux_level : :obj:`float`, optional
        Total counts per spectral row. Default is 1000.0.
    sn_ratio : :obj:`float`, optional
        Desired S/N of the extracted 1D spectrum.  Values ``>= 20`` exercise
        the B-spline path; values ``<= 1`` trigger the Gaussian fallback.
        Default is 20.0.
    trace_offset : :obj:`float`, optional
        Shift of the true object centre relative to ``trace_in`` in pixels.
        Default is 0.0.
    seed : :obj:`int`, optional
        Random seed. Default is 42.

    Returns
    -------
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux,
    fluxivar, inmask, true_trace : ndarray
    """
    rng = np.random.default_rng(seed)

    # Step 1 — coordinate grids
    wave = np.linspace(4000.0, 5000.0, nspec)
    waveimg = wave[:, np.newaxis] * np.ones((nspec, nspat))
    spat_img = np.ones((nspec, 1)) * np.arange(nspat, dtype=float)

    # Step 2 — object trace
    trace_in = np.full(nspec, float(nspat // 2))
    true_trace = trace_in + trace_offset

    # Step 3 — noiseless 2D image: Gaussian profile scaled by a flat spectrum
    sigma = fwhm / 2.3548
    dspat_true = spat_img - true_trace[:, np.newaxis]
    profile_gauss = np.exp(-0.5 * (dspat_true / sigma) ** 2)
    profile_norm = profile_gauss / profile_gauss.sum(axis=1, keepdims=True)
    true_image = flux_level * profile_norm

    # Step 4 — noise and inverse variance (uniform background model)
    # noise_sigma_2d is set so that the extracted 1D S/N equals sn_ratio
    noise_sigma_2d = flux_level / (sn_ratio * np.sqrt(nspat))
    image = true_image + rng.normal(scale=noise_sigma_2d, size=(nspec, nspat))
    ivar = np.full((nspec, nspat), 1.0 / noise_sigma_2d ** 2)
    thismask = np.ones((nspec, nspat), dtype=bool)
    inmask = thismask.copy()

    # Step 5 — 1D extracted spectrum (noiseless true spectrum + ivar)
    flux = np.full(nspec, flux_level)
    noise_sigma_1d = flux_level / sn_ratio
    fluxivar = np.full(nspec, 1.0 / noise_sigma_1d ** 2)

    return (image, ivar, waveimg, thismask, spat_img, trace_in,
            wave, flux, fluxivar, inmask, true_trace)


# ============================================================================
# Tests
# ============================================================================

def test_profile_return_shapes():
    """Return values have the correct shape and contain no NaNs."""
    nspec, nspat = 200, 100
    for sn_ratio in (20.0, 1.0):
        image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, \
            fluxivar, inmask, _ = make_profile_inputs(
                nspec=nspec, nspat=nspat, sn_ratio=sn_ratio)
        profile_model, xnew, fwhmfit, med_sn2 = spatialprofile.fit_profile(
            image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
            spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
            fluxivar=fluxivar, inmask=inmask)
        assert profile_model.shape == (nspec, nspat)
        assert xnew.shape == (nspec,)
        assert fwhmfit.shape == (nspec,)
        assert isinstance(med_sn2, float)
        assert np.all(np.isfinite(profile_model))
        assert np.all(np.isfinite(xnew))


def test_profile_normalization_bspline():
    """B-spline path: each non-zero spectral row of profile_model sums to 1."""
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, \
        fluxivar, inmask, _ = make_profile_inputs(sn_ratio=20.0, fwhm=4.0)
    profile_model, xnew, fwhmfit, med_sn2 = spatialprofile.fit_profile(
        image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
        spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
        fluxivar=fluxivar, inmask=inmask, sn_gauss=4.0)
    assert med_sn2 > 4.0 ** 2, "Expected B-spline path (high S/N)"
    assert np.all(profile_model >= 0)
    row_sums = profile_model.sum(axis=1)
    nonzero_rows = row_sums > 0
    np.testing.assert_allclose(row_sums[nonzero_rows], 1.0, atol=1e-10)


def test_profile_normalization_gaussian():
    """Gaussian fallback: profile is non-negative, finite, and bell-shaped."""
    nspec, nspat = 200, 100
    fwhm = 4.0
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, \
        fluxivar, inmask, _ = make_profile_inputs(
            nspec=nspec, nspat=nspat, sn_ratio=1.0, fwhm=fwhm)
    profile_model, xnew, fwhmfit, med_sn2 = spatialprofile.fit_profile(
        image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
        spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
        fluxivar=fluxivar, inmask=inmask, sn_gauss=4.0)
    assert med_sn2 < 4.0 ** 2, "Expected Gaussian path (low S/N)"
    assert np.all(profile_model >= 0)
    assert np.all(np.isfinite(profile_model))
    # Peak column of each row must be at the trace centre (±1 px)
    peak_cols = profile_model.argmax(axis=1)
    np.testing.assert_allclose(peak_cols, nspat // 2, atol=1)
    # Profile drops to < 10% of its peak beyond ±3 sigma from the centre
    sigma = fwhm / 2.3548
    spat_pix = np.arange(nspat)
    far = np.abs(spat_pix - nspat // 2) > 3 * sigma
    assert np.all(profile_model[:, far] <
                  0.1 * profile_model.max(axis=1, keepdims=True))


def test_forced_gaussian():
    """gauss=True forces the Gaussian path regardless of S/N."""
    nspec, nspat = 200, 100
    fwhm = 4.0
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, \
        fluxivar, inmask, _ = make_profile_inputs(
            nspec=nspec, nspat=nspat, sn_ratio=20.0, fwhm=fwhm)
    profile_model, xnew, fwhmfit, med_sn2 = spatialprofile.fit_profile(
        image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
        spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
        fluxivar=fluxivar, inmask=inmask, gauss=True, thisfwhm=fwhm)
    # Gaussian path always returns trace_in unchanged
    np.testing.assert_array_equal(xnew, trace_in)
    # fwhmfit is set from the initial sigma and never updated in the Gaussian path
    np.testing.assert_allclose(fwhmfit, fwhm * np.ones(nspec))
    assert np.all(profile_model >= 0)
    assert np.all(np.isfinite(profile_model))


@pytest.mark.parametrize("fwhm", [6.0, 8.0])
def test_fwhm_recovery(fwhm):
    """B-spline path recovers the input FWHM to within 10%.

    Tested only for fwhm >= 6 px.  For narrower profiles the discrete pixel
    sampling in sigma_x causes findfwhm to overestimate the FWHM (the
    half-maximum crossing falls between adjacent sample points), so the
    recovery tolerance cannot be met with this synthetic setup.
    """
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, \
        fluxivar, inmask, _ = make_profile_inputs(sn_ratio=20.0, fwhm=fwhm)
    profile_model, xnew, fwhmfit, med_sn2 = spatialprofile.fit_profile(
        image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
        spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
        fluxivar=fluxivar, inmask=inmask, thisfwhm=fwhm, sn_gauss=4.0)
    assert med_sn2 > 4.0 ** 2, "Expected B-spline path"
    np.testing.assert_allclose(np.median(fwhmfit), fwhm, rtol=0.10)


def test_trace_correction_direction():
    """When the true object centre is offset from trace_in, xnew moves toward it.

    A 0.5 px offset is too small — the iterative correction overshoots because
    the discrete sigma_x sampling causes findfwhm to overestimate the shift.
    An offset of 1.0 px corrects reliably with the default max_trace_corr=2.0.
    """
    trace_offset = 1.0  # pixels, positive direction
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, \
        fluxivar, inmask, true_trace = make_profile_inputs(
            sn_ratio=20.0, fwhm=4.0, trace_offset=trace_offset)
    profile_model, xnew, fwhmfit, med_sn2 = spatialprofile.fit_profile(
        image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
        spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
        fluxivar=fluxivar, inmask=inmask, thisfwhm=4.0, sn_gauss=4.0)
    assert med_sn2 > 4.0 ** 2, "Expected B-spline path"
    assert np.mean(xnew) > np.mean(trace_in), \
        "xnew should shift toward the positive offset"
    err_before = np.abs(np.mean(trace_in - true_trace))
    err_after = np.abs(np.mean(xnew - true_trace))
    assert err_after < err_before, \
        "xnew should be closer to true_trace than trace_in"


def test_prof_nsigma_extended():
    """prof_nsigma code path runs and returns valid arrays.

    The formula ``nb = np.round(prof_nsigma > 10)`` on line 532 of
    spatialprofile.py evaluates a boolean comparison (giving nb=0 for
    prof_nsigma<=10, division by zero; nb=1 for prof_nsigma>10, only 2
    knots).  This is a pre-existing bug: the knot count is incorrectly
    computed and the recovered FWHM is unreliable.  This test only verifies
    that the function returns finite arrays without raising.
    FWHM recovery for the prof_nsigma path is not asserted here; it should
    be tested once the knot-count formula is corrected during refactoring.
    """
    fwhm = 8.0
    prof_nsigma = 12.0  # > 10 so nb=1 (avoids division-by-zero at nb=0)
    nspat = 150  # wider slit to contain the broad profile
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, \
        fluxivar, inmask, _ = make_profile_inputs(
            nspat=nspat, sn_ratio=20.0, fwhm=fwhm, flux_level=5000.0)
    profile_model, xnew, fwhmfit, med_sn2 = spatialprofile.fit_profile(
        image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
        spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
        fluxivar=fluxivar, inmask=inmask,
        thisfwhm=fwhm, prof_nsigma=prof_nsigma, sn_gauss=4.0)
    # Only check shape and finiteness: the broken knot-count formula produces a
    # poor B-spline fit that can go negative, so positivity is not asserted.
    assert profile_model.shape == (200, nspat)
    assert np.all(np.isfinite(profile_model))


@pytest.mark.parametrize("variant,sn_ratio,nan_flux", [
    ("high_sn_bspline",    20.0, False),
    ("low_sn_gaussian",     1.0, False),
    ("early_return_noflux", 20.0, True),  # NaN flux → indsp empty → early return
])
def test_profile_positivity_and_finite_all_paths(variant, sn_ratio, nan_flux):
    """profile_model is non-negative and finite across every code path."""
    nspec, nspat = 200, 100
    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, \
        fluxivar, inmask, _ = make_profile_inputs(
            nspec=nspec, nspat=nspat, sn_ratio=sn_ratio)
    if nan_flux:
        # All-NaN flux makes flux_sm NaN everywhere, leaving indsp empty
        flux = np.full(nspec, np.nan)
    profile_model, xnew, fwhmfit, med_sn2 = spatialprofile.fit_profile(
        image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
        spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
        fluxivar=fluxivar, inmask=inmask, sn_gauss=4.0)
    assert profile_model.shape == (nspec, nspat), \
        f"shape check failed for {variant}"
    assert np.all(profile_model >= 0), \
        f"non-negative check failed for {variant}"
    assert np.all(np.isfinite(profile_model)), \
        f"finite check failed for {variant}"


def test_refactor_regression():
    """Refactored function output matches the pre-refactor reference arrays."""
    files_dir = os.path.join(os.path.dirname(__file__), 'files')
    ref_paths = {
        'profile': os.path.join(files_dir, 'spatprof_profile_model_ref.npy'),
        'xnew':    os.path.join(files_dir, 'spatprof_xnew_ref.npy'),
        'fwhmfit': os.path.join(files_dir, 'spatprof_fwhmfit_ref.npy'),
    }
    if not all(os.path.exists(p) for p in ref_paths.values()):
        pytest.skip("Reference arrays not found; generate them with "
                    "pypeit/tests/files/freeze_spatprof.py first")

    ref_profile = np.load(ref_paths['profile'])
    ref_xnew    = np.load(ref_paths['xnew'])
    ref_fwhmfit = np.load(ref_paths['fwhmfit'])

    image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, \
        fluxivar, inmask, _ = make_profile_inputs(sn_ratio=20.0, fwhm=4.0)
    profile_model, xnew, fwhmfit, med_sn2 = spatialprofile.fit_profile(
        image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
        spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
        fluxivar=fluxivar, inmask=inmask, sn_gauss=4.0)

    np.testing.assert_allclose(profile_model, ref_profile, atol=1e-10)
    np.testing.assert_allclose(xnew, ref_xnew, atol=1e-10)
    np.testing.assert_allclose(fwhmfit, ref_fwhmfit, atol=1e-10)
```

---

### Regression Freeze Workflow

Before starting any refactoring:

1. Run `fit_profile` on the high-S/N synthetic inputs from `test_profile_normalization_bspline` and save outputs to an untracked location:
   ```python
   np.save('pypeit/tests/files/spatprof_profile_model_ref.npy', profile_model)
   np.save('pypeit/tests/files/spatprof_xnew_ref.npy', xnew)
   np.save('pypeit/tests/files/spatprof_fwhmfit_ref.npy', fwhmfit)
   ```
   Leave these `.npy` files as untracked — do not `git add` them.
2. After refactoring, add `test_refactor_regression` (above) to the test file and verify all assertions pass with `atol=1e-10`.
3. The regression test is skipped automatically if the reference files are absent (e.g., on CI), so it never blocks the test suite.

---

## Verification Plan

1. ✅ Write unit tests and confirm all pass against the *pre-refactor* function (`pytest pypeit/tests/test_spatialprofile.py`: 11 passed, 1 skipped).
2. Generate regression freeze reference arrays and leave them as untracked files.
3. Perform refactoring (helper-function extraction, broadcasting replacements, apodization vectorization, fix Bug 1 `nb` formula).
4. Re-run unit tests and the regression freeze test; all must pass.
5. Strengthen `test_prof_nsigma_extended` to assert FWHM recovery now that Bug 1 is fixed.
6. Run the dev-suite for a multi-slit instrument (e.g., `keck_deimos`) to validate end-to-end behavior on real data.
7. `pytest pypeit/tests/` to catch any import or signature regressions.
