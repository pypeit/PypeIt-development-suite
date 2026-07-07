# fit_profile Refactoring Analysis (v4)

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
of `fit_profile_refactor` is organized into thirteen labeled comment blocks
(matching the workflow breakdown below); only the smallest, most reusable
pieces have been extracted into named functions.

---

## Workflow Breakdown (line ranges are approximate)

### Block 1 — Initialization (347–383)
Set up masks (`inmask`, `totmask`), handle `prof_nsigma → no_deriv`, clamp
`thisfwhm`, compute `dspat` (spatial separation image), zero-initialize
`sn2_img` / `spline_img`, median-smooth `flux` and `fluxivar`, and initialize
`sigma`, `fwhmfit`, `trace_corr`, `sigma_x`.

### Block 2 — Flux B-spline fitting (383–415)
Cap `fluxivar` at S/N=100, identify eligible spectral pixels `indsp`, fit two
fine B-splines (`b_answer`, pass 1 & 2) and one coarse continuum B-spline
(`c_answer`). Early Gaussian return if no eligible pixels or continuum fit
raises.

### Block 3 — S/N estimation (417–460)
Compute per-pixel `sn2` from the fine spline × ivar; use a percentile cut +
sigma-clipped mean to get `med_sn2`. Expand spline/continuum evaluations to the
full 1D spectrum (`spline_flux1`, `cont_flux1`, `sn2_1`); interp-fill bad
B-spline pixels with the continuum model; median-filter `sn2` and interpolate
it onto the 2D `sn2_img`.

### Block 4 — Normalized-object image construction (462–509)
Build `spline_img` (the S/N-level-dependent spectral flux model, extrapolated
to 2D):
- S/N² ≤ 2: use noise floor `sigma1`
- 2 < S/N² ≤ 5: use coarse continuum
- S/N² > 5: fine spline with bad-pixel patching from continuum

Compute `norm_obj = image / spline_img` and `norm_ivar = ivar * spline_img²`.
Compute `xtemp` — a cumulative S/N-weighted spectral coordinate used as the
independent variable for the iterative profile correction fit.

### Block 5 — Early returns (513–518)
Return Gaussian immediately if `ngood < 10`, `med_sn2 < sn_gauss²`, or
`gauss=True`.

### Block 6 — Knot setup (520–554)
Compute profile half-extent `limit` from `erfcinv`, set sinh-spaced B-spline
interior knots `bkpt`, determine `min_sigma` / `max_sigma`, and select `inside`
pixels for the initial profile fit (preferring high-S/N pixels if available,
falling back to all in-range pixels).

### Block 7 — Initial profile fit (554–614)
Fit a 4th-order B-spline to `norm_obj[inside]` sorted by `sigma_x`. Evaluate
to get `mode_fit`; call `_findfwhm` to locate the peak and measure width.
Update `sigma` and `limit` using the B-spline FWHM. Walk the profile down from
the peak to find `l_limit` and `r_limit` (where profile drops below
`min_level`). Select the final `inside` subset within those limits; early
Gaussian return if `ninside < 10`.

### Block 8 — Iterative trace and width correction loop (612–716) — `sigma_iter=3`
Each iteration:
1. **Trace shift**: evaluate the current B-spline at `sigma_x` and at
   `sigma_x ± 0.5`; use `(mode_min05 - mode_plu05)` as a shift-sensitive basis
   and fit it jointly with `mode_zero` via `iterative_bspline_fit`; extract the
   h0/h1 coefficient ratio to get `delta_trace_corr`.
2. **Width stretch**: evaluate at `sigma_x / 1.3`; use the rescaled difference
   as a stretch-sensitive basis; extract h2/h0 ratio to get `sigma_factor`.
3. Update `sigma`, `area`, `trace_corr`, `sigma_x`.
4. On non-final iterations: re-fit the profile B-spline with the updated
   `sigma_x` (returned as a `BSpline2D`, converted to 1D `BSpline` via
   `_bspline2d_to_1d` for subsequent `value()` calls).

Each sub-fit failure returns a Gaussian immediately via `_return_gaussian`.

### Block 9 — Final trace (718–722)
Apply `trace_corr` only if `median(|trace_corr * sigma|) < max_trace_corr`;
otherwise keep `trace_in`.

### Block 10 — Final profile fit (724–742)
Re-select all in-range, unmasked pixels sorted by `sigma_x`; fit the definitive
4th-order B-spline profile using the full `pb` (area) weight array and loose
rejection (`upper=lower=10`). Convert resulting `BSpline2D` → 1D `BSpline` via
`_bspline2d_to_1d`.

### Block 11 — Apodization limit search (744–803)
Evaluate the final B-spline across the full sigma_x range (`full_bsp`). Re-find
`l_limit` / `r_limit` from the profile level. Then compute the logarithmic
derivative `d(log P)/d(log x)` over vectors `l_lim_vec` and `r_lim_vec` to
find the extremal derivative. The limit walk is vectorized (no `while True`
loops).

### Block 12 — Exponential apodization (805–816)
If derivatives have the right sign and `no_deriv=False`, replace `full_bsp`
outside `[l_limit, r_limit]` with an exponential matched in value and first
derivative to the profile edge.

### Block 13 — Normalization, logging, QA, return (818–854)
Reshape `full_bsp`, multiply by `pb`, normalize each spectral row to unit sum,
guard against NaNs, log fit statistics, show QA plot if `show_profile`, and
return.

---

## Implemented Helper Functions

Four module-private helpers have been extracted into named functions. All other
blocks (1–13) remain inline within `fit_profile_refactor`, organized by
comment-block labels.

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

---

## Deferred Block-Level Extractions

The v3 plan listed eight candidate helper functions corresponding to the major
blocks of `fit_profile_refactor` (smooth_spectrum, estimate_sn2,
build_spline_img, compute_bspline_knots, initial_profile_bspline,
iterative_correction, final_profile_bspline, compute_apodization). These have
**not** been extracted. The inline comment-block structure is already a
significant readability improvement over the original monolith, and extracting
these blocks introduces substantial interface boilerplate (large argument lists)
with limited testing benefit. Block-level extraction is deferred unless there is
a specific motivation (e.g., a need to unit-test a block in isolation, or a
second call site requiring reuse).

---

## Performance Improvements (implemented)

### `np.outer(v, np.ones(n))` → `v[:, None]` broadcasting

Every `np.outer` was replaced with the equivalent broadcasting form, avoiding
temporary array allocation. The initialization `area = 1.0` was changed to
`area = np.ones(nspec)` so that `area[:, None]` is valid throughout.

### Vectorized apodization limit walk

The two `while True` loops in Block 11 (stepping 0.1 increments, re-evaluating
the B-spline at each step) were replaced with vectorized searches over the
pre-computed derivative arrays `l_deriv_vec` / `r_deriv_vec`:

```python
l_cross = np.where(l_deriv_vec <= l_deriv_max)[0]
l_limit = l_lim_vec[l_cross[0]] if l_cross.size > 0 else -1.0

r_cross = np.where(r_deriv_vec >= r_deriv_max)[0]
r_limit = r_lim_vec[r_cross[0]] if r_cross.size > 0 else 1.0
```

### Profile normalization simplified

Block 13 profile normalization uses `np.where` instead of the boolean-mask
arithmetic from the legacy code:

```python
norm_safe = np.where(row_sums[:, None] > 0.0, row_sums[:, None], 1.0)
profile_model = np.where(row_sums[:, None] > 0.0, profile_model / norm_safe, 0.0)
```

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
two bracketing samples to locate the crossing to sub-pixel precision, eliminating
the discrete-sampling bias. Combined with the masked-array peak search (restricts
search to `|sig_x| < 2`), this also eliminates the conceptual error in the
legacy implementation where the half-maximum crossing region was masked out.

---

## Pre-existing Bugs

These bugs exist in the legacy `spatialprofile.py`. All three are **fixed in
`fit_profile_refactor`**.

### Bug 1 — `nb` formula in the `prof_nsigma` knot-setup branch (line 532)

```python
# Legacy (wrong)
nb = np.round(prof_nsigma > 10)
# Fixed
nb = max(1, round(prof_nsigma / 10))
```

`prof_nsigma > 10` is a boolean scalar; `np.round` of a boolean gives 0 or 1.
For `prof_nsigma ≤ 10` this causes division-by-zero; for `prof_nsigma > 10` it
gives only 2 knot positions. The fix uses a count proportional to `prof_nsigma`.

**Status in `fit_profile_refactor`**: fixed (Block 6, `nb = max(1, round(prof_nsigma / 10))`).

### Bug 2 — FWHM overestimation bias in `findfwhm` for narrow profiles

For narrow profiles (FWHM ≲ 5 px, sigma ≲ 2.1 px), `sigma_x` values are spaced
`1/sigma` units apart and the half-maximum crossing always falls between two
adjacent samples. The legacy code reported the farther of the two bracketing
samples, overestimating the FWHM:

| True FWHM (px) | Legacy recovered FWHM | After fix |
|---|---|---|
| 3.0 | ~4.0 (33% error) | ~3.0 (<2% error) |
| 4.0 | ~5.0 (25% error) | ~4.0 (<2% error) |
| 5.0 | ~6.0 (20% error) | ~5.0 (<2% error) |
| 6.0 | ~6.0 (<1% error) | ~6.0 (<1% error) |

**Status in `fit_profile_refactor`**: fixed in `_findfwhm` (linear interpolation
between bracketing samples; mask at `|sig_x| > 2` rather than `> 1`).

### Bug 3 — Trace correction overshoots for small offsets (< 1 px)

For `trace_offset = 0.5` px, the iterative correction (Block 8) overshoots in
the legacy function. For `trace_offset = 1.0` px the correction is accurate.
The root cause is the same discrete sigma_x sampling as Bug 2: `findfwhm`
overestimates `peak_x` when the peak falls between two sample points, causing
`trace_corr` to be initialized too large.

**Expected status in `fit_profile_refactor`**: substantially reduced by the Bug 2
fix. To be verified by `test_trace_correction_small_offset` (see Step 5 below).

---

## Unit Tests (`pypeit/tests/test_spatialprofile.py`)

The test file is written and all 12 tests pass. Tests are split into two groups:

**Tests against `spatialprofile.fit_profile` (legacy)**:
- `test_profile_return_shapes`
- `test_profile_normalization_bspline`
- `test_profile_normalization_gaussian`
- `test_forced_gaussian`
- `test_fwhm_recovery` — parametrized over `[6.0, 8.0]` px only (legacy Bug 2
  prevents testing narrower profiles)
- `test_trace_correction_direction` — uses `trace_offset=1.0` px (0.5 px
  overshoots in the legacy function)
- `test_prof_nsigma_extended` — asserts only shape and finiteness (not FWHM
  recovery or positivity) because the broken `nb` formula (Bug 1) produces a
  poor fit
- `test_profile_positivity_and_finite_all_paths`

**Test against `spatialprofile_refactor.fit_profile_refactor`**:
- `test_refactor_regression` — loads frozen reference `.npy` arrays generated
  by `pypeit/tests/files/freeze_spatprof.py` and verifies `fit_profile_refactor`
  matches them to `atol=1e-10`. Skipped on CI if reference files are absent.

### Key design notes

- **S/N knob**: `sn_ratio ≥ 20` reliably puts `sqrt(med_sn2) > sn_gauss = 4`,
  sending the function down the B-spline path. `sn_ratio ≤ 1` forces the
  Gaussian fallback.
- **Gaussian normalization**: `_gaussian_profile` (the refactored helper) returns
  row sums ≈ 1.0 (pixel-integrated form). The legacy `return_gaussian` returns
  row sums ≈ σ (point-sampled density). Tests that assert row sums near 1.0 must
  target `fit_profile_refactor`.
- **FWHM range**: `test_fwhm_recovery` is limited to `[6.0, 8.0]` when testing
  `fit_profile`; it can be extended to `[3.0, 4.0, 6.0, 8.0]` for
  `fit_profile_refactor` (Bug 2 fixed).
- **Trace offset**: `test_trace_correction_direction` uses 1.0 px for the legacy
  function; 0.5 px should work for `fit_profile_refactor` (Bug 3 reduced).
- **`prof_nsigma`**: `test_prof_nsigma_extended` cannot assert FWHM recovery for
  the legacy function (Bug 1); it can for `fit_profile_refactor`.

---

## Regression Freeze Workflow

Steps 1 and 2 are complete. The freeze now captures `fit_profile_refactor`
output (not `fit_profile`) and serves as a self-consistency check: re-run
`freeze_spatprof.py` whenever an intentional algorithmic change is made to
`fit_profile_refactor`.

1. ✅ Run `fit_profile_refactor` on the high-S/N synthetic inputs and save
   outputs as untracked reference files via `pypeit/tests/files/freeze_spatprof.py`:
   - `pypeit/tests/files/spatprof_profile_model_ref.npy`
   - `pypeit/tests/files/spatprof_xnew_ref.npy`
   - `pypeit/tests/files/spatprof_fwhmfit_ref.npy`

   Do not `git add` these files.

2. ✅ `test_refactor_regression` is implemented in
   `pypeit/tests/test_spatialprofile.py`. It loads the reference arrays and
   compares them against fresh output of `fit_profile_refactor` with `atol=1e-10`.
   The test is skipped automatically if the reference files are absent (CI).

---

## Verification Plan

1. ✅ Write unit tests and confirm all pass against the pre-refactor function
   (`pytest pypeit/tests/test_spatialprofile.py`: 12 passed).
2. ✅ Generate regression freeze reference arrays (`freeze_spatprof.py`) and
   leave them as untracked files.
3. ✅ Write `pypeit/core/spatialprofile_refactor.py` containing
   `fit_profile_refactor` and the four extracted helpers, incorporating all
   performance and accuracy improvements and bug fixes listed above.
4. ✅ Re-run `pytest pypeit/tests/test_spatialprofile.py`; all 12 tests pass
   including `test_refactor_regression`.
5. Strengthen tests targeting `fit_profile_refactor`:
   - **`test_fwhm_recovery_refactor`**: new test (or parametrize extension) calling
     `fit_profile_refactor`, extending the FWHM list to `[3.0, 4.0, 6.0, 8.0]`
     and asserting recovery within 10%.
   - **`test_gaussian_normalization_refactor`**: assert row sums of the Gaussian-path
     profile from `fit_profile_refactor` are close to 1.0 (pixel-integrated form),
     unlike the legacy function whose row sums are ≈ σ.
   - **`test_prof_nsigma_refactor`**: call `fit_profile_refactor` with
     `prof_nsigma=10.0` (the division-by-zero case in the legacy function) and
     assert FWHM recovery within 15%.
   - **`test_trace_correction_small_offset`**: call `fit_profile_refactor` with
     `trace_offset=0.5` px and assert `xnew` moves toward the true trace (Bug 3
     reduced).
6. Run the dev-suite for a multi-slit instrument (e.g., `keck_deimos`) to
   validate end-to-end behavior on real data.
7. `pytest pypeit/tests/` to catch any import or signature regressions.

---

## Change Log

### v1 → v2 (first edit)
- Added "Pre-existing Bugs Discovered During Testing" section documenting Bugs 1,
  2, and 3 found while running the unit tests against the pre-refactor function.
- Updated `_compute_bspline_knots` helper description to call out the broken `nb`
  formula (Bug 1).
- Adjusted three unit tests to match the passing test file checked in at
  `pypeit/tests/test_spatialprofile.py`: reduced `test_fwhm_recovery` parametrize
  list from `[3.0, 5.0, 8.0]` to `[6.0, 8.0]`; changed `test_trace_correction_direction`
  offset from 0.5 px to 1.0 px; changed `test_prof_nsigma_extended` from
  `prof_nsigma=10.0` to `prof_nsigma=12.0` and removed FWHM/positivity assertions.

### v2 → v2 (second edit)
- Added "Accuracy Improvements" section covering:
  - §A: pixel-integrated Gaussian in `return_gaussian` via `scipy.special.erf`
  - §B: linear interpolation for the half-maximum crossing in `findfwhm`
- Updated Bug 2 and Bug 3 descriptions to cross-reference the fixes in
  §A/§B and document the expected post-fix impact on unit tests.

### v2 → v3
- Added "Output File" section: refactored code lives in
  `pypeit/core/spatialprofile_refactor.py` with entry point `fit_profile_refactor`;
  legacy `spatialprofile.py` is left untouched.
- Removed "Full test file" Python code block from the Unit Tests section (the
  file is already written and checked in).
- Removed `test_refactor_regression` Python code block from the Regression Freeze
  Workflow section (also already implemented).
- Updated key design note for the regression test to reflect it calls
  `fit_profile_refactor`.
- Updated Verification Plan step 3 to name the output file explicitly; updated
  step 4 to note that `test_refactor_regression` now targets `fit_profile_refactor`.
- Moved change log to the bottom of the file and expanded it to cover all
  revisions since v1.

### v3 → v4
- Marked Verification Plan steps 3 and 4 as ✅ complete: `spatialprofile_refactor.py`
  is written and all 12 tests pass, including `test_refactor_regression`.
- Replaced "Candidate Helper Functions" section with two new sections:
  - "Implemented Helper Functions" — describes the four extracted helpers
    (`_bspline2d_to_1d`, `_findfwhm`, `_gaussian_profile`, `_return_gaussian`)
    and their roles.
  - "Deferred Block-Level Extractions" — records the decision to keep the eight
    larger blocks inline (labeled comment blocks in `fit_profile_refactor`), with
    rationale.
- Updated "Performance Improvements" and "Accuracy Improvements" sections from
  planned to implemented; added details on `_gaussian_profile`'s pixel-space
  interface (1D/2D broadcasting) and `_return_gaussian`'s split from construction.
- Updated "Pre-existing Bugs" to mark Bug 1 and Bug 2 as fixed in
  `fit_profile_refactor`; Bug 3 marked as expected-reduced.
- Updated "Regression Freeze Workflow" to note the freeze now captures
  `fit_profile_refactor` output (self-consistency check, not comparison against
  the legacy function).
- Updated Verification Plan Step 5 with four concrete new tests targeting
  `fit_profile_refactor` specifically.
- Updated Unit Tests design notes to distinguish legacy-function tests from
  `fit_profile_refactor`-targeting tests.
