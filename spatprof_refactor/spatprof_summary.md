# `spatprof_refactor` Branch — PR Summary

## Overview

This branch refactors `pypeit/core/spatialprofile.py`, which contains the
spatial profile fitting code used for optimal extraction.  The module is a
direct port of the IDL LOWREDUX routine `long_gprofile.pro` that had never been
touched after its initial translation.  The goals were: (1) make the code
readable and maintainable; (2) fix pre-existing bugs; (3) improve numerical
accuracy; (4) wire diagnostic QA output into the standard `run_pypeit`
pipeline; and (5) document the algorithm.

---

## Files changed

| File | Change |
|------|--------|
| `pypeit/core/spatialprofile.py` | Full rewrite (primary change) |
| `pypeit/core/skysub.py` | Updated call site; QA wiring; `spec_img` |
| `pypeit/extraction.py` | Updated QA path parameters |
| `pypeit/qa.py` | New `spat_profile_qa` filename case |
| `pypeit/utils.py` | `growth_lim` extension |
| `pypeit/tests/test_spatialprofile.py` | New — comprehensive unit test suite |
| `doc/spatprof.rst` | New — algorithm documentation page |
| `doc/index.rst` | Registered new page in toctree |
| `doc/skysub.rst` | Cross-reference to new page |
| `doc/scripts/make_spatprof_figures.py` | New — figure generation script |
| `doc/figures/spatprof_example_bspline.png` | New — example QA figure |
| `doc/figures/spatprof_example_gaussian.png` | New — example QA figure |

---

## `spatialprofile.py` — rewrite

### Structure

The 593-line monolith `fit_profile` has been restructured around seven
module-private helpers, each covering a coherent sub-task.  The helpers, in
calling order:

| Helper | Encapsulates |
|--------|-------------|
| `_findfwhm` | Sub-pixel half-maximum crossing (used in Blocks 7 and 11) |
| `_gaussian_profile` | Pixel-integrated Gaussian profile (replaces point-sampled form) |
| `_fit_spectrum_and_normalize` | Spectral B-spline fitting, S/N estimation, normalized-object construction (Blocks 2–5) |
| `_profile_coordinates_and_model_sampling` | B-spline knot setup (Block 6) |
| `_apodize_profile` | Apodization limit search and exponential tail fitting (Blocks 11–12) |
| `_bin_profile` | Binned percentile statistics for QA |
| `_fit_profile_qa` | Full diagnostic figure |

The `fit_profile` body itself follows the original thirteen-block structure for
readability.

### Accuracy improvements

**Pixel-integrated Gaussian** (`_gaussian_profile`): the legacy code evaluated
the Gaussian density at pixel centers (`exp(-0.5*sigma_x²)/sqrt(2π)`).  The
new code integrates across each pixel using the error function, which matters
significantly for narrow profiles: the peak is underestimated by ~10% at
FWHM = 3 px with point-sampling; the integrated form reduces this to < 1%.
Row sums of the returned profile are ~1.0, matching the B-spline normalization
path.

**Sub-pixel FWHM** (`_findfwhm`): the legacy function returned the first
discrete sample below half-peak, systematically overestimating the FWHM by one
pixel for narrow profiles (33% error at FWHM = 3 px).  The new function
interpolates between the two bracketing samples to locate the half-maximum
crossing to sub-pixel precision (< 2% error for FWHM ≥ 3 px).

**`sig2fwhm` constant**: `2.3548` → `np.sqrt(8*np.log(2))` throughout.

### Bug fixes

**Bug 1 — `nb` formula (`prof_nsigma` branch)**: the legacy code computed
`nb = np.round(prof_nsigma > 10)` — a boolean cast that yielded 0 or 1
regardless of `prof_nsigma`, causing division-by-zero for `prof_nsigma ≤ 10`
and only two knot positions for larger values.  Fixed to
`nb = max(1, round(prof_nsigma / 10))`.

**Bug 2 — FWHM overestimation**: see `_findfwhm` above.

### Performance

The `np.outer(v, np.ones(n))` pattern (which allocated a temporary
`np.ones(nspat)` array before forming the outer product) was replaced with
`v[:, None]` broadcasting throughout.  The two `while True` loops in the
apodization limit walk were replaced with vectorized searches over the
pre-computed derivative arrays.  Benchmarks on three synthetic image sizes show
~7–20% speedup and ~7–12% reduction in peak heap allocation relative to the
legacy function.

### New `spec_img` parameter

`fit_profile` accepts an optional `spec_img` parameter: an integer image where
`spec_img[i, j] = i` (the spectral row index).  When provided, the function
uses it directly instead of calling `np.where(totmask)` to derive row indices
internally, avoiding a redundant full-image traversal.  The call site in
`skysub.py` now constructs `spec_img_sub` via `np.broadcast_to` and passes it
in explicitly.

---

## QA wiring

Spatial profile QA figures are now written automatically to `QA/PNGs/` for
each extracted object on the final profile-fitting iteration.  Earlier
iterations produce intermediate, unconverged fits and are not written.

**Filename convention** follows the standard `qa.set_qa_filename` factory:
```
QA/PNGs/{basename}_{det}_S{slit:04d}_O{obj:03d}_spat_prof.png
```
A new `spat_profile_qa` case was added to `qa.set_qa_filename`, along with an
`obj: int | None = None` keyword argument for cases where multiple objects
share a slit.

**Interface change**: the single `qa_path` parameter on
`local_skysub_extract` and `ech_local_skysub_extract` was replaced by three
parameters: `qa_outdir`, `qa_basename`, and `qa_detname`, following the pattern
established by `get_findobj_qa_filename` in `find_objects.py`.  Both extraction
wrappers in `extraction.py` were updated accordingly.

**Bug fix — trace coordinate offset**: the image panels in `_fit_profile_qa`
were displayed with `extent = [0, nspat, ...]` (sub-image local coordinates)
while the trace arrays `trace_in` and `xnew` are in absolute spatial pixel
coordinates.  For objects whose local-sky sub-image does not start at column 0,
the trace overlay was displaced from the flux peak.  Fixed by computing
`spat_offset = float(spat_img[0, 0])` in `fit_profile` and passing it to
`_fit_profile_qa`, which uses it to set `extent = [spat_offset, spat_offset + nspat, ...]`.

---

## `utils.py` — `growth_lim` extension

`growth_lim` accepts two new string values for the `midpoint` keyword:

- `'median'` (new default, preserves previous behavior): midpoint is the
  sample median.
- `'center'`: midpoint is the center of the empirical data range.  Used by
  `_fit_profile_qa` to set symmetric y-axis limits around zero in the spatial
  profile panel.

A `PypeItError` guard was added for invalid `midpoint` values.

---

## Unit tests (`pypeit/tests/test_spatialprofile.py`)

The module had no tests before this branch.  The new test file contains
44 tests covering:

- Return shapes and finiteness (both B-spline and Gaussian paths, and the
  early-return NaN-flux path)
- Row normalization for the B-spline path
- Profile positivity and bell-shape for the Gaussian path
- `gauss=True` forced Gaussian
- FWHM recovery to within 10% for three profile widths (parametrized)
- Trace correction direction and magnitude
- `prof_nsigma` extended-object mode
- `_findfwhm` sub-pixel accuracy and edge cases
- `_gaussian_profile` normalization and pixel-integration accuracy
- `_fit_spectrum_and_normalize` (success path, early-return path, S/N
  estimation)
- `_profile_coordinates_and_model_sampling` (standard and `prof_nsigma` modes)
- `_apodize_profile` (B-spline path and `no_deriv` path)
- Refactor regression against saved reference arrays (skipped automatically
  on CI if reference files are absent)
- Fourier-flux + curved trace end-to-end test

---

## Documentation (`doc/spatprof.rst`)

A new documentation page describes the spatial profile fitting algorithm,
covering: the S/N decision (B-spline vs. Gaussian), the spectral B-spline
fitting, the iterative trace and width correction loop, apodization, the
Gaussian fallback, profile normalization, QA output, and key parameters.  Two
example figures are included (B-spline path at S/N = 20 and Gaussian fallback
at S/N = 2), generated from the same synthetic data construction used in
`test_fourier_flux_and_curved_trace`.  The page is registered in the Sphinx
toctree and cross-referenced from `doc/skysub.rst`.
