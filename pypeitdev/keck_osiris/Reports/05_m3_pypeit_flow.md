# 05 - Milestone M3: rectification inside PypeIt and spec1d output

Date: 2026-10-11.  Author: Claude (Opus 5.5) with JXP.
Task: Prompt 1 of `claude_prompts/keck_osiris/M3_pypeit_flow.md` (Stage 3 of
`keck_osiris_implementation.md`).  Environment: laptop (Apple M5 Max),
conda env `pypeit14`, PypeIt branch `osiris_port`, dev-suite branch
`dev_osiris`.  No git operations were performed.

## 1. Summary

| Acceptance item | Target | Measured | Status |
|---|---|---|---|
| `run_pypeit keck_osiris_hn3_100.pypeit` | completes | exit 0, 2 min 45 s | pass |
| spec1d spectra | ~1000, one per lenslet | **3040** spectra = 3040 distinct lenslets (Hn3 is narrowband: 1024 slices x 3 lenslets, section 3.2) | pass (see deviation 1) |
| `MASKDEF_ID` | lenslet id | `row * 51 + col` of the Keck 66 x 51 grid (35-3325, all distinct); `MASKDEF_OBJNAME = L<col>_<row>` | pass |
| Vacuum wavelengths | flagged in header | `AIRORVAC = 'vac'`, `WAVESRC = 'Keck coeffs'`; OH lines match the vacuum `OH_NIRES` list to -0.01 A (air would be off by 4.4 A) | pass |
| Quasar position, collapsed spec1d vs DRP cube | match | peak (row 30, col 23) in both; centroid (30.139, 22.898) vs (30.146, 22.901) | pass |
| OH-line tilt over a block | < 0.05 px | **0.073 px** median (0.03-0.15 px, systematic by detector region) | **fail** (section 5.2; Q&A 1) |
| Per-lenslet OH centroid scatter after the Keck shifts | < 0.1 px | **0.034 px** median | pass |
| Weighted `lstsq` quasar spaxel vs `drp` (task 8) | < 1 % | **+0.03 %** (median spaxel 0.09 %) | pass; default switched to `lstsq` |
| Existing unit tests | pass | `pypeit/tests`: **889 passed** (`test_pkgdata.py` excluded, section 6) | pass |

## 2. What was built

### 2.1 PypeIt repo (`/Users/xavier/Projects/PypeIt/PypeIt`)

OSIRIS-specific code:

| File | Change |
|---|---|
| `pypeit/spectrographs/keck_osiris.py` | `has_rectification = True`; `rectify` (section 3.1); lenslet bookkeeping (`set_lenslet_state`, `lenslet_layout`, `lenslet_wavelengths`, `frame_expansion`, `rect_matrix`, `unorient_image`); Fiber hooks `get_fiber_blocks`, `get_block_slit_edges`, `identify_fibers_in_block`, `get_fiber_metadata`, `load_sky_layout`, `get_science_fiber_layout_indices`, `ifu_sky_wcs`, `get_ifu_datacube_meta`, `get_wcs`, `get_datacube_bins`; `get_fixed_tilts`; `get_fixed_wavecalib`; `get_fiber_wavelengths`; `subheader_for_spec` (vacuum flag); new defaults (`wavelengths.method = 'fixed'`, L.A.Cosmic `sigclip = objlim = 5`, `trim_edge = [0, 0]`, `no_poly`, `slitedges.pad = 0`); `bpm` tolerates a `DQ` of another shape |
| `pypeit/core/lenslet_rectify.py` | least-squares core rewritten (`_lstsq_core`, `_normal_blocks`, `_shift_blocks`) with the **sandwich variance**; `rectify_lstsq(..., var=)`, `lstsq_variance`, `background_variance` (flux-independent weights); `cosmic_mask_drp2016` (moved from the M2 vet test); Keck wavelength helpers `keck_order`, `keck_pixels`, `keck_wavelengths`; lenslet image `lenslet_layout`, `block_reference`, `lenslet_image` |
| `pypeit/scripts/osiris_rectify.py` | flux-independent weights + sandwich variance for `lstsq`; `drp` variance from `lstsq_variance` |
| `pypeit/par/pypeitpar.py` | `WavelengthSolutionPar.method` option **`fixed`**; new `CalibrationsPar.osiris_rect_crmask` (`lacosmic` / `drp2016` / `none`); `osiris_rect_method` default **`lstsq`**; descriptions updated (`osiris_rectmat_dir` expands environment variables and falls back to `<raw dir>/calib`) |
| `doc/pypeit_par.rst`, `doc/help/pypeit_osiris_rectify.rst` | regenerated (`build_par_rst.py`, `write_script_help.py`) |

Shared PypeIt code (all generic and gated; other spectrographs unchanged):

| File | Change | Gate |
|---|---|---|
| `pypeit/spectrographs/spectrograph.py` | base attributes/hooks: `has_rectification = False`, `runtime_par`/`set_runtime_par`, `rectify` (NotImplementedError), `get_fixed_tilts` (returns None), `get_fixed_wavecalib` (raises), `get_fiber_wavelengths` (returns None) | defaults are no-ops |
| `pypeit/images/rawimage.py` | new step `rectify` (after `subtract_dark`, before mosaic/flexure/flat/CR); after it, `build_ivar` adds no shot noise or noise floor (already propagated), the PypeItImage gets `shot_noise = False`, and saturation is not flagged on the rectified image | `spectrograph.has_rectification` / `steps['rectify']` |
| `pypeit/images/combineimage.py` | same shot-noise/noise-floor/saturation rule for stacks of rectified images | `'rectify' in process_steps` |
| `pypeit/pypeit.py` | `spectrograph.set_runtime_par(self.par)` (the rectify step needs the `CalibrationsPar`) | stores only |
| `pypeit/calibrations.py` | `get_tilts`: uses `spectrograph.get_fixed_tilts(slits)` when not None (no tilt frame traced); `get_wv_calib`: new `get_fixed_wv_calib` branch for `method = 'fixed'` (no arc fitting, no `msarc` required) | hook returns None / method |
| `pypeit/find_objects.py` | `FiberFindObjects`: per-fiber wavelengths from `get_fiber_wavelengths` cached as offsets from the block solution; no fiber sky model for an A-B reduction when `joint_fit = False`; `global_skysub`/`create_skymask` overrides return "no model / no mask" in that case; optional per-block `box_radius` | hook / `bkg_redux and not joint_fit` / key present |
| `pypeit/extraction.py` | `FiberExtract`: applies the cached per-fiber offsets; boxcar extraction on the column window of each aperture (identical result, ~50x faster; unit-tested) | offsets present / always (pure speed-up) |
| `pypeit/specobj.py` | internal `wave_fiber_offset_AA` | - |

New unit tests: `pypeit/tests/test_keck_osiris.py` (+7: inversion round trip on the
shipped coefficients, broadband and narrowband layouts, fiber hooks,
identity tilts, fixed WaveCalib + per-lenslet wavelengths, generic `rectify`
step with a stand-in transform), `pypeit/tests/test_lenslet_rectify.py` (+5:
sandwich variance incl. Monte Carlo, background variance, DRP strip rule,
synthetic inversion, lenslet image), `pypeit/tests/test_fiber_extract_window.py`
(new, 4: window boxcar == full-image boxcar).

### 2.2 Dev suite

| File | Change |
|---|---|
| `vet_tests/test_keck_osiris_flow.py` | new: OH-line tilt test (scatter, strict tilt as `xfail`, regression guard), vacuum check, spec1d lenslet test, quasar-position test |
| `vet_tests/test_keck_osiris_rectify.py` | `lstsq` vs `drp` quasar now asserted < 1 % (task 8); strip rule imported from PypeIt |
| `pypeit_files/keck_osiris_hn3_100.pypeit` | `[calibrations] osiris_rectmat_dir = ${PYPEIT_DEV}/RAW_DATA/keck_osiris/Hn3_100/calib` |
| `pypeitdev/keck_osiris/Reports/figures/make_m3_figures.py` | `m3_oh_tilts.png`, `m3_quasar.png`, `m3_lenslet_image.png` |

## 3. Design of the flow

### 3.1 The `rectify` step

`RawImage.process` -> ... -> `orient` -> `build_dark`/`subtract_dark` ->
**`rectify`** -> (flexure, flat: off) -> `build_ivar` -> PypeItImage
(`mask_cr` stays off).  `RawImage.rectify` passes the oriented image, its
full inverse variance (read noise + dark + processing + Poisson + noise
floor), the flux-independent variance, the read-noise variance, the BPM, the
headers, the `ProcessImagesPar` and the file name to
`KeckOSIRISSpectrograph.rectify`, which:

1. un-orients the arrays (inverse of `orient_image`: flip spectral axis,
   transpose);
2. subtracts the per-column median of rows 2032-2044
   (`subtract_reference_rows`) from the frame;
3. flags raw saturation (counts >= 0.8 x 65535 ADU) and cosmic rays
   (`osiris_rect_crmask = lacosmic`: PypeIt's L.A.Cosmic with `sigclip = 5`,
   `objlim = 5`) into the raw BPM;
4. selects (latest matrix of the filter/scale on or before the date) and
   caches the matrix;
5. solves with `osiris_rect_method` (`lstsq` default, `drp` available);
   weights = 1 / (read noise + dark + processing variance + per-column
   median background);
6. computes the variance with the least-squares machinery for **either**
   method: sandwich `(A^T W A)^-1 A^T W V W A (A^T W A)^-1` with the true
   variance `V`, and separately with `V` = read noise only;
7. maps the slices to the **lenslet image** (3.2) and flips the spectral axis.

It returns `image`, `rn2img` (read-noise part), `proc_var` (the rest) and
`bpm`; the step sets `datasec_img = 1`, `dark = 0`.  13 s per frame
(L.A.Cosmic 3.4 s, solve + two variances 7 s).

Because PypeIt processes A and B separately and subtracts the processed
images (`bkg_id`), the reference rows and cosmic rays are handled on each
raw frame rather than on the raw A-B difference; the solvers are linear, so
`rect(A) - rect(B)` is the rectified difference up to the masks.
L.A.Cosmic on a single raw frame flags ~4000 pixels, only 58 of them on the
brightest 0.5 % (OH) pixels of the sky frame (0.25 % of those pixels); the
2016-DRP strip rule on the same single frame flags 2239 OH pixels, so the
fallback was not needed (it is available as `osiris_rect_crmask =
drp2016`).

### 3.2 Lenslet image and blocks

A narrowband slice holds three (sometimes four) short spectra of different
lenslets head to tail, so "one row per slice" cannot carry one wavelength
solution per spectrum.  The rectified frame is therefore re-laid out as a
**lenslet image**: one spatial pixel per lenslet, `pos = 64 * bidx +
slice % 64`, where `bidx` indexes the Keck lenslet columns present for the
filter.  Each lenslet's spectrum is copied from its slice over a window of
raw columns covering the filter band (DRP grid) extended halfway into the
gaps to its neighbours (464-474 px for Hn3); everything else is masked.
Broadband filters give the identity (`pos = slice`, 19 blocks, 1216
pixels); Hn3 gives 51 blocks (3264 pixels) and 3040 lenslets with Keck
coefficients (the DRP extras clipped onto row 65 are dropped).

Blocks = block-slits with fixed edges `64 b - 0.5 ... 64 (b+1) - 0.5`
(`get_block_slit_edges`); one `SpecObj` per lenslet, `BOX_R_PIX = 0.5`
(one pixel; per-block `box_radius`).  The filter and date needed by the
hooks (which receive no header) are recorded on the spectrograph instance by
`rectify`/`get_fixed_wavecalib` (`set_lenslet_state`, the same mechanism as
`dispname` for echelle spectrographs).

### 3.3 Tilts and wavelengths

- `get_fixed_tilts`: identity (`tilt = spec / (nspec - 1)`; Legendre
  coefficients 0.5, 0.5 and zero spatial terms); written as `Tilts_*.fits`.
- `get_fixed_wavecalib`: one `WaveFit` per block = the Keck solution of the
  block's reference lenslet (row 32), obtained by inverting `x(lambda) =
  poly(lambda * order / 3 / frac - 2200)` on a 0.05 nm grid (`frac` from
  `DTMP7`; 1 with a warning when missing), fit with a 7th-order Legendre
  polynomial (rms 1e-7 px); constant FWHM map (2 px).
- **Finding:** within a block the lenslets are staggered by **1.92 px per
  lenslet row** (~120 px over the block) and their dispersions differ (a
  pure shift leaves up to 0.6 px).  The "small per-lenslet shifts" of the
  design (Q&A 21) do not exist, so `get_fiber_wavelengths` returns the
  exact Keck solution of every lenslet; `FiberFindObjects` stores it as a
  per-pixel offset from the block solution and `FiberExtract` adds it to
  `BOX_WAVE`/`OPT_WAVE` (verified: spec1d wavelengths equal the Keck
  per-lenslet solution to 0.17 A, which is the heliocentric correction,
  factor 1.0000102).  The spec2d wavelength image is exact only for the
  reference lenslets (Q&A 2).

### 3.4 Sky and extraction (task 5)

With the A-B pair no fiber sky model is fit (`joint_fit = False`), the
slit-by-slit sky model of the un-subtracted frame is skipped, and
`BOX_COUNTS_SKY = 0`.  Without a flat, `FiberExtract` only does the boxcar
(no optimal extraction; `_refine_fiber_traces` is a no-op).  Multi-pixel
assumptions found and fixed:

1. `BOX_R_PIX` from the inter-fiber spacing (wider apertures at gaps and
   block ends) -> optional per-block `box_radius` (0.5 for OSIRIS);
2. the bkg-redux slit-by-slit sky fit masked every fiber and **rejected all
   51 slits, so no spec1d was written** -> `FiberFindObjects.global_skysub`
   returns a zero model in that configuration;
3. `create_skymask` over 3040 objects took 167 s for a mask that is not
   used -> skipped in that configuration;
4. `FiberExtract` ran full-image operations per fiber (3040 x 80 ms) ->
   boxcar on the aperture's column window.  Run time 9 min 48 s -> 2 min
   45 s (profile: extraction 249 s -> ~5 s, skymask 167 s -> 0).

One lenslet (of 3040) has no good pixel (its window is fully masked) and an
all-masked spectrum.

## 4. Task 8: production `lstsq` weighting

Weights = flux-independent (read noise + per-column median background);
variance = sandwich with the true pixel variance (diagonal; the
off-diagonal blocks of `N^-1` are built with the forward/backward Schur
complements, `G_{r,s} = -S_r^-1 U_r G_{r+1,s}` and `G_{r,s} = -T_r^-1
U_{r-1}^T G_{r-1,s}`, within +/-2 bands; the result is identical to +/-6
bands at float32 precision).  Checks on the Hn3 frame (columns 1000-1031,
200 noise realisations): Monte-Carlo variance / sandwich = 0.99 (5-95 %:
0.84-1.16); `diag(N^-1)` would be wrong by 5.5 %.

`vet_tests/test_keck_osiris_rectify.py` (5 passed):

```
drp:   3024 spaxels with S/N > 10; median |frac| = 0.0074; 90th pct = 0.0640; quasar ratio = 1.0004
lstsq: 3024 spaxels with S/N > 10; median |frac| = 0.0077; 90th pct = 0.0660; quasar ratio = 1.0008
lstsq vs drp: median |frac| = 0.0009; quasar ratio = 1.0003
```

The quasar spaxel agrees with `drp` to **0.03 %** (M2: -2.1 % with Poisson
weights), so `osiris_rect_method` now defaults to `lstsq` (Q&A 5).

Note on the reference rows: the script (comparison with the DRP) subtracts
them from the A-B difference, as the DRP did; the PypeIt flow subtracts them
from each frame (M2 Q&A 4).  The two are statistically equivalent (the
per-column difference is 0.66 ADU rms, the noise of a 13-row median), but
the per-frame order is a different noise realisation from the DRP's and
raises the median spaxel difference to the DRP cube from 0.74 % to 2.0 %
(S/N 10-30 spaxels).

## 5. Verification

### 5.1 run_pypeit

```
$ cd <scratch>/REDUX_OUT/keck_osiris/Hn3_100
$ cp $PYPEIT_DEV/pypeit_files/keck_osiris_hn3_100.pypeit .
$ PYPEIT_DEV=... run_pypeit keck_osiris_hn3_100.pypeit -o
...  Raw-frame masks: 0 saturated and 4068 cosmic-ray pixels (lacosmic).
...  Rectification of s160321_a002010.fits took 13.3 s; 3040 lenslets.
...  Keck wavelength solution: 51 blocks, max inversion-fit rms 9.58e-08 px.
...  Background frame subtracted and joint_fit is False: no fiber sky model is fit.
...  Applied per-fiber wavelength offsets to 3040 fibers
...  Execution time: 2m 44.01s          (exit 0)
Science/spec1d_s160321_a002010-3C273_OSIRIS_20160321T101742.422.fits
Science/spec2d_s160321_a002010-3C273_OSIRIS_20160321T101742.422.fits
```

Calibrations written: `Arc`, `Tiltimg`, `Slits` (51 block-slits), `Tilts`
(identity), `WaveCalib` (fixed).

spec1d check (`SpecObjs.from_fitsfile`):

```
nobj 3040, unique MASKDEF_ID 3040, min/max id 35 3325
AIRORVAC vac | WAVESRC Keck coeffs | PYPELINE Fiber
good pixels per lenslet: median 466, 1 lenslet with none
brightest lenslet 1808 = L23_35 (row 35 -> flipped row 30, col 23)
```

### 5.2 OH-line tilt test (design 8, test 3)

`vet_tests/test_keck_osiris_flow.py` processes the sky frame 011 through
`RawImage` (rectification included), measures the OH centroids of every
lenslet with `arc.detect_lines` (Gaussian fits), predicts each line in every
lenslet of a block from its Keck solution, and fits the residuals versus the
lenslet row:

```
OH tilt test: 649 (block, line) pairs in 47 blocks; median |tilt| = 0.073 px,
90th pct = 0.271 px; median scatter = 0.034 px, 90th pct = 0.158 px;
median offset = -0.025 px
```

The scatter passes; the tilt does not.  Diagnosis (`figures/m3_oh_tilts.png`):

| Run | median abs. tilt | signed median tilt, blocks 3-18 / 19-34 / 35-50 | scatter |
|---|---|---|---|
| frame 011, `lstsq` | 0.073 | +0.133 / +0.034 / +0.072 | 0.034 |
| frame 011, `drp` | 0.066 | +0.129 / +0.029 / +0.070 | 0.032 |
| frame 011, `lstsq`, no `DTMP7` term | 0.073 | +0.133 / +0.034 / +0.072 | 0.034 |
| frame 010, `lstsq` | 0.080 | +0.149 / +0.036 / +0.071 | 0.032 |

The tilt is a smooth, signed function of the detector region (largest for
the left third of the detector), identical for both solvers and both
frames, with a 0.03 px per-lenslet scatter and no saw-tooth.  The `lstsq`
solver does not mix detector columns, so line centroids along a rectified
row are the raw detector positions of that lenslet's light: the residual
trend (0.002 px per lenslet row, i.e. 0.0005 nm) is in the Keck per-lenslet
coefficients, not in the rectification.  A slice/lenslet mis-mapping or a
`hilo` error would give ~2-px offsets (one lenslet row) or a saw-tooth
(the regression test guards against both: median |tilt| < 0.2 px, offset <
0.1 px).  The strict 0.05-px test is kept as `xfail` with this reason
(Q&A 1).

### 5.3 Vacuum wavelengths and the `DTMP7` term

Bright OH lines of the block reference lenslets of the sky frame vs the
vacuum `OH_NIRES` list (625 isolated lines):

| Grating-expansion term | median (data - list) | MAD |
|---|---|---|
| `DTMP7 = 55.34 K` applied (frame 011, `frac = 1.000025`) | **+0.398 A** | 0.097 A |
| none (`frac = 1`) | **-0.009 A** | 0.096 A |

The wavelengths are vacuum (air would be off by 4.4 A).  The DRP
temperature term, implemented as specified, makes the 2016 solution worse
by 0.4 A (0.19 px) for this frame (Q&A 4).  The science frame 010 has no
`DTMP7`, so its spec1d is uncorrected (warning in the log).

### 5.4 Quasar position (collapsed spec1d vs DRP cube)

Boxcar counts summed over 15960-16740 A per lenslet, Keck I row flip,
compared with the DRP cube summed over the same range
(`figures/m3_quasar.png`):

```
3C 273: PypeIt peak (30, 23) centroid [30.139 22.898];
        DRP    peak (30, 23) centroid [30.146 22.901]
```

The quasar spectrum agrees in shape and wavelength (absorption features
aligned once the DRP's 0 K `DTMP7` scale is applied to its grid); PypeIt is
~3 % high at the red end, as in M2 (the DRP cube has its "Correct
Dispersion" shift and 2016 cosmic-ray input).

### 5.5 Tests

```
$ pytest pypeit/tests --ignore=pypeit/tests/test_pkgdata.py
889 passed in 85.87s
$ pytest pypeit/tests/test_keck_osiris.py         20 passed
$ pytest pypeit/tests/test_lenslet_rectify.py     18 passed
$ pytest pypeit/tests/test_fiber_extract_window.py 4 passed
$ PYPEIT_DEV=... pytest vet_tests/test_keck_osiris_rectify.py
5 passed
$ PYPEIT_DEV=... pytest vet_tests/test_keck_osiris_flow.py --redux_out <scratch>/REDUX_OUT
5 passed, 1 xfailed            (xfail = strict 0.05-px tilt, section 5.2)
$ PYPEIT_DEV=... pytest unit_tests/test_load_images.py -k keck_osiris
1 passed
```

`test_pkgdata.py` was excluded: it downloads from the GitHub branch
`osiris_port`, which does not exist on GitHub, and one of its tests deletes
the tracked file `pypeit/data/tests/cache_test.txt` (Reports 03, 04).

## 6. Deviations from the design document

1. **One spectrum per lenslet, not per slice**: the narrowband Hn3 spec1d
   has 3040 spectra (the prompt expected ~1000; broadband gives up to
   1216).  The lenslet image (3.2) replaces "slice index + 0.5" positions.
2. **Per-lenslet wavelengths are full solutions**, not scalar shifts (3.3);
   the new generic hook is `get_fiber_wavelengths` (absolute wavelengths,
   converted to offsets by the Fiber pypeline).
3. **Tilts are the identity as specified, but OH lines are tilted by 1.92
   px per lenslet row in the rectified frame**; the spec2d wavelength image
   is exact only for the block reference lenslets (Q&A 2).
4. **Cosmic rays and reference rows per raw frame**, not on the raw A-B
   difference (3.1).
5. **`osiris_rect_method` default is now `lstsq`** (task 8 passed).
6. **No fiber sky model with A-B** (`joint_fit = False`; Davies OH scaling
   is M4); spec1d has boxcar spectra only (no flat).
7. The fixed WaveCalib takes the header of the first arc (else trace) frame
   of the calibration group; per-lenslet wavelengths use the science frame's
   `DTMP7`.
8. `osiris_rectmat_dir` expands environment variables and falls back to
   `<raw dir>/calib` and `<raw dir>`.
9. `pypeit/wavetilts.py` is imported inside `get_fixed_tilts` (module-level
   import is circular: `wavetilts -> slittrace -> specobj -> flexure ->
   onespec -> spectrograph registry`).
10. The spec1d wavelengths are heliocentric (PypeIt default `refframe`).

## 7. Not done / next (M4)

- OH refinement of the Keck solution (would absorb the 0.03-0.15 px block
  tilts and the `DTMP7` offset), Davies OH-scaled sky subtraction, the
  generic `pypeit_fiber_cube` (the hooks it needs exist but the cube was not
  built), Kbb/050 (needs the Keck matrix).
- The dev-suite `reduce` test of `keck_osiris/Hn3_100` was not run through
  `pypeit_test`; the equivalent `run_pypeit` was run by hand (5.1).
- Documentation and release notes (M5).
