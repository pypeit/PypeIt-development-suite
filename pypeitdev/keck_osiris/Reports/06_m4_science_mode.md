# 06 - Milestone M4: science-mode reduction of NGC 3245 (Kbb/050)

Date: 2026-10-12.  Author: Claude (Opus 5.5) with JXP.
Task: Prompt 1 of `claude_prompts/keck_osiris/M4_science_mode.md` (Stage 4 of
`keck_osiris_implementation.md`).  Environment: laptop, conda env
`pypeit14`, PypeIt branch `osiris_port`, dev-suite branch `dev_osiris`.
No git operations were performed.  All reductions ran under the session
scratchpad (`.../scratchpad/REDUX_DS`, `REDUX_HN3`), not in the repos.

## 1. Summary

| Acceptance item | Target | Measured | Status |
|---|---|---|---|
| Matrix selection by validity epoch (task 11) | 2024-02-03 selects `s240524_c006`; 2016 Hn3 selects `s160318_c005` | both, unit tests `test_matrix_epochs`, `test_select_matrix_epoch`, `test_matrix_from_gzip` | pass |
| NGC 3245 reduced end to end | completes | dev-suite `keck_osiris/Kbb_050` reduce: PASSED, 39 min 31 s for 50 reduced frames (15 NGC 3245 A-B, 31 standards, run concurrently with other jobs) | pass |
| Continuum image | physical | smooth bulge, peak at lenslet (row 33, col 23), well inside the field (`m4_ngc3245.png`) | pass |
| Emission-line map | physical | H2 1-0 S(1) concentrated on the nucleus (3x3 flux 5745 counts vs map noise 970; nuclear EW +0.39 A) but with a horizontal step at flipped lenslet row ~37 (section 5.4) | partial (Q&A 3) |
| OH residuals, Davies scaling vs plain A-B | reduced (quantify) | median |integrated OH-line residual| / sky-line flux: 600-s NGC 3245 frames **0.053 -> 0.037** (12 of 13 improve, 1 equal); sky-sky pair **0.082 -> 0.030**; Hn3 pair line-pixel rms **225 -> 62** counts.  The two 44-s frames get **worse** (0.040 -> 0.068, 0.074 -> 0.081) | pass on science exposures; short frames open (Q&A 2) |
| Telluric-corrected A0V | flat to a few percent | HD101060 (`OS.20240203.59407.32`): corrected/model **1.5 %** rms in 50-A bins (median 1.0 %, max 4.5 %), 4.9 % per pixel (`m4_telluric.png`) | pass |
| Dev-suite Kbb_050 reduce + afterburn | pass | reduce, `pypeit_fiber_cube` (NGC3245), `pypeit_fiber_cube` (HD101060 point source), `pypeit_sensfunc`, `pypeit_tellfit`: **5/5 PASSED** | pass |
| Kbb vet tests | pass | `vet_tests/test_keck_osiris_kbb.py`: **3 passed** (spec1d, cube structure, telluric flatness) | pass |
| OH-line tilt test, refined solution, no `xfail` | < 0.05 px | **Hn3 0.023 px**, **Kbb (ISSKY=1 frame `OS.20240203.43383.79`) 0.020 px** (Keck coefficients alone: 0.073 and 0.106 px) | pass |
| Keck-derived tilts vs per-lenslet Keck solutions | < 0.05 px rms | **Kbb 0.044 px**, **Hn3 0.048 px** rms (robust 0.021 / 0.022 px) | pass (close to the limit, Q&A 6) |
| `DTMP7` re-test (task 9) | decide | +0.40 A shift on both epochs; term now off by default (`osiris_grating_temp = False`) | done |
| Existing + new PypeIt unit tests | pass | `pypeit/tests`: **900 passed** (`test_pkgdata.py` excluded) | pass |

## 2. What was built

### 2.1 PypeIt repo (`/Users/xavier/Projects/PypeIt/PypeIt`)

OSIRIS-specific:

| File | Change |
|---|---|
| `pypeit/core/lenslet_rectify.py` | `RECTMAT_URL`, `load_matrix_epochs`, `matrix_epoch`; `select_matrix` by validity epoch (file date inside the epoch, latest-before fallback with a warning, error naming the epoch and the tkserver directory); `.fits.gz` accepted; `lenslet_band_mask`; `lenslet_block_tilts` (2D Legendre per block); broadband lenslet windows = filter band |
| `pypeit/data/spectrographs/keck_osiris/rectmat_epochs.ecsv` | new: the 17 validity epochs copied from http://tkserver.keck.hawaii.edu/osiris/ |
| `pypeit/data/arc_lines/lists/OH_OSIRIS_lines.dat` | new: the OSIRIS DRP K-band OH list (`ohlines_kbb_all.dat`, vacuum; NaN-intensity blend components flagged `NIST = 0`) |
| `pypeit/spectrographs/keck_osiris.py` | `get_fixed_tilts` returns Keck-derived tilts (`TILT_ORDER = (6, 4)`, uses the raw header); `frame_expansion(apply=None)` follows `osiris_grating_temp`; `fiber_cube_setup` (lenslet state from the spec1d header, scale = `SSCALE`, `wcs_ref = 'origin'`); block edges at -0.5/+0.5 px (first/last lenslet were dropped); defaults: `fiber_sky_refine = zeropoint`, `osiris_oh_scaling = global`, `OH_OSIRIS` list for K filters, `sensfunc`/`telluric` `polyorder = 7` |

Generic shared code (all gated on a parameter or a spectrograph hook; to be
flagged to the PR reviewers, Q&A 7):

| File | Change |
|---|---|
| `pypeit/par/pypeitpar.py` | new `[calibrations][wavelengths] fiber_sky_refine` (`zeropoint`/`block`), `[reduce][skysub] osiris_oh_scaling` (`global`/`block`/`fiber`), `[calibrations] osiris_grating_temp`; `doc/pypeit_par.rst` regenerated (verified current) |
| `pypeit/core/wavecal/wvutils.py` | `isolated_lines`, `fiber_sky_line_shifts`, `stacked_sky_line_shift`, `fiber_sky_refinement`, `fiber_shift_qa` |
| `pypeit/core/skysub.py` | `OH_GROUP_EDGES`/`OH_GROUP_NAMES` (Davies 2007 groups as in `scaledskysub_000.pro`), `running_median`, `oh_scaled_sky` |
| `pypeit/find_objects.py` | `FiberFindObjects`: `_refine_fiber_sky_lines` (when `fiber_sky_refine` is set), `_fiber_oh_scaling` (when `osiris_oh_scaling` is set, A-B frames only, not standards), `_fiber_column_spectra`, `_stamp_fiber_spectra`, `_fiber_sky_image`; new attribute `bkg_redux_sciimg` |
| `pypeit/pypeit_steps.py`, `pypeit/exposure.py`, `pypeit/scripts/reduce_by_step.py` | pass the un-subtracted A frame (`bkg_redux_sciimg`) to object finding |
| `pypeit/images/rawimage.py`, `pypeit/images/buildimage.py` | `process(rectify=None)`; bias and dark frames are not rectified (so darks are subtracted from raw frames before the rectification) |
| `pypeit/calibrations.py` | `get_tilts` passes the header of a raw frame of the calibration group to an overridden `get_fixed_tilts` |
| `pypeit/spectrographs/spectrograph.py` | `get_fixed_tilts(slits, raw_header=None)`; new hook `fiber_cube_setup(header)` (returns cube defaults) |
| `pypeit/core/datacube.py` | `wcs_ref` (`center`/`origin`), `fiber_extraction_prefix` (boxcar fallback), `extract_fiber_point_source` (single-object `_ps` spec1d for sensfunc/tellfit) |
| `pypeit/scripts/fiber_cube.py` | new generic `pypeit_fiber_cube` (`--point_source`, `--no_cube`); entry point in `pyproject.toml`, `__init__.py`, `doc/help/pypeit_fiber_cube.rst` |
| `pypeit/scripts/binospec_ifu_cube.py` | thin alias of `FiberCube` (Binospec only, 0.27" default; same options) |
| `pypeit/qa.py`, `pypeit/specobj.py` | QA names `fiber_sky_refine`/`fiber_oh_scaling`; internal attribute `wave_sky_shift_pix` (not in the datamodel) |

Tests: `pypeit/tests/test_lenslet_rectify.py` (epoch selection, gzip),
`test_keck_osiris.py` (defaults, Keck-derived tilts, **tilt-based
wavelength image vs per-lenslet Keck solutions** for Kbb 2024 and Hn3 2016,
dark not rectified, `fiber_cube_setup`, cube grid/WCS),
new `test_fiber_sky.py` (synthetic refinement and Davies scaling).

### 2.2 Dev suite

| File | Change |
|---|---|
| `pypeit_files/keck_osiris_kbb_050.pypeit` | `osiris_rectmat_dir`; darks in ITIME calibration groups 1-3; only the 600-s pair 42705/43383 carries `arc,tilt,trace` |
| `test_scripts/pypeit_tests.py`, `test_scripts/test_setups.py` | new `PypeItFiberCubeTest` (afterburn) and `_fiber_cube` setups (NGC 3245 cube; HD101060 point source), sensfunc and telluric entries for `Kbb_050` |
| `sensfunc_files/keck_osiris_kbb_050.sens`, `tellfit_files/keck_osiris_kbb_050.tell` | IR algorithm, A0, V = 8.78; `objmodel = star` |
| `vet_tests/test_keck_osiris_flow.py` | tilt test on Hn3 *and* Kbb on the OH-refined solution, `xfail` removed; refinement repeatability test |
| `vet_tests/test_keck_osiris_kbb.py` | new: spec1d (lenslets, vacuum, Davies model in the spec2d), cube structure/WCS, telluric flatness |
| `pypeitdev/keck_osiris/Reports/figures/make_m4_figures.py` | figures `m4_oh_tilts`, `m4_tilt_model`, `m4_oh_scaling`, `m4_ngc3245`, `m4_telluric` |

## 3. Design notes

### 3.1 Matrix selection (task 11, done first)

Keck builds one matrix per (filter, scale) and validity epoch; the file
date is the scan date and can be after the observations
(`s240524_c006___infl_Kbb_050` serves 2024-02-03, epoch 230701-250310).
`select_matrix`: (1) epoch containing MJD-OBS (nearest epoch with a warning
in a gap), (2) latest matrix of the filter/scale dated inside the epoch,
(3) otherwise latest-before with a warning, (4) otherwise an error that
lists the directory and names the tkserver subdirectory.  Recorded in design
Section 7 and in the M5 documentation task.

### 3.2 Keck-derived tilts (task 8)

For every block, the tilt of pixel (x, lenslet) is the pixel of the block's
reference lenslet at the same wavelength, `t = lambda_ref^-1(lambda_p(x))`,
fit with a 2D Legendre (6, 4) in PypeIt's normalisation, using the pixels
inside the filter band.  The spec2d wavelength image `wave_ref(tilts)`
reproduces the per-lenslet Keck solutions to 0.044 px rms (Kbb) and
0.048 px (Hn3), robust 0.021-0.022 px; p99 0.09/0.18 px; a handful of
lenslets with peculiar coefficients reach 1-4 px (`m4_tilt_model.png`).
The per-lenslet offsets of `get_fiber_wavelengths` now carry only this
residual, so the spec1d wavelengths stay exact.  Broadband lenslet windows
are restricted to the filter band (out-of-band pixels masked), which also
removes the filter-edge steps that confused line finding.

### 3.3 Darks and A-B (task 2)

Paired frames are A-B through `bkg_id` without darks (Q&A 23).  Darks are
processed unrectified (`buildimage` passes `rectify=False` for
bias/dark), so a dark of the same ITIME can be subtracted from an
unpaired raw frame before its rectification: on `OS.20240203.42705.45`
with the three 600-s darks the processing steps are `apply_gain, trim,
orient, subtract_dark, rectify` and the rectified dark contribution is
0.7 e- median (90th pct 18.6 e-).  The dev-suite file keeps the darks in
calibration groups 1-3 (600, 1.5, 8.85 s); this night has no unpaired
science frame, so they are not used.

### 3.4 OH refinement (task 4) and `DTMP7` (task 9)

`fiber_sky_refine = zeropoint`: the isolated lines (no listed neighbour,
blend components included, within 3 FWHM) of `OH_OSIRIS` (K) or
`OH_NIRES` (J/H) are centroided in every lenslet of the un-subtracted A
frame; the median offset is applied as a pixel zero point.  Lenslets with
< 3 lines get the frame median; frames with < 10 measured lenslets the
shift of the median-stacked spectrum; QA
`QA/PNGs/*_fiber_sky_refine.png`.

| Frame type | Result |
|---|---|
| NGC 3245 600 s | 1012-1013/1207 lenslets measured, 32 isolated lines (median 13 per lenslet), median shift +0.536 px, std 0.054 px, error 0.021 px |
| NGC 3245 44 s, HD101060 3/8.85 s | stack shift +0.53 to +0.63 px |
| HD89239 1.5 s | no lines: Keck solution kept (+1.5 A off; Q&A 4) |
| Hn3 2016 (dev-suite run) | 3039/3040 measured, median +0.003 px, std 0.060 px |

The per-lenslet shifts are reproducible between the science and the sky
frame (r = 0.971 Hn3, 0.898 Kbb; std of the difference 0.014 / 0.024 px),
i.e. structure of the Keck coefficients, not noise (`m4_oh_tilts.png`).

`DTMP7` (median data - list over all lenslets, `scratchpad/m4/dtmp7b.py`):

| Frame | `DTMP7` | with term | without |
|---|---|---|---|
| Kbb `OS.20240203.42705.45` (A) | 54.245 K | +1.917 A | +1.519 A |
| Kbb `OS.20240203.43383.79` (B) | 54.248 K | +1.913 A | +1.514 A |
| Hn3 `s160321_a002011` | 55.339 K | +0.409 A | -0.000 A |

The term adds +0.40 A at both epochs, moving the lines away from the
vacuum list; it is now off by default (`osiris_grating_temp = False`), and
the 2024 epoch offset (+1.5 A = +0.54 px) is absorbed by the refinement.

### 3.5 Davies OH scaling (task 3)

`skysub.oh_scaled_sky` follows `scaledskysub_000.pro`: OH lines of the B
spectrum (continuum removed by a 41-px running median, pixels within 4 px
of a listed line) are grouped by vibrational transition; one factor `r_g -
1` per group is fit (weighted least squares, MAD-based 3-sigma clipping,
pixels where the residual exceeds half the sky line excluded) and only the
line spectrum is rescaled (DRP `scale_fitted_lines_only`); factors that are
not 3-sigma significant, or that do not lower the residual of their lines,
are zeroed; the K tail group takes the factor of 9-7.  `FiberFindObjects`
calls it on the boxcar spectra of A-B and of B = A - (A-B) per lenslet
(pooling `global` by default) and stamps the model back onto the lenslet
image, where it becomes the sky model of the A-B frame (spec2d `SKYMODEL`).
Fitted factors on NGC 3245 are -3 % to +3 % for 600-s pairs and -8 to -10 %
for the 44-s pairs; on the Hn3 pair +26 / +27 %.

## 4. Datacube (task 5) and telluric (task 6)

`pypeit_fiber_cube` reads `PYP_SPEC`, calls `fiber_cube_setup` (OSIRIS:
`set_lenslet_state` from `SFILTER`, MJD, `SSCALE`; M3-10 respected), then
`build_cube_from_spec1d`.  NGC 3245 cube (`cube_OS.20240203.42705.45-...`):
(1474, 66, 21) = 19 lenslet columns x 64 rows + padding, 0.05" pixels,
2.82 A steps over 19651-23808 A, `CRPIX = (12, 32)` = layout origin (DRP
reference spaxel), `CRVAL` = header RA/DEC, rotation from the PA hooks;
continuum peak inside the field (vet test).  No DRP cube of this night is
available (`Kbb_050/ref/` is empty), so no direct comparison was made.
`pypeit_binospec_ifu_cube` is now a thin alias (same options, Binospec only).

Telluric: `pypeit_fiber_cube --point_source 0.4 --no_cube` sums the
lenslets within 0.4" of the star (184 lenslets) into a single-object
`_ps` spec1d; `pypeit_sensfunc` (IR, A0, V = 8.78) and `pypeit_tellfit
--objmodel star` run on it with `TellPCA_3000_26000_R10000.fits` (6.3 MB,
already in the PypeIt data tree; no download) and a polynomial order 7 for
the Kbb filter ripples.  The 8.85-s frames `56998.98` and `57040.80`
originally chosen for the dev suite contain **no star** (brightest lenslet
9-13 counts): the dev suite now uses `OS.20240203.59407.32` (nodded pair
with 59447).  Corrected/model: 1.5 % rms in 50-A bins, 1.0 % median
absolute deviation, 4.5 % max (telluric transmission > 0.5).  The 3-s
frames 56832 and 57720 give 2.7 % and 2.5 %.

## 5. Verification (commands and outputs)

### 5.1 Dev-suite reduce and afterburn

```
PYPEIT_DEV=$PWD ./pypeit_test reduce afterburn -s keck_osiris/Kbb_050 -o <scratch>/REDUX_DS
  keck_osiris/Kbb_050 pypeit                       PASSED  0:39:31
  keck_osiris/Kbb_050 pypeit_fiber_cube (NGC3245)  PASSED  0:00:13
  keck_osiris/Kbb_050 pypeit_fiber_cube (HD101060) PASSED  0:00:04
  keck_osiris/Kbb_050 pypeit_sensfunc              PASSED  0:00:04
  keck_osiris/Kbb_050 pypeit_tellfit               PASSED  0:00:24
--- PYPEIT DEVELOPMENT SUITE PASSED 5/5 TESTS ---
```
After switching the standard to 59407.32 the afterburn phase was re-run
(`./pypeit_test afterburn ...`): 4/4 PASSED.  Hn3:
`./pypeit_test reduce -s keck_osiris/Hn3_100`: PASSED in 3 min 32 s
(3040 spectra; OH scaling +26 %/+27 %, line-pixel rms 225 -> 62 counts).

### 5.2 Vet tests

```
pytest -s vet_tests/test_keck_osiris_kbb.py --redux_out <scratch>/REDUX_DS
  test_spec1d_lenslets  1207 lenslet spectra, AIRORVAC = vac;
                        Davies model on 0.061 of the pixels, rms 14.2 counts   PASSED
  test_cube_structure   cube (1474, 66, 21); CRPIX = 12.0, 32.0; peak (33, 10) PASSED
  test_telluric_flat    50-A binned rms 0.015, max |dev| 0.045, 79 bins        PASSED
pytest -s vet_tests/test_keck_osiris_flow.py  (tilt tests, both setups)
  Hn3_100 Keck: median |tilt| 0.073 px; refined on s160321_a002010: 0.023 px
  Kbb_050 Keck: median |tilt| 0.106 px; refined on OS.20240203.42705.45: 0.020 px
  (on the ISSKY=1 frame OS.20240203.43383.79; scatter 0.030/0.039 px)
  9 passed, 2 deselected
pytest vet_tests/test_keck_osiris_flow.py -k "spec1d or quasar" --redux_out <scratch>/REDUX_HN3
  3040 spectra; 3C 273 centroid (30.140, 22.898) vs DRP (30.146, 22.901)  2 passed
pytest vet_tests/test_keck_osiris_rectify.py                               5 passed
```

### 5.3 OH residuals (`make_m4_figures.py ... ohscale`, `m4_oh_scaling.png`)

Metric: for every isolated `OH_OSIRIS` line, the residual (continuum
removed) integrated over +-3 px and summed over lenslets, divided by the
sky-line flux of the B frame; median of |.| over lines.

| Data | A-B | A-B-Davies |
|---|---|---|
| Sky-sky pair 43383 - 45436 (r-1 = +7.6 % / +6.0 %) | 0.082 | 0.030 |
| NGC 3245, 13 x 600 s (median) | 0.053 | 0.037 |
| NGC 3245 42266 (44 s) | 0.040 | 0.068 |
| NGC 3245 42460 (44 s) | 0.074 | 0.081 |

### 5.4 NGC 3245 maps (`m4_ngc3245.png`)

Stack of the 13 600-s spec1d files on a common grid: smooth K continuum
peaked at lenslet (row 33, col 23); CO 2-0 bandhead depth 0.08 nearly
uniform (old bulge); H2 1-0 S(1) (redshifted to 21311 A) detected on the
nucleus (EW +0.39 A, 3x3 flux 5745 counts vs map noise 970), Br gamma not
detected (EW +0.05 A).  Both line maps show a horizontal step at flipped
lenslet row ~37 that is not physical.  Hypothesis (not verified): within a
block the lenslet stagger (1.92 px/row) moves the line and its side bands
across a raw-detector boundary (e.g. column 1024 or a readout-channel
edge) at that row, and an additive offset there biases the continuum
subtraction.  Q&A 3.

### 5.5 Unit tests

```
pytest -q pypeit/tests --ignore=pypeit/tests/test_pkgdata.py
900 passed, 9959 warnings in 113.56s
pytest pypeit/tests/test_keck_osiris.py test_lenslet_rectify.py test_fiber_sky.py
49 passed
```
`test_pkgdata.py` is excluded because it downloads from a GitHub branch that
does not exist and deletes a tracked file.

## 6. Deviations from the design document

1. The Davies scaling runs on the 1D lenslet spectra inside
   `FiberFindObjects` (the lenslet image has one column per lenslet) and is
   stamped back as a 2D sky model; pooling defaults to `global` (as the
   DRP) rather than per lenslet or block.  It is skipped for standards.
2. In A-B fiber frames the spec1d `*_COUNTS_SKY` columns hold the sky of the
   un-subtracted frame (PypeIt convention; zero here); the Davies model is
   in the spec2d `SKYMODEL` only.  The vet test checks the spec2d.
3. `DTMP7` is off by default (design 6.1 kept it until this re-test).
4. The telluric/sensfunc path needs a point-source extraction from the
   lenslet spectra (`--point_source`), not foreseen in the design.
5. The emission-line and continuum maps were built from the stack of the
   600-s spec1d files (the cube script builds one cube per exposure; no
   cube co-add yet).
6. Only one A/B pair carries `arc,tilt,trace` in the dev-suite file
   (combining 23 rectified frames for fixed calibrations cost ~15 min).

## 7. Open items / next (M5)

- Davies scaling worsens the 44-s frames (Q&A 2).
- Step in the line maps (Q&A 3).
- No OH refinement on the 1.5-s HD89239 frames (Q&A 4).
- Residual cosmic rays missed by L.A.Cosmic on single frames (visible as
  isolated hot lenslets in the line maps).
- `get_fixed_tilts` refits the 19-51 block polynomials for every frame
  (168 calls in the Kbb run); caching per configuration would save time.
- Shared-code changes above to flag for the PR reviewers; user
  documentation (tkserver matrices, new parameters, `pypeit_fiber_cube`,
  telluric recipe) in M5.
