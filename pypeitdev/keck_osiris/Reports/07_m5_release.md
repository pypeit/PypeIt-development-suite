# 07 - Milestone M5: hardening, documentation and release

Date: 2026-10-13.  Author: Claude (Opus 5.5) with JXP.
Task: Prompt 1 of `claude_prompts/keck_osiris/M5_release.md` (Stage 5 of
`keck_osiris_implementation.md`).  Environment: laptop, conda env
`pypeit14`, PypeIt branch `osiris_port`, dev-suite branch `dev_osiris`.
No git state was changed (read-only `git status/diff/log/archive` only).
All reductions ran under the session scratchpad
(`.../scratchpad/REDUX_M5A`, `REDUX_M5B`, `REDUX_KN3`); `REDUX_DS` is the M4
reduction used as the "before" reference.

## 1. Summary

| Acceptance item | Target | Measured | Status |
|---|---|---|---|
| Dev-suite reduce + afterburn, all OSIRIS setups | pass | `./pypeit_test reduce afterburn vet -s keck_osiris/Hn3_100 keck_osiris/Kbb_050 keck_osiris/Kn3_035`: **7/7 PASSED** (Hn3 3 min 38 s, Kn3 5 min 33 s, Kbb 38 min 02 s; NGC 3245 and HD101060 fiber cubes, sensfunc, tellfit) | pass |
| Dev-suite vet, OSIRIS | pass | **21/21 passed** (`test_keck_osiris_rectify` 5, `_flow` 11, `_kbb` 3, `_kn3` 2); the `vet` phase also runs the other instruments' vet tests, which fail for lack of their reductions (96 failed / 2 errors, none OSIRIS) | pass |
| `cd doc; make htmlonly` | no new warnings | **0 warnings** (baseline tree at HEAD: 2, both from M3/M4 docstrings, now fixed) | pass |
| PR description | ready | section 8 and `pypeitdev/keck_osiris/pr_description.md` | done |
| Davies scaling no longer degrades any NGC 3245 frame | no frame worse | 44-s frames **0.0402 -> 0.0402** and **0.0735 -> 0.0735** (M4: 0.068, 0.081); **no frame worse** on either metric definition; 600-s median 0.0528 -> **0.0369** (M4 0.0362) | pass; 600-s "keep improvement": 10/13 improve, 3 equal (Q&A 1) |
| Line-map step fixed or documented | evidence | not instrumental: fixed +-8 A window vs H2 rotation of +-230 km/s (step ratio 3.6 -> 1.5 with +-25 A), plus an intrinsic H2 edge; no step in a sky-sky pair (21 vs 32 counts x A, scatter 86) nor in the sky continuum; per-block step row independent of the raw column of the line (1011-1273) | documented (Q&A 2) |
| 1.5-s standards carry the science OH refinement | ~0.1 px | HD89239 vs refined 600-s NGC 3245: **+0.538 px -> +0.001 px** median, per-lenslet std 0.014 px, p95 0.027 px (each of the 4 frames; vs each 600-s frame -0.001..+0.003) | pass |
| Existing + new unit tests | pass | `pytest pypeit/tests --ignore=pypeit/tests/test_pkgdata.py`: **908 passed** (M4: 900) | pass |
| Another KOA mode (task 1) | run | **Kn3/0.035"** (2024-08-13, nova V1405 Cas, 2 A/B pairs): reduces unchanged, 3040 lenslets, all >= 394 usable pixels, Br gamma EW 27-29 A; new dev-suite setup + vet test | pass (tilt rms Q&A 3) |
| Tilt cache / residual CRs (task 10) | cache; quantify | Kbb: `get_fixed_tilts` fits **1x** + 167 cache hits (M4: 168 fits, 0.6 s each); 44064: 8-sigma residual spikes **260 -> 3** events, 5-sigma 558 -> 302 | pass |

## 2. What was built

### 2.1 PypeIt repo (`/Users/xavier/Projects/PypeIt/PypeIt`)

OSIRIS-specific:

| File | Change |
|---|---|
| `pypeit/spectrographs/keck_osiris.py` | `supported = True`; module docstring; `_tilt_cache` (Keck-derived tilt coefficients per filter / coefficient file / expansion / nslits, shared by instances); `_set_state_from_header` names missing `SFILTER`/`MJD-OBS`; `rect_matrix` error for a missing directory says how to set `osiris_rectmat_dir`; default `fiber_cr_nsig = 8`; table comment links the doc page |
| `pypeit/core/lenslet_rectify.py` | `select_matrix`: missing-directory error with the parameter and URL; the no-match error lists matrices *and* other files (first 20) with the naming rule |
| `doc/spectrographs/keck_osiris.rst` | new user page (overview, supported data, data requirements, matrices and epoch rule, typing, `bkg_id`, darks/ITIME, processing, parameter table, outputs, `pypeit_fiber_cube`, telluric recipe, limitations, references) |
| `doc/tutorials/keck_osiris_howto.rst` | new tutorial (Kbb/050 night end to end) |
| `doc/spectrographs/spectrographs.rst`, `doc/tutorials/tutorials.rst`, `doc/include/spectrographs_table.rst`, `doc/pypeit_par.rst` | links; regenerated table and parameters (`build_spectbl_rst.py`, `build_par_rst.py`; datamodel and script-help regeneration produced no change other than a version string, reverted by hand) |
| `doc/releases/2.1.0dev.rst` | Documentation, Functionality, Instrument-specific (Keck/OSIRIS), Script and Testing entries (`CHANGES.rst` is deprecated; not edited) |

Generic shared code (gated on parameters; flagged to the reviewers):

| File | Change |
|---|---|
| `pypeit/core/skysub.py` | `oh_scaled_sky(guard='integrated', guard_lines, guard_halfwidth)`: subset-search guard on the integrated line residual (`'rms'` keeps the M4 rule); `stats['int_before'/'int_after']`; new `oh_line_residuals`; new `fiber_cosmic_rays` |
| `pypeit/find_objects.py` | `FiberFindObjects`: isolated guard lines passed to the OH scaling; `osiris_oh_min_exptime` floor; `_arc_fiber_sky_shifts` (arc-frame fallback, cached per file); `_fiber_residual_crs` (second CR pass, A-B science frames); `global_skysub` docstring fixed (doc warning) |
| `pypeit/par/pypeitpar.py` | new `[reduce][skysub] osiris_oh_min_exptime` and `fiber_cr_nsig` (`add-parameter` procedure); descriptions of `osiris_oh_scaling`, `fiber_sky_refine`, `osiris_rectmat_dir` (epoch rule) updated |
| `pypeit/core/wavecal/wvutils.py` | `fiber_sky_line_shifts` silences the per-line fit warnings (Kbb log: 2688 -> 0 `OptimizeWarning`s) |
| tests | `test_keck_osiris.py` (+5: truncated-header hooks, missing cards, era boundary, matrix errors, tilt cache; `supported` assertion), `test_fiber_sky.py` (+3: `oh_line_residuals`, guard vs M4 rule, `fiber_cosmic_rays`) |

### 2.2 Dev suite

| File | Change |
|---|---|
| `RAW_DATA/keck_osiris/Kn3_035/` | new: 4 KOA frames (2 A/B pairs, 300 s), `calib/s240523_c002___infl_Kn3_035.fits.gz` (tkserver, epoch 230701-250310, gzip verified), `koa_file_list.csv` |
| `pypeit_files/keck_osiris_kn3_035.pypeit` | new (A/B via `bkg_id`, `arc,tilt,trace` on the first pair) |
| `test_scripts/setups.py` | `Kn3_035` registered (Hn3_100 and Kbb_050 were already registered with their afterburn tests in M4) |
| `unit_tests/test_load_images.py` | Kn3 frame added (exptime 299.481 s) |
| `vet_tests/test_keck_osiris_kn3.py` | new |
| `pypeitdev/keck_osiris/Reports/figures/make_m5_figures.py` | figures `m5_oh_guard`, `m5_std_shift`, `m5_linemap`, `m5_crs`, `m5_kn3` |
| `pypeitdev/keck_osiris/pr_description.md` | PR text (both PRs) |
| design doc, implementation plan | M5 status (decision row "M5", sections 4.8, 6.2, 8; Stage 5 status) |

## 3. Tasks

### 3.1 Task 1: another KOA mode (Kn3/0.035")

KOA TAP (`koa_osiris`, public rows only) lists 207 Kn3 frames >= 300 s at
0.035/0.050" since 2016; most nights lack sky frames.  Chosen: 2024-08-13,
V1405 Cas (nova 2021), Kn3/0.035", three 300-s object/sky pairs; two pairs
downloaded (`getKOA`, 4 x 37.9 MB) and the epoch matrix from
`http://tkserver.keck.hawaii.edu/osiris/data/230701-250310/` (25.8 MB gz,
`gzip -t` ok, used gzipped).  The B frames are 1" away, so the nova appears
negative in the upper part of each A-B image (expected).

`run_pypeit` (5 min 33 s): 3040 lenslet spectra per frame (51 blocks of the
narrowband layout); every lenslet keeps >= 394 of ~462 window pixels
(`m5_kn3.png`); OH refinement 2918/3040 lenslets (median 5 lines,
-0.24 px); the one Davies group in band (9-7) is rejected as insignificant;
192 residual CR pixels.  Nova continuum peaks at lenslet (34, 28) and
(46, 26) (dithered); Br gamma EW 26.6 / 28.6 A; `pypeit_fiber_cube` builds a
(468, 68, 53) cube and the point-source spectrum.  **No lenslet-mapping
edge case appeared**: the narrowband layout and band windows of M3/M4 hold.
The Keck-derived tilts reach 0.065 px rms (robust 0.025 px) on the 2024 Kn3
coefficients, dominated by 40 interior lenslets with peculiar coefficients
(0.035 px without them; Hn3 2016: 0.048, Kbb 2024: 0.044); those lenslets
keep their exact solution in the spec1d (Q&A 3).

### 3.2 Task 2: robustness

Unit tests: the 2016 truncated header passes through every header hook
(`frame_expansion` = 1 without `DTMP7`, fixed wavelengths and tilts, cube
setup, all pypeit-file keys); headers without `SFILTER`/`MJD-OBS` raise a
`PypeItError` naming the card; the H2RG era starts at MJD 57388 inclusive;
missing matrix directory -> error naming `osiris_rectmat_dir` and the URL;
wrong matrices -> error listing the matrices and the misnamed file, epoch
230701-250310 and tkserver.  `pypeit_setup` on the 2015 Kn5 frames warns
that both predate the H2RG (existing `vet_assigned_ftypes`).

### 3.3 Task 7: Davies guard

The M4 guard (pixel rms of the line pixels) accepted factors that increase
the integrated residual of weak lines.  M5: for every pool, all subsets of
the significant group factors are tried and the one with the lowest median
over the isolated lines of |integrated residual| / sky-line flux is kept;
none unless the median decreases (exactly the Report 06 metric, measured on
the pipeline's per-lenslet refined wavelengths).  A first M5 version
(per-group median, then per-group sum) still dropped or flipped small
factors; the subset rule is the one adopted.  `osiris_oh_min_exptime`
(default None) is available but not needed.

Metric (`make_m5_figures.py ... ohscale`; M4 definition on the spec2d
`WAVEIMG`, and the same with the spec1d per-lenslet refined wavelengths):

| Frame | t (s) | A-B | Davies M4 | Davies M5 | refined: A-B | M4 | M5 |
|---|---|---|---|---|---|---|---|
| 42266.81 | 44 | 0.0402 | 0.0675 | **0.0402** | 0.0309 | 0.0597 | 0.0309 |
| 42460.10 | 44 | 0.0735 | 0.0808 | **0.0735** | 0.0693 | 0.0689 | 0.0693 |
| 42705.45 | 599 | 0.0528 | 0.0347 | 0.0437 | 0.0491 | 0.0380 | 0.0424 |
| 44064.64 | 599 | 0.0387 | 0.0329 | 0.0387 (=) | 0.0322 | 0.0354 | 0.0322 |
| 44757.73 | 599 | 0.0475 | 0.0391 | 0.0391 | 0.0486 | 0.0440 | 0.0441 |
| 46116.63 | 599 | 0.0501 | 0.0358 | 0.0358 | 0.0507 | 0.0377 | 0.0377 |
| 46816.35 | 599 | 0.0634 | 0.0480 | 0.0511 | 0.0579 | 0.0369 | 0.0443 |
| 48172.98 | 599 | 0.0838 | 0.0464 | 0.0464 | 0.0882 | 0.0515 | 0.0514 |
| 48868.72 | 599 | 0.0821 | 0.0366 | 0.0366 | 0.0839 | 0.0455 | 0.0455 |
| 50226.88 | 599 | 0.0369 | 0.0373 | 0.0369 (=) | 0.0395 | 0.0396 | 0.0392 |
| 50920.99 | 599 | 0.1332 | 0.0542 | 0.0542 | 0.1324 | 0.0537 | 0.0538 |
| 52278.78 | 599 | 0.0313 | 0.0281 | 0.0304 | 0.0306 | 0.0274 | 0.0296 |
| 52965.80 | 599 | 0.0258 | 0.0236 | 0.0258 (=) | 0.0300 | 0.0302 | 0.0300 |
| 54326.00 | 599 | 0.0770 | 0.0362 | 0.0362 | 0.0854 | 0.0424 | 0.0423 |
| 55001.79 | 599 | 0.0638 | 0.0343 | 0.0342 | 0.0713 | 0.0408 | 0.0408 |

No frame is worse with M5 on either definition.  The 600-s median goes
0.0528 -> 0.0369 (M4 0.0362).  On the M4 definition 10/13 600-s frames
improve and 3 are unchanged (44064, 50226, 52965); on the refined
wavelengths 11/13 improve.  For 44064 and 52965 the M4 "improvement" was not
real on the refined metric (M4 made them slightly *worse*: 0.0322 ->
0.0354, 0.0300 -> 0.0302); the WAVEIMG metric uses the block wavelengths
without the +0.54 px per-lenslet refinement.  M5 gives smaller gains than
M4 on 42705, 46816 and 52278, where the subset search drops one group
factor.  The "keep their improvement" clause is therefore met for 10 (or
11) of 13 frames, not 13; see Q&A 1.

### 3.4 Task 8: line-map step

Evidence (`make_m5_figures.py ... linemap`, `m5_linemap.png`; stack of the
13 M4 spec1d files):

1. Column-1024 / readout-channel hypothesis rejected: in the central blocks
   (cols 18-27) the step sits at flipped row 35-36 in every block while the
   raw detector column of H2 at that row runs from 1011 to 1273 (+29 px per
   block); a fixed raw-column boundary would give a slanted pattern.
   Lenslet row 28-29 (crow) sits at raw rows ~1100-1140.
2. Instrument-only data show nothing: a sky-sky A-B pair (43383 - 45436)
   measured with the same on/side windows gives 32 vs 21 counts x A on
   either side (row-to-row scatter 86); the galaxy pair 42705 - 43383 jumps
   120 -> 587.  The sky continuum of a single sky frame changes by < 3 %
   across those rows; raw-frame band medians show no feature.  No
   readout-channel structure is visible in raw darks or sky frames at the
   1e-3 ADU/s level.
3. Kinematics: H2 1-0 S(1) velocity per lenslet row goes from +230 km/s
   (flipped rows >= 37) through ~0 near the nucleus to -200 km/s (rows
   <= 22) -- the rotating nuclear gas of NGC 3245.  The +-8 A window
   (+-113 km/s) loses the line exactly at the step rows.  With +-25 A the
   contrast across the step drops from 3.6 to 1.5.
4. The remaining edge (H2 brighter on the nucleus side, sharp at flipped
   row ~36) has no counterpart in the continuum or in any sky-only test:
   an intrinsic feature of the H2 emission (e.g. obscuration by the
   nuclear dust disk), not something PypeIt can fix.

No DRP comparison: IDL is not installed (`which idl`: not found).  Noted
also: the M4 maps stack 13 frames dithered by up to 8 lenslet columns
without registration.  Documented in the user page (Limitations).

### 3.5 Task 9: OH refinement of very short frames

When neither the lenslets nor the stacked spectrum of a frame give a line,
`FiberFindObjects` measures the shifts on `Arc_<calib_key>.fits` (in the
dev-suite file the combined 600-s A/B pair 42705/43383), extracted along the
same traces, and applies them with a warning (`median +0.537 px; 1207
fibers shifted`).  Fallback order: own lenslets -> own stack -> arc.  The
3/8.85-s HD101060 frames still use their own stack (+0.51/+0.54 px).

| HD89239 frame (1.5 s) | M4 offset vs refined 600-s frames | M5 |
|---|---|---|
| 41472.23, 41635.43, 41774.15, 41802.20 | +0.538 px (std 0.045, p95 0.589) | **+0.001 px** (std 0.014, p95 0.027) |

(per-lenslet median over pixels of (lambda_std - lambda_sci)/dispersion,
heliocentric correction removed; range over the 13 science frames
-0.001..+0.003 px).

### 3.6 Task 10: tilt cache and residual cosmic rays

Tilt cache: one fit per configuration (Kbb 0.6 s, Hn3 ~2 s); Kbb run
1 fit + 167 hits (M4: 168 fits), Hn3 1 + 3, Kn3 1 + 7.

Second CR pass (`fiber_cr_nsig = 8`): a high-pass spike > 8 sigma that is
> 3x the residual of the two neighbouring lenslets on either side at the
same wavelength (scaled up by the continuum ratio, so bright-fiber
absorption lines survive; unit test) is flagged with the CR bit (added to
existing flags); standards are skipped.  On 44064 (same detector run on M4
and M5 spec2d): 8-sigma spikes 260 events in 223 lenslets -> 3; 5-sigma
558 -> 302 (the rest are mostly lower-amplitude, partly real structure).
About 70 of the 260 are not cosmic rays but the raw column at spectral row
80 (raw column ~1967) across many lenslets; masking them is harmless.
Per frame: 96-340 pixels (A-B, 600 s), 101-103 (44 s), 192 (Kn3, 300 s).
Frame-combination rejection was not used (frames are reduced singly).

## 4. Documentation (task 3)

Baseline: `git archive HEAD` exported to the scratchpad, `make htmlonly`
with `PYTHONPATH` on the copy: **2 warnings** (`FiberFindObjects.global_skysub`
field list; undefined label `pypeit_fiber_cube` in the Binospec alias
docstring).  After M5 (`cd doc; make htmlonly`, full `-aE` rebuild):
**0 warnings** (`sphinx_warnings.out` empty; the label now exists in
`keck_osiris.rst`, docstring fixed).

## 5. Google Drive staging (task 4; nothing uploaded)

Copy to the dev-suite drive under `RAW_DATA/keck_osiris/`:

| Directory | Files | Size |
|---|---|---|
| `Hn3_100/` | `s160321_a002010.fits`, `s160321_a002011.fits` | 75.7 MB |
| `Hn3_100/calib/` | `s160318_c005___infl_Hn3_100.fits` | 159.4 MB |
| `Hn3_100/ref/` | `s160321_a002010_Hn3_100_ref.fits` (DRP cube, rectify vet test) | 12.5 MB |
| `Kbb_050/` | 73 `OS.20240203.*.fits` + `koa_file_list.txt`, `koa_20240203_all_spec_frames.csv` | 2765 MB |
| `Kbb_050/calib/` | `s240524_c006___infl_Kbb_050.fits` (or the 32-MB `.fits.gz` from tkserver, Q&A 4) | 159.4 MB |
| `Kn3_035/` | `OS.20240813.{45703.90,46043.26,46751.14,47096.41}.fits`, `koa_file_list.csv` | 151.5 MB |
| `Kn3_035/calib/` | `s240523_c002___infl_Kn3_035.fits.gz` | 25.8 MB |
| `Kn5_035/` | `s150531_a025001.fits`, `s150531_a025002.fits` (era-error tests) | 75.6 MB |

Total ~3.4 GB (~3.3 GB with the gzipped Kbb matrix).  Do not copy the empty
`Kbb_050/ref/`, nor `Kn5_035/calib/s150901_c008___infl_Kn5_035.fits` (a
254-byte stub, not a matrix) or `Kn5_035/ref/` (not used by any test).
Pypeit/afterburn inputs are in git: `pypeit_files/keck_osiris_{hn3_100,
kbb_050,kn3_035}.pypeit`, `sensfunc_files/keck_osiris_kbb_050.sens`,
`tellfit_files/keck_osiris_kbb_050.tell`.  Registered: `setups.py`
(Hn3_100, Kbb_050, Kn3_035), `test_setups.py` (Kbb fiber cubes, sensfunc,
telluric; M4).

## 6. Verification commands

```
PYPEIT_DEV=$PWD ./pypeit_test reduce afterburn vet -s keck_osiris/Hn3_100 \
    keck_osiris/Kbb_050 keck_osiris/Kn3_035 -t 2 -o <scratch>/REDUX_M5B
  1 PASSED keck_osiris/Hn3_100 pypeit                      0:03:38
  2 PASSED keck_osiris/Kn3_035 pypeit                      0:05:33
  3 PASSED keck_osiris/Kbb_050 pypeit                      0:38:02
  4 PASSED keck_osiris/Kbb_050 pypeit_fiber_cube (NGC3245) 0:00:14
  5 PASSED keck_osiris/Kbb_050 pypeit_fiber_cube (HD101060) 0:00:04
  6 PASSED keck_osiris/Kbb_050 pypeit_sensfunc             0:00:27
  7 PASSED keck_osiris/Kbb_050 pypeit_tellfit              0:01:31
  --- PYPEIT DEVELOPMENT SUITE PASSED 7/7 TESTS ---
  vet phase (whole vet_tests dir): 96 failed, 23 passed, 2 errors -- all
  failures are other instruments without reductions in REDUX_M5B
pytest -s vet_tests/test_keck_osiris_{flow,kbb,kn3,rectify}.py --redux_out <scratch>/REDUX_M5B
  OH tilt: Hn3 0.073 -> 0.023 px; Kbb 0.106 -> 0.020 px (refined)
  3C 273 centroid (30.140, 22.898) vs DRP (30.146, 22.901)
  Kbb: 1207 spectra, vac; cube (1474, 66, 21); telluric 50-A rms 0.015
  Kn3: 3040 spectra x 2, min 394 good px; Br gamma EW 26.6 / 28.6 A
  rectify: drp 0.74 %, lstsq 0.77 % median vs DRP cube
  21 passed
pytest -q pypeit/tests --ignore=pypeit/tests/test_pkgdata.py   908 passed
cd doc; make htmlonly                                             0 warnings
```
`test_pkgdata.py` is excluded because it downloads from a GitHub branch
that does not exist and deletes a tracked file.  An earlier full run
(`REDUX_M5A`, before the final guard) also passed 6/6 (Hn3, Kbb + 4
afterburn).

## 7. Deviations from the design document

1. Davies guard: subset search on the integrated metric (always on);
   the exposure-time floor exists but defaults to None.
2. Short-frame OH fallback uses the `Arc` calibration image, not a
   "median of the night's science frames" (standards are reduced before
   science frames, so no science shifts exist yet).  Own stacked spectrum
   still has priority over the arc.
3. Line-map step: no code change; diagnosis astrophysical; no DRP run.
4. Kn3/035 Keck-derived tilts exceed the 0.05 px rms of decision M4-13
   (0.065 px, robust 0.025).
5. The dev-suite `vet` phase is not restricted by `-s`; OSIRIS vet results
   were also run directly.

## 8. Pull-request text

The full text (PypeIt PR and dev-suite PR, cross-referenced, crediting
Lyke et al. 2017 and Lockhart et al. 2019, shared-code changes of M3, M4
and M5 itemised for the reviewers, ending with the required attribution
line) is in `pypeitdev/keck_osiris/pr_description.md`.  Its PypeIt summary:

> Adds `keck_osiris`, the Keck I AO lenslet IFU (H2RG era), as a supported
> spectrograph.  Raw frames are rectified with the Keck matrices
> (`pypeit/core/lenslet_rectify.py`: DRP-solver port and an exact
> least-squares solver; matrix selection by validity epoch) into a lenslet
> image reduced by the `Fiber` pypeline with one spectrum per lenslet;
> Keck per-lenslet wavelengths refined on OH lines, A-B sky with Davies OH
> scaling, a second cosmic-ray pass, `pypeit_fiber_cube` (cubes and
> point-source spectra) and the IR telluric tools.  Shared-code changes
> for review: `runtime_par`, the `FiberExtract` column-window boxcar, the
> `rectify` step, fixed tilts/wavelength hooks (+ `raw_header`),
> `bkg_redux_sciimg`, `fiber_sky_refine`, `osiris_oh_scaling` (+ guard,
> `osiris_oh_min_exptime`), `rectify=False` for bias/dark, `fiber_cube_setup`
> and `pypeit_fiber_cube`, the arc-frame OH fallback, `fiber_cr_nsig`.
> Credits: Keck OSIRIS DRP (Lyke et al. 2017; Lockhart et al. 2019),
> Davies (2007).

## 9. Figures

`m5_oh_guard.png` (OH metric per frame), `m5_std_shift.png` (HD89239
offsets M4/M5), `m5_linemap.png` (H2 maps +-8/+-25 A, H2 velocity per row,
sky-sky test), `m5_crs.png` (A-B residual image of 44064, M4/M5 masks),
`m5_kn3.png` (Kn3 continuum images, nova spectrum).  Regenerate with
`PYPEIT_DEV=... python make_m5_figures.py <REDUX_M5B>/keck_osiris/Kbb_050
<REDUX_DS>/keck_osiris/Kbb_050 <REDUX_M5B>/keck_osiris/Kn3_035`.
