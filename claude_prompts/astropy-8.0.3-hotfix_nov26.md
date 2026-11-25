# Plan: handle FITS BLANK pixels after astropy 8.0.3 (fixes the 23 Nov dev-suite failures)

## Context

`westfall-release-23nov26` failed 4 tests (3 `subaru_focas` reductions and `vet_tests/test_wavelengths.py::test_keck_hires`). They passed on 8 Nov with identical PypeIt and dev-suite code. Your local runs in `$ENVS_HOME/pypeit-main` (pre-update) and `$ENVS_HOME/pypeit-main-upd` (post-update) reproduce all of them exactly. Both are editable installs of the same checkout.

**Root cause (confirmed): astropy 8.0.3 (released 17 Nov) changed how pseudo-unsigned FITS images with a `BLANK` card are read.**
- Raw FOCAS and HIRES files have `BITPIX=16, BZERO=32768, BLANK=-32768`.
- astropy ≤ 8.0.2 (via `ImageHDU._convert_pseudo_integer`, now removed) returned `uint16` and ignored BLANK, so BLANK pixels read as **0 counts**.
- astropy 8.0.3 returns `float32` with BLANK pixels set to **NaN**.
- PypeIt has no guard for non-finite raw pixels.

**FOCAS.** The bias `RAW_DATA/subaru_focas/300B_None/FCSA00249756.fits` has one BLANK pixel, in amplifier 1's overscan.
- `procimg.subtract_overscan` computes `np.median(overscan)`, which returns NaN. All 256 columns of amplifier 1 in the bias then become NaN.
- Those NaNs propagate into the trace and science images. The trace image's NaN columns hide the leftmost slit edge, which leads to three different downstream errors from one cause.

**HIRES** (`J1723+2243_W241_RED_C5_ECH_-0.15_XD_0.90_2x2`). The raw arc has 2,048 / 12,669 / 26,553 BLANK pixels on its three chips. On each chip, one column is entirely BLANK, and the rest are saturated ThAr line cores.
- **Pre-update:** these pixels read as 0, so they appear as about −1957 after overscan subtraction, and are left unmasked.
- **Post-update:** they become BPM with value 0 in the mosaic (35,126 pixels).
- **Effect:** edges and slits are identical. Orders 24–36 change, and order 32's RMS goes from 0.083 to 0.259, above the threshold of 0.21.

**Decisions (yours):**
- Honor BLANK, i.e. treat those pixels as bad, rather than restoring the old read behavior.
- Flag these pixels with the **existing `IS_NAN` bit** ("Pixel value is undefined") together with `BPM`. Adding a new `ImageBitMask` bit was rejected: as-is it would make older `spec2d` and calibration files unreadable, a breaking change.
- No datamodel change, so this is **not breaking**.
- Leave `gen_kube_devsuite` unchanged.
- Work on the existing `ds_failures` branch in both repos.

## Update during implementation (24 Nov 2026): BLANK policy changed

The base `Spectrograph.get_rawimage` (`pypeit/spectrographs/spectrograph.py:1641`) already opened raw files with `ignore_blank=True` (added in 2023, commit 7b9304e53, as part of a Keck/NIRSPEC pre-upgrade fix). The 19 spectrographs with custom readers did not, so only they were affected by astropy 8.0.3. Given that, you chose to **ignore BLANK everywhere**:

- `pypeit/io.py:fits_open` now defaults to `ignore_blank=True` (callers can override). Its `ignore_missing_end` retry no longer raises a duplicate-keyword `TypeError`.
- Steps 1–4 (`raw_undef` tracking, BPM + IS_NAN flags, a NaN-safe overscan, IS_NAN through `build_mask` and `CombineImage`) are kept as a safety net for genuinely non-finite raw data. BLANK pixels are no longer flagged.
- Consequence: results under astropy 8.0.3 match pre-8.0.3 exactly. Under astropy 8.0.3, the FOCAS 300B_None bias is bit-identical to the pre-update product, and the HIRES J1723 arc has 0 undefined pixels. **Step 5 (the HIRES threshold review) should no longer be needed** if your `pypeit_test` runs confirm it.
- New tests: `pypeit/tests/test_rawimage_undefined.py`, plus additions to `pypeit/tests/test_procimg.py`.

## Implementation (existing `ds_failures` branch in both repos)

### 0. Save this plan (first action after approval)
Copy this plan verbatim to `/Users/westfall/Work/packages/pypeit-main/PypeIt-development-suite/claude_prompts/astropy-8.0.3-hotfix_nov26.md` (dev-suite repo, `ds_failures` branch). Don't commit it unless you ask.

**Testing split:** you run every `pypeit_test` execution (dev-suite reduce, vet and Nautilus runs). These runs don't gate the code changes: all of steps 1–4, 6 and 7 are implemented without waiting for them. I run only `pytest pypeit/tests` and the documentation build, and I give you the exact `pypeit_test` commands to run.

All PypeIt changes go on `ds_failures` in `pypeit-main/pypeit`. Any dev-suite change (the `vet_tests` threshold in step 5) goes on `ds_failures` in `pypeit-main/PypeIt-development-suite`.

### 1. Track undefined raw pixels in `RawImage` (`pypeit/images/rawimage.py`)
- In `__init__`, right after `spectrograph.get_rawimage` (around lines 129–160):
  - Build `self.raw_undef = np.logical_not(np.isfinite(self.rawimage))`, with the same `expand_dims` handling.
  - Set those pixels to a finite placeholder (0, the pre-8.0.3 value) in `rawimage`/`image`.
  - Log a warning with the count when any are found.
- Carry `raw_undef` through the same geometry steps as `datasec_img`/`_bpm`:
  - `trim()` (around line 1350, `procimg.trim_frame`)
  - `orient()` (around line 906, `spectrograph.orient_image`)
  - `build_mosaic()` (around lines 1386–1440): use the same `build_image_mosaic` transform and rounding as `_bpm`, and treat pixels outside the detector footprint as not undefined
- In `process()`:
  - After the processed-shape BPM is built (around lines 612–637), OR `raw_undef` into `self._bpm`. Every consumer that selects `'BPM'` by name then excludes these pixels without further changes: edgetrace (4 sites), flatfield (4), `wavecalib`/`wavetilts` (`['BPM','CR']`), and the lacosmic bpm.
  - After the `PypeItImage` is constructed (around lines 700–725), run `pypeitImage.update_mask('IS_NAN', indx=raw_undef)`.

### 2. Keep `IS_NAN` through mask rebuilds (`pypeit/images/pypeitimage.py`)
- `build_mask(from_scratch=True)` (around lines 590–600) currently saves and restores only `BPM` and `CR`. Restore the existing `IS_NAN` flags in the same way, then continue to OR in `~isfinite(image)`. Update the docstring ("all bits are turned off except for BPM, CR and IS_NAN").
- Risk to verify: a pixel that was NaN earlier and later became finite would now keep `IS_NAN`. Check that nothing in the pipeline relies on `IS_NAN` being cleared. Search for consumers of `IS_NAN`, `bad_flags` and `flagged(` (mostly "any bit" selections, which are unaffected because these pixels are already BPM).

### 3. Propagate through combining (`pypeit/images/combineimage.py`, around lines 160–310)
- While stacking, also collect `rawImage.select_flag(flag='IS_NAN')` per frame.
- After `comb.build_mask(...)`, set `comb.update_mask('IS_NAN', indx=np.all(isnan_stack, axis=0))`, i.e. the pixel is undefined in every input frame. Those pixels are already `~gpm`, so they are BPM/STCKMASK.
- Pixels undefined in only some frames are already excluded per frame through `gpm_stack`, so the other frames fill them.

### 4. Exclude bad pixels from the overscan statistics (`pypeit/core/procimg.py:subtract_overscan`, around lines 710–850)
- Add an optional `bpm` argument (default None, which keeps current behavior). Compute the per-amplifier overscan statistics from good pixels only, for the `median`, `savgol`, `chebyshev`/`polynomial` and `odd_even` paths (e.g. `np.nanmedian` on a copy with bad pixels set to NaN).
- Pass `raw_undef` from `RawImage.subtract_overscan` (around line 1084). Check the other callers (`grep -rn "subtract_overscan(" pypeit`); the default leaves them unchanged.

### 5. Evaluate HIRES under the fix (after you've run it; not gating)
- You rerun J1723 `-0.15` in `pypeit-main-upd`: `./pypeit_test reduce -i keck_hires -s J1723+2243_W241_RED_C5_ECH_-0.15_XD_0.90_2x2`. I then analyze the products. The BLANK pixels are now BPM and IS_NAN, much like post-update, so expect order 32's RMS to stay near 0.259.
- Compare order 32 pre- vs. post-fix (`pypeit_chk_wavecalib`, WaveCalib QA): which lines were identified or rejected, and whether masked saturated cores bias the line centroids.
- Bring the outcome back to you for a decision:
  - (a) update the J1723 RMS threshold in dev-suite `vet_tests/test_wavelengths.py`, with a justification comment
  - (b) improve how the arc handles masked or saturated line cores
- Also run the other `keck_hires` setups and an echelle/mosaic sample.

### 6. Tests, scope check, documentation
- **Unit tests** (`pypeit/tests`, plain functions, fixed seeds). They must pass under astropy 8.0.0 and 8.0.3: inject the NaN directly into the array, or open the file with `ignore_blank=False`.
  - `subtract_overscan` with a NaN or flagged pixel in the overscan gives a finite result equal to the median of the good pixels (`median` and `savgol` methods).
  - `RawImage` with a synthetic `uint16` FITS file (`BZERO=32768`, a `BLANK` pixel in the overscan and in the data region), using a minimal spectrograph or a monkeypatched `get_rawimage`. Check that the result has no NaNs, the data-region pixel has `BPM` and `IS_NAN` set, and the overscan pixel has no effect.
  - `build_mask(from_scratch=True)` preserves `IS_NAN`.
  - `CombineImage`: a pixel undefined in all frames is `IS_NAN`; a pixel undefined in one frame is not.
- **Scope check:** scan the headers of `$PYPEIT_DEV/RAW_DATA/**` for `BLANK` with integer `BITPIX` (a read-only scan, which I run) and give you the list of affected setups, so you can run their `reduce` tests in `pypeit-main-upd`.
### 7. Hotfix release notes (procedure in `doc/dev/development.rst`, "Hotfixes", lines 377–401)
- This is a patch release, 3.0.1, branched from `release`. `ds_failures` is already based on `release` @ c1798aa8a.
- **Create `doc/releases/3.0.1.rst`**, following the style of `doc/releases/2.0.1.rst`: a `Version 3.0.1` title with a `=` underline, then bullets. Describe:
  - astropy ≥ 8.0.3 converts FITS `BLANK` pixels in pseudo-unsigned raw images to NaN.
  - PypeIt now flags undefined raw pixels with `BPM` and `IS_NAN`, and excludes them from the overscan statistics.
  - This fixes NaN-poisoned amplifiers (e.g. Subaru/FOCAS) and how saturated pixels stored as BLANK are handled (e.g. Keck/HIRES).
  - Mention any HIRES vet-threshold or arc-handling outcome from step 5.
  - No datamodel or breaking changes.
- **Edit `doc/whatsnew.rst`:** add `.. include:: releases/3.0.1.rst` followed by `----` above the `releases/3.0.0.rst` include.
- **No development doc to rename on `release`.** `3.1.0dev.rst` exists only on `develop`. Per the procedure (lines 366–371), the hotfix entry must also be added to `develop`'s `doc/releases/3.1.0dev.rst` when the fix is merged forward to `develop`. Note this in the PR description.
- Don't use the `update-changelog` skill for this, since it targets the in-development doc.

## Critical files
- `pypeit/images/rawimage.py`: `__init__`, `trim`, `orient`, `build_mosaic`, `subtract_overscan`, `process`
- `pypeit/images/pypeitimage.py`: `build_mask`
- `pypeit/images/combineimage.py`: mask stacking and combined-mask construction
- `pypeit/core/procimg.py`: `subtract_overscan`
- `pypeit/tests/`: new or extended tests
- `doc/releases/3.0.1.rst` (new) and `doc/whatsnew.rst`
- Dev-suite `vet_tests/test_wavelengths.py`: only if step 5 settles on option (a)
- Not changed: `pypeit/images/imagebitmask.py`, where no new bit is added

## Verification
**Run by me:**
- `pytest pypeit/tests` in both environments.
- `cd doc; make htmlonly`.

**Run by you** (`pypeit_test`; I analyze the outputs afterward):
- **FOCAS:** in `pypeit-main-upd` (astropy 8.0.3), `./pypeit_test reduce -i subaru_focas` passes all three setups:
  - 300B_None finds 3 left edges again
  - `Bias_A_0_DET01.fits` has no NaNs
  - results match the pre-update run, apart from negligible overscan-median differences
- **FOCAS on old astropy:** the same command passes in `pypeit-main` (astropy 8.0.0).
- **HIRES:** J1723 is resolved through step 5. The BLANK pixels show `BPM` and `IS_NAN` in the Arc full mask, and `vet` `test_keck_hires` passes in `pypeit-main-upd`.
- **Unit tests** (run by me): `pytest pypeit/tests` passes in both environments.
- **Docs** (run by me): `cd doc; make htmlonly` builds without new warnings, and 3.0.1 appears in "What's New".
- **Final check** (run by you): a full Nautilus `all` run with `-p ds_failures -d ds_failures` passes 335/335.
