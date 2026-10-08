# Develop the SOAR/Triplespec spectrograph

## Goals

We wish to develop the SOAR/Triplespec spectrograph for PypeIt.  

## Claude

### Skills

Consider using the skills in PypeIt-development-suite/.claude/skills/

## Coding

Here are guidelines for coding: 

- Use Python where possible
- The main spectrograph class is in PypeIt/pypeit/spectrographs/soar_tspec.py
- Add inline comments to explain the effort
- Reuse existing code when possible
- Place import statements at the top of the file.
- Include a description of inputs/outputs in the doc string of all methods

## Testing

If you need to test the code: 

- use the files in PypeIt-development-suite/pypeitdev/soar_tspec
- The PypeIt file is PypeIt-development-suite/pypeitdev/soar_tspec/TripleSpec/soar_tspec_A.pypeit
- run on the "pypeit14" conda environment.

## Prep

1. Begin by running PypeIt on the PypeIt file  PypeIt-development-suite/pypeitdev/soar_tspec/TripleSpec/soar_tspec_A.pypeit.  Inspect the QA files generated and assess the current state of the pipeline.  Be sure to log below with a new Subsection; include any Skills used and the commands used to perform the task.  Make recommendations for the next stage of development.  

2. I have modified the PypeIt file to use the science frames for wavelength calibration.  Run PypeIt again and inspect the QA files.  Be sure to log below with a new Subsection.  Make recommendations for the next stage of development.  

3. I removed the existing Calibration files. Run PypeIt again and inspect the QA files.  Be sure to log below with a new Subsection.  Make recommendations for the next stage of development.  



## Development

### Wavelength Calibration

The current wavelength calibration fails only Order 3 of the 5 orders (see the Logs for details).

1. Try multiple methods of wavelength calibration to capture order 3.   You may need to remove the WaveCalib_A_0_DET01.fits each time you run PypeIt and try another method.  Here are some additional things to consider:

- Use the Holy Grail algorithm to bootstrap a solution from the OH lines without any archive. Run it on the science-frame arc we just built. This both (a) tests whether `OH_NIRES` covers SOAR's range for orders 3–7, and (b) produces an initial per-order solution we can promote into a new `soar_triplespec.fits` reid archive.
- Verify the line list. Compare `OH_NIRES` vs `OH_triplespec` coverage against SOAR's order wavelengths; pick whichever matches the detected OH lines better.

2. Using the solution you found in your previous effort, generate a new reidentify file for the soar_triplespec spectrograph.  Check to see if you can then increase the cc_thresh parameter back to its default.

### Frame typing

PypeIt uses images of a range of frame types to calibrate and reduce the data.  Two of these are "dark" and "lampoffflats".  These are used to calibrate the detector and remove the thermal emission from the telescope and dome.   

1. The test data set includes an image that is named with DFLAT but has a 'target' of DFLAT_OFF.  Currently we are treating it as a "dark" frmae.  Please:

- Consider if we should type it as a "lampoffflats" frame instead
- Remind me what will be different in the image processing if we have a "lampoffflats" frame instead of a "dark" frame
- Do not modify any code (yet)
- Log your findings below.
- Include your recommendations

2. I have modified the PypeIt file to type the DFLAT_OFF image as a "lampoffflats" frame.  Please:

- Confirm that by default the "lampoffflats" frame is not used for calibration.
- Modify the soar_triplespec.rst doc to describe how to use the "lampoffflats" frame for calibration, if the user desires
- Run PypeIt again and inspect the QA files.  
- Propose any other additional changes to the code or docs 
- Log your findings below.

## Docs

Here are some guidelines when generating the docs:

- Use .rst files
- Use :ref: tags to link to other pages (instead of paths and :doc:)  

1. Generate a soar_triplespec.rst doc in doc/spectrographs.  Model it after the other files in that folder.  Update any other files in doc/

2. Generate a tutorial for soar_triplespec in doc/tutorials.  Model it after the other files in that folder.  Use outputs in the 
PypeIt-development-suite/pypeitdev/soar_tspec/TripleSpec/ folder as needed.  As always, Log your work.  Run make_html to generate the docs and inspect the outputs.

3. Make the following modifications for the new docs:

- Use :ref: instead of :doc: for internal links.  Do so in soar_triplespec_howto.rst and soar_triplespec.rst.  Replace any instances of :doc: with :ref:.

## Finishing up

1. Follow standard PypeIt development procedure to indicate the changes and double check we have included proper doc strings for any new methods.

2. Set up for Dev-Suite testing by following what already exists for p200_tspec:

- Copy the raw files necessary for the test from PypeIt-development-suite/pypeitdev/soar_tspec/TripleSpec/ to 
PypeIt-development-suite/RAW_DATA/soar_tspec/TSPEC/
- Copy the PypeIt file into PypeIt-development-suite/pypeit_files/ and rename it accorinding the naming convention there
- Edit the test files so that this test is run automatically

3. From your work on this spectrograph, update any of the Skills as you see fit.

## Prompts

1. Read this doc.  Perform the first task under Prep.
2. Read this doc.  Perform the 2nd task under Prep.
3. Read this doc.  Perform the 3rd task under Prep
4. Read this doc.  Perform the 1st task under Development/Wavelength Calibration
5. Read this doc.  Perform the 2nd task under Development/Wavelength Calibration
6. Read this doc.  Remove all the files in Calibrations/ and do a fresh run.  Inspect the outputs and report your findings in the Logs.
7. Read this doc.  Execute the first request in the Docs Section.
8. Read this doc.  Execute the 2nd request in the Docs Section.
9. Read this doc.  Execute the 1st request in the Finishing up Section.
10. Read this doc.  Execute the 2nd request in the Finishing up Section.
11. Read this doc.  Execute the 3rd request in the Finishing up Section.
12. Read this doc.  Execute the 3rd request in the Docs Section.
13. Read this doc.  Perform the 1st task under Development/Frame typing.
14. Read this doc.  Perform the 2nd task under Development/Frame typing.

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-05-30 (First end-to-end run of `soar_tspec`; pipeline reaches wavelength calibration then fails on the P200 reid archive)

**Goal:** Prep task 1 — run PypeIt on `soar_tspec_A.pypeit`, inspect QA, assess the current state, and recommend next steps.

**Environment / setup**
- Ran in the `pypeit14` conda env, which uses the editable install at `PypeIt/pypeit` on branch `soar_triplespec`.
- Verified `soar_tspec` is registered in `pypeit/spectrographs/__init__.py` and `SOARTelescopePar` exists in `telescopes.py`.
- Working dir: `PypeIt-development-suite/pypeitdev/soar_tspec/TripleSpec/`.

**Command used**
```bash
conda activate pypeit14
run_pypeit soar_tspec_A.pypeit -o > run_log.txt 2>&1
```

**Skills:** No skill tool was formally invoked for this run, but I followed the approach of the `diagnose-reduction` skill (read `run_log.txt`, inspect `Calibrations/`, trace the failure to its calibration step). The `wavelength-calibration` skill is the relevant one for the recommended next stage.

**Code fixes required just to run** (both are stale `msgs`→`log`/`PypeItError` migration leftovers, not algorithmic changes):
1. `PypeIt/pypeit/spectrographs/soar_tspec.py` — the file was written against the old `from pypeit import msgs` API. Changed to `from pypeit import log` + `from pypeit import PypeItError`; `msgs.error(...)` → `raise PypeItError(...)`, `msgs.info(...)` → `log.info(...)`.
2. `PypeIt/pypeit/edgetrace.py` (lines 901–903) — **pre-existing core bug on this branch**: three `msgs.info(...)` calls remained although the module no longer imports `msgs`. Changed to `log.info(...)`. This is unrelated to SOAR and would break edge tracing for *any* fixed-format echelle on this branch.

**What worked (assessment of current state)**
- Frame typing from the `.pypeit` file is correct (arc/tilt, dark, pixelflat/trace, science, standard).
- Calibrations built successfully: `Arc`, `Dark`, `Tiltimg`, `Edges`, `Slits` all written to `Calibrations/`.
- **Slit/order tracing is in good shape.** All 5 orders are found; `order_spat_pos` and `orders` (7,6,5,4,3) are well tuned — `match_order` reports order-matching separations of ≈ 0.00 for every order (max |sep| = 0.003). The `auto_pca=False` + polynomial-fit + `det_buffer`/`minimum_slit_length` choices in `default_pypeit_par` are working.
- Arc central spectra extracted for all 5 slits; measured arc-line FWHM ≈ 2.0–2.2 px.

**What failed (the single blocker)**
- **Wavelength calibration fails for every order.** `default_pypeit_par` reuses the Palomar archive `reid_arxiv = 'p200_triplespec.fits'`. Reidentification cross-correlation is below threshold on all orders:
  - order 7: shift=-260.6, cc=0.247 → failed
  - order 6: shift=-0.8, cc=0.223 → failed
  - order 5: shift=+340.9, cc=0.276 → failed
  - order 4: shift=-399.3, cc=0.237 → failed
  - order 3: shift=-282.3, cc=0.247 → failed
- Result: 0 of 5 orders calibrated → `wavecalib.py: "No successful 2D Wavelength fits. Cannot proceed."` The run aborts here, so **no Science output and no QA HTML are produced** (only the calibration FITS exist; QA dir is empty).
- Root cause: SOAR/TripleSpec ≠ Palomar/TripleSpec. They share a name and use OH sky lines, but the dispersion/spectral coverage and detector geometry differ enough that the P200 template does not register. SOAR needs its **own** wavelength reidentification archive.

**Recommendations for the next stage**
1. **Build a SOAR-specific wavelength solution** (highest priority; unblocks everything downstream). Use the `wavelength-calibration` skill / `pypeit_identify` on the SOAR arc (the OH sky lines in the `SPEC_ARC` frames), order by order, to construct a new reid archive `soar_triplespec.fits` and point `reid_arxiv` at it. The arc-line FWHM (~2.1 px) and order trace are already reliable, so identification should be tractable.
   - As an interim sanity check, try `method='holy-grail'` (full-template-free) on a single order to confirm the OH line list (`OH_NIRES`) matches the SOAR coverage before committing to a template.
2. **Confirm `OH_NIRES` is the right line list** for SOAR's wavelength range / orders 3–7; swap to `OH_triplespec` if its coverage matches better.
3. After wavelength calib works, re-run to generate tilts QA, science extraction, and the full QA HTML, then assess sky subtraction / extraction and the standard-star / sensfunc path.
4. **Upstream cleanup:** the two `msgs` fixes above should be folded into the branch; the `edgetrace.py` one in particular is a genuine core regression worth flagging in the PR.
5. Once a reduction completes end-to-end, register the setup in the dev suite (`add-devsuite-setup` skill) so it becomes a regression test.

**Files touched this run:** `soar_tspec.py`, `edgetrace.py` (code fixes); `run_log.txt` (output, git-ignored).

### 2026-05-30 (Prep task 2 — science frames as the arc source; same wavecalib failure, confirms the P200 archive is the sole blocker)

**Goal:** Prep task 2 — the `.pypeit` file was modified (by JXP) to drive wavelength calibration off the science frames instead of the dedicated arc, then re-run and inspect QA.

**`.pypeit` change (made by user):** the dedicated `SPEC_ARC_..._0103.fits` line is now commented out, and the two 300 s science frames are re-typed `arc,tilt,science`. Rationale (correct for NIR): a 2 s arc is sparse, whereas a 300 s science exposure is rich in OH sky emission lines — the natural wavelength reference for the near-IR.

**Command used**
```bash
conda activate pypeit14
run_pypeit soar_tspec_A.pypeit -o > run_log2.txt 2>&1
```

**What changed vs. run 1**
- The frame-typing change took effect: `Calibrations/Arc_A_0_DET01.fits` was rebuilt from the science frames, central arc spectra extracted for all 5 slits, OH lines detected, measured FWHM ≈ 2.0–2.2 px. So the OH-from-science approach is mechanically sound — the input spectrum is fine.

**What did NOT change (the blocker is unmoved)**
- Reidentification still fails on **every** order against `reid_arxiv = 'p200_triplespec.fits'`, with cross-correlation values **identical to run 1** to the decimal (e.g. order 4: shift=-399.261, cc=0.237; order 3: shift=-282.251, cc=0.247). The identical numbers are themselves informative: the SOAR "arc" and "science" OH spectra are essentially the same OH-sky spectrum, so swapping the frame source cannot change the cross-correlation — the mismatch is entirely between SOAR's OH spectrum and the **P200 template**, not between arc and science frames.
- Same terminal error: `wavecalib.py: "No successful 2D Wavelength fits. Cannot proceed."` → no Science output, **QA dir still empty** (the crash precedes QA HTML/PNG generation).

**Conclusion:** Using the science frames for wavelength calibration is the right long-term choice (and should stay), but it does not unblock the reduction on its own. The single remaining blocker is unchanged from run 1: **there is no SOAR-appropriate reid archive.** No further `.pypeit`/frame-typing tweak will help until a SOAR wavelength solution exists.

> **⚠️ Correction (discovered in Prep task 3):** the analysis above is WRONG, and the bold-faced cc/shift "identical numbers" were misread. Run 2 never actually rebuilt the Arc from the science frames — it silently **reused the stale `Arc_A_0_DET01.fits` cached from run 1** (built from the dedicated 2 s arc). PypeIt reuses existing calibration FITS regardless of `-o` or `.pypeit` edits. The numbers were identical because it was literally the same cached file. The P200 archive is **not** the blocker; the stale arc was. See the task-3 entry below.

**Recommendations for the next stage** (this maps directly onto the new `## Development → Wavelength Calibration` section):
1. **Switch to the Holy-Grail algorithm** (`method = 'holy-grail'`, template-free) to bootstrap a solution from the OH lines without any archive. Run it on the science-frame arc we just built. This both (a) tests whether `OH_NIRES` covers SOAR's range for orders 3–7, and (b) produces an initial per-order solution we can promote into a new `soar_triplespec.fits` reid archive.
2. **Verify the line list.** Compare `OH_NIRES` vs `OH_triplespec` coverage against SOAR's order wavelengths; pick whichever matches the detected OH lines better.
3. Once Holy-Grail (or `pypeit_identify`) yields good per-order fits, archive them as `soar_triplespec.fits`, point `reid_arxiv` there, switch `method` back to `reidentify`, and re-run to finally generate tilts/extraction/QA.

**Files touched this run:** none (run only); `run_log2.txt` (output, git-ignored). The science-frame `.pypeit` edit was made by the user.

### 2026-05-30 (Prep task 3 — removed stale calibrations → FIRST full end-to-end reduction; the P200 archive works after all)

**Goal:** Prep task 3 — the user removed the cached `Arc`/`Tiltimg` calibration files; re-run, inspect QA, recommend.

**Command used**
```bash
conda activate pypeit14
run_pypeit soar_tspec_A.pypeit -o > run_log3.txt 2>&1   # EXIT CODE: 0  🎉
```

**The key discovery (overturns the runs 1–2 conclusion):**
With the cached `Arc` gone, PypeIt **rebuilt the Arc by combining the two science frames** (`...0116.fits` + `...0117.fits`) — the log shows `Combining frames: ... Preparing a Arc calibration frame ... process: SPEC_2025cy_..._0116/0117`. Reidentification against the **same** `p200_triplespec.fits` archive now succeeds for 4 of 5 orders:

| Order | shift (px) | cc | 1D RMS (px) |
|------:|-----------:|------:|------------:|
| 7 | −16.2 | 0.909 | 0.102 |
| 6 | −15.4 | 0.938 | 0.073 |
| 5 | −15.9 | 0.913 | 0.105 |
| 4 | −15.9 | 0.957 | 0.065 |
| 3 | −15.6 | 0.660 | — (failed, cc<thresh) |

Note the now-consistent ~−16 px shift across all orders (vs the nonsensical −260/−400 px scatter of runs 1–2). **So the P200 archive is appropriate for SOAR/TripleSpec after all** — the real blocker in runs 1–2 was a *stale cached Arc built from the sparse 2 s dedicated arc frame*, which PypeIt reused because `-o` does **not** force calibration rebuilds and the `.pypeit` frame-typing edit alone didn't invalidate the cache. Deleting `Calibrations/` was the fix.

**Lesson (important for the workflow):** after changing frame typing / calibration grouping in a `.pypeit` file, you must delete the affected files in `Calibrations/` (or use `pypeit_clean`/`--remove_calibs`); `run_pypeit -o` only overwrites *science* products, not calibrations.

**Full pipeline now runs end-to-end.** Products written:
- `Calibrations/`: `Arc`, `Tiltimg`, `WaveCalib`, `Edges`, `Slits`, `Dark`, plus `Flat` (illum/pixel) — all present.
- `Science/`: `spec2d` + `spec1d` (FITS + .txt) for both science frames (0116, 0117) and both standards (0120, 0121).
- `QA/`: full `MF_A.html` + per-exposure HTML and ~60 PNGs (Arc 1D/2D fits, FWHM, tilts, object profiles/traces, illum-flat fine-corr).

**QA assessment (inspected the PNGs):**
- **2D wavelength solution** (`Arc_2dfit_global`): clean, 4 orders spanning ≈ 8000–18700 Å (0.8–1.87 µm; z/J/H into K), `norder_coeff=6, nspec_coeff=4`, **RMS = 0.107 Ang·Order#**. Order 3 (reddest, K-band) is absent.
- **1D fits** (e.g. order 7): dozens of OH lines identified, linear solution, RMS 0.10 px, residuals within ±0.3 px. Flux correctly restricted to pixels ≳1024 (matches `spec_min_max` for order 7).
- **Tilts**: well-behaved, fit RMS 0.018–0.057 px (RMS/FWHM 0.008–0.027) across orders.
- **Object finding / extraction**: target 2025cy found as a single clean object at **SNR ≈ 330** (order 7); sky subtraction and boxcar/optimal extraction succeed on all orders.

**Remaining issues / open questions:**
1. **Order 3 (K-band) wavelength calibration fails** (cc=0.66 < threshold; `redo_echelle_orders` couldn't recover it). It is currently filled only by the global 2D fit / extrapolation and flagged bad. This is the natural target for Development task 2 (line-list coverage).
2. Need to verify the absolute wavelength solution against known OH lines (esp. the extrapolated red end) and check flexure (currently `spec_method='skip'`).
3. The standards reduced fine, but the sensfunc/telluric/flux path is untested.

**Recommendations for the next stage** (the `## Development → Wavelength Calibration` section is now partly moot — the archive works — so reprioritize):
1. **Down-weight the "build a new archive / Holy-Grail" task.** The P200 archive already gives RMS ~0.1 px on 4/5 orders. A native SOAR archive is now a *nice-to-have* (mainly to capture order 3 and remove the P200 dependence), not a blocker.
2. **Focus on order 3 (Development task 2):** compare `OH_NIRES` vs `OH_triplespec` line-list coverage in the K-band order-3 range; the reddest order likely has few/poorly-matched OH lines. Try `OH_triplespec`, and/or relax `rms_thresh_frac_fwhm` / `cc_thresh` for that order, or run `pypeit_identify` on order 3 to anchor it and fold the result into a SOAR archive.
3. **Bake the OH-from-science requirement into the spectrograph defaults / docs** so users don't hit the stale-arc trap — e.g. document that the dedicated 2 s arc is not used and OH sky lines from science frames drive wavecalib.
4. Once order 3 is solid, build `soar_triplespec.fits` from this run's `WaveCalib` (covers all 5 orders), then proceed to sensfunc/telluric and register the setup as a dev-suite regression test (`add-devsuite-setup`).

**Files touched this run:** none (run only); `run_log3.txt`, `Calibrations/`, `Science/`, `QA/` (all git-ignored outputs).

### 2026-05-30 (Development/Wavelength task 1 — recovered order 3; settled on `cc_thresh=0.60` with OH_NIRES)

**Goal:** Development → Wavelength Calibration task 1 — try multiple wavelength-calibration methods to capture the one failing order (order 3, K band), removing `WaveCalib_A_0_DET01.fits` between runs.

**Diagnosis first.** Loaded the run-3 `WaveCalib` and confirmed order 3 had **no solution** (masked). The reason: order 3 reidentifies against the P200 arxiv with **cc = 0.660**, sitting just under the default **`cc_thresh = 0.70`**, so it is dropped before fitting. The per-order shift (−15.6 px) is consistent with the other orders, i.e. the cross-correlation is geometrically fine — only the correlation *strength* is marginal (K band has higher thermal background and a sparser OH forest). Also checked the line lists: `OH_NIRES` has **118** lines redward of 19000 Å (out to **24980 Å**), vs `OH_triplespec`'s **67** (out to **23535 Å**) — so the line list is *not* the limitation for order 3; OH_NIRES covers it better.

**Method:** in-place experiments via a `[calibrations][[wavelengths]]` block in the `.pypeit`, deleting `WaveCalib_*` before each `run_pypeit -o`. (`.pypeit` backed up to `.bak` and restored at the end; the working WaveCalib was saved to `/tmp` first.)

| Exp | Change | Result |
|----|--------|--------|
| **1** | `cc_thresh = 0.60` (keep reidentify + P200 + OH_NIRES) | ✅ **All 5 orders.** Order 3: rms **0.108 px**, **78** lines, 18806–24681 Å. Bluer orders unchanged. |
| **2** | `method = holy-grail` (no arxiv), OH_NIRES | ❌ **Crashes** in `wavecalib.py:redo_echelle_orders` → `AttributeError: 'BuildWaveCalib' object has no attribute 'arcspec_arxiv'` (the bad-order-redo path assumes a reid arxiv exists). Even before the crash, per-order fits were erratic (good orders ~0.10 px, others 0.88–3.1 px). Not viable. |
| **3** | `cc_thresh = 0.60` + `lamps = OH_triplespec` | ✅ All 5 orders, order 3 rms **0.054 px** but only **52** lines, and OH_triplespec stops at 23535 Å so the reddest ~1000 Å of order 3 is **extrapolated**, not anchored. |

**Decision — Exp 1.** `cc_thresh = 0.60` is the minimal, robust fix: it admits order 3 while every order's final reidentify RMS stays ~0.06–0.11 px, and OH_NIRES anchors order 3 across its *full* wavelength range with real lines (Exp 3's lower RMS is partly an artifact of fewer points + extrapolation). Holy-Grail is both buggy here and less reliable, so a bespoke SOAR arxiv is **not** needed — the P200 arxiv + relaxed threshold suffices.

**Code change (permanent):** added to `soar_tspec.py:default_pypeit_par` (with an explanatory inline comment), so users get this for free:
```python
par['calibrations']['wavelengths']['cc_thresh'] = 0.60
```
Kept `lamps = ['OH_NIRES']` and `method = 'reidentify'` as-is.

**Verification (final clean run with the new default, no `.pypeit` overrides):** `run_pypeit` exits 0; **all 5 orders solved**, no "Masking … bad orders (RMS)" warning. Per-order: ord7 0.102 px (55 ln), ord6 0.073 px (102), ord5 0.105 px (95), ord4 0.065 px (100), ord3 0.108 px (78). 2D global fit now spans **~8000–24700 Å (0.80–2.47 µm)** across all 5 orders, RMS 0.123 Ang·Order# (`Arc_2dfit_global_A_0_DET01.png`). Full Science + QA regenerated.

**Side finding (worth an upstream issue/PR):** `redo_echelle_orders` is not safe for `method='holy-grail'` — it unconditionally references `self.arcspec_arxiv`, which only exists on the reidentify path. Any echelle run that uses holy-grail and has a bad order will crash there.

**Recommendations for the next stage:**
1. Wavelength calibration is now complete for all 5 orders — close this Development task.
2. Sanity-check the **absolute** wavelength scale on a few known OH lines (esp. order 3's red end) using `pypeit_chk_wavecalib` / the 1D-fit QA before trusting the K-band tail.
3. Proceed down the afterburn path: sensfunc from HD 116699 (A0 V telluric standard), telluric correction (`TellPCA_3000_26000`), then 1D coadd of the ABBA science pairs.
4. Optionally build a native `soar_triplespec.fits` arxiv from this 5-order `WaveCalib` to drop the P200 dependence (low priority — current setup works).
5. Register the setup as a dev-suite regression test (`add-devsuite-setup`).

**Files touched this run:** `PypeIt/pypeit/spectrographs/soar_tspec.py` (added `cc_thresh=0.60` + comment); `run_exp1.txt`/`run_exp2.txt`/`run_exp3.txt`/`run_final.txt`, `Calibrations/`, `Science/`, `QA/` (git-ignored outputs). The `.pypeit` was modified during experiments and **restored** to its task-3 state.

### 2026-05-30 (Development/Wavelength task 2 — built native `soar_triplespec.fits` arxiv; restored `cc_thresh` to default)

**Goal:** Development → Wavelength Calibration task 2 — promote the 5-order solution from task 1 into a native SOAR reid arxiv, then check whether `cc_thresh` can go back to its default (0.70).

**Skill:** invoked the `wavelength-calibration` skill. It documents `pypeit_compile_wvarxiv`, but that tool ingests per-order `wvarxiv` files emitted by `pypeit_identify`; I instead already had a complete automatic `WaveCalib_A_0_DET01.fits`, so I built the arxiv directly from it (same end format).

**How the arxiv was built.** Inspected the existing `p200_triplespec.fits`: it is a single `BinTableHDU` with one row per echelle order and columns `wave[2048]`, `flux[2048]`, `order`. My `WaveCalib` carries exactly these per order (`wave_soln`, `spec`, `ech_order`), all length-2048 and monotonic. I wrote a 5-row table (orders 7,6,5,4,3 — blue→red, matching the P200 ordering) to:
```
PypeIt/pypeit/data/arc_lines/reid_arxiv/soar_triplespec.fits
```
(`flux` = the per-order OH-sky arc spectrum, `wave` = the task-1 wavelength solution.)

**Spectrograph change.** In `soar_tspec.py:default_pypeit_par`:
- `reid_arxiv`: `'p200_triplespec.fits'` → `'soar_triplespec.fits'`
- **removed** the `cc_thresh = 0.60` override (back to the default 0.70); replaced the comment block with one documenting the arxiv's bootstrap provenance.

**Verification (final run, default params, no `.pypeit` overrides):** `run_pypeit` exits 0. Every order now self-reidentifies essentially perfectly against the native arxiv:

| Order | shift (px) | cc | 1D RMS (px) | nlines |
|------:|-----------:|----:|------------:|-------:|
| 7 | 0.000 | 1.000 | 0.099 | 55 |
| 6 | 0.000 | 1.000 | 0.070 | 101 |
| 5 | −0.000 | 1.000 | 0.105 | 95 |
| 4 | −0.000 | 1.000 | 0.065 | 100 |
| 3 | −0.000 | 1.000 | 0.108 | 78 |

No "Masking … bad orders" warning; 2D global fit RMS 0.547 Ang·Order#, coverage ~8000–24700 Å. **So yes — with the native arxiv, `cc_thresh` returns to its default and all 5 orders (including the previously-failing order 3) pass.** Confirmed `WavelengthSolutionPar` default `cc_thresh = 0.7`.

**Important caveat (honest assessment).** The arxiv was built from *this exact* arc, so cc = 1.000 / shift = 0.000 is a **self-reference** — it proves the arxiv format/wiring is correct and that the default threshold is no longer the limiter, but it is **not** an independent validation. The real test is a *different* SOAR/TripleSpec dataset (another night / target), where the proper checks are a high-but-not-unity cc, a small consistent shift, and per-order RMS ~0.1 px. Until then, treat the restored default `cc_thresh` as validated only for setup A.

**Recommendations for the next stage:**
1. Re-reduce an **independent** SOAR/TripleSpec dataset against `soar_triplespec.fits` to confirm cross-dataset reidentification and lock in the default `cc_thresh`. (If a second night isn't available yet, keep an eye on cc when one is.)
2. Consider seeding the arxiv with **2+ exposures** (more rows / S-N) once more data exist, so it generalizes better than a single-arc template.
3. The arxiv is a binary file under `pypeit/data/arc_lines/reid_arxiv/` — make sure it is committed with the branch (it is *not* git-ignored there, unlike dev-suite outputs) and noted in the PR.
4. Proceed to the afterburn path (sensfunc → telluric → coadd) now that wavecal is fully solved.
5. Register the dev-suite setup (`add-devsuite-setup`) + add to the changelog.

**Files touched this run:** `PypeIt/pypeit/spectrographs/soar_tspec.py` (reid_arxiv → native, cc_thresh override removed); **new** `PypeIt/pypeit/data/arc_lines/reid_arxiv/soar_triplespec.fits` (committed data file); `run_arxiv.txt`, `Calibrations/`, `Science/`, `QA/` (git-ignored outputs).

### 2026-05-30 (Prompt 6 — clean-slate validation: wiped Calibrations/ and re-reduced from scratch)

**Goal:** delete everything in `Calibrations/` and run fresh, to prove the full pipeline rebuilds end-to-end from raw frames with the current code (native `soar_triplespec.fits` arxiv, default `cc_thresh`, no `.pypeit` overrides).

**Command used**
```bash
rm -rf Calibrations/* Science/* QA/*    # full clean slate (also cleared stale Science/QA)
conda activate pypeit14
run_pypeit soar_tspec_A.pypeit -o > run_fresh.txt 2>&1   # EXIT 0, no ERROR/WARNING
```

**Result — clean reduction, no errors or warnings.** Every calibration was rebuilt from the raw frames: `Arc, Dark, Edges, Flat, Slits, Tiltimg, Tilts, WaveCalib`. All 5 orders reidentified against the native arxiv at **cc = 1.000, shift = 0.000**; 2D fit RMS 0.547 Ang·Order#; tilt RMS 0.018–0.057 px. No "bad orders" masking. Products: `spec1d`+`spec2d` for both science (0116, 0117) and both standards (0120, 0121); QA = 5 HTML + 72 PNGs.

**Extraction quality (spec1d for 0116).** The science target 2025cy (OBJ0697, spat_frac 0.70) is extracted on all 5 orders with per-order `wv_rms` matching the WaveCalib (0.099/0.070/0.105/0.065/0.108 px):

| Order | S/N | opt_fwhm (px) |
|------:|----:|--------------:|
| 7 | 28.3 | 1.03 |
| 6 | 30.3 | 0.99 |
| 5 | 16.6 | 0.95 |
| 4 | 21.0 | 0.89 |
| 3 | 4.0 | 0.59 |

Order 3 (K band) has the lowest S/N (4.0), as expected for the faintest/highest-background order.

**Minor finding.** A second object OBJ0432 (spat_frac 0.43) is also extracted on every order but is essentially **noise** (S/N from −1.8 to +5.5, `opt_fwhm` pinned at the 4.291 px ceiling on several orders). It appears only because `default_pypeit_par` sets `reduce.findobj.maxnumber_sci = 2`. The inline comment there actually says "Slit is narrow so allow one object per order" — which argues for **1**, not 2. Recommend revisiting `maxnumber_sci` (set to 1, or raise the SNR threshold) so spurious second objects aren't carried through to coadd/flux. Not a blocker.

**Conclusion.** The SOAR/TripleSpec reduction is now fully reproducible from a clean slate end-to-end. Wavelength calibration (the long-standing blocker) is solved for all 5 orders. The pipeline is ready for the afterburn path.

**Recommendations for the next stage:**
1. Sensfunc from HD 116699 (A0 V) → telluric (`TellPCA_3000_26000`) → 1D coadd of the ABBA science pairs.
2. Tidy `maxnumber_sci` (see minor finding) to avoid spurious objects.
3. Validate the native arxiv on an **independent** SOAR dataset (still the one outstanding wavelength caveat from the previous entry).
4. Commit the staged code/data changes (`soar_tspec.py`, `edgetrace.py`, `reid_arxiv/soar_triplespec.fits`), update `CHANGES.rst`, and register a dev-suite setup.

**Files touched this run:** none (run only); `run_fresh.txt`, `Calibrations/`, `Science/`, `QA/` (git-ignored outputs).

### 2026-05-30 (Docs task 1 — authored `doc/spectrographs/soar_triplespec.rst` and refreshed the auto-generated doc tables)

**Goal:** Docs task 1 — create `soar_triplespec.rst` under `doc/spectrographs/`, modeled on the existing files, and update any other affected `doc/` files.

**Model chosen.** Surveyed `doc/spectrographs/`; `keck_nires.rst` is the closest analog (fixed-format NIR cross-dispersed echelle that wavelength-calibrates off OH sky lines in the science frames), with `soar_goodman.rst` for the SOAR-house style. Wrote the new page in the same structure: *Overview → PypeIt File (setup command + data-block example) → Calibrations (Wavelength Calibration, Flat Fielding) → Additional Reading.*

**New file:** `PypeIt/doc/spectrographs/soar_triplespec.rst`. Content captures everything learned in the earlier log entries:
- Echelle pipeline, 5 fixed orders (7→3), ~0.8–2.47 µm.
- `pypeit_setup -s soar_tspec -b -c A` plus a representative data block.
- Frame typing (dome-on = `pixelflat,trace`, dome-off = `dark`, science = `science`).
- **Wavelength calibration**: OH sky lines from the *science* frames (type them `arc,tilt`, omit the short arc); native `soar_triplespec.fits` arxiv + `OH_NIRES`; why OH_NIRES anchors order 3.
- A `.. warning::` documenting the stale-cached-calibration trap (the run-2 lesson): `run_pypeit -o` does not rebuild `Arc_*`/`WaveCalib_*`, so delete them after re-typing frames.

**Other `doc/` files updated:**
- `doc/spectrographs/spectrographs.rst` — added `soar_triplespec` to the instrument toctree (after `soar_goodman`). *(hand-edited)*
- Regenerated the three auto-built tables that key off the spectrograph registry, by running the same scripts the `apirst` Make target uses (no dev-suite data needed):
  - `scripts/build_spectbl_rst.py` → `include/spectrographs_table.rst` now has the `soar_tspec` row (Echelle, Supported=True).
  - `scripts/build_detector_table.py` → `include/inst_detector_table.rst` now has the `soar_tspec` detector row.
  - `scripts/build_par_rst.py` → `pypeit_par.rst` now has the "SOAR TSPEC (`soar_tspec`)" parameter section (reflects the `reid_arxiv = soar_triplespec.fits` default).

**Validation.** Parsed the new `.rst` with docutils — no structural errors (the only messages are "unknown role" for `:ref:`/`:doc:`/`:meth:`, which are Sphinx-only roles also used by `keck_nires.rst`; they resolve in a real build). Avoided two easy pitfalls: (a) dropped the `|micron|`/`|Angstrom|` substitutions since they are not defined in `include/links.rst`, using plain text; (b) made "Additional Reading" plain `:doc:` bullets rather than a `toctree`, to avoid "document already in a toctree" warnings.

**Not done (out of scope / needs full build):** the API stub `doc/api/pypeit.spectrographs.soar_tspec.rst` is produced by `sphinx-apidoc` during `make apirst`; it will appear on the next full doc build (`cd doc; make html`, which needs internet + `PYPEIT_DEV`). A complete `make html` was not run here.

**Recommendations for the next stage:**
1. On the next full `make html`, confirm `soar_triplespec` renders and the API stub is generated; fix any cross-ref warnings.
2. Add a CHANGES.rst / release-notes entry (the `update-changelog` skill) for the new spectrograph + arxiv.
3. Proceed to the afterburn path (sensfunc → telluric → coadd) and the dev-suite registration.

**Files touched this run:** `PypeIt/doc/spectrographs/soar_triplespec.rst` (new); `PypeIt/doc/spectrographs/spectrographs.rst`, `PypeIt/doc/include/spectrographs_table.rst`, `PypeIt/doc/include/inst_detector_table.rst`, `PypeIt/doc/pypeit_par.rst` (updated/regenerated).

### 2026-05-30 (Docs task 2 — authored the SOAR/TripleSpec tutorial and built the HTML)

**Goal:** Docs task 2 — write a tutorial in `doc/tutorials/` modeled on the existing ones, using the real reduction outputs in `pypeitdev/soar_tspec/TripleSpec/`, then run the HTML build and inspect.

**Model chosen.** `nires_howto.rst` (Keck/NIRES) — the same kind of fixed-format NIR echelle that wavelength-calibrates off OH sky lines. Mirrored its structure: *Overview → Setup (with data block + wavecal-frame and A-B notes) → Core Processing (Order Edges, Wavelength Calibration, Field Flattening, Object Finding) → Outputs (2D/1D) → Next Steps.*

**New file:** `PypeIt/doc/tutorials/soar_triplespec_howto.rst`, populated with **real numbers and products from the prompt-6 clean run** (not invented):
- The actual `pypeit_chk_wavecalib` table (5 orders, RMS 0.065–0.108 px, ~8100–24700 Å).
- The actual standard-star (HD 116699) `spec1d` `.txt` summary (S/N 110–150 per order).
- Real `pypeit_show_2dspec` / `pypeit_show_1dspec` / `pypeit_chk_*` commands with the real filenames.
- The OH-from-science-frames rationale, the stale-calibration `.. warning::`, and the A-B differencing setup.

**Figures.** Copied four auto-generated QA PNGs from `TripleSpec/QA/PNGs/` into the (git-tracked) `doc/figures/` with descriptive names, and embedded them:
- `soar_triplespec_arc1d_order3.png` (order-3 1D wave fit)
- `soar_triplespec_2dwave_global.png` (global 2D solution, all 5 orders)
- `soar_triplespec_tilts2d_order3.png` (order-3 2D tilts)
- `soar_triplespec_obj_prof.png` (2025cy detected in order 7, S/N ~28)

**Wiring:** added the tutorial to the `doc/tutorials/tutorials.rst` toctree (after Shane Kast), and linked it from `doc/spectrographs/soar_triplespec.rst` (Additional Reading).

**HTML build + inspection.** A full `make html` needs network (`apirst` does `wget`) and the dev suite (`examples`), so I ran the offline equivalent: regenerated the API stubs with `sphinx-apidoc` (creating the missing `api/pypeit.spectrographs.soar_tspec.rst`) and then `sphinx-build -b html . _build/html`.
- First pass surfaced **two warnings on my files**, both fixed:
  1. `undefined label: 'telluric'` → the correct label is `telluric_correction` (also confirmed `sensitivity_function`, `fluxing`, `coadd1d/2d` labels).
  2. `api/...soar_tspec.rst not in any toctree` → incremental `sphinx-apidoc` does not overwrite the existing `api/pypeit.spectrographs.rst` index (a real `make html` runs `make clean` first, which would regenerate it); added the `soar_tspec` entry to that index by hand to match.
- Second pass: **build exit 0, zero soar/triplespec warnings.** Verified the rendered pages (`spectrographs/soar_triplespec.html`, `tutorials/soar_triplespec_howto.html`) and that all four figures were copied into `_build/html/_images/`.

**Recommendations for the next stage:**
1. On a real `cd doc; make html` (clean build with network + `PYPEIT_DEV`), re-confirm — the api index regenerates automatically there.
2. Afterburn path (sensfunc → telluric → coadd) is the substantive next step; the tutorial's "Next Steps" section is the placeholder to expand once that's run.
3. Commit the docs + figures with the code/data changes; add a `CHANGES.rst` entry.

**Files touched this run:** `PypeIt/doc/tutorials/soar_triplespec_howto.rst` (new); `PypeIt/doc/tutorials/tutorials.rst`, `PypeIt/doc/spectrographs/soar_triplespec.rst`, `PypeIt/doc/api/pypeit.spectrographs.rst` (updated); `PypeIt/doc/figures/soar_triplespec_{arc1d_order3,2dwave_global,tilts2d_order3,obj_prof}.png` (new, git-tracked); plus the `sphinx-apidoc`-generated `api/*.rst` stubs and the untracked `doc/_build/` output.

### 2026-05-31 (Finishing up task 1 — release notes + docstring audit, standard PypeIt dev procedure)

**Goal:** Finishing-up task 1 — follow standard PypeIt development procedure to record the changes, and verify docstrings for any new methods.

**State at start.** The prior sessions' work was already committed to the `soar_triplespec` branch across three commits (`msgs`, `docs`, `mo docs`); the working tree was clean. So this step is the changelog + docstring audit, not new code.

**Docstring audit.** Ran an AST check over `pypeit/spectrographs/soar_tspec.py`: all 14 methods (`init_meta`, `compound_meta`, `configuration_keys`, `raw_header_cards`, `get_detector_par`, `default_pypeit_par`, `pypeit_file_keys`, `check_frame_type`, `bpm`, `norders`, `order_spat_pos`, `orders`, `spec_min_max`, `order_platescale`) have docstrings. No *new* methods were added during this development (the work was parameter defaults + a reid archive + docs), so nothing was undocumented.

**Changelog (used the `update-changelog` skill).** `CHANGES.rst` is deprecated; entries now go in `doc/releases/`. Added bullets to the in-development `doc/releases/2.1.0dev.rst`:
- *Instrument-specific Updates*: new `soar_tspec` echelle support (5 orders 7--3, ~0.8--2.47 um; OH-sky wavecal + native `soar_triplespec.fits` arxiv).
- *Documentation Updates*: new spectrograph page + tutorial (with `:ref:` links).
- To make the `:ref:`<soar_triplespec>`` link resolve, added a `.. _soar_triplespec:` label atop `doc/spectrographs/soar_triplespec.rst` (matching the `gemini_gmos` convention; the tutorial already had `.. _soar_triplespec_howto:`).

**Not added:** a *Bug Fixes* bullet for the `edgetrace.py` stale-`msgs` fix — the user declined that entry (treated as an internal branch regression rather than a user-facing release fix), so it was intentionally left out of the release notes.

**Verification.** `sphinx-build -b html` exits 0 with **zero soar/triplespec/2.1.0dev warnings** after adding the label.

**Remaining (standard procedure, not yet done):** the changelog edit to `doc/releases/2.1.0dev.rst` and the `soar_triplespec.rst` label are **uncommitted**; commit them with the branch. Then the usual PR steps (push branch, open PR against `develop`) plus dev-suite registration via `add-devsuite-setup`. The afterburn path (sensfunc -> telluric -> coadd) and independent-dataset wavelength validation remain open from earlier entries.

**Files touched this run:** `PypeIt/doc/releases/2.1.0dev.rst` (release notes), `PypeIt/doc/spectrographs/soar_triplespec.rst` (added ref label).


### 2026-05-31 (Finishing up task 2 — registered the SOAR/TripleSpec dev-suite reduce test; PASSED 1/1)

**Goal:** Finishing-up task 2 — set up dev-suite testing for `soar_tspec`, mirroring `p200_tspec`.

**Skill:** used the `add-devsuite-setup` skill.

**Path note.** `$PYPEIT_DEV` was stale (pointed at a missing `.../PypeIt-Codes/...`). The live dev suite is the working-tree repo `/home/xavier/Projects/PypeIt/PypeIt-development-suite`, with `RAW_DATA` symlinked to `/media/xavier/SamsungT7/RAW_DATA`. Set `PYPEIT_DEV` to the working-tree path for the test run.

**What "mirror p200_tspec" actually means.** Studied the existing `p200_tspec` wiring:
- raw data at `RAW_DATA/p200_tspec/TSPEC/`, input file `pypeit_files/p200_tspec_tspec.pypeit`;
- registered only in `test_scripts/setups.py` `all_setups` (NOT in `test_setups.py` directly — that module imports `all_setups` from `setups.py` and auto-derives `_reduce_setups` from it). The `supported_instruments` list referenced in the `test_setups.py` docstring no longer exists as code; instruments come from `all_setups` keys.
- `p200_tspec` has no `test_load_images.py` entry, so (mirroring it) I added none.
- Filename convention: `pypeit_file_name(instr, setup)` = `{instr}_{setup.lower()}.pypeit`, and the harness rewrites the `path` line to the real RAW_DATA dir at runtime (`fix_pypeit_file_directory`).

**Changes made (this repo, PypeIt-development-suite):**
1. **Raw data** — created `RAW_DATA/soar_tspec/TSPEC/` and copied the 6 frames the reduction uses (2 dome flats on/off, 2 science, 2 standard) from `pypeitdev/soar_tspec/TripleSpec/`. (The commented-out dedicated arc is not used.) These live on the external-drive symlink and go to Google Drive, not git.
2. **Input file** — wrote `pypeit_files/soar_tspec_tspec.pypeit` (modeled on the working file; dropped the commented arc line; set `path ../soar_tspec/TSPEC`).
3. **Registration** — added `'soar_tspec': ['TSPEC']` to `all_setups` in `test_scripts/setups.py` (alphabetically, after `soar_goodman_blue`). This auto-enrolls it in the `reduce` phase, so the test runs automatically.

**Verification.**
- `./pypeit_test list` now shows `soar_tspec / TSPEC`.
- `./pypeit_test reduce -s soar_tspec/TSPEC` → **PASSED 1/1** (~3:54 run, ~1.5 GB peak RAM, no error messages; full QA written). Output under `REDUX_OUT/soar_tspec/TSPEC/`.

**Files to commit (left to the user, who handles git):** `test_scripts/setups.py`, `pypeit_files/soar_tspec_tspec.pypeit`. Raw data → Google Drive. `RAW_DATA/`, `REDUX_OUT/`, and `pypeitdev/` scratch are git-ignored / not for commit.

**Recommendations / remaining:**
1. Upload the 6 staged raw frames to the dev-suite Google Drive so CI/other users get them.
2. Optional: add a `unit_tests/test_load_images.py` entry for `soar_tspec` (p200_tspec lacks one, but it is the "supported instrument" convention).
3. Finishing-up task 3 (update Skills) is the remaining prompt.


### 2026-05-31 (Finishing up task 3 — updated four Skills from lessons learned)

**Goal:** Finishing-up task 3 — update the `.claude/skills/` as I see fit based on the SOAR/TripleSpec work. Each edit below comes from a concrete problem hit during this development, not speculation.

1. **`diagnose-reduction`** — added a "Common causes" bullet on **stale cached
   calibrations**: `run_pypeit -o` overwrites only science products and reuses
   existing `Calibrations/*.fits` regardless of `-o` or `.pypeit` edits; delete
   the affected files after re-typing frames. Tell-tale = two runs with
   byte-identical calibration results after changing inputs. (This is exactly
   what made runs 1-2 look like an archive problem when it was a stale arc.)

2. **`wavelength-calibration`** — three additions:
   - NIR with no usable arc lamp → use OH sky lines from the *science* frames
     (type them `arc,tilt,science`, give standards the same `calib`).
   - A code recipe to build an echelle `reid_arxiv` **directly from a
     `WaveCalib`** (one `BinTableHDU` row per order: `wave`, `flux`, `order`),
     since `pypeit_compile_wvarxiv` only ingests `pypeit_identify` outputs;
     plus the note that a self-built arxiv self-reidentifies at cc~1 and needs
     independent validation.
   - Tips: lower `cc_thresh` (~0.6) to recover a single marginal reddest order;
     prefer a line list that covers the full range with real lines (OH_NIRES).

3. **`add-devsuite-setup`** — corrected two inaccuracies and added detail:
   - `all_setups` lives in **`test_scripts/setups.py`** (not `test_setups.py`,
     which imports it and auto-derives reduce coverage); there is no longer a
     `supported_instruments` list in code.
   - Documented the lower-cased `.pypeit` filename rule
     (`f'{instr}_{setup.lower()}.pypeit'`) and that the `path` line is rewritten
     at run time; noted `_raw_data_dirs` for non-standard layouts; noted the
     load-image unit test is not universal (mirror the analog instrument).

4. **`new-spectrograph`** — added a note to use the current
   `from pypeit import log` / `PypeItError` API (not the deprecated `msgs`),
   since copying an older module drags in `from pypeit import msgs`, which now
   fails at import (the original `soar_tspec.py` failure).

**Verification.** Edits are documentation-only (Markdown skill files); reviewed each for accuracy against the code I touched this project. No build needed.

**This completes all prompts in the doc (Prep 1-3, Wavelength 1-2, fresh run, Docs 1-2, Finishing-up 1-3).** Open follow-ups carried in earlier entries: independent-dataset wavelength validation, the afterburn path (sensfunc -> telluric -> coadd), uploading the staged raw frames to the dev-suite Google Drive, and committing the staged changes (left to the user).

**Files touched this run:** `.claude/skills/{diagnose-reduction,wavelength-calibration,add-devsuite-setup,new-spectrograph}/SKILL.md`.

### 2026-06-09 (Docs task 3 — converted all `:doc:` internal links to `:ref:` in the SOAR docs)

**Goal:** Docs task 3 — replace every `:doc:` internal link with `:ref:` in both `soar_triplespec_howto.rst` and `soar_triplespec.rst`.

**Method.** Grepped both files for `:doc:` (7 instances total), then for each target document found its top-level label (`grep '^.. _'`) so the `:ref:` points at a real label rather than a path. Mapping used:

| `:doc:` path | `:ref:` label |
|---|---|
| `../A-B_differencing` | `a-b_differencing` |
| `../tutorials/soar_triplespec_howto` | `soar_triplespec_howto` |
| `../calibrations/wave_calib` | `wave_calib` |
| `../calibrations/flat` | `flat` |
| `../spectrographs/soar_triplespec` | `soar_triplespec` |
| `../running` | `run-pypeit` (section label; `running.rst` has no doc-level label) |
| `../out_spec1D` | `spec-1d-output` |

**Edits.**
- `doc/spectrographs/soar_triplespec.rst` (4): the A-B_differencing inline link + the three Additional Reading bullets. Custom link text ("-- a full worked example reduction.") preserved; bare `:ref:`label`` renders the target's section title, matching the surrounding `:ref:` usage already in the file.
- `doc/tutorials/soar_triplespec_howto.rst` (3): the instrument-notes link, the `../running` link, and the `../out_spec1D` link.

**Verification.** `grep ':doc:'` on both files now returns nothing. All seven `:ref:` labels were confirmed to exist in their target `.rst` files. Edits are doc-only; no build run (the offline `sphinx-build` path was validated for these files in Docs task 2).

**Files touched this run:** `PypeIt/doc/spectrographs/soar_triplespec.rst`, `PypeIt/doc/tutorials/soar_triplespec_howto.rst`.

### 2026-06-09 (Frame-typing task 1 — analysis: `DFLAT_OFF` should be `lampoffflats`, not `dark` (no code changed))

**Goal:** Frame-typing task 1 — decide whether the `DFLAT_OFF` image (currently typed `dark`) ought to be a `lampoffflats` frame, explain how the processing differs, and log findings/recommendations. **No code modified** (as instructed).

**The frame in question.** `SPEC_DFLAT_04-02-2025_0088.fits`. Header: `OBSTYPE='LAMPFLAT'`, `OBJECT='DFLAT_OFF'`, `EXPTIME=2.0` s. It is the **lamp-off companion of the dome flat** `SPEC_DFLAT_04-02-2025_0042.fits` (`OBSTYPE='LAMPFLAT'`, `OBJECT='DFLAT_ON'`, `EXPTIME=2.0` s, typed `pixelflat,trace`). Same instrument config, same 2.0 s exposure — a matched dome-flat ON/OFF pair, which is the textbook definition of a lamp-off flat, **not** a dark.

**Current typing logic** (`soar_tspec.py:check_frame_type`):
- `dark` → `idname=='LAMPFLAT' & target=='DFLAT_OFF'`  ← the frame lands here now.
- `pixelflat,trace` → `idname=='LAMPFLAT' & target=='DFLAT_ON'`.
There is no `lampoffflats` branch at present.

**What the `dark` typing actually does today: nothing.** `default_pypeit_par` calls `reset_all_processimages_par(use_darkimage=False)`, and the generated `.par` shows `use_darkimage = False` in all 15 frame-process blocks and `True` in none. So PypeIt builds a `Dark` master from the `DFLAT_OFF` frame and then **never subtracts it** from any frame. The lamp-off exposure currently contributes nothing to the reduction.

**How `dark` vs `lampoffflats` differ in image processing:**

| | `dark` | `lampoffflats` |
|---|---|---|
| Master built | `Dark_*` (combined dark) | folded into the flat build |
| When/where applied | only if a frame's `process.use_darkimage=True`; subtracted from that frame **scaled by the exposure-time ratio** (dark current model) | subtracted **directly, 1:1 (no scaling)** from the lamp-on pixel/illum flat *before* normalization — `calibrations.py:897-913`, `pixel_flat = pixel_flat.sub(lampoff_flat)` |
| Purpose | remove dark-current accumulation (needs a long, dark exposure to be meaningful) | remove the **thermal/dome glow + dark + bias structure imprinted in the flat itself** at the flat's exposure time |
| Effect here | none (`use_darkimage=False`) | the 2.0 s lamp-off is subtracted from the 2.0 s lamp-on flat — exposure-matched, so the direct subtraction is **exact** |

In short: a 2 s frame is a useless *dark* (captures essentially no dark current and is the wrong tool for scaling onto 300 s science), but it is an ideal *lamp-off flat* — it records exactly the additive background (telescope/dome thermal + bias + dark) present in the 2 s lamp-on dome flat, and `lampoffflats` removes that background from the flat field at matched exposure with no scaling.

**Precedent in the codebase.** `keck_nires` — the closest analog (fixed-format NIR echelle, OH-sky wavecal) — types its lamp-off dome flats as `lampoffflats` and explicitly returns **no** `dark` frames (`keck_nires.py:457-466`). `lbt_luci`, `keck_mosfire`, and `vlt_xshooter` follow the same pattern. SOAR/TripleSpec fits this NIR paradigm exactly.

**Recommendation — yes, retype `DFLAT_OFF` as `lampoffflats`.** It is the physically correct type and it actually *uses* the frame (today it's dead weight). When we make the change (next task), it will require two coordinated edits, **plus** confirming the calib grouping:
1. `soar_tspec.py:check_frame_type`: replace the `dark` branch with a `lampoffflats` branch (`idname=='LAMPFLAT' & target=='DFLAT_OFF'`); the `dark` branch can then return all-`False` (mirroring NIRES).
2. The `.pypeit` data block: change the `0088` frametype from `dark` → `lampoffflats`.
3. **Calib grouping:** the lamp-off must share the pixelflat's `calib` ID so `find_frame_files('lampoffflats', calib_ID=...)` matches it — both are `calib = 0` here, so this already holds.
4. **Dev-suite knock-on:** the registered test file `PypeIt-development-suite/pypeit_files/soar_tspec_tspec.pypeit` types this frame `dark` too; it must be updated in lockstep so the regression test stays consistent.

**Minor caveats (not blockers).** Only one `DFLAT_OFF` exists, so there's no combine/S-N benefit; subtracting a single lamp-off adds a little read noise to the flat, but it removes a systematic background and is the standard NIR approach — worth it. `darkframe` defaults (`exprng=[0,None]`) become moot once nothing is typed `dark`.

**Files touched this run:** none (analysis only). Inspected `soar_tspec.py`, `calibrations.py`, `keck_nires.py`, the `.pypeit`/`.par`/`.log`, and the two DFLAT FITS headers.

### 2026-06-11 (Frame-typing task 2 — empirically confirmed `lampoffflats` IS subtracted by default; effect negligible; kept as lampoffflats)

**Goal:** Frame-typing task 2 — the `check_frame_type` typing of `DFLAT_OFF` was changed (by JXP) from `dark` to `lampoffflats`. Confirm whether the `lampoffflats` frame is used by default, run PypeIt and inspect, and propose changes.

**Setup observed at start.** `soar_tspec.py:check_frame_type` now auto-types `DFLAT_OFF` (`OBSTYPE='LAMPFLAT'`, `target='DFLAT_OFF'`) as `lampoffflats` (was `dark`). The working `soar_tspec_A.pypeit` data block still typed `0088` as `dark` (explicit frametype overrides auto-typing), so I retyped that row to `lampoffflats` to exercise the new path.

**Headline finding — the lamp-off flat IS used by default (premise was incorrect).** A fresh `pypeit_setup` auto-types `0088` as `lampoffflats` and assigns it `calib = 0`, the **same calibration group as the `pixelflat,trace`**; the `.calib` association file lists `lampoffflats → Flat_A_0_DET01.fits`. A clean `run_pypeit` (deleted `Calibrations/*`, exit 0) then logged the subtraction in **two** places:
- `calibrations.py:get_slits:1096 — Subtracting lamp off flats using files:` (trace/edges image)
- `calibrations.py:get_flats:904 — Subtracting lamp off flats using files:` (pixel flat)

So PypeIt subtracts a `lampoffflats` frame from the trace and pixelflat **unconditionally** whenever it is present in their calib group. **There is no enable/disable parameter** (grep of the codebase: no `use_lampoff`/equivalent; the subtraction in `calibrations.py` is gated only on `len(raw_lampoff_files) > 0`). This matches the core docs (`doc/calibrations/flat_fielding.rst`: "The `lampoffflats` frames are *always* subtracted … If distinct frames are desired … we currently advise users to simply not use the lampofflats") and `keck_nires.rst`. The only lever is the **frametype** (`lampoffflats` → used; `dark`/omitted → not used).

**Impact is negligible (per-order science S/N, frame 0116, vs the prompt-6 run with no lamp-off subtraction):**

| Order | S/N (no lamp-off) | S/N (lamp-off) | wv_rms |
|------:|------------------:|---------------:|-------:|
| 7 | 28.3 | 28.3 | 0.099 |
| 6 | 30.3 | 30.2 | 0.070 |
| 5 | 16.6 | 16.6 | 0.105 |
| 4 | 21.0 | 21.0 | 0.065 |
| 3 | 4.0 | 4.4 | 0.108 |

All 5 orders still solved, full Science+QA written, no new errors. The lamp-off subtraction changes the final products at the noise level only — consistent with "users report no significant improvement."

**Decision (JXP).** Keep `DFLAT_OFF` typed as `lampoffflats` even though it is therefore subtracted by default; the negligible effect is acceptable. JXP will update the docs to discuss this. (I initially proposed reverting the auto-typing to `dark` so the example would match the existing doc's "not subtracted by default" wording, and asked which default was wanted — JXP chose to keep `lampoffflats`.) **No code was changed this run**; `check_frame_type` keeps the `lampoffflats` branch. The example `soar_tspec_A.pypeit` now types `0088` as `lampoffflats` to match the auto-typing.

**Note for the docs update (JXP is handling):** the current `doc/spectrographs/soar_triplespec.rst` "Flat Fielding" section says *"By default, the lamp-off flats are **not** subtracted as users report not significant improvement."* — this is now **inaccurate**: with the `lampoffflats` typing they **are** subtracted by default (negligibly). The accurate statement is that the lamp-off flats are subtracted from the trace/pixelflat by default, with negligible effect, and a user who wants to skip the subtraction should re-type the frame as `dark` (or omit it), since there is no on/off parameter.

**Other proposed changes (optional, not made):**
1. If skipping the lamp-off without retyping is ever wanted, that would require a *core* PypeIt feature (an opt-out flag on `lampoffflatsframe`); out of scope here.
2. The dev-suite raw set / `pypeit_files/soar_tspec_tspec.pypeit` types this frame `dark`; if we want the regression test to exercise the lamp-off path, retype it `lampoffflats` there too (the reduction passes either way given the negligible effect).

**Files touched this run:** `PypeIt-development-suite/pypeitdev/soar_tspec/TripleSpec/soar_tspec_A.pypeit` (typed `0088` as `lampoffflats`); `run_lampoff.txt` + regenerated `Calibrations/`/`Science/`/`QA/` (git-ignored). No PypeIt source or doc files modified (JXP owns the doc edit).
