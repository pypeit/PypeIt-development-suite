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

The current wavelength calibration fails because the P200 reid archive is not appropriate for SOAR/Triplespec.  We need to build a new reid archive for SOAR/Triplespec.

1. Try using the Holy Grail algorithm to build a wavelength solution.  You may try 

2. Confirm OH_NIRES is the right line list for SOAR's wavelength range / orders 3–7; swap to OH_triplespec if its coverage matches better.


## Prompts

1. Read this doc.  Perform the first task under Prep.
2. Read this doc.  Perform the 2nd task under Prep.
3. Read this doc.  Perform the 3rd task under Prep

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