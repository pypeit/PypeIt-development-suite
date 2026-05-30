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

## Prompts

1. Read this doc.  Perform the first task under Prep.

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