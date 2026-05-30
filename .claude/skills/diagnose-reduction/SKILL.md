---
name: diagnose-reduction
description: Triage a failed or suspect run_pypeit reduction using the run log, QA HTML, and the pypeit_chk_* inspection tools (edges, flats, tilts, wavecalib, scattlight, noise). Use when a reduction errors, produces bad spectra, or the user asks why a calibration/extraction looks wrong.
---

# Diagnose a PypeIt reduction

Lives in `PypeIt-development-suite/`; operates on a reduction directory produced
by `run_pypeit` (often under this repo's `REDUX_OUT/`).

## Reference material

- Inspection scripts in `PypeIt/pypeit/scripts/` (installed as `pypeit_*`):
  `chk_edges`, `chk_flats`, `chk_tilts`, `chk_wavecalib`, `chk_scattlight`,
  `chk_alignments`, `chk_noise_1dspec`, `chk_noise_2dspec`, `chk_for_calibs`,
  `parse_slits`, `print_bpm`, `qa_html`.
- Display tools: `pypeit_show_2dspec`, `pypeit_show_1dspec`,
  `pypeit_show_wvcalib`, `pypeit_chk_flexure`.
- User docs: https://pypeit.readthedocs.io/ (QA pages and per-step guides).

## Procedure

1. **Read the log first.** Find the `*.log` in the redux dir (and the terminal
   output). Identify the failing step and the traceback. `[INFO]` lines mark
   each stage; the last successful stage tells you where to look.

2. **Check that frames were typed and calibrations grouped correctly.** Inspect
   the `.pypeit` file's data block and `pypeit_chk_for_calibs`. Many "failures"
   are mis-typed frames or a missing arc/flat in a calibration group.

3. **Inspect the relevant calibration**, working upstream from the failure:
   - Slit edges: `pypeit_chk_edges Calibrations/Edges_*.fits.gz`
   - Flats: `pypeit_chk_flats Calibrations/Flat_*.fits`
   - Wavelength: `pypeit_chk_wavecalib Calibrations/WaveCalib_*.fits` /
     `pypeit_show_wvcalib`
   - Tilts: `pypeit_chk_tilts Calibrations/Tilts_*.fits`
   - Scattered light: `pypeit_chk_scattlight ...`

4. **Inspect the 2D/1D products:**
   - `pypeit_show_2dspec Science/spec2d_*.fits` (look at sky-subtraction
     residuals, object traces, masking).
   - `pypeit_show_1dspec Science/spec1d_*.fits`.
   - Noise vetting: `pypeit_chk_noise_2dspec`, `pypeit_chk_noise_1dspec`.

5. **Open the QA HTML** (`QA/MF_*.html` and the science QA) for an overview of
   each step's diagnostic plots.

## Common causes to consider

- Wrong/insufficient calibration frames or mis-set `configuration_keys`.
- Wavelength solution failure → revisit with the `wavelength-calibration` skill.
- Parameter mistuning → adjust via the `.pypeit` file (see `add-parameter` for
  what each parameter does).
- A genuine code bug → reproduce minimally and check against `develop`.

## Verify the fix

- Re-run with `run_pypeit <file>.pypeit -o` and confirm the step now passes.
- If reproducible from dev-suite data, capture it as a setup (see
  `add-devsuite-setup`) and run it (see `run-dev-suite`).
