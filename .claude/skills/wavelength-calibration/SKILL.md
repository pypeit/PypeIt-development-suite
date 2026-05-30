---
name: wavelength-calibration
description: Build, validate, and archive a PypeIt wavelength solution — interactively identify arc lines with pypeit_identify, manage line lists, and construct/archive a reusable arc template. Use when a new setup needs a wavelength solution, an existing solution is poor, or the user wants to add a template to the PypeIt archive.
---

# Wavelength calibration and templates

Lives in `PypeIt-development-suite/`; operates on the sibling `PypeIt/` repo and
on dev-suite arc data.

## Reference material

- Interactive tool: `PypeIt/pypeit/scripts/identify.py` (`pypeit_identify`).
- Core algorithms & line lists: `PypeIt/pypeit/core/wavecal/`.
- Archive/template tools: `pypeit_compile_wvarxiv`, `pypeit_show_arxiv`,
  `pypeit_install_wvarxiv`, `pypeit_install_linelist`.
- Docs: the wavelength-calibration and "construct template" pages at
  https://pypeit.readthedocs.io/ (see `PypeIt/doc/calibrations/`).

## Procedure

1. **Get an arc spectrum.** Run `run_pypeit` far enough to produce the
   `Calibrations/Arc_*` / `WaveCalib_*` products, or use `pypeit_identify`
   directly on the arc frame for the setup.

2. **Identify lines interactively:**
   ```console
   pypeit_identify Arc_A_1_DET01.fits Slits_A_1_DET01.fits
   ```
   Mark known arc lines against the appropriate line list, fit, and check the
   RMS / residuals. Save the solution and (optionally) the new line IDs.

3. **Tune the automatic solution.** If a future automatic run should reproduce
   this, set the relevant `WavelengthSolutionPar` options (method, `lamps`,
   `reid_arxiv`, sigdetect, fwhm, etc.) in the `.pypeit` file — see the
   `add-parameter` skill for what each does.

4. **Build / archive a template** so other setups of the instrument can reuse
   it: follow the "construct template" doc, then `pypeit_compile_wvarxiv` and
   place/install it in the wavelength archive
   (`pypeit/data/arc_lines/reid_arxiv/`). Inspect with `pypeit_show_arxiv`.

5. **New or updated line list**: install via `pypeit_install_linelist`; keep
   lists in `pypeit/data/arc_lines/lists/`.

## Verify

- Re-run wavelength calibration non-interactively and confirm low RMS and a
  sensible dispersion solution (`pypeit_chk_wavecalib`,
  `pypeit_show_wvcalib`).
- Confirm sky/arc line positions land at expected wavelengths in
  `pypeit_show_2dspec`/`pypeit_show_1dspec`.
- If this supports a new instrument, register a dev-suite setup
  (`add-devsuite-setup`) and run it (`run-dev-suite`).
