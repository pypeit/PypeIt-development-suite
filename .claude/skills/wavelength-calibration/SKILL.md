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

   - **Near-IR with no usable arc lamp:** use the OH sky-emission lines from the
     *science* frames instead. Type the science exposures `arc,tilt,science` in
     the `.pypeit` file and comment out any short dedicated arc — a long science
     exposure has far more OH signal, especially in the reddest order. Give the
     standard frames the same `calib` value so they inherit the solution. (See
     `keck_nires`, `soar_tspec`.)

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

   - `pypeit_compile_wvarxiv` expects per-order `wvarxiv` files written by
     `pypeit_identify`. If you instead already have a good automatic
     `WaveCalib_*.fits` (e.g. an echelle solution bootstrapped from a related
     instrument's arxiv), you can build the reid arxiv **directly from it** — it
     is just a `BinTableHDU` with one row per order and columns `wave[Nspec]`,
     `flux[Nspec]`, `order`. Read the `WaveCalib` (`wavecalib.WaveCalib.from_file`)
     and write each order's `wave_soln`, `spec`, and `ech_order`:
     ```python
     from pypeit import wavecalib
     from astropy.table import Table
     wc = wavecalib.WaveCalib.from_file('Calibrations/WaveCalib_A_0_DET01.fits')
     n = len(wc.wv_fits[0].wave_soln)
     t = Table(names=('wave','flux','order'),
               dtype=(f'({n},)>f8', f'({n},)>f8', '>i8'))
     for w in wc.wv_fits:
         t.add_row([w.wave_soln, w.spec, int(w.ech_order)])
     t.write('pypeit/data/arc_lines/reid_arxiv/<instr>.fits', format='fits')
     ```
     Mirror an existing arxiv's row order (e.g. `p200_triplespec.fits`). After
     this, the same arc self-reidentifies at cc ~ 1, so a previously-needed
     ``cc_thresh`` reduction can be reverted to the default — but validate on an
     *independent* dataset, since a self-built arxiv matches its own arc trivially.

   - **A single marginal order** (e.g. the reddest, low-S/N echelle order)
     failing reidentification at cc just under the default ``cc_thresh`` (0.7) is
     often recoverable by lowering ``cc_thresh`` (e.g. to 0.6); check the
     resulting per-order RMS is consistent with the others before trusting it.

5. **New or updated line list**: install via `pypeit_install_linelist`; keep
   lists in `pypeit/data/arc_lines/lists/`. For echelle, prefer the line list
   that covers the *full* wavelength range with real lines (e.g. `OH_NIRES`
   reaches ~24980 A) over one that stops short and forces extrapolation in the
   reddest order.

## Verify

- Re-run wavelength calibration non-interactively and confirm low RMS and a
  sensible dispersion solution (`pypeit_chk_wavecalib`,
  `pypeit_show_wvcalib`).
- Confirm sky/arc line positions land at expected wavelengths in
  `pypeit_show_2dspec`/`pypeit_show_1dspec`.
- If this supports a new instrument, register a dev-suite setup
  (`add-devsuite-setup`) and run it (`run-dev-suite`).
