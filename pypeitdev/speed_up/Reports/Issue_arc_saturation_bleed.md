# Arc line detection keeps the edges of saturated/bleed plateaus as "lines" (Keck/DEIMOS)

## Summary

In bright Keck/DEIMOS arc frames, `pypeit.core.arc.detect_lines` accepts
features that are the **edges or flat tops of extended saturated regions**,
not arc lines.  Their Gaussian centroids mark the plateau edge (offset
inward by a few pixels), not a line wavelength.  The existing
`nonlinear_counts` cut does not catch them, because the plateau sits well
below the non-linear level.  They are never used in a final wavelength fit,
but they do enter the line lists that `full_template` identification works
from, and whether they are present changes the identification of a few
slits.

## What the features are

Example: `keck_deimos/600ZD_M_6500`, `Arc_A_0_MSC04.fits`, slit 4.
- Around spectral pixels 6505–6574, the 2D arc reads **~39 900 counts right
  across the slit for ~60 spectral pixels**.
- Line cores inside that region reach **79 503 counts**, the maximum value in
  the frame, i.e. hard saturation.
- The region then drops to ~1 300 counts within one pixel.

This looks like charge bleeding from the hard-saturated cores along the CCD
columns, which on DEIMOS run along the spectral direction.  That reading has
not been verified further.

Within a 7- or 11-pixel fit window, the extracted arc spectrum then shows a
vertical step or a flat top, not a Gaussian:
- `curve_fit` places the steepest part of a Gaussian on the step, or matches
  the curvature of a slightly convex top.
- It returns a σ that depends on the window geometry; with `fwhm = 4`,
  σ ≈ 1.69 or 2.16 px on almost every one of them.
- The centroid lands 2–4 px inside the plateau edge.

## How many

Over the 370 arc spectra of the dev-suite wavelength calibrations (38
WaveCalib files):
- 46 such features pass every cut of `detect_lines`: 45 in
  `keck_deimos/600ZD_M_6500` (MSC01–04) and 1 in `keck_deimos/1200G_M_5500`;
- 45 of the 46 have raw peaks above 30 000 counts;
- 39 of those 45 sit on a vertical edge (within ±4 px of the peak pixel,
  adjacent pixels jump by 88–98% of the peak); the other 6 sit on a flat top;
- **none is used in a final wavelength fit**; the nearest final-fit line is
  ≥ 7 px away, median 26 px.

## Why it matters

- **Identification.**  In the speed-up work (PR "C") an alternative line fit
  that did not converge on these features dropped them.  Losing one or two of
  them flipped the `full_template` identification of several DEIMOS 600ZD
  slits.  For example, MSC01 slit 2 moved by ~30 px at the blue end, with
  4359.56 Å placed on the line develop calls 4377.35 Å, while the slit RMS
  changed only from 0.588 to 0.617 px.  The identification therefore depends
  on artefacts, in either direction.
- **FWHM measurement.**  Features with artefact σ can bias the per-slit FWHM
  that `autoid.measure_fwhm`/`map_fwhm` estimate from the detected lines.

## Why the existing cuts miss them

- `detect_lines` rejects lines with `tampl_true >= nonlinear_counts`.  The
  plateau is at ~39 900 counts, below DEIMOS `nonlinear_counts`
  (0.95 × saturation × gain).
- For `full_template`, `BuildWaveCalib` deliberately keeps saturated pixels in
  the arc image ("We want even the saturated lines in full_template for the
  cross-correlation"), so the 2D `msarc.image < nonlinear_counts` mask is not
  applied there.  It would not cover the bleed plateau anyway.

## Cuts tried on the 1D spectrum (none is clean)

These shape tests were applied to the fit windows of the
`curve_fit` + median-width-floor line list.  Each catches the 46 only by also
removing real lines:

| cut | of the 46 removed | final-fit lines removed |
|---|---:|---:|
| flat top: 3rd-highest / max pixel > 0.98 | 38 | 4 314 |
| vertical edge: max adjacent-pixel jump > 2.5× that of a Gaussian of the measured FWHM | 46 | 954 |
| poor fit: RMS residual / amplitude > 0.15 | 46 | 157 |
| raw > 30 000 and residual > 0.15 | 45 | 141 |
| raw > 30 000, flat top > 0.995 and jump > 0.8 | 33 | 35 |

Wide-slit DEIMOS arc lines are themselves flat-topped (FWHM 7–11 px), and
bright blends fit a single Gaussian poorly, so in 1D they resemble these
features.

## Proposal

Handle it in 2D, before the arc spectrum is extracted:
1. Flag pixels at or above `nonlinear_counts` in the arc frame (the
   saturated line cores).
2. Grow that flag along the bleed direction (the spectral direction on
   DEIMOS) over the contiguous plateau.  For example, extend it while the
   value stays within a few percent of the plateau level, or by the
   detector's known bleed length.
3. Carry the flag into the extracted arc spectrum (e.g. via `arccen_bpm`) so
   that `detect_lines` does not accept peaks inside or at the edge of a
   flagged region, while `full_template` can still use the flux for the
   cross-correlation if that is wanted.

Open questions:
- Is the plateau charge bleeding, or a property of the DEIMOS amplifiers?
- Should the flagged regions also be excluded from the `full_template`
  cross-correlation, or only from line detection?
- Which other instruments show this?  In the dev suite only DEIMOS does.

## Reproduction

Scripts (PypeIt development suite, `pypeitdev/speed_up/scripts/arcfit/`):
- `report_fig_saturated.py --list` lists the features, with a figure;
- `arcfit_revert_survey.py` and `arcfit_revert_summary.py` give the counts
  above;
- `arcfit_revert_jump.py` and `arcfit_revert_resid.py` test the 1D cuts.

They run on the stored WaveCalib files of a dev-suite reduction.
