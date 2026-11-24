# PR C (`speed_up_arcfit`) vs `develop`: full dev-suite comparison

*2026-11-18.  Model: Opus 5.5.  Companion to `Speed_Up_Arcfit_Report.md`
(the design and performance report of PR C).*

> **Note (2026-11-24).**  This comparison is of the *vectorized* version of
> PR C.  PR C has since been reverted to develop's per-line `curve_fit`,
> keeping the σ floor, the archive-sky FWHM cache and the §7 fixes (C.md Q20;
> `Speed_Up_Arcfit_Report.md` §0–§8).  The floor's effects described here
> carry over; the effects of the vectorized fit itself (e.g. dropped plateau
> edges, gained lines on which `curve_fit` gave up) do not.  The revert is
> validated in `Speed_Up_Arcfit_Report.md` §5.

## Summary

The dev suite ran end to end on `develop` (90c41e1a9) and on `develop` plus
PR C as it stands after Q16 ("Bxa").  PR C's code has these parts:
- the vectorized arc-line fit, B, is the default (`fast_line_fit = True`);
- the σ floor (`min_frac_med_width = 0.5`) is applied in `detect_lines`;
- the synthetic cross-correlation arc in `get_xcorr_arc` keeps the legacy
  fit (A).

| | develop | PR C |
|---|---|---|
| Reduction / afterburner tests | 332 / 334 | 332 / 334 (the same two GNIRS coadd2d/tellfit failures) |
| Vet tests | 109 passed / 6 failed | 108 passed / 7 failed |
| Dev-suite unit tests | 168 passed / 1 failed (ginga, headless) | 167 / 2 (+ `test_add_missing_slits`) |
| Slits/orders paired with a solution on both branches | 2 780 | 2 780 |
| Lines in the final fits (Nlin, sum) | 149 525 | 148 460 (−0.7%) |
| Median RMS (px) | 0.142 | 0.140 |
| Median IDs_Wave_cov | 94.1% | 94.0% |
| Median measured FWHM (px) | 3.990 | 4.004 |

**Wavelength solutions** (max |Δλ| / dispersion over the common line range):

| max shift | slits |
|---|---|
| < 0.01 px | 2 138 (77%) |
| < 0.1 px | 2 556 (92%) |
| < 0.5 px | 2 718 (98%) |
| ≥ 1 px | 37, in 24 setups |

Control runs on three setups show no run-to-run noise: re-running develop
moved 0 of 75 slits.  Every difference here is therefore caused by the code.

**The four regressions you asked me to trace (Q17) all have a root cause
outside PR C's line fit.**  PR C changes which lines are kept, and that
exposes a weakness that was already in develop.

| Problem | Root cause | Already in develop? | Proposed fix | Tested |
|---|---|---|---|---|
| `test_add_missing_slits` fails | A rejects a narrow Sobel spike only because `curve_fit` returned σ = −0.66.  B keeps it, so an extra intermediate edge appears and the test's hard-coded indices point at different edges.  **The final slits are unchanged.** | the rejection was luck | Make the test find edges by position, not by index | passes on both branches |
| GNIRS 32_SB_SXD order 8: blue-end wavelengths not monotonic, so `test_load_spec1d` fails | The order-8 FWHM, extrapolated to a row outside the order, goes 3.73 → 4.16 px because the σ floor removes 7 spikes.  Tilt tracing has a cliff at 4.03 px, where the bluest tilt line is lost and the 4th-order tilt fit extrapolates over 210 rows. | yes: develop sits 0.3 px from the cliff | `spec_order = 2` for the order-8 tilts (or a guard against extrapolation) | tilt error 62 → 1.3 px |
| MMIRS nep_as1_MOS slit 18: RMS 0.13 → 23 px | The template `mmt_mmirs_J_zJ.fits` holds three identical copies of the spectrum, so the global cross-correlation has two exactly tied peaks.  A difference at 1e-5 decides which copy is used, and one choice puts the red end of copy 1 at the slit's blue end.  (`HK_zJ` and `K3000_Kspec` are tripled too.) | yes: a tie | Store each template once | slit 18 fixed; develop unchanged on all 31 MMIRS slits |
| LDSS3 VPH-RED_mos slit 21: IDs_Wave_cov 102% → 1.7% | A 27-px box slit.  **Develop's solution is also wrong**: it does not increase monotonically and starts about 1 100 Å too red.  Its 102% coverage hides this. | yes: both fail | None in PR C.  A follow-up could add an RMS + monotonicity check to `full_template`. | — |

The rest of this report covers:
- §1: the suite-wide metrics;
- §2: the wavelength solutions;
- §3: the arc-line measurements;
- §4: the arc-line fits;
- §5: each trace, with a figure of the problem and of the proposed fix;
- §6: the remaining items and a recommendation.

## How the comparison was made

- **Code.**  Two git worktrees in `/mnt/tank/Astronomy/PypeIt/worktrees`:
  - `develop` at 90c41e1a9;
  - `arcfit`: the same base plus the squashed `speed_up_arcfit` changes.

  Each run selects its code with `PYTHONPATH`.
- **Runs.**  `./pypeit_test all` at `-t 1`, one suite at a time, each in a
  memory-capped systemd scope:
  - develop → `REDUX_OUT_PR/develop`;
  - PR C → `REDUX_OUT_PR/arcfitBxa`.
- **Controls.**  Three setups (GMOS GN_E2V_MULTI_R400_600,
  GS_HAM_MULTI_R400_700, MOSFIRE J_multi), re-run as:
  - `develop2`: develop again;
  - `arcfitA`: PR C with the legacy fit forced.
- **Metrics.**  `scripts/arcfit/arcfit_pr_metrics.py` reads every
  `WaveCalib_*.fits` and pairs slits by `spat_id`, or by order for echelle.
  - Nlin, RMS and measured_FWHM come from each `WaveFit`.
  - IDs_Wave_cov(%) is computed as `pypeit_chk_wavecalib` does.
  - The solution change is the maximum |λ_PR − λ_dev| over the range covered
    by both line lists, divided by the dispersion.
  - The stored `tcent` arrays give the detected-line comparison.
- **Outputs.**
  - Tables: `REDUX_OUT_PR/cmp_bxa/pr_metrics.{csv,npz}`.
  - Figures: `figures/pr_compare_bxa/` (`arcfit_pr_figs.py`) and
    `figures/q17/` (`arcfit_q17_figs.py`).
  - Trace harnesses: `scripts/arcfit/q17/`.
  - Harness outputs: `REDUX_OUT_PR/q17/`.

## 1. Metrics

![RMS](figures/pr_compare_bxa/rms.png)

*Figure 1.  RMS of the final fit per slit, develop (x) vs PR C (y), coloured
by how far the solution moved.*  Where the solution did not move, the RMS
changes only by the small amount that comes from the centroids.  Where the
solution moved, the RMS spreads both ways.  Of the 37 slits that moved by
≥ 1 px, the RMS improves for 20 and worsens for 17.

![Nlin and coverage](figures/pr_compare_bxa/nlin_cov.png)

*Figure 2.  Left: change in Nlin per slit.  Right: change in IDs_Wave_cov,
clipped at ±20%.*  Nlin is unchanged on 2 136 of the 2 780 slits.  Where it
changes, it usually falls by 1–5 lines.  These are lines that the σ floor
removes, or plateau edges that B no longer passes to `curve_fit` (Q14).  Over
the suite, Nlin falls by 0.7%.  IDs_Wave_cov is unchanged (< 0.05%) on 2 645
slits.  71 slits lose ≥ 1% and 17 gain ≥ 1%.  Eight slits lose ≥ 20%:
- the flips of §5.3–5.4;
- UVES 860 B_1 order 11 and UVES 760 G_6 order 0;
- DEIMOS 1200G_M_5500 B_0 DET02 slit 12;
- three FIRE echelle orders (1, 2, 6).

The FIRE orders and the DEIMOS slit lose lines at the ends of their range
while the solution moves by < 0.75 px.  For them, IDs_Wave_cov measures the
span of the lines that were kept, not the quality of the solution.

![FWHM](figures/pr_compare_bxa/fwhm.png)

*Figure 3.  measured_FWHM per slit.*  The median shifts by +0.004%, but
213 slits change by > 1%, and almost all of them become wider.  This is the
σ floor at work.  The narrowest detections (σ < 0.5 × the median σ of the
spectrum) used to pull `measure_fwhm`/`map_fwhm` down, and now they are
gone.  This matters in §5.2: the FWHM is passed on to tilt tracing and to
line detection for identification, so a wider FWHM changes which lines are
found downstream.

## 2. Wavelength solutions

![Δλ](figures/pr_compare_bxa/dlambda.png)

*Figure 4.  Left: cumulative distribution of the solution change.  Right: per
instrument, the slits that are unchanged (< 0.01 px), moved 0.01–0.5 px, or
moved ≥ 0.5 px.*

- **77% of slits** reproduce develop to better than 0.01 px.
- **92%** reproduce it to better than 0.1 px.
- Most of the 0.01–0.5 px band is centroid noise propagating into the fit
  (§3).
- **The changes cluster by method.**
  - `holy-grail` and `reidentify` setups with many lines barely move.
  - `full_template` multislit setups and echelle orders with few lines can
    jump.  There a different line list, or a different FWHM, changes a
    cross-correlation or a pattern match, and the solution flips (§5).

The 37 slits that moved by ≥ 1 px:

| setup | key | slit | Nlin dev → PR | RMS dev → PR (px) | cov dev → PR (%) | max shift (px) | note |
|---|---|---|---|---|---|---|---|
| vlt_uves_red/860 | B_1_MSC01 | 11 | 48 → 29 | 0.54 → 0.46 | 85 → 27 | 1.5e5 | order lost at one end; also failed under B |
| magellan_ldss3/VPH-RED_mos | A_1_DET01 | 21 | 14 → 10 | 0.71 → 0.29 | 102 → 1.7 | 2 615 | §5.4: wrong on both |
| mmt_mmirs/nep_as1_MOS | A_0_DET01 | 18 | 101 → 96 | 0.13 → 23.1 | 94 → 46 | 2 290 | §5.3: template tie |
| gemini_gmos/GN_HAM_B480 | A_0_MSC01 | 0 | 16 → 10 | 1.07 → 0.25 | 91 → 82 | 1 130 | gain |
| keck_mosfire/J_multi | A_0_DET01 | 2 | 16 → 12 | 0.57 → 0.95 | 71 → 61 | 189 | OH, flat-topped lines |
| gtc_osiris/R1000RMOS | A_0_DET02 | 2 | 16 → 28 | 1.78 → 3.71 | 67 → 70 | 139 | poor on both |
| magellan_ldss3/VPH-RED_mos | A_1_DET01 | 0 | 22 → 22 | 0.74 → 0.86 | 76 → 73 | 36 | box slit |
| keck_deimos/600ZD_M_6500 | A_0_MSC02 | 20 | 44 → 43 | 21.2 → 17.2 | −58 → −59 | 12 | junk on both |
| apf_levy/bin2 | A_0_DET01 | 8 | 57 → 47 | 0.05 → 0.10 | 91 → 91 | 11 | |
| shane_hamspec/Hamilton_e2v | A_0_DET01 | 78 | 17 → 10 | 0.70 → 0.24 | 95 → 85 | 11 | |
| apf_levy/bin2 | A_0_DET01 | 16 | 45 → 39 | 0.29 → 0.13 | 96 → 96 | 9.4 | |
| apf_levy/bin2 | A_0_DET01 | 41 | 35 → 33 | 0.05 → 0.21 | 89 → 94 | 8.7 | |
| keck_lris_red_orig/multi_600_5000_d500 | A_0_DET01 | 15 | 35 → 34 | 0.04 → 0.05 | 97 → 97 | 8.5 | |
| keck_mosfire/J_multi | A_0_DET01 | 0 | 20 → 24 | 0.78 → 0.75 | 80 → 60 | 6.2 | |
| keck_hires/HS1700+6416 (RED B2) | A_0_MSC01 | 53 | 9 → 6 | 0.30 → 0.035 | 67 → 65 | 3.8 | |
| keck_kcwi/large_bl | D_0_DET01 | 16 | 78 → 84 | 0.22 → 0.36 | 88 → 88 | 3.7 | |
| keck_kcwi/large_bl | D_0_DET01 | 10 | 84 → 80 | 0.33 → 0.23 | 86 → 87 | 3.4 | |
| keck_hires/J1723+2243 (RED C5) | A_0_MSC01 | 30 | 31 → 27 | 0.49 → 0.27 | 94 → 95 | 3.4 | |
| magellan_ldss3/VPH-RED_mos | A_1_DET01 | 9 | 61 → 62 | 0.20 → 0.23 | 100 → 96 | 3.2 | |
| keck_lris_red/multi_1200_9000_d680_1x2 | A_0_DET01 | 0 | 29 → 32 | 0.08 → 0.22 | 91 → 100 | 2.7 | |
| mmt_mmirs/nep_as1_MOS | A_0_DET01 | 15 | 56 → 59 | 1.69 → 1.59 | 78 → 79 | 2.5 | |
| keck_lris_red/multi_900_5500_d560 | A_0_DET01 | 2 | 33 → 33 | 0.09 → 0.03 | 90 → 87 | 2.5 | |
| magellan_ldss3/VPH-RED_mos | A_1_DET01 | 5 | 31 → 31 | 0.89 → 0.92 | 64 → 54 | 2.5 | box slit |
| keck_hires/J0100+2802 (RED C1) | A_0_MSC01 | 20 | 56 → 56 | 0.24 → 0.26 | 95 → 93 | 2.4 | |
| apf_levy/bin2 | A_0_DET01 | 54 | 27 → 30 | 0.04 → 0.20 | 96 → 95 | 2.3 | |
| gemini_gmos/GS_HAM_MULTI_R400_700 | A_0_MSC01 | 3 | 50 → 49 | 1.26 → 1.14 | 98 → 99 | 2.3 | |
| vlt_sinfoni/K_0.8 | A_1_DET01 | 30 | 20 → 19 | 0.15 → 0.18 | 93 → 93 | 2.1 | |
| mmt_mmirs/nep_as1_MOS | A_0_DET01 | 2 | 75 → 48 | 1.39 → 1.37 | 84 → 85 | 2.0 | |
| mmt_binospec_ifu/G270 | A_0_DET01 | 15 | 24 → 26 | 0.10 → 0.21 | 97 → 97 | 1.8 | |
| mmt_binospec_ifu/G270_Sky | A_0_DET01 | 15 | 24 → 26 | 0.10 → 0.21 | 97 → 97 | 1.8 | |
| gemini_gmos/GS_HAM_B480_550 | A_0_MSC01 | 2 | 80 → 78 | 0.30 → 0.26 | 100 → 90 | 1.6 | |
| vlt_sinfoni/K_0.8 | A_1_DET01 | 16 | 21 → 19 | 0.17 → 0.08 | 93 → 93 | 1.4 | |
| keck_lris_red_mark4/multi_1200_9000_d680 | A_all_DET01 | 11 | 26 → 27 | 0.07 → 0.08 | 92 → 90 | 1.4 | |
| vlt_uves_red/760 | G_6_MSC01 | 0 | 17 → 15 | 0.42 → 0.38 | 59 → −18 | 1.4 | |
| vlt_uves_red/860 | B_1_MSC01 | 18 | 25 → 15 | 0.80 → 0.12 | 92 → 91 | 1.3 | |
| vlt_xshooter/VIS_1x1 | A_0_DET01 | 11 | 42 → 41 | 0.75 → 0.70 | 100 → 99 | 1.2 | |
| gemini_gmos/GS_HAM_MULTI_R400_700 | A_0_MSC01 | 40 | 71 → 58 | 1.11 → 1.07 | 99 → 99 | 1.0 | |

![RMS vs Δλ](figures/pr_compare_bxa/rms_vs_dlambda.png)

*Figure 5.  RMS ratio (PR / develop) against the solution change, with the
outliers labelled.*  The three big jumps at the top right and right of the
plot are traced in §5.  Most slits that moved by 1–10 px are echelle orders
with 6–60 lines, where one line more or fewer moves the edge of the
polynomial.

Two remarks on what IDs_Wave_cov and RMS measure:
- **A good metric value can hide a wrong solution.**  LDSS3 slit 21 has
  102% coverage on develop and is wrong.  MMIRS slit 18 under PR C keeps 46%
  coverage while its RMS is 23 px.
- **Multislit `full_template` has no safeguard against either case**
  (§5.5).

## 3. Arc-line measurements

![tcent](figures/pr_compare_bxa/tcent.png)

*Figure 6.  Left: detected lines lost and gained per slit.  Right:
per-slit median and 99th-percentile centroid change of the lines that both
runs detect.*

- **Detections.**
  - 223 636 → 218 838 detected lines.
  - 217 788 are matched between runs: 5 848 lost and 1 050 gained.
  - The losses are mostly the σ floor (sub-resolution spikes) and plateau
    edges that B rejects instead of passing to `curve_fit`.
  - The gains are mainly lines on which `curve_fit` gives up (a
    `RuntimeError`, the −999 return, or a negative σ as in §5.1) but the
    analytic fit converges.  This was measured on the stored arcs in Task 10,
    not re-measured line by line here.
- **Centroids.**  For the matched lines the per-slit median centroid change
  is 10⁻⁷–10⁻⁶ px, i.e. B and `curve_fit` converge to the same minimum.  The
  99th percentile reaches 10⁻³–10⁻¹ px on blended or asymmetric lines, where
  the two fits stop at different points of a flat χ² surface.
- **These are the inputs of every downstream step**: identification, FWHM
  map, tilts and flexure.  Where the downstream step is continuous in its
  inputs (a least-squares fit to many lines), the solutions agree to
  < 0.01 px.  Where it is discontinuous, small input changes can flip the
  result (§5):
  - an argmax over a cross-correlation;
  - a threshold on the FWHM;
  - a hard-coded index.

## 4. Arc-line fits

The final 1D wavelength fits differ in only three ways:

1. **Line lists that differ by 1–5 lines.**  The σ floor removes narrow
   features, and B handles failed fits differently.  On well-sampled setups
   this changes the RMS by a few per cent and the solution by < 0.1 px.
2. **FWHM-driven changes.**  `rms_thresh_frac_fwhm`, the match tolerance in
   `reidentify` and `detect_lines` all scale with the measured FWHM.  A wider
   FWHM can:
   - admit two detections 1.5–2.7 px apart that are both matched to the same
     wavelength, as in GNIRS order 8 (8772.5 Å and 9317.7 Å each appear
     twice in PR C's fit);
   - lose a line.
3. **Discontinuous flips** in `full_template` and `reidentify`, where the
   correlation peak or the anchor lines change (§5.3–5.4).

The vet tolerances change in two places:
- **`test_gmos`:** GS_HAM_R150_869's RMS rises from 0.234 to 0.272 px against
  a 0.25 px tolerance, with a solution shift of 0.1 px.  This is case 1, and
  it is identical under B.
- **`test_keck_hires`:** fails on both branches, on different setups
  (develop J1723+2243, PR C J0306+1853).  This is case 1 on echelle orders
  already close to the tolerance.

## 5. The traced regressions

### 5.1 `test_maskdesign_matching::test_add_missing_slits` (DEIMOS 1200G_Cooper)

![maskdesign](figures/q17/q17_maskdesign.png)

*Figure 7.  Left: the collapsed, rectified Sobel image around spatial pixel
299, the B fit, and the identical Gaussian that `curve_fit` returns with a
negative σ.  Right: edge positions at the PCA reference row; the PR has one
extra edge at index 3.*

- **What happens.**
  1. `EdgeTraceSet.peak_refine` calls `trace.peak_trace`, which detects
     edge peaks in the Sobel image with `detect_lines`.
  2. A narrow positive spike (σ = 0.66 px, FWHM_nominal 3 px) sits at the
     foot of a deep negative trough at 303 px.
  3. `curve_fit` converges on it with **σ = −0.66 px**: the Gaussian depends
     only on σ², so the sign is arbitrary.  `detect_lines` then drops it on
     `twid > 0`.
  4. B returns |σ| and keeps it.
- **Result.**  B finds 62 traces instead of 61.  Every index after 3 shifts
  by one, so `edge_fit[..., 7]`, `[9]` and `[10]` are different edges.
- **The σ floor does not help.**  The spike's σ is 0.51 × the median σ (the
  floor is 0.5), and `peak_trace` disables the floor anyway
  (`min_frac_med_width=0.`).
- **Effect on the reduction: none.**  With `auto=True`, the final DEIMOS slits
  are identical on both branches (31 slits, the same `maskdef_id`, edges
  within 0.47 px).  The extra edge is removed by the later
  sync/maskdesign steps.  MOSFIRE J_multi gives the same 21 traces on both
  branches.
- **Proposed fix (dev suite, not PypeIt).**  Find the four edges in the test
  by position (the nearest trace to 458/1736/496/561, and 349/2032/973/1100
  for MOSFIRE), not by index, and keep the same position assertions before
  and after `maskdesign_matching()`.  The rewritten test
  (`scripts/arcfit/q17/tmm_pos.py`) **passes on develop and on PR C**.

### 5.2 GNIRS 32_SB_SXD order 8 (`test_specutils::test_load_spec1d`)

![GNIRS FWHM](figures/q17/q17_gnirs_fwhm.png)

*Figure 8.  Order 8's FWHM map.  Grey: lines kept by both branches.  Blue
triangles: seven sub-resolution spikes kept by A only.  Lines: the
`fwhm_spec_order = 1` fits.  The order is only on the detector for rows
513–1021, yet its FWHM is evaluated at row 511.*

![GNIRS tilts](figures/q17/q17_gnirs_tilts.png)

*Figure 9.
(a) Maximum |tilt − row| along the order centre against the FWHM passed to
the tilt tracing, for `spec_order` 4 (the GNIRS default) and 2.
(b) The order-8 tilt residual for develop, PR C, and PR C with
`spec_order = 2`.
(c) The resulting spec1d wavelengths at the blue end of order 8.*

The chain (each step reproduced in isolation, `scripts/arcfit/q17/gn6.py`,
`gn7.py`):

1. **The wavelength solution is not the problem.**  The 2D echelle fit is
   monotonic and almost identical on both branches.
2. **The tilts are the problem.**  On PR C, the order-8 tilt at row 520 is
   577 instead of 520, so the wavelength image is evaluated about 57 rows too
   red there.  That is why the spec1d starts at 8 433 Å instead of 8 238 Å
   and wobbles (14 non-increasing pixels).
3. **Which fit traces the tilts (A or B) makes no difference.**  With the
   same FWHM, A and B give identical tilts.
4. **The tilts change because `calibrations.py` passes the order's measured
   FWHM from WaveCalib to `BuildWaveTilts`.**
   - That FWHM is `map_fwhm(...).eval(nspec//2, 0.5)`, i.e. row 511, which
     is off the order: an extrapolation of a linear fit.
   - On develop, A keeps seven spikes with FWHM 0.6–1.6 px at rows 597–831.
     They pull the fit down to 3.73 px at row 511.
   - On PR C the σ floor removes the spikes, and the fit gives 4.16 px.
   - Over rows 600–1000 the fitted FWHM is 3.86–4.43 px on develop and
     4.22–4.46 px on PR C (Fig. 8).  So the PR's value is the more honest
     one, and develop's FWHM is above the cliff over most of the order.
5. **The tilt tracing has a cliff.**
   - At FWHM ≤ 4.02 px the bluest tilt line (row 552) is traced and used.
   - At ≥ 4.05 px a neighbouring feature at row 559 is detected instead and
     rejected in `trace_tilts`.  The next usable line is at row 723.
   - The 4th-order spectral tilt fit then extrapolates over rows 513–723 and
     curls by up to 62 px.
   - Develop's 3.73 px is only 0.3 px below the cliff.
6. **Proposed fix (tested).**  Set `spec_order = 2` for the order-8 tilts:
   - the cliff drops from 52–63 px to ≤ 1.3 px for every FWHM from 3.5 to
     4.5 px;
   - at 3.73 px it is 0.14 px, against develop's 0.43 px.

   In config form, for GNIRS 32_SB_SXD:

   ```ini
   [calibrations]
     [[tilts]]
       spec_order = 2, 4, 4, 4, 4, 4
   ```

   A general alternative, for a follow-up: cap the spectral order of the tilt
   fit when the usable lines cover only part of the order.
7. **Separate issue, for a follow-up.**  `map_fwhm` is evaluated at
   `nspec//2` even for orders that do not cover that row
   (`wavecalib.py:649`).  Evaluating at the order's own spectral midpoint
   would be more meaningful.  It does **not** fix this case on its own:
   the order-8 midpoint (row 767) gives 4.10 px on develop and 4.32 px on
   PR C, both above the cliff, so develop would fail too.

### 5.3 MMIRS nep_as1_MOS slit 18 (RMS 0.13 → 23 px)

![MMIRS](figures/q17/q17_mmirs.png)

*Figure 10.
(a) The `full_template` template `mmt_mmirs_J_zJ.fits`: its wavelength vector
runs 9954 → 15066 Å three times, and the three 2040-px segments are
bit-identical.
(b) Slit 18's global cross-correlation with the template on both branches.
(c) The solutions: develop; PR C; PR C with a one-copy template; develop with
a one-copy template.*

- **Where it diverges.**  A harness that re-runs `full_template` on the
  stored arc spectra (`scripts/arcfit/q17/ft.py`) reproduces the pipeline
  exactly: develop RMS 0.129, PR C 23.11, B 0.114.  Instrumenting it shows
  the divergence in the very first, global `xcorr_shift`.  Develop picks a
  lag of +1793.6 px and PR C picks −246.4 px.
- **Why the two lags compete.**
  - Because the template repeats, the correlation function has peaks of
    *identical* height (0.253480) at lags 5872.6 and 7912.6.  The two
    branches' correlation functions differ by at most 3 × 10⁻⁷.
  - `xcorr_shift` chooses between them on the continuum-subtracted peak
    height from `detect_lines(nfind=1)`:
    - develop: 0.137161 vs 0.137154;
    - PR C: 0.136362 vs 0.136351.
  - The tie is broken by noise at 10⁻⁵, and PR C with A forced everywhere
    falls on the bad side too.
- **Why one copy is bad.**  Slit 18 reaches bluer (9 335 Å) than the
  template (9 954 Å).  With one choice the template window is zero-padded at
  the blue end.  With the other, the window reaches into the previous copy, so
  the slit's first ~250 px are matched to 14 500–15 066 Å.  The first
  `reidentify` snippet then anchors 28 wrong lines, and the 4th-order fit
  folds over.
- **The other MMIRS templates.**  `mmt_mmirs_HK_zJ.fits` (3 × 998 px) and
  `mmt_mmirs_K3000_Kspec.fits` (3 × 2040 px) are tripled the same way.  A
  scan of all 254 `reid_arxiv` templates finds no other tiled copies; the
  other non-monotonic files are multi-slit concatenations or single repeated
  pixels.
- **Proposed fix (data only).**  Rewrite the three MMIRS templates with one
  copy each.  Truncating at the first seam in the harness:
  - **fixes slit 18 on PR C** (RMS 0.114, 110 lines, the same solution as
    develop);
  - raises the global correlation from 0.25 to 0.40;
  - leaves **develop unchanged on all 31 MMIRS slits** in nep_as1_MOS,
    J_zJ, HK_zJ and K_K (< 0.01 px).
  - On PR C, the only nep_as1_MOS slits still > 1 px from develop are 2, 5,
    15 and 24, all ordinary line-measurement differences (≤ 11 px, RMS
    within ±0.15 px).

  The proposed files would be written by
  `wave[:2040], flux[:2040]` (`[:998]` for HK_zJ).

### 5.4 LDSS3 VPH-RED_mos A_1 slit 21 (IDs_Wave_cov 102% → 1.7%)

![LDSS3](figures/q17/q17_ldss3.png)

*Figure 11.
(a) Wavelength at pixel 0 for every slit on the mask.  Stars mark the four
27-px box slits (0, 1, 5, 21); the others are 38 px.
(b) Slit 21's solution on both branches against its neighbours.*

- **Not a regression.**  Slit 21 is a box slit, with broader lines
  (FWHM 8.1 px against 6.4–6.7 px for the science slits).
  - **Develop's solution is already wrong.**  It turns over at the red end
    and starts about 1 100 Å redder than slits 20 and 22.
  - It fits only 14 lines, against 59–71 on the science slits.
  - The 102% coverage comes from a polynomial that spans more than the
    detector.
  - PR C's solution is wrong in a different way.
- **Where it fails.**  In the harness, on both branches the second
  `reidentify` snippet identifies 0 lines.  The first identifies 3–4 with a
  shift of +328 px (develop) or −208 px (PR C): both are spurious matches
  with cc ≈ 0.22.
- **The other box slits.**  Slits 0, 1 and 5 are also poor on develop (RMS
  0.57–0.89 px, 16–35 lines).
- **No PR C change proposed.**  A box slit is not a spectrum the
  `full_template` method can calibrate reliably.  Either flag such slits, or
  add the safeguard of §5.5.

### 5.5 The missing safeguard in multislit `full_template`

Neither MMIRS slit 18 (RMS 23 px) nor LDSS3 slit 21 (not monotonic) was
flagged as BADWVCALIB.
- The `rms_thresh_frac_fwhm` check runs only in the echelle redo path
  (`wavecalib.py:793`) and before the echelle 2D fit.  Multislit
  `full_template` solutions are accepted as they are.
- There is no monotonicity check anywhere.

Checking all 179 `full_template` multislit setups in the suite:

| | develop | PR C |
|---|---|---|
| `full_template` slits | 1 391 | 1 394 |
| RMS > `rms_thresh_frac_fwhm` × FWHM | 49 | 41 |
| non-monotonic solution | 6 | 8 |
| either | 54 (3.9%) | 47 (3.4%) |

The non-monotonic slits are:
- **on both branches:** DEIMOS 1200G_M_5500 B_0 DET02 slit 3 and 600ZD MSC02
  slit 20; LRIS-blue slitless; LDSS3 slit 21; NOT ALFOSC grism10 and
  grism11;
- **PR C only:** MMIRS slit 18, and GMOS GS_HAM_B480_550 MSC01 slit 0
  (develop: no solution; PR C: a 5-line fit with −0.06% coverage).

Adding the checks would flag 54 develop slits.  That is a behaviour change
in its own right and needs its own review, so I recommend a separate PR.

## 6. Remaining items and recommendation

**Also changed, not traced in depth:**
- **UVES red 860 B_1 order 11:** cov 85% → 27%, 48 → 29 lines.  It already
  failed the same way under B.  This is an echelle order with lines lost at
  one end; the mechanism is probably that of §5.2 (FWHM or line-list
  change → reidentify anchor change).
- **`test_gmos`:** R150_869 RMS 0.272 against a 0.25 px tolerance, with a
  0.1 px shift.  This is a real but small line-measurement change.  Either
  relax the tolerance to 0.30 or accept the failure as documented.

**Recommended course:**
1. **In this PR:**
   - replace the three tripled MMIRS templates (§5.3);
   - set `spec_order = 2` for GNIRS 32_SB_SXD order 8 tilts (§5.2);
   - update `test_add_missing_slits` to find edges by position (§5.1);
   - re-run the affected setups (MMIRS × 6, GNIRS 32_SB_SXD) and the vet
     tests.
2. **Follow-ups:**
   - the multislit `full_template` RMS and monotonicity safeguard (§5.5);
   - evaluating `map_fwhm` at each order's own spectral midpoint (§5.2, item 7);
   - a guard against extrapolation in the tilt fit.

## 7. Fixes applied (Q18 = (a), 2026-11-18)

Three fixes went into PR C, plus the `test_gmos` tolerance:

| Fix | Where | Result |
|---|---|---|
| One copy per MMIRS template (`J_zJ`, `HK_zJ`, `K3000_Kspec`) | PypeIt `pypeit/data/arc_lines/reid_arxiv/` | nep_as1_MOS slit 18: RMS 23.1 → **0.114 px**, 110 lines (develop 0.129 px, 111 lines); within 0.08 Å of develop.  The other 27 slits and the 5 longslit setups are unchanged from the unfixed PR C run. |
| GNIRS 32_SB_SXD tilts `spec_order = [2, 4, 4, 4, 4, 4]` | PypeIt `gemini_gnirs.py` | Order 8 is monotonic in both standard-star spec1d files (unfixed PR C: 14 and 1 non-monotonic pixels), spanning 8 233–9 466 Å (develop 8 238–9 463 Å).  `test_load_spec1d` passes. |
| `test_add_missing_slits` finds edges by position | dev suite `unit_tests/test_maskdesign_matching.py` | Passes (3/3 in the module). |
| `test_gmos` R150_869 tolerance 0.25 → 0.30 px | dev suite `vet_tests/test_wavelengths.py` | Passes (RMS 0.272 px). |

**Re-run.**  The 6 MMIRS setups and GNIRS 32_SB_SXD were re-reduced with the
fixes (`scripts/arcfit/arcfit_q18_fix.sh`, output
`REDUX_OUT_PR/arcfitBxaFix`).  The vet tests were run on
`REDUX_OUT_PR/bxafix_merged`, which links every other setup to the unfixed
PR C run.
- Reductions: 10/12 pass.  The two failures (GNIRS `pypeit_coadd_2dspec`
  failed, `pypeit_tellfit` skipped) are the same on develop.
- Vet tests: 107 passed, 7 failed.  Five failures are the same as on
  develop (`test_specutils` × 4, `test_keck_hires`).  The other two are
  artefacts of re-running the vet tests on an already-vetted tree, not of
  the fixes:
  - `test_spat_flexure_science`: an earlier vet test (Kast red `redoslit`)
    re-reduces Kast red 600_5000_d46 into the same directory without
    spatial flexure.  It ran during the unfixed run's vet pass, after
    `test_spat_flexure_science` had passed, so the linked spec2d files no
    longer carry a flexure value.  None of the five setups this test reads
    was re-run.
  - `test_show_2dspec`: ginga failed to start (display timeout), like
    `test_view_fits_proc` in the earlier runs.
- PR C's three regressions in the vet and unit tests
  (`test_add_missing_slits`, `test_load_spec1d`, `test_gmos`) now pass.
