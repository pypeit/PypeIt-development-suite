# Report 04 — Wavelength-calibration deep dive (Shane/Hamspec)

**Date:** 2026-06-17
**Author:** JXP and Claude
**Task:** Development / Wavelength Calibration #1 — examine the provided
wavelength reference files, characterize the ThAr solution, and lay out the
path to a PypeIt echelle archive.
**Scripts:** `pypeitdev/shane_hamspec/scripts/` (run with the `pypeit14` env)

This is the blocker identified in [Report 03](03_first_reduction.md): the
`echelle` wavelength step needs the archive files returned by
`get_echelle_angle_files()` (`shane_hamspec_angle_fits.fits`,
`shane_hamspec_composite_arc.fits`), which do not exist yet. We have **two
independent legacy reductions** of the same ThAr exposure to build them from.

---

## 1. Inputs examined

In `pypeitdev/shane_hamspec/`:

| File | What it is |
|------|-----------|
| `IRAF/ecTHAR20240125.rsfbo_1` | IRAF `ecidentify` **database** (text): 1355 identified ThAr lines + the global 2D Chebyshev solution. |
| `IRAF/THAR20240125.rsfbo_1.fits` | Extracted ThAr arc, **103 orders × 4096 pix** (equispec / pixel WCS). |
| `IRAF/THAR20240125.wrsfbo_1.fits` | Same, **MULTISPE** WCS — per-order wavelength solutions in 285 `WAT2` cards. |
| `Arc_01_fit.idl` | XIDL (HIRES/ESI-style) save: `guess_ordr`, `sv_aspec` (extracted arc), `sv_lines` (IDed lines), `all_arcfit` (per-order 1D fits). |
| `Hamilton_ThArTi.pdf` | Screenshots of the IRAF solution (human reference). |

**Key discovery:** PypeIt already ships a reader for the exact XIDL schema —
`pypeit.core.wavecal.templates.xidl_esihires()` — which evaluates the
per-order fits and returns `(orders, wave_vacuum, spec)` directly. This is the
single most useful fact for building the archive.

## 2. Scripts written (`scripts/`)

- **`parse_iraf_ecthar.py`** — parses the `ecidentify` database into a
  structured array of features `(aperture, order, pixel, λ_fit, λ_ref,
  weight)` plus the 2D-fit header/coefficients; `summary()` prints stats.
- **`read_xidl_arc.py`** — wraps `templates.xidl_esihires()` to load the XIDL
  `(orders, wave, spec)` and pulls the per-order RMS from `all_arcfit`.
- **`characterize_wavesoln.py`** — combines both, prints a comparison, and
  writes the figures below.

## 3. What the solution looks like

### Geometry / physics
- Cross-dispersed echelle; dispersion along columns (`specaxis=1`, confirmed
  in Report 01). **4096 pixels per order.**
- **Grating constant** `m·λ ≈ 5.71 × 10⁵ Å` (IRAF: median 570 663, fit
  constant 570 757; XIDL center: 571 145 — agree to ~0.1%). This is the single
  number that ties physical order number to wavelength.

### Coverage (Figure 1)
- **IRAF**: 102 orders, **physical order m = 58 → 159**, **λ ≈ 3577 – 9788 Å
  (air)**, 1355 lines (median 11 lines/order, range 1–37).
- **XIDL**: 64 orders, **m = 65 → 146**, **λ ≈ 3858 – 8899 Å (vacuum)**,
  median 10 lines/order (4–21).
- The two use the **same physical order numbering** (Figure 1: XIDL points sit
  inside the IRAF line spans). IRAF simply solved **more** orders — it adds the
  reddest (58–64) and bluest (147–159) that XIDL left out. So IRAF is the more
  complete catalog; XIDL is the more convenient (ready-made extracted spectra +
  vacuum wavelengths).

### Accuracy
- **IRAF** global 2D Chebyshev (xorder=5 in pixel, yorder=4 in order):
  residual **RMS 13 mÅ**, max 55 mÅ — i.e. ≈ 0.01 px, excellent.
- **XIDL** per-order 1D fits (CHEBY/POLY, order 2–5): **RMS ≈ 1 mÅ** median
  (≈ 35 mÅ worst order). Both are sub-pixel and mutually consistent.

### Dispersion (Figure 2)
- Central dispersion runs ≈ **0.026 Å/pix in the blue (m≈146)** up to
  ≈ **0.057 Å/pix in the red (m≈65)** — it grows toward the red as λ/m rises
  (a single order spans more Å). Line density peaks in the mid/red orders
  (per-order ThAr counts up to ~37; see Figure 2).

### Example order (Figure 3)
- A mid-format XIDL order (m=106, ~5310–5460 Å) showing the extracted ThAr
  spectrum on its vacuum-wavelength grid — clean, well-resolved lines; the
  per-order product a PypeIt archive needs.

## 4. What PypeIt needs, and how these map to it

`shane_hamspec` uses `wavelengths['method'] = 'echelle'` with
`get_echelle_angle_files()` → two archive files:

1. **`shane_hamspec_composite_arc.fits`** — a **reference ThAr spectrum per
   order with its wavelength solution**, used for reidentification by
   cross-correlation. **Directly buildable** from either reduction:
   - *XIDL path (easiest):* `templates.xidl_esihires('Arc_01_fit.idl')` already
     returns `(orders 65–146, wave_vac, spec)` — feed straight into the
     PypeIt archive builder. Covers 64 orders.
   - *IRAF path (most complete):* combine `THAR…rsfbo_1.fits` (raw extracted
     spectra) with the `ecTHAR` 2D solution (or the `wrsfbo` MULTISPE WCS) to
     get all 102 orders. More work to parse the MULTISPE/2D fit.
2. **`shane_hamspec_angle_fits.fits`** — maps an **order's spatial (cross-
   dispersion) position to the cross-disperser angle** so PypeIt can predict
   which physical orders land where for a given `xdangle` (GTILTRAW). **This is
   the gap:** we have data at **only one** cross-disperser setting
   (GTILTRAW = 9194). We can anchor that single setting (associate the traced
   orders with m = 58…159 via `m·λ`), but we cannot *fit* an order-vs-angle
   relation from one point. See Q&A.

The reidentification itself (matching the 171 traced orders from Report 03 to
physical m and assigning λ) is what failed at `build_wv_calib: Finding the
echelle orders`; it needs both files above.

## 5. Recommended path (for the next task)
1. **Build the composite arc** from the **XIDL** solution first (one call to
   `xidl_esihires`), since it returns vacuum λ + spectra ready to go; optionally
   extend to all 102 orders later using the IRAF catalog. Air→vacuum is already
   handled by `xidl_esihires` (XIDL/IRAF wavelengths are **air**).
2. **Decide the order/angle strategy** (Q&A): since there is a single XD
   setting, either (a) treat Hamspec as effectively **fixed-format** for this
   setup and provide a single-setting `angle_fits` anchored by `m·λ`, or (b)
   plan to add exposures at other GTILTRAW values to fit a real order-vs-angle
   relation.
3. Cross-check the new archive by re-running the Report-03 reduction and
   inspecting the `WaveCalib` QA (per-order RMS should land near the ~1–13 mÅ
   the legacy fits achieve).

## Figures
| File | Description |
|------|-------------|
| `04_fig1_coverage.png` | Wavelength coverage vs physical order (IRAF spans + XIDL points). |
| `04_fig2_dispersion_lines.png` | Dispersion (Å/pix) and ThAr lines/order vs order. |
| `04_fig3_example_order.png` | Example extracted XIDL ThAr order on its vacuum λ grid. |

---

## Q&A (questions for JXP)

17. **`angle_fits` with a single XD setting.** We only have data at one
    cross-disperser angle (GTILTRAW = 9194). Should I (a) build a
    single-setting / effectively fixed-format `angle_fits` anchored by `m·λ`
    (good enough if observers always use this setting), or (b) hold off until
    we have ThAr at additional GTILTRAW values to fit a real order-vs-angle
    relation? Do you know if Hamspec observers vary the cross-disperser?
18. **Which reduction to base the composite arc on?** XIDL is turnkey (64
    orders, vacuum, via `xidl_esihires`) but misses the 7 reddest + 13 bluest
    orders; IRAF is complete (102 orders) but needs MULTISPE/2D parsing. Build
    from XIDL now and extend with IRAF later, or invest in the full IRAF set up
    front?
19. **Order numbering.** Both legacy reductions number the physical orders the
    same way (m = 58…159), with `m·λ ≈ 5.71e5 Å`. Can you confirm this is the
    true spectral order number (not an instrument-specific offset) so we record
    the correct `m` in the PypeIt archive?
20. **`keck_nires` shift code.** The prompt mentions the dev-suite code used to
    derive the NIRES J-band shift "to solve for the shift" — should the
    Hamspec archive include an analogous per-setup pixel-shift solve, or is the
    single-setting cross-correlation sufficient?
