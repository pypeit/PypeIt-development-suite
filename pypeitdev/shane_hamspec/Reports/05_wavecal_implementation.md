# Report 05 — Wavelength-calibration implementation (NIRES-style)

**Date:** 2026-06-17
**Author:** JXP and Claude
**Task:** Development / Wavelength Calibration #3 — implement a wavelength
calibration for Shane/Hamspec and see how far it gets.
**Scripts:** `pypeitdev/shane_hamspec/scripts/` (`pypeit14` env)

Builds on [Report 04](04_wavelength_deep_dive.md). The plan there (Q&A 20) was
the Keck/NIRES recipe: a single composite-arc `reid_arxiv` + cross-correlation
reidentification. This report implements it and reports the outcome.

---

## 1. What was implemented

1. **Built the `reid_arxiv` composite arc** `shane_hamspec.fits`
   (`scripts/build_reid_arxiv.py`) from the XIDL solution via
   `templates.xidl_esihires` — 64 orders (m = 65–146), columns
   `wave`/`flux`/`order`, same structure as `keck_nires.fits`. Written into
   `pypeit/data/arc_lines/reid_arxiv/`.
2. **Switched the wavecal recipe** in `shane_hamspec.default_pypeit_par` to the
   NIRES route: `method='reidentify'`, `reid_arxiv='shane_hamspec.fits'`
   (keeping the echelle refit params).
3. **Re-mapped the angle meta (Q&A 21):** `echangle ← GTILTRAW` (grating tilt,
   sets wavelength) and `xdangle ← DHEITRAW` (dewar height / cross-disperser,
   sets order placement); added `echangle` to `configuration_keys` so both the
   grating tilt and dewar height define a configuration.
4. **Disabled the leftover debug `embed` (Q&A 22)** in `wavecalib.py`
   (`method='echelle'` branch) so that path can't halt a reduction.

`pypeit_setup` still yields **one configuration**; the reduction was launched
with `pypeit_test reduce -i shane_hamspec`.

## 2. How far it got

Big step forward from Report 03: the reduction now runs the **whole
calibration chain and enters wavelength calibration** — edges, slits, flats,
tiltimg, arc extraction, then `reidentify` actually executes (the 404 and the
`embed` are both gone). The `WaveCalib` step begins cross-correlating each
observed order against the 64 archive orders.

**But the reidentification does not match.** Across the cross-correlations that
ran (before I stopped it — see §3), the correlation coefficient was:

```
n = 124 attempts,  cc:  min 0.115,  median 0.172,  max 0.382
matches with cc ≥ cc_thresh (0.6):  0
```

Every arxiv order is rejected (`cc < cc_thresh`) for every observed order — no
wavelength solution is produced.

## 3. Why it fails — the archive and the raw data are *different
configurations* (key finding)

`scripts/check_arxiv_config.py` compares the instrument setting of the legacy
ThAr (which the archive is built from) to the dev-suite raw arcs:

| Source | DATE | GTILTRAW (echangle) | DHEITRAW (xdangle) | DHEITVAL | PLATENAM |
|--------|------|:------:|:------:|:------:|:------:|
| **Archive** (IRAF/XIDL ThAr) | 2024-08-08 | **9308** | **3780** | 19385 | 960:5.0 |
| **Raw data** d2084/d2085 | 2014-07-25 | **9194** | **3970** | 20359 | 640:5.0 |

The legacy wavelength solution was taken at a **different grating tilt** (9308
vs 9194 → a different wavelength zero-point per order) **and a different dewar
height** (3780 vs 3970 → the cross-disperser is in a different place, so the
orders land at different spatial/spectral positions), a decade apart. So the
archive's per-order ThAr spectra simply do not line up with the observed ones,
and the cross-correlation cc stays at the noise level (~0.17). This — not a
coding problem — is why reidentification fails.

This is exactly the multi-setting issue anticipated in Report 04 §6 / Q&A 17:
the dewar height (cross-disperser) **does** vary between data sets, and here the
only legacy solution we have is for a setting that no raw dev-suite frame
matches.

## 4. Secondary issues found
- **Speed.** With `ech_fixed_format=False`, every observed order is
  cross-correlated against **all 64** archive orders (≈ 2 s each). With the
  **171** traced orders from Report 03 that is ~64×171 ≈ 11 000 correlations
  (hours). A fixed-format setup (pre-assigned `order_spat_pos` + `orders`, as
  NIRES does) would compare each order to only its counterpart — but it needs a
  reliable order count, and edge tracing currently over-finds 171 (Q&A 15).
- **Over-tracing.** The 171 traced orders (vs ~100 real) means many spurious
  slits that can never match; this both slows and pollutes the reidentification.

## 5. Status & recommended next steps
- The **machinery is in place and runs**: composite arc, NIRES `reidentify`
  recipe, corrected angle meta, no debug halt. What's missing is a wavelength
  template **for the configuration of the dev-suite data**.
- **The cleanest fix is data:** either (a) add a ThAr arc taken at the raw
  data's setting (GTILTRAW≈9194, DHEITRAW≈3970) and build the archive from it,
  or (b) replace the dev-suite raw set with frames matching the legacy
  solution's setting (GTILTRAW≈9308, DHEITRAW≈3780, plate 960:5.0). See Q&A 23.
- **Or bootstrap a solution for the 2014 arcs** with `pypeit_identify` on
  d2084/d2085 (seeded by `m·λ ≈ 5.71e5` and the line list), then archive that.
  More effort but self-contained. See Q&A 24.
- Separately, **reduce the order over-tracing to ~the true count** and consider
  a **fixed-format** setup (`order_spat_pos`, `orders`) to make
  reidentification fast and unambiguous (Q&A 25).

---

## Q&A (questions for JXP)

23. **Config mismatch (blocker).** The legacy ThAr solution
    (GTILTRAW=9308, DHEITRAW=3780, plate 960:5.0; 2024) does **not** match the
    dev-suite raw data (GTILTRAW=9194, DHEITRAW=3970, plate 640:5.0; 2014). Can
    you provide a ThAr arc + solution **at the 2014 setting**, or swap the
    dev-suite raw data to the 2024 setting so the archive applies?
24. **Bootstrap instead?** Alternatively, should I derive a fresh solution for
    the existing 2014 arcs (d2084/d2085) with `pypeit_identify` (seeded by the
    legacy `m·λ` and order numbering) and archive *that*? This keeps the
    current dev-suite data but is more hands-on.
25. **Fixed-format + order count.** To make reidentification fast and robust,
    should we move to `ech_fixed_format=True` with explicit `order_spat_pos`
    and `orders` (like NIRES), once we settle the true order count (the edge
    tracer currently finds 171; Q&A 15)?
