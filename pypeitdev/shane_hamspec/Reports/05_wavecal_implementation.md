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

> **CORRECTION (see §5).** §3 below originally concluded the archive and data
> were *incompatible* configurations. A direct per-order cross-correlation
> (§5, added after JXP's pushback) shows that is **wrong**: the orders match
> very well (cc up to 0.85, ~3–4 px shift, no stretch). The header settings
> differ, but the echellograms are nearly identical. The real reason the full
> `reidentify` failed is method/tuning + order assignment, not the data. Read
> §3 as "what the headers say," §5 as "what the spectra actually show."

## 3. Header configurations differ (but see §5)

`scripts/check_arxiv_config.py` compares the instrument setting of the legacy
ThAr (which the archive is built from) to the dev-suite raw arcs:

| Source | DATE | GTILTRAW (echangle) | DHEITRAW (xdangle) | DHEITVAL | PLATENAM |
|--------|------|:------:|:------:|:------:|:------:|
| **Archive** (IRAF/XIDL ThAr) | 2024-08-08 | **9308** | **3780** | 19385 | 960:5.0 |
| **Raw data** d2084/d2085 | 2014-07-25 | **9194** | **3970** | 20359 | 640:5.0 |

The legacy solution was taken at a different grating tilt (9308 vs 9194) and
dewar height (3780 vs 3970), a decade apart. My initial inference was that this
makes the per-order spectra non-overlapping — **but §5 shows that inference is
wrong**: the resulting echellograms differ by only ~3–4 px. The header deltas
are real but do not prevent matching.

## 4. Secondary issues found
- **Speed.** With `ech_fixed_format=False`, every observed order is
  cross-correlated against **all 64** archive orders (≈ 2 s each). With the
  **171** traced orders from Report 03 that is ~64×171 ≈ 11 000 correlations
  (hours). A fixed-format setup (pre-assigned `order_spat_pos` + `orders`, as
  NIRES does) would compare each order to only its counterpart — but it needs a
  reliable order count, and edge tracing currently over-finds 171 (Q&A 15).
- **Over-tracing.** The 171 traced orders (vs ~100 real) means many spurious
  slits that can never match; this both slows and pollutes the reidentification.

## 5. UPDATE — per-order cross-correlation: the orders **do** match

Prompted by JXP (the §3 "incompatible" inference was suspect), I cross-
correlated individual **red-side** orders of the *observed* extracted arc
(rebuilt from the `Slits`/`Arc` calibrations) against the **XIDL** archive,
using `wvutils.xcorr_shift` (fast scan) then `wvutils.xcorr_shift_stretch`
(shift + stretch refine). Module: `scripts/xcorr_orders.py`.

**Result — strong, consistent matches:**

| XIDL order m | best observed order # | m + index | refined cc | shift (pix) | stretch |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 80  | 120 | 200 | 0.70 | 3.6 | 0.9998 |
| 90  | 110 | 200 | 0.68 | 3.3 | 1.0000 |
| 110 |  90 | 200 | 0.57 | 4.6 | 0.9988 |
| 120 |  80 | 200 | 0.58 | 3.3 | 0.9997 |
| 130 |  70 | 200 | 0.85 | 3.0 | 1.0000 |

![Best red-side order match](05_fig1_xcorr_best.png)

Three conclusions:

1. **The archive is usable for this data.** Every probed order matches with
   cc 0.57–0.85, a tiny **~3–4 px shift**, and **no stretch** (≈1.000). The
   overlay (XIDL m=80 vs observed order #120) shows the ThAr lines landing on
   top of each other. The header config delta (§3) produces only a few-pixel
   offset — *not* an incompatibility. **§3's original conclusion is retracted.**
2. **The order mapping is exact and linear:** physical order
   **m = 200 − (observed order index)** across the whole probed range. This is
   the order-identification relation the pipeline was missing — it lets us
   assign physical `m` (and hence the right archive order) to each traced slit
   directly.
3. **So why did the full `reidentify` fail (cc≈0.17)?** Not the data. The
   `autoid.reidentify` line-pattern cc (with `reid_cont_sub=False`,
   `sigdetect`, `percent_ceil`, and all 64 arxiv orders tried per slit, diluted
   by the 171 spurious traces) scores much lower than the direct spectrum
   cross-correlation here (0.7–0.85). The fix is method/tuning + order
   assignment, not new data.

### Revised recommended path
1. **Keep the XIDL archive** — it works.
2. **Assign orders up front** using `m = 200 − index` (equivalently a
   fixed-format `order_spat_pos`/`orders` once the tracing is cleaned), so each
   observed order is reidentified against its single correct archive order.
3. **Tune `reidentify`**: set `reid_cont_sub = True`, relax/curate `cc_thresh`
   to ~0.5, and revisit `sigdetect`/`percent_ceil` so the line-pattern cc
   reflects the real ~0.7 correlation.
4. **Trim the order over-tracing** (171 → ~100) so spurious slits stop
   polluting/slowing the match.
5. Then extend the archive XIDL→IRAF (64→102 orders) as planned (Report 04).

---

## Q&A (questions for JXP)

*(Q23–25 from the previous version are withdrawn: Q23's "config mismatch
blocker" is resolved by §5 — the orders match. The relevant questions now:)*

26. **Order mapping.** §5 finds `m = 200 − (traced order index)` for this
    setup. Is hard-coding/encoding this relation (or the equivalent
    `order_spat_pos`) acceptable, given the dewar height can move the orders
    (then the constant 200 would shift)? Should we instead derive the offset
    per-frame from the cross-correlation (robust to dewar moves)?
27. **`reidentify` tuning vs fixed-format.** Prefer I (a) tune the existing
    non-fixed `reidentify` (`reid_cont_sub=True`, `cc_thresh≈0.5`) and let it
    find orders, or (b) implement `ech_fixed_format=True` with
    `order_spat_pos`/`orders` for a deterministic per-order match? (b) is
    faster and more robust but needs the cleaned order count.
28. **Order over-tracing.** Still want the tracer trimmed from 171 to the true
    count (~100)? That helps both reidentification and the final 2D solution
    (was Q15).
