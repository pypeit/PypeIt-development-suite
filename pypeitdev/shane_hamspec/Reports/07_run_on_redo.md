# Report 07 — Running PypeIt on the Cooke dataset (a different detector era)

**Date:** 2026-06-27
**Author:** JXP and Claude
**Data:** `pypeitdev/shane_hamspec/Cooke_data/redo_from_scratch/`
**Task:** WaveCal #13 — build an archival wavelength template from the Cooke
XIDL outputs and run PypeIt on that dataset.
**Outputs:** `Cooke_data/redo_from_scratch/PypeIt/` (setup file + Cooke
archive); code in `scripts/build_angle_arxiv.py` (`cooke` source).

---

## TL;DR

The Cooke wavelength solution is excellent and I built an archival template from
it. But running PypeIt on the Cooke *raw data* is **blocked by a fundamental
fact: the Cooke 2010 data is a different instrument era / detector** than the
2014 dev-suite data the `shane_hamspec` class was built for. Supporting it is a
substantial, separate piece of work (a new detector definition + header-era
handling) that also cannot share the current 2014 test's archive. **A direction
decision is needed (Q&A 41–43) before the reduction can run.**

## 1. Archival template generated (done)

Built a composite-arc + angle_fits archive from the Cooke `Arc_01_fit.idl` at
its **native binning** (no resampling), via `build_angle_arxiv.py`
(`ARC_SOURCE='cooke'`, `specbin=1`):

- `Cooke_data/redo_from_scratch/PypeIt/shane_hamspec_composite_arc.fits`
- `Cooke_data/redo_from_scratch/PypeIt/shane_hamspec_angle_fits.fits`
- 83 orders (m=65–147), 80 with a real arc, nspec=2056, **~36 lines/order** —
  the rich solution from Report 06.

It is written to the Cooke folder (not PypeIt's `reid_arxiv/`) so the **2014
dev-suite test keeps its working XIDL archive**. (The two cannot share one
archive — see §3.)

## 2. `pypeit_setup` ran — with era warnings (PypeIt/ written)

`pypeit_setup -s shane_hamspec -r Raw/ -c all -d PypeIt` produced
`PypeIt/shane_hamspec_A/shane_hamspec_A.pypeit`, but with problems that reveal
the era mismatch:

- **Header keys "Bad": `mjd`, `airmass`, `xdangle`.** Our `compound_meta`
  reads `DATE` for the MJD, but the 2010 frames use **`DATE-OBS`**; and
  `xdangle←DHEITVAL`, but 2010 frames have **no `DHEITVAL`** (only `DHEITRAW`).
- **Flats untyped.** `check_frame_type` routes pixelflat/trace by the 2014
  plate names (`800:6.0`/`640:5.0`), but the Cooke plate is **`800:2.5`**, so
  the trace/pixel flats (d111–138) are left uncommented-out / untyped.
- Arcs (d139–141, ThAr) and the science/standard frames typed OK.

## 3. The blocker: the Cooke data is a *Loral 2Kx2K* detector

The decisive finding. Comparing the two datasets:

| | **2014 dev-suite** (RAW_DATA) | **2010 Cooke** (redo_from_scratch) |
|---|---|---|
| sensor | e2v CCD203-82, **4Kx4K** | **Loral 2Kx2K** (`DSENSOR`) |
| raw frame | 4136 × 4160 | **1485 × 2080** |
| amplifiers | **2** (split in columns) | **1** (`AMPSROW=1`) |
| spectral pix | 4096 | ~2048 (+32 overscan) |
| date card | `DATE` | `DATE-OBS` |
| cross-disp meta | `DHEITVAL` (µm) | `DHEITVAL` **absent**; only `DHEITRAW` |
| aperture plate | 800:6.0 / 640:5.0 | 800:2.5 |
| epoch | 2014-07 | 2010-11 |

`shane_hamspec.get_detector_par` returns a 2-amp datasec
`['[:, 1:2048]', '[:, 2049:4096]']` — **invalid for a 2080-column, single-amp
frame** (the second amp section runs off the array). So `run_pypeit` cannot
process the Cooke images; it would fail in image processing before reaching
wavelength calibration. **The reduction was therefore not run.**

(Aside: the XIDL `hamspec_setup` labels **ECH angle = DHEITRAW** and
**XD angle = GTILTRAW** — the *opposite* of our current
`echangle←GTILTRAW`, `xdangle←DHEITRAW` meta mapping. Worth revisiting; it has
not mattered so far because the reidentification self-corrects, but it should be
set correctly.)

## 4. What running on the Cooke data would require

Supporting the Cooke (2010, Loral) era in `shane_hamspec` means:
1. **A second detector definition** (Loral 2Kx2K: single amp, ~2048+32 overscan,
   its own gain/RON — values not in hand), selected by frame size in
   `get_detector_par`.
2. **Header-era handling:** `mjd`/`airmass` from `DATE-OBS` when `DATE` is
   absent; `xdangle` from `DHEITRAW` when `DHEITVAL` is absent.
3. **Plate-based frame typing** that also recognizes the 2010 plate (`800:2.5`).
4. A **binning-matched archive** (the native-2048 Cooke template built in §1),
   which **conflicts** with the 4096 XIDL archive the 2014 test needs — they
   can't both live under the one `get_echelle_angle_files()` name.

This is essentially adding a second instrument era — a real chunk of work, and
the archive conflict means a decision about which dataset is canonical.

## 5. Status
- Cooke archival template: **built** (§1).
- Cooke `.pypeit`: **generated** (§2), but with untyped flats.
- Reduction: **not run** — blocked on the Loral detector (§3).
- Production 2014 archive + dev-suite test: **untouched / still passing**.

---

## Q&A — round 8 (new, for JXP)

41. **Which dataset is the canonical dev-suite — 2014 (e2v 4Kx4K) or 2010 Cooke
    (Loral 2Kx2K)?** They are different detectors and cannot share one archive
    or one detector definition. If we *move* to Cooke, the 2014 setup/test is
    retired; if we *add* Cooke, `shane_hamspec` becomes multi-era.
42. **OK to add a Loral 2Kx2K detector era to `shane_hamspec`?** I can add a
    frame-size-based `get_detector_par` branch + `DATE-OBS`/`DHEITRAW` header
    handling + `800:2.5` plate typing. I need the **Loral 2Kx2K gain and read
    noise** (not in the headers) — do you have them, or shall I estimate?
43. **`echangle`/`xdangle` mapping.** XIDL's `hamspec_setup` uses
    **ECH=DHEITRAW, XD=GTILTRAW** — the reverse of what I implemented (Q21).
    Should I flip our meta to match XIDL?
