# 10 — Bright-object (forced slit-center) extraction test

**Date:** 2026-08-20
**Task:** Development / Finishing up #3 — "When a very bright object is in the
slit, the object finding fails. Let's explore how best to fix this. Do a test
run and then ask me questions below."

## Purpose

First end-to-end test of the forced slit-center extraction implemented for
Extraction #1 (`[reduce][findobj] force_center_obj = True`, sky subtraction
off via `skip_skysub = True` + `no_local_sky = True`).  That code had only a
synthetic unit test; the real-data verification was blocked because the
9-frame Cooke push to `RAW_DATA/shane_hamspec/Hamilton/` was never mirrored
to Google Drive (Q&A 48) and is absent from this machine.

## Test setup

- **Data:** the seven **2014 e2v** frames in `RAW_DATA/shane_hamspec/`
  (d2000–d2100).  The 300 s **HD188209** (O9.5Iab) standard is the very
  bright star.  The Cooke/Loral frames (the ~7 px orders where object
  finding actually failed) are **not available** — see Caveats.
- **Run:** `REDUX_OUT/shane_hamspec/2014_brightobj/shane_hamspec_A/`
  (`pypeit_setup` → 1 clean config, correct frame typing; then `run_pypeit`,
  `pypeit14b` env, ~17 min, exit 0).
- **Wavelength archive:** the production `reid_arxiv` is now the **Cooke**
  (Loral 2048-px) archive and cannot calibrate the 4096-px e2v frames.  I
  rebuilt the **XIDL**-based archive with `build_angle_arxiv.py`
  (`ARC_SOURCE='xidl'`, anchors at the 2014 flipped angles
  `ECH_ANCHOR=3970` (DHEITRAW) / `XD_ANCHOR=9194` (GTILTRAW)) into the run
  directory, which PypeIt prefers over the package copy (bare-filename
  lookup resolves against the CWD first).  The package archive was not
  touched.  Note: the script's `from pypeit.core.wavecal import templates`
  now needs a shim to `pypeitdev/wavecal/templates.py` (the module moved out
  of the PypeIt package).

## Results

### Wavelength calibration — much better than before

**63/99 orders solved (m=65–146), median per-order RMS 0.278 px** — the
*blue orders now solve* (m=101–146 at RMS 0.06–0.75 px).  Compare the last
2014 run (WaveCal #9): 27 orders, red half only (m=65–97), same median RMS.
The order identification is confirmed correct: m·λ = 5.711–5.712×10⁵ Å,
constant to ~0.02% across all 63 orders.  Nothing in our code changed the
matching; the branch merged `origin/develop` since WaveCal #11, so the
develop-side echelle wavecal improvements appear to have fixed the
blue-order matching failure we fought in Q36/Q37 (which is now moot for the
e2v era).

The default gate `rms_thresh = 0.3 px` (`RMS/FWHM = 0.1`) then masks 30 of
the 63 solved orders for the reduction, leaving **33 extracted orders**
(the masked 30 have RMS 0.31–1.04 px).

### Forced extraction — works end-to-end

`FindObjects` logged "Skipping global sky sub as per user request" and
"Forced one object at the center of each of the 33 orders
(force_center_obj=True)"; the reduction then wrote **spec1d + spec2d** for
HD188209 — the spec1d that peak-detection object finding never produced.

- 33 SpecObjs (one per good order), all with **boxcar and optimal**
  extractions; wavelength coverage 3858–8147 Å (vac).
- Median per-pixel S/N **89–185** per order (optimal slightly above boxcar).
- Boxcar radius ≈ 5.5–6.6 px (half the order width, as designed); the fitted
  object FWHM ≈ 3.5 px.
- The forced center trace is within **−1.1 to +1.0 px** of the star's true
  flux-weighted centroid (≤9% of the ~12 px order width, a smooth trend
  with order — spectrograph distortion), so aperture losses are negligible.

![Extracted orders](10_fig1_extracted_orders.png)

![Spatial profiles vs forced trace](10_fig2_spatial_profiles.png)

### A key observation

On the **e2v** era the star does **not** fill the order: FWHM ≈ 3.5 px in
~12 px-wide orders (Fig 2 shows a clear peak with sky pixels either side).
The "smashed profile has no peak" failure that motivated the forced
extraction is specific to the **Loral** era's ~7 px orders (Report 09).
Consequences: (a) on e2v data the normal `ech_objfind` would likely work and
would center the trace on the star; (b) sky pixels *do* exist on the e2v
orders, so sky subtraction is not impossible there; (c) the forced-center
mechanism nevertheless works fine on e2v (≲1 px offset, full-order boxcar).

## Caveats

- The **actual failing case** — the Loral ~7 px orders where the star fills
  the slit — is still **unverified**: the 9 Cooke frames exist on no drive
  reachable from this machine and were never mirrored to Google Drive
  (Q&A 48).  This run verifies the mechanism, not the original failure.
- The dev-suite `pypeit_test reduce -i shane_hamspec` could not be used (its
  `.pypeit` points at the missing `Hamilton/` Cooke frames); this was a
  manual `run_pypeit` on the 2014 set.

## Q&A

New questions 49–53 are in `claude_prompts/shane_hamspec.md`
(Development / Finishing up § Q&A).

---

## UPDATE (2026-08-24, Finishing-up #4): Q&A 49–53 implemented; Loral case verified

JXP answered Q49–53 and re-downloaded the full Cooke dataset to
`pypeitdev/shane_hamspec/Cooke_data/`.  Actions taken:

- **Cooke data restored (Q53).**  All 60 raw frames are in
  `Cooke_data/redo_from_scratch/Raw/`; the 9 dev-suite frames (arcs
  d139/140, trace d133–135, pixelflat d111–113, standard d144) were copied
  back to `RAW_DATA/shane_hamspec/Hamilton/`.  (Still needs mirroring to
  Google Drive.)
- **Q49/Q50 (no change):** `force_center_obj=True` stays unconditional; the
  trace stays at the geometric order center.
- **Q51 — RMS gate raised:** `rms_thresh_frac_fwhm` 0.1 → 0.25 (~0.75 px at
  the ~3 px arc FWHM) in `shane_hamspec.default_pypeit_par`.
- **Q52 — per-era archives shipped.**  The Cooke archive was renamed
  `shane_hamspec_loral_{angle_fits,composite_arc}.fits` and a new e2v
  archive (XIDL `Arc_01_fit.idl`, anchors ECH=3970/XD=9194) was built and
  shipped as `shane_hamspec_e2v_*`.  A new compound meta `detector`
  (`Loral2Kx2K`/`e2v4Kx4K` from the raw frame size, mirroring
  `get_detector_par`) selects the era in `get_echelle_angle_files`, which
  now takes an optional `meta_dict` (the core `wavecalib.py` call site
  passes the arc-frame metadata; the base class and the keck_hires,
  keck_nirspec, and vlt_uves overrides accept-and-ignore the new kwarg).
  `build_angle_arxiv.py` gained an era CLI (`build_angle_arxiv.py
  [loral|e2v]`); its `loral` output matches the shipped archive to
  float-noise (≤1e-11 Å).  The scripts also now fall back to the dev-suite
  `pypeitdev/wavecal/templates.py` when the (removed) package module is
  absent.

### Verification

**Loral / dev-suite (the original failing case):**
`pypeit_test reduce -i shane_hamspec` **PASSES 1/1 (~5 min)** with the
shipped loral archive selected via the `detector` meta.  Wavecal unchanged
(19 orders @ 0.147 px, all below even the old gate).  **HR 7373 — the star
that fills the ~6 px Loral orders — is now extracted**: 19 SpecObjs, boxcar
+ optimal, boxcar spans the full order (BOX_R≈2.9 px), FWHM≈2.0–2.5 px,
median per-pixel S/N ≈ 42–63 (boxcar) / 38–101 (optimal), 4757–5700 Å.
This closes the Q47/Extraction-#1 loop end-to-end.

![Loral spatial profiles](10_fig3_loral_profiles.png)

![Loral extracted orders](10_fig4_loral_extracted.png)

**e2v / 2014 (fresh calibrations, shipped archives only):** run completes;
the e2v archive is auto-selected; **64/99 orders solved (m=65–146, median
RMS 0.284 px)** and, with the raised gate, **all 64 are extracted** (vs 33
before), optimal S/N 100–186/pix, 3858–8899 Å.

Note: a first re-run reused the on-disk calibrations (run_pypeit `-o`
overwrites science outputs but reuses existing calibration files), silently
bypassing both changes — the definitive run required deleting
`Calibrations/` first.
