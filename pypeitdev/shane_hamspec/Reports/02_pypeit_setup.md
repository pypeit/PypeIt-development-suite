# Report 02 — `pypeit_setup` on the Shane/Hamspec data set

**Date:** 2026-06-16
**Author:** JXP and Claude
**Data:** `PypeIt-development-suite/RAW_DATA/shane_hamspec/` (7 frames)
**Task:** Prep #2 — run `pypeit_setup` on the main data set and report findings.

Builds on [Report 01](01_raw_data_inventory.md) and the answered Q&A in the
prompt (`claude_prompts/shane_hamspec.md`). Key answers applied here: the class
stays **`shane_hamspec`**; **pixelflats = wide plate (800:6.0)**, the
narrow-plate (640:5.0) flats are **traceflats**; **HD188209 is the standard**.

---

## 1. Command run

```bash
conda activate pypeit14
pypeit_setup -s shane_hamspec \
  -r .../RAW_DATA/shane_hamspec/ \
  -c all -d .../pypeitdev/shane_hamspec
```

It **completed (exit 0)** — the Report-01 `msgs`→`logging` fix is what allows
the CLI to import at all. It wrote one `.pypeit` (+ `.calib`) per detected
configuration: `shane_hamspec_A` … `shane_hamspec_E`.

---

## 2. What `pypeit_setup` produced

It found **5 configurations** and typed the frames as follows:

| Config | decker (PLATENAM) | filter1 (DFILTNAM) | file | OBJECT | PypeIt frametype |
|:------:|:-----------------:|:------------------:|------|--------|------------------|
| A | 800:6.0 | BG13 | d2000 | pixflat    | pixelflat,illumflat,trace |
| B | 800:6.0 | BG12 | d2040 | pixflat    | pixelflat,illumflat,trace |
| C | 640:5.0 | Open | d2084 | arcs       | arc,tilt |
| C | 640:5.0 | Open | d2100 | HD188209   | **None (commented out!)** |
| D | 640:5.0 | BG12 | d2085 | arcs       | arc,tilt |
| D | 640:5.0 | BG12 | d2095 | blazeflats | pixelflat,illumflat,trace |
| E | 640:5.0 | BG13 | d2090 | blazeflats | pixelflat,illumflat,trace |

All frames also produced a warning: `Bad Header key (airmass)`.

**The data set is not reducible as-is.** No single configuration contains a
complete calibration set (arc **and** trace-flat) together with the standard.
There are four distinct, related problems, all in
`pypeit/spectrographs/shane_hamspec.py`.

---

## 3. Problems found (with root causes)

### Problem 1 — Over-split into 5 configurations *(highest priority)*
`configuration_keys()` returns `['decker', 'filter1', 'xdangle', 'binning']`.
Both `decker` (PLATENAM aperture plate) and `filter1` (DFILTNAM blocking
filter) vary across frames that all belong to the **same** spectroscopic setup:

- `decker` differs *by design* — pixelflats use the **wide** plate (800:6.0)
  while arcs/traceflats/standard use the **science** plate (640:5.0). Splitting
  on it permanently separates the pixelflats from everything else.
- `filter1` is just the order-blocking filter: arcs were taken `Open`, flats
  `BG12`/`BG13`, the standard `Open`. Splitting on it scatters arcs, flats and
  the standard into different configs.
- `xdangle` (GTILTRAW) is **constant (9194)** here, and `binning` is fixed
  (`1,1`), so those two alone would yield **a single configuration** — which is
  what we want.

There is also a self-contradiction in the source: the comment directly above
the return reads *"decker is not included because arcs are often taken with a
0.5" slit"*, yet `'decker'` **is** in the returned list. The comment describes
the intended behavior; the code does the opposite.

**Recommended fix:** `configuration_keys()` → `['xdangle', 'binning']` (drop
`decker` and `filter1`). That collapses the data to one config containing all
calibrations + the standard. (HIRES, the template this class was copied from,
keys on cross-disperser angle + binning, not decker/filter.)

### Problem 2 — Standard/science never typed (HD188209 dropped)
`d2100` (the standard) was typed `None` and **commented out**. Root cause is a
**case mismatch** in `lamps(fitstbl, 'off')`: it tests the lamp value against
the lowercase string `'off'`, but the header card is `LAMPPOS = 'Off'`
(capital O). So `('Off' == 'off')` is `False`, `check_frame_type` returns False
for both `science` and `standard`, and the frame is discarded. (The `arcs`
test matches `'Thorium-Argon'` and the `flat` test matches `'PolarQuartz'`
exactly, which is why those frames *do* type correctly.)

**Recommended fix:** make the off-test case-insensitive, e.g. compare
`str(val).lower()` against `{'off','none'}`.

Secondary issue once that is fixed: `d2100` is 300 s, but
`standardframe['exprng'] = [1, 61]`, so it would type as **science**, not
**standard** (science exprng is `[61, None]`). Since HD188209 is the intended
standard, we need a way to mark it as such — e.g. widen the standard exprng
and/or distinguish standard from science some other way (a bright-star /
object-name heuristic, or simply let the user re-type it). Flagged in Q&A.

### Problem 3 — pixelflats vs traceflats not distinguished
Per the answered Q&A, the **wide-plate (800:6.0)** flats are the true
**pixelflats**, and the **narrow-plate (640:5.0)** "blazeflats" are
**traceflats** (trace/illumflat only — *not* pixelflat). Currently
`check_frame_type` types **every** `PolarQuartz` frame as
`pixelflat,illumflat,trace`, because it keys only on lamp status. So the
blazeflats are wrongly tagged `pixelflat`, and (because of Problem 1) the
pixelflats sit in their own configs anyway.

**Recommended fix:** route on `PLATENAM`/decker inside `check_frame_type`:
`pixelflat` ← `PolarQuartz` AND wide plate (800:6.0); `trace,illumflat` ←
`PolarQuartz` AND narrow plate (640:5.0).

### Problem 4 — `airmass` missing (warning only)
`init_meta` maps `airmass` to header card `AIRMASS`, which **does not exist** in
these files (there is `HA`, `RA`/`DEC`, `DATE`, but no AIRMASS/elevation card).
Hence the `Bad Header key (airmass)` warning on every frame. It is non-fatal
(setup continues) but airmass is needed downstream (extinction, telluric). It
can be computed via a `compound_meta` from RA/DEC + MJD + the Shane telescope
location.

---

## 4. Net effect / what a correct setup should look like

With Problems 1–3 fixed, `pypeit_setup` should yield **one configuration**
containing:

| frametype | files |
|-----------|-------|
| pixelflat (+illumflat) | d2000, d2040  (wide plate 800:6.0) |
| trace, illumflat       | d2090, d2095  (narrow plate 640:5.0 "blazeflats") |
| arc, tilt              | d2084, d2085 |
| standard               | d2100  (HD188209) |

That is a complete, reducible echelle setup (modulo the standard-vs-science
typing decision in Problem 2).

---

## 5. Carry-forward / next steps
- These four fixes belong to the **Development** phase (editing
  `shane_hamspec.py`): `configuration_keys`, `lamps` case-insensitivity,
  `check_frame_type` plate-based flat routing, and an airmass `compound_meta`.
- After the fixes, re-run `pypeit_setup -c all` and confirm a single config;
  that single `.pypeit` becomes the basis for the first `run_pypeit` attempt and
  for the eventual dev-suite `pypeit_files/shane_hamspec_*.pypeit`.
- Generated scratch configs (`shane_hamspec_A…E/`) live under
  `pypeitdev/shane_hamspec/` and are git-ignored; they can be deleted/regenerated.

---

## Q&A (questions for JXP)

1. **Configuration keys:** OK to set `configuration_keys()` →
   `['xdangle', 'binning']` (drop `decker` and `filter1`) so all frames land in
   one configuration? (This matches the existing in-code comment's intent.)
2. **Standard vs science typing:** HD188209 is 300 s, which collides with the
   science exptime range. How should we tag it as `standard` automatically —
   widen `standardframe['exprng']`, use an object-name/bright-star heuristic, or
   just rely on the user editing the `.pypeit` frametype? (For the dev-suite
   test there are no faint science frames yet, so treating the lone hot star as
   the standard is fine.)
3. **Airmass:** OK to compute `airmass` via `compound_meta` from RA/DEC + MJD +
   Shane location (since no `AIRMASS` card exists), or leave it `None` for now?
4. **illumflat source:** should the **pixelflats** (wide plate) or the
   **traceflats** (narrow plate) be used for `illumflat`? Currently both carry
   the `illumflat` tag; for an echelle the slit-illumination usually comes from
   the science-plate (trace) flats. Please confirm.

---

## UPDATE — 2026-06-17 (Prep #3: fixes implemented + airmass computed)

JXP implemented the fixes for Problems 1–3; this update records the result and
the **airmass** work (Problem 4), which was assigned to Claude.

### Resolution of Problems 1–3 (per the answered Q&A)
- **Config keys (Q6):** `configuration_keys()` is now `['xdangle','binning']`.
- **Standard exprng (Q7):** `standardframe['exprng']` widened to `[1, 301]`.
- **Plate-based flat typing (Q9):** `check_frame_type` now types the wide
  plate (800:6.0) as `pixelflat` and the narrow plate (640:5.0) as
  `trace,illumflat` — confirming **illumflat comes from the narrow
  (science) plate**, as expected for echelle.

### Airmass (Problem 4 / Q8) — implemented by Claude
- The raw headers have **no `AIRMASS` card**, so airmass is computed with
  `astropy` from the pointing (`RA`/`DEC`), the readout time (`DATE`), and the
  Shane/Lick observatory location.
- The calculation lives in a **reusable common helper**,
  `pypeit.core.meta.airmass(ra, dec, obstime, longitude, latitude, elevation)`,
  so any spectrograph can call it from `compound_meta`. It builds an
  `EarthLocation`, transforms the target `SkyCoord` into the local `AltAz`
  frame at `obstime`, and returns `sec z` (NaN if below the horizon) — i.e. it
  follows the astropy observation-planning example.
- `shane_hamspec.compound_meta` now calls it, passing the
  `ShaneTelescopePar` longitude/latitude/elevation.

### Two follow-on bugs found while verifying (fixed by Claude)
Re-running `pypeit_setup` after JXP's edits initially **crashed**; two small
issues in the new typing code had to be fixed to get a working setup (both
flagged in the prompt Q&A, items 10–11):
1. **`KeyError: 'PLATENAM'`** — `check_frame_type` referenced the raw header
   card `fitstbl['PLATENAM']`, but the metadata table is keyed by the PypeIt
   meta name. Changed to **`fitstbl['decker']`** (decker ← PLATENAM).
2. **Lamp-off case** — `check_frame_type` calls `self.lamps(fitstbl, 'Off')`,
   but `lamps()` only handled lowercase `'off'` (and compared values against
   `'off'`), so it fell through to `raise ValueError`. Made the off-test
   **case-insensitive** (matches `'off'`/`'none'` regardless of case), which
   keeps the `'Off'` call working and matches the `LAMPPOS='Off'` header.

### Verified result — one configuration, correct types, airmass populated
`pypeit_setup -s shane_hamspec -c all` now exits 0 with **1 configuration**:

| file | frametype | decker | airmass |
|------|-----------|:------:|:-------:|
| d2084, d2085 | arc,tilt        | 640:5.0 | ~1.0000 |
| d2090, d2095 | illumflat,trace | 640:5.0 | ~1.0000 |
| d2000, d2040 | pixelflat       | 800:6.0 | ~1.0000 |
| d2100 (HD188209) | science     | 640:5.0 | **1.2394** |

Calibrations sit at the zenith (telescope parked; DEC ≈ Lick latitude) → airmass
≈ 1.0; the standard HD188209 gives a realistic **1.24**.

### One residual item (new Q&A 12)
`pypeit_setup` warns *"Some frames are assigned both science and standard
types; choosing the most likely type"* and tags **d2100 as `science`**, not
`standard`. With `standard=[1,301]` and `science=[61,None]`, the 300 s standard
falls in **both** windows. For the dev-suite (no faint science yet) this hot
star should be the **standard**; how would you like it disambiguated? (e.g.
narrow `science['exprng']`, or rely on the user setting the frametype in the
`.pypeit`.)
