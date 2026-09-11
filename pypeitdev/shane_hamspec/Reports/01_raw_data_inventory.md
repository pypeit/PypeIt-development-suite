# Report 01 — Lick/Hamspec Raw Data Inventory & First Look

**Date:** 2026-06-16
**Author:** JXP and Claude
**Data:** `PypeIt-development-suite/RAW_DATA/lick_hamspec/`
**Task:** Prep #1 — view the raw data with `pypeit_view_fits` and report findings.

> **Report-naming convention adopted for this project** (used for all
> subsequent reports): `NN_short-title.md`, where `NN` is a zero-padded,
> monotonically increasing index. Associated figures share the report's index
> as a prefix: `NN_figM_description.png`. All reports live in
> `pypeitdev/lick_hamspec/Reports/`.

---

## 1. How the data were viewed

The task asks to use `pypeit_view_fits`. That script drives the **ginga** GUI
viewer over an RC connection, which is not available in this headless
development session (the command returns immediately with no display). I
therefore "viewed" the frames programmatically with `astropy.io.fits` +
`matplotlib` using a ZScale stretch (the same interval ginga uses by default),
and saved PNGs alongside this report (see Figures below). The numeric
characterization (statistics, order counts, overscan levels) comes from the
same script.

**Blocker found & fixed to get this far:** importing any spectrograph on this
`hamspec` branch failed because `pypeit/spectrographs/shane_hamspec.py` still
used the old `from pypeit import msgs` API. The branch has migrated the whole
package to the stdlib-`logging` convention (`from pypeit import log`;
`raise PypeItError(...)`), and `shane_hamspec.py` was the *only* file left
behind — so it broke `from pypeit.spectrographs import *`, and hence every
PypeIt CLI script (`pypeit_view_fits` included). I migrated it:
`msgs.error(...)` → `raise PypeItError(...)`, `msgs.warn(...)` →
`log.warning(...)`. Imports now succeed. See Q&A item 1.

---

## 2. File inventory (7 frames, all from UT 2014-07-25)

| File | OBJECT | Frame type (PypeIt) | EXPTIME | LAMPPOS | DFILTNAM | PLATENAM | TEMPDET |
|------|--------|---------------------|--------:|---------|----------|----------|--------:|
| d2000.fits | pixflat    | pixelflat (flood)  |   3 s | PolarQuartz   | BG13 | 800:6.0 | −80.2 |
| d2040.fits | pixflat    | pixelflat (flood)  |  22 s | PolarQuartz   | BG12 | 800:6.0 | −79.9 |
| d2084.fits | arcs       | arc/tilt           |   3 s | Thorium-Argon | Open | 640:5.0 | −80.4 |
| d2085.fits | arcs       | arc/tilt           |  20 s | Thorium-Argon | BG12 | 640:5.0 | −75.1 |
| d2090.fits | blazeflats | trace/illum (slit) |   4 s | PolarQuartz   | BG13 | 640:5.0 | −79.0 |
| d2095.fits | blazeflats | trace/illum (slit) |  25 s | PolarQuartz   | BG12 | 640:5.0 | −79.9 |
| d2100.fits | HD188209   | science/standard   | 300 s | Off           | Open | 640:5.0 | −79.1 |

HD188209 is a bright O9.5Iab star — i.e. this is effectively a **standard /
telluric** exposure rather than a faint science target. There is **no bias**
and **no dark** frame in this set, and no on-sky science target other than the
hot star. (A full reduction will need biases or an overscan-only bias model.)

### Key header structure
- Single `PrimaryHDU`, `uint16` (BZERO=32768), **NAXIS1=4160 × NAXIS2=4136**.
- `DATASEC = '[1:4096,1:4136]'` → 4096 science columns + **64 overscan
  columns** (4097–4160). `AMPSROW=2`, `AMPSCOL=1`: **two amplifiers split along
  columns** (left/right). `DSENSOR = 'e2v CCD203-82 4kx4k thin'`.
- Time: `DATE` is the UT of readout (ISOT). Midpoint cards (`MP-MID`) also exist.
- No `BINX/BINY/CCDSUM/BINNING` card present — these data are unbinned
  (`1,1`); binning will have to default rather than be read (matches the
  existing `shane_hamspec` `binning` default of `'1,1'`).
- Useful typing/config cards: `LAMPPOS` (lamp), `DFILTNAM` (order-blocking
  filter), `PLATENAM` (aperture/decker plate), `GTILTRAW` (grating tilt,
  constant `9194` here), `OBJECT` (human label literally encodes intent:
  `pixflat` / `blazeflats` / `arcs`).

---

## 3. Detector image: counts & overscan

ZScale, sigma-clipped statistics over the 4096 science columns; overscan over
the 64 trailing columns:

| File | type | min | max | median(data) | median(overscan) |
|------|------|----:|----:|-------------:|-----------------:|
| d2000 | pixflat    | 925 | 51324 | 1563 | 954 |
| d2040 | pixflat    | 923 | 48090 | 1574 | 954 |
| d2084 | arc        | 922 | 65535 |  958 | 954 |
| d2085 | arc        | 920 | 65535 |  957 | 954 |
| d2090 | blazeflat  | 925 | 50277 | 1502 | 954 |
| d2095 | blazeflat  | 924 | 52202 | 1382 | 954 |
| d2100 | HD188209   | 923 | 28717 | 1173 | 954 |

- **Bias level ≈ 954 ADU**, stable across all frames and both amps (amp1
  median ≈ 962, amp2 ≈ 947 on d2095 → a small ~15 ADU amp-to-amp offset, so
  per-amp overscan subtraction is appropriate).
- The two arcs **saturate** (65535 = the `uint16`/`saturation` ceiling in
  `shane_hamspec`): the brightest ThAr lines are saturated even at 3 s. Arc
  exposure choice / saturated-line rejection matters for wavelength work.

---

## 4. Echelle format (the important part)

The dispersion runs **along columns (NAXIS1)** and orders are stacked **along
rows (NAXIS2)** — consistent with the existing `shane_hamspec` detector
(`specaxis=1`, `specflip=False`, `spatflip=True`). Confirmed visually:

- **Figure 4** (blazeflat d2095): dozens of horizontal echelle orders, brightest
  near the center and falling off top/bottom following the blaze + cross-disperser
  throughput.
- **Figure 5** (cross-dispersion cut): a peak-find at col=2048 gives **~50–90
  resolvable orders** with a **median order separation ≈ 16 px**, spanning rows
  ≈ 430–2673. Orders are **densely packed** (separations as small as ~8 px in
  the crowded center) — edge tracing will be demanding. The orders also **shift
  in row with column** (order curvature/tilt), so peaks at col=1000/2048/3000 do
  not line up.
- The format shows **two illuminated bands** (a bright main band, rows ≈
  1500–2700, and a fainter lower band, rows ≈ 400–1000, with a gap ≈ 1000–1500).
  This needs confirmation — it may be the full cross-dispersed format with a
  throughput dip, or vignetting. Flagged in Q&A.

### pixflat vs blazeflat — both are `PolarQuartz`, distinguished by the plate
This is the single most important data-organization finding:
- **pixflats** (`OBJECT=pixflat`, `PLATENAM=800:6.0`) are taken with a **wide
  aperture plate** that floods the detector — orders are washed out into a
  smooth illumination (**Figures 1–3**: no order structure, just a broad blaze
  hump). These are true **pixel-to-pixel flats**.
- **blazeflats** (`OBJECT=blazeflats`, `PLATENAM=640:5.0`, same plate as
  science/arc) are taken through the **science slit** and **show the orders**
  (Figure 4). These are the frames to use for **slit/order edge tracing** and
  the blaze/illumination shape.

The current `shane_hamspec.check_frame_type` types *all* `PolarQuartz` frames
together as `pixelflat`+`trace`+`illumflat` purely on lamp status — it cannot
tell these two apart. For `lick_hamspec` we will likely need to use `PLATENAM`
(or `OBJECT`) to route the wide-plate frames to `pixelflat` and the
narrow-plate frames to `trace`/`illumflat`. (See Q&A item 3.)

---

## 5. Implications for `lick_hamspec` development (carry-forward notes)

1. **Reuse `shane_hamspec` as the template.** Geometry (4160×4136, 2 amps split
   in columns, 64 overscan cols, `specaxis=1`, `spatflip=True`,
   `pypeline='Echelle'`, `ech_fixed_format=False`) all match this data. Note
   `shane_hamspec.get_echelle_angle_files()` already references
   `lick_hamspec_angle_fits.fits` / `lick_hamspec_composite_arc.fits` — there is
   clear intent that lick_hamspec is the "real" name. Worth clarifying whether
   `lick_hamspec` should *replace* or *coexist with* `shane_hamspec` (Q&A 2).
2. **Frame typing** must separate wide-plate pixflats from narrow-plate
   blazeflats (§4). Also: arcs taken with `Open` *and* `BG12` filters; science
   `Off` lamp.
3. **No bias/dark** in this set → rely on overscan; set `use_biasimage=False`
   unless biases are added to the dev-suite setup later.
4. **Saturated arc lines** → enable saturated-pixel masking for wavecal; the
   short 3 s arc (d2084) is a useful unsaturated companion to the 20 s arc.
5. **Wavelength calibration** (the hard part, per the prompt) has supporting
   material in this folder: `Arc_01_fit.idl` (XIDL output), `Hamilton_ThArTi.pdf`
   (IRAF screenshots), and the `IRAF/` outputs (`ecTHAR*`, `THAR*.fits`,
   `*.wrsfbo*`). These are the next reports' subject.

---

## Figures

| File | Description |
|------|-------------|
| `01_fig1_overview.png` | ZScale overview of pixflat / arc / science (full frame). |
| `01_fig2_pixflat_zoom.png` | Pixflat zooms — smooth flood, no order structure. |
| `01_fig3_pixflat_cuts.png` | Pixflat row/column cuts — broad blaze hump. |
| `01_fig4_blazeflat.png` | Blazeflat full + center zoom — echelle orders clearly resolved. |
| `01_fig5_blazeflat_xd_cuts.png` | Cross-dispersion cuts — order count/spacing & two-band format. |

---

## Q&A (questions for JXP)

1. **`shane_hamspec` msgs→logging migration:** I made the minimal fix to unbreak
   package imports (it blocked *all* CLI scripts). OK to keep this in the PR for
   the lick_hamspec work, or should it be split into its own commit/PR?
2. **`shane_hamspec` vs `lick_hamspec`:** the Hamspec spectrograph physically
   lives at Lick (Shane 3 m). Should `lick_hamspec` *replace* `shane_hamspec`
   (rename + deprecate), or be a separate spectrograph class? The existing
   `get_echelle_angle_files()` already points at `lick_hamspec_*` files.
3. **Frame typing of flats:** confirm that pixflats are the wide-plate
   (`PLATENAM=800:6.0`) flood frames and blazeflats are the narrow
   (`640:5.0`) slit frames, and that we should type them via `PLATENAM`/`OBJECT`
   (the lamp alone can't distinguish them).
4. **Two illuminated bands** in the blazeflat (§4): is the fainter lower band a
   real part of the cross-dispersed format, or scattered light / vignetting we
   should mask?
5. **Standard vs science:** HD188209 (O9.5Iab) is the only on-sky frame — is the
   intent to use it as a telluric/standard, with science frames to be added to
   the dev-suite setup later?
