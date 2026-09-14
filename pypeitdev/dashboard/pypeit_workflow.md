# The PypeIt Workflow

**Version:** 0.2
**Date:** 2026-06-12
**Author:** JXP and Claude

**Changelog**
- 0.1 (2026-06-12): Initial draft from the PypeIt docs.
- 0.2 (2026-06-12): Added Section 6, a concrete worked example from a real
  `run_pypeit` reduction of the Shane Kast blue 600/4310 dev-suite dataset
  (observed file listing, runtime, and QA observations).

---

## Purpose

This document records our working understanding of the end-to-end PypeIt
data-reduction workflow.  It is intended to guide the design and development of
the **PypeIt Dashboard** by capturing (a) the sequence of user actions, (b) the
scripts and files involved at each stage, (c) the data products generated, and
(d) the inspection / quality-assessment (QA) touch-points where a user most
needs visibility.  The dashboard's job is to make this workflow legible and
navigable: to show where a reduction is, what it has produced, and where it may
have gone wrong.

This is a living document.  It draws primarily on `PypeIt/doc/running.rst`,
`PypeIt/doc/setup.rst`, the instrument tutorials in `PypeIt/doc/tutorials/`
(especially the Shane Kast HOWTO), `PypeIt/doc/outputs.rst`, `PypeIt/doc/qa.rst`,
and `PypeIt/doc/scripts.rst`.

---

## 1. High-level picture

PypeIt converts **raw spectroscopic frames** from a supported spectrograph into
**calibrated 2D spectral images** and **1D extracted spectra**, optionally
followed by **flux calibration and coaddition**.  A core design principle is the
separation of *spectrograph-specific* code from *generally applicable* slit-based
reduction code; instrument quirks live in `Spectrograph` subclasses, while the
reduction algorithms are shared.

Three primary reduction paths exist, selected automatically per spectrograph:

| Path        | Used for                                   | Key differences |
|-------------|--------------------------------------------|-----------------|
| `MultiSlit` | Long-slit and multi-object slit (MOS)      | Object finding/extraction done independently per slit |
| `Echelle`   | Cross-dispersed echelle spectrographs      | Object traces matched across orders |
| `IFU`       | Slicer-based integral field units          | Different calibration order; different sky-sub & flat-fielding |

The typical user journey has **four phases**:

1. **Setup** — organize raw data and generate a PypeIt reduction file.
2. **Reduction** — run the core pipeline (`run_pypeit`).
3. **Inspection / QA** — review calibrations and reduced spectra.
4. **Further processing** — fluxing, coaddition, telluric, collation.

The dashboard will be focused on the second and third phases of the PypeIt workflow, i.e. the reduction of science frames and the inspection of the calibrations and reduced spectra.

---

## 2. Phase 1 — Setup

The single most important artifact in any PypeIt run is the **PypeIt file**
(`*.pypeit`); a correct setup is the precondition for a successful reduction.

This section reviews that effort, although an existing GUI tool is 
already available to do this.

### 2.1 Organize raw data
- All raw frames for a run are placed in a single folder, organized in the
  dev-suite as `RAW_DATA/<instrument>/<setup>/`.
- Only one copy (compressed *or* uncompressed) of each frame should be present.

### 2.2 `pypeit_obslog` (optional)
- Lists the raw files and the metadata PypeIt pulls from their headers.
- `pypeit_obslog <spectrograph> -k` prints the metadata-key → header-card map.
- Useful as a pre-flight check; not strictly required.

### 2.3 Instrument configurations ("setups")
- PypeIt groups frames into unique **instrument configurations** using a defined
  subset of header metadata (the `configuration_keys` for that spectrograph).
- Configuration-defining metadata typically includes: `dispname` (disperser),
  `dispangle`/`cenwave` (central wavelength), `decker` (slit mask), `dichroic`,
  and `binning`; the full set is instrument-specific.
- Each unique configuration is labeled with a capital letter: **A**, **B**, ...
- **Constraint:** each `.pypeit` file should contain data from exactly one
  configuration.

### 2.4 `pypeit_setup`
The automated tool that types frames and groups them by configuration.

- **First execution (no `-c`):**
  `pypeit_setup -r <raw_path>/<prefix> -s <spectrograph>`
  - `-s`: spectrograph (camera-specific, e.g. `keck_lris_blue`).
  - `-r`: full path + filename root prefix of the raw FITS files.
  - Produces, in `setup_files/`:
    - `<spectrograph>.sorted` — configurations and the frames in each (review
      only; **edits here are ignored** by PypeIt).
    - `<spectrograph>.obslog` — the default obslog listing.
- **Second execution (with `-c`):**
  `pypeit_setup -r <raw_path>/<prefix> -s <spectrograph> -c A`
  - `-c A` / `-c A,C` / `-c all` selects which configuration(s) to write.
  - Creates a sub-folder per configuration (e.g. `keck_lris_blue_A/`) containing
    the `.pypeit` file (`keck_lris_blue_A.pypeit`).
- **Useful options:**
  - `-b`: add `calib`, `comb_id`, `bkg_id` columns (calibration grouping, frame
    combination, A-B background sequences). Added by default for most near-IR
    spectrographs.
  - `-m`: add a column for manual extraction input.

### 2.5 The PypeIt file
The `.pypeit` file has the following structure and is the primary thing a user
edits by hand:

- **Parameter block** — user-level `PypeItPar` overrides (reduction parameters).
- **Setup block** — the instrument-configuration metadata (between `###...` and
  `#--...` lines). This (not the raw FITS headers) is what fixes the
  configuration at run time.
- **Data block** — the table of frames, their `frametype`, and relevant
  metadata. Optional columns: `calib`, `comb_id`, `bkg_id`, manual.

**Common manual edits before running:**
- Fix mis-assigned `frametype`s (e.g. standard stars are often mistyped as
  science; arcs/flats may need correction, especially for Kast).
- Remove or fix any entries with `None` frametype — **`run_pypeit` will crash**
  otherwise.
- Add/remove frames; set science↔calibration associations (calibration groups).
- Set frame-combination groups (2D combine) and background/AB-offset groups.

**Dashboard relevance:** parsing and presenting the `.pypeit` file (frame table,
frametypes, configuration, parameter overrides) is a natural dashboard feature.
Flagging `None` frametypes and missing calibrations pre-run would be high value.

---

## 3. Phase 2 — Reduction (`run_pypeit`)

### 3.1 Invocation
From within the configuration sub-folder:
```bash
cd keck_lris_blue_A
run_pypeit keck_lris_blue_A.pypeit -o
```

Key options:
- `-o` / `--overwrite`: overwrite existing outputs (recommended in most cases;
  does **not** rebuild existing processed calibrations).
- `-m` / `--do_not_use_calibs`: ignore calibrations on disk (slow; rarely
  wanted).
- `-c`: reduce calibrations only (required if the `.pypeit` file has only
  calibration frames).
- `-s`: developer debugging mode (floods the screen with plots).
- Re-reducing: prefer deleting files in `Calibrations/` over using `-m`.

### 3.2 The general workflow

1. **Initialization** — parse the `.pypeit` file, consolidate with header
   metadata (pypeit-file metadata wins), set up the output directory structure.
2. **Reduce standard-star frames** (if any) — their sole role in `run_pypeit` is
   to provide an object-tracing crutch for faint science objects. (Fluxing
   itself happens later, in separate scripts.)
3. **Reduce science frames** — the main reduction; associated standard-star
   traces are used as a tracing crutch when available.

### 3.3 Reduction loop structure
Both reduction steps differ only in the tracing crutch (slit edges for
standards; standard-star trace for science when available). Beware
differential atmospheric refraction whenever a tracing crutch is used.

Nested loops, outermost → innermost:
1. calibration groups,
2. frame-combination groups (or individual frames),
3. detectors or detector **mosaics**.

For each detector / mosaic:
- All frames undergo **image processing** (`image_proc`) first.
- Process the relevant **calibration** frames.
- First-pass **global sky model** estimate + **object finding**.
- (Instruments with slit-mask metadata, e.g. DEIMOS) match detected objects to
  expected slit-mask positions; flag undetected and serendipitous objects.
- A final loop over detectors/mosaics performs spectral **extraction**.

### 3.4 Per-path alterations
- `Echelle`: object traces matched across orders (vs. independent for MultiSlit).
- `IFU`: calibrations done in a slightly different order; sky-subtraction and
  flat-fielding handled differently.

### 3.5 Calibration products (written to `Calibrations/`)
Generated and reused across runs. The principal calibration frames and their
inspection tools:

| Calibration | File (example)              | Inspect with |
|-------------|-----------------------------|--------------|
| Bias        | `Bias_A_0_DET01.fits`       | `ginga` |
| Arc         | `Arc_A_0_DET01.fits`        | `ginga` |
| Slit edges  | `Edges_A_0_DET01.fits.gz`   | `pypeit_chk_edges` |
| Wave (1D)   | `WaveCalib_A_*.fits`        | `pypeit_chk_wavecalib`, `pypeit_show_wvcalib` |
| Tilts (2D)  | `Tilts_A_0_DET01.fits.gz`   | `pypeit_chk_tilts` |
| Flat        | `Flat_A_0_DET01.fits`       | `pypeit_chk_flats`, `pypeit_show_pixflat` |
| Alignments  | (IFU)                       | `pypeit_chk_alignments` |

Most calibration and output files subclass `pypeit.datamodel.DataContainer`,
enforcing a strict on-disk datamodel and a common IO interface.

---

## 4. Phase 3 — Outputs and inspection

### 4.1 Directory structure
Run from `${RDXDIR}/<PYP_SPEC>_<SETUP>/`, `run_pypeit` produces:

- `Calibrations/` — all calibration frames (see above).
- `Science/` — reduced science and standard frames: `spec2d_*` and `spec1d_*`.
- `QA/` — quality-assessment output: `PNGs/` plus per-setup HTML
  (e.g. `MF_A.html`; the PNGs are more reliable than the HTML).

### 4.2 Primary science products
- **`spec2d_*.fits`** — calibrated 2D spectral image (sky-subtracted), one per
  science/standard frame. Inspect with `pypeit_show_2dspec` (ginga).
- **`spec1d_*.fits`** — 1D extracted spectra. Inspect with `pypeit_show_1dspec`.
- **Bitmasks** — used throughout to mark untrustworthy pixels and record *why*
  pixels were flagged (`out_masks`).
- Processing **history** is logged in FITS headers.

### 4.3 QA figures (in `QA/PNGs/`)
Fixed-format PNGs generated during the run, by setup / detector / (sometimes)
slit. The QA the user most needs to review:

**Calibration QA**
- *Wavelength 1D fit* (`Arc_1dfit_*`): arc lines IDed in green spanning the full
  range; target RMS < 0.1 px; random residuals.
- *Wavelength tilts 2D* (`Arc_tilts_2d_*`): arc-line centroids vs. spatial
  position; few rejections; RMS < 0.1.
- *Echelle order prediction* (`Edges_*_orders_qa`): measured order widths/gaps
  vs. predicted polynomial; identifies missing/outlier orders.

**Exposure QA** (per science exposure)
- *Spatial flexure* (`*_spat_flex_corr`): slit-edge shift correction.
- *Spectral flexure* (`flex_corr_*`, `flex_sky_*`): correlation-lag fit and sky
  lines vs. archived sky spectrum.
- *Object finding* (`*_obj_prof.png`): S/N-collapsed spatial profile; detected
  objects vs. SNR threshold — key for confirming detections / threshold tuning.
- *Object tracing* (`*_obj_trace.png`): spatial peak vs. spectral pixel — key
  for catching traces that have "gone off the rails" (e.g. atmospheric
  dispersion at high airmass).

### 4.4 Inspection scripts (the `pypeit_chk_*` / `pypeit_show_*` family)
These are the primary CLI tools a user runs to vet a reduction. The dashboard
should ideally either wrap or replicate their views:

- Calibrations: `pypeit_chk_edges`, `pypeit_chk_tilts`, `pypeit_chk_flats`,
  `pypeit_chk_wavecalib`, `pypeit_show_wvcalib`, `pypeit_chk_alignments`,
  `pypeit_show_pixflat`.
- Pre-run sanity: `pypeit_chk_for_calibs`, `pypeit_chk_calibs`.
- Spectra: `pypeit_show_2dspec`, `pypeit_show_1dspec`, `pypeit_parse_slits`.
- Flexure: `pypeit_chk_flexure`.
- General FITS view: `pypeit_view_fits`.

Many of these launch a `ginga` RC viewer; spec1d uses an interactive GUI.

---

## 5. Phase 4 — Further processing

Post-reduction tools (run after `run_pypeit`), each driven by its own input file
and producing its own data products:

| Step             | Script                   | Input file   | Notes |
|------------------|--------------------------|--------------|-------|
| Sensitivity func | `pypeit_sensfunc`        | `.sens`      | From a reduced standard `spec1d`; outputs sensfunc FITS + QA. |
| Flux calibration | `pypeit_flux_calib`      | `.flux`      | **Edits the `spec1d` file in place.** `pypeit_flux_setup` helps build inputs. |
| 1D coadd         | `pypeit_coadd_1dspec`    | `.coadd1d`   | Produces a `OneSpec` file. |
| Collation        | `pypeit_collate_1d`      | —            | Groups/coadds by object; output identical to coadd1d. |
| 2D coadd         | `pypeit_coadd_2dspec`    | `.coadd2d`   | Combined `spec2d`. |
| 3D coadd (IFU)   | `pypeit_coadd_datacube`  | —            | Datacube. |
| Telluric         | `pypeit_tellfit`         | `tellfit`    | Telluric correction. |
| Quick-Look       | `pypeit_ql`              | —            | Fast, on-the-fly reduction. |

The dev suite mirrors these with per-setup input-file directories
(`coadd1d_files/`, `coadd2d_files/`, `fluxing_files/`, `sensfunc_files/`,
`tellfit_files/`, `flexure_files/`).

---

## 6. Worked example — Shane Kast blue 600/4310 (real run)

This section records a concrete, end-to-end `run_pypeit` reduction of one of the
simplest spectrographs in PypeIt, to ground the abstract workflow above in real
files and timings.

- **Spectrograph / path:** `shane_kast_blue` (long-slit → `MultiSlit`).
- **Dataset:** dev-suite `RAW_DATA/shane_kast_blue/600_4310_d55` (disperser
  `600/4310`, dichroic `d55`, single detector `DET01`, single slit `S0175`).
- **Reduction dir:** `REDUX_OUT/shane_kast_blue/600_4310_d55/shane_kast_blue_A`.
- **Command:** `run_pypeit shane_kast_blue_A.pypeit -o`.
- **Runtime:** ~1m 8s (single calibration group, exit code 0).

### 6.1 Input frames (from the `.pypeit` data block)
A single configuration (Setup A), one calibration group (`calib = 0`):

| frametype                   | count | files |
|-----------------------------|-------|-------|
| `arc,tilt`                  | 1     | b1 |
| `bias`                      | 10    | b14–b23 |
| `pixelflat,illumflat,trace` | 11    | b3–b13 |
| `science`                   | 2     | b27, b28 (J1217p3905) |
| `standard`                  | 1     | b24 (Feige 66) |

Note the single `arc,tilt` frame serves both wavelength calibration and tilt
tracing, and the dome flats triple as pixelflat/illumflat/trace. The standard
(Feige 66) is correctly typed here — unusually easy for Kast.

### 6.2 Output directory structure (what `run_pypeit` produced)
```
shane_kast_blue_A/
├── Calibrations/
│   ├── Arc_A_0_DET01.fits          # processed arc image
│   ├── Bias_A_0_DET01.fits         # master bias
│   ├── Edges_A_0_DET01.fits.gz     # slit-edge traces
│   ├── Flat_A_0_DET01.fits         # pixel/illum flat
│   ├── Slits_A_0_DET01.fits.gz     # finalized slit definitions
│   ├── Tiltimg_A_0_DET01.fits      # processed tilt image
│   ├── Tilts_A_0_DET01.fits        # 2D tilt solution
│   └── WaveCalib_A_0_DET01.fits    # wavelength solution
├── Science/
│   ├── spec1d_b24-Feige66_*.fits   # + .txt summary, one per frame
│   ├── spec1d_b27-J1217p3905_*.fits
│   ├── spec1d_b28-J1217p3905_*.fits
│   ├── spec2d_b24-Feige66_*.fits
│   ├── spec2d_b27-J1217p3905_*.fits
│   └── spec2d_b28-J1217p3905_*.fits
├── QA/
│   ├── MF_A.html                   # master QA index (per setup)
│   ├── <frame>.html                # per-exposure QA index (×3)
│   └── PNGs/                        # all QA figures (see 6.3)
├── shane_kast_blue_A.pypeit        # the input file
├── shane_kast_blue_A.calib         # calibration association record
├── shane_kast_blue_A_state.json    # reduction state / status
├── shane_kast_blue_A.status.log    # per-frame status table
├── shane_kast_blue_A.log           # full run log
└── shane_kast_blue_A_UTC_<date>.par  # the full, resolved parameter set
```

Observations relevant to the dashboard:
- Calibration file naming encodes `<type>_<setup>_<calibgrp>_<detector>`, e.g.
  `Arc_A_0_DET01` — directly parseable for a product browser/status grid.
- `spec1d`/`spec2d` names encode the raw frame, target, instrument, and UTC.
- Beyond the three doc'd directories, the run dropped several **state/log
  artifacts** (`*_state.json`, `*.status.log`, `*.calib`, `*.par`). The
  `*_state.json` / `*.status.log` pair is especially promising as a
  machine-readable source for the dashboard's run-status view (worth a closer
  look in the design phase — see Section 8).

### 6.3 QA figures produced (`QA/PNGs/`)
Calibration QA (per setup/detector/slit):
- `Arc_1dfit_A_0_DET01_S0175.png` — 1D wavelength fit.
- `Arc_FWHMfit_A_0_DET01_S0175.png` — arc-line FWHM fit.
- `Arc_tilts_2d_A_0_DET01_S0175.png`, `Arc_tilts_spat_*`, `Arc_tilts_spec_*` —
  tilt solution.
- `Spatillum_FineCorr_A_0_DET01_S0175_spatillum_finecorr.png` — illumination
  fine correction.

Per-exposure QA (one set per science/standard frame):
- `*_obj_prof.png` — object-finding S/N profile.
- `*_obj_trace.png` — object trace centroid fit.
- `*_spec_flex_corr.png`, `*_spec_flex_sky.png` — spectral flexure (both
  `global` and `local` variants written).

### 6.4 QA assessment (this run was healthy)
Inspected the key figures directly:
- **Wavelength 1D fit** (`Arc_1dfit`): RMS = **0.074 px** (target < 0.1 ✓),
  dispersion ΔÅ = 1.032 / pix; green-IDed arc lines span the full spectral
  range; residuals scatter randomly about zero.
- **Tilts 2D** (`Arc_tilts_2d`): RMS = **0.019** (fit orders spat=3, spec=5);
  only a small fraction of centroids rejected (red).
- **Object finding** (`obj_prof`, b27): a single clean object on slit 175, SNR
  peak ≈ **142** well above the `SNR_THRESH = 10` threshold; "1 Good Objects".
- **Object trace** (`obj_trace`, b27): smooth, monotonic centroid trace that the
  fit follows closely across all 2048 spectral pixels — no "off the rails"
  behavior.

These four figures, with their headline numbers (RMS, SNR, object count) and the
"what you hope to see" guidance, are exactly the at-a-glance signals the
dashboard's QA view should surface. Notably, the key quality metrics (wavelength
RMS, tilt RMS, SNR threshold, object count) are rendered into the PNG **titles /
legends** rather than being separately tabulated — so the dashboard likely needs
to pull these numbers from the data products (e.g. `WaveCalib`, `spec1d`
headers/datamodels) rather than scraping figures.

----

## 7. Dashboard implications (initial notes)

These are early observations to refine in the design document.

- **Phase-oriented layout.** Reduce → Inspect/QA → Further
  processing. At any moment the dashboard should answer: *which phase is this
  run in, and what has it produced?*
- **The `.pypeit` file is the hub.** Parse and display the parameter/setup/data
  blocks; validate (flag `None` frametypes, missing calibrations, suspicious
  frametype assignments) *before* the user runs.
- **Run monitoring.** `run_pypeit` is verbose and long-running; surfacing
  progress through the nested loop (calib group → comb group → detector/mosaic →
  extraction) and tailing the log would be valuable.  The code also produces a state file and a status log that can be used to monitor the progress.
- **Product browser.** Tree of `Calibrations/`, `Science/`, `QA/` with status
  (present/missing/stale) and one-click launch of the matching `pypeit_chk_*` /
  `pypeit_show_*` tool.
- **QA gallery.** The `QA/PNGs/` are fixed-format and the richest at-a-glance
  signal of reduction health (wavelength RMS, tilts, object find/trace). A
  thumbnail gallery with the "what you hope to see" guidance inline would lower
  the expertise barrier.
- **Path awareness.** MultiSlit vs. Echelle vs. IFU change which products and QA
  exist; the dashboard should adapt its views per spectrograph path.
- **Reuse.** Wherever possible, reuse existing PypeIt code (`DataContainer` IO,
  the metadata/`PypeItMetaData` machinery, the existing `pypeit_chk_*` logic)
  rather than reimplementing parsing.

---

## 8. Open questions / to investigate

- The exact mapping from `run_pypeit` runtime state to an observable progress
  signal (log parsing vs. an instrumentation hook). The Kast run (Section 6)
  revealed `*_state.json` and `*.status.log` artifacts that look like a
  machine-readable status source worth investigating before resorting to log
  scraping.
- Where the headline QA metrics (wavelength/tilt RMS, object SNR/count) live in
  the data products so the dashboard can read them directly rather than parsing
  the QA PNGs (which embed them only in titles/legends).
- Whether the dashboard should *launch* reductions or only *observe* them. -- Answer: both
- How calibration reuse / staleness is best detected from disk.
- IFU-specific products and QA (datacubes, alignments) need a closer look.
- Multi-configuration / multi-setup projects: how to present several `.pypeit`
  runs side by side. -- Answer:  we will focus on a single run at a time.

---

## References (PypeIt docs)

- `doc/running.rst` — core executable and workflow.
- `doc/setup.rst`, `doc/pypeit_file.rst` — setup and the PypeIt file.
- `doc/tutorials/` — instrument HOWTOs (Kast long-slit; DEIMOS MOS; MOSFIRE NIR;
  GNIRS echelle; coadd2d).
- `doc/outputs.rst`, `doc/out_spec1D.rst`, `doc/out_spec2D.rst`,
  `doc/out_masks.rst` — output products and datamodels.
- `doc/qa.rst` — QA figures.
- `doc/calibrations/` — per-calibration detail.
- `doc/scripts.rst` — CLI scripts, including the `pypeit_chk_*` / `pypeit_show_*`
  inspection family.
- `doc/fluxing.rst`, `doc/coadd1d.rst`, `doc/coadd2d.rst`, `doc/coadd3d.rst`,
  `doc/telluric.rst`, `doc/collate1d.rst`, `doc/quicklook.rst` — further
  processing.
