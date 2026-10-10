# Keck/OSIRIS in PypeIt: design document (DRAFT 1)

Date: 2026-10-10.  Authors: JXP and Claude (Fable 5.1).
Status: draft for discussion.  Decisions recorded here come from the Q&A
in `claude_prompts/keck_osiris.md` (questions 1-17); open points are
listed in Section 10 and mirrored as new Q&A items.

Companion report: `Reports/01_test_data_inventory.md` (file inventory,
header tables, measured geometry, rectification-matrix format).

---

## 1. Instrument and data summary

OSIRIS is the Keck I adaptive-optics near-infrared **lenslet** integral
field spectrograph (Larkin et al. 2006).  Facts that drive the design:

| Property | Value | Source |
|---|---|---|
| Spectrograph detector (since 2016-01-01) | Teledyne H2RG, 2048 x 2048, gain 2.15 e-/ADU, RN 22 e- (CDS) to 7.8 e- (MCDS-32), dark < 0.025 e-/s, 32 readout channels | specs.html, manual App. A |
| Plate scales | 0.020, 0.035, 0.050, 0.100 arcsec/lenslet (header `SSCALE`) | manual 2.3 |
| Filters | 5 broadband (Zbb, Jbb, Hbb, Kbb, Kcb) + 18 narrowband (Zn4, Jn1-4, Hn1-5, Kn1-5, Kc3-5); header `SFILTER` | scale_filter.html |
| Orders / dispersion | Z 6th, J 5th, H 4th, K 3rd; 0.141 / 0.169 / 0.212 / 0.282 nm/px raw | manual Table 2-3 |
| Resolution | R ~ 3800 (2800-4500 across the field); ~2 px FWHM | manual 2.4 |
| Lenslet grid | 51 columns x 66 rows partially illuminated; broadband 16x64 complete (19 columns stored = 1216 slices); narrowband up to 48x64 (3072 spectra) | manual 2.2, App. E |
| Detector layout | lenslet grid rotated 3.6 deg (tan = 1/16); spectra horizontal, 2 px FWHM, **2 rows apart**, staggered 29 px in dispersion; pattern repeats every 32 rows (16 spectra x 64 bands) | manual 2.1; measured in Report 01 |
| Raw FITS (post-2016) | ext0 `SCI` float32 in ADU/s/coadd; ext1 `VAR`; ext2 `DQ` byte (9 = good); ITIME in ms; no flats/arcs taken by observers | Report 01 |
| Calibrations supplied by Keck | rectification matrices per (filter, scale, epoch): `sYYMMDD_cNNN___infl_<filt>_<scale>.fits`, ~160 MB; global wavelength coefficients per epoch (`osiris_wave_coeffs_*.fits`, 66x51x4) | OsirisDRP `data/`, App. E |
| Observer calibrations | darks (`SFILTER='Drk'`), sky frames (`ISSKY=1`, unreliable), A0V tellurics | manual 4.3 |
| Wavelengths | vacuum; DRP solution has 0.3-0.6 A scale/epoch offsets vs OH lines | manual 4.3.2 |

The consequence that shapes everything: **individual lenslet spectra
cannot be isolated by slit tracing** (neighbours 2 px apart with a 2 px
PSF).  The Keck DRP solves a per-column linear inverse problem using a
measured PSF model (the rectification matrix).  PypeIt will adopt that
model.

## 2. Decisions (from Q&A 1-17)

| # | Decision |
|---|---|
| 1 | Port the DRP rectification (PSF deconvolution) to numpy; run PypeIt's **Fiber** pypeline on the rectified frame (one row per lenslet). |
| 2 | Test data: OsirisDRP Hn3/100 (2016) pair + matrix as the first setup; a KOA science-mode set (Kbb/Kn3 at 35/50 mas) to follow (Section 9). |
| 3 | Support the H2RG era only (MJD >= 57388); raise a clear error otherwise.  2015 frames are kept only as a header-era test. |
| 4 | `get_rawimage` returns DN (`SCI x ITIME[s] x COADDS`); BPM = (`DQ != 9`) OR Keck hot mask OR Keck dead mask. |
| 5 | Frame types: dark (`SFILTER=='Drk'`), science/sky/standard from `OBSTYPE`/`ISSKY`; pairing via `comb_id`/`bkg_id`; `configuration_keys = ['dispname','decker']`. |
| 6 | Wavelengths: Keck coefficients (epoch-matched) as the template, refined on OH sky lines where available; vacuum. |
| 7 | Rectification matrices supplied by the user through a directory parameter with automatic selection by filter/scale/date; test matrices live beside the raw data in the dev suite. |
| 8 | Names: `keck_osiris.py`, `KeckOSIRISSpectrograph`, `name='keck_osiris'`, `camera='OSIRIS'`. |
| 9 | Design document first; code in phases afterwards. |
| 10 | Reports in `pypeitdev/keck_osiris/Reports/`, log in the prompt file. |
| 11 | Kn5/035 matrix unobtainable; 2015 set is header-only; search KOA for data. |
| 12 | Header tolerance rules (ITIME ms/s, truncated headers, missing weather and `INSTRUME`). |
| 13 | VAR used only as an optional weight for the solver; PypeIt computes its own variance. |
| 14 | All on-sky frames typed `science` (plus `sky` when `ISSKY=1`); user sets `bkg_id`. |
| 15 | Sky: A-B subtraction **and** Davies (2007) OH-scaled sky subtraction. |
| 16 | Acceptance test: rectified spectra agree with the DRP reference cube. |
| 17 | This outline. |

## 3. Processing flow

```
raw SCI/VAR/DQ (ADU/s)                       Keck calibrations
      |                                      - rectification matrix (filter, scale, epoch)
      v                                      - wave coeffs (epoch), OH line list
[A] get_rawimage: -> DN, exptime, BPM        - static hot/dead masks (pypeit/data)
      |
      v
[B] RawImage.process: dark subtraction (optional), BPM, cosmic-ray mask
      |                 (no bias, no overscan, no flat)
      v
[C] A-B / A-sky frame subtraction (PypeIt bkg_id machinery, raw-frame level)
      |
      v
[D] Rectification: (2048 x 2048 detector) -> (2048 spec x 1216 slices)
      |            per-column linear solve with the influence matrix;
      |            output image, ivar and mask; "detector" image for the
      |            rest of PypeIt is this rectified frame
      v
[E] Fiber pypeline calibrations:
      slits  = block-slits defined geometrically (64 slices per block)
      tilts  = identity (spectra are already straight)
      wvcalib= Keck coefficients -> WaveCalib, optional OH refinement
      flats  = none (uniform fiber profile) or a Keck white-light scan
      |
      v
[F] FiberFindObjects: one SpecObj per lenslet; sky: residual OH scaling
    (Davies 2007) on the rectified 2D frame or on the 1D spectra
      |
      v
[G] FiberExtract: boxcar + optimal (profile = influence-matrix row sum
    collapsed to 1 row, i.e. essentially unit-width) -> spec1d
      |
      v
[H] Datacube: spec1d -> (lambda, y, x) cube via the lenslet (row, col)
    map and the SSCALE/ROTPOSN WCS; telluric/flux with PypeIt's IR tools
```

Where the rectification runs is the main architectural choice; see
Section 5.3.

## 4. Class design: `pypeit/spectrographs/keck_osiris.py`

Modelled on `keck_kcwi.py` for structure and `mmt_binospec.py` for the
Fiber hooks.  Logging via `from pypeit import log` and `PypeItError`.

### 4.1 Class attributes

```python
class KeckOSIRISSpectrograph(spectrograph.Spectrograph):
    ndet = 1
    name = 'keck_osiris'
    telescope = telescopes.KeckTelescopePar()
    camera = 'OSIRIS'
    header_name = 'OSIRIS'          # label only; not used for matching
    url = 'https://www2.keck.hawaii.edu/inst/osiris/'
    pypeline = 'Fiber'
    supported = False               # until a dev-suite reduce test passes
    comment = 'Lenslet IFU; H2RG era (>= 2016) only; requires Keck rectification matrices'
    # Geometry constants (manual App. E, verified on data)
    NSLICE_BB = 1216                # 19 columns x 64 rows stored in the matrix
    NROW_BAND = 32                  # detector rows per lenslet-column band
    MJD_H2RG = 57388.0              # 2016-01-01
```

### 4.2 Metadata (`init_meta`, `compound_meta`)

| PypeIt key | Source | Notes |
|---|---|---|
| `ra`, `dec` | `RA`, `DEC` (deg) | present in all frames |
| `target` | `TARGNAME` | `OBJECT` is often empty |
| `dispname` | `SFILTER` | configuration key |
| `decker` | `SSCALE` | configuration key; string `'0.100'`, compared as string |
| `binning` | default `'1,1'` | |
| `mjd` | `MJD-OBS` | |
| `exptime` | compound: `ITIME x COADDS`, ITIME/1000 if `ITIME > 2000` | TRUITIME never present |
| `airmass` | `AIRMASS` | |
| `dithoff`/`ra_off`/`dec_off` | `RAOFF`, `DECOFF` | informational only |
| `posang` | compound: `ROTPOSN - INSTANGL` (+ `ROTREFAN`) | `PA_SPEC` may be absent; Keck I cube is flipped (DRP HISTORY) |
| `idname` | `OBSTYPE` (first letter) | `'a'`stro, `'s'`tar, `'c'`alib |
| `issky` | `ISSKY`, `required=False`, default 0 | |
| `instrument` | compound: `CURRINST`, fallback `INSTRUME` | 2016 files lack INSTRUME |
| `sampmode`, `numreads` | `SAMPMODE`, `NUMREADS` | sets read noise |
| `pressure`, `temperature`, `humidity` | `WXPRESS`, `WXOUTTMP`, `WXOUTHUM` with KCWI-style defaults + warning | required by the IFU subheader |
| `parangle` | `PARANG` (deg -> rad) | |
| `tmatemp` | `DTMP7`, `required=False` | grating-temperature wavelength correction |
| `slitwid` | compound: `SSCALE` in degrees | IFU subheader; unused otherwise |

`configuration_keys()` -> `['dispname', 'decker']`.
`raw_header_cards()` -> `['SFILTER', 'SSCALE', 'SAMPMODE', 'NUMREADS',
'ITIME', 'COADDS', 'ROTPOSN', 'INSTANGL']`.
`pypeit_file_keys()` adds `['issky', 'idname', 'calib', 'comb_id', 'bkg_id']`.

### 4.3 Frame typing (`check_frame_type`)

| Type | Rule |
|---|---|
| `dark` | `dispname == 'Drk'` |
| `science` | `idname == 'a'` and `dispname != 'Drk'` (all on-sky frames; pairing via `bkg_id`) |
| `sky` | as science **and** `issky == 1` (helps the user fill `bkg_id`) |
| `standard` | `idname == 's'` |
| `trace` | same as science (needed so `IFUCalibrations.get_slits` has a frame to build block-slits from; see 5.4) |
| `arc`, `tilt` | same as science (sky frame if available) - the OH lines are the wavelength reference for the optional refinement |
| `bias`, `pixelflat`, `illumflat`, `align`, `pinhole`, `scattlight` | never |

An era check in `get_rawimage` raises `PypeItError` for `mjd < 57388`.

### 4.4 Detector (`get_detector_par`)

```python
dict(det=1, binning='1,1', dataext=0, specaxis=1, specflip=False,
     spatflip=False, platescale=<SSCALE from header; 0.05 default>,
     darkcurr=90.0,            # e-/pix/hr (0.025 e-/s)
     saturation=65535., nonlinear=0.80,   # 80 % full well, per DRP practice
     mincounts=-1e10, numamplifiers=1,
     gain=np.atleast_1d(2.15),
     ronoise=np.atleast_1d(rn),  # from SAMPMODE/NUMREADS: CDS 22; MCDS n:
                                 # 12.2 (4), 8.1 (16), 7.8 (32); UTR ~ 8
     datasec=np.atleast_1d('[1:2048,1:2048]'), oscansec=None)
```

`specaxis = 1` because the raw frame has wavelength along columns (axis
1 in numpy); PypeIt will transpose to (spec, spat).  After rectification
the "spatial" axis is the slice index.

### 4.5 Raw image reading (`get_rawimage`)

1. Open the file; read `SCI`, `VAR` (optional), `DQ` (optional).
2. Era check; compute `exptime`; multiply `SCI` by `exptime` to get DN.
3. Build `rawdatasec_img = 1` everywhere, `oscansec_img = 0`.
4. Stash `VAR` and `DQ` on the returned `hdu` (or as attributes) so the
   rectification step and `bpm()` can use them without re-reading.

### 4.6 Bad-pixel mask (`bpm`)

`DQ != 9` OR `badpixelmask20170902_sigma50` OR `bpm_dead_4pixbufgood`
(both shipped as `pypeit/data/static_calibs/keck_osiris/*.fits.gz`,
~100 kB).  The mask is used by the solver (pixels excluded from the
residual) and by PypeIt.

### 4.7 Fiber-pypeline hooks

| Method | Returns | OSIRIS implementation |
|---|---|---|
| `get_fiber_blocks(det)` | list of dicts (`block_id`, `nfibers`, `type`, `fiber_positions`, `fiber_names`, `fiber_ids`, `min_pix`, `max_pix`) | 19 blocks (broadband) of 64 slices each; positions = slice index + 0.5 within the rectified frame; names `L<col>_<row>`; type `'science'` for all (no sky fibers) |
| `get_block_slit_edges(traceimg, det)` | `(left, right)` arrays `(nspec, nblocks)` | constant edges at `64*k - 0.5` and `64*(k+1) - 0.5`; no cross-correlation needed because the rectified geometry is exact |
| `identify_fibers_in_block(det, block_idx, detected_positions)` | ids/names/types aligned with detections | nearest-integer match to slice index |
| `get_fiber_metadata(det, slit_spat_ids, slit_centers)` | ids/names/types per slit | same mapping, applied after extraction |
| `load_sky_layout()` | `(x_arcsec, y_arcsec)` for every lenslet | from the slice-to-lenslet map (`kbb_2016_slice_to_lenslet.txt` generalised to all filters via the DRP `assembcube` rules) times `SSCALE` |
| `get_science_fiber_layout_indices(det, fiber_ids, fiber_types)` | indices into the layout | identity on lenslet id |
| `ifu_sky_wcs(raw_hdr, scale)` | `(SkyCoord, CD)` | `RA`/`DEC`, position angle from `ROTPOSN`/`INSTANGL`, Keck I flip |
| `get_wcs`, `get_datacube_bins` | 3-D WCS and bins | as KCWI but with lenslet grid (66 x 51 or the filter's subset) |

Narrowband filters: the DRP packs up to three lenslet-column blocks head
to tail along the same 2048-px slice.  Phase 1 supports broadband only;
narrowband needs the split rule from `assembcube_000.pro` (column
offsets 0, +16, +32 and the `sp > 831` remainder) applied when building
`get_fiber_blocks` and the layout.  Hn3/100 (our only test set) is
narrowband, so the split must be implemented early (Section 10).

### 4.8 Parameters (`default_pypeit_par`)

```python
par['rdx']['...']                                   # nothing special
par['calibrations']['darkframe']['exprng'] = [None, None]
par['calibrations']['biasframe']   -> never used
par['calibrations']['slitedges']   -> bypassed by get_block_slit_edges
par['calibrations']['tilts']       -> see 5.4 (identity tilts)
par['calibrations']['wavelengths']['method'] = 'full_template'   # phase 2
par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES'] (or an OSIRIS OH list)
par['calibrations']['wavelengths']['reference'] = 'sky'
par['calibrations']['wavelengths']['fwhm'] = 2.0
par['calibrations']['flatfield']['method'] = 'skip'
par['calibrations']['osiris_rectmat_dir'] = None   # NEW: directory of rect matrices
par['calibrations']['osiris_rect_niter'] = 40      # NEW
par['calibrations']['osiris_rect_relax'] = 1.0     # NEW (DRP 'relaxation')
par['scienceframe']['process'] : use_biasimage=False, use_overscan=False,
      use_darkimage=True (optional), use_pixelflat=False, use_illumflat=False,
      mask_cr=False (the DRP recommends against LACosmic on >=2016 data;
      cosmic rays are removed by frame combination)
par['reduce']['findobj']  -> fiber mode (no peak finding)
par['reduce']['skysub']['joint_fit'] = False
par['reduce']['skysub']['osiris_oh_scaling'] = True   # NEW (Davies 2007), phase 2
par['reduce']['extraction']['boxcar_radius'] = 0.5 * platescale (one slice)
par['reduce']['cube']['slit_spec'] = False
par['flexure']['spec_method'] = 'skip'
par['sensfunc']['algorithm'] = 'IR'
par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'
```

New parameters go through the `add-parameter` skill (names are
proposals; see Q&A).

## 5. The rectification step

### 5.1 Model

For detector column `x` (one wavelength sample), raw pixel rows `j`,
slices `s` with bottom row `b_s` (matrix ext0[:,0]) and influence
`B[s, l, x]` (ext2, `l = 0..15`):

    F[j, x] = sum_s  B[s, j - b_s, x] * c[s, x]  + noise,     DQ[j, x] == 9

Each column is an independent sparse linear system: 2048 equations, 1216
unknowns (1017 effective), at most 16 nonzeros per unknown.  The system
is well posed (each slice's PSF overlaps only its +/- 2-px neighbours
within the same 32-row band) but ill-conditioned where two PSFs overlap
strongly, which is why the DRP damps and smooths.

### 5.2 Solver options

| Option | Pros | Cons |
|---|---|---|
| **(a) Faithful port of `spatrectif_000.c`** (damped blame iteration: 15 sharpened iterations with +/-64-slice coupling, then plain back-projection; 5-tap spectral smoothing at iterations 12-17; divide by 1.28; 40 iterations) | reproduces the DRP exactly, so the reference cube is a direct test; vectorises cleanly over all 2048 columns at once (arrays of shape `(nslice, 16, ncol)`) | hand-tuned schedule with no noise model; the 1.28 factor and smoothing are empirical |
| (b) Weighted linear least squares per column (`scipy.sparse.linalg.lsqr`/`lsmr`, weights from VAR or PypeIt variance, non-negativity optional) | principled; propagates ivar; converges in ~20 iterations | will not match the DRP cube bit-for-bit; needs regularisation where PSFs overlap |
| (c) Dense per-band solve: within each 32-row band only ~16-19 slices contribute, so each column splits into 64 independent 32 x 19 systems solvable with `numpy.linalg.lstsq` in one batched call | fastest and exact; natural ivar propagation | assumes no cross-band leakage (true by construction of the matrix, `BASESIZE = 15 < 32`) |

**Recommendation:** implement (a) first as the validation reference,
then (c) as the production solver, keeping both behind a parameter.
Both are pure numpy; the per-column work is ~1216 x 16 multiply-adds per
iteration, i.e. ~40 x 2048 x 20k = 1.6 Gflop for (a), seconds in numpy.

### 5.3 Where it runs

Two candidate hooks:

1. **Inside `get_rawimage`** (spectrograph-level).  `RawImage` then sees a
   `(2048 spec, 1216 spat)` image from the start.  Simple, but dark
   subtraction, BPM and A-B happen *after* rectification, and the mask
   cannot be applied to the raw residual.  Rejected.
2. **A new `RawImage` processing step `rectify`** placed after `trim`,
   `orient`, dark subtraction and BPM, and before cosmic-ray masking,
   invoked only when `spectrograph.has_rectification` is True.  The
   step calls `spectrograph.rectify(image, ivar, bpm, hdr)` which loads
   the matrix (cached per filter/scale/epoch) and returns the rectified
   image, ivar and mask; `datasec_img`/`rn2img`/`proc_var` are
   re-shaped accordingly.  A-B subtraction in PypeIt operates on
   processed `PypeItImage` objects, so it works on rectified frames;
   because the solver is linear this is equivalent to rectifying the
   difference (up to the DQ mask).

**Recommendation:** option 2.  It touches `rawimage.py` once (a generic
"spectrograph-provided transform" step, which other lenslet or
image-slicer instruments could reuse) and keeps all OSIRIS knowledge in
`keck_osiris.py` plus a core module `pypeit/core/lenslet_rectify.py`
holding the solvers.

### 5.4 Calibrations on the rectified frame

- **Slits**: `IFUCalibrations.get_slits` uses `get_block_slit_edges`
  when the spectrograph defines it.  It still requires a `trace` frame,
  hence science frames are also typed `trace`.  Edges are the fixed
  block boundaries.
- **Tilts**: spectra are exactly horizontal after rectification, so the
  tilt model is the identity.  `Calibrations.get_tilts` needs an
  `mstilt` frame and a wavelength solution; the cleanest path is a
  `WaveTilts` object built directly by the spectrograph (a new
  `Spectrograph.get_fixed_tilts(slits)` hook, or reuse of the existing
  `tilts` par with `['calibrations']['tilts']['method'] = 'identity'`
  if such an option is added).  Phase-1 fallback: let PypeIt fit tilts
  on OH lines of the rectified sky frame (they will come out ~0).
- **Wavelengths**: see Section 6.
- **Flats**: none.  `FiberExtract` falls back to a uniform profile when
  no flat image exists, and each lenslet occupies exactly one rectified
  row, so boxcar with radius 0.5 row is the natural extraction; optimal
  extraction adds nothing in phase 1.
- **Scattered light, pattern noise, non-linearity**: off.

## 6. Wavelength calibration

### 6.1 Keck coefficients as the template

`osiris_wave_coeffs_<epoch>.fits` is a `(66 row, 51 col, 4)` cube of
polynomial coefficients giving **pixel as a function of wavelength**:

    x(lambda) = poly(lambda', coeffs[row, col]),
    lambda' = lambda * order/3 / frac_expan(DTMP7) - 2200 nm

(3rd-order K-band reference; `frac_expan` from the aluminium CTE table in
`assembcube_000.pro`; `order` from the filter's minimum wavelength).  For
each slice we therefore have an exact, invertible wavelength solution
along the 2048-px rectified row.  Implementation:

1. Ship the 11 epoch files (small, ~55 kB each) and the CTE table in
   `pypeit/data/spectrographs/keck_osiris/`; pick the epoch from
   `MJD-OBS` like `calibrations.xml` does.
2. Build a `WaveCalib` object directly: for every slit (block) a
   `WaveFit` whose `pypeitfit` is a polynomial fit of lambda(x) derived
   by inverting the Keck relation on a fine grid (the inverse is smooth
   and monotonic), with `wave_soln` on the 2048 pixels.  Per-lenslet
   offsets within a block are small but nonzero; since the Fiber
   pypeline evaluates wavelengths per fiber from the block solution plus
   tilts, the per-lenslet term is encoded either as a per-fiber shift
   (there is already `_refine_fiber_wavelengths`/`_solve_fiber_shift`
   machinery in `FiberFindObjects`) or by making every lenslet its own
   "slit" (1216 slits; heavier but exact).  **Proposal:** one slit per
   block, per-lenslet shifts from the Keck coefficients, consistent with
   how Binospec handles fibers.
3. The science header's `DTMP7` applies the grating-temperature scaling;
   when absent (truncated header) use the epoch default (no scaling) and
   warn.

### 6.2 OH refinement (phase 2)

With `reference = 'sky'` and an OH line list (`ohlines_kbb_all.dat` for
K; PypeIt's `OH_NIRES`/`OH_FIRE_Echelle` lists for J/H), run the
existing fiber shift solver on the A frame (before sky subtraction) to
measure a zero-point (and optionally linear) correction per lenslet; the
DRP documents 0.3-0.6 A offsets depending on scale and epoch.  Where
fewer than 3 lines are available (long-K, narrowband) keep the Keck
solution.  Report the per-lenslet shift statistics in QA.

### 6.3 Vacuum

Keck solutions are vacuum.  Set the PypeIt wavelength frame accordingly
(no air-to-vacuum conversion) and record it in the spec1d header.

## 7. Calibration-file handling

| File | Size | Location | Selection |
|---|---|---|---|
| Rectification matrix | ~160 MB per (filter, scale, epoch) | user directory via `osiris_rectmat_dir` (new parameter); dev-suite copies in `RAW_DATA/keck_osiris/<setup>/calib/` | parse `sYYMMDD_cNNN___infl_<filt>_<scale>.fits`; choose the latest matrix with date <= MJD-OBS matching filter and scale; error listing what was found otherwise |
| Wavelength coefficients (11 epochs) + CTE table | 0.6 MB total | `pypeit/data/spectrographs/keck_osiris/` | epoch table copied from `calibrations.xml` |
| Hot/dead pixel masks | ~100 kB gz | `pypeit/data/static_calibs/keck_osiris/` | always |
| Slice-to-lenslet maps | text | `pypeit/data/spectrographs/keck_osiris/` | per filter family (bb / nb) |
| OH line lists | text | `pypeit/data/arc_lines/lists/` (new `OH_OSIRIS_K.dat` from `ohlines_kbb_all.dat`) | by band |

The matrix loader caches the parsed `(hilo, effective, influence)`
tuple in memory for the run (one matrix per configuration).  A later
convenience could add a PypeIt-cache download of the common-mode
matrices if Keck agrees to host them with stable URLs.

## 8. Tests and validation

1. **Unit tests** (`pypeit/tests/test_keck_osiris.py`, no large files):
   metadata parsing from synthetic headers (both eras, truncated
   header), frame typing, exptime unit logic, matrix filename parsing,
   wavelength inversion round-trip on the shipped coefficients, and a
   solver test on a synthetic 64-row band with known input.
2. **Rectification acceptance test** (dev-suite "vet" test): rectify
   `s160321_a002010 - s160321_a002011` with
   `s160318_c005___infl_Hn3_100.fits`, resample each slice onto the DRP
   wavelength grid and compare with `s160321_a002010_Hn3_100_ref.fits`
   (spaxel-by-spaxel, after the narrowband split).  Target: median
   absolute fractional difference < 2 % on spaxels with S/N > 10;
   quasar spaxel flux within 5 %.
3. **Dev-suite reduce test** `keck_osiris/Hn3_100` (pypeit file with the
   pair, `bkg_id` set), then a Kbb/050 science set from KOA (Section 9)
   with its A0V telluric for the full flow including cube and telluric
   correction.
4. **Load-images unit test** entry in `unit_tests/test_load_images.py`.
5. **Header-era test**: the 2015 Kn5 files are read by `pypeit_setup`
   and produce a clear "pre-2016 data not supported" error.

## 9. Data search (KOA)

A TAP query against `koa_osiris` (public spectrograph frames, ITIME >=
60 s, since 2016) shows the most used science modes are Kn3/020, Kbb/020,
Kn3/050, Kn3/035, Kbb/050 and Kbb/035 (3.5k-11k frames each).  Nights
with many Kbb or Kn3 frames at 35/50 mas on one target, with sky frames
flagged, include:

| Night | Target | Mode | Frames | Total itime (s) | Sky flagged | Note |
|---|---|---|---|---|---|---|
| 2024-02-03 | NGC 3245 | Kbb/050 | 23 | 11580 | 8 | plus HD101060 (A0V) same night: a complete science + telluric set on the current detector |
| 2024-07-16 | PGC 70520 | Kbb/050 | 22 | 10920 | 8 | plus LHS2803A, 36 Her same night |
| 2023-07-25/26 | tt001 / tt011 / tt025 | Kbb/050 | 23-31 | ~11000 each | 9-13 | Keck "tt" target names; galaxies/nuclei |
| 2017-06-18 | PG 1411+442 | Kbb/035 | 20 | 7620 | 1 | quasar; 35 mas |
| 2021-07-28 | HR 8799 | Kbb/035 | 191 | 9720 | 0 | bright planet host, many short frames |
| 2019-06-04 | ZTF19aalypsi etc. | Kn3/050 | 26-108 | 100-430 | 0 | transients, 4 s frames |

**Recommendation:** adopt the 2024-02-03 NGC 3245 + HD101060 set
(Kbb/050) as the science-mode dev-suite setup, together with the darks
from that night (`SFILTER='Drk'`), and request the matching Kbb/050
rectification matrix (2023-2024 epoch) from Keck.  Download via PyKOA or
the KOA TAP `filehand` column once JXP confirms.

## 10. Open issues and phases

### Open issues (new Q&A 18-22)

1. Narrowband support is needed for the only test set we have (Hn3);
   implement the DRP's three-block split in phase 1 or obtain a
   broadband set first?
2. Where to put the rectification: new generic `RawImage` step (proposed)
   vs a spectrograph-only hook.
3. Tilts: add an "identity tilts" option to PypeIt vs fit OH lines.
4. Wavelength bookkeeping: one slit per 64-lenslet block with per-lenslet
   shifts (proposed) vs one slit per lenslet.
5. KOA dataset choice (Section 9) and who requests the matrix from Keck.

### Phases

| Phase | Deliverable | Depends on |
|---|---|---|
| 1 | `keck_osiris.py` metadata, typing, detector, `get_rawimage`, BPM, era check; static files; unit tests; `pypeit_setup` runs on both test sets | - |
| 2 | `pypeit/core/lenslet_rectify.py` (solver a and c), matrix loader/selection, `RawImage.rectify` step; vet test against the DRP cube | Q18, Q19 |
| 3 | Fiber hooks (blocks, ids, layout, WCS), identity tilts, Keck-coefficient `WaveCalib`; full run to spec1d on Hn3/100 | Q20, Q21 |
| 4 | Datacube, A-B and Davies OH-scaled sky subtraction, OH wavelength refinement; Kbb/050 KOA set; dev-suite reduce test; docs and changelog | Q22, matrix from Keck |
| 5 | Narrowband generalisation, telluric/flux calibration, `supported = True` | 4 |
