# Keck/OSIRIS in PypeIt: design document (DRAFT 2)

Date: 2026-10-10.  Authors: JXP and Claude (Fable 5.1).
Status: draft 2, accepted by JXP on 2026-10-10 with the answers to Q&A
18-22 folded in.  Decisions come from the Q&A in
`claude_prompts/keck_osiris/keck_osiris.md` (questions 1-25).  The staged build plan
is in `keck_osiris_implementation.md`.

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
| 4 | `get_rawimage` returns DN (`SCI x exptime[s]`, with `exptime = TRUITIME x COADDS` when `TRUITIME` is present, else `ITIME x COADDS`; M1 Q&A 2); BPM = (`DQ != 9`) OR Keck hot mask OR Keck dead mask. |
| 5 | Frame types: dark (`SFILTER=='Drk'`), science/sky from `OBSTYPE`/`ISSKY`, standard from `OBSTYPE=='s'` **or the target name** (M1 Q&A 1); pairing via `comb_id`/`bkg_id`; `configuration_keys = ['dispname','decker']`. |
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
| 18 | Implement the narrowband three-block split in phase 1 (needed for the Hn3 acceptance test); Kbb/050 is the first "supported" mode. |
| 19 | Rectification runs as a new generic `RawImage` step (`rectify`) calling a spectrograph-provided transform; solvers live in `pypeit/core/lenslet_rectify.py`. |
| 20 | Identity tilts via a new `Spectrograph.get_fixed_tilts` hook; **plus** an OH-line tilt measurement on the rectified sky frame as a test that the rectification is working (Section 8, test 3). |
| 21 | One slit per 64-lenslet block with per-lenslet wavelength shifts from the Keck coefficients. |
| 23 | No dark subtraction when A-B (or A-sky) subtraction is used; darks only for unpaired frames, grouped by ITIME. |
| 24 | A generic `pypeit_fiber_cube` script dispatching on the spectrograph's fiber hooks; the Binospec scripts become thin aliases; no OSIRIS-specific cube script. |
| 25 | Keep all 73 KOA frames for now; trim the dev-suite copy later. |
| M1-1 | Standards are typed from the target name (spectral-type tag such as `(A0V)`, or a bright-star catalogue name with exptime <= 60 s); science and standard are exclusive. |
| M1-2 | `exptime = TRUITIME x COADDS` when `TRUITIME` (actual per-coadd integration, seconds) is present; `ITIME` rules otherwise. |
| M1-3 | In the dev-suite pypeit files, `ISSKY=1` frames are typed `sky` only (not `science`), so they serve as backgrounds and are not reduced on their own; `pypeit_setup` still types them `science,sky`. |
| M1-4 | `specflip = True`; the M3 `rectify` step un-orients the image (transpose + flip), solves in the raw orientation, and flips the spectral axis of the output (M2 finding). |
| M2-1 | Production `lstsq` uses flux-independent weights (read noise + smoothed sky/background variance, not the object's shot noise) and propagates the true pixel variance with the sandwich formula; `drp` remains the default until this is in place. |
| M2-2 | Cosmic rays are flagged on the raw (A-B) frame *before* `rectify` and passed in the BPM: PypeIt's LACosmic with conservative settings first, the 2016-DRP strip rule as a fallback if OH lines are eaten. |
| M2-3 | The rectified `ivar` is always computed from the least-squares normal matrix, whichever solver produced the image; the DRP noise rule is kept in `rectify_drp` for fidelity only. |
| M2-4 | Missing `DTMP7`: no grating-temperature correction (warning), never the DRP's 0 K extrapolation; the per-column reference-row correction (rows 2032-2044) is applied to every raw frame in the `rectify` step, before A-B. |
| M3-1 | The OH-line tilt criterion (< 0.05 px over a block) is accepted open at M3 (0.073 px median, a smooth trend of the Keck per-lenslet coefficients, not a rectification error); it is closed in M4 by the OH refinement (per-lenslet zero point or per-block linear term) and re-run on the refined solution. |
| M3-2 | Lenslets within a block are staggered by 1.92 px per lenslet row; every lenslet carries its exact Keck solution through the `get_fiber_wavelengths` hook (spec1d exact).  In M4 `get_fixed_tilts` returns Keck-derived tilts (2D Legendre per block, 0.02-0.03 px rms) so the spec2d wavelength image and 2D sky models are exact too; per-lenslet offsets keep the residual. |
| M3-3 | Narrowband filters: the rectified frame is re-laid out as a *lenslet image* with one spatial pixel per lenslet (one block per Keck lenslet column: 19 broadband, 51 narrowband); one spec1d spectrum per lenslet (3040 for Hn3), `MASKDEF_ID = row * 51 + col`. |
| M3-4 | The `DTMP7` grating-temperature term is kept (it worsens the 2016 Hn3 solution by +0.4 A); it is re-tested on the 2024 Kbb/050 frames in M4 and, if it also degrades them, made switchable by a parameter or absorbed by the OH refinement. |
| M3-5 | `osiris_rect_method` defaults to `lstsq` (flux-independent weights, sandwich variance; quasar within 0.03 % of `drp`). |
| M3-6 | Reference-row correction and L.A.Cosmic run on each raw frame (A and B separately), not on the raw A-B difference; the DRP strip rule is an option only. |
| M3-7 | The generic shared-code additions of M3 (the `rectify` step, base hooks `rectify`/`get_fixed_tilts`/`get_fixed_wavecalib`/`get_fiber_wavelengths`/`set_runtime_par`, `method = 'fixed'`, per-fiber wavelengths and the column-window boxcar in the Fiber pypeline) are accepted; the `runtime_par` mechanism and the `FiberExtract` speed-up are flagged for the PR reviewers. |
| M4-1 | Rectification matrices are selected by Keck **validity epoch** (`rectmat_epochs.ecsv`, from http://tkserver.keck.hawaii.edu/osiris/): a matrix whose file date lies inside the epoch of MJD-OBS, else latest-before with a warning; `.fits.gz` accepted. |
| M4-2 | `get_fixed_tilts` returns Keck-derived tilts (2D Legendre (6, 4) per block, fit inside the filter band); residual to the per-lenslet Keck solutions 0.044 px rms Kbb, 0.048 px Hn3 (0.021 px robust); the per-lenslet offsets carry the rest. The base hook gains an optional `raw_header`. |
| M4-3 | Broadband lenslet windows are the filter band (DRP cube range); pixels outside it are masked in the lenslet image. |
| M4-4 | Per-lenslet OH zero points on every science frame (`[calibrations][wavelengths] fiber_sky_refine = zeropoint`, list `OH_OSIRIS` in K, `OH_NIRES` otherwise); unmeasured lenslets get the frame median, very short frames the shift of the median-stacked spectrum.  Closes the M3 tilt criterion (0.023 px Hn3, 0.020 px Kbb). |
| M4-5 | The `DTMP7` term moves the OH lines by +0.40 A on both epochs (2016 Hn3, 2024 Kbb): new `osiris_grating_temp` parameter, default False; the refinement absorbs the remaining epoch offset (+1.5 A on 2024 Kbb). |
| M4-6 | Davies (2007) OH scaling: `[reduce][skysub] osiris_oh_scaling` (`global` default for OSIRIS, `block`, `fiber`), generic `skysub.oh_scaled_sky`, applied in `FiberFindObjects` on A-B frames (the un-subtracted A frame is passed to object finding); not applied to standards. |
| M4-7 | Darks are never rectified (`buildimage` passes `rectify=False` for bias/dark frames), so they can be subtracted from unpaired raw frames; the dev-suite file groups darks by ITIME in calibration groups 1-3. |
| M4-8 | Generic `pypeit_fiber_cube` (hook `fiber_cube_setup`, WCS referenced to the layout origin, optional point-source extraction to a single-object spec1d); `pypeit_binospec_ifu_cube` is a thin alias. |
| M4-9 | Telluric/flux with the IR algorithm on the point-source spectrum of the A0V star (`pypeit_fiber_cube --point_source 0.4`; dev suite: HD101060 frame `OS.20240203.59407.32`, 1.5 % rms flatness); instrument-response polynomial order 7 (`sensfunc` and `telluric` `polyorder`) for the Kbb filter ripples. |
| M4-10 | Davies OH scaling guard (M5): each group factor is judged on the *integrated* line residual (Report 06, 5.3) rather than the line-pixel rms, and/or the scaling is skipped below a minimum sky-line S/N (e.g. `EXPTIME < 100 s`); `global` stays the default.  Motivation: the scaling improves 12/13 600-s frames but worsens the two 44-s frames. |
| M4-11 | The non-physical horizontal step at flipped lenslet row ~37 in the H2 and Br gamma maps (continuum and CO index are fine) is diagnosed in M5: map the raw detector column of the line per lenslet (hypothesis: the 1.92 px/row stagger carries the line or a side band across column 1024 / a readout channel and an additive offset biases the continuum subtraction) and compare with a DRP reduction of one frame. |
| M4-12 | OH refinement for frames too short to show sky lines (1.5-s HD89239; no line even in the stack): use a science frame of the same night and configuration as the "arc", i.e. apply its per-lenslet OH shifts (fallback: the median shift of the night's science frames), with a warning; never leave the raw Keck solution (+1.5 A off on this epoch).  M5. |
| M4-13 | Accepted as is: the Keck-derived tilt rms close to the 0.05 px limit (dominated by a few lenslets with peculiar coefficients, which keep their exact per-lenslet solution; unit test asserts rms < 0.05, robust < 0.03 px); the telluric standard frame `OS.20240203.59407.32`; and the M4 shared-code changes (listed in M5 task 6 for the reviewers). |
| 22 | KOA 2024-02-03 night (NGC 3245 + HD101060/HD89239 A0V + darks, Kbb/050) adopted as the science-mode setup; frames downloaded to `RAW_DATA/keck_osiris/Kbb_050/`; Kbb/050 matrix to be requested from Keck (`Reports/02_koa_dataset_and_matrix_request.md`). |

## 3. Processing flow

```
raw SCI/VAR/DQ (ADU/s)                       Keck calibrations
      |                                      - rectification matrix (filter, scale, epoch)
      v                                      - wave coeffs (epoch), OH line list
[A] get_rawimage: -> DN, exptime, BPM        - static hot/dead masks (pypeit/data)
      |
      v
[B] RawImage.process: dark subtraction (optional), BPM; inside the
      |               `rectify` step, per raw frame: reference-row
      |               correction, L.A.Cosmic (no bias, overscan, flat)
      v
[C] A-B / A-sky frame subtraction (PypeIt bkg_id machinery, raw-frame level)
      |
      v
[D] Rectification: (2048 x 2048 detector) -> (2048 spec x nslice)
      |            per-column linear solve with the influence matrix;
      |            output image, ivar and mask, re-laid out as a
      |            *lenslet image* (one spatial pixel per lenslet;
      |            identity for broadband); "detector" image for the
      |            rest of PypeIt is this rectified frame
      v
[E] Fiber pypeline calibrations:
      slits  = block-slits defined geometrically (64 slices per block)
      tilts  = identity (M3); Keck-derived 2D Legendre per block (M4)
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
| `exptime` | compound: `TRUITIME x COADDS` when `TRUITIME` is present (s per coadd; 2024 headers); else `ITIME x COADDS`, ITIME/1000 when `MJD-OBS >= 57388` (fallback: `ITIME > 2000`) | `TRUITIME` absent in 2015/2016; `ITIME` overstates short 2024 frames by up to 36 % (M1 Q&A 2) |
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
| `science` | `idname == 'a'` and `dispname != 'Drk'` and not `standard` (all on-sky frames; pairing via `bkg_id`) |
| `sky` | on-sky **and** `issky == 1` (helps the user fill `bkg_id`); in the dev-suite pypeit files the `science` type is removed from these frames so they are only used as backgrounds (M1 Q&A 3) |
| `standard` | `idname == 's'` **or** target name with a spectral-type tag (`HD101060(A0V)`) or a bright-star catalogue name (HD, HIP, HR, BD, ...) with `exptime <= 60 s`; exclusive with `science` (M1 Q&A 1; `OBSTYPE` is `'a'` for every KOA frame) |
| `trace` | same as science (needed so `IFUCalibrations.get_slits` has a frame to build block-slits from; see 5.4) |
| `arc`, `tilt` | same as science (sky frame if available) - the OH lines are the wavelength reference for the optional refinement |
| `bias`, `pixelflat`, `illumflat`, `align`, `pinhole`, `scattlight` | never |

An era check in `get_rawimage` raises `PypeItError` for `mjd < 57388`.

### 4.4 Detector (`get_detector_par`)

```python
dict(det=1, binning='1,1', dataext=0, specaxis=1,
     specflip=True,   # Keck coefficients: wavelength decreases with column
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
2. Era check; compute `exptime` (`TRUITIME x COADDS` if present, else `ITIME`); multiply `SCI` by `exptime` to get DN.
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
| `get_fiber_blocks(det)` | list of dicts (`block_id`, `nfibers`, `type`, `fiber_positions`, `fiber_names`, `fiber_ids`, `min_pix`, `max_pix`) | one block per Keck lenslet column present for the filter (19 broadband, 51 narrowband) of up to 64 lenslets each; positions = lenslet pixel + 0.5 within the lenslet image (`pos = 64 * block + slice % 64`); names `L<col>_<row>`, ids `row * 51 + col`; type `'science'` for all (no sky fibers) |
| `get_block_slit_edges(traceimg, det)` | `(left, right)` arrays `(nspec, nblocks)` | constant edges at `64*k - 0.5` and `64*(k+1) - 0.5`; no cross-correlation needed because the rectified geometry is exact |
| `identify_fibers_in_block(det, block_idx, detected_positions)` | ids/names/types aligned with detections | nearest-integer match to slice index |
| `get_fiber_metadata(det, slit_spat_ids, slit_centers)` | ids/names/types per slit | same mapping, applied after extraction |
| `load_sky_layout()` | `(x_arcsec, y_arcsec)` for every lenslet | from the slice-to-lenslet map (`kbb_2016_slice_to_lenslet.txt` generalised to all filters via the DRP `assembcube` rules) times `SSCALE` |
| `get_science_fiber_layout_indices(det, fiber_ids, fiber_types)` | indices into the layout | identity on lenslet id |
| `ifu_sky_wcs(raw_hdr, scale)` | `(SkyCoord, CD)` | `RA`/`DEC`, position angle from `ROTPOSN`/`INSTANGL`, Keck I flip |
| `get_wcs`, `get_datacube_bins` | 3-D WCS and bins | as KCWI but with lenslet grid (66 x 51 or the filter's subset) |

Narrowband filters: the DRP packs up to three lenslet-column blocks head
to tail along the same 2048-px slice.  M3 implements the split rule from
`assembcube_000.pro` (column offsets 0, +16, +32 and the `sp > 831`
remainder) by re-laying the rectified frame out as a *lenslet image*:
one spatial pixel per lenslet, one block per Keck lenslet column, so
every spectrum has one wavelength solution and one sky position (M3
Q&A 3).  Broadband filters reduce to the identity (1216 pixels).  The
hooks need the filter/date, recorded on the spectrograph by `rectify` or
`set_lenslet_state` (cube scripts must call the latter).

### 4.8 Parameters (`default_pypeit_par`)

```python
par['rdx']['...']                                   # nothing special
par['calibrations']['darkframe']['exprng'] = [None, None]
par['calibrations']['biasframe']   -> never used
par['calibrations']['slitedges']   -> bypassed by get_block_slit_edges
par['calibrations']['tilts']       -> see 5.4 (identity tilts)
par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']  # 'OH_OSIRIS' in K (M4)
par['calibrations']['wavelengths']['reference'] = 'sky'
par['calibrations']['wavelengths']['fiber_sky_refine'] = 'zeropoint'  # NEW (M4)
par['calibrations']['osiris_grating_temp'] = False  # NEW (M4): DTMP7 term off
par['calibrations']['wavelengths']['fwhm'] = 2.0
par['calibrations']['flatfield']['method'] = 'skip'
par['calibrations']['osiris_rectmat_dir'] = None   # NEW: directory of rect matrices
par['calibrations']['osiris_rect_method'] = 'lstsq' # NEW: 'drp' | 'lstsq'; default lstsq (M3)
par['calibrations']['osiris_rect_niter'] = 40      # NEW
par['calibrations']['osiris_rect_relax'] = 1.0     # NEW (DRP 'relaxation')
par['calibrations']['osiris_rect_use_var'] = False # NEW: weight 'lstsq' by VAR (M2)
par['calibrations']['osiris_rect_crmask'] = 'lacosmic'  # NEW: CR flagging per raw frame (M3)
par['calibrations']['wavelengths']['method'] = 'fixed'  # NEW option: spectrograph-provided WaveCalib (M3)
par['scienceframe']['process'] : use_biasimage=False, use_overscan=False,
      use_darkimage=False (no dark subtraction for A-B pairs, Q&A 23;
      darks grouped by ITIME only for unpaired frames), use_pixelflat=False,
      use_illumflat=False,
      mask_cr=False (post-rectification LACosmic is useless: cosmic rays
      are flagged on the raw A-B frame inside the `rectify` step and
      passed in the BPM, M2 Q&A 2 / Section 5.3)
par['reduce']['findobj']  -> fiber mode (no peak finding)
par['reduce']['skysub']['joint_fit'] = False
par['reduce']['skysub']['osiris_oh_scaling'] = 'global'  # NEW (Davies 2007, M4);
      # M5: factors judged on the integrated line residual and/or skipped
      # below a minimum sky-line S/N (short frames), decision M4-10
par['reduce']['extraction']['boxcar_radius'] = 0.5 * platescale (one slice)
par['reduce']['cube']['slit_spec'] = False
par['flexure']['spec_method'] = 'skip'
par['sensfunc']['algorithm'] = 'IR'
par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'
par['sensfunc']['polyorder'] = 7 ; par['telluric']['polyorder'] = 7   # M4
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
| (c) Batched per-band solve: within each 32-row band only ~16-19 slices contribute, but adjacent bands are coupled through the slice profiles, so each column is an exact **block-tridiagonal** system (64 blocks of <= 19 unknowns) solved in one batched call | fast and exact; natural ivar propagation from the normal matrix | independent bands are wrong by up to 5x at band edges (M2 finding); Poisson weighting biases bright sources (see below) |

**Recommendation:** implement (a) first as the validation reference,
then (c) as the production solver, keeping both behind a parameter
(`osiris_rect_method`, M2).

**Weighting and variance (M2 Q&A 1, 3).**  With Poisson + read-noise
weights, (c) under-estimates bright point sources by ~2 % (6 % at the
red end of Hn3) because the Keck PSF model leaves residuals of up to
~40 sigma at the core and down-weighting the core moves flux.  The
production solver therefore uses *flux-independent* weights (read noise
plus a smoothed sky/background variance, not the object's own shot
noise) and propagates the true pixel variance with the sandwich formula
`(A^T W A)^-1 A^T W V W A (A^T W A)^-1` (diagonal only, same
block-tridiagonal machinery).  The rectified `ivar` is always taken from
this least-squares machinery, whichever solver produced the image: the
DRP noise rule (linear sum of `fblame * sigma`, not divided by 1.28) is
0.53x the least-squares sigma and is kept in `rectify_drp` only for
fidelity.  `drp` stays the default method until the weighted solver is
in place (M3).
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

   Pre-processing inside the step (M2 Q&A 2, 4):

   - *Orientation*: both solvers work in the raw orientation and return
     `(2048, nslice)` in raw column order; the step un-orients the
     image (transpose + flip, undoing `specflip`), solves, and flips the
     spectral axis of the result rather than transforming the matrix.
   - *Reference rows*: subtract the per-column median of rows 2032-2044
     (`lenslet_rectify.subtract_reference_rows`, the DRP's H2RG
     `adjchan` behaviour) from every raw frame, before A-B; otherwise
     the solver assigns the column offsets to lenslets.
   - *Cosmic rays*: flag them on the raw (A-B) frame before solving and
     pass them in the BPM (masked pixels are handled gracefully by both
     solvers; after rectification a hit is spread over neighbouring
     slices and cannot be found).  Use PypeIt's LACosmic on the raw
     difference with conservative settings; fall back to the 2016-DRP
     strip rule (ported in the M2 vet test) if LACosmic eats OH lines.
     The post-rectification `mask_cr` stays off.

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
  *M3 status:* `get_fixed_tilts` exists and returns the identity; but
  lenslets within a block are staggered by 1.92 px per lenslet row, so
  in M4 the hook returns Keck-derived tilts (a 2D Legendre per block,
  0.02-0.03 px rms) and the per-lenslet offsets carry only the residual
  (M3 Q&A 2).  This makes the spec2d wavelength image exact and lets a
  2D sky model (Davies scaling) work without special cases.
  *M4 status:* done.  The tilts map every pixel of a block onto the
  spectral pixel of the reference lenslet at the same wavelength (fit
  inside the filter band, Legendre orders (6, 4)); the wavelength image
  reproduces the per-lenslet Keck solutions to 0.044 px rms (Kbb) and
  0.048 px (Hn3), 0.021 px robust; a handful of lenslets with peculiar
  coefficients dominate the tail and keep their exact solution through
  the per-lenslet offsets.
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
   offsets within a block are **not** small: OH lines run at 1.92 px
   per lenslet row (~120 px across a 64-lenslet block) and a pure shift
   leaves up to 0.6 px (M3 finding).  One slit per block is kept, but
   every lenslet receives its exact Keck solution through the
   `Spectrograph.get_fiber_wavelengths` hook (applied as per-fiber
   wavelength arrays by the Fiber pypeline), so the spec1d wavelengths
   are exact.  From M4 the block tilts are Keck-derived (Section 5.4)
   so the spec2d wavelength image is exact as well.
3. The science header's `DTMP7` applies the grating-temperature scaling;
   when absent (truncated header) use the epoch default (no scaling) and
   warn.  Do **not** copy the DRP, which reads a missing `DTMP7` as 0 K
   and extrapolates the CTE table (a ~1.5 px shift in the 2016
   reference cube; M2 Q&A 4).  On the 2016 Hn3 frame the term *worsens*
   the OH-line agreement (+0.40 A with it, -0.01 A without); it is kept
   for now and re-tested on the 2024 Kbb/050 frames in M4, after which
   it becomes switchable or is absorbed by the OH refinement (M3 Q&A
   4).  *M4 decision:* on the 2024 Kbb/050 frames the term also moves the
   OH lines by +0.40 A (+1.92 A with, +1.52 A without, against the
   vacuum `OH_OSIRIS` list), so it is now off by default
   (`osiris_grating_temp = False`) and the per-lenslet OH refinement
   absorbs the remaining epoch offset.  `resample_to_drp_grid` reproduces that
   behaviour only for the comparison with DRP cubes.

### 6.2 OH refinement (phase 2)

With `reference = 'sky'` and an OH line list (`ohlines_kbb_all.dat` for
K; PypeIt's `OH_NIRES`/`OH_FIRE_Echelle` lists for J/H), run the
existing fiber shift solver on the A frame (before sky subtraction) to
measure a zero-point (and optionally linear) correction per lenslet; the
DRP documents 0.3-0.6 A offsets depending on scale and epoch.  Where
fewer than 3 lines are available (long-K, narrowband) keep the Keck
solution.  Report the per-lenslet shift statistics in QA.

*M4 status:* implemented as a generic Fiber-pypeline step
(`[calibrations][wavelengths] fiber_sky_refine`; `wvutils.fiber_sky_refinement`):
isolated lines of the list (no neighbour within three FWHM) are centroided
in each lenslet of the un-subtracted science frame and the median offset
from the prediction is applied as a pixel zero point (or, `block`, a linear
fit across each block).  Lenslets with < 3 lines get the frame median;
frames with < 10 measured lenslets (short standards) the shift of the
median-stacked spectrum.  The per-lenslet shifts are reproducible between
frames (r = 0.97 Hn3, 0.90 Kbb; 0.014-0.024 px difference), i.e. they are
real structure of the Keck coefficients; QA `QA/PNGs/*_fiber_sky_refine.png`.
Frames with no measurable line at all (the 1.5-s HD89239 standards)
currently keep the Keck solution (+1.5 A off); M5 makes a science frame
of the same night/configuration the "arc" for such frames (its
per-lenslet shifts, else the night median) with a warning (M4-12).

### 6.3 Vacuum

Keck solutions are vacuum.  Set the PypeIt wavelength frame accordingly
(no air-to-vacuum conversion) and record it in the spec1d header.

## 7. Calibration-file handling

| File | Size | Location | Selection |
|---|---|---|---|
| Rectification matrix | ~160 MB per (filter, scale, epoch); 30-90 MB gzipped | **Source:** Keck's public listing by validity epoch, http://tkserver.keck.hawaii.edu/osiris/ (plain HTTP; files are `.fits.gz`).  User directory via `osiris_rectmat_dir`; dev-suite copies in `RAW_DATA/keck_osiris/<setup>/calib/` | parse `sYYMMDD_cNNN___infl_<filt>_<scale>.fits[.gz]`; select by **validity epoch** (shipped `rectmat_epochs.ecsv` copied from the tkserver table): the epoch containing MJD-OBS, then a matrix of that filter/scale whose file date lies inside the epoch (the file date can be *after* the observation, e.g. `s240524` serves 2024-02-03); fall back to latest-before with a warning; error listing what was found otherwise (M4 task 11) |
| Wavelength coefficients (11 epochs) + CTE table | 0.6 MB total | `pypeit/data/spectrographs/keck_osiris/` | epoch table copied from `calibrations.xml` |
| Hot/dead pixel masks | ~100 kB gz | `pypeit/data/static_calibs/keck_osiris/` | always |
| Slice-to-lenslet maps | text | `pypeit/data/spectrographs/keck_osiris/` | per filter family (bb / nb) |
| OH line lists | text | `pypeit/data/arc_lines/lists/` (new `OH_OSIRIS_K.dat` from `ohlines_kbb_all.dat`) | by band |

The matrix loader caches the parsed `(hilo, effective, influence)`
tuple in memory for the run (one matrix per configuration).  A later
convenience could add a PypeIt-cache download of the common-mode
matrices from the tkserver URLs above (they are stable and public, so
the hosting question in Report 02 is moot).

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
3. **OH-line tilt test (rectification check, requested in Q&A 20)**:
   run PypeIt's standard tilt tracing (`tracewave`) on the rectified
   sky-only frame (frame B of the Hn3 pair, and an `ISSKY=1` NGC 3245
   frame) with the OH lines as the tilt source.  If the rectification is
   correct every OH line must be a straight, vertical feature across all
   64 lenslets of a block: the fitted tilt slope must be < 0.05 px over
   the block and the per-lenslet line-centroid scatter < 0.1 px after
   removing the Keck per-lenslet shift.  Any residual tilt or saw-tooth
   pattern points at a wrong `hilo` offset, a slice/lenslet mis-mapping
   or an unconverged solver.  This test is part of the dev-suite vet
   tests and its QA plot (line centroid vs lenslet row) is kept in the
   Reports folder.  *M3 result:* scatter 0.034 px (passes); slope
   0.073 px median (fails), a smooth signed trend of the Keck
   coefficients identical for both solvers and frames.  Accepted open
   with the strict assertion as `xfail` plus a regression guard (median
   < 0.2 px, offset < 0.1 px); to be closed in M4 on the OH-refined
   solution (M3 Q&A 1).  *M4 result:* on the solution refined on the
   paired science frame the median tilt is 0.023 px (Hn3, sky frame 011)
   and 0.020 px (Kbb, `ISSKY=1` frame `OS.20240203.43383.79`); the
   `xfail` is removed.
4. **Dev-suite reduce tests** `keck_osiris/Hn3_100` (pypeit file with
   the pair, `bkg_id` set) and `keck_osiris/Kbb_050` (NGC 3245 with its
   A0V telluric and darks; needs the Keck matrix) for the full flow
   including cube and telluric correction.
5. **Load-images unit test** entry in `unit_tests/test_load_images.py`.
6. **Header-era test**: the 2015 Kn5 files are read by `pypeit_setup`
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

**Decision (Q&A 22):** the 2024-02-03 NGC 3245 + HD101060 set
(Kbb/050) is the science-mode dev-suite setup.  The night's spectrograph
frames (117 in total) were listed through the KOA TAP service; 73 were
downloaded through the KOA `getKOA` endpoint into
`RAW_DATA/keck_osiris/Kbb_050/` (file list kept beside the data):

| Target | Frames | ITIME | ISSKY=1 | Role |
|---|---|---|---|---|
| NGC 3245 | 19 + 4 | 600 s (19), 45 s (4) | 8 | science, 20-min A-B pairs |
| HD101060 (A0V) | 24 | 10 s (18), 3 s (6) | 3 | telluric standard, nodded |
| HD89239 (A0V) | 7 | 2 s | 3 | second telluric standard |
| darks | 15 | 600, 10, 3, 2 s (matched to the above) | - | `SFILTER='Drk'`, 20 mas scale |

The 2024 headers are complete (526 cards; `INSTRUME`, `ISSKY`,
`PA_SPEC`, `DTMP7`, `WXPRESS` all present; `DQ` all 9).  The matching
Kbb/050 rectification matrix (epoch 2023-2024) must be requested from
Keck; the draft request is in `Reports/02_koa_dataset_and_matrix_request.md`.

## 10. Open issues and phases

### Resolved (Q&A 18-22, 2026-10-10)

1. Narrowband split in phase 1; Kbb/050 first supported mode.
2. Rectification as a generic `RawImage` step with solvers in
   `pypeit/core/lenslet_rectify.py`.
3. Identity tilts through a `get_fixed_tilts` hook, verified by the
   OH-line tilt test.
4. One slit per block, per-lenslet shifts.
5. KOA 2024-02-03 set downloaded; matrix request drafted.

### Resolved (Q&A 23-25, 2026-10-10)

- No dark subtraction with A-B pairs (Q23).
- Generic `pypeit_fiber_cube` script; Binospec names kept as aliases (Q24).
- All 73 KOA frames kept for now; dev-suite copy trimmed later (Q25).

The implementation is driven by the milestone prompts in
`claude_prompts/keck_osiris/M1_*.md` ... `M5_*.md`.

### Phases

| Phase | Deliverable | Depends on |
|---|---|---|
| 1 | `keck_osiris.py` metadata, typing, detector, `get_rawimage`, BPM, era check; static files; unit tests; `pypeit_setup` runs on both test sets | - |
| 2 | `pypeit/core/lenslet_rectify.py` (solver a and c), matrix loader/selection, `RawImage.rectify` step; vet test against the DRP cube | Q18, Q19 |
| 3 | Fiber hooks (blocks, ids, layout, WCS), identity tilts, Keck-coefficient `WaveCalib`; full run to spec1d on Hn3/100 | Q20, Q21 |
| 4 | Datacube, A-B and Davies OH-scaled sky subtraction, OH wavelength refinement; Kbb/050 KOA set; dev-suite reduce test; docs and changelog | Q22, matrix from Keck |
| 5 | Narrowband generalisation, telluric/flux calibration, `supported = True` | 4 |
