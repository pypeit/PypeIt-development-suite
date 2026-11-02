# HostSub_GP + PypeIt: current workflow and integration options

- **Version:** 0.1
- **Date:** 2026-10-18
- **Author:** LK and Claude
- **Status:** Draft for discussion. No code has been written.

## 0. Sources read

- `dev_hostsub/refs/HostSub_GP/` (HEAD `d456a49`, "update: likelihood
  evaluation in GP"): README, `docs/cookbook.rst`,
  `docs/reference/input_file.rst`, `pyproject.toml`,
  `src/hostsub_gp/{inputfiles.py, spec_proc.py, _utils/_par.py}`,
  `src/hostsub_gp/scripts/{hostsub_setup.py, hostsub.py, scriptbase.py}`.
- Paper: Liu & Miller 2025, [arXiv:2508.15278](https://arxiv.org/abs/2508.15278)
  (`dev_hostsub/refs/liu_miller_2025.pdf`). Not re-read in detail for this
  doc. The integration questions below are about plumbing, not the method.
- PypeIt on this fork (`hostsub_gp`, the same commit as `develop` `f555483dd`,
  `git describe` = `2.0.1-1356`): `CLAUDE.md`, `doc/running.rst`,
  `pypeit/par/pypeitpar.py` (`PypeItPar`, `ReduxPar`, `SkySubPar`),
  `pypeit/inputfiles.py`, `pypeit/outputfiles.py`, `pypeit/pypeit.py`,
  `pypeit/scripts/run_pypeit.py`, `pypeit/scripts/coadd_2dspec.py`.
  Where this doc says "release", it means the local `release` branch
  (2.0.1, 2026-09-08). Where it says "develop", it means `develop`
  (2026-10-08).

---

## 1. The current two-pipeline workflow, end to end

```
 raw/ ──run_pypeit──► Science/spec2d_*.fits, spec1d_*.fits   (PypeIt)
                          │
          hostsub_setup pypeit -c X.pypeit --target T          (HostSub_GP)
                          ▼
                 HostSub_T/hostsub.txt  ◄── user hand-edits ra/dec/z/seeing/raw_dir
                          │
          cd HostSub_T; hostsub hostsub.txt [--coadd2d ...]    (HostSub_GP)
                          ▼
   ../Science_hostsub/spec2d_*_{preproc,rect,hostsub}.fits, HostSub_T/QA/*.pdf,
   HostSub_T/{ps1,ls,sdss}_cutout/  [+ spec1d_<coadd>_hostsub.fits if --coadd2d]
                          │
          pypeit_coadd_2dspec X.coadd2d (on *_hostsub.fits)   (PypeIt)
                          ▼
                 Science_coadd/ spec1d/spec2d of the transient
```

### 1.1 Step 1: `run_pypeit` (PypeIt)

A standard MultiSlit longslit reduction. The cookbook recommends two
settings in the `.pypeit` file:

| Setting | PypeIt default | Recommended | Why |
|---|---|---|---|
| `[calibrations][[tilts]] tracethresh` | 20.0 | 5.0 | More arc lines go into the tilt fit, so the 2D wavelength solution is better. HostSub's rectification depends on `waveimg`/`tilts`. |
| `[reduce][[skysub]] no_local_sky` | False | True | PypeIt's local sky step flags bright host pixels as CRs/outliers. HostSub models the host+sky itself. |

There are also implicit requirements on the reduction:

- **A standard star is effectively required.** It is used for the trace
  when the science object is not found, and for `spat_resln` (seeing) when
  the user does not set it (`FWHM` x pixel scale of the highest-S2N std
  object). It is also used in `hostsub_setup` to pick the science object
  (see 1.2).
- **A science object should be found** in each exposure. Its `TRACE_SPAT`
  defines the "distance from trace" coordinate. With `no_local_sky=True`,
  finding faint curved traces can fail. The cookbook warns about this, and
  the code falls back to the standard trace.
- **Only `keck_lris_{blue,red,red_mark4}`, `mmt_binospec`, `not_alfosc`
  and `vlt_fors2` are supported.** `SpecData.from_pypeit` hard-codes the
  detector, pixel scale, PA and slit-width parsing per `PYP_SPEC`, and
  raises `NotImplementedError` for anything else. For Binospec, ALFOSC and
  FORS2 it reads the **raw** header (PA, rotator, mask), so it needs
  `raw_dir`.
- Products used: `Science/spec2d_*.fits` (sciimg, ivarraw, waveimg, tilts,
  bpmmask, wavesol, primary header `FILENAME`, `PYP_SPEC`, `RA`, `DEC`,
  `BINNING`, `ROTPOSN`, `SLITNAME`, `DECKER`) and `Science/spec1d_*.fits`
  (SpecObj `NAME`, `TRACE_SPAT`, `S2N`, `FWHM`).

### 1.2 Step 2: `hostsub_setup pypeit -c X.pypeit --target T [T2 ...]`

Run from the PypeIt reduction directory.

1. `PypeItFile.from_file(X.pypeit)` reads the data block. The first row with
   `frametype == 'standard'` gives the standard's `target` name.
2. Glob `./Science/spec1d*[0-9].fits`. Science files = those whose
   *filename contains* `T`. Standard files = those whose filename contains
   the std target name. This is substring matching on output filenames, not
   a metadata query. It is fragile (e.g. `T=SN1` matches `SN10`; target
   names with spaces or `+/-` are sanitized in PypeIt filenames).
3. Std object = the highest-`S2N` SpecObj. Its `NAME`
   (`SPATxxxx-SLITyyyy-DETzz`) is split into spat/slit/det ints.
4. Per science file: the object on the same slit+det with the spatial pixel
   nearest the standard's. Its `NAME` goes in `objid`.
5. Write `HostSub_T/hostsub.txt`: a templated `[hostsub]` config block
   (section 3) plus a `hostsub read ... hostsub end` data block with
   `path ../Science/` and columns `filename | frametype | objid`. Filenames
   are **spec1d** names. The script later derives spec2d names by string
   replacement.

### 1.3 Step 3: hand-edit `hostsub.txt`

At minimum the user sets `ra`, `dec` (header coords can be wrong, e.g.
LRIS), `spat_resln`, `raw_dir`, `[[host_prior]] survey/filters`, and
`[[host_emission]] z`. The template leaves `ra`/`dec`/`z` commented out.

### 1.4 Step 4: `hostsub hostsub.txt [--overwrite] [--skip_seeing_match] [--skip_gp] [--coadd2d] [--write_spec2d]`

Run from inside `HostSub_T/`. QA and cutouts go to the CWD.

1. Parse with `HostSubInput.from_file` (a subclass of PypeIt's `InputFile`;
   the config block is plain ConfigObj, *not* a ParSet, so there is no
   validation and every value is a string cast ad hoc with
   `Float`/`Int`/`Digitize`).
2. Per science exposure: `SpecData.from_pypeit(...)` loads spec2d+spec1d,
   builds a signed "distance from trace" image (arcsec) and a preliminary
   global sky, and rectifies onto a common (spat, wave) grid. It caches
   `Science_hostsub/spec2d_*_preproc.fits` and `*_rect.fits`. All exposures
   share the first exposure's rectification grid.
3. If `--coadd2d`: `SpecData.coadd2d` coadds the rectified frames (HostSub's
   own coadd, not PypeIt's) and writes a CR mask back into each `_rect`.
4. Model (`_model_host_subtraction`), per exposure or once on the coadd:
   - `build_host_prior`: **downloads archival cutouts** (PS1/LS/SDSS) and
     projects the host image onto the slit (needs network).
     `host_image.py` also calls astrometry.net in some paths.
   - `construct_spec_wrapper`: optional host emission-line detection
     (`z`, `p_value`, `kernel_wid`) and adaptive batching (`batch_2d`).
   - Seeing match: fit `dseeing` (bounded by `dseeing_lower/upper`), then
     convolve the spectrum or the prior.
   - `model_host`: tinygp/JAX GP fit (1D x 2D kernel; `params_init_*`,
     `params_limit_*`).
   - QA: `QA/<base>_{host_prior,raw,conv,host_profile_prior,pred,host_profile_pred,sci}.pdf`
     and `QA/<base>_sci.txt` (HostSub's own 1D extraction).
   - `--skip_gp` stops after the "classic" extraction.
5. Write-back:
   - Without `--coadd2d` (or with `--write_spec2d`):
     `update_pypeit_skymodel` writes `Science_hostsub/spec2d_*_hostsub.fits`.
     This is a copy of the PypeIt spec2d with `skymodel` replaced by
     (preliminary global sky + post-rectification global sky + GP local
     host/sky model) and `bpmmask` OR-ed with OFFSLITS outside the modeled
     region and with the CR mask. **Caveat in the code:** in the per-exposure
     mode, the *last* exposure's model is used for every rewrite (comment:
     "intentionally ... expected to be coadded").
   - With `--coadd2d`: also writes
     `Science_hostsub/spec1d_<coadd_basename>_hostsub.fits` (+ `.txt`) by
     building a fresh `SpecObj` and running PypeIt's
     `extract.extract_boxcar` / `extract.fit_profile` /
     `skysub.optimal_bkpts` / `skysub.skyoptimal` on the rectified coadd
     (`pypeit_local_sky` sub-dict controls this).
6. Doc/code mismatch: the cookbook says outputs go to `Science/` with
   suffix `_hostsub.fits`. The code writes to a sibling `Science_hostsub/`.

### 1.5 Step 5: `pypeit_coadd_2dspec` on the `*_hostsub.fits` files

The user writes a `.coadd2d` file listing the `_hostsub` spec2d files.
Because the replaced skymodel already removes the host, object finding in
the coadd sees the transient. Offsets/weights from the original spec1d are
often unusable, so the cookbook recommends `[coadd2d] offsets = 0,0,...`.
PypeIt's coadd2d finds and extracts the transient (with PypeIt's local sky
on by default in coadd2d, which partly reverses the `no_local_sky` choice).

### 1.6 Pain points the integration should fix

1. Two config languages (`.pypeit` ParSet vs untyped ConfigObj) and two CLIs.
2. Target/standard/object association by filename substring and
   `NAME` string parsing, instead of `PypeItMetaData` / `outputfiles`
   helpers.
3. Hand-editing a generated file, so it is not reproducible from the
   `.pypeit` file alone.
4. Spectrograph-specific geometry (pixel scale, PA, slit width, det) is
   hard-coded in HostSub, duplicating what the PypeIt `Spectrograph`
   classes know.
5. Output paths depend on CWD (`QA/`, cutouts), and there is a doc/code
   mismatch on the output directory.

---

## 2. PypeIt APIs imported by HostSub_GP, and their stability on develop

HostSub_GP says it was tested on PypeIt 1.17. Every PypeIt import is lazy
(inside functions) except `hostsub_gp/inputfiles.py`, which imports
`pypeit.inputfiles.InputFile` at module level. So `import hostsub_gp`
itself does not need PypeIt, but both scripts do.

Legend: **OK** = present with a compatible signature on develop.
**Drift** = present, but semantics/format changed in a way that may matter.
**BROKEN** = the call as written fails on develop.

| # | Import / call site (HostSub_GP) | PypeIt API | Used for | develop status |
|---|---|---|---|---|
| 1 | `inputfiles.py` (top level) | `pypeit.inputfiles.InputFile` (subclassed: `data_block`, `flavor`, `setup_required`, `datablock_required`, `required_columns`; `.from_file`, `.config`, `.data`, `.filenames`) | Parse `hostsub.txt` | **OK.** Class attributes and `from_file(input_file, vet, preserve_comments)` are unchanged. `inputfiles.py` had 30 commits since 2025 (+619/-182 vs release), mostly new helpers. A private subclass is still a moderate coupling risk. |
| 2 | `hostsub_setup.py` | `pypeit.inputfiles.PypeItFile.from_file`, `.data`, `.file_paths` | Read the `.pypeit` data block, raw dir | **OK.** |
| 3 | `hostsub_setup.py`, `spec_proc.py` | `pypeit.specobjs.SpecObjs.from_fitsfile(file, det=...)`, iteration, `obj['S2N']`, `obj['NAME']`, `obj['FWHM']`, `obj['TRACE_SPAT']`, `.name_indices(name)` | Choose std/sci objects, trace, seeing | **OK** for the signatures. **Drift:** `NAME` is parsed by string split (`SPATxxxx-SLITyyyy-DETzz`). That is unchanged for MultiSlit on develop, but mosaics (`MSCnn`) and the new Fiber naming differ. `FWHM` is in pixels (`SPAT_FWHM` is arcsec), so HostSub's `FWHM x pixel_scale` is correct. It could just read `SPAT_FWHM`. |
| 4 | `spec_proc.py` (`from_pypeit`, `update_pypeit_skymodel`) | `pypeit.spec2dobj.Spec2DObj.from_file(file, detname)`; attributes `sciimg`, `ivarraw`, `waveimg`, `tilts`, `bpmmask.mask`, `wavesol['measured_fwhm']`, `skymodel` (write) | Load the 2D frame, write back the new sky | **OK.** 6 commits since 2025, +13/-1 vs release. **Drift:** `BUNIT` moved from the per-det header to the primary header when gain was applied. HostSub copies the primary header from the original file, so this is harmless. |
| 5 | `spec_proc.py` | `pypeit.spec2dobj.AllSpec2DObj()`, `[det] = ...`, `.write_to_fits(outfile, pri_hdr=...)` | Write `*_hostsub.fits` | **OK.** The signature is unchanged. Only the one det is written, so other detectors in the original file are dropped. |
| 6 | `spec_proc.py` (hard-coded) | Bit values `CR = 2**1`, `OFFSLITS = 2**4` of `ImageBitMask` | OR new masks into `bpmmask` | **OK** (bit order unchanged), but it should use `ImageBitMask().turn_on(...)` instead of literals. |
| 7 | `scripts/hostsub.py` (`--coadd2d`) | `pypeit.coadd2d.CoAdd2D.default_basename(spec2d_files)` | Name the coadd spec1d | **OK** (static method, unchanged). `coadd2d.py` changed heavily (+1062/-710), but not this method. |
| 8 | `spec_proc.py` (`write_pypeit_spec1d`) | `pypeit.specobj.SpecObj(PYPELINE, DET, OBJTYPE, SLITID)`, `.set_name()`, attributes `BOX_*`, `OPT_*`, `S2N`, `FWHM`, `SPAT_FWHM`, `FWHMFIT`, `TRACE_SPAT`, `SPAT_PIXPOS`, `SPAT_FRACPOS`, `RA`, `DEC`, `DETECTOR`, **`BOX_RADIUS`** | Build the coadd spec1d | **BROKEN on release and develop.** `BOX_RADIUS` was renamed to `BOX_R_PIX` after 1.17.0 (absent in 1.18.0). `DataContainer.__setattr__` raises `AttributeError` for unknown attributes. Only the `--coadd2d` path is affected. |
| 9 | `spec_proc.py` | `pypeit.specobjs.SpecObjs([sobj]).write_to_fits(subheader, outfile, overwrite)`, `.write_info(txt, pypeline)` | Write the coadd spec1d | **OK.** |
| 10 | `spec_proc.py` | `pypeit.images.detector_container.DetectorContainer.from_hdu` | Copy detector metadata into the new SpecObj | **OK** (`DataContainer.from_hdu`). |
| 11 | `spec_proc.py` (`_extract_pypeit_optimal`) | `pypeit.core.extract.extract_boxcar(imgminsky, ivar, mask, waveimg, skyimg, sobj)` | Boxcar on the coadd | **BROKEN on develop.** Commit `db453e9b4` "extract without specobj" (2026-02-26) changed the signature to `extract_boxcar(box_radius, trace_spat, imgminsky, ivar, mask, waveimg, skyimg, ...)` and no longer takes a SpecObj. The SpecObj wrapper is now `SpecObj.extract_boxcar(...)`. OK on release. |
| 12 | same | `pypeit.core.extract.fit_profile(...)` | Profile fit on the coadd | **BROKEN on develop.** It moved to `pypeit.core.spatialprofile.fit_profile` (same commit). The keyword names look unchanged. OK on release. |
| 13 | `spec_proc.py` (`_fit_pypeit_local_sky`) | `pypeit.core.skysub.optimal_bkpts(bkpts_optimal, bsp_min, piximg, sampmask, skyimage=)`, `skysub.skyoptimal(piximg, data, ivar, oprof, sigrej, npoly, spatial_img, fullbkpt)` | Residual local sky on the coadd | **OK.** Signatures are unchanged, though `skysub.py` had 36 commits since 2025. Errors are caught by a broad `except`, so failures are silent. |

Non-PypeIt fragility worth noting: `spec_proc.py` imports
`jax._src.typing` (a private JAX module). HostSub_GP requires Python
>= 3.12, but PypeIt allows 3.11.

**Summary.** The core path (rows 1–7) works on develop. The `--coadd2d`
spec1d writer (rows 8, 11, 12) is broken on develop, and partly broken
even on release. Any integration should either (a) pin and patch upstream
HostSub_GP, or (b) re-implement the thin PypeIt-facing layer (load, trace,
write-back, extraction) inside PypeIt and call HostSub_GP only for the
modeling (`SpecData` rectification → `SpecModel` → GP). Option (b) removes
rows 1–3 and 8–13 from HostSub's surface.

---

## 3. `hostsub.txt` parameters: scope and proposed `.pypeit` home

Scope: **T** = per-target (differs between transients in one `.pypeit`
file). **G** = global (instrument/reduction choice, the same for every
target). **T/G** = usually global, but sometimes tuned per target.

"Proposed home" assumes a new `HostSubPar` registered in `PypeItPar` as
`[hostsub]`. This is required in any case. `PypeItPar.from_dict`
whitelists top-level keys and **raises `ValueError` on an unknown
`[hostsub]` block**, so even approach C needs PypeItPar to know about it.
Per-target values need a home that ParSets do not have today; see the
options after the table.

| hostsub.txt key | Code default | Type | Scope | Proposed .pypeit home |
|---|---|---|---|---|
| *(new)* enable flag | n/a | bool | G | `[hostsub] use_gp = False` (the master switch; False ⇒ output unchanged) |
| `ra`, `dec` | None → header | float (deg) | **T** | Per-target (see options) |
| `[[host_emission]] z` | None | float | **T** | Per-target |
| `[[host_emission]] z_err` | None | float | **T** | Per-target |
| `spat_resln` | None → std FWHM | float (arcsec) | **T** (seeing per night/target) | Per-target, with a global fallback `[hostsub] spat_resln` |
| `sky_offset` | 0.0 | float (arcsec) | T | Per-target |
| `mask_offset` | 0.0 | float (arcsec) | T | Per-target |
| `host_wid` / `host_region` | `[-5, 5]` | float / 2-list (arcsec) | T/G | `[hostsub] host_region`; per-target override |
| `sky_region` | `[-5, 5]` (template: `-10, 10`) | 2-list (arcsec) | T/G | `[hostsub] sky_region` |
| `mask_wid` | 2.0 (template 2.5) | float (x seeing) | G | `[hostsub] mask_wid` |
| `slit_len` | None (template 60) | float (arcsec) | G (instrument) | Derive from `SlitTraceSet` / spectrograph. Keep only as an override. |
| `slit_trim` | `[5, 5]` | 2-int | G | `[hostsub] slit_trim` |
| `spec_range` | None | 2-float (Å) | G | `[hostsub] spec_range` |
| `spec_resln` | None → `wavesol measured_fwhm` | float | G | Drop (derived from `WaveCalib`) or keep as an override |
| `batch_2d` | `[2, 128]` (template 256) | 2-int | G | `[hostsub] batch_2d` |
| `sigma_clip` | 5.0 | float | G | `[hostsub] sigma_clip` |
| `raw_dir` | None (template: first `.pypeit` path) | str | G | **Drop.** It is already in the `.pypeit` data block (`file_paths`). Raw-header items should come from `fitstbl` / spectrograph metadata. |
| `[[host_prior]] survey` | `PS1` | str ∈ {PS1, LS, SDSS} | G (T if dec < −30 ⇒ no PS1) | `[hostsub][[host_prior]] survey` |
| `[[host_prior]] filters` | `grizy` | str | G | `[hostsub][[host_prior]] filters` |
| `[[host_emission]] find_host_emission` | True | bool | G | `[hostsub][[host_emission]] find_host_emission` |
| `[[host_emission]] p_value` | 0.05 | float | G | `[hostsub][[host_emission]] p_value` |
| `[[host_emission]] kernel_wid` | None | float | G | `[hostsub][[host_emission]] kernel_wid` |
| `[[seeing_match]] dseeing` | None (0 ⇒ None) | float | T | Per-target |
| `[[seeing_match]] dseeing_lower`, `dseeing_upper` | 0.0, 1.0 | float | G | `[hostsub][[seeing_match]] ...` |
| `params_init_1d`, `params_init_2d` (dicts: `log_amp`, `log_scale`, `mean`, `log_amp_line`, `scale_line`) | package defaults | dict of float/list | G (expert) | `[hostsub][[gp]] params_init_1d = ...` (ConfigObj subsection; expert-only) |
| `params_limit_1d`, `params_limit_2d` (`<p>_lower`, `<p>_upper`) | none | dict | G (expert) | `[hostsub][[gp]] params_limit_*` |
| `pypeit_local_sky` (`enabled`, `bsp`, `sigrej`, `npoly`, `bkpts_optimal`) | True, 0.6, 3.5, 1, True | dict | G | Reuse the existing `[reduce][[skysub]]` keys (`bspline_spacing`, `sky_sigrej`, ...) instead of a parallel set |
| CLI `--coadd2d` | False | bool | G | `[hostsub] coadd2d` |
| CLI `--skip_seeing_match` | False | bool | G | `[hostsub][[seeing_match]] skip` |
| CLI `--skip_gp` | False | bool | G | Covered by `use_gp = False` (or a `mode = classic`) |
| CLI `--write_spec2d` | False | bool | G | `[hostsub] write_spec2d` |
| CLI `--overwrite` | False | bool | G | Reuse `run_pypeit -o` / `[rdx]`-level overwrite |
| Data-block `frametype` | n/a | str | n/a | Already in the `.pypeit` data block |
| Data-block `objid` | auto (nearest std on same slit) | str | **T**, per frame | Auto by default. Optional override via a new data-block column, e.g. `hostsub_objid`, or reuse `manual` |
| *(implicit)* target selection (`--target`) | n/a | str | T | The `.pypeit` `target` column. Which targets to process: per-target section names or a `[hostsub] targets` list |

**Options for per-target values** (ra, dec, z, z_err, spat_resln,
offsets, dseeing):

- **P1: data-block columns** (`hostsub_ra`, `hostsub_dec`, `hostsub_z`,
  ...). Per frame, so they are repeated for every exposure of a target.
  This fits PypeIt's existing user-column pattern
  (`PypeItMetaData.set_user_added_columns`: `manual`, `shift`). It makes
  the data block wide and puts science metadata in a calibration-oriented
  table.
- **P2: named subsections** `[hostsub][[targets]][[[SN2026abc]]] ra = ...`.
  This is natural ConfigObj, but `ParSet` has no "dict of child ParSets
  keyed by arbitrary names". It needs a small new container (precedent:
  none exact; closest is how `Coadd2DPar`/`Coadd3DPar` carry
  per-file lists).
- **P3: per-target parameter lists** `targets = A, B`, `ra = 1.0, 2.0`.
  The current ParSet machinery handles this with no new container, but it
  is error-prone for users.

I lean to **P2**. Values are per target, not per frame, and the file
stays readable. `target_matches` (already in `inputfiles.py`) gives a
robust name match. See Q4.

---

## 4. Evaluation of approaches A, B, C

Constraints from Goals:

- **K1:** HostSub_GP (+ JAX/tinygp) remains optional.
- **K2:** With hostsub disabled, PypeIt output is unchanged.
- **K3:** MultiSlit longslit only.
- Plus the fact that HostSub_GP is a **post-processor over a target's set
  of exposures + standard** and needs **network access** (cutouts).

| Criterion | A. In-pipeline step in `run_pypeit` | B. Afterburner at end of `run_pypeit` when `[hostsub]` set | C. New `pypeit_hostsub` script reading `HostSubPar` from `.pypeit` |
|---|---|---|---|
| Fits HostSub's data model (all exposures of a target + std, after extraction) | **Poor.** `reduce_exposure` works one exposure (or comb group) at a time. HostSub needs every exposure of the target plus the std spec1d, and shares a rectification grid across them. That forces restructuring of `reduce_all`. | **Good.** After the science loop in `reduce_all`, all spec1d/spec2d exist. | **Good.** It runs on finished outputs, like `coadd_2dspec` / `setup_datacube`. |
| K1 optional deps | Hard. The import guard sits inside the core loop. | OK. Lazy import inside the afterburner, with a clear error if `use_gp=True` and the package is missing. | **Best.** Only the script imports it. `run_pypeit` never does. |
| K2 output unchanged when off | Risky. It touches `reduce_exposure`/`extraction` code paths, and any refactor can change bits. | Good if gated by `use_gp`. The only shared change is the `PypeItPar` schema (new section, default off). | **Best.** `run_pypeit` code is untouched. Only `PypeItPar` gains a default-off section. |
| K3 scope check | Spread across pypelines. | One gate (`spectrograph.pypeline == 'MultiSlit'`). | One gate in the script. |
| Re-running only the host model (iterating on z, regions, GP params) | Requires re-reducing. | Requires re-running `run_pypeit` (calibs are reused, but science is re-reduced unless output is reused). This is the main UX cost. | **Cheap.** Re-run just the script, and the `_rect` cache is reused. This matches how users actually iterate on HostSub. |
| Network / failure isolation | A cutout-download failure would kill the reduction. | It would fail at the end of a long run. It can be made non-fatal, but then it fails silently. | Isolated. A failure leaves the reduction intact. |
| Coadd2d hand-off | n/a | Could call coadd2d afterwards, but that duplicates `coadd_2dspec`. | Can write the `.coadd2d` file (like `pypeit_setup_coadd2d`) or call coadd2d directly. |
| Single config source (`.pypeit`) | Yes | Yes | Yes (reads the same `.pypeit`) |
| Precedent in PypeIt | None (no post-extraction multi-exposure step inside `run_pypeit`) | Weak. `[rdx] fluxcalib/coadd1d` flags exist but only drive *separate* scripts. | **Strong.** `coadd_1dspec`, `coadd_2dspec`, `sensfunc`, `flux_calib`, and `setup_datacube` (which already uses `outputfiles.existing_spec2d_files(pypeit_file, target, ...)`) |
| Implementation size | Large | Medium | Medium |
| Reviewability upstream (PR to pypeit) | Low | Medium | **High.** Self-contained script, ParSet, docs, tests. |

---

## 5. Recommendation

**C (a `pypeit_hostsub` script reading `[hostsub]` from the `.pypeit`
file). Optionally, a thin B later as a convenience that only calls C's
entry point.**

Reasons:

1. **It matches what HostSub_GP is.** It operates on a target's full
   set of extracted exposures plus a standard, after `run_pypeit`. C keeps
   it a post-processor, just as PypeIt already treats coadding, fluxing
   and datacube building.
2. **It is the strongest way to meet K1 and K2.** `run_pypeit` is not
   modified. The only change seen by users who do not set `use_gp` is a
   new, default-off `[hostsub]` section in `PypeItPar`. With `use_gp=False`
   the reduction is byte-identical; this can be tested by diffing a
   dev-suite spec2d before and after.
3. **It fits the real workflow.** Users iterate on `z`, regions, seeing
   and GP limits many times per target. Re-running a script on cached
   `_rect` files is minutes. Re-entering `run_pypeit` (B) is not.
4. **It isolates network and JAX failures** from the reduction.
5. **One config source is still achieved.** The `.pypeit` file holds both
   the reduction and the `[hostsub]` block. `hostsub_setup` and
   `hostsub.txt` disappear. Target/std/object association uses
   `PypeItMetaData`, `inputfiles.matching_science_rows`, and
   `outputfiles.existing_spec2d_files` instead of filename substrings.

Implied design directions (to be detailed in `hostsub_design.md`):

- Keep the PypeIt-facing layer (loading spec2d/spec1d, choosing the
  trace, deriving pixel scale/PA/slit width from the `Spectrograph` and
  `fitstbl`, writing back the skymodel, extracting) **in PypeIt**. Call
  HostSub_GP only for rectification → `SpecModel` → GP. This avoids the
  broken rows in §2 and the hard-coded instrument list.
- Write products to `Science/` (or a `Science_hostsub/` sibling; see Q6)
  as `DataContainer`-backed files, with QA in the standard `QA/` tree.
- Optionally write the `.coadd2d` file for the user, or run coadd2d
  directly (Q7).

---

## Clarifications

Please answer inline.

1. **Local environment.** Which Python environment should I use for later
   tasks? The default `python` (miniconda base, 3.14) has neither an
   importable `pypeit` nor JAX/tinygp, and HostSub_GP needs Python
   >= 3.12. Should I create a dedicated conda env with
   `pip install -e ".[dev]"` plus HostSub_GP, or is there an existing
   one? `PYPEIT_DEV` is unset. Where is the dev-suite `RAW_DATA`?
2. **Upstream HostSub_GP.** Can we change HostSub_GP itself (e.g. PR fixes
   for `BOX_RADIUS`, `extract_boxcar`, `fit_profile`; expose a
   pipeline-agnostic API that takes arrays instead of PypeIt filenames)?
   Or must PypeIt treat it as a frozen dependency? Is Chang Liu in the
   loop?
3. **Where the PypeIt-facing layer lives.** Do you agree with moving
   loading, trace selection, geometry and write-back into PypeIt and
   calling HostSub_GP only for modeling? Or should we call
   `SpecData.from_pypeit` / `update_pypeit_skymodel` as-is and accept the
   four-instrument limit?
4. **Per-target parameters.** P1 (data-block columns), P2 (named
   `[[targets]][[[NAME]]]` subsections), or P3 (parallel lists)? I lean
   to P2.
5. **Parameter naming.** Is the master switch `use_gp` (as in the design
   prompt), or `[hostsub] process = True`, or something else? Should
   `--skip_gp` / "classic" mode be exposed at all?
6. **Outputs.** Should `_hostsub` spec2d files go in `Science/` (cookbook)
   or `Science_hostsub/` (code)? Should the intermediate `_rect`/`_preproc`
   caches be kept as user-visible products?
7. **Coadd step.** Should `pypeit_hostsub` (i) stop after writing
   `_hostsub` spec2d files, (ii) also write a ready-to-run `.coadd2d` file,
   or (iii) run PypeIt's coadd2d itself? Is HostSub's own `--coadd2d` mode
   (model on the coadd, then its own spec1d writer) in scope at all?
8. **Per-exposure rewrite caveat.** In non-coadd mode, HostSub rewrites
   every exposure with the *last* exposure's model. Is that intended
   science behaviour we must preserve, or a bug we should fix (a per-exposure
   model, or a model fit on the coadd and projected back to each frame)?
9. **Recommended reduction settings.** Should `pypeit_hostsub` check
   (warn) or enforce `tracethresh` and `no_local_sky`? Or should
   `use_gp = True` set them automatically, which would change
   `run_pypeit` output only when hostsub is on?
10. **Network.** Archival cutouts and astrometry.net need internet. Is
    a pre-downloaded cutout directory (offline mode) required, e.g. for
    the dev suite and CI?
11. **Instruments for validation.** Which spectrograph and dataset in
    `RAW_DATA` should be the reference (LRIS, Binospec longslit, ALFOSC,
    FORS2)? Binospec on this fork is now `Fiber` for the IFU. Is
    Binospec *longslit* still MultiSlit and in scope?
12. **Base branch.** This fork is on `develop` (2.0.1-1356). Should the
    integration target `develop`? That implies fixing rows 11–12 of §2,
    or not using that HostSub path.

---

## Changelog

- **0.1 (2026-10-18), LK and Claude:** First draft. It covers the current
  workflow, the PypeIt API audit against develop, the parameter mapping,
  the A/B/C evaluation, and the recommendation for C.
