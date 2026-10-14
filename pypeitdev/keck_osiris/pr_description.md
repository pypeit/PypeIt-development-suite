# Pull requests: Keck/OSIRIS support

Two pull requests, to be opened together and cross-referenced.

---

## PypeIt: Add Keck/OSIRIS (lenslet IFU) support

**Base:** `develop`  **Branch:** `osiris_port`
**Dev-suite PR:** pypeit/PypeIt-development-suite#XXX (data, pypeit files,
dev-suite registrations and vet tests; must be merged with this one)

### Summary

Adds `keck_osiris`, the Keck I adaptive-optics near-infrared lenslet
integral-field spectrograph (H2RG detector era, 2016 onwards), as a
supported spectrograph.  OSIRIS records up to ~3000 interleaved lenslet
spectra only two pixels apart, so the spectra cannot be traced; as in the
Keck OSIRIS data reduction pipeline, the raw frames are *rectified* with the
Keck rectification matrices (the measured detector footprint of every
lenslet), and the resulting lenslet image (one spatial pixel per lenslet) is
reduced by the `Fiber` pypeline with one spectrum per lenslet.

- New module `pypeit/core/lenslet_rectify.py`: matrix I/O (`.fits` and
  `.fits.gz`), matrix selection by Keck validity epoch (shipped table
  `rectmat_epochs.ecsv`, from http://tkserver.keck.hawaii.edu/osiris/), a
  vectorised port of the DRP solver (`spatrectif_000.c`) and an exact
  block-tridiagonal least-squares solver (default) with flux-independent
  weights and exact variance propagation, the narrowband lenslet layout,
  the Keck wavelength relation and the Keck-derived tilts.
- New spectrograph `pypeit/spectrographs/keck_osiris.py`: metadata for all
  header eras (including the truncated 2016 headers), frame typing
  (standards from the target name), DN conversion, static bad-pixel masks,
  the rectification hook, block-slits, lenslet identities
  (`MASKDEF_ID = row * 51 + col`), Keck-derived tilts and wavelengths
  (vacuum, epoch-matched coefficients shipped in `pypeit/data`), cube setup,
  and the defaults (A-B without darks, per-lenslet OH refinement, Davies OH
  scaling, second cosmic-ray pass, IR telluric settings).
- Sky: per-lenslet OH-line wavelength refinement of every science frame;
  frames with no measurable sky line (1.5-s standards) adopt the shifts
  measured on the arc frame of their calibration group.  Davies (2007) OH
  scaling of A-B frames, with each group factor kept only if it lowers the
  integrated line residual.  A second cosmic-ray pass on the A-B lenslet
  spectra catches the hits L.A.Cosmic misses on single raw frames.
- New generic script `pypeit_fiber_cube` (datacubes and point-source
  spectra from `Fiber`-pypeline spec1d files); `pypeit_binospec_ifu_cube`
  becomes a thin alias.  Development script `pypeit_osiris_rectify`.
- Documentation: `doc/spectrographs/keck_osiris.rst` (data requirements,
  obtaining matrices and the epoch rule, the `bkg_id` convention, darks,
  parameters, cubes, telluric recipe, limitations), the tutorial
  `doc/tutorials/keck_osiris_howto.rst`, release notes, regenerated
  parameter docs and spectrograph table.

### Shared-code changes for the reviewers

All are generic and inactive unless a spectrograph hook or a parameter
enables them; existing unit tests pass unchanged.  Please review in
particular:

*From M3 (decision M3-7):*

1. **`runtime_par` mechanism**: `Spectrograph.set_runtime_par(par)`, called
   by `PypeIt.__init__`, gives a spectrograph the run parameters (OSIRIS
   needs the matrix directory and solver inside `RawImage.process`, which
   only receives the processing parameters).
2. **Column-window boxcar speed-up in `FiberExtract`**: the boxcar of a
   fiber only touches the columns of its aperture (`_fiber_columns`)
   instead of the full image; identical results, needed for ~3000 fibers.
3. **`rectify` step of `RawImage.process`** (runs when
   `Spectrograph.has_rectification` is True): replaces the image, the
   variance components and the BPM by the spectrograph's transform
   (`Spectrograph.rectify`).
4. **Fixed-calibration hooks**: `Spectrograph.get_fixed_tilts`
   (`Calibrations.get_tilts` uses it when overridden) and
   `get_fixed_wavecalib` with the new wavelength `method = 'fixed'`;
   per-fiber wavelengths `get_fiber_wavelengths` applied as offsets by the
   Fiber pypeline.

*From M4 (Report 06, section 2):*

5. `raw_header` argument of `get_fixed_tilts` (header of a raw frame of the
   calibration group, passed by `Calibrations.get_tilts`).
6. The un-subtracted frame of an A-B reduction (`bkg_redux_sciimg`) is
   passed to `FindObjects` (`pypeit_steps`, `exposure`, `reduce_by_step`).
7. `[calibrations][wavelengths] fiber_sky_refine` (wvutils +
   `FiberFindObjects._refine_fiber_sky_lines`) and
   `[reduce][skysub] osiris_oh_scaling` (`skysub.oh_scaled_sky` +
   `FiberFindObjects._fiber_oh_scaling`).
8. `rectify=False` for bias and dark frames (`buildimage`/`rawimage`), so
   darks are subtracted from raw frames before their rectification.
9. `Spectrograph.fiber_cube_setup` hook, `wcs_ref` and
   `extract_fiber_point_source` in `core/datacube.py`, and
   `pypeit_fiber_cube` with `pypeit_binospec_ifu_cube` as an alias.

*From M5 (this milestone):*

10. `skysub.oh_scaled_sky` gains an integrated-residual guard (`guard`,
    `guard_lines`; helper `skysub.oh_line_residuals`) and the new parameter
    `[reduce][skysub] osiris_oh_min_exptime` (default None).
11. `FiberFindObjects._arc_fiber_sky_shifts`: arc-frame fallback of the
    sky-line refinement (reads the `Arc_<calib_key>.fits` calibration image
    next to the `WaveCalib`; cached per file).
12. `skysub.fiber_cosmic_rays` + `FiberFindObjects._fiber_residual_crs`,
    parameter `[reduce][skysub] fiber_cr_nsig` (default None; 8 for
    OSIRIS): adds CR flags to the A-B science image.
13. `wvutils.fiber_sky_line_shifts` silences the per-line fit warnings
    (thousands per frame).

### Testing

- `pytest pypeit/tests --ignore=pypeit/tests/test_pkgdata.py`: 908 passed
  (`test_pkgdata.py` downloads from a GitHub branch that does not exist and
  deletes a tracked file).
- Dev suite: `./pypeit_test reduce afterburn vet -s keck_osiris/Hn3_100
  keck_osiris/Kbb_050 keck_osiris/Kn3_035`: 7/7 reduce and afterburn tests
  PASSED (Hn3 3.6 min, Kn3 5.5 min, Kbb 38 min; fiber cubes, sensfunc,
  telluric); the 21 OSIRIS vet tests (`test_keck_osiris_rectify.py`,
  `test_keck_osiris_flow.py`, `test_keck_osiris_kbb.py`,
  `test_keck_osiris_kn3.py`) pass; the `keck_osiris` load-image unit test
  passes.  (The `vet` phase runs the whole `vet_tests` directory; tests of
  other instruments fail for lack of their reductions in that output
  directory.)
- `cd doc; make htmlonly`: no warnings (the two warnings of the previous
  tree, from M3/M4 docstrings, are fixed).

### Credits

The rectification algorithm, the lenslet geometry, the wavelength
coefficients, the OH line list and the Davies OH scaling follow the Keck
OSIRIS data reduction pipeline (https://github.com/Keck-DataReductionPipelines/OsirisDRP;
Lyke et al. 2017, ascl:1710.021; Lockhart et al. 2019, AJ 157, 75; OH
scaling: Davies 2007, MNRAS 375, 1099).  The rectification matrices are
provided by W. M. Keck Observatory.  Developed by J. X. Prochaska in
collaboration with Claude (Anthropic's AI assistant).

---

## PypeIt-development-suite: Keck/OSIRIS setups and tests

**PypeIt PR:** pypeit/PypeIt#YYYY (must be merged with this one)

- `RAW_DATA/keck_osiris/` (Google Drive; see the staging list in
  `pypeitdev/keck_osiris/Reports/07_m5_release.md`): `Hn3_100` (2016 Hn3
  pair, matrix and DRP reference cube), `Kbb_050` (KOA 2024-02-03: NGC 3245,
  two A0V standards, darks; 73 frames for now, to be trimmed; matrix),
  `Kn3_035` (KOA 2024-08-13: nova V1405 Cas, two A/B pairs; matrix), and
  `Kn5_035` (2015 header-era frames for the error test).
- Pypeit files `keck_osiris_hn3_100.pypeit`, `keck_osiris_kbb_050.pypeit`,
  `keck_osiris_kn3_035.pypeit`; sensfunc and telluric inputs
  `keck_osiris_kbb_050.sens`/`.tell`.
- `test_scripts/setups.py`, `test_setups.py`, `pypeit_tests.py`: the three
  setups for `reduce`; `pypeit_fiber_cube` (new `PypeItFiberCubeTest`),
  `sensfunc` and `telluric` afterburn tests for `Kbb_050`.
- Vet tests: `test_keck_osiris_rectify.py` (rectification versus the DRP
  cube), `test_keck_osiris_flow.py` (OH-line tilt test on Hn3 and Kbb, OH
  refinement repeatability, quasar position), `test_keck_osiris_kbb.py`
  (spec1d, cube, telluric flatness), `test_keck_osiris_kn3.py` (narrowband
  layout, nova); `unit_tests/test_load_images.py` entry.
- Development record: design document, implementation plan, milestone
  prompts and Reports 01-07 in `pypeitdev/keck_osiris/` and
  `claude_prompts/keck_osiris/`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
