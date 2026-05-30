---
name: new-spectrograph
description: Add a new spectrograph (or new mode/arm of an existing one) to PypeIt — scaffold and register a Spectrograph subclass, define metadata/detectors/frame-typing/default parameters, and wire up the required dev-suite test. Use when the user wants to support a new instrument, telescope, grating mode, or detector in PypeIt.
---

# Add a new spectrograph to PypeIt

This skill lives in the `PypeIt-development-suite/` repo but operates on the
sibling `PypeIt/` repo (the main code) and on this repo (for the required
tests). Confirm both repos are present as sibling clones before starting.

## Reference material

- Authoritative guide: `PypeIt/doc/dev/new_spectrograph.rst`
- Existing classes to copy from: `PypeIt/pypeit/spectrographs/` (47 instruments).
  Start from the closest analog — same pipeline path and similar detector/header
  format. Good templates:
  - MultiSlit with arms: `keck_lris.py` (shared base class + per-arm subclasses)
  - MultiSlit single detector: `keck_deimos.py` (most-documented methods)
  - Echelle (fixed-format): `vlt_xshooter.py`
  - Near-IR (calibrations partly disabled): `gemini_gnirs.py`
  - SlicerIFU: `keck_kcwi.py`
- Base class & required abstractions: `pypeit/spectrographs/spectrograph.py`
- Example PR (SOAR/Goodman): https://github.com/pypeit/PypeIt/pull/1179

## Procedure

Work in the `PypeIt/` repo on a branch off `develop`.

1. **Create the module** `pypeit/spectrographs/{telescope}_{spectrograph}.py`
   (e.g. `soar_goodman.py`). Subclass `Spectrograph` (or an existing base class
   if adding an arm/mode). Set the `name` attribute (None only for base classes)
   and the `telescope`. **`pypeline` must be `'MultiSlit'`, `'Echelle'`, or
   `'SlicerIFU'`.**

2. **Register the module** by adding its name to `__all__` in
   `pypeit/spectrographs/__init__.py`, keeping alphabetical order.

3. **Telescope**: add a `Telescope` object in `pypeit/telescopes.py` if the
   telescope is new, and add the telescope name to `valid_telescopes` in
   `pypeit/par/pypeitpar.py`.

4. **Implement the required methods** (compare to the template class):
   - `default_pypeit_par()` — instrument default parameters (see the
     `add-parameter` skill if you need new parameters).
   - `bpm()` — default bad-pixel mask.
   - `init_meta()` / `compound_meta()` — map raw FITS header cards to PypeIt
     metadata keys (see `doc/metadata.rst`).
   - `configuration_keys()` — metadata keys defining a unique configuration.
   - `check_frame_type()` — frame typing (arc, bias, dome, science, …) from
     metadata.
   - `get_detector_par()` — detector geometry/gain/readnoise, one per detector.
   - `get_rawimage()` — read the raw data (only override if the default in the
     base class is insufficient).
   - For fixed-format echelle, add the order-format methods (see
     `vlt_xshooter.py`).

5. **Wavelength solution**: most setups need an arc solution/template. Use the
   `wavelength-calibration` skill (`pypeit_identify`, then archive a template).

## Required tests (to be "supported")

Per `new_spectrograph.rst`, a supported spectrograph requires, in **this**
(`PypeIt-development-suite/`) repo:

- A full pipeline run for each grating/mode — add raw data and a `.pypeit` file
  and register the setup. Use the `add-devsuite-setup` skill.
- A load-images unit test entry in `unit_tests/test_load_images.py`.

## Docs to update (in `PypeIt/`)

- Note the addition in the latest release doc under `doc/releases/`
  ("Instrument-specific" section).
- Optionally add an instrument doc under `doc/spectrographs/` and link it from
  `doc/spectrographs/spectrographs.rst`.
- Add or point to a tutorial in `doc/tutorials/`.

## Verify

- `python -c "from pypeit.spectrographs.util import load_spectrograph; load_spectrograph('{name}')"`
- `pypeit_setup -s {name} -r <rawdir>` correctly types frames and identifies the
  configuration.
- Run the new dev-suite setup (see the `run-dev-suite` skill).
