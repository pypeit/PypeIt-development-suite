# Report 02: `pypeit_setup` on the INT/IDS test data

**Date:** 2026-08-24
**Author:** Claude (with JXP)
**Branch:** `int_ids` (PypeIt repo, off `develop`)
**Environment:** `pypeit14b` conda env (editable install of this checkout)
**Data:** `pypeitdev/int_ids/INT_20201218/` (arc, bias, flat, science;
EEV10 detector, R1200B grating, cenwave 4502 Å, 1×1 binning)

## Summary

`pypeit_setup` runs **cleanly end-to-end** on the test data set with the new
`int_ids_eev10` spectrograph: all four frames are ingested, correctly
frame-typed, grouped into a single configuration (Setup A), and a vetted
`.pypeit` file is written. The log contains **zero warnings or errors**.

## Commands run

```bash
# First pass: inspect the automated grouping (.sorted file)
pypeit_setup -s int_ids_eev10 \
    -r $PYPEIT_DEV/pypeitdev/int_ids/INT_20201218/

# Second pass: write the .pypeit file for configuration A
pypeit_setup -s int_ids_eev10 \
    -r $PYPEIT_DEV/pypeitdev/int_ids/INT_20201218/ -c A
```

Both `int_ids_eev10` and `int_ids_redplus2` appear in the script's list of
valid spectrograph identifiers.

## Results

One unique configuration was found:

```
Setup A
    dispname: R1200B
     cenwave: 4502.182
    detector: EEV10
     binning: 1,1
```

The data table (from the generated `int_ids_eev10_A.pypeit`):

| filename | frametype | ra | dec | exptime | airmass | lampstat01 |
|---|---|---|---|---|---|---|
| arc.fits | arc,tilt | 356.6472 | -4.3084 | 90.01 | 1.672 | CuAr CuAr CuNe CuN |
| bias.fits | bias | 0.0 | 0.0 | 0.0 | 1.012 | off |
| flat.fits | pixelflat,illumflat,trace | 0.0 | 0.0 | 60.01 | 1.515 | W |
| science.fits | science | 356.6472 | -4.3084 | 2500.01 | 1.355 | off |

Everything is as expected:

- Frame typing (via `IMAGETYP` + `AGARCLMP`) is correct for all four frames.
- RA/DEC are correctly converted from the sexagesimal header strings to
  decimal degrees (bias/flat report the parked-telescope 0/0 values, which
  is what their headers contain).
- The configuration keys (`dispname`, `cenwave`, `detector`, `binning`)
  produce a single setup, and the `calib` group assignment is uniform.
- A `.calib` association file is written and the `.pypeit` file passes
  PypeIt's vetting step.

## Airmass check

The IDS headers provide `AIRMASS` ("effective mean airmass") plus
`AMSTART`/`AMEND`. As a cross-check, the airmass was recomputed with astropy
(`SkyCoord` → `AltAz` at the Roque de los Muchachos `EarthLocation`) for the
science frame:

| epoch | astropy sec(z) | header |
|---|---|---|
| exposure start (`MJD-OBS`) | 1.2938 | `AMSTART` = 1.2930 |
| mid-exposure (start + EXPTIME/2) | 1.3528 | `AIRMASS` = 1.3549 |

Agreement is better than 0.2%, confirming both the header values and our
site/coordinate/time metadata handling.

## Recommended changes

1. **Airmass robustness**: rather than trusting the `AIRMASS` card, compute
   the mid-exposure airmass from `ra`, `dec`, `mjd`, and the telescope
   `EarthLocation` with astropy (as validated above), falling back to the
   header card when coordinates are unavailable (e.g. bias/flat frames with
   RA=DEC=0). This logic is generic and belongs in a common module (e.g.
   `pypeit.core.meta` or similar) so any spectrograph with missing or
   unreliable airmass cards can reuse it as a `compound_meta` helper.

2. **Dev-suite registration** (needed before the instrument can be declared
   supported): stage the raw data as `RAW_DATA/int_ids_eev10/INT_20201218/`,
   add a `pypeit_files/int_ids_eev10_int_20201218.pypeit`, register the setup
   in `test_scripts/test_setups.py`, and add an entry to
   `unit_tests/test_load_images.py`.

3. **Documentation**: note the new instrument in the current release doc
   (`doc/releases/`), and add `doc/spectrographs/int_ids.rst` linked from
   `spectrographs.rst`.

4. **Frame-typing coverage**: the current logic handles the four frame types
   in this data set. Darks, sky flats, and standard-star frames have not been
   seen yet; the 120 s science/standard exposure-time split (Q&A #2) should
   be revisited when a standard-star observation is available.

5. **`supported = True`** for `int_ids_eev10` once the dev-suite reduction
   passes (prompt 6).

## Next steps

Generate the dev-suite `.pypeit` file and run the full reduction via
`pypeit_test` (prompt 6), which will exercise slit tracing on the windowed
366-pixel-wide detector, the R1200B `full_template` wavelength solution, and
extraction.
