# Report 03: First full reduction via `pypeit_test`

**Date:** 2026-08-25
**Author:** Claude (with JXP)
**Branch:** `int_ids` (PypeIt repo, off `develop`)
**Environment:** `pypeit14b` conda env; `PYPEIT_DEV` pointed at this repo
**Data:** EEV10 / R1200B / cenwave 4502 Å (arc, bias, flat, science)

## Summary

The first end-to-end reduction of INT/IDS data **passed the dev-suite
`reduce` test on the first attempt**: `PYPEIT DEVELOPMENT SUITE PASSED 1/1
TESTS` in 3m07s, producing calibrations, QA plots, and extracted 1D/2D
science spectra with a 0.23-pixel wavelength RMS.

## Dev-suite registration

The setup was registered following the standard conventions:

- Raw frames staged as `RAW_DATA/int_ids_eev10/INT_20201218/` (arc, bias,
  flat, science; copied from `pypeitdev/int_ids/INT_20201218/`). These are
  git-ignored and still need to be uploaded to the shared Google Drive.
- New reduction input `pypeit_files/int_ids_eev10_int_20201218.pypeit`
  (filename lower-cased per the harness convention; the `path` line is
  rewritten at run time).
- `test_scripts/setups.py`: added `'int_ids_eev10': ['INT_20201218']` to
  `all_setups`, which auto-derives the default `reduce` coverage.
- `unit_tests/test_load_images.py`: new `test_load_int_ids` (verified:
  1 passed).

## Command and result

```bash
cd $PYPEIT_DEV
PATH=<pypeit14b>/bin:$PATH ./pypeit_test reduce -s int_ids_eev10/INT_20201218
```

```
 1 passed/ 0 failed/ 0 skipped PASSED  int_ids_eev10/INT_20201218 pypeit
--- PYPEIT DEVELOPMENT SUITE PASSED 1/1 TESTS  ---
Total Time: 0:03:07  |  Disk usage: 0.157 GiB
```

## How far it went: all the way

Every reduction stage completed:

| Stage | Outcome |
|---|---|
| Bias | `Bias_A_0_DET01.fits` built from the single bias frame |
| Slit edges | 1 slit traced and synced across the 366-pixel window (`bound_detector=True`) |
| Wavelengths | `full_template` with `int_ids_R1200B.fits`: **RMS = 0.226 pixels** (~0.11 Å) |
| Tilts | RMS = 0.066 pixels (RMS/FWHM = 0.023) |
| Flat field | Pixel flat + spatial illumination built (illumflat falls back to the pixelflat bspline — expected with a single tungsten flat) |
| Object finding | 1 object found and traced at spatial pixel 172.9 in both passes |
| Sky subtraction + extraction | Boxcar and optimal extraction completed |
| QA | Full HTML + PNG QA generated (arc fits, tilts, object profile/trace) |

Extracted-spectrum sanity checks (`spec1d` file):

- Wavelength coverage 3487–5412 Å, matching the R1200B template range.
- Object FWHM 1.88 pixels ≈ 0.75" — a plausible seeing value, and an
  independent confirmation of the 0.40"/pix EEV10 plate scale.
- Median S/N ≈ 3.7 per pixel (faint Gaia source, 2500 s) — consistent with
  expectations for this target.
- The object sits at spatial pixel 172.9, matching the aperture
  (168–178) of the original IRAF reduction of this data set.

## Warnings (all benign)

- `No raw dark/scattlight frames found ... Continuing` — none were taken;
  darks are disabled in the instrument defaults.
- `illumflat has no spatial bspline fit - using the pixelflat` — only one
  tungsten flat exists, used for both pixelflat and illumflat.
- A handful of numpy `RankWarning`/`RuntimeWarning` and pydantic serializer
  notices from deep in the fitting stack — none affect the products.

There were **no errors** and no instrument-specific warnings.

## Remaining items

1. **Google Drive**: upload `RAW_DATA/int_ids_eev10/INT_20201218/` so other
   developers (and Nautilus runs) can execute the test.
2. **`supported = True`** for `int_ids_eev10` now that the dev-suite
   reduction passes (awaiting go-ahead — Q&A #10).
3. **Docs**: release-notes entry and a `doc/spectrographs/int_ids.rst` page.
4. **Afterburn coverage**: no standard star was observed in this data set,
   so fluxing/sensfunc tests are not possible yet; a future data set with a
   standard would also let us validate the 120 s science/standard split.
