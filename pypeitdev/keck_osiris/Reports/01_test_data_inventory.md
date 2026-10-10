# 01 - Keck/OSIRIS test-data inventory and header verification

Date: 2026-10-10.  Author: Claude (Fable 5.1) with JXP.
Task: Design task 2 of `claude_prompts/keck_osiris.md` (act on Q&A answers
1-10; download the OsirisDRP test data; verify header assumptions).

## 1. What was downloaded

Source: `OsirisDRP/tests/drptestbones/map_file_urls.txt` (Dropbox links).
Files are untrusted downloads and were placed in their own directories.
`pypeit_setup` globs a raw directory non-recursively, so the `calib/` and
`ref/` sub-directories are invisible to it.

| Path under `RAW_DATA/keck_osiris/` | Size | Content | Status |
|---|---|---|---|
| `Hn3_100/s160321_a002010.fits` | 37.8 MB | 3C 273, Hn3, 0.100", 240 s, 2016-03-21 (H2RG) | OK |
| `Hn3_100/s160321_a002011.fits` | 37.8 MB | same setup, next frame; quasar absent (sky/off frame) | OK |
| `Hn3_100/calib/s160318_c005___infl_Hn3_100.fits` | 159 MB | rectification matrix, Hn3/100, 2016-03-18 | OK |
| `Hn3_100/ref/s160321_a002010_Hn3_100_ref.fits` | 12.5 MB | DRP reference cube (411 x 66 x 51) | OK |
| `Kn5_035/s150531_a025001.fits` | 37.8 MB | V5558 Sgr, Kn5, 0.035", 600 s, 2015-05-31 (Hawaii-2), ISSKY=0 | OK |
| `Kn5_035/s150531_a025002.fits` | 37.8 MB | same, ISSKY=1 | OK |
| `Kn5_035/ref/s150531_a025002_Kn5_035_ref.fits` | 14.1 MB | DRP reference cube (465 x 66 x 51) | OK |
| `Kn5_035/calib/s150901_c008___infl_Kn5_035.fits` | - | rectification matrix, Kn5/035 | **unavailable**: both the Dropbox link and the tkserver mirror return an HTML page |

The 2017 Kbb/035 matrix referenced in `OsirisDRP/data/calibrations.xml`
lives on a Google Drive folder ("Test Data/Trace/") with no public link
in the repo.

## 2. Raw file structure (both eras)

| Item | 2015 (Hawaii-2) | 2016 (H2RG) |
|---|---|---|
| ext0 | PRIMARY, float32 2048x2048, DN/s (BUNIT absent) | `SCI`, float32, `BUNIT='ADU per second per coadd'` |
| ext1 | float32 2048x2048, no EXTNAME (DRP "IntFrame" = noise) | `VAR`, `BUNIT='ADU squared per coadd'` |
| ext2 | uint8, values {0, 9}; ~1400 pixels = 0 | `DQ`, `BUNIT='bitmap'`; **all pixels = 9** (no bad-pixel map applied) |
| DQ bit meaning (ext2 COMMENTs, 2016) | - | bit0 bad pixel (map), bit1 raw saturation, bit2 read-delta saturation, bit3 max delta dropped in UTR, bit4 min delta dropped |
| ext1/ext2 headers | 7 cards only | full copies of the primary header plus `DATASEC='[1:2048,1:2048]'`, `ITIME0`, `COADDS0`, `READS0`, `GROUPS0` |

Note on units: the "9 = good" convention (bits 0 and 3 set) is what the
DRP tests for; PypeIt should treat `DQ != 9` as bad.  The VAR values do
not follow a simple Poisson/(t*gain) relation (ratio VAR/SCI ~ 7e-4 for
bright pixels vs 1.9e-3 expected for t = 240 s, g = 2.15), so VAR should
be treated as an opaque weight, not as PypeIt's variance.

## 3. Header cards: manual vs reality

Cards checked on all four raw frames.  "trunc" = `s160321_a002010`,
whose primary header has only 259 cards (vs 384 in the next frame); the
DRP reference cube records the warning "Some crucial keywords missing
from header".  This is the OSIRIS "global server hiccup" the DRP cookbook
warns about and PypeIt must tolerate it.

| Card | 2015 value | 2016 value | 2016 trunc | Comment |
|---|---|---|---|---|
| `INSTRUME` | `'OSIRIS'` | **absent** | absent | use `CURRINST` (= `'OSIRIS'` in all) as fallback |
| `TELESCOP` | `''` | `'Keck I'` | `'Keck I'` | not reliable for identification |
| `INSTR` | `'spec'` | absent | absent | spec vs imager flag only in 2015 |
| `OBSTYPE` | `'astro'` | `'a'` | `'a'` | match on first letter |
| `ISSKY` | 0 / 1 | 0 | **absent** | optional, default 0 |
| `SFILTER` | `'Kn5'` | `'Hn3'` | `'Hn3'` | dispname |
| `SSCALE` | `'0.035'` (str) | `'0.100'` (str) | `'0.100'` | decker; string, parse to float |
| `SS1NAME`/`SS2NAME`/`SFWNAME`/`SLMNAME` | present | present | absent | redundant with SFILTER/SSCALE |
| `ITIME` | 600.0 (**s**) | 240000 (**ms**) | 240000 | units changed with the 2016 upgrade |
| `TRUITIME` | absent | absent | absent | manual lists it; not in these files |
| `RDITIME` | 19.354 (s between reads) | 238.9954 (~ total) | 238.9954 | meaning changed too |
| `COADDS` | 1 | 1 | 1 | |
| `NUMREADS` / `SAMPMODE` | 32 / 1 (UTR) | 162 / 3 (MCDS) | 162 / 3 | sets the read noise |
| `GAIN` | 0.23 | absent | absent | 2015 value is not e-/ADU; ignore, hard-code 2.15 |
| `SATURATE` | 33000 | 65535 | 65535 | |
| `DATE-OBS`, `UTC`, `MJD-OBS` | present | present | present | MJD-OBS is the clean choice |
| `RA`, `DEC` | deg (float) | deg (float) | deg | `EQUINOX = 2000` |
| `AIRMASS`, `PARANG`, `EL` | present | present | present | |
| `ROTPOSN`, `INSTANGL`, `ROTREFAN` | 312, 312, 0 | 312, 312, 0 | 312, 312, 0 | |
| `PA_SPEC` | 0.0 | 0 | **absent** | used by DRP for WCS; derive from ROTPOSN/INSTANGL instead |
| `PONAME` | `'ospec'` | `'ospec'` | `'ospec'` | |
| `RAOFF`, `DECOFF` | ~0 | ~0 | 0 | tiny; not a usable dither record |
| `WXPRESS`, `WXOUTTMP`, `WXOUTHUM` | 618.3, 3.0, 8.4 | **absent** | absent | needed by the IFU subheader; need defaults |
| `DTMP7` (TMA temperature, used by DRP wavelength scaling) | 59.2 | 55.3 | **absent** | |
| `DATAFILE` | `'s150531_a025001'` | `'s160321_a002010.fits'` | same | |
| `FRAMENUM` | 1 | 10 | 10 | |

Implications for the design (feeding the Q&A):

- `exptime` must be `ITIME * COADDS`, with ITIME divided by 1000 when the
  value is in ms (detect by `ITIME > 2000` or by MJD >= 57388).
- `ISSKY`, `PA_SPEC`, `DTMP7`, `WX*` and `INSTRUME` must all be optional
  with defaults; a truncated header must not abort a reduction.
- The two 2016 frames are both `ISSKY = 0`, yet frame 011 contains no
  quasar (peak 0.13 vs 80.8 DN/s in frame 010) and the DRP uses it as the
  "Subtract Frame" calibration.  Automatic science/sky typing from the
  header is therefore unreliable; the pairing must come from the PypeIt
  file (`bkg_id`), which PypeIt already supports.
- The 2015 pair is the reverse: the DRP reduces the `ISSKY = 1` frame and
  subtracts the `ISSKY = 0` frame (test_emission_line looks at OH lines).

## 4. Geometry measured on the data

| Quantity | Measured | Manual |
|---|---|---|
| Spectral PSF across the slice (rect matrix, Hn3/100, slice 600) | FWHM 2 px, peak fraction 0.37-0.39 | ~2 px |
| Spectra tilt along dispersion (rect matrix centroid, 1412 px) | 0.00 deg | spectra aligned to rows |
| Horizontal stagger between spectra 2 rows apart (raw frame cross-correlation, rows 1010 vs 1012) | 29 px | "should be 32, ~29 because of anamorphism" |
| Spacing of illuminated spectra at one column (point source, 100 mas) | 32 px (rows 842, 870, 902, 933, 965, 997, 1030) | 64 lenslet columns x 32 px = 2048 |
| Illuminated columns for one slice (Hn3, narrowband) | 434-1846 (1413 px) | 3 x ~390 px narrowband spectra head-to-tail |
| Column sum of one slice profile at one wavelength | 0.79 | DRP divides the solution by 1.28 (= 1/0.78) |

So the raw frame is a comb of horizontal spectra, one every 2 rows with a
2-px FWHM, each shifted 29 px along the dispersion from its vertical
neighbour, with the pattern repeating every 32 rows (16 spectra per
band, 64 bands).  This confirms the Q1 decision: the Keck PSF model
(rectification matrix) is required to separate them.

## 5. Rectification matrix format (verified)

`s160318_c005___infl_Hn3_100.fits`:

| Ext | Shape / dtype | Meaning (from `spatrectif_000.pro` and the C code) |
|---|---|---|
| 0 | (1216, 2) int16 | `hilo`: bottom and top detector row of the 16-row window of each slice; slice 0 = rows 2015-2029, then -32 rows per slice within a lenslet column |
| 1 | (1216,) int16 | `effective`: 1 if the slice is illuminated (1017 of 1216 for this matrix) |
| 2 | (1216, 16, 2048) float32 | influence (PSF) cube: for slice s, column x, the fraction of that lenslet's flux landing in rows hilo[s,0] + 0..15 |
| header | `BASESIZE=15`, `WTLIMIT=0.01`, `SLICE=14`, `SHIFT=-26` | parameters of `mkrecmatrx` |

Slice index to lenslet: slice = 64 * column + row (19 columns x 64 rows
for broadband; the DRP's `assembcube` maps the 1216 narrowband slices to
three lenslet-column blocks offset by 0, +16 and +32 columns, plus a few
extra spectra for `sp > 831`).  The Kbb mapping is tabulated in
`OsirisDRP/testing_scripts/trace_spectra/kbb_2016_slice_to_lenslet.txt`.

## 6. The DRP extraction algorithm (for the design document)

From `spatrectif_000.c` (366 lines).  Let `F[j,x]` be the raw frame,
`B[s,l,x]` the influence cube, `b[s]` the bottom row of slice `s`, and
`c[s,x]` the unknown lenslet spectra.  The forward model is
`M[j,x] = sum_s B[s, j-b[s], x] * c[s,x]`, i.e. each detector column `x`
is an independent linear system of 2048 pixels in 1216 unknowns.  The
solver, per column, is an iterative "blame" redistribution:

1. Precompute two normalised kernels per (s, x): `blame = B^3 / sum(B^3)`
   (sharpened, early iterations) and `fblame = B / sum(B)` (final).
2. Start with `c = 0`.  For each iteration (default 40 for science,
   `relaxation` parameter):
   - residual `r = F - M` on pixels with `DQ == 9`, else 0;
   - iterations 0-14: accumulate `t[s] = relax * sum_l r[b+l] * blame`
     and apply it, with a +/-0.4 (iters < 10) then -0.2 (iters 10-14)
     coupling to the slices +/-64 (the spatial neighbours);
   - iterations >= 15: `c += relax * sum_l r[b+l] * fblame`;
   - iterations 12-17: smooth `c` along the dispersion with a
     [0.1, 0.2, 0.4, 0.2, 0.1] kernel.
3. Divide by 1.28; propagate noise as `sum_l fblame * N`; mark the output
   pixel good if more than half the kernel weight falls on `DQ == 9`
   pixels.

This is a damped Landweber / Gauss-Seidel iteration with a hand-tuned
schedule.  A numpy port can vectorise over all columns at once with a
banded sparse matrix per column (each column has <= 16 x 1216 nonzeros)
or simply reproduce the loops with array ops over (s, x).  A least-squares
alternative (`scipy.sparse.linalg.lsqr` per column with non-negativity)
would be a cleaner formulation and can be validated against the DRP
reference cubes we now have.

## 7. Open items raised by this inventory (see Q&A 11-17)

1. No Kn5/035 matrix is obtainable, so the 2015 set cannot be rectified.
2. Header tolerance rules (ITIME units, truncated headers, missing
   weather cards) need agreement.
3. VAR/DQ semantics: whether to use VAR at all.
4. Science/sky pairing must be manual (`bkg_id`) for these data.
5. A science-mode (Kbb/Kn3, 35/50 mas) H2RG dataset with its matrix is
   still needed; KOA is the obvious source.
