# 02 - KOA science-mode dataset (2024-02-03, Kbb/050) and matrix request

Date: 2026-10-10.  Author: Claude (Fable 5.1) with JXP.
Task: Design task 4 of `claude_prompts/keck_osiris.md` (Q&A 22).

## 1. Why this night

The KOA TAP query of public OSIRIS spectrograph frames since 2016
(design document, Section 9) showed 2024-02-03 as a complete set on the
current detector in the most common science mode: NGC 3245 in Kbb at
0.050" with sky frames flagged, two A0V telluric standards and a full
dark sequence.

## 2. Download

- Source: KOA TAP (`koa_osiris` table) for the listing;
  `https://koa.ipac.caltech.edu/cgi-bin/getKOA/nph-getKOA?filehand=...`
  for the files (public data, no login).
- Destination: `RAW_DATA/keck_osiris/Kbb_050/` with the selected file
  list (`koa_file_list.txt`) and the full night listing
  (`koa_20240203_all_spec_frames.csv`) kept beside the data.
- Selection: all Kbb/050 frames of NGC 3245, HD101060 and HD89239, plus
  darks whose ITIME matches those frames (600, 10 and 2 s; there are
  no 45 s or 3 s darks).  73 of the night's 117 spectrograph frames;
  2.6 GB, all verified as valid FITS after download.
- A probe frame (`OS.20240203.56800.59.fits`, HD101060, 3 s) confirmed
  the 2024 header is complete: 526 cards, `INSTRUME='OSIRIS'`,
  `ISSKY=1` on the nod-off frames, `PA_SPEC=175`, `DTMP7=54.19`,
  `WXPRESS=617.57`, `BUNIT='ADU per second per coadd'`, `DQ` all 9.

## 3. Night listing

| KOAID | UT | Target | Filter | Scale | ITIME (s) | ISSKY | Fetched |
|---|---|---|---|---|---|---|---|
| OS.20240203.11652.80.fits | 03:14:12 | unknown | Drk | 0.020 | 600 | 0 | yes |
| OS.20240203.12309.99.fits | 03:25:09 | unknown | Drk | 0.020 | 600 | 0 | yes |
| OS.20240203.12966.89.fits | 03:36:06 | unknown | Drk | 0.020 | 600 | 0 | yes |
| OS.20240203.13624.32.fits | 03:47:04 | unknown | Drk | 0.020 | 600 | 0 | yes |
| OS.20240203.14281.72.fits | 03:58:01 | unknown | Drk | 0.020 | 600 | 0 | yes |
| OS.20240203.14938.60.fits | 04:08:58 | unknown | Drk | 0.020 | 2 | 0 | yes |
| OS.20240203.14952.37.fits | 04:09:12 | unknown | Drk | 0.020 | 2 | 0 | yes |
| OS.20240203.14967.16.fits | 04:09:27 | unknown | Drk | 0.020 | 2 | 0 | yes |
| OS.20240203.14980.42.fits | 04:09:40 | unknown | Drk | 0.020 | 2 | 0 | yes |
| OS.20240203.14993.68.fits | 04:09:53 | unknown | Drk | 0.020 | 2 | 0 | yes |
| OS.20240203.15006.94.fits | 04:10:06 | unknown | Drk | 0.020 | 5 | 0 |  |
| OS.20240203.15023.26.fits | 04:10:23 | unknown | Drk | 0.020 | 5 | 0 |  |
| OS.20240203.15041.11.fits | 04:10:41 | unknown | Drk | 0.020 | 5 | 0 |  |
| OS.20240203.15057.43.fits | 04:10:57 | unknown | Drk | 0.020 | 5 | 0 |  |
| OS.20240203.15074.77.fits | 04:11:14 | unknown | Drk | 0.020 | 5 | 0 |  |
| OS.20240203.15091.90.fits | 04:11:31 | unknown | Drk | 0.020 | 7 | 0 |  |
| OS.20240203.15110.47.fits | 04:11:50 | unknown | Drk | 0.020 | 7 | 0 |  |
| OS.20240203.15128.83.fits | 04:12:08 | unknown | Drk | 0.020 | 7 | 0 |  |
| OS.20240203.15147.70.fits | 04:12:27 | unknown | Drk | 0.020 | 7 | 0 |  |
| OS.20240203.15168.10.fits | 04:12:48 | unknown | Drk | 0.020 | 7 | 0 |  |
| OS.20240203.15186.97.fits | 04:13:06 | unknown | Drk | 0.020 | 10 | 0 | yes |
| OS.20240203.15212.47.fits | 04:13:32 | unknown | Drk | 0.020 | 10 | 0 | yes |
| OS.20240203.15238.99.fits | 04:13:58 | unknown | Drk | 0.020 | 10 | 0 | yes |
| OS.20240203.15265.51.fits | 04:14:25 | unknown | Drk | 0.020 | 10 | 0 | yes |
| OS.20240203.15292.30.fits | 04:14:52 | unknown | Drk | 0.020 | 10 | 0 | yes |
| OS.20240203.15318.55.fits | 04:15:18 | unknown | Drk | 0.020 | 15 | 0 |  |
| OS.20240203.15350.17.fits | 04:15:50 | unknown | Drk | 0.020 | 15 | 0 |  |
| OS.20240203.15381.79.fits | 04:16:21 | unknown | Drk | 0.020 | 15 | 0 |  |
| OS.20240203.15413.92.fits | 04:16:53 | unknown | Drk | 0.020 | 15 | 0 |  |
| OS.20240203.15445.30.fits | 04:17:25 | unknown | Drk | 0.020 | 15 | 0 |  |
| OS.20240203.15476.65.fits | 04:17:56 | unknown | Drk | 0.020 | 20 | 0 |  |
| OS.20240203.15513.37.fits | 04:18:33 | unknown | Drk | 0.020 | 20 | 0 |  |
| OS.20240203.15550.90.fits | 04:19:10 | unknown | Drk | 0.020 | 20 | 0 |  |
| OS.20240203.15584.77.fits | 04:19:44 | unknown | Drk | 0.020 | 20 | 0 |  |
| OS.20240203.15621.30.fits | 04:20:21 | unknown | Drk | 0.020 | 20 | 0 |  |
| OS.20240203.15656.74.fits | 04:20:56 | unknown | Drk | 0.020 | 25 | 0 |  |
| OS.20240203.15697.54.fits | 04:21:37 | unknown | Drk | 0.020 | 25 | 0 |  |
| OS.20240203.15737.32.fits | 04:22:17 | unknown | Drk | 0.020 | 25 | 0 |  |
| OS.20240203.15777.10.fits | 04:22:57 | unknown | Drk | 0.020 | 25 | 0 |  |
| OS.20240203.15816.89.fits | 04:23:36 | unknown | Drk | 0.020 | 25 | 0 |  |
| OS.20240203.15856.67.fits | 04:24:16 | unknown | Drk | 0.020 | 30 | 0 |  |
| OS.20240203.15903.80.fits | 04:25:03 | unknown | Drk | 0.020 | 30 | 0 |  |
| OS.20240203.15950.10.fits | 04:25:50 | unknown | Drk | 0.020 | 30 | 0 |  |
| OS.20240203.15995.94.fits | 04:26:35 | unknown | Drk | 0.020 | 30 | 0 |  |
| OS.20240203.16042.92.fits | 04:27:22 | unknown | Drk | 0.020 | 30 | 0 |  |
| OS.20240203.16088.83.fits | 04:28:08 | unknown | Drk | 0.020 | 120 | 0 |  |
| OS.20240203.16230.62.fits | 04:30:30 | unknown | Drk | 0.020 | 120 | 0 |  |
| OS.20240203.16373.93.fits | 04:32:53 | unknown | Drk | 0.020 | 120 | 0 |  |
| OS.20240203.16514.69.fits | 04:35:14 | unknown | Drk | 0.020 | 120 | 0 |  |
| OS.20240203.16656.47.fits | 04:37:36 | unknown | Drk | 0.020 | 120 | 0 |  |
| OS.20240203.16798.76.fits | 04:39:58 | unknown | Drk | 0.020 | 180 | 0 |  |
| OS.20240203.17014.49.fits | 04:43:34 | unknown | Drk | 0.020 | 180 | 0 |  |
| OS.20240203.17228.69.fits | 04:47:08 | unknown | Drk | 0.020 | 180 | 0 |  |
| OS.20240203.17443.92.fits | 04:50:43 | unknown | Drk | 0.020 | 180 | 0 |  |
| OS.20240203.17659.65.fits | 04:54:19 | unknown | Drk | 0.020 | 180 | 0 |  |
| OS.20240203.41442.14.fits | 11:30:42 | HD89239(A0V) | Kbb | 0.050 | 2 | 1 | yes |
| OS.20240203.41472.23.fits | 11:31:12 | HD89239(A0V) | Kbb | 0.050 | 2 | 0 | yes |
| OS.20240203.41606.87.fits | 11:33:26 | HD89239(A0V) | Kbb | 0.050 | 2 | 1 | yes |
| OS.20240203.41635.43.fits | 11:33:55 | HD89239(A0V) | Kbb | 0.050 | 2 | 0 | yes |
| OS.20240203.41743.55.fits | 11:35:43 | HD89239(A0V) | Kbb | 0.050 | 2 | 1 | yes |
| OS.20240203.41774.15.fits | 11:36:14 | HD89239(A0V) | Kbb | 0.050 | 2 | 0 | yes |
| OS.20240203.41802.20.fits | 11:36:42 | HD89239(A0V) | Kbb | 0.050 | 2 | 0 | yes |
| OS.20240203.42266.81.fits | 11:44:26 | NGC3245 | Kbb | 0.050 | 45 | 0 | yes |
| OS.20240203.42350.96.fits | 11:45:50 | NGC3245 | Kbb | 0.050 | 45 | 1 | yes |
| OS.20240203.42460.10.fits | 11:47:40 | NGC3245 | Kbb | 0.050 | 45 | 0 | yes |
| OS.20240203.42543.75.fits | 11:49:03 | NGC3245 | Kbb | 0.050 | 45 | 1 | yes |
| OS.20240203.42705.45.fits | 11:51:45 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.43383.79.fits | 12:03:03 | NGC3245 | Kbb | 0.050 | 600 | 1 | yes |
| OS.20240203.44064.64.fits | 12:14:24 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.44757.73.fits | 12:25:57 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.45436.78.fits | 12:37:16 | NGC3245 | Kbb | 0.050 | 600 | 1 | yes |
| OS.20240203.46116.63.fits | 12:48:36 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.46816.35.fits | 13:00:16 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.47493.63.fits | 13:11:33 | NGC3245 | Kbb | 0.050 | 600 | 1 | yes |
| OS.20240203.48172.98.fits | 13:22:52 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.48868.72.fits | 13:34:28 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.49547.50.fits | 13:45:47 | NGC3245 | Kbb | 0.050 | 600 | 1 | yes |
| OS.20240203.50226.88.fits | 13:57:06 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.50920.99.fits | 14:08:40 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.51598.82.fits | 14:19:58 | NGC3245 | Kbb | 0.050 | 600 | 1 | yes |
| OS.20240203.52278.78.fits | 14:31:18 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.52965.80.fits | 14:42:45 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.53647.19.fits | 14:54:07 | NGC3245 | Kbb | 0.050 | 600 | 1 | yes |
| OS.20240203.54326.00.fits | 15:05:26 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.55001.79.fits | 15:16:41 | NGC3245 | Kbb | 0.050 | 600 | 0 | yes |
| OS.20240203.56018.23.fits | 15:33:38 | NGC3245-PSF1 | Kbb | 0.050 | 30 | 0 |  |
| OS.20240203.56075.86.fits | 15:34:35 | NGC3245-PSF1 | Kbb | 0.050 | 30 | 1 |  |
| OS.20240203.56170.72.fits | 15:36:10 | NGC3245-PSF1 | Kbb | 0.050 | 120 | 0 |  |
| OS.20240203.56321.68.fits | 15:38:41 | NGC3245-PSF1 | Kbb | 0.050 | 120 | 1 |  |
| OS.20240203.56800.59.fits | 15:46:40 | HD101060(A0V) | Kbb | 0.050 | 3 | 1 | yes |
| OS.20240203.56832.21.fits | 15:47:12 | HD101060(A0V) | Kbb | 0.050 | 3 | 0 | yes |
| OS.20240203.56959.20.fits | 15:49:19 | HD101060(A0V) | Kbb | 0.050 | 10 | 1 | yes |
| OS.20240203.56998.98.fits | 15:49:58 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.57040.80.fits | 15:50:40 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.57146.88.fits | 15:52:26 | HD101060(A0V) | Kbb | 0.050 | 3 | 1 | yes |
| OS.20240203.57176.46.fits | 15:52:56 | HD101060(A0V) | Kbb | 0.050 | 3 | 0 | yes |
| OS.20240203.57690.54.fits | 16:01:30 | HD101060(A0V) | Kbb | 0.050 | 3 | 1 | yes |
| OS.20240203.57720.12.fits | 16:02:00 | HD101060(A0V) | Kbb | 0.050 | 3 | 0 | yes |
| OS.20240203.57823.65.fits | 16:03:43 | HD101060(A0V) | Kbb | 0.050 | 10 | 1 | yes |
| OS.20240203.58063.88.fits | 16:07:43 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58103.67.fits | 16:08:23 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58208.25.fits | 16:10:08 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58250.10.fits | 16:10:50 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58347.54.fits | 16:12:27 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58388.34.fits | 16:13:08 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58534.71.fits | 16:15:34 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58576.53.fits | 16:16:16 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58774.41.fits | 16:19:34 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58807.56.fits | 16:20:07 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58918.23.fits | 16:21:58 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.58960.50.fits | 16:22:40 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.59126.31.fits | 16:25:26 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.59166.90.fits | 16:26:06 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.59407.32.fits | 16:30:07 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.59447.61.fits | 16:30:47 | HD101060(A0V) | Kbb | 0.050 | 10 | 0 | yes |
| OS.20240203.59546.55.fits | 16:32:26 | HD101060(A0V) | Kbb | 0.050 | 20 | 0 | yes |
| OS.20240203.59597.56.fits | 16:33:17 | HD101060(A0V) | Kbb | 0.050 | 20 | 0 | yes |

## 4. Draft request to Keck for the rectification matrix

To: Jim Lyke (Keck OSIRIS instrument scientist), cc Tuan Do.
Subject: OSIRIS Kbb/050 rectification matrix for 2024-02-03 (PypeIt support)

Dear Jim,

We are adding Keck/OSIRIS support to PypeIt (the Python spectroscopic
reduction pipeline, https://pypeit.readthedocs.io).  Our approach keeps
the OSIRIS DRP's spatial rectification step (we are porting
`spatrectif_000` to numpy and will validate it against DRP reference
cubes) and uses the Keck rectification matrices as a required input.

For our first science-mode test we are using public KOA data from
2024-02-03 (NGC 3245 and the A0V standard HD101060, Kbb, 0.050" scale,
Keck I).  Could you send us, or point us to, the rectification matrix
appropriate for that date and mode (the file would be named like
`sYYMMDD_cNNN___infl_Kbb_050.fits`, epoch 2023-2024)?  If the matrices
for the other common modes (Kbb/035, Kn3/035, Kn3/050) are easy to
share, we would gladly include them as well.

We would also like to confirm two things with you: (1) that the
Dropbox/tkserver links in the OsirisDRP repository are no longer
maintained (the Kn5/035 test matrix `s150901_c008___infl_Kn5_035.fits`
is not downloadable), and (2) whether Keck would be willing to host the
current matrices at a stable URL so PypeIt could fetch them on demand.

We will of course credit the OSIRIS DRP (Lyke et al. 2017; Lockhart et
al. 2019) in the PypeIt documentation and in any publication.

Many thanks,
Xavier Prochaska

## 5. Follow-ups

- When the matrix arrives: put it in `RAW_DATA/keck_osiris/Kbb_050/calib/`
  and run the OH-line tilt test (design document, Section 8, test 3).
  *2026-10-12: the request was not needed.*  Keck publishes every
  matrix by validity epoch at http://tkserver.keck.hawaii.edu/osiris/ (HTTP only);
  the 2023-07-01 to 2025-03-10 epoch's Kbb/050 file,
  http://tkserver.keck.hawaii.edu/osiris/data/230701-250310/s240524_c006___infl_Kbb_050.fits.gz,
  was downloaded, gunzipped and placed in `Kbb_050/calib/`.  Its file
  date (2024-05-24) is later than the observations, which the M2
  latest-before selection rule does not handle (M4 task 11).
- Reduce the NGC 3245 A-B pairs with the DRP as well, if IDL is
  available, to obtain a second reference cube for the vet test.
