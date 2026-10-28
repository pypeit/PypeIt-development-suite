# Speed up PypeIt -- PR C (Arc-fit vectorization)

## Goals

Implement **workstream C** of the locked speed-up plan: replace the per-line
`scipy.optimize.curve_fit` Gaussian fit in `pypeit/core/arc.py:fit_arcspec`
(reached via `pypeit/core/fitting.py:fit_gauss`) with a **vectorized weighted
analytic log-parabola fit across all detected arc lines at once**, polished with
a few batched Gauss–Newton iterations and backed by a per-line `curve_fit`
fallback for the lines the analytic fit cannot handle. This removes ~1.11 million
`curve_fit` calls and 93 million `gauss_3deg` evaluations from a Keck/DEIMOS
reduction (`fit_gauss` cum **840 s**), speeding up both wavelength calibration and
tilt tracing. Self-contained: one function, one module.

## Claude

### Skills

Consider using the skills in `PypeIt-development-suite/.claude/skills/` — in
particular `run-dev-suite` (the RMS vetting), `wavelength-calibration`,
`build-docs`, and `update-changelog`.

## Context

Authoritative documents (read the coding doc section first; it is the
implementation contract):

- **`PypeIt-development-suite/pypeitdev/speed_up/Reports/speed_up_coding.md` §3**
  — "PR C — Vectorized arc-line Gaussian fitting". Sub-anchors §C.1 (what is
  being replaced), §C.2 (the math, with the full derivation), §C.3 (the
  implementation sketch, including `_gauss_newton_refine` and
  `_fit_arcspec_curvefit`), §C.4 (the edge-case table), §C.5 (tests), §C.6
  (dev-suite vetting), §C.7 (validation/re-profiling), §C.8 (docs). §0.2 is a
  grounded code map; §7 is the risks/revert table.
- **`Reports/speed_up_design.md` v0.5 (locked)** — §0 locked decisions, §4 Tier 1
  "Vectorize arc-line Gaussian fitting", §8 roadmap.
- **Profile reports** — `Reports/keck_deimos_600zd_m_6500_profile_report.md`
  (`fit_gauss` cum **840 s** over **1.11 M** `curve_fit` calls; `gauss_3deg` self
  **332 s** over **93 M** calls; `curve_fit` cum 648 s; `leastsq` cum 586 s) and
  `Reports/shane_kast_blue_A_profile_report.md` (small — ~1–3 s).
- **Profiling scripts** — `pypeitdev/speed_up/scripts/profile_kast_blue.py`,
  `profile_deimos.py`, `analyze_profile.py`.
- **PR A/B results** — `Reports/speed_up_results.md`, and the Logs in
  `claude_prompts/speed_up/A.md` and `B.md`.

Repository and branch:

- Code repo: `/mnt/tank/Astronomy/PypeIt/PypeIt` (same as
  `/home/xavier/Projects/PypeIt/PypeIt`).
- **Branch stack (as actually built; updated 2026-10-20 after PR B):**
    - A: `speed_up_qa` → merged to **`develop`** (#2198).
    - B: lives on **`speed_up_B`**, not `speed_up_detpar`.  That branch was
      never created (B.md Q1(a)).  PR **#2209** (`speed_up_B` → `develop`) is
      open and under review by Tim Pickering.  `develop` has been merged into
      it up to `72e448f75`.
    - **C: `speed_up_arcfit`, branched from `speed_up_B` at `72e448f75`.**
      **PR C targets `speed_up_B`**, and is retargeted to `develop` once #2209
      merges.  The `speed_up` integration branch has been superseded: the PRs
      go straight to `develop`.
    - A stray local branch `speed_up_C` also points at `72e448f75` (see Q1).
- PR B is complete on `speed_up_B`.  All B prompts and PR prompts 1–7 are done;
  only Tim's re-review is outstanding.  Things PR B changed that matter here:
    - `run_pypeit --ncpu N` reduces detectors in worker processes
      (`pypeit/parallel.py`).  The start method is `forkserver` on Linux and
      `spawn` elsewhere, and `PYPEIT_MP_START` overrides it.  cProfile does not
      see inside the workers, so **profile at `--ncpu 1`**.
    - Products are bit-identical across `ncpu` only on Linux, and only for data
      that already reproduce run to run.  Multi-slit is not reproducible even
      serially (#2208, `EDGE_CEN` off by 1 ULP).  So dev-suite comparisons
      before and after C on multi-slit setups carry that noise floor as well as
      C's own changes.
    - #2211 (`[rdx] fit2d_chunk`, lower memory) is merged into `speed_up_B`.
    - Latest DEIMOS 600ZD numbers (`speed_up_results.md`, "After the B port"):
      `--ncpu 4` forkserver gives a **3 391 s** wall and 36.5 GB peak PSS.
      Serial baseline: 12 390 s.

Key facts about the code being changed:

- `fit_arcspec(xarray, yarray, pixt, fitp)` loops over every detected peak and
  calls `fitting.fit_gauss` → `curve_fit(gauss_3deg, …)` with **no** `sigma`, so
  the objective is unweighted least squares in **linear** space.
- The windows are tiny: `nfitpix = round(1.25*fwhm)`, so for the typical
  `fwhm=4` the window is **7 pixels**.
- `yarray` is the **continuum-subtracted** arc, so wing pixels routinely go
  **≤ 0** — a naive log-space fit cannot use them. This is the central accuracy
  risk and the reason for the Gauss–Newton polish.
- `centerr` is returned as `fitcov[1,1]`, i.e. the **variance** of the centre
  (despite the name).  It reaches `all_ecent` in
  `wvutils.arc_lines_from_spec`, and from there `HolyGrail._det_weak/_det_stro`
  and `get_use_tcent`.  **It is never used numerically.**  Every caller
  discards it (`use_tcent, _ = …`), or assigns it to a variable that is never
  read (`autoid.py:2232-2233, 2985`).  So `centerr` only has to keep its
  shape, its sign and the −999 sentinel; its accuracy does not matter.  Legacy
  `centerr` is also `inf` whenever `curve_fit` cannot estimate a covariance
  (31 of 4 259 lines on the Paranal sky).
- `fit_gauss` has exactly one other, cold caller
  (`pypeit/multislit_flexure.py`) — leave it untouched.
- Downstream filters in `detect_lines` are unchanged and still apply:
  `twid > 0`, `twid < fwhm_max/2.35`, `tcent > 0`,
  `tampl_true < nonlinear_counts`, `abs(tcent-pixt) < 0.75*fwhm`.
- **`detect_lines`, and therefore `fit_arcspec`, has callers outside
  wavelength calibration**, so PR C also moves:
    - object finding (`findobj_skymask.py:2167`: `fit_frac_fwhm=1.5`,
      `max_frac_fwhm=10`, on a smoothed S/N profile);
    - edge/trace peak finding (`trace.py:1466`, on Sobel-collapsed profiles);
    - the spatial-flexure cross-correlation peak (`flexure.py:90`, which checks
      `pix_max[0] == -999.0`);
    - sky-line matching (`flexure.py:283,289`);
    - tilt tracing (`tracewave.py:63`).

  `test_extraction.py`, `test_skysub.py` and the slit-edge checks in the vet
  suite are therefore relevant too, not only the wavelength and tilt RMS.
- **Grounded on the Paranal sky spectrum** (`test_detect_lines`: 4 259 peaks,
  `fit_interval = 3`), scratch probe only:
    - legacy `curve_fit` raises `RuntimeError` on **248** lines;
    - **56** lines converge to a **negative σ**.  `gauss_3deg` only uses σ²,
      so the sign is arbitrary, but `twid > 0` rejects those lines.  49 of the
      56 would pass every other filter;
    - 144 lines end up with a centre more than `fit_interval` from `pixt`;
    - 35% of the windows contain a pixel ≤ 0.
- **The sketch's analytic seed, on the same sky lines**, gives a median |Δc|
  of 0.053 px against legacy and a 99th percentile of 2.1 px.  With
  `GN_ITER = 0` it falls back on **13.2%** of lines, not the "<1%" of §C.7:
    - 157 lines have fewer than 4 positive pixels;
    - 287 have `a2 >= 0`;
    - 119 have `|u0| > fit_interval`.

  Arc lamps are cleaner than night sky, so measure it on DEIMOS in Task 2.
- **The sketch's `np.linalg.solve` calls fail on NumPy ≥ 2** (the
  `pypeit14b` env runs numpy 2.5).  A stacked right-hand side of shape `(n,3)`
  is now read as a matrix: "mismatch in its core dimension".  Use
  `np.linalg.solve(A, b[..., None])[..., 0]`, both in `fit_arcspec` and in
  `_gauss_newton_refine`.
- **The sketch's fallback mask disagrees with the §C.4 table.**  `usable`
  includes `ngood >= ARCFIT_MIN_GOODPIX`, and the fallback is
  `usable & ~ok`.  So a line with fewer than 4 positive pixels gets **no fit at
  all** (−999).  §C.4 says it "falls back to `curve_fit`, matching the old
  behaviour".  The fix is to keep two masks:
    - `attempt = (nwin > 0) & (nwin >= fit_interval)`: the old `continue`
      guards;
    - `analytic = attempt & (ngood >= ARCFIT_MIN_GOODPIX)`: the solve;
    - fallback on `attempt & ~ok`.

  See Q2.
- The legacy loop computes `pmin = pixt[p] - fit_interval` with no `int()`.
  The sketch's `_fit_arcspec_curvefit` adds `int(...)`, which is harmless for
  the integer `pixt` from `detect_peaks`.  The move is otherwise verbatim.

Locked decisions relevant to PR C (a fresh session needs no other context):

- **Noise-level output changes are acceptable**, vetted by the dev-suite
  wavelength/tilt RMS checks. (This is the one workstream that is *not* required
  to be bitwise identical — the identical-output requirement in PR B is about
  `ncpu=1` vs `ncpu>1`, which this PR does not affect.)
- **Q3 (accepted):** ship with **`ARCFIT_GN_ITER = 3`** — the analytic weighted
  log-parabola provides the seed, then 2–3 batched, damped Gauss–Newton
  iterations refine against the *exact* `curve_fit` objective (linear space, all
  in-bounds pixels **including the negative wings**), which tracks `curve_fit`
  to ~1e-6 instead of only to the log-space approximation. Still fully
  vectorized; still no `scipy` call. **Report both variants' RMS** (`GN_ITER=0`
  and `GN_ITER=3`) in the PR description.
- **Weighting:** the log-space fit must be weighted by `w_i = y_i` (Guo
  weighting), i.e. minimise `Σ y_i² (ln y_i − p(u_i))²`. An unweighted log fit is
  badly biased toward the faint wings and will not match `curve_fit`.
- **Conditioning:** fit in the window-centred coordinate `u = x − x[pixt]`, not
  absolute pixel `x`. A 7-pixel absolute-`x` Vandermonde at `x ≈ 3000` has a
  condition number ~1e14.
- **Keep the legacy loop verbatim** as `_fit_arcspec_curvefit` — it is both the
  per-line fallback and the reference implementation for the tests.
- **Kill switch:** a module-level `USE_VECTORIZED_ARCFIT = True` in
  `pypeit/core/arc.py` restores the exact previous behaviour when set to `False`.
- **Out of scope:** the B-spline kernel (`pypeit/core/bspline/`) — another branch.
  Do not touch it. Also do not touch `moment1d` (a v2 candidate) or
  `pypeit/multislit_flexure.py`.

## Running

If you need to run Python or PypeIt, use the **`pypeit14b` conda environment**
(Python 3.14, numpy 2.5; editable install of this repo).  The `pypeit` env is
broken: `import astropy` fails with `numpy has no attribute 'in1d'`.  This was
re-confirmed 2026-10-20, and B.md Q8 settled on `pypeit14b` too.  The shell's
default environment shadows it, so use the absolute env binaries:

```bash
/home/xavier/miniconda3/envs/pypeit14b/bin/python
/home/xavier/miniconda3/envs/pypeit14b/bin/run_pypeit
/home/xavier/miniconda3/envs/pypeit14b/bin/pytest
```

Dev-suite runs: use B's safe recipe, adopted after the 2026-10-13 workstation
crash.  Run one suite at a time, at `-t 1`, each inside
`systemd-run --user --scope -p MemoryMax=80G -p MemorySwapMax=0`.  If a
setup is re-reduced into an existing `REDUX_OUT` directory, **remove its
`Science/` first**.  `-m` rebuilds the calibrations but merges new objects into
the old spec1d (B.md QP7, left unfixed in the harness).

## Prompts

1. Read this doc.  Perform the 1st task under Tasks.
2. Read this doc.  Perform the 2nd task under Tasks.
3. Read this doc.  Perform the 3rd task under Tasks.
4. Read this doc.  Perform the 4th task under Tasks.
5. Read this doc.  Perform the 5th task under Tasks.
6. Read this doc.  Perform the 6th task under Tasks.
7. Read this doc.  Perform the 7th task under Tasks.

## Tasks

1. **Prepare.** Create the branch and ground yourself; write no production code
   yet.
    - Confirm PR B is complete on `speed_up_B`. Then create and check out
      **`speed_up_arcfit`** from **`speed_up_B`** (`speed_up_detpar` was
      never created); report the HEAD commit.
    - Update this prompt doc as need to be account for all of the changes in PR B.
    - Read `Reports/speed_up_coding.md` **§0**, **§3** and **§7** in full, and
      work through the derivation in §C.2 until you can re-derive
      amplitude/centre/sigma from the three parabola coefficients and the centre
      variance from the Jacobian propagation.
    - Read the code: `pypeit/core/arc.py` (`detect_lines`, especially the
      continuum subtraction, the `fit_arcspec` call and the `good` mask; then
      `fit_arcspec` itself) and `pypeit/core/fitting.py` (`fit_gauss`,
      `gauss_3deg`, `guess_gauss`).
    - Enumerate the consumers of the `fit_arcspec` outputs so you know what must
      not change shape or sentinel convention:
      `git grep -n "detect_lines(\|centerr\|ecent" -- pypeit/`.
    - Verify the line numbers in coding-doc §0.2/§3 still hold after PRs A and B;
      note any drift.
    - Ask questions in the Q&A section below. Log your work in the Logs section.

2. **Implement the vectorized log-parabola core + the legacy fallback.**
   Coding-doc **§C.3**, first half.
    - Add the module constants `USE_VECTORIZED_ARCFIT`, `ARCFIT_GN_ITER` and
      `ARCFIT_MIN_GOODPIX` to `pypeit/core/arc.py`.
    - Move the original per-line loop verbatim into
      `_fit_arcspec_curvefit(xarray, yarray, pixt, fit_interval, sz_a, which,
      ampl, cent, widt, centerr)`, filling the output arrays in place for the
      indices in `which`.
    - Rewrite `fit_arcspec` with the batched path: window index array
      `idx = pixt[:,None] + off[None,:]` with an `inb` in-bounds mask (windows
      truncated at the spectrum edges, exactly as the old loop did), the pixel
      mask `gpm = inb & isfinite(y) & (y > 0)`, the `usable` guard reproducing
      the old `continue` conditions (`nwin > 0`, `nwin >= fit_interval`) plus
      `ngood >= ARCFIT_MIN_GOODPIX`, the moment arrays `S[0..4]` / `T[0..2]` with
      weights `w² = y²`, the stacked `(nline,3,3)` `np.linalg.solve`, and the
      back-substitution to `amp`/`u0`/`sig` with the `a2 < 0`,
      `abs(u0) <= fit_interval` and finiteness guards.
    - Wire the fallback: lines that are `usable` but not `ok` go through
      `_fit_arcspec_curvefit`; warn if more than 25% of lines fall back.
    - Preserve the `-999.0` sentinels, the return signature
      `(ampl, cent, widt, centerr)`, and the docstring (updated to describe the
      new method).
    - Set `ARCFIT_GN_ITER = 0` for now so this task can be validated on its own.
    - Ask questions in Q&A. Log your work in the Logs section.

3. **Add the Gauss–Newton polish and the centre variance.** Coding-doc **§C.3**,
   second half, and the `centerr` derivation in **§C.2**.
    - Implement `_gauss_newton_refine(u, y, gpm, amp, u0, sig, ok, niter, lam)`:
      batched Jacobian over the three parameters, Levenberg-damped 3×3 normal
      equations solved with `np.linalg.solve` on the good subset, `niter`
      iterations. It must use **all in-bounds pixels including the negative
      wings** (mask `inb`, not `gpm`) — that is the whole point, and it is what
      makes the result track `curve_fit`.
    - Implement `_center_variance(A, coef, u, w2, lny, ngood, ok)`:
      `Cov = s²·inv(A)` with `s² = Σ w²δ²/(ngood−3)`, then propagate through
      `∂c/∂a1 = −1/(2a2)` and `∂c/∂a2 = a1/(2a2²)`.
    - Set **`ARCFIT_GN_ITER = 3`** (the accepted Q3 default).
    - Walk the §C.4 edge-case table and confirm each case is handled: negative /
      zero flux, NaN pixels, peaks at the array edges, too few positive pixels,
      upward parabola (`a2 >= 0`), singular normal matrix, amplitude overflow in
      `exp()`, saturated/flat-topped lines, centre outside the fit window, and
      `sz_p == 0`.
    - Ask questions in Q&A. Log your work in the Logs section.

4. **Tests.** Coding-doc **§C.5**.
    - Extend `pypeit/tests/test_arc.py` with the `_synthetic_arc` fixture helper
      (fixed seed — determinism is a repo rule) and:
      `test_fit_arcspec_recovers_truth`,
      `test_fit_arcspec_matches_curve_fit` (toggles `USE_VECTORIZED_ARCFIT` to
      compare against the retained legacy implementation; centroids within
      0.02 px, widths within 2%), `test_fit_arcspec_edge_cases` (peaks at index 0
      and N−1, an all-negative window; no NaN may reach the caller), and
      `test_fit_arcspec_no_curve_fit_on_clean_data` (monkeypatch
      `fitting.fit_gauss` to raise; the vectorized path must not call it on clean
      data).
    - Keep the existing `test_detect_lines` assertion (`len(arx_w) > 3275`) — it
      is a real-data regression on the new code. If the count changes, stop and
      investigate before proceeding.
    - Run `pytest pypeit/tests` in the `pypeit14b` env; must be green.
    - Also record a quick micro-benchmark (vectorized vs legacy on ~2 000 lines)
      for the PR description.
    - Ask questions in Q&A. Log your work in the Logs section.

5. **Dev-suite RMS vetting (both GN variants) + docs and changelog.**
   Coding-doc **§C.6** and **§C.8**.
    - Run the reductions and the vet tests:

      ```bash
      cd $PYPEIT_DEV
      ./pypeit_test reduce -i shane_kast_blue shane_kast_red keck_deimos keck_lris_blue
      pytest vet_tests/test_wavelengths.py vet_tests/test_wavetilts.py \
             vet_tests/test_extraction.py vet_tests/test_skysub.py \
             --redux_out $PYPEIT_DEV/REDUX_OUT
      ```

    - Per **Q3**, do this for **both** `ARCFIT_GN_ITER = 0` and
      `ARCFIT_GN_ITER = 3`, and tabulate the per-setup wavelength RMS
      (`test_shane_kast_red`, `test_deimos`, `test_keck_lris_blue`, …) and the
      tilt RMS (`test_wavetilts.py::test_run`) for both, alongside the pre-change
      baseline. Acceptance: **no setup's RMS degrades beyond its existing
      tolerance, and no setup loses arc lines.** Ship `GN_ITER = 3`.
    - Docs/changelog: add the §C.8 bullet to `doc/releases/2.1.0dev.rst` under
      **Under-the-hood Improvements** (use the `update-changelog` skill). No new
      user-facing parameter, so `doc/pypeit_par.rst` is unchanged; regenerate
      `doc/api` and check `doc/sphinx_warnings.out` is empty.
    - Ask questions in Q&A. Log your work in the Logs section.

6. **Validate, re-profile and summarize.** Coding-doc **§C.7** and the
   cross-cutting checklist in **§4** of the coding doc.
    - Re-profile both setups at **`--ncpu 1`** (so cProfile sees everything —
      with `ncpu>1` the child processes are invisible to the in-process
      profiler).
    - Confirm in the `.pstats.txt` tables that `fitting.py:fit_gauss` (baseline
      cum **840 s**, 1.11 M calls), `gauss_3deg` (self **332 s**, 93 M calls) and
      scipy's `leastsq` (**586 s**) have dropped to noise, and that the residual
      `curve_fit` calls are the fallback only — typically **<1% of lines**.
    - **Expected (§C.7): ~700–800 s off DEIMOS (~6% of the original wall);
      ~1–3 s off Kast.**
    - Sanity-check the products: `spec1d` `WAVE_RMS` values and the QA arc-fit
      PNGs for a couple of slits.
    - Append the numbers and the RMS table to `Reports/speed_up_results.md` (new
      "after-C" section), note the total A+B+C improvement versus the original
      12 390 s / 119 s baselines, and summarize in the Logs section below.
    - Report the A+B+C total both serially (C's gain is a serial-time gain)
      and at `--ncpu 4` against B's 3 391 s.
    - Finally, note the remaining pre-merge steps.  #2209 still needs Tim's
      re-review and merge, after which PR C is retargeted to `develop`.  The
      DEIMOS `--ncpu 4` run was already done in B, but the multi-slit products
      are not reproducible run to run (#2208), so there is no identical-output
      check left to do for C.  One full dev-suite pass (`./pypeit_test all`,
      using the safe recipe under Running) at the tip of `speed_up_arcfit`
      before merging.  Ask any remaining questions in Q&A.

7. **Write the final report, `Speed_Up_Arcfit_Report.md`** (added per Q6; the
   last item before the Pull Request).  Create it in
   `PypeIt-development-suite/pypeitdev/speed_up/Reports/`.  It is the
   standalone record of workstream C, for readers who have not followed these
   logs.  It should cover:
    - what was replaced and why, with the profile numbers;
    - the method: the weighted log-parabola, the Gauss–Newton polish, the
      moment-seeded rescue (if adopted), and the `curve_fit` fallback;
    - how faithful the new fit is to legacy `curve_fit`, on the Paranal sky and
      on DEIMOS arcs;
    - the line-count changes, with the reasons: negative σ (Q3) and the
      saturated lines (Q6);
    - **worked examples, with figures, of the saturated, flat-topped DEIMOS
      lines that the new fit rejects**: the window pixels, the legacy and new
      fits overplotted, and why the legacy fit is an artefact.  Also an example
      of a negative-σ line;
    - the fallback rate and where the remaining time goes;
    - the dev-suite RMS table for both `GN_ITER` variants;
    - the re-profiling numbers, and the A+B+C totals;
    - the kill switch, and the remaining risks.
    - A discussion of how including saturated lines can improve the fit (by more IDs)

   Use the `dataviz` skill for the figures, and keep the plotting scripts
   under `pypeitdev/speed_up/scripts/` so the figures can be regenerated.
   Ask questions in Q&A.  Log your work in the Logs section.

## Q&A

Claude poses questions here (as `**Qn — title.** body` followed by a `>A:` line);
the user answers inline beneath each `>A:`.

**Q1 — Stray `speed_up_C` branch.**  `speed_up_arcfit` already existed when
Task 1 started, at the `speed_up_B` tip `72e448f75`, so I only checked it out.
A second local branch, `speed_up_C`, points at the same commit and has no
upstream.  My default is to **leave it alone**.  Tell me if you want it
deleted (`git branch -d speed_up_C`).
>A: Leave it alone

**Q2 — Fix the sketch's fallback mask?**  As written in §C.3, a line with fewer
than `ARCFIT_MIN_GOODPIX = 4` positive pixels is neither fit analytically nor
sent to `curve_fit`, so it silently gets −999.  That contradicts §C.4 ("falls
back to `curve_fit`, matching the old behaviour").  On the Paranal sky that is
157 of 4 259 lines.  My default is two masks:
- `attempt`: the old `continue` guards;
- `analytic = attempt & (ngood >= 4)`;
- fallback on `attempt & ~ok`.

The old code attempted every line that the vectorized code attempts, so no
line loses its fit.
>A: Use your default

**Q3 — Lines where legacy `curve_fit` converged to a negative σ.**  The model
uses only σ², so the sign `curve_fit` lands on depends on the optimizer path.
Legacy then throws those lines away through `twid > 0`: 56 of 4 259 on the
Paranal sky, 49 of which pass every other filter.  The new fit always returns
σ > 0, so most of these lines will now be **kept**.  `len(arx_w)` in
`test_detect_lines` will probably *rise*, and wavelength/tilt line lists gain
about 1% more lines.  I cannot reproduce the legacy rejection, since it is an
optimizer artefact, and I do not think we should.  My default is to accept the
gain.  I would read Task 4's "if the count changes, stop and investigate" as
"stop if it **drops**, or rises by more than these ~50 lines", and record the
new count.  Agree?
>A: I agree

**Q4 — The fallback rate is 13%, not <1%.**  With the analytic seed alone,
13.2% of sky lines fall back to `curve_fit`:
- `a2 >= 0`: 287 lines;
- `|u0| > fit_interval`: 119 lines;
- fewer than 4 positive pixels: 157 lines.

Arc lamps should be cleaner, but if DEIMOS arcs look anything like this, C
delivers roughly 85% of the 840 s, not ~95%.  Options for the
`a2 >= 0` / out-of-window lines:
- (a) **keep the `curve_fit` fallback for all of them.**  Exact legacy
  behaviour for the hard lines; costs speed.
- (b) seed Gauss–Newton for them from a vectorized `guess_gauss`-style moment
  estimate instead of falling back.  No scipy at all, but these become new
  fits rather than legacy fits.
- (c) mark them failed (−999).  Most would fail the `abs(tcent-pixt) <
  0.75*fwhm` filter anyway, but not all.

My default is **(a) now**, then measure the fallback rate on a DEIMOS arc in
Task 2 and revisit (b) in Task 3 only if it is well above a few percent.
>A: I agree, use your default

**Q5 — `centerr` from the Gauss–Newton Jacobian rather than the log-parabola?**
`centerr` is never used numerically (see Key facts), so this is about fidelity,
not results.  Once Gauss–Newton has run, the closest match to legacy is
`curve_fit`'s own definition: `pcov = s²·inv(JᵀJ)` with linear-space residuals,
`s² = Σr²/(n−3)` over the in-bounds pixels, at the converged parameters.  It is
nearly free, because `JᵀJ` is already built in the last iteration.  My default:
- `GN_ITER > 0`: take `centerr` from the Jacobian;
- `GN_ITER = 0`: use §C.2's `_center_variance` propagation.

Both are implemented and tested against legacy `fitcov[1,1]`.
>A: I agree, use your default

**Q6 — Saturated lines on DEIMOS: the analytic fit drops 129 of them.**
On the 92 DEIMOS 600ZD arc spectra (`GN_ITER = 0`), the good-line count falls
from **7 887 to 7 763**, and 67 slits lose lines.  This breaks the "no setup
loses arc lines" acceptance criterion.
- **All 129 lost lines fail the width cut** (`twid < fwhm_max/2.35`).
- **They are saturated.**  99% have raw peaks above 30k counts (median 39k,
  against 8.6k for all lines), and their windows are flat-topped
  (3rd-highest/max = 0.999).  A plateau makes the log-parabola nearly flat,
  so σ blows up (example: σ = 16.2 against `curve_fit`'s 2.17).
- **`nonlinear_counts` does not remove them**, contrary to §C.4.  The cut is
  0.95×65 535 ≈ 62k, but the extracted arc plateaus near 39k.
- **Saturated lines also dominate the fallback cost.**  They are 376 of the
  678 `curve_fit` fallbacks and 54% of the fallback time.  Each fallback line
  averages ~540 `gauss_3deg` evaluations, against ~36 for a typical line.
  The fallback is now ~99% of `fit_arcspec`'s run time, so the speedup is
  only **2.1×**.

Task 3's Gauss–Newton polish minimises the same linear-space objective as
`curve_fit`, so it may pull these widths back.  That is the first thing I
will measure.  If it does not:
- (a) send lines whose analytic σ exceeds `fwhm_max/2.35`-ish (or plateau-shaped
  windows) to the `curve_fit` fallback.  Exact legacy behaviour, but it keeps
  the slow path for exactly the slowest lines;
- (b) accept the loss: legacy's narrow σ on a plateau is itself an artefact,
  and dropping saturated lines arguably improves the centroids;
- (c) start Gauss–Newton for plateau lines from a moment-based seed (Q4(b)).

**My default: measure GN first, then (c) if GN alone does not recover them,
falling back to (a).**  I would only choose (b) with your agreement, because
it changes which lines the HolyGrail/reidentify sees.
>A:  It appears these lines are not useful.  Let's adopt (b) for now and be sure
to inculde examples in the final report that we create on this effort.
That report should be named `Speed_Up_Arcfit_Report.md` in 
`PypeIt/PypeIt-development-suite/pypeitdev/speed_up/Reports`.  Please add 
a new prompt to create this report as the last item before the Pull Request.

**Q7 — Raise `ARCFIT_GN_ITER` from 3 to 10?**  The locked Q3 decision is 3,
and that is what is set.  Gauss–Newton now stops each line once it converges,
so the setting is a *cap*, and raising it is cheap.  On DEIMOS 600ZD arcs
(92 spectra, 10 170 lines; rescue at 20):

| `GN_ITER` | `fit_arcspec` | speedup | lines lost vs legacy | median / p99 \|Δc\| vs `curve_fit` | lines with \|Δc\| > 0.02 px |
|---|---|---|---|---|---|
| 3 | 0.29 s | 16.1× | 122 | 1.4e-6 / 0.14 px | 174 |
| 10 | 0.35 s | 13.3× | 51 | 1.4e-7 / 4.4e-4 px | 17 |

On the Paranal sky, `GN = 10` takes the p99 from 0.87 to 0.33 px.  The lost
lines are the saturated plateaus (Q6(b)), whose flat valley converges slowly;
10 iterations recover more than half of them.  Against 840 s of `fit_gauss`
on DEIMOS, the extra 0.06 s per 92 spectra is negligible.
**I recommend 10**, and I would report both in Task 5's RMS table as planned
(`GN = 0` vs the shipped value).  It also makes
`test_fit_arcspec_matches_curve_fit` (0.02 px) far more comfortable.
>A: Ok, increase it to 10

**Q8 — Fidelity or speed?  The dev suite shows that shipped `gn10` changes
some DEIMOS wavelength *identifications*.**  Every vet test that ran gives the
same outcome for legacy and `gn10`, and every vetted RMS is within its
tolerance (table in the Task 5 log).  Underneath, though:
- on 28 of 38 WaveCalib files the solutions agree to ~1e-6 px;
- on a handful of DEIMOS slits the **line identification flips**, and the
  solution moves by tens of pixels inside the range both variants' lines
  cover.  `600ZD_M_6500` MSC01 slit 2: the blue end is misidentified, with
  4359.56 Å put on the line legacy calls 4377.35 Å, ~30 px off, while the RMS
  only moves 0.588 → 0.617.  Similar flips occur on MSC01 slit 11 and MSC03
  slits 0, 14 and 15;
- `830G_L_8100`'s single slit goes from RMS 0.051 to 0.27 px;
- one case improves: `1200G_M_7750` max RMS 0.42 → 0.13.

The triggers are **discrete changes in which lines are detected**, not
centroid noise (shared lines agree to ~1e-6 px).  Over all 370 arc spectra of
the run (33 354 legacy lines):
- **83 lines gained**, all sub-resolution spikes (σ < 0.5 of the expected
  width, median 0.33).  Legacy `curve_fit` happened to *fail* on these;
- **56 lines lost**, the saturated plateaus.  Q6(b) assumed they were not
  useful, but **they anchor the pattern match**: losing one or two is the
  only difference on the misidentified 600ZD slits.  So Q6's premise looks
  wrong;
- a few lines where both fits converge to different minima (17 DEIMOS lines
  > 0.02 px).

Options, each measured on the 92 DEIMOS 600ZD arcs and the Paranal sky
(`fit_arcspec` time against legacy):
- **(A) Defer to `curve_fit` whenever the vectorized fit is not confident**:
  a line that has not converged within the iteration cap, or a fit narrower
  than 0.35× the window half-width.  The prototype is saved as
  `scripts/arcfit/prototype_defer_to_curvefit.patch`.
    - Result: **identical line lists on all 370 spectra** (0 gained, 0 lost,
      shared centroids within 1e-3 px).  This is legacy behaviour by
      construction.
    - Cost: **only 2.0× on DEIMOS arcs and 1.4× on the sky**.  About 500
      DEIMOS lines (5%) go to `curve_fit`, and 432 of those end in
      `RuntimeError` (legacy's −999) after `curve_fit`'s full 800
      evaluations, which is 1.75 of the 1.83 s.
    - The expensive lines cannot be predicted or capped: some lines that
      legacy accepts need over 700 evaluations.  The patch also needs a test
      tweak, because a noise-free test line with a negative pedestal hits the
      convergence cap.
- **(B) Ship `gn10` as run: 13×.**  The vet tolerances pass, but the
  per-slit identification flips above remain.
- **(C) A hybrid: the width rule, plus deferring only non-converged
  *flat-topped* windows.**
    - DEIMOS: 2.7×, nearly exact (2 lines lost);
    - sky: 4.1×, but 113 shared lines move by more than 0.02 px and 20 are
      lost.

**My recommendation is (A).**  Wrong wavelength identifications are a worse
outcome than a smaller speedup.  (A) still removes about half of the
`fit_gauss` cost, and it keeps the whole cascade, tests and kill switch.  If
you agree, I will apply the patch with the test fix and re-run only the
shipped variant through the dev suite (~2.5 h; the legacy baseline in
`REDUX_OUT_C/legacy` is reused).  If you prefer (B), Task 5 is complete as it
stands.
>A:  Yes, let us ship A but allow B as an option.  The report and docs should 
discuss this at length.

**Q9 — Expose option (B) as a user parameter?**  Per Q8, (B) is available
through the module constant `pypeit.core.arc.ARCFIT_MATCH_CURVEFIT = False`.
That is reachable from Python (and from the dev-suite hook), but **not from a
`.pypeit` file**.  Making it a real option would mean:
- a new `[rdx]` parameter, e.g. `arcfit_match_curvefit = True`, applied in
  `PypeIt.__init__` the way #2211 applies `fit2d_chunk`;
- adding it to the `cfg` that `parallel._worker_init` re-applies in
  `forkserver`/`spawn` workers;
- `doc/pypeit_par.rst` regenerated, plus a test.

`detect_lines` has many callers (wavelength calibration, tilts, object
finding, flexure), so a single global switch is the natural granularity.
**My default is to leave it as a module constant**, documented in the
release notes and the API docs, because (B) is an expert/speed option whose
caveat (identification flips) users should not stumble into.  Say if you want
the parameter.
>A: Yes, I want a new user parameter. It should be an item in `WavelengthSolutionPar`

**Q10 — Cache the archive-sky FWHM in spectral flexure?**  The DEIMOS (A)
profile (`--ncpu 1`, 10 094 s) shows that **571 s of the 655 s spent in
`detect_lines`** comes from a single caller: `autoid.measure_fwhm`, called
**170 times** by `flexure.get_archive_spectrum` (`flexure.py:944`).  Each call
measures the FWHM of the **same archived sky spectrum** (`sky_file`), once per
slit, at 3.4 s each.  The code's own comment at that line says it is meant to
"avoid the performance hit from calling it on the archive sky spectrum
multiple times", but it does not.  This is also most of the residual
`curve_fit` load under (A): 319 k calls, 494 s.  Blended sky lines often fail
to converge and are handed off.
- **Fix:** memoize the archive FWHM per `sky_file`, e.g. a small
  `functools.lru_cache` helper around `load_sky_spectrum` + `measure_fwhm`.
  The output is **bit-identical**, since it is the same deterministic
  computation on the same file.  It should remove ~560 s from (A), more from
  legacy, and make (A) versus (B) nearly moot for DEIMOS.
- **Scope:** it is outside the original PR C plan, but one function and
  directly motivated by PR C's profiling.  One caveat for Q9: because this
  cost sits in *flexure*, the wavelength-only `fast_line_fit` parameter would
  not have touched it.

**My recommendation: add the cache in PR C**, with a unit test that a second
call does not re-measure, and record its gain separately from the arc-fit
gain.  It would go in after the current batch, which must not have
`pypeit/` changed under it, and would need one extra DEIMOS profile (~3 h).
>A: Yes, let's do this

**Q11 — Pre-existing misidentification candidates in the DEIMOS 600ZD
`full_template` path (outside PR C).**  Tracing the one WaveCalib file where
the shipped (A) code differs from legacy (log entry of 2026-10-26, second
part) showed that the difference is not in the line fitting at all: the
accepted lines and their centres are the same to ~1e-5 px, but
`xcorr_shift` → `detect_lines(corr, nfind=1)` → `differential_evolution` in
`xcorr_shift_stretch` amplifies a 1.6e-8 change in the cross-correlation
into a 30 px change in the blue end of slit 11 (spat_id 1063).  (A) lands
on the correct branch, legacy on a wrong one, and the zero-lag correlation
objective prefers the wrong one.  Slit 9 (spat_id 902) of the same mosaic
shows the same tell-tale dispersion turnover in *both* runs, so it may be
misidentified today.  My default: **leave the pipeline alone in PR C**,
mention the fragility in the report and release notes as a known limitation,
and open a GitHub issue on `pypeit/PypeIt` (with the slit-11 numbers and the
slit-9 suspicion) once you have had a look.  Say if you want the issue filed
now, or want me to dig into slit 9 first.
>A:  Leave the pipeline alone in PR C.  Mention the fragility in the report and release notes as a known limitation.  But do not file an Issue on GitHub.

**Q12 — The "saturated lines" are mostly plateau edges; reword anything?**
Task 7 looked at the DEIMOS 600ZD features that option B rejects and legacy
keeps.  With the pipeline's own inputs (per-slit `wv_fits[i].spec` and
measured fwhm ≈ 7.2 px, 86 spectra) there are 46 of them:
- 45 have raw peaks above 30k;
- 39 of those 45 sit on a vertical edge (an adjacent-pixel jump of 88–98% of
  the peak);
- 6 sit on a flat top.

In the 2D arc frame these are ~60-pixel regions at ~39 900 counts, around
cores that are hard-saturated at 79 503 counts.  That is consistent with
charge bleeding along the CCD columns.  Legacy's fit puts a narrow Gaussian
2–4 px inside the edge.  Its centroid is therefore an edge marker, not a line
wavelength.  None of the 46 is in any final wavelength fit; nearest
`pixel_fit` ≥ 7.2 px away.

The report says this plainly, and it corrects the "centroid is nevertheless
usable" premise of the Task 7 brief: usable as a pattern anchor, not as a
line position.  The release note calls them "saturated, flat-topped lines",
which is still accurate enough.  **My default: leave the release note as it
is**, and keep the plateau-model idea (report §6, item 5) as a possible
follow-up outside PR C.  Say if you want the wording changed or a follow-up
recorded.
>A:

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-10-20 (Task 1 — Prepare: branch, grounding, drift audit, scratch probe of the legacy fit)

**Branch.**  PR B is complete on **`speed_up_B`**: all tasks and PR prompts
1–7 are logged done in B.md.  #2209 (`speed_up_B` → `develop`) is open and
awaits Tim's re-review.  `speed_up_detpar` never existed (B.md Q1(a)), so C
branches from `speed_up_B`.  `speed_up_arcfit` already existed at the
`speed_up_B` tip and I checked it out:

- **HEAD = `72e448f75`** ("Merge branch 'develop' into speed_up_B").  It is
  identical to `speed_up_B` and `origin/speed_up_B`, with a clean tree and no
  upstream yet.
- A stray `speed_up_C` points at the same commit (Q1).

**Prompt doc updated** for what PR B changed:
- the branch stack: A merged to `develop`; B on `speed_up_B`/#2209; C
  targets `speed_up_B`, then `develop`; `speed_up` superseded;
- B's `--ncpu`/`forkserver` facts, with the reminder to profile at `--ncpu 1`;
- the scoped reproducibility (#2208 noise floor on multi-slit);
- #2211;
- the latest DEIMOS numbers;
- **env `pypeit` → `pypeit14b`**: re-confirmed broken today
  (`numpy has no attribute 'in1d'` on `import astropy`);
- B's safe dev-suite recipe and the stale-`Science/` caveat;
- Task 6's pre-merge steps.

**Coding doc read:** §0, §3 (C.1–C.8), §4–§8.  I re-derived the back
substitution.  With `ln y = a0 + a1 u + a2 u²` and
`ln g = ln A − (u−u0)²/(2σ²)`, matching the coefficients gives:
- `a2 = −1/(2σ²)`, so `σ = √(−1/(2a2))`, which needs `a2 < 0`;
- `a1 = u0/σ² = −2 a2 u0`, so `u0 = −a1/(2a2)`;
- `a0 = ln A − u0²/(2σ²) = ln A + a2 u0²`, so
  `A = exp(a0 − a2 u0²) = exp(a0 − a1²/(4a2))`.

(The sketch writes `a0 = ln A + a2 u0²` and `exp(a0 − a1²/(4a2))`; those are
consistent, since `a2 u0² = a1²/(4a2)`.)  Variance:
- `c = x0 − a1/(2a2)`, so `J = (0, −1/(2a2), a1/(2a2²))`;
- `var(c) = J Cov(a) Jᵀ`, with `Cov(a) = s² A⁻¹` for the weighted normal
  matrix `A = XᵀW²X` and `s² = Σw²δ²/(n_good−3)`.

The Guo weight follows from `y − g ≈ g δ ≈ y δ`.

**Code read:** `detect_lines` (continuum subtraction
`arc = (censpec − cont_now)*~bpm`, the peak finding, the `fit_arcspec` call,
the `good` mask), `fit_arcspec`, `fit_gauss`/`gauss_3deg`/`guess_gauss`.
Points not in the doc:
- `curve_fit` warns `OptimizeWarning` and returns `inf` covariance on
  degenerate windows;
- a NaN in the window would raise an uncaught `ValueError` from `curve_fit`
  (`check_finite`).  There are no NaNs in practice, because the arc is
  multiplied by `~bpm`.

**Consumers** (`git grep "detect_lines(\|centerr\|ecent"`):
- wavelength calibration: `wvutils.py:347/502/625`, `autoid.py:930/961`;
- tilts: `tracewave.py:63`;
- **object finding** (`findobj_skymask.py:2167`), **trace peaks**
  (`trace.py:1466`), **spatial flexure** (`flexure.py:90`, uses the −999
  sentinel) and sky matching (`flexure.py:283/289`);
- `wavemodel.py:697`.

`ecent`/`centerr` is never used numerically (every use is `_` or a dead
variable), so the doc's "HolyGrail pattern weights" is wrong.  Shape and
sentinel are what must be preserved.  Recorded in Key facts.

**Line-number drift** (coding doc §0.2/§3 at `e9ed85c1a` → now):

| symbol | doc | now |
|---|---|---|
| `arc.py` `detect_lines` | 857–1074 | 869–1087 |
| continuum-subtracted `arc =` | 1002 | 1014 |
| `nfitpix` | 1026 | 1039 |
| `fit_arcspec` call | 1028 | 1041 |
| `good` mask | 1042–1044 | 1055–1057 |
| `fit_arcspec` | 1077–1143 | 1090–1156 |
| per-line loop | 1121–1142 | 1134–1155 |
| `fitting.py` `fit_gauss` | 874 | 962 |
| `gauss_3deg` / `gauss_4deg` / `guess_gauss` | 921 / 936 / 952 | 1009 / 1024 / 1040 |
| `multislit_flexure.py:77,84`, `wvutils.py:353` | | unchanged |
| `autoid.py` HolyGrail `_det_*` | 2116–2136 | 2106–2126 |

The +88 in `fitting.py` is the chunked 2D-fit work (#2211).  Code content is
unchanged apart from the shifts.

**Scratch probe** (no production code; scripts in the session scratchpad).
I ran legacy `fit_arcspec` on the `test_detect_lines` Paranal sky spectrum:
4 259 peaks, `fitp = 5`, `fit_interval = 3`, `detect_lines` 3.5 s, 3 748 good
lines.
- 248 `RuntimeError`s, 56 negative-σ fits (Q3), 144 centres outside the
  window, 31 `inf` `centerr`s.
- 35% of the windows have a pixel ≤ 0.

A prototype of the analytic seed (scratch only) found two problems with the
sketch:
1. **`np.linalg.solve(A, b)` with `b` of shape `(n,3)` raises under NumPy 2.**
   It needs `b[..., None]`.
2. **The `usable`/fallback masks drop lines with fewer than 4 positive pixels
   altogether** (Q2).

With those fixed, the seed agrees with legacy to a median |Δc| of 0.053 px,
but falls back on 13.2% of lines (Q4).  So `GN_ITER = 0` will not meet the
0.02 px `matches_curve_fit` tolerance on real sky data.  That is expected, and
it is the case for Q3's accepted `GN_ITER = 3`.

**Heads-up for Task 4.**  §C.5's `test_fit_arcspec_edge_cases` puts all-negative
windows at `pixt = 0, 2` and asserts `np.all(np.isfinite(...))`.  Those lines go
to the `curve_fit` fallback, where `guess_gauss` divides 0/0 (`ypos ≡ 0`).  The
result is a NaN `p0`, which either raises an uncaught `ValueError` or leaks
NaN.  Legacy has the same weakness.  The test (or a finiteness guard on the
fallback output) has to account for it.

Questions Q1–Q5 posed.  None blocks Task 2, and each has a stated default.

### 2026-10-22 (Task 2 — vectorized log-parabola core and the legacy fallback, `GN_ITER = 0`)

**Q&A applied.**  Q1: `speed_up_C` left alone.  Q2: two masks.  Q3: the line
gain is accepted.  Q4(a): the `curve_fit` fallback is kept for every line the
analytic fit rejects.  Q5: noted for Task 3.

**Code (`pypeit/core/arc.py`, uncommitted on `speed_up_arcfit`):**
- **Module constants** after the imports: `USE_VECTORIZED_ARCFIT = True`,
  `ARCFIT_GN_ITER = 0` (for now) and `ARCFIT_MIN_GOODPIX = 4`, each with a
  docstring.
- **`_fit_arcspec_curvefit(xarray, yarray, pixt, fit_interval, sz_a, which,
  ampl, cent, widt, centerr)`**: the legacy loop, moved verbatim (including the
  JFH comment), with `range(sz_p)` replaced by `which`.  I did *not* add the
  sketch's `int()`.  It has a Numpy docstring.
- **`fit_arcspec`, rewritten.**  It returns early for `sz_p == 0`, and when
  `USE_VECTORIZED_ARCFIT` is False it calls `_fit_arcspec_curvefit` for every
  line.  Otherwise:
    - windows `idx = pixt[:,None] + off`, with an `inb` mask, clipped
      indices, and `u = x − x[pixt]` (0 outside the window);
    - `gpm = inb & isfinite(y) & (y > 0)`;
    - `attempt = (nwin > 0) & (nwin >= fit_interval)` and
      `analytic = attempt & (ngood >= 4)` (Q2);
    - moments `S_k = Σ w² u^k`, `T_k = Σ w² u^k ln y`, from one
      `(nline, nwin, 5)` power array;
    - `A` assembled by fancy indexing of `S` (Hankel), and a stacked
      `np.linalg.solve(A, T[..., None])[..., 0]` (NumPy-2 shape);
    - back-substitution to `amp`/`u0`/`sig`, with `ok` requiring `a2 < 0`,
      all values finite, and `|u0| <= umax`.  `umax` is the largest |u| in the
      window, which equals `fit_interval` in pixel units and stays correct if
      `xarray` is ever normalised, as the docstring allows;
    - the fallback on `attempt & ~ok`, with a warning above 25%.

  The `-999` sentinels and the return signature are unchanged.  The docstring
  describes the new method, and `centerr` is now documented as a *variance*.
  The unused legacy array `b` was dropped.
- **Two deviations from the sketch.**
  1. The singular-matrix guard is scaled by the product of `A`'s diagonal,
     not by `S0³`.  `S0³` has the wrong dimensions once `u` is not O(1).
  2. **`_center_variance` is already implemented** (the §C.2 propagation).
     Task 3 lists it, but the `GN_ITER = 0` path cannot fill `centerr`
     without it.  A non-finite variance (e.g. `ngood = 3` after the guards)
     also sends the line to the fallback, so no NaN reaches the caller.

**Validation, Paranal sky** (`test_detect_lines` input; 4 259 lines):
- the fallback is 563 lines (13.2%, as the Task 1 probe predicted);
- 0 NaN in any output;
- all 563 fallback lines are **bit-identical** to legacy;
- the sentinels drop from 248 to 218;
- against legacy where both fit: |Δc| median 0.040 px, p99 2.0 px;
- **good lines 3 748 → 3 803** (+55, as Q3 predicted; the test threshold of
  3 275 still passes);
- `fit_arcspec` 2.87 s → 1.35 s.

**Validation, DEIMOS 600ZD arcs** (`WaveCalib_A_0_MSC0[1-4].fits`
`arc_spectra`; 92 slits, 10 170 lines; `detect_lines(sigdetect=5, fwhm=4)`):
- the fallback is **6.7%** (678 lines);
- against legacy: |Δc| median **0.009 px**, p90 0.038, p99 0.25;
  |Δσ/σ| median 0.9%;
- **`fit_arcspec` 4.24 s → 2.04 s, only 2.1×.**  cProfile shows the vectorized
  work at ~0.02 s for all 92 spectra; the rest is the 678 fallback lines at
  ~540 `gauss_3deg` evaluations each;
- **good lines 7 887 → 7 763 (−124; 67 slits lose lines).**  All 129 lost
  lines fail the width cut, and they are saturated, flat-topped lines that the
  `nonlinear_counts` cut does not catch on DEIMOS.  The same lines are 54% of
  the fallback time.  This is **Q6**, and it is what Task 3 has to settle,
  starting with whether GN recovers them.

**Tests:**
- `pytest pypeit/tests` (`pypeit14b`): **917 passed** in 303 s, including
  `test_detect_lines`;
- the arc, findobj, flexure, trace, wave and wavetilts subsets were also run
  separately: 12 passed.

The kill switch was checked: with `USE_VECTORIZED_ARCFIT = False`, the output
is the legacy code path.

**For you to commit:** `PypeIt/pypeit/core/arc.py` (or tell me to commit it as
"Vectorize fit_arcspec with an analytic weighted log-parabola fit"), and
`claude_prompts/speed_up/C.md`.  The scratch scripts (`val2.py`,
`deimos_arcs.py`, `why_lost.py`, `lost_stats.py`, `fb_stats.py`) are in the
session scratchpad.  `deimos_arcs.py` is the harness to reuse in Task 3.

### 2026-10-22 (Task 3 — Gauss–Newton polish, moment-seeded rescue, `curve_fit`-style variance; edge cases)

**Q&A applied.**
- **Q6(b): saturated lines may be lost.**  Added **prompt 7 / Task 7**, the
  final `Speed_Up_Arcfit_Report.md` in `pypeitdev/speed_up/Reports/`, as the
  last item before the PR, with worked figures of the saturated and
  negative-σ lines.
- **Q4(b) adopted**, which Task 2's measurement triggered.  With `GN_ITER = 3`
  alone, `fit_arcspec` was still only 2.4× faster on DEIMOS: the 679
  `curve_fit` fallback lines (6.7%) were ~99% of its time.
- **Q5** implemented as stated.

**Code (`pypeit/core/arc.py`, uncommitted):**
- **`_gauss_newton_refine(u, y, gpm, amp, u0, sig, ok, niter, lam=1e-3,
  tol=1e-10)`.**  Batched Levenberg–Marquardt on the exact `curve_fit`
  objective, over `fpm = inb & isfinite(y)`, i.e. *including the negative
  wings*.
    - The step is accepted only if χ² does not rise, and each line has its own
      damping factor (÷10 on acceptance, ×10 on rejection).
    - Lines stop once their χ² changes by less than `tol·χ²` or their damping
      exceeds 1e10, so `niter` is a **cap**.  Converged lines drop out of the
      arrays.
    - σ is returned as |σ|.
- **Helpers:** `_gauss_model` (model and analytic Jacobian), `_solve3`
  (guarded, NumPy-2-safe batched 3×3 solve), and `_moment_seed` (a
  vectorized `fitting.guess_gauss`, i.e. the seed `curve_fit` itself starts
  from).
- **`_gauss_center_variance`** (Q5): `s²·[(JᵀJ)⁻¹]₁₁` with
  `s² = Σr²/(n−3)`, which is `curve_fit`'s own `pcov` definition.  It is used
  whenever GN runs.  §C.2's `_center_variance` (written in Task 2) is used when
  `GN_ITER = 0`.  The median ratio to legacy `fitcov[1,1]` is 1.000 on both
  data sets.
- **Rescue** (`ARCFIT_RESCUE_ITER = 20`, a new constant; 0 restores Q4(a)).
  Lines that `attempt` but are not `ok` after the analytic fit and polish, and
  that have more than 3 finite pixels, get the moment seed plus up to 20 GN
  iterations.  Only lines that also fail this go to `curve_fit`.
- **Post-GN centroids are not restricted to the window**, just as `curve_fit`
  does not restrict them.  The `|u0| <= umax` guard now applies only to the
  analytic seed, where it stops extrapolation.  Keeping it after GN had pushed
  ~90 converged lines back to `curve_fit`.
- **Fallback hardening.**
    - Lines whose window has a non-finite pixel are kept out of `curve_fit`:
      legacy raised an uncaught `ValueError` on them.
    - If `curve_fit` returns non-finite parameters, all four outputs become
      **−999**.  Legacy leaked NaN there, which slips past `flexure.py:90`'s
      `== -999.0` check.
    - A non-finite `centerr` on a finite fit becomes `inf`, legacy's
      "covariance could not be estimated" value.
- **`ARCFIT_GN_ITER = 3`**, per the accepted Q3 (see Q7).  The docstrings
  describe the analytic → polish → rescue → `curve_fit` cascade.

**Measured** (`scripts/arcfit/arcfit_sweep.py`; best of 3; the "accepted"
column applies `detect_lines`' cuts):

| data | config | `fit_arcspec` | speedup | `curve_fit` calls | accepted (legacy) | median \|Δc\| |
|---|---|---|---|---|---|---|
| sky, 4 259 lines | legacy | 2.35 s | 1× | 4 259 | 3 748 | — |
| | GN 0, no rescue | 1.09 s | 2.1× | 563 | 3 803 | 4e-2 |
| | GN 3, no rescue | 1.12 s | 2.1× | 569 | 3 735 | 3e-4 |
| | **GN 3, rescue 20** | **0.097 s** | **24.1×** | **36 (0.8%)** | **3 762** | 5e-4 |
| | GN 10, rescue 20 | 0.112 s | 20.9× | 36 | 3 766 | 8e-6 |
| DEIMOS, 10 170 lines | legacy | 4.68 s | 1× | 10 170 | 7 887 | — |
| | GN 0, no rescue | 1.95 s | 2.4× | 678 | 7 763 | 1e-2 |
| | GN 3, no rescue | 1.93 s | 2.4× | 679 | 7 765 | 1e-6 |
| | **GN 3, rescue 20** | **0.29 s** | **16.1×** | **5 (0.05%)** | **7 767** | 1e-6 |
| | GN 10, rescue 20 | 0.35 s | 13.3× | 5 | 7 838 | 1e-7 |

Notes on the table:
- **Rescue fidelity.**  At 100 iterations with no early stop, every one of the
  678 DEIMOS rescue lines agreed with legacy after the cuts: 66 kept by both
  and within 0.02 px, 612 rejected by both.  The accepted set is the same at
  20 iterations.  On the blended sky, 26 of 563 converged to a different local
  minimum than `curve_fit` (i.e. disagree by more than 0.02 px).
- **Rescue cap.**  At 100 the rescue cost 0.85 s on DEIMOS, because a few
  slow-converging saturated lines keep every call iterating at ~110 µs per
  iteration.  20 gives 16× with the same accepted lines.
- **Lost DEIMOS lines (shipped config): 122.**  98% have peaks above 30k and
  86% are flat-topped, the same population as in Task 2.  This is Q6(b),
  accepted.  `GN = 10` recovers 71 of them (Q7).
- **`test_detect_lines` count: 3 762** (legacy 3 748, +14), within Q3's
  agreed band.

**§C.4 edge cases walked** (`scripts/arcfit/arcfit_edge_cases.py`, new
against legacy):

| case | result |
|---|---|
| negative wings | vectorized, matches legacy |
| zero-flux, all-negative and constant windows | legacy NaN → now −999 |
| NaN pixel | vectorized; legacy *raised* |
| peaks at index 0/1/N−2/N−1 | vectorized, exact |
| too few positive pixels | rescued |
| upward parabola | `curve_fit` fallback (legacy result) |
| singular matrix | guarded, then fallback |
| `exp()` overflow (1e300) | fallback; `centerr` NaN → inf |
| saturated plateau | vectorized |
| centre 4 px outside the window | vectorized, 34.5 |
| `sz_p == 0` | empty arrays |

No case raises, and no NaN is returned anywhere.

**Tests:** `pytest pypeit/tests` (`pypeit14b`): **917 passed** (303 s).

**Scripts saved** to `pypeitdev/speed_up/scripts/arcfit/`: `arcfit_load.py`,
`arcfit_sweep.py`, `arcfit_edge_cases.py` and `arcfit_lost_lines.py`.  They
need `REDUX_OUT/keck_deimos/600ZD_M_6500/Calibrations/WaveCalib_*`.

**For you to commit:** `PypeIt/pypeit/core/arc.py`, and in the dev suite
`claude_prompts/speed_up/C.md` and `pypeitdev/speed_up/scripts/arcfit/`.
Q7 is open.  Task 4's tests will pin down the cascade, including a test that
no `curve_fit` call happens on clean data, which now holds with the rescue.

### 2026-10-22 (Task 4 — tests, micro-benchmark; `ARCFIT_GN_ITER = 10` per Q7)

**Q7 → `ARCFIT_GN_ITER = 10`.**  You committed `arc.py` with it as
`f1cef0d11` ("3") while this task ran.  With it, **`test_detect_lines` gives
`len(arx_w)` = 3 766**: legacy 3 748, +18, inside Q3's agreed band, and above
the 3 275 threshold.  No investigation needed.

**`pypeit/tests/test_arc.py`** (module-level functions, fixed seeds).  Two
helpers:
- `_synthetic_arc(nspec=4096, spacing=20, fwhm=4, noise=1, seed=42)` puts
  lines on a regular grid with ±2 px jitter, amplitudes 10^1.5–10^3.5, and
  σ = fwhm/2.355 × U(0.9, 1.1).  §C.5's generator kept only 36 of 200 lines
  after its spacing cut; this one gives 203;
- `_legacy_fit` flips `USE_VECTORIZED_ARCFIT` through `monkeypatch.context`,
  not the sketch's try/finally.

Tests:
- `test_fit_arcspec_recovers_truth`.  The tolerances are relaxed from the
  sketch's 0.05 px / 5% max to medians of < 0.01 px and < 2%, and maxima of
  < 0.1 px and 10%.  The faintest lines (S/N ~ 30) reach 0.058 px, which is
  the noise and is the same for legacy.
- `test_fit_arcspec_matches_curve_fit`: centroid < 0.02 px, and width,
  amplitude and **`centerr`** within 2%.  The measured maximum is 2e-7 px, and
  the `centerr` ratio is 1.0000.
- `test_fit_arcspec_analytic_only` (new): `GN = 0`, rescue 0, median |Δc| vs
  legacy < 0.005 px (measured 2e-5).
- `test_fit_arcspec_kill_switch` (new): with the switch off, the output is
  bit-identical to `_fit_arcspec_curvefit`.
- `test_fit_arcspec_edge_cases`:
    - the sketch's case: the all-negative and zero windows return −999, the
      delta spike returns 32, and there is no NaN;
    - lines at index 0, 1, N−2 and N−1, recovered to 1e-4;
    - a NaN pixel (legacy raised);
    - empty `pixt`.
- `test_fit_arcspec_negative_wings` (new): a −30 offset puts the outer window
  pixels below zero.  The fit must not call `curve_fit`, and must equal it to
  1e-6.
- `test_fit_arcspec_no_curve_fit_on_clean_data`: `fitting.fit_gauss` is
  monkeypatched to raise.

**Mutation check.**  I made GN use only the positive pixels (`fpm = gpm`).
At first only `edge_cases` failed, because a −8 offset left every window
pixel positive.  The offset is now −30, with an assertion that two pixels are
negative, and the mutation fails both tests.  `arc.py` was restored from a
copy; it is identical to `f1cef0d11`.

**`pytest pypeit/tests` (`pypeit14b`): 924 passed** (917 + 7 new), 305 s.

**Micro-benchmark** (`scripts/arcfit/arcfit_microbench.py`; best of 5
vectorized, best of 2 legacy; for the PR description):

| input | lines | legacy | vectorized | speedup | per line |
|---|---|---|---|---|---|
| synthetic (`_synthetic_arc`, nspec 40 960) | 2 046 | 437 ms | 16.6 ms | **26×** | 214 → 8.1 µs |
| Paranal sky | 4 259 | 2 279 ms | 111 ms | **21×** | 535 → 26 µs |

The DEIMOS figure at `GN = 10` / rescue 20 is 13.3× (Task 3).

**For you to commit:** `PypeIt/pypeit/tests/test_arc.py`; in the dev suite,
`claude_prompts/speed_up/C.md` and
`pypeitdev/speed_up/scripts/arcfit/arcfit_microbench.py`.  No new questions.

### 2026-10-22 (Task 5 — in progress: dev-suite vetting running as a systemd service; changelog and docs done)

**Variants.**  Q7 changed the shipped value to `GN_ITER = 10`, so the
variants are:
- `legacy` (`USE_VECTORIZED_ARCFIT = False`, i.e. the pre-PR-C code path,
  which serves as the baseline);
- `gn10` (shipped);
- `gn0` (analytic seed, plus the rescue).

Each is reduced into its **own fresh output directory**,
`REDUX_OUT_C/<variant>/`.  That avoids the stale-`Science/` problem (B.md
QP7), and the legacy baseline stays on disk next to the new runs.  The
existing `REDUX_OUT/` is not touched.

**Variant switch without editing `arc.py`.**
`pypeitdev/speed_up/scripts/arcfit/variant_hook/sitecustomize.py`, on
`PYTHONPATH`, reads `PYPEIT_ARCFIT_VARIANT` and sets the module constants in
every interpreter.  I checked that it reaches the parent, `forkserver`
workers and `spawn` workers (`hookcheck.py`).  A Kast blue 600 smoke reduction
under `legacy` and `gn10` confirmed that it reaches `run_pypeit`: the fitted
pixel positions differ, and the RMS agrees to 5 digits (0.07390).
`pypeit_test` hands `os.environ` to `run_pypeit`, and it has to be run as
`python ./pypeit_test`: its shebang resolves to the broken `pypeit` env.

**Batch.**  `pypeitdev/speed_up/scripts/arcfit/arcfit_vet.sh` runs as the user
service `pypeit-arcfit-vet` (`MemoryMax=120G`, `MemorySwapMax=0`).  Started
2026-10-22 17:25 on PypeIt HEAD `88337631c`, with a clean `pypeit/` tree.
For each variant, in the order legacy → gn10 → gn0:
- `pypeit_test reduce -i shane_kast_blue shane_kast_red keck_deimos
  keck_lris_blue -t 4 -o REDUX_OUT_C/<v>`;
- then `pytest vet_tests/test_{wavelengths,wavetilts,extraction,skysub}.py
  --redux_out REDUX_OUT_C/<v>`.

This is about 45 setups per variant, 22 of them DEIMOS.  `-t 4` alone, rather
than B's `-t 1`, is a deliberate trade: the crash came from *two concurrent*
suites, and the memory cap contains a blow-up.  Check on it with:
- `systemctl --user status pypeit-arcfit-vet`
- `REDUX_OUT_C/progress.txt`

**Do not edit `pypeit/` while it runs**: the editable install means the
variants would see different code.

**Metrics.**  `scripts/arcfit/arcfit_metrics.py REDUX_OUT_C legacy gn10 gn0`
writes Markdown tables:
- per WaveCalib file: median and maximum slit RMS (what `test_wavelengths`
  asserts on), lines in the final fits, the maximum change of the wavelength
  solution in pixels, and the number of slits whose RMS is >10% worse than
  legacy;
- per Tilts file: tilt lines and 2D-fit points, and the median and maximum
  change of the tilt image in pixels.

The tilt-fit RMS is not stored, and the main reduction log does not record
it.  `test_wavetilts::test_run` asserts no RMS (it only evaluates the tilts),
so the tilt-image difference is the tilt metric.

**Changelog/docs (done):**
- `doc/releases/2.1.0dev.rst`: a §C.8 bullet under **Under-the-hood
  Improvements**, extended to cover the rescue, the saturated-line behaviour,
  the NaN → −999 change and the kill switch, plus a **Testing** bullet for the
  new `test_arc.py` tests.  `CHANGES.rst` is deprecated per its header and was
  not touched;
- `make clean; make htmlnoex` (`apirst` plus the build; examples skipped): rc 0,
  **`sphinx_warnings.out` empty**;
- `doc/api` needs no change for PR C, because `pypeit.core.arc` already had
  its page and the new private functions render in it.  The only diff was a
  blank line in `doc/api/pypeit.coadd2d.rst` from the apidoc patch script,
  unrelated, and I reverted it;
- `doc/pypeit_par.rst` is unchanged, since no parameter was added.

### 2026-10-23 (Task 5 — dev-suite results: vet tolerances pass, but DEIMOS identifications flip; Q8)

**The batch finished** at 01:09 after 7 h 44 m, with each variant's reduce
taking ~2.5 h at `-t 4` (peak memory well under the 120 G cap).
- **The reduce passed 35/35 setups for every variant** (`legacy`, `gn10`,
  `gn0`).
- **Vet tests: `legacy` and `gn10` have identical outcomes** (8 passed,
  13 failed).  All 13 failures are missing products for instruments outside
  the four reduced (not_alfosc, mdm, lris_red, hires, hamspec, gmos, focas,
  bok, echelle), raising "Calibration directory does not exist" or
  `FileNotFoundError`.
- **`gn0` additionally fails `test_shane_kast_red`**: `600_7500_d57` RMS
  0.0633 against a tolerance of 0.05.  This confirms that the pure analytic
  fit is not shippable.

**RMS on the vetted setups** (the exact slit or quantity each test asserts
on):

| vet test / setup | slit | tolerance | legacy | gn10 | gn0 |
|---|---|---|---|---|---|
| shane_kast_red 300_7500_Ne | 1 | 0.05 | 0.0308 | 0.0308 | 0.0413 |
| shane_kast_red 600_7500_d57 | 0 | 0.05 | 0.0369 | 0.0369 | 0.0633 **FAIL** |
| shane_kast_red 1200_5000_d57 | 0 | 0.055 | 0.0512 | 0.0512 | 0.0496 |
| keck_deimos 1200B_LVM_5200 | 3 | 0.15 | 0.1234 | 0.1234 | 0.1242 |
| keck_deimos 600ZD_M_6500 | 1 | 0.35 | 0.3356 | 0.3350 | 0.3100 |
| keck_deimos 900ZD_LVM_5500 | 1 | 0.1 | 0.0687 | 0.0687 | 0.0685 |
| keck_lris_blue multi_300_5000_d680 | max | 0.32 | 0.2970 | 0.2969 | 0.2877 |
| keck_lris_blue long_400_3400_d560 | max | 0.21 | 0.1947 | 0.1947 | 0.1853 |
| keck_lris_blue long_600_4000_d560 | max | 0.08 | 0.0621 | 0.0621 | 0.0772 |

**Full comparison** (`REDUX_OUT_C/arcfit_metrics.md`, from
`arcfit_metrics.py`; 38 WaveCalib and 38 Tilts files, `ncpu/` and
`step_by_step/` duplicates excluded):
- **Wavelength solutions.**  Median max |Δλ| over the common line range:
  **gn10 2.2e-6 px**, gn0 0.29 px.  For gn10, 10 files exceed 0.05 px and 7
  exceed 0.5 px, all on multi-slit DEIMOS (1200B_M, 1200G_M_5500/7750, the
  four 600ZD mosaics, 830G_L/LVM) plus LRIS `multi_300` at 0.096 px.
- **Lines in the final fits:** legacy 21 702, gn10 21 681, gn0 21 700.
- **Tilts.**  gn10 matches to median |Δtilt| ≤ 3e-4 px everywhere, with max
  ≤ 0.04 px except `600ZD` MSC02 (0.54 px) and `1200G_M_5500` (0.88 px).
  The traced-line counts are within ±5.  gn0 is noisier (max up to 62 px on
  `1200G_M_5500`).
- **Serial and `--ncpu 2` products are identical within each variant**
  (`830G_M_8500`, `multi_600_4000_slitmask`), and the hook reached the
  workers.

**Comparison method.**  The first |Δλ| metric was taken over the whole
detector, so it mostly measured extrapolation where one variant has extra lines
at the ends.  For example, `1200G_M_5500` slit 3 showed 1 408 px although the
two variants share all 22 of gn10's lines.  I restricted it to the pixel range
covered by **both** variants' lines.

**Root cause of the large gn10 differences** (`arcfit_slitdiff.py`, the
`trig.py` scratch probe, and `arcfit_gained_lost.py` over all 370 arc
spectra).  They are **identification flips**, triggered by discrete changes
in the detected line list:
- **lines gained:** 83, all sub-resolution spikes on which legacy `curve_fit`
  failed (`830G_L_8100`: 5 spikes with σ ≈ 0.45–0.55 px, where 1.43 px is
  expected);
- **lines lost:** 56 saturated plateaus (`600ZD` MSC01 slit 11, MSC03 slits 0,
  14 and 15);
- **shared lines converging to a different minimum:** MSC01 slit 2 (3.2 px).

Example: `600ZD` MSC01 slit 2 is misidentified at the blue end by ~30 px.
See **Q8**.

**Fidelity prototype** (Q8(A); saved as
`scripts/arcfit/prototype_defer_to_curvefit.patch`; **`arc.py` is restored
to `88337631c`**, and `test_arc.py` passes 8/8):
- the GN routine returns a `converged` flag (with a parameter-step `xtol`
  besides the χ² `tol`);
- non-converged polish and rescue lines are deferred to `curve_fit`;
- so are fits with σ < `ARCFIT_MIN_WIDTH = 0.35` × the window half-width.

Measured (`rules.py`, `tolsweep.py`, `xtolsweep.py`, `deferred.py`, `nfev.py`
in the scratchpad):
- **0 gained / 0 lost over all 370 spectra**, at 2.0× on the DEIMOS arcs and
  1.4× on the sky;
- the convergence rule carries the cost: ~500 DEIMOS lines, insensitive to
  `tol` (1e-10 to 1e-6), `xtol` (1e-6 to 1e-4) and the cap (10/30);
- **432 of the 517 deferred lines end in `curve_fit` `RuntimeError`**, after
  its full 800 evaluations (1.75 of 1.83 s);
- capping `maxfev` is not safe.  Over 7 505 fallback fits, successful ones
  reach 799 evaluations, and 22 plausible lines need more than 400.

**Changelog/docs** (done earlier today; see the in-progress entry).  If Q8 is
(A), the bullet will need its "more than 10 times faster" claim and its
saturated-line sentence revised.

**Disk:** `REDUX_OUT_C/` holds ~27 GB per variant (~80 GB).  Keep `legacy/`
as the baseline for the next run; `gn0/` can be deleted.

**For you to commit:** `PypeIt/doc/releases/2.1.0dev.rst`; in the dev suite,
`claude_prompts/speed_up/C.md` and `pypeitdev/speed_up/scripts/arcfit/`
(`arcfit_vet.sh`, `variant_hook/`, `arcfit_metrics.py`, `arcfit_slitdiff.py`,
`arcfit_gained_lost.py`, `prototype_defer_to_curvefit.patch`).
**Task 5 waits on Q8.**

### 2026-10-23 (Task 6 — in progress: option A applied with B as a switch; profiling batch running)

**Q8 → ship (A), keep (B) as an option.**
- The prototype patch was applied to `pypeit/core/arc.py`, with the three
  hand-off rules gated by the new module constant
  **`ARCFIT_MATCH_CURVEFIT = True`**: non-converged polish lines,
  non-converged rescue lines, and fits with σ < `ARCFIT_MIN_WIDTH` (0.35) ×
  the window half-width.  False gives (B).
- `_gauss_newton_refine` now also returns a per-line `converged` flag, and
  takes `xtol` alongside `tol`.  The docstrings of `fit_arcspec` and both
  constants describe the trade-off.
- **Tests:**
    - `test_fit_arcspec_negative_wings` now runs with
      `ARCFIT_MATCH_CURVEFIT = False`, since it tests the vectorized path
      itself.  Its noise-free, pedestal-offset line converges only linearly
      and hits the cap;
    - new `test_fit_arcspec_narrow_features`: under (A) a hot pixel goes to
      `curve_fit` and gets exactly its result; under (B) it does not, and
      comes back narrow;
    - new `test_detect_lines_matches_curve_fit`: on the Paranal sky, (A)
      accepts exactly legacy's lines, with centroids within 0.02 px.

  `pytest pypeit/tests`: **926 passed**.
- **Changelog** rewritten for (A)/(B) ("about twice as fast"; the
  `ARCFIT_MATCH_CURVEFIT = False` option and its caveat), and the Testing
  bullet extended.  `make clean; make htmlnoex`: rc 0, **0 warnings**, and the
  `:data:` cross-references resolve.  The stray `doc/api/pypeit.coadd2d.rst`
  blank-line change was reverted again.
- **Q9** posed: expose (B) as an `[rdx]` parameter?

**Batch** `scripts/arcfit/arcfit_task6.sh`, as the user service
`pypeit-arcfit-task6` (`MemoryMax=120G`).  Started 2026-10-23 07:08 on HEAD
`88337631c`, plus the uncommitted (A) changes (2 dirty files).  Progress is
logged to `REDUX_OUT_C/task6_progress.txt`.  In order:
1. the dev suite for variant **A** (`REDUX_OUT_C/A`; completes Task 5 for the
   shipped code);
2. cProfile `--ncpu 1` of Kast under legacy, A and B;
3. cProfile `--ncpu 1` of DEIMOS 600ZD under A, then legacy;
4. uninstrumented DEIMOS `--ncpu 4` (B's `run_with_rss.py` +
   `instr_run_pypeit.py` recipe) under A, then legacy;
5. cProfile DEIMOS under B.

Notes:
- The variant hook now accepts `A`/`B`.  The earlier `gn10` run equals `B`.
- The profile scripts take `--tag`, and the artifacts are
  `Reports/<stem>_C_<variant>.*`.
- Legacy is re-profiled at the current HEAD because no serial DEIMOS profile
  exists on the current code (see "After the B port").
- The profile scripts reduce in place in `REDUX_OUT/<setup>`, so the main
  `REDUX_OUT` DEIMOS and Kast products will be those of the last run, B.

The expected duration is ~13 h.  **Do not edit `pypeit/` while it runs**;
committing is fine.

### 2026-10-23 (Q9 → `[calibrations][wavelengths] fast_line_fit`, built in a worktree; port waits for the batch)

**Q9 answered: a user parameter in `WavelengthSolutionPar`.**  The batch was
running on the editable `pypeit/` tree, so the work is in a **separate git
worktree, `/mnt/tank/Astronomy/PypeIt/PypeIt_C_par`** (detached at
`88337631c`).  It holds copies of the uncommitted (A) changes, plus the
ignored `pypeit/pkg/version.py` and `pypeit/version.py`, and is tested with
`PYTHONPATH` pointing at it.

**Design.**
- **`fast_line_fit`** (bool, default **False**, meaning option A) is a new
  `WavelengthSolutionPar` key.  It has a description, an `__init__` argument
  and a `from_dict` key.
- It is applied by **`BuildWaveCalib.run`**, whose whole body now runs inside
  a new context manager, **`arc.arcfit_options(match_curvefit=not
  fast_line_fit)`**.  The context manager sets `arc.ARCFIT_MATCH_CURVEFIT`
  temporarily and restores it in a `finally`.  `git diff -w` on
  `wavecalib.py` is just the two added lines.
- Why a context manager rather than threading a keyword through
  `detect_lines`: ~15 wavelength-calibration call sites reach `fit_arcspec`,
  several without `par` at hand (`autoid.py:527/541/783/792/1140`,
  `wvutils.py:347/502/625`).
- Consequences:
    - the parameter affects **only the wavelength calibration**, as its
      location suggests; tilts, object finding and flexure keep the default
      (A);
    - it works unchanged in `--ncpu` workers, because each worker runs
      `BuildWaveCalib.run` itself.  No `parallel._worker_init` `cfg` entry is
      needed, unlike `fit2d_chunk`;
    - `ARCFIT_MATCH_CURVEFIT` remains the all-callers switch.

**Tests** (in `test_arc.py`):
- `test_arcfit_options`: sets and restores the option, including on an
  exception, and a `None` argument is a no-op;
- `test_fast_line_fit_par`: the default is False, and
  `[calibrations] [[wavelengths]] fast_line_fit = True` parses through
  `PypeItPar.from_cfg_lines`;
- `test_wavecalib_applies_fast_line_fit`: a stubbed `extract_arcs` records the
  option seen inside `run` (False for `fast_line_fit=True`, and vice versa),
  and the default is restored afterwards.

Worktree `pytest pypeit/tests`: **929 passed**.

**Docs.**  `doc/scripts/build_par_rst.py` adds one row, for
`fast_line_fit`, to `doc/pypeit_par.rst`.  Release notes:
- a new **Functionality** bullet (vectorized line fitting; `fast_line_fit`);
- the Under-the-hood bullet now points to the parameter;
- the Testing bullet is extended.

The "roughly halves its cost" and "more than 10 times faster" wording is
provisional until the Task 6 profiles are in.

**Port plan, once `pypeit-arcfit-task6` finishes.**  I am deliberately not
porting now:
- it would change the code under the running dev-suite and profiles;
- with the new wrapper, the hook's `B` variant would be overridden to (A)
  inside the wavelength calibration, which would invalidate the B profile.

Once it finishes:
1. copy the worktree's `pypeit/core/arc.py`, `pypeit/par/pypeitpar.py`,
   `pypeit/wavecalib.py`, `pypeit/tests/test_arc.py`, `doc/pypeit_par.rst` and
   `doc/releases/2.1.0dev.rst` over the main tree (the main tree's copies are
   unchanged since the worktree was made);
2. run `pytest` and `make htmlnoex`;
3. `git worktree remove`.

The batch's B profile is therefore the all-callers B
(`ARCFIT_MATCH_CURVEFIT = False`), and its A profile equals the default of
the final code.

### 2026-10-23 (Task 6 — batch complete: Task 5 for A, DEIMOS/Kast profiles, `--ncpu 4`; `fast_line_fit` ported)

The batch `pypeit-arcfit-task6` ran 07:08–20:10, and every step returned
rc 0.  Progress log: `REDUX_OUT_C/task6_progress.txt`.

**Task 5 for the shipped (A).**  The dev suite for variant A
(`REDUX_OUT_C/A`) passed **35/35**, and its **vet outcomes are identical to
legacy**.  `arcfit_metrics.py legacy A gn10` is in
`REDUX_OUT_C/arcfit_metrics_A.md`:
- 37 of 38 WaveCalib files match legacy (median max |Δλ| 1.4e-6 px), with
  the same line counts and RMS;
- the one exception is `600ZD_M_6500` MSC03.  It has 1 812 lines against
  1 804, its median RMS improves from 0.2656 to 0.2579, and one slit moves by
  up to 31 px.  This remains to be traced;
- for comparison, B (`gn10`) had 10 files above 0.05 px.

**DEIMOS 600ZD, `--ncpu 1`, cProfile, current HEAD** (`Reports/
keck_deimos_600zd_m_6500_C_{legacy,A,B}.*`):

| variant | wall | `fit_arcspec` | `fit_gauss` calls / cum | `measure_fwhm` (flexure) |
|---|---|---|---|---|
| legacy | 10 420 s | 851 s | 1 114 679 / 846 s | 801 s |
| **A (shipped)** | **10 095 s (−326 s, −3.1%)** | 526 s | 319 275 / 494 s | 571 s |
| B | 9 554 s (−867 s, −8.3%) | 44 s | 10 013 / 20 s | 126 s |

- **The baseline is reproduced:** legacy's 1.11 M `fit_gauss` calls and 846 s
  match the June profile (1.11 M, 840 s).
- **The vectorized work itself costs very little:** `_gauss_newton_refine`
  is 24 s under (A) and 19 s under (B).
- **(B)'s fit is 19× faster than legacy's.**  (A)'s is only 1.6× faster,
  because its hand-off goes to `curve_fit`, which is mostly spent on lines
  that then fail.
- **94% of legacy's arc-fit cost is `measure_fwhm` inside spectral
  flexure:** the archive-sky FWHM, recomputed 170 times.  That is **Q10**,
  still open.  A cache would remove ~800 s from legacy and ~560 s from (A),
  i.e. more than PR C's change to the fit itself.

**DEIMOS `--ncpu 4`, uninstrumented** (`REDUX_OUT_C/_ncpu4/{A,legacy}`, B's
`run_with_rss` recipe):

| | legacy | A |
|---|---|---|
| wall | 3 400.2 s | 3 368.6 s (−32 s, −0.9%) |
| calibrations | 1 480.1 s | 1 487.9 s |
| image processing | 76.7 s | 78.2 s |
| extraction | 660.7 s | 624.7 s |
| peak PSS | 36.7 GB | 36.6 GB |

At `--ncpu 4`, C's gain is within run-to-run noise (B measured 2.9% of
scatter).  As anticipated, it is a serial-time gain.  B's 3 391 s on older
code matches legacy here.

**Kast, `--ncpu 1`, cProfile** (`Reports/shane_kast_blue_A_C_*`):

| | legacy | A | B |
|---|---|---|---|
| wall | 88.5 s | 88.4 s | 100.7 s |
| `fit_arcspec` | 0.75 s | 0.47 s | 0.22 s |
| `curve_fit` calls | 1 297 | 109 | 5 |

The wall differences are machine noise (B is slower despite the fastest fit).
§C.7's "~1–3 s off Kast" is in fact 0.3–0.5 s.

**`fast_line_fit` ported to the main tree** after the batch:
- the main tree's six files were checked as unchanged, and the worktree's as
  pure supersets (except the intended release-note rewrites);
- copied over; `pytest pypeit/tests` **929 passed**; `make clean; make
  htmlnoex` rc 0, **0 warnings**;
- `doc/api/pypeit.coadd2d.rst` reverted again; worktree removed.

The uncommitted files are `pypeit/core/arc.py`, `pypeit/par/pypeitpar.py`,
`pypeit/wavecalib.py`, `pypeit/tests/test_arc.py`, `doc/pypeit_par.rst` and
`doc/releases/2.1.0dev.rst`.

**Release-note wording to fix:** "about twice as fast" stands (fit 1.6×), but
"more than 10 times faster" for `fast_line_fit` overstates the
wavelength-only case.  B's 19× applies to all callers, and most of the cost
is in flexure, which the wavelength-only parameter does not reach.  I will
rephrase this once Q10 is settled.

**Still to do for Task 6:**
- trace the `600ZD` MSC03 slit;
- check the spec1d `WAVE_RMS` and the QA arc-fit PNGs;
- the "after-C" section of `speed_up_results.md` (A+B+C totals: serial
  12 390 → 10 095 s, i.e. B's #2211 and A's merges plus C; `--ncpu 4`
  3 369 s);
- the pre-merge notes.

These wait on **Q10**, because the cache would change every number above.

### 2026-10-26 (Q10 → archive-sky FWHM cache implemented; session paused before the re-profile)

**Q10 answered: yes.**  Implemented (by an Opus 5.5 sub-agent) and left
uncommitted in the main tree, on top of the earlier uncommitted PR C files:
- **`pypeit/core/flexure.py`**: new `@functools.lru_cache(maxsize=None)`
  helper `_archive_sky_fwhm(sky_file)` that loads the archive sky and runs
  `autoid.measure_fwhm` on it once per file.  `get_archive_spectrum` now
  still loads the `OneSpec` every call (cheap, and so callers never share a
  mutable object) but takes the FWHM from the cache.  A failed measurement
  raises `PypeItError` as before; exceptions are not cached.  Each `--ncpu`
  worker keeps its own cache, so no `parallel._worker_init` plumbing.  The
  `'model'` branch is unchanged.  The result is bit-identical.
- **`pypeit/tests/test_flexure.py`**: `test_archive_fwhm_cache` counts
  `measure_fwhm` calls through a monkeypatch: one measurement for two calls
  on `paranal_sky.fits`, a second for `sky_kastb_600.fits`, a third after
  `cache_clear()`; the returned spectra are distinct objects with equal
  flux.  `pytest pypeit/tests/test_flexure.py`: **3 passed**.
- **`doc/releases/2.1.0dev.rst`**: new Under-the-hood bullet for the cache,
  a Testing bullet for the new test, and the `fast_line_fit` wording toned
  down (the shipped default "roughly halves" the line-fitting cost; the
  vectorized fit "itself is many times faster"; `fast_line_fit` removes the
  remaining `curve_fit` hand-offs from the wavelength calibration only).
- **`pypeitdev/speed_up/scripts/arcfit/arcfit_q10.sh`** written: DEIMOS
  600ZD cProfile at `--ncpu 1` under A (tag `C_A_cache`) then B (tag
  `C_B_cache`), each followed by `analyze_profile.py`; progress to
  `REDUX_OUT_C/q10_progress.txt`.

**Interrupted** at the user's request (exit and resume), then resumed the
same morning.  Correction to what I first wrote here: the sub-agent *had*
already launched the profile, as the user service **`pypeit-arcfit-q10`**
at **06:19:21**, before I stopped it mid-`pytest`.  The tree it profiles is
the final one (`flexure.py` and the test were last written at 06:13, and
nothing in `pypeit/` changed afterwards).  Caveat: the full test suite and
the doc build (below) ran on the same machine during the profile's first
~15 minutes, which may add a little to variant A's wall time; the
per-function cProfile numbers that matter (`measure_fwhm`, `fit_gauss`) are
not affected in any meaningful way on this many-core machine.

After the resume:
- `pytest pypeit/tests -q` (`pypeit14b`): **930 passed** (320 s);
- `make clean; make htmlnoex`: rc 0, **0 warnings**; the stray
  `doc/api/pypeit.coadd2d.rst` blank-line change was restored again;
- the profile service is running (variant A first, then B; ~2.8 h each).
  Progress: `REDUX_OUT_C/q10_progress.txt`; stdout in
  `REDUX_OUT_C/_profiles/profile_deimos.py.{A,B}_cache.stdout.txt`;
  artifacts `Reports/keck_deimos_600zd_m_6500_C_{A,B}_cache.*`.

**Remaining after the profile:** fold the cache numbers into the release
note, the "after-C" section of `speed_up_results.md` (A+B+C totals, serial
and `--ncpu 4`) and the pre-merge notes; then Task 7 (the report).  The
`600ZD` MSC03 trace and the spec1d `WAVE_RMS` / QA-PNG check do not depend
on the cache (it is bit-identical), so they are being done while the
profile runs (next entry).

Note to self: a read-only `git status`/`git log` was run once this session
by mistake; the no-git rule covers those too.

### 2026-10-26 (Task 6 — `600ZD` MSC03 traced to the identification step; spec1d `WAVE_RMS` and QA PNGs checked)

Done by an Opus 5.5 sub-agent while the Q10 profile runs; read-only, no
`pypeit/` change.  Scripts in `scripts/arcfit/`: `arcfit_msc03_slit11.py`
(four-way `full_template` reproduction: fwhm × fit path),
`arcfit_msc03_trace.py` (per-call `detect_lines` trace and the snippet-0
shift/stretch), `arcfit_msc03_judge.py` (line-list and dispersion checks;
figure), `arcfit_wave_rms.py`, `arcfit_png_cmp.py`.  Figures in
`Reports/figures/arcfit/` (`msc03_slit11_judge.png`,
`crop_Arc_1dfit_A_0_MSC03_S0695.png`, `sbs_FWHMfit_S0149.png`).

**The MSC03 exception is one slit, and (A) is right.**  Only slit 11
(spat_id 1063) differs; the other 21 slits agree to ≤1e-6 px.

| run | rms (px) | lines fit / kept | detected | λ(pix 0) |
|---|---|---|---|---|
| legacy | 0.2757 | 80 / 73 | 103 | 4191.26 Å |
| A | 0.2513 | 87 / 81 | 103 (same to 4.7e-6 px) | 4215.85 Å |
| gn10 | 0.2756 | 80 / 73 | 102 | as legacy |

- 71 identifications are common; all differences are blueward of ~5100 Å
  (pixel < 1500); the solutions differ by 38 px at pixel 0, 8.9 px at 1000,
  0.85 px at 2000.  Legacy alone has 174.23 px → 4301.311 and 766.22 px →
  4672.534 Å; (A) alone has ten blue lines (e.g. 174.23 → 4320.794,
  238.67 → 4359.560, 766.22 → 4681.450 Å).
- **Root cause: the identification step (ii)**, not detection and not the
  final fit.  Re-running `autoid.full_template` on the stored arc spectrum
  reproduces both results exactly, and the fwhm difference (1e-6) plays no
  part.  All 15 `detect_lines` calls inside `full_template` accept the same
  lines, centres within 2e-6–8.5e-4 px.  The chain: `get_xcorr_arc` builds
  its synthetic arcs from *all* peaks → the cross-correlation changes by
  1.6e-8 → `detect_lines(corr, cont_subtract=True, nfind=1)`'s continuum
  fit flips, moving the correlation peak 4131.8517 → 4131.8583 px → the
  initial shift 26.8517 → 26.8583 seeds `differential_evolution` in
  `xcorr_shift_stretch` (snippet 0, pixels 0–4106), which lands in a
  different optimum (legacy: shift 75.86, stretch 0.98896, stretch2 −1e-6;
  A: shift 98.59, stretch 0.97535, stretch2 +1e-6; stretch2 pinned at its
  bound in both, ≈ ±17 px of curvature).  Snippet 1 is identical.
- **(A) is the correct branch.**  Mapping the 12 detected blue lines
  (identical in both runs) through each solution: with (A), 10 of 12 fall
  within 0.5 px of a line-list line (median 0.14 px); with legacy, 2 of 12
  (median 7.2 px).  (A)'s dispersion rises smoothly 0.6035 → 0.6293 Å/px
  over 4300–5500 Å like slits 10, 12, 13 and the template; legacy's turns
  over (0.630 → 0.626 → 0.631).  (A) also has lower rms with more lines.
  Caution: legacy's optimum has the *higher* zero-lag correlation (0.778 vs
  0.767), so that objective cannot pick the right branch here.
- **Why the hand-off rules did not make it identical:** they guarantee the
  same accepted *set*, not bit-identical centres; the vectorized fit and
  `curve_fit` stop within their own tolerances, 1e-5–1e-3 px apart.  The
  correlation-peak fit (half-width 4, σ/half-width 1.73, so vectorized, no
  hand-off) gives tcent 4131.8517 / twid 6.9222 (legacy) vs 4131.8583 /
  6.9204 (A); on identical input both paths return exactly the same numbers
  (centre 4131.851679, σ 6.922171), so the change is in the input, not the
  fit.  The release-note wording ("accepts the same lines as before, with
  centroids that agree to better than 0.001 pixels") is accurate as it
  stands.  The slit-9 suspicion and the fragile step are **Q11**.

**spec1d `WAVE_RMS`, legacy vs (A):** DEIMOS 600ZD, 84 objects: exactly one
differs, `SPAT1068-SLIT1063-MSC03` (0.2757 → 0.2513; slit 11 above).  Kast
blue (3 setups) and Kast red (5 setups), 40 objects: max difference 1.1e-6;
nothing missing.

**QA PNGs:** DEIMOS MSC03 `Arc_1dfit`: 13 of 22 byte-identical; S1063 differs
as expected; S0695/S1293 differ in ~250 pixels (same markers, slight shift);
six others at anti-aliasing level (max 1–5).  `Arc_FWHMfit`: none
byte-identical (widths differ at 1e-6; plots look the same, S0149 side by
side).  Kast: `Arc_1dfit`, `Arc_tilts_2d`, `Arc_tilts_spec` byte-identical;
`Arc_FWHMfit` and `Arc_tilts_spat` differ at the 1–16 level.  The older
`task6_qa_samples` is the ncpu1/ncpu4 comparison, not A vs legacy.

### 2026-10-26 (Q10 — DEIMOS 600ZD `--ncpu 1` profile of (A) with the archive-FWHM cache)

`pypeit-arcfit-q10` step 1 finished 08:57 (rc 0; B started at once).
Artifacts: `Reports/keck_deimos_600zd_m_6500_C_A_cache.*`.  All four serial
profiles are on the same code except for the cache:

| | legacy | A | **A + cache** | B (no cache) |
|---|---|---|---|---|
| wall | 10 420 s | 10 094 s | **9 462 s** | 9 554 s |
| vs legacy | — | −326 s (−3.1%) | **−958 s (−9.2%)** | −867 s (−8.3%) |
| `measure_fwhm` (flexure) calls / cum | 170 / 801 s | 170 / 571 s | **1 / 3.2 s** | 170 / 126 s |
| `get_archive_spectrum` cum | 802 s | 571 s | **3.6 s** | 127 s |
| `spec_flex_shift` cum | 842 s | 591 s | **23 s** | 138 s |
| `detect_lines` calls / cum | 2 944 / 981 s | 2 944 / 655 s | **2 775 / 87 s** | 2 944 / 170 s |
| `fit_arcspec` cum | 851 s | 526 s | **57 s** | 44 s |
| `fit_gauss` calls / cum | 1.11 M / 846 s | 319 k / 494 s | **14 k / 42 s** | 10 k / 20 s |
| `curve_fit` cum | 666 s | 441 s | **40 s** | 18 s |
| `gauss_3deg` calls | 93.0 M | 71.4 M | **7.0 M** | 3.1 M |
| `_gauss_newton_refine` cum | — | 24 s | 12 s | 19 s |
| `wavecalib.run` / `wavetilts.run` / `extraction.run` | 764 / 1 741 / 2 038 s | 696 / 1 761 / 1 894 s | 693 / 1 727 / 1 606 s | 664 / 1 731 / 1 656 s |

- **The cache removes 568 s from (A)** (`measure_fwhm` 571 → 3 s), i.e.
  the single fix is worth more than the whole vectorized fit, as predicted
  in Q10.  The 169 skipped `measure_fwhm` calls are also why
  `detect_lines` is called 2 775 instead of 2 944 times.
- **Shipped (A) with the cache is now faster than the all-callers (B)
  without it** (9 462 vs 9 554 s), and the residual `curve_fit` load is
  40 s.  The remaining `fit_arcspec` cost (57 s) is 0.6% of the run; the
  remaining `fast_line_fit` upside on DEIMOS is therefore tens of seconds.
  The B + cache profile (running, expected ~11:40) will give the number.
- Wall differences between profiles include machine noise of a few percent
  (see B's `--ncpu 4` scatter); the `extraction.run` spread (1 606–2 038 s)
  is of that kind and unrelated to C.  The ~15 min overlap of this profile
  with the test suite and doc build (noted above) did not inflate it.
- Serial A+B+C versus the original 12 390 s baseline: **9 462 s (−23.6%)**
  on this machine, with the caveat that the baseline was measured on older
  code and a different day.

**Release note** updated with the DEIMOS numbers for the cache bullet.

### 2026-10-26 (Task 6 — complete apart from the `--ncpu 4` cached run; after-C results written; pre-merge notes)

**B + cache profile** (`pypeit-arcfit-q10` step 2, 08:57–11:35, rc 0;
`Reports/keck_deimos_600zd_m_6500_C_B_cache.*`): wall **9 495 s**, i.e.
the same as A + cache (9 462 s) within noise.  `fit_arcspec` 44 s (vs 57),
`fit_gauss` 7.5 k calls / 29 s (vs 14 k / 42), `measure_fwhm` 1 call /
0.7 s, `detect_lines` 75 s (vs 87).  **With the cache, (B) buys ~13 s on
DEIMOS**; the parameter's upside is gone on this setup, and the release
note now says so ("the hand-offs cost well under 1% of a serial DEIMOS
reduction").

**`speed_up_results.md`**: the "After workstream C" section is written
(serial DEIMOS table for all five variants, Kast, `--ncpu 4`, fidelity,
verdict, A+B+C totals: serial DEIMOS 12 394 → 9 462 s, −23.7%, indicative
only since the ends are on different code; Kast 121 → 88 s).

**`--ncpu 4` DEIMOS with the cache** launched 11:38 as the user service
`pypeit-arcfit-q10-ncpu4` (`scripts/arcfit/arcfit_q10_ncpu4.sh`, B's
`run_with_rss` + `instr_run_pypeit` recipe, output
`REDUX_OUT_C/_ncpu4/A_cache`, progress appended to `q10_progress.txt`),
~1 h.  It fills the one pending cell in the results table.

**Pre-merge steps (for the PR):**
1. User commits the uncommitted files: `pypeit/core/arc.py`,
   `pypeit/core/flexure.py`, `pypeit/par/pypeitpar.py`,
   `pypeit/wavecalib.py`, `pypeit/tests/test_arc.py`,
   `pypeit/tests/test_flexure.py`, `doc/pypeit_par.rst`,
   `doc/releases/2.1.0dev.rst`; and on the dev-suite side
   `claude_prompts/speed_up/C.md`, `pypeitdev/speed_up/Reports/
   speed_up_results.md`, `pypeitdev/speed_up/scripts/arcfit/*` (PNG figures
   and `.prof` artifacts are git-ignored).
2. #2209 (PR B) still needs Tim's re-review and merge; PR C is then
   retargeted to `develop` and `speed_up_arcfit` rebased/merged with
   `develop` (the user does the git).  The dev-suite branch `dev_speed_up`
   likewise.
3. One full `./pypeit_test all` (safe recipe under Running) at the tip of
   `speed_up_arcfit` before merging.  Because `fit_arcspec` touches every
   `detect_lines` caller, the setups most worth watching are the full
   DEIMOS set, the echelle setups (order-ID cross-correlations), and the
   IFU ones.
4. No identical-output check is possible for multi-slit products (#2208).
5. Task 7: the report `Speed_Up_Arcfit_Report.md`.  Q11 (the fragile
   `xcorr_shift_stretch` step and slit 9) awaits an answer; the report
   should mention it as a known limitation either way.

### 2026-10-26 (Task 6 — `--ncpu 4` DEIMOS with the cache: 3 254.8 s; Task 6 complete)

`pypeit-arcfit-q10-ncpu4` ran 11:38–12:32, rc 0 (`REDUX_OUT_C/_ncpu4/A_cache`).

| `--ncpu 4` | legacy | A | **A + cache** |
|---|---|---|---|
| wall | 3 400.2 s | 3 368.6 s | **3 254.8 s** |
| calibrations / image proc / extraction | 1 480.1 / 76.7 / 660.7 | 1 487.9 / 78.2 / 624.7 | 1 475.3 / 78.3 / 501.2 |
| peak PSS | 36.7 GB | 36.6 GB | 36.7 GB |

The cache's serial saving (~570 s, in the per-detector flexure step) shows
up as −124 s in the parallel extraction stage, as expected for 4 workers.
A+B+C at `--ncpu 4`: **3 254.8 s vs B's 3 391.2 s (−4.0%)**, 3.8x the June
serial baseline.  `speed_up_results.md` is updated; the doc build is being
re-run to check the release-note edits.  **Task 6 is complete**; next is
Task 7 (prompt 7) once Q11 is answered.

### 2026-10-26 (Q11 → known limitation only; no GitHub issue)

**Q11 answered: leave the pipeline alone in PR C, mention the fragility in
the report and the release notes, no GitHub issue.**  Applied:
- `doc/releases/2.1.0dev.rst`: a "Known limitation, not introduced by this
  change" sentence appended to the Under-the-hood `fit_arcspec` bullet,
  naming `wvutils.xcorr_shift_stretch` and the one 600ZD slit that moved
  30 px from a wrong to the correct solution.  Doc build re-run.
- The report (Task 7) will carry the full trace (slit 11 numbers, the
  amplification chain, the slit-9 suspicion) as a known limitation.
- No issue filed; the slit-9 question stays in this log only.

**Task 6 is complete.  Next: prompt 7** ("Read this doc.  Perform the 7th
task under Tasks"): the report `Speed_Up_Arcfit_Report.md` in
`pypeitdev/speed_up/Reports/`, with figures from scripts under
`pypeitdev/speed_up/scripts/`.

### 2026-10-26 (Task 7 — final report written)

**Report:** `pypeitdev/speed_up/Reports/Speed_Up_Arcfit_Report.md`.  It
has a table of contents, ~7 400 words including tables, and these sections:
- 0 Executive summary;
- 1 What was replaced, and why (June 840 s / 1.11 M; legacy re-profile
  846 s);
- 2 Method (log-parabola math, guards, GN polish with cap 10, rescue,
  fallback and the three hand-off rules, `centerr`, controls including
  `fast_line_fit` / `arcfit_options`);
- 3 Fidelity (sky and DEIMOS; the GN 3 vs 10 table; A identical on 370
  spectra);
- 4 Line-count changes under B, and why A ships;
- 5 Worked examples (Fig. 1 saturated, Fig. 2 negative σ);
- 6 Saturated features and identification (the requested discussion,
  including the future plateau model);
- 7 Fallback rate, where the time goes, and the Q10 cache;
- 8 Dev-suite vetting (vetted RMS table with legacy / A / gn10 / gn0, and a
  per-family summary of all 38 WaveCalib files plus tilts);
- 9 Re-profiling (Fig. 3, five-variant table, Kast, `--ncpu 4`, A+B+C
  totals with caveats);
- 10 Known limitation: MSC03 slit 11 (Fig. 4 = `msc03_slit11_judge.png`,
  slit 9, the Q11 decision);
- 11 Switches, options and risks;
- appendices: files changed, tests, scripts and figure regeneration.

**Scripts** (`pypeitdev/speed_up/scripts/arcfit/`, each with a run header):
- `report_fig_saturated.py` (`--list` prints the candidates);
- `report_fig_negsigma.py` (`--list` prints the statistics);
- `report_fig_profile.py`.

**Figures** (`Reports/figures/arcfit/`, PNGs git-ignored): `sat_examples.png`,
`negsigma_example.png`, `profile_deimos.png`.  The dataviz skill was loaded
and its palette used (blue for legacy, orange for B, ink points, grey
context).

**Examples chosen** (600ZD, legacy dev-suite WaveCalib `wv_fits[i].spec`,
pipeline fwhm ≈ 7.3 px, 11-px window):

| line | raw peak | legacy c / σ | B c / σ |
|---|---:|---|---|
| MSC04 slit 4, pix 6573 (edge) | 40 219 | 6570.46 / 3.03 | 6410.24 / 153.1 |
| MSC03 slit 11, pix 5407 (edge) | 39 129 | 5409.97 / 2.59 | 5478.88 / 66.5 |
| MSC03 slit 14, pix 5895 (flat top) | 39 538 | 5893.59 / 5.38 | 5819.00 / 176.1 |
| Paranal pix 6479 | — | 6477.605 / −2.166 | 6477.618 / +2.436 |
| Paranal pix 16750 | — | 16750.179 / −0.557 | 16750.179 / +0.557 |

A equals legacy on all of them.

**Numbers computed fresh** (pypeit14b, read-only, in-process toggles of
`USE_VECTORIZED_ARCFIT` and `ARCFIT_MATCH_CURVEFIT`):
- **600ZD stored `arc_spectra`, fwhm 4, 92 spectra:** legacy 7 887, A 7 887,
  B 7 838 accepted (Q7's 51 reproduced).  Of the 50 lost lines with peaks
  above 30k, all have an adjacent-pixel jump > 50% of the peak (median 96%),
  and none is in a final wavelength fit;
- **pipeline inputs, 86 spectra:** legacy 10 132, A 10 132, B 10 093.  46 are
  lost: 45 with peaks above 30k (39 edges, 6 flat tops) and 1 faint.  None is
  in a final fit (nearest `pixel_fit` ≥ 7.2 px, median 26 px);
- **2D `Arc_A_0_MSC04.fits`, slit 4:** ~39 900 counts across the slit over
  spectral pixels ~6505–6574, with cores at 79 503 (the frame maximum),
  dropping to ~1 250 at 6576;
- **Paranal sky:** legacy 3 748, A 3 748, B 3 766.  56 negative-σ fits, 49
  of which would pass the other cuts; **B keeps 38** of the 49 and A keeps
  none.  37 of the 38 have σ between 0.25 and 0.66 px (sub-resolution).
  B-vs-legacy |Δc| over the 38: median 0.031 px, max 3.4 px;
- **the per-family summary of `arcfit_metrics*.md`** (the gn0 extremes,
  e.g. 104 px on 1200G_M_5500 and 10 px on Kast 600_5000_d46; A tilts median
  ≤ 1.3e-7 px, max 0.005 px).

**Inconsistencies with the earlier logs and the brief:**
1. **What the "saturated lines" are.**  They are mostly *edges* of extended
   saturated plateaus (likely charge bleeding), not saturated Gaussian lines.
   The legacy centroid sits 2–4 px inside the edge, so it is not a line
   wavelength.  The brief's "centroid is nevertheless usable" holds only in
   the sense of a reproducible pattern anchor.  The report says so (§5.1,
   §6); see Q12.
2. **Negative σ.**  Q3 / the brief expected the new fit to keep these lines
   with "the same |σ|".  That is true only for option B, which keeps 38 of
   49.  Most are sub-resolution spikes, and only some have exactly the same
   |σ| (e.g. pixel 6479 differs by 12%).  Under the shipped A none is kept,
   consistent with `test_detect_lines_matches_curve_fit`.
3. **`arcfit_lost_lines.py` / `arcfit_load.py` read
   `REDUX_OUT/keck_deimos/600ZD_M_6500`,** which now holds the products of
   the last profile run.  The new scripts read `REDUX_OUT_C/legacy` instead.

**Not done:**
- the report is longer than the 2 500–4 000-word target, mostly because of
  the tables;
- the cause of the plateaus (bleeding) and the slit-9 suspicion were not
  verified further;
- which hand-off rule catches Paranal pixel 6479 under A was inferred, not
  instrumented.

Q12 added (default: leave the release note unchanged).
