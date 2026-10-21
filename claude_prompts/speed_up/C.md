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

## Q&A

Claude poses questions here (as `**Qn — title.** body` followed by a `>A:` line);
the user answers inline beneath each `>A:`.

**Q1 — Stray `speed_up_C` branch.**  `speed_up_arcfit` already existed when
Task 1 started, at the `speed_up_B` tip `72e448f75`, so I only checked it out.
A second local branch, `speed_up_C`, points at the same commit and has no
upstream.  My default is to **leave it alone**.  Tell me if you want it
deleted (`git branch -d speed_up_C`).
>A:

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
>A:

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
>A:

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
>A:

**Q5 — `centerr` from the Gauss–Newton Jacobian rather than the log-parabola?**
`centerr` is never used numerically (see Key facts), so this is about fidelity,
not results.  Once Gauss–Newton has run, the closest match to legacy is
`curve_fit`'s own definition: `pcov = s²·inv(JᵀJ)` with linear-space residuals,
`s² = Σr²/(n−3)` over the in-bounds pixels, at the converged parameters.  It is
nearly free, because `JᵀJ` is already built in the last iteration.  My default:
- `GN_ITER > 0`: take `centerr` from the Jacobian;
- `GN_ITER = 0`: use §C.2's `_center_variance` propagation.

Both are implemented and tested against legacy `fitcov[1,1]`.
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
