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
- **Branch stack:** A branched from `speed_up` as `speed_up_qa`; B from
  `speed_up_qa` as `speed_up_detpar`; **C branches from `speed_up_detpar` as
  `speed_up_arcfit`**. **PR C targets `speed_up_detpar`** (retarget as the lower
  PRs merge). `speed_up` itself eventually merges to `develop`.
- PR B must be complete on `speed_up_detpar` before starting.

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
  (despite the name); it is consumed only as `all_ecent` in
  `wvutils.arc_lines_from_spec` → `HolyGrail` pattern weights.
- `fit_gauss` has exactly one other, cold caller
  (`pypeit/multislit_flexure.py`) — leave it untouched.
- Downstream filters in `detect_lines` are unchanged and still apply:
  `twid > 0`, `twid < fwhm_max/2.35`, `tcent > 0`,
  `tampl_true < nonlinear_counts`, `abs(tcent-pixt) < 0.75*fwhm`.

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

If you need to run Python or PypeIt, use the **`pypeit` conda environment**.
The shell's default environment shadows it, so use the absolute env binaries:

```bash
/home/xavier/miniconda3/envs/pypeit/bin/python
/home/xavier/miniconda3/envs/pypeit/bin/run_pypeit
/home/xavier/miniconda3/envs/pypeit/bin/pytest
```

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
    - Confirm PR B is complete on `speed_up_detpar`. Then create and check out
      **`speed_up_arcfit`** from **`speed_up_detpar`**; report the HEAD commit.
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
    - Run `pytest pypeit/tests` in the `pypeit` env; must be green.
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
    - Finally, note the remaining pre-merge steps for the stack: the manual
      DEIMOS `--ncpu 4` identical-output run (per Q6) and one full dev-suite pass
      (`./pypeit_test all -t 2`) at the tip of `speed_up` before merging to
      `develop`. Ask any remaining questions in Q&A.

## Q&A

Claude poses questions here (as `**Qn — title.** body` followed by a `>A:` line);
the user answers inline beneath each `>A:`.

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs
