# Speed up PypeIt -- Design doc

## Goals

We wish to generate a design document that will be used to guide a refactor of PypeIt to speed it up.

## Claude

### Skills

Consider using the skills in PypeIt-development-suite/.claude/skills/

## Context

Here is information that will be useful to you:

   - Read the docs in PypeIt/doc/
   - Identify the key components of the PypeIt workflow
   - Read the the `pypeit_workflow.md` document in `PypeIt-development-suite/pypeitdev/dashboard`
   - Log your work in the Logs section below

## Running

If you need to run Python, use the `pypeit` environment

## Prompts

### Profile

1. Read this doc.  Perform the first task under Profile
2. Read this doc.  Perform the 2nd task under Profile
3. Read this doc.  Perform the 3rd task under Profile
4. Read this doc.  Perform the 4th task under Profile
5. Read this doc.  Perform the 5th task under Profile


## Profile

1. Now let's profile the current PypeIt on a test dataset.  Use the dataset in PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55 and run the reduction in PypeIt-development-suite/REDUX_OUT/shane_kast_blue/600_4310_d55/shane_kast_blue_A where the PypeIt file is shane_kast_blue_A.pypeit.  Do:

    - Be aware this is one of the simplest spectrographs in PypeIt
    - Prepare to profile this command 
    - Ask questions in the Q&A section below
    - Log your work in the Logs section below

2. I have responded to your questions.  Read them and then:
    - Prepare to profile this command 
    - Ask any additional questions in the Q&A section below
    - Log your work in the Logs section below

3. Read my additional answers in the Q&A section.  Then, do it:
    - Generate the code necessary to do your profiling in `PypeIt-development-suite/pypeitdev/speed_up/scripts`
    - Run it
    - Write a report in `PypeIt-development-suite/pypeitdev/speed_up/Reports`
    - Log your work

4. Ok, that report is excellent.  Now, prepare to profile the data in `PypeIt-development-suite/REDUX_OUT/keck_deimos/600ZD_M_6500`.  If you have any questions, put them in the Q&A/keck_deimos section. Log your work.

5. I have answered your questions.  Please read them and proceed to profile `keck_deimos`.  Write a new report in `PypeIt-development-suite/pypeitdev/speed_up/Reports`.

## Plan

1. Given what you have learned so far, develop a plan to speed up PypeIt under these constraints:
    - KISS: Keep it simple stupid!  
        - Avoid any significant refactoring
        - Make use of multi-CPUs
        - Identify places to use more `numpy` array manipultions/math
    - The b-spline code is being worked on in another branch, so ignore it

Then start a new doc named `PypeIt-development-suite/pypeitdev/speed_up/Reports/speed_up_design.md` for us to iterate on.   Ask question in Q&A/Planing below. Log your work


## Q&A

### Planning

Questions from Claude on the speed-up plan (`Reports/speed_up_design.md`, v0.1).
Please answer inline beneath each "A:". B-spline is excluded throughout (other
branch).

**Q1 — Parallelism mechanism, default, and new parameter.** I propose
parallelizing the existing per-detector loops in `pypeit/exposure.py` with a
`concurrent.futures.ProcessPoolExecutor` (stdlib, no new dep), gated by a new
parameter `par['rdx']['ncpu']` (name negotiable: `ncpu` vs `nproc`), **default
`1`** so current behavior is byte-for-byte preserved and >1 is strictly opt-in.
OK with stdlib multiprocessing + an opt-in `ncpu` (default 1)? Or do you want a
smarter default (e.g. `min(n_detectors, os.cpu_count()-1)`) given the heavy
per-detector RAM?
A:

**Q2 — Granularity to target first.** I recommend **detector/mosaic-level**
parallelism first (biggest win for DEIMOS/LRIS/MOSFIRE; finest memory
granularity; the loops are already isolated). Coarser levels (calibration
groups, frame/comb-id loops) are Tier 3 / later. Agree to start at the detector
level only?
A:

**Q3 — Numerical-identicality tolerance.** Vectorizing the arc-line Gaussian fit
(replacing per-line `scipy.curve_fit` with an analytic weighted log-parabola fit)
will change wavelength/tilt results at the **noise level** (not bitwise
identical). Is "within tolerance, vetted by the dev-suite wavelength/tilt RMS
checks" acceptable, or must outputs remain numerically identical (which would
rule this change out)?
A:

**Q4 — Determinism under multiprocessing.** The dev suite requires deterministic
output. The parallel path recombines results by detector key (order-independent)
and seeds are per-call, so I expect `ncpu>1` to match `ncpu=1`. Do you want a
hard requirement that parallel and serial outputs be identical (I'll add a test),
or is "deterministic for a fixed `ncpu`" sufficient?
A:

**Q5 — Validation scope.** I plan to validate on the two profiled setups
(shane_kast_blue for the single-detector/vectorization path, keck_deimos for the
multi-detector/parallel path) plus the relevant `vet_tests`. Want any additional
instrument/setup in the first validation pass (e.g. a keck_lris or mosfire setup)?
A:

**Q6 — Early cheap-wins PR.** The QA cheap wins (force matplotlib `Agg` for
headless QA; option to skip/defer/parallelize PNG writes) are independent and
low-risk. Want these as a **separate early PR** ahead of the parallelism/vectorize
work, or bundled into one effort?
A:

### shane_kast_blue

Questions from Claude for the Profile task (please answer inline beneath each
"A:"). Profiling (task 2) will proceed once these are settled.

**Q1 — `speed_up` branch import is broken; OK with my workaround?**
On this branch `import pypeit` fails: `pypeit/__init__.py:8` does
`from .pkg.version import version`, but the setuptools-scm generated version file
still sits at the *old* path `pypeit/version.py` while `pyproject.toml:189` now
points to `pypeit/pkg/version.py` (the "rename pkg dir" commit moved the import
target but the env was never reinstalled, so the file at the new path was never
written). `run_pypeit` cannot start at all in this state. To unblock, I copied
`pypeit/version.py` → `pypeit/pkg/version.py` (both are git-ignored, so no commit
impact). Do you want a proper fix instead — `pip install -e .` to regenerate, or
add `setuptools_scm` to the `pypeit` env — or is the copy fine for the profiling
work?
A: Can you try agani?  If you need to run PypeIt, use the `pypeit` environment.

**Q2 — Profiling method?** I propose **cProfile** as the primary tool
(deterministic, function-level call graph; saved as a `.prof` plus a pstats
top-by-cumulative-time table), complemented by a **per-step wall-clock timeline**
parsed from the run `.log` (PypeIt timestamps each calibration/reduction step).
Optionally I can add **py-spy** sampling for a flame graph if you want one (needs
install). OK with cProfile + log-timeline?
A: Use **cProfile**

**Q3 — Cold or warm run?** The redux dir already contains built `Calibrations/`
(and an `Intermediate/`). For a speed-up baseline I plan a **clean cold run**
(remove the existing `Calibrations/`, `Science/`, `QA/`, `Intermediate/` first,
then run with `-o`) so calibration-building time is included in the profile.
Confirm cold, or do you also want a warm (calibs reused, science-only) timing for
comparison?
A: Yes, a **clean cold run**

**Q4 — Repeats?** Single timed run, or average over N (e.g. 3) for stability? On
this small dataset (~1m runtime) a single run is usually enough; I'll note
run-to-run variance caveats either way.
A: Single run.

**Q5 — Report granularity?** I plan to map the cProfile hotspots onto the
`pypeit_workflow.md` phases (initialization → calibrations [bias/arc/tiltimg/
slits/wave/tilts/flat] → standard-star reduction → science reduction → QA), so
the report answers "where does the wall-clock go, by phase and by function." Is
that the breakdown you want?
A:  That sounds good.

### Additional questions (task 2)

No blocking questions — the five answers above fully settle the plan. Two minor
decisions I'll take by default unless you say otherwise:

**Q6 — Cold-run cleanup scope.** To make the run cold I'll remove the output
products `Calibrations/`, `Science/`, `QA/`, `Intermediate/` and the stale
`shane_kast_blue_A_state.json`. I'll **move aside** (not delete) the prior
`shane_kast_blue_A.log` and the old `_UTC_*.par` files into a `pre_profile_bak/`
subfolder so nothing is lost, then run with `-o`. The `.pypeit` input and the raw
data are untouched. OK?
A: You can remove the log and .par files too.

**Q7 — cProfile completeness.** Confirmed the shane_kast_blue reduction is
**single-process**: there is no `multiprocessing`/`concurrent.futures`/`joblib`
use in `pypeit/` core and no `nproc` parameter, so cProfile captures the entire
workload in one process (no child-process blind spots). Noting this as a fact,
not a question — no action needed from you.
A: ok

### keck_deimos

Questions from Claude for the **keck_deimos 600ZD_M_6500** Profile task (4th task
under Profile). This is a prepare-only step; the profiling run will proceed once
these are settled. Please answer inline beneath each "A:".

This dataset is a much heavier target than shane_kast_blue: it is **MultiSlit
with a slitmask (many slits per detector)** and keck_deimos defaults to
**4 detector mosaics** `[(1,5),(2,6),(3,7),(4,8)]` (the 8 chips paired), with
`detnum` left commented out in the `.pypeit` so the full 4-mosaic reduction runs
by default. Expect the run (and the cProfile overhead on top) to be substantially
longer than Kast — plausibly tens of minutes.

**Q1 — Detector/mosaic scope?** Profile (a) the **full 4-mosaic cold run** (the
default; a true real-world baseline, but long under cProfile), or (b) **restrict
to a single mosaic**, e.g. `detnum = (3,7)` as in the commented QL line, for a
faster, still-representative profile? I recommend **(a) full** for an honest
baseline and will report a per-mosaic breakdown, but I'll flag the runtime. Your
call.

A: Yes, the full run.

**Q2 — Run location & cold scope.** Unlike Kast, the `.pypeit` sits *directly* in
`REDUX_OUT/keck_deimos/600ZD_M_6500/` (no `_A` subdir). I'll run in place there;
PypeIt will create `Calibrations/`, `Science/`, `QA/` in that dir. The dir is
already essentially clean (only the `.pypeit` and a stale `.test.log`), so the
cold-clean is nearly a no-op — I'll remove any products/`_state.json`/`.calib`/
old `.log`/`.par` if a prior partial run exists, then run with `-o`. The raw data
and `.pypeit` are untouched. OK?

A: ok

**Q3 — Reuse the Kast methodology?** Same approach as the Kast report: in-process
`cProfile` driver, pstats cumulative + tottime tables, a per-step wall-clock
timeline parsed from the timestamped run log, and a report mapping hotspots onto
the `pypeit_workflow.md` phases — **extended with a per-mosaic breakdown** so we
can see how cost scales with detector count and slit count. New scripts go in
`pypeitdev/speed_up/scripts/` (e.g. `profile_deimos.py`, reusing
`analyze_profile.py`) and the report in `Reports/`. Confirm.

A: Yes, same methodology.

**Q4 — Single science frame, single cold run?** Use the `.pypeit` as-is: one
science frame (`d1010_0056`; `0057` is commented out), one cold run, no repeats —
same as Kast. Confirm.

A: confirmed

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-06-26 (Prepared to profile shane_kast_blue_A; unblocked broken env)

Performed the **1st task under Profile** (prepare + ask questions; the actual
profiling run is task 2). Re-grounded myself in the workflow via
`pypeitdev/dashboard/pypeit_workflow.md` (esp. §6, the same Kast 600/4310 dataset)
and inspected the target reduction.

Environment / dataset reconnaissance:
- The `pypeit` conda env (`/home/xavier/miniconda3/envs/pypeit`) is the editable
  install of *this* working repo, currently on branch `speed_up` (HEAD
  `e9ed85c1a`). The shell's default env is `base`, which has no `run_pypeit`, so
  I use the absolute env binary `…/envs/pypeit/bin/{python,run_pypeit}` (matches
  the conda-env-shadowing memory).
- **Found and fixed a blocker:** `import pypeit` failed with
  `ModuleNotFoundError: No module named 'pypeit.pkg.version'`. Cause:
  `pypeit/__init__.py:8` imports `from .pkg.version import version`, and
  `pyproject.toml:189` sets `version_file = "pypeit/pkg/version.py"`, but the
  setuptools-scm generated file still lives at the *old* path `pypeit/version.py`
  — the "rename pkg dir" commit moved the import target without a reinstall to
  rewrite the generated file. `setuptools_scm` isn't installed in the runtime
  env, so I unblocked by copying `pypeit/version.py` → `pypeit/pkg/version.py`
  (both git-ignored). `run_pypeit -h` and `import pypeit` (v1.17.5.dev1018+
  g42ed28072) now work. Raised this as Q1 — a proper reinstall may be the right
  fix.
- Confirmed dataset: `RAW_DATA/shane_kast_blue/600_4310_d55` (25 frames:
  1 arc/tilt, 10 bias, 11 dome flats, 2 science, 1 standard), one config/calib
  group/detector/slit — matches workflow §6. The redux dir already holds built
  `Calibrations/`, `Science/`, `QA/`, `Intermediate/`, a prior `.log`, and a
  `_state.json` (so a cold run needs them cleared first — Q3).
- Verified the `Reports` output dir exists at `pypeitdev/speed_up/Reports`.

Profiling plan staged for task 2 (pending Q&A answers): run
`run_pypeit shane_kast_blue_A.pypeit -o` from the redux dir under cProfile
(`python -m cProfile -o <out>.prof <run_pypeit> …`), produce a pstats
cumulative-time table + a per-step wall-clock timeline parsed from the run
`.log`, map hotspots onto the `pypeit_workflow.md` phases, examine the QA PNGs,
and write the report to `pypeitdev/speed_up/Reports`.

Posed five questions in the Q&A section (env-fix workaround, profiling method,
cold-vs-warm run, repeats, report granularity).

### 2026-06-27 (Read Q&A answers; staged profiling driver; no blockers remain)

Performed the **2nd task under Profile** (read answers, prepare to profile, ask
any additional questions, log). No actual profiling run yet — that is task 3.

Answers incorporated:
- **Q1 (broken import):** Re-tried in the `pypeit` env and `import pypeit`
  succeeds — `1.17.5.dev1018+g42ed28072`. The generated `pypeit/pkg/version.py`
  is present (my earlier git-ignored copy), so the env is usable as-is; a
  `pip install -e .` would be the clean long-term fix but is not needed to
  profile. `which run_pypeit` is empty in the default shell (base env), so I'll
  keep using the absolute env binaries.
- **Q2:** cProfile (primary), plus a per-step wall-clock timeline from the run
  `.log`.
- **Q3:** Clean cold run with `-o`.
- **Q4:** Single run.
- **Q5:** Report mapped onto `pypeit_workflow.md` phases.

Preparation done:
- Re-confirmed the dataset/redux dir: 1 arc/tilt, 10 bias, 11 dome flats, 2
  science, 1 standard; the redux dir holds prior `Calibrations/`, `Science/`,
  `QA/`, `Intermediate/`, a `.log`, a `_state.json`, and several `_UTC_*.par`
  files (cold run must clear these — see Q6).
- Verified the `run_pypeit` console wrapper calls
  `RunPypeIt.entry_point()` → `RunPypeIt.main(RunPypeIt.parse_args())`, and that
  `parse_args(options)` accepts an explicit args list and resolves cwd-defaults
  via `Path.cwd()`. This lets me drive the run under `cProfile` in-process.
- **Confirmed single-process** (Q7): no `multiprocessing`/`concurrent.futures`/
  `joblib`/fork use in `pypeit/` core and no `nproc` param, so cProfile sees the
  whole workload — no child-process blind spots.
- **Staged the profiling driver** (not yet executed) at
  `scratchpad/profile_kast_blue.py`: it `chdir`s to the redux dir, runs
  `RunPypeIt.main(parse_args([pypeit_file, "-o"]))` wrapped in `cProfile.Profile`,
  dumps `shane_kast_blue_A.prof` into `pypeitdev/speed_up/Reports`, and prints a
  top-40 cumulative-time table. This will be finalized/run in task 3 (along with
  the log-timeline parse, phase mapping, QA inspection, and the written report).

Raised two minor follow-ups in Q&A: Q6 (cold-cleanup scope — propose moving the
old `.log`/`.par` aside to `pre_profile_bak/` rather than deleting) and Q7 (a
confirmation, not a question). Neither blocks task 3.

### 2026-06-27 (Task 3 — profiled the cold run; wrote report; B-spline is the headline)

Performed the **3rd task under Profile**: per Q6 (delete the old log/.par too)
generated the profiling code, ran it cold, examined QA, and wrote the report.

Code (in `pypeitdev/speed_up/scripts/`):
- `profile_kast_blue.py` — wipes the redux dir cold (`Calibrations/`, `Science/`,
  `QA/`, `Intermediate/`, `*_state.json`, `*.calib`, prior `.log`, old
  `*_UTC_*.par`; raw + `.pypeit` untouched), then runs
  `RunPypeIt.main(parse_args([file, "-o"]))` in-process under `cProfile`, dumping
  `Reports/shane_kast_blue_A.prof` + a `runmeta.json` (wall time) + a copy of the
  timestamped run log.
- `analyze_profile.py` — pstats tables (top-40 cumulative + tottime →
  `.pstats.txt`) and a per-phase wall-clock timeline parsed from the timestamped
  log → `.timeline.txt`.

Run: cold, single, return code 0. **PypeIt internal exec time 1m 58.72 s
(≈119 s)**; cProfile wall 121.2 s; 50.9 M function calls. Healthy outputs
(3 spec1d + 3 spec2d, 20 QA PNGs across arc/tilt/wave/flat/trace/profile/flexure).

Key findings (full report:
`pypeitdev/speed_up/Reports/shane_kast_blue_A_profile_report.md`):
- **Headline: the pure-Python B-spline kernel (`pypeit/core/bspline/`) is ≈39 s
  of self-time, ~32 % of the whole run** — `intrv` 10.7 s, `cholesky_band` 8.7 s,
  `bsplvn` 6.8 s, `solution_arrays` 4.6 s, `cholesky_solve` 3.9 s,
  `bspline_model` 2.3 s. It is the shared engine of global sky, local sky, and
  extraction, so optimizing it (Cython/Numba/vectorize) helps *every* path.
- **Important nuance:** this pure-Python kernel is the status quo on `release`
  too (both import `bspline/util.py`; neither wires in a compiled path) — *not* a
  refactor regression. A stale leftover `pypeit/bspline/_bspline.*.so` (no longer
  used) shows these kernels were historically Cython-compiled → re-compiling is
  feasible/precedented.
- Phase split (cumtime, mapped to `pypeit_workflow.md` §3): extraction /
  local-sky 40.5 s (34 %), global sky 17.2 s (14 %), build calib images 13.6 s
  (11 %), slit edges 8.1 s, per-frame image proc 7.2 s (lacosmic+CR ≈ 6 s),
  object trace/`moment1d` 8.2 s, tilts 4.4 s, wave 2.2 s, flat 1.4 s,
  QA/plotting ≈ 12 s (10 %), init ≈ 1 s.
- **QA plotting ≈ 12 s** and the default matplotlib backend is the interactive
  **`qtagg`** even though QA is headless PNGs → forcing `Agg` is an easy win.
- Ranked opportunities in the report: (1) B-spline kernel, (2) cut redundant
  B-spline fits, (3) headless/optional/parallel QA, (4) CR rejection, (5) NumPy
  allocation churn, (6) detector/frame parallelism (structural; helps big
  instruments more than single-detector Kast). Items 1–3 ≈ 45 % of runtime.

Note: the log-timeline labels are approximate (the "sky subtraction" aggregate
folds in the local-sky+extract block); the cProfile tables are authoritative.

### 2026-06-28 (Task 4 — prepared to profile keck_deimos 600ZD_M_6500; posed questions)

Performed the **4th task under Profile** (prepare + ask questions; the actual
keck_deimos profiling run is a later step). Reconnaissance of the new, much
heavier target.

Dataset / redux recon (`REDUX_OUT/keck_deimos/600ZD_M_6500/`):
- The `.pypeit` (`keck_deimos_600zd_m_6500.pypeit`) sits **directly** in the
  redux dir (no `_A` subdir, unlike Kast). The dir is otherwise clean — only the
  `.pypeit` and a stale `.test.log` (records a Ctrl+C from a prior aborted run) —
  so the cold-clean here is essentially a no-op.
- Frames (raw under `RAW_DATA/keck_deimos/600ZD_M_6500/`, gzipped, large —
  42–124 MB each): 2 bias, 2 arc/tilt, 1 illum/pixel/trace flat, 1 science
  (`d1010_0056`, 1800 s; `0057` commented out). Several blue calib frames are
  commented out (pipeline can't handle that lamp combo yet).
- `.pypeit` overrides: wavelength `lamps = NeI, ArI, CdI, KrI, XeI, HgI, ZnI`;
  `detnum` is **commented out** so the keck_deimos default applies.

Why this is a different beast from Kast (from `pypeit/spectrographs/keck_deimos.py`):
- keck_deimos has **`ndet = 8`** and defaults to **4 mosaics**
  `par['rdx']['detnum'] = [(1,5),(2,6),(3,7),(4,8)]` (`allowed_mosaics`,
  `default_mosaic`). So the default run reduces 4 mosaicked detectors, each
  pairing 2 chips.
- It is **MultiSlit with a slitmask** → many slits per detector, vs Kast's single
  slit. The B-spline-heavy sky/extract work (the Kast headline, ~32 % of runtime)
  scales with slits × detectors, so I expect it to dominate even more here.
- Still single-process (the Q7 finding — no multiprocessing/`nproc` in core — is
  global), so cProfile will capture the whole workload; but the wall-clock and
  cProfile overhead will be far larger (plausibly tens of minutes).

Posed four questions in **Q&A/keck_deimos**: (Q1) detector scope — full 4-mosaic
baseline vs a single-mosaic `(3,7)` fast profile; (Q2) run-in-place + cold-clean
scope; (Q3) reuse the Kast methodology, extended with a per-mosaic breakdown; (Q4)
single science frame / single cold run. No code run, no profiling yet — awaiting
answers.

### 2026-06-29 (Task 5 — profiled keck_deimos cold full 4-mosaic run; wrote report)

Performed the **5th task under Profile**: read the keck_deimos Q&A answers (full
run, run-in-place cold, same Kast methodology + per-mosaic breakdown, single
science frame) and profiled `keck_deimos_600zd_m_6500`. (The cold cProfile run
itself had been launched on 2026-06-28 — artifacts dated Jun 28 08:49 — but we
were cut off before the analysis/report; this entry finalizes it.)

Code: added `pypeitdev/speed_up/scripts/profile_deimos.py` (cold clean + in-process
`cProfile` driver, same shape as `profile_kast_blue.py`) and extended
`analyze_profile.py` with `--stem` + a per-mosaic wall-clock breakdown
(`DET_MARKER` on the "Calibrating detector (a, b)" log lines).

Run: cold, single, **rc 0**, healthy (1 spec1d + 1 spec2d, calibrations for all 4
mosaics, 971 QA PNGs; the 8 log "error" hits are benign `Max centroid error: None`
INFO lines). **PypeIt wall-clock 12 393.8 s ≈ 3 h 26 m** (~104× Kast's 119 s);
4.10 billion function calls.

Key findings (full report:
`pypeitdev/speed_up/Reports/keck_deimos_600zd_m_6500_profile_report.md`):
- **Headline confirms Kast at scale:** the pure-Python B-spline kernel
  (`pypeit/core/bspline/`) is again the dominant engine — **≈2 750 s self-time,
  ~22 % of the run** (`cholesky_band` 856 s, `solution_arrays` 654 s,
  `cholesky_solve` 437 s, `intrv` 312 s, `bsplvn` 252 s, `bspline_model` 238 s),
  reached via `bspline_profile` → `local_skysub_extract` (cum 2 952 s) and
  `global_skysub` (cum 1 677 s).
- **Phase split (wall-clock):** sky subtraction 4 284 s (35 %), flat 2 173 s
  (18 %), slit edges 1 896 s (15 %), extraction 1 097 s (9 %), tilts 779 s, wave
  771 s, object finding 728 s, image combine 299 s. Sky+extract ≈ 44 %, calib
  building ≈ 48 %.
- **New at DEIMOS scale (negligible on Kast):** (a) mosaic resampling via
  scipy.ndimage — `geometric_transform` 462 s, `correlate` 577 s,
  `build_image_mosaic` cum 662 s; (b) arc-line fitting — **1.11 M** `curve_fit`
  calls, `fit_gauss` cum 840 s, `gauss_3deg` 332 s over 93 M calls; (c)
  `moment1d` 212 s over 1.05 M calls; (d) flexure `spec_flexure_slit` cum 825 s;
  (e) QA matplotlib text-metrics cum 894 s + 971-PNG encode 276 s.
- **Per-mosaic caveat:** calibration building is even (~1 300–1 600 s/mosaic), but
  the marker-based split lumps the whole post-calib science loop (~6 200 s) into
  the (4,8) interval — the science half is shared, not (4,8)-specific. Noted in
  the report; finer per-mosaic science attribution needs extra instrumentation.
- **Ranked opportunities:** (1) B-spline kernel [helps every instrument],
  (2) mosaic resampling, (3) arc-line `curve_fit` vectorization, (4) headless/
  parallel QA, (5) `moment1d`, (6) flexure, (7) NumPy/masked-array churn,
  (8) detector/mosaic + frame parallelism (structural — 4 independent mosaics,
  the obvious DEIMOS/Echelle/IFU win to design in).

This completes the Profile section (Kast + DEIMOS). Next up: the **Design
document** prompts (Preamble + Initialization).

### 2026-06-29 (Plan task 1 — drafted KISS speed-up plan; new doc + questions)

Performed the **1st task under Plan**: developed a speed-up plan under the KISS
constraints (no significant refactor, use multi-CPU, more numpy vectorization,
**B-spline excluded** — other branch), started the iteration doc, and posed
questions.

Grounding (beyond the two profiles): re-read `pypeit_workflow.md` §3 (reduction
loop structure) and inspected the actual code paths to confirm the plan is
KISS-feasible:
- **Reduction is already factored for parallelism.** On this branch `pypeit.py`'s
  `reduce_all` → `reduce_calibID` → `exposure.reduce_exposure` are module-level
  functions taking explicit args, and the per-detector work lives in
  `exposure.process_exposure` / `objfind_exposure` / `extract_exposure`, each
  looping `for det in detectors` and returning per-detector dicts (calling
  `pypeit_steps.{process_one_det, objfind_one, extract_det}`). The loop bodies are
  independent → a process-pool `map` over detectors is a near-drop-in change. This
  is the single biggest, lowest-risk win (DEIMOS = 4 independent mosaics ⇒ ~4× on
  the per-detector bulk; single-detector Kast gains nothing, as expected).
- **Top vectorization target confirmed:** `fit_arcspec` (`core/arc.py:1138`) fits
  every arc line with `scipy.curve_fit` via `fitting.fit_gauss` in a Python loop →
  the 1.11 M `curve_fit` calls / 93 M `gauss_3deg` evals seen in the DEIMOS
  profile. A 3-param Gaussian has a closed form (weighted log-parabola), batchable
  across all lines with numpy — localized, helps wavecal + tilts.
- **Secondary vectorize candidate:** `moment1d` (`core/moment.py`) — 1.05 M calls
  via `trace.masked_centroid`/`follow_centroid`; numpy-internal already, cost is
  call volume from row-by-row trace following ⇒ more invasive, flagged not
  committed.
- **Cheap wins:** force matplotlib `Agg` for headless QA (run uses interactive
  `qtagg`); parallelize/defer the 971 QA PNG writes; trim obvious numpy churn.

Wrote the draft `pypeitdev/speed_up/Reports/speed_up_design.md` (v0.1): purpose,
constraints, profile-derived cost classes, a ranked plan (Tier 1 = detector-level
multiprocessing + arc-fit vectorization; Tier 2 = cheap QA/numpy wins + moment1d;
Tier 3 = coarser calib-group/frame parallelism), validation strategy (regress
Kast + DEIMOS, vet_tests, re-profile, determinism check), out-of-scope list, and
open questions.

Posed six questions in **Q&A/Planning**: (Q1) multiprocessing mechanism + new
`ncpu`/`nproc` param + default; (Q2) start at detector-level granularity; (Q3)
tolerance for non-identical numerical output from vectorized fits; (Q4)
determinism requirement under multiprocessing; (Q5) validation scope/instruments;
(Q6) whether the cheap QA wins ship as a separate early PR. Awaiting answers
before any code.