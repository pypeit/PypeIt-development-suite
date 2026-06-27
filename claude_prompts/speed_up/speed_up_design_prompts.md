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
2. Read this doc.  Perform the 3rd task under Profile

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

## Design document

Guidelines for the design document which will be named pypeit_dashboard_design.md and will be stored in PypeIt-development-suite/pypeitdev/dashboard/.  Keep in mind:

- When performing any of the design, be sure to refer to the pypeit_workflow.md document
- You are encouraged to suggest your own design ideas 
- This document will be used to guide the development of the dashboard
- It will not include specific code recommendations, but do keep in mind we will use Python and PyQt6 

### Prep

1. Start the design document by including a preamble of what it is for.  Title that section "Preamble".

   - Add any other information you think is relevant
   - Add a version number to the document (0.1)
   - Add a date to the document (today's date)
   - Add a author to the document (JXP and Claude)

### Initialization 

The following will set the requirements and design for the Initialization phase, i.e. how to launch the Dashboard and what to display when the Dashboard is launched.

Here are a few basic requirements to get us started:

- The Dashboard will launch by running a python script
- The user will launch it from a folder that includes a PypeIt file
- The default view will show the "state" of the PypeIt reduction as a formatted table with color-coded information
- If a state file is present in the folder, load that file to get the state
- Otherwise, get the state as done in the pypeit_status.py script

1. Generate a draft of the initialization design

   - Make its first section titled "Initialization"
   - Include a subsection titled "Requirements"
   - Include sub-subsections to "Requirements" titled "Pending" "Implemented", "Deferred"
   - Examine the state.py module in PypeIt/pypeit/state.py to understand the structure of the state
   - Include the basic requirements above
   - Refine them with proposed design ideas on color, formatting, etc.
   - Emphasize ease of readability and clarity
   - Include a list of your own based on your understanding of the PypeIt workflow
   - If you have any questions, ask them in the Clarifications section below
   - Update the version number to the document (0.1)
   - Log your work in the Logs section below

2. I have answered the questions in the Clarifications section below.  Please:

   - Note that I have hand edited the design document (slightly).  In particular, I removed a few of the Requirements that you proposed. Reread it
   - Update the document to reflect the answers.  
   - Check for consistency throughout the document
   - Update the version number 
   - Log your work in the Logs section below



## Q&A

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