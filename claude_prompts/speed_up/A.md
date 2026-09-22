# Speed up PypeIt -- PR A (QA cheap wins + `ncpu` plumbing)

## Goals

Implement **workstream A** of the locked speed-up plan: force the headless
matplotlib `Agg` backend for reductions, introduce the `par['rdx']['ncpu']`
parameter and the matching `run_pypeit --ncpu` command-line flag (the CLI
overrides the parameter), and write the QA PNGs concurrently under that same
`ncpu` control. This is the smallest and lowest-risk of the three stacked PRs,
and it is where the `ncpu` plumbing that PRs B and C build on lands. Target:
recover a chunk of the ~0.9–1.2 ks QA budget on Keck/DEIMOS and ~10 s on
shane_kast_blue.

## Claude

### Skills

Consider using the skills in `PypeIt-development-suite/.claude/skills/` — in
particular `add-parameter` (for `[rdx] ncpu`), `build-docs`, `update-changelog`,
and `run-dev-suite`.

## Context

Authoritative documents (read the coding doc section first; it is the
implementation contract):

- **`PypeIt-development-suite/pypeitdev/speed_up/Reports/speed_up_coding.md` §1**
  — "PR A — QA cheap wins (+ the `ncpu` plumbing)". Sub-anchors §A.1 (parameter),
  §A.2 (CLI flag), §A.3 (`Agg`), §A.4 (deferred QA writes), §A.5 (call-site
  conversions), §A.6 (tests), §A.7 (validation/re-profiling), §A.8 (docs).
  §0.2 of that doc is a grounded code map (file → symbol → line numbers on the
  `speed_up` branch at `e9ed85c1a`).
- **`Reports/speed_up_design.md` v0.5 (locked)** — the design this implements;
  §0 lists all locked decisions, §8 the roadmap.
- **Profile reports** — `Reports/shane_kast_blue_A_profile_report.md` and
  `Reports/keck_deimos_600zd_m_6500_profile_report.md` (the QA numbers this PR
  targets: matplotlib `_get_text_metrics_with_cache_impl` cum **894 s**, PNG
  `ImagingEncoder.encode` **276 s** over **971 PNGs** on DEIMOS; ~12 s over 20
  PNGs on Kast; the default backend is the *interactive* `qtagg`).
- **Profiling scripts** — `pypeitdev/speed_up/scripts/profile_kast_blue.py`,
  `profile_deimos.py`, `analyze_profile.py`.

Repository and branch:

- Code repo: `/mnt/tank/Astronomy/PypeIt/PypeIt` (same as
  `/home/xavier/Projects/PypeIt/PypeIt`).
- **Branch stack:** A branches from `speed_up` as **`speed_up_qa`**; B will
  branch from `speed_up_qa` as `speed_up_detpar`; C from `speed_up_detpar` as
  `speed_up_arcfit`. **PR A targets `speed_up`.**
- The working tree may currently be on some other branch. Create and check out
  `speed_up_qa` from `speed_up` before editing anything.

Locked decisions relevant to PR A (a fresh session needs no other context):

- **`ncpu` default is 1, and `ncpu=1` must remain the literal current serial
  path.** No behavioral change for a default run.
- **CLI overrides the parameter**: `run_pypeit --ncpu N` wins over
  `[rdx] ncpu` in the `.pypeit` file.
- **A thread pool is sufficient for the QA PNGs** (render/encode bound).
  Figures are created and closed **only on the main thread**; only
  `Figure.savefig` is handed to a worker thread.
- **No skip-QA switch.** No new parameter beyond `ncpu`.
- **Q1 (accepted):** forcing `Agg` respects `run_pypeit -s/--show` —
  `matplotlib.use('Agg', force=True)` is applied `if not args.show:`. No new
  parameter is introduced by this guard.
- **Q2 (accepted):** if the QA thread pool turns out to be GIL-bound (matplotlib's
  FreeType text-metrics path may not release the GIL; PIL's PNG encode does),
  **accept the partial win** — keep `Agg` + the deferred-save machinery (free at
  `ncpu=1`) and revisit only if QA is still a top-5 hotspot after PR B. Do **not**
  escalate QA to a process pool in this PR.
- **Q5 (accepted, forward-looking):** in PR B, QA will be written serially inside
  detector worker processes (`parallel._worker_init` calls `qa.init_qa_pool(1)`).
  Write `qa.init_qa_pool` in this PR so it can be safely called a second time to
  *reset* the pool in a forked child.
- **Out of scope:** the B-spline kernel (`pypeit/core/bspline/`) — another branch.
  Do not touch it.
- Do not convert the coadd / sensfunc / telluric / `PdfPages` QA writers; only
  the high-volume per-slit reduction writers listed in coding-doc §A.5.

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
    - In `/mnt/tank/Astronomy/PypeIt/PypeIt`, create and check out
      **`speed_up_qa`** from **`speed_up`**. Confirm the tree is clean first and
      report the HEAD commit you branched from.
    - Read `Reports/speed_up_coding.md` **§0** (scope, branch stack, code map)
      and **§1** in full.
    - Read the code you are about to change: `pypeit/par/pypeitpar.py`
      (`ReduxPar`, ~2758–2900), `pypeit/scripts/run_pypeit.py` (all of it),
      `pypeit/pypeit.py` (`PypeIt.__init__`, `build_qa`, `calib_all`,
      `reduce_all`), and `pypeit/qa.py` (imports, `set_qa_filename`, and the
      `arc_tilts_*_qa` / `spec_flexure_qa` / `spat_flexure_qa` writers).
    - Verify the line numbers quoted in coding-doc §0.2/§1 still hold; note any
      drift.
    - Confirm the environment works: `run_pypeit -h` in the `pypeit` env.
    - Ask questions in the Q&A section below. Log your work in the Logs section.

2. **Add the `ncpu` parameter and the `--ncpu` CLI flag.** Coding-doc **§A.1**
   and **§A.2**.
    - `pypeit/par/pypeitpar.py`: add `ncpu` to the `ReduxPar.__init__` signature,
      add the `defaults`/`dtypes`/`descr` block (default **1**, dtype `int`), add
      `'ncpu'` to `parkeys` in `from_dict`, and add the `ncpu >= 1` check to
      `validate`.
    - `pypeit/scripts/run_pypeit.py`: add `--ncpu` (type `int`, default `None`)
      to `get_parser`, and pass `ncpu=args.ncpu` to the `PypeIt(...)`
      instantiation in `main`.
    - `pypeit/pypeit.py`: add `ncpu=None` to `PypeIt.__init__`, apply the
      override **immediately after the existing `redux_path` override and before
      the `.par` file is dumped**, and document it in the class docstring.
    - Add the `test_ncpu_default_and_override` test to
      `pypeit/tests/test_pypeitpar.py` and run it.
    - Regenerate the parameter docs (`add-parameter` / `build-docs` skill).
    - Ask questions in Q&A. Log your work in the Logs section.

3. **Force `Agg` and add the deferred-QA machinery.** Coding-doc **§A.3** and
   **§A.4**.
    - `pypeit/scripts/run_pypeit.py`: `matplotlib.use('Agg', force=True)` at the
      top of `main`, guarded by `if not args.show:`, **before** `from pypeit
      import pypeit` (the backend must be selected before any pypeit module
      imports `pyplot`).
    - `pypeit/qa.py`: add the module-level `_QA_POOL` / `_QA_PENDING` /
      `_QA_MAX_PENDING` state and the `init_qa_pool(ncpu)`,
      `save_figure(fig, outfile, show=False, close=True, **kwargs)` and
      `flush_qa()` functions exactly as sketched in §A.4. `init_qa_pool` must be
      safe to call repeatedly (it is used to *reset* the pool in a forked child
      in PR B).
    - Wire the lifecycle: create the pool in `PypeIt.__init__` (after the `ncpu`
      override), and drain with `qa.flush_qa()` at the end of `calib_all`, at the
      end of `reduce_all`, in `RunPypeIt.main` before `pypeIt.build_qa()`, at the
      end of `exposure.reduce_exposure`, and at the end of
      `pypeit_steps.calib_one`.
    - Do **not** convert any QA call site yet — that is task 4.
    - Ask questions in Q&A. Log your work in the Logs section.

4. **Convert the high-volume QA call sites to `qa.save_figure`.** Coding-doc
   **§A.5**, in the priority order given there:
    1. `pypeit/qa.py`: `arc_tilts_2d_qa`, `arc_tilts_spec_qa`,
       `arc_tilts_spat_qa` (3 PNGs per slit — the biggest single block).
    2. `pypeit/core/wavecal/autoid.py`: `arc_fit_qa`, `arc_fwhm_qa`.
    3. `pypeit/flatfield.py`: `spatillum_finecorr_qa` and the second writer.
    4. `pypeit/core/findobj_skymask.py`: the object trace/profile QA writers.
    5. `pypeit/qa.py`: `spec_flexure_qa`, `spat_flexure_qa`.
    - Obey the three conversion rules in §A.5 without exception: **(a)** capture
      the `Figure` object (convert `plt.figure()` + `plt.subplot(gs[...])` sites
      to bind `fig`); **(b)** every `plt.close('all')` must become
      `plt.close(fig)` — `autoid.py` has three of them and they would destroy
      queued figures; **(c)** all layout calls (`tight_layout`,
      `subplots_adjust`) must run *before* the `save_figure` hand-off, and `dpi`
      must be passed explicitly.
    - Before finishing, `git grep -n "close('all')"` over the reduction path and
      confirm none remain in a converted module.
    - Ask questions in Q&A. Log your work in the Logs section.

5. **Tests, docs and changelog.** Coding-doc **§A.6** and **§A.8**.
    - Add the QA tests to `pypeit/tests/test_qa.py`:
      `test_save_figure_serial`, `test_save_figure_threaded_matches_serial`
      (compare **decoded pixel arrays** via PIL, not raw bytes — matplotlib
      embeds a `Software` PNG chunk), and `test_flush_qa_reraises`. Every test
      must restore `qa.init_qa_pool(1)` before returning.
    - Run the CI-safe unit tests: `pytest pypeit/tests` in the `pypeit` env.
      They must be green.
    - Docs: create the "Running on multiple CPUs" subsection in
      `doc/running.rst` (state honestly that in this PR `ncpu` only affects QA
      figure writing); note the `Agg` backend and concurrent PNG writes in
      `doc/qa.rst`; regenerate `doc/pypeit_par.rst` and
      `doc/help/run_pypeit.rst`; check `doc/sphinx_warnings.out` is empty.
    - Changelog: add the bullets from §A.8 to
      `doc/releases/2.1.0dev.rst` under **Functionality/Performance Improvements
      and Additions** and **Testing** (use the `update-changelog` skill).
    - Ask questions in Q&A. Log your work in the Logs section.

6. **Validate and re-profile; summarize.** Coding-doc **§A.7**.
    - Add an `--ncpu N` option to `pypeitdev/speed_up/scripts/profile_kast_blue.py`
      and `profile_deimos.py`, suffixing the artifact stems with `_ncpu{N}` so the
      existing baselines are not overwritten.
    - Re-profile `shane_kast_blue` at `--ncpu 1` and `--ncpu 4` (~2 min each) and
      `keck_deimos 600ZD_M_6500` at `--ncpu 1` and `--ncpu 4` (~3.5 h each), then
      run `analyze_profile.py` on each.
    - **Expected (from §A.7):** `Agg` alone ≈ 1–4 s off Kast (1–3%) and 100–400 s
      off DEIMOS (1–3%); `Agg` + threaded writes at `--ncpu 4` a further 2–3× on
      the PNG render/encode budget, i.e. up to ~300–500 s more on DEIMOS and
      ~4–6 s on Kast. **Combined PR A target: DEIMOS 12 390 s → ~11 500–11 900 s
      (4–7%); Kast 119 s → ~108–114 s (4–9%).**
    - Confirm the QA PNG *count* is unchanged versus the baselines (971 on
      DEIMOS, 20 on Kast) and spot-check a few PNGs render correctly.
    - Per **Q2**, if the threaded writes give little over `Agg` alone, record the
      measurement and **keep the machinery** — do not escalate to a process pool.
    - Write the before/after numbers into a short
      `Reports/speed_up_results.md` (create it; sections: baseline, after-A) and
      summarize in the Logs section below. Ask any remaining questions in Q&A.

7. **Conflicts.** I have issued a Pull Request on GitHub.  In doing so, I find that there are conflicts
    with with the `develop` branch.  Please resolve the conflicts.  Use Opus if you can.  Log your work.

8. **Q2.** I have answered your question.  Please read it and proceed.  
    Use Opus if you can.  Log your work.  Add to your Memory -- JXP performs all
    git operations.

## Q&A

Claude poses questions here (as `**Qn — title.** body` followed by a `>A:` line);
the user answers inline beneath each `>A:`.

**Q1 — The `pypeit` conda env is broken; use `pypeit14b` instead?** This doc
says to run in the `pypeit` env, but it currently cannot run `run_pypeit` at
all: it holds numpy 2.4.6 + astropy 7.0.2, and astropy 7.0.2 references the
removed `np.in1d` during import (`AttributeError: module 'numpy' has no
attribute 'in1d'`). The `pypeit14b` env (numpy 2.5.0, astropy 8.0.0) is an
editable install of this same working repo — so it sees the `speed_up_qa`
checkout — and `run_pypeit -h` works there. Options: (a) use **`pypeit14b`**
for all PR-A work, including the task-6 re-profiling (note: the June baselines
were profiled in the old `pypeit` env, so a fresh `--ncpu 1` baseline in task 6
is needed anyway); (b) you repair the `pypeit` env (e.g. upgrade astropy to
≥8.0) and I keep using it. I won't modify your env myself. Recommendation:
**(a)** — use `pypeit14b`.
>A:  Yes, use pypeit14b for all work.

**Q2 — Task-6 finding: worker-thread *rendering* is unsafe; I changed the
deferred design to render-on-main-thread + encode-in-workers. OK, or revert to
fully serial writes?** The DEIMOS `--ncpu 4` validation run crashed within
minutes: matplotlib's mathtext parser is process-global and **not
thread-safe**, and rendering figures in worker threads races the main thread's
figure layout (`ParseFatalException: Unknown symbol: \mathdefault` on log-axis
tick labels). Kast passed at `ncpu=4` only because its 20 figures make the
race unlikely; DEIMOS's 971 make it near-certain. Fix, already implemented and
committed (`c39818589`): `qa.save_figure` renders on the main thread
(`fig.canvas.draw()` at the requested dpi) and defers only the **PIL PNG
encode** — the thread-safe, GIL-releasing part, i.e. exactly the "encode half"
your accepted planning answer Q2 anticipated. Verified crash-free (40-figure
log-axis stress at `ncpu=8`), pixel-identical to serial, all tests green;
measured wins: DEIMOS −222.9 s, Kast −5.9 s at `--ncpu 4`. Alternative if you
prefer maximum conservatism: drop threading entirely (deferred writes stay
serial; `ncpu` remains plumbing-only until PR B). Recommendation: keep the
implemented encode-in-workers design.
>A: Let's follow your recommendation.

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-09-09 (Task 1 — branched `speed_up_qa`; grounded; found §A.3 already landed + broken env)

Performed the **1st task under Tasks** (prepare; no production code written).

Branch:
- Working tree was clean (on `int_ids`). Created and checked out **`speed_up_qa`**
  from `speed_up` at HEAD **`5120b01ea`** ("Merge pull request #2193 …
  binospec_ifu_comb_group").
- Note: `speed_up` has moved well past the coding doc's grounding commit
  `e9ed85c1a` — it has since absorbed `develop`, the `int_ids` PR (#2188), the
  Binospec IFU comb-group work (#2193), and the mosaic-weights hotfix (#2189).

Read coding-doc §0 + §1 in full, then the target code on the new branch. Two
**major findings**:

1. **§A.3 (force `Agg`) has already landed** on `speed_up` as commit
   `bbfe2675a` (2026-08-19, "Force headless backend for noninteractive
   reductions"): a shared `scriptbase.configure_matplotlib(show)` helper that
   does `matplotlib.use('Agg', force=True)` unless `--show` — exactly the
   accepted Q1 behavior — called at the top of `RunPypeIt.main`
   (`run_pypeit.py:70`) and also wired into `run_to_calibstep` and
   `reduce_by_step`, with unit tests (`test_scriptbase.py`). Scope impact:
   task 3 reduces to §A.4 only (the deferred-QA machinery); the commit shipped
   **without** any doc/release-notes updates (`doc/qa.rst`, `doc/running.rst`,
   `doc/releases/2.1.0dev.rst` say nothing about `Agg`), so task 5 absorbs that
   documentation. The task-6 "Agg alone" saving is already in any fresh
   baseline, so re-profiling must compare against a **new** `--ncpu 1` run, not
   the June numbers.
2. **The `pypeit` conda env cannot run `run_pypeit`** (numpy 2.4.6 +
   astropy 7.0.2; astropy trips on the removed `np.in1d` at import). The
   `pypeit14b` env (numpy 2.5.0 / astropy 8.0.0, editable install of this same
   repo) works: `run_pypeit -h` OK. Raised as **Q1** in Q&A (recommend
   `pypeit14b`); I did not modify either env.

Line-number verification vs coding-doc §0.2/§1 (current branch):
- `pypeit/qa.py` — **matches exactly**: imports 6–19 (`pyplot` at 13),
  `set_qa_filename` 27–141, `arc_tilts_2d_qa`/`_spec_`/`_spat_` at 544/610/693
  with `plt.savefig` at 600/684/750 (each followed by bare `plt.close()` +
  `plt.rcdefaults()`), `spec_flexure_qa` savefig 843/914, `spat_flexure_qa`
  (fig-based) savefig 1098. `arc_tilts_2d_qa` already binds
  `fig, ax = plt.subplots(...)` (line 578) — conversion rule (a) is easy there.
- `pypeit/core/wavecal/autoid.py` — `arc_fit_qa` at 43 with `plt.close('all')`
  at **67 and 200** and savefigs at 152/199 (matches). `arc_fwhm_qa` (207) has
  **grown**: `plt.close('all')` at 242, but its savefig is now at **288**
  (multi-spatial FWHM plotting added) — the doc's "~245" end is stale.
- `pypeit/par/pypeitpar.py` — `ReduxPar` shifted **+12**: class 2770, `__init__`
  2782–2784 (signature ends `chk_version=None`), `quicklook` block ends ~2814,
  `from_dict`/`parkeys` 2885/2889–2890, `validate` 2902.
- `pypeit/pypeit.py` — `PypeIt.__init__` 61–145; `redux_path` override at 87–88
  and `.par` dump at 92–94 (matches); `build_qa` 157–166; `calib_all` 168–210
  (serial det loop 200–207, `print_end_time()` at 210); `reduce_all` 212–251
  (ends `print_end_time()` at 251). `qa` already imported at line 18.
- `pypeit/scripts/run_pypeit.py` — `get_parser` 38–65 (`-c/--calib_only` at
  62–63, matches); `main` now 67–105: `configure_matplotlib` at 70,
  `PypeIt(...)` at 90–93, `build_qa()` at 103.
- `pypeit/exposure.py` — small drift (+1 to +36): `adjust_for_slitmask` 22,
  `process_exposure` 85, `findobj_on_exposure` 141, `extract_exposure` 267,
  `reduce_exposure` 351, `save_exposure` 504.
- `pypeit/pypeit_steps.py` — larger drift: `calib_one` **166**,
  `process_one_det` **253**, `extract_det` **730**.
- `pypeit/flatfield.py` — larger drift: `spatillum_finecorr_qa` **2423**,
  savefigs **2532 and 2600**.
- `pypeit/core/findobj_skymask.py` — savefigs 1456/1621, each followed by a
  `plt.close('all')` (1457, 1622) that task 4's rule (b) must also convert.
- Confirmed still true on current HEAD: **no
  multiprocessing/ThreadPool/ncpu/nproc anywhere in `pypeit/`** (PR A introduces
  the first concurrency), and `spat_flexure_qa` is the only already-fig-based
  writer of the set.

Environment: `run_pypeit -h` **fails** in `pypeit` (see Q1); works in
`pypeit14b`; both envs import this repo editable
(`pypeit 2.0.2.dev584+gde3919738`).

Next (task 2, after Q1 is answered): add `ncpu` to `ReduxPar` + the `--ncpu`
CLI flag + the `PypeIt.__init__` override, with the `test_ncpu_default_and_override`
test and regenerated parameter docs.

### 2026-09-09 (Task 2 — `[rdx] ncpu` parameter + `run_pypeit --ncpu` flag; committed)

Performed the **2nd task under Tasks** (coding-doc §A.1 + §A.2), on branch
`speed_up_qa`. Q1 (env) is still unanswered; since only `pypeit14b` can run
anything, I used it for the test/verification runs — consistent with the Q1
recommendation. No env was modified.

Code changes (two commits, per the coding doc's §6 commit sequence):

- **`79d3a55c0` — "Add [rdx] ncpu parameter"**
  - `pypeit/par/pypeitpar.py` (`ReduxPar`): `ncpu=None` added to the `__init__`
    signature; `defaults['ncpu'] = 1`, `dtypes['ncpu'] = int`, and the §A.1
    description block (placed after the `quicklook` block); `'ncpu'` added to
    `parkeys` in `from_dict`; `validate` now raises `ValueError` for
    `ncpu < 1`.
  - `pypeit/tests/test_pypeitpar.py`: added `test_ncpu_default_and_override`
    (default is 1; `from_dict` round-trip gives 4; `ncpu=0` raises).
  - `doc/pypeit_par.rst` regenerated via `doc/scripts/build_par_rst.py`.
    **Note:** the diff is larger than one row because the generated file was
    *stale on the branch* — an unrelated `lamps_wvrng` parameter had landed in
    code without a doc regen, and the long `ncpu` description reflows the
    ReduxPar table borders. The regenerated file now matches the code.
- **`cfeb4007b` — "Add run_pypeit --ncpu; plumb through PypeIt.__init__"**
  - `pypeit/scripts/run_pypeit.py`: `--ncpu` (type `int`, default `None`) added
    to `get_parser` after `-c/--calib_only`; `ncpu=args.ncpu` passed to
    `PypeIt(...)` in `main`.
  - `pypeit/pypeit.py`: `ncpu=None` added to `PypeIt.__init__`; the override
    (`self.par['rdx']['ncpu'] = ncpu` when not None) applied immediately after
    the `redux_path` override and **before** the `_UTC_*.par` dump, so the
    recorded parameter file reflects what ran; `ncpu` documented in the class
    docstring `Args:` block.

Verification (all in `pypeit14b`):
- `pytest pypeit/tests/test_pypeitpar.py` — **34 passed** (incl. the new test).
- `PypeItPar().to_config()` shows `ncpu = 1` under `[rdx]`.
- Full config-machinery parse: `PypeItPar.from_cfg_lines(... 'ncpu = 4')`
  yields `int` 4; `ncpu = 0` is rejected with "ncpu must be a positive
  integer."
- `run_pypeit -h` shows the new `--ncpu` option with its help text.
- `doc/pypeit_par.rst` line 678 carries the new `ncpu` row.

Deferred to task 5 (as planned there): regenerating `doc/help/run_pypeit.rst`
and the release-notes bullets. No new questions — **Q1 in Q&A still awaits an
answer** (matters most for task 6 profiling; recommendation stands:
`pypeit14b`).

### 2026-09-09 (Task 3 — deferred-QA machinery added and wired; §A.3 was already on the branch)

Performed the **3rd task under Tasks** (coding-doc §A.3 + §A.4) on
`speed_up_qa`. **Q1 answered: use `pypeit14b` for all work** — done so here.

- **§A.3 (force `Agg`): no code change needed.** As found in task 1, it already
  landed on `speed_up` as commit `bbfe2675a`
  (`scriptbase.configure_matplotlib(args.show)` at `run_pypeit.py:70`, exactly
  the accepted `--show`-guarded behavior). Its docs/changelog still land in
  task 5.
- **§A.4 machinery** (commit **`31f822c9f`**, "Defer QA figure writes to an
  optional thread pool"):
  - `pypeit/qa.py`: added `_QA_POOL` / `_QA_PENDING` / `_QA_MAX_PENDING` (=16)
    and `init_qa_pool(ncpu)`, `save_figure(fig, outfile, show, close,
    **kwargs)`, `flush_qa()` exactly as sketched in §A.4 (thread cap
    `min(ncpu, 8)`; figures created/closed only on the main thread; only
    `Figure.savefig` goes to workers; `init_qa_pool` re-callable for the PR-B
    forked-child reset; `flush_qa` re-raises worker exceptions on the main
    thread). `qa.py` deliberately imports nothing from `pypeit` (module warns
    it must not import `log`), so no circular-import risk.
  - Lifecycle wiring: pool created in `PypeIt.__init__` immediately after the
    `--ncpu` override (`qa.init_qa_pool(self.par['rdx']['ncpu'])`); drains
    added at the end of `PypeIt.calib_all` and `PypeIt.reduce_all` (both before
    `print_end_time()`), in `RunPypeIt.main` before `pypeIt.build_qa()` (HTML
    never races PNGs), at the end of `exposure.reduce_exposure` (before its
    return), and at the end of `pypeit_steps.calib_one` (before `return
    caliBrate`) — the last two bound in-flight figures to ~one detector's
    worth. `exposure.py` and `pypeit_steps.py` gained `from pypeit import qa`.
  - **No QA call site converted yet** (that is task 4), so at any `ncpu` the
    pool currently sits idle — behavior is provably unchanged.

Verification (in `pypeit14b`):
- Functional smoke test: serial path writes in-line with `_QA_POOL is None`;
  threaded path (`ncpu=4`, 20 figures > `_QA_MAX_PENDING`, exercising the
  auto-flush) produced 20 PNGs **pixel-identical** to the serial reference
  (PIL-decoded arrays); no matplotlib figures left open after `flush_qa`;
  repeated `init_qa_pool(4/1/2)` re-inits clean; `save_figure(fig, None)`
  closes without writing.
- `pytest pypeit/tests/test_pypeitpar.py pypeit/tests/test_scriptbase.py` —
  **36 passed**; `run_pypeit -h` OK.

No new questions. Next: task 4 — convert the high-volume per-slit QA call
sites to `qa.save_figure` (rules a/b/c), including the two extra
`plt.close('all')` in `findobj_skymask.py` found in task 1.

### 2026-09-11 (Task 4 — converted all §A.5 QA writers; found+fixed a pyplot-state contamination bug)

Performed the **4th task under Tasks** (coding-doc §A.5) on `speed_up_qa`;
commit **`92feaf418`** ("Convert the per-slit QA writers to qa.save_figure").
No new Q&A answers were pending (Q1 already applied: all runs in `pypeit14b`).

Conversions (all five §A.5 priorities, rules a/b/c obeyed at every site):
1. `pypeit/qa.py` — `arc_tilts_2d_qa` / `arc_tilts_spec_qa` /
   `arc_tilts_spat_qa`: identical tails collapsed to
   `save_figure(fig, outfile, show=show_QA, dpi=400)`; `arc_tilts_spec_qa`
   now binds `fig = plt.figure(...)` (rule a).
2. `pypeit/core/wavecal/autoid.py` — `arc_fit_qa` (both the `ids_only` and the
   full-panel exits, dpi 800/400 preserved) and `arc_fwhm_qa` via
   `qa.save_figure(fig, outfile, show=…)`; added `from pypeit import qa`
   (no circular import — `qa.py` imports nothing from pypeit). The two
   `plt.close('all')` at the *starts* of these functions (lines 67, 242) were
   **removed** rather than converted — they close nothing the function owns and
   would destroy queued figures; the one at the end (line 200) became the
   per-figure close inside `save_figure`.
3. `pypeit/flatfield.py` — `spatillum_finecorr_qa` and `detector_structure_qa`
   (the "second writer"): bound `fig`, converted tails (dpi=400,
   `show=outfile is None` preserving show-vs-save semantics). The
   "Saved QA" log line became **"Saving QA" before the hand-off** — under
   deferred writes the file is not yet on disk at log time.
4. `pypeit/core/findobj_skymask.py` — `objfind_QA` and `objtrace_QA`,
   including their trailing `plt.close('all')` calls (rule b).
5. `pypeit/qa.py` — `spec_flexure_qa` (both per-slit figures: correlation +
   sky-lines; original writers pass **no dpi** and that was preserved) and
   `spat_flexure_qa` (kept its `debug → plt.show()` branch; save path via
   `save_figure(fig, outfile, dpi=200)`).

**Bug found by the serial-vs-threaded integration test, and fixed:**
`arc_tilts_spat_qa` PNGs differed between `ncpu=1` and `ncpu=4` when followed
by `objfind_QA`. Cause: `objfind_QA` plots via bare pyplot state (`plt.plot`
with no figure of its own); with closes deferred, `gcf()` returned the still
*pending* `arc_tilts_spat_qa` figure and `objfind_QA` drew into it before the
worker rendered it. (The old code was accidentally protected by every writer
closing everything.) Two-part fix, both in commit `92feaf418`:
- `qa.save_figure` now **deregisters a queued figure from pyplot immediately**
  (`plt.close(fig)` before submit, when `close=True`): the `Figure` object
  stays alive via `_QA_PENDING` and Agg renders a closed figure identically in
  the worker (verified byte-identical). Pending figures are thereby invisible
  to any later pyplot-state code — a general defense, not just for
  `objfind_QA`.
- `objfind_QA` now creates its own figure (`fig = plt.figure()`) instead of
  drawing into whatever is current.

Verification (in `pypeit14b`):
- Integration test: `arc_tilts_{2d,spec,spat}_qa` + `objfind_QA` rendered in
  one sequence at `ncpu=1` vs `ncpu=4` → all four PNGs **pixel-identical**
  (PIL-decoded); no figures leaked. (First run flagged `objfind` differing —
  a test artifact: the random SNR vector was regenerated between runs; fixed
  data → identical.)
- `git grep "close('all')"` over the four converted modules → **none remain**;
  the only remaining `savefig` in them is `autoid.py:358` `pp.savefig` — a
  `PdfPages` writer, explicitly out of scope.
- `pytest test_wvcalib/test_flatfield/test_fiberflatfield/test_wavetilts/
  test_wave/test_wavemodel` — **14 passed**.

Note for task 5: add a regression test for the pyplot-contamination fix (a
pyplot-state plot after a queued figure must not alter the queued PNG).
Coadd/sensfunc/telluric/PdfPages writers untouched, per scope.

### 2026-09-11 (Task 5 — tests green (696), docs + changelog written; warnings file empty)

Performed the **5th task under Tasks** (coding-doc §A.6 + §A.8) on
`speed_up_qa`; commit **`8b6eb3c6f`** ("Add QA-write tests; document ncpu and
the Agg backend"). No new Q&A answers pending (only Q1, already applied —
everything run in `pypeit14b`).

Tests (`pypeit/tests/test_qa.py` — the file already existed with
`test_get_dimen`; extended, not replaced):
- `test_save_figure_serial` — in-line write at `ncpu=1`, figure closed.
- `test_save_figure_threaded_matches_serial` — 8 deferred writes at `ncpu=4`
  pixel-identical to the serial reference (PIL-decoded arrays, per §A.6 — raw
  bytes differ via matplotlib's `Software` PNG chunk); fixed RNG seed.
- `test_save_figure_pending_isolated_from_pyplot` — the **task-4 regression
  test**: a bare `plt.plot` + `plt.close('all')` issued while a figure is
  queued must not alter the queued PNG.
- `test_flush_qa_reraises` — a failed background write surfaces on the main
  thread.
- Every test restores `qa.init_qa_pool(1)`; module forces `Agg`.
- **Full CI-safe suite: `pytest pypeit/tests` → 696 passed** (~5.3 min).

Docs (§A.8):
- `doc/running.rst`: new subsection **"Running on multiple CPUs"**
  (`.. _run-pypeit-ncpu:`) — states honestly that in this PR `ncpu` only
  affects QA figure writing, and warns that future detector-level use scales
  peak memory.
- `doc/qa.rst`: paragraph on the forced `Agg` backend (with the `--show`
  exception) and concurrent, output-identical PNG writes; links to
  `run-pypeit-ncpu`.
- `doc/releases/2.1.0dev.rst`: two bullets under
  **Functionality/Performance Improvements and Additions** (the `Agg` change —
  crediting that it covers `run_to_calibstep`/`reduce_by_step` too, previously
  undocumented from commit `bbfe2675a` — and the `ncpu` parameter + `--ncpu`
  flag), plus one under **Testing** for the new QA tests. `CHANGES.rst`
  untouched (its own header deprecates it, per coding-doc §A.8).
- Regenerated `doc/help/run_pypeit.rst` (now shows `--ncpu`);
  `doc/pypeit_par.rst` unchanged since its task-2 regeneration.
- Doc build (`make htmlnoex` then `htmlonly`, sphinx from `pypeit14b`,
  `PYPEIT_DEV` set): **build succeeded; `doc/sphinx_warnings.out` is empty.**
  Two warning classes were fixed along the way: (1) my `qa.py` docstrings
  referenced `matplotlib.figure.Figure(.savefig)`_ targets missing from
  `doc/include/links.rst` — added both targets there; (2) a **pre-existing**
  warning unrelated to PR A: the `int_ids` spectrograph API page (from PR
  #2188) was never registered in `doc/api/pypeit.spectrographs.rst` — added
  the toctree line + the generated `pypeit.spectrographs.int_ids.rst`.

Housekeeping notes:
- The apirst regeneration also showed ~10 *other* stale generated docs on the
  branch (`doc/help/pypeit_obslog.rst` [-G removal], `pypeit_ql`, `setup`,
  `trace_edges`, `view_fits`, `chk_for_calibs`, `cache_github_data`,
  `doc/include/{inst_detector_table,spectrographs_table}.rst`,
  `doc/api/pypeit.coadd2d.rst`). **Restored, not committed** — unrelated to
  PR A; worth a separate docs-regen commit on `speed_up` at some point.
- The full pytest run left untracked test outputs under
  `pypeit/data/tests/{REDUX_OUT_TEST,shane_kast_blue_A}/` — left untouched
  (pre-existing suite behavior, not git-ignored).

No new questions. Next: task 6 — add `--ncpu` to the two profiling scripts,
re-profile Kast + DEIMOS at `--ncpu 1` and `--ncpu 4` (fresh `ncpu=1` baselines;
the June profiles predate the `Agg` commit and used the old env), and write
`Reports/speed_up_results.md`.

### 2026-09-21 (Task 6 — re-profiled; found+fixed a matplotlib thread-safety crash; results doc written)

Performed the **6th task under Tasks** (coding-doc §A.7). No new Q&A answers
were pending at start (Q1 applied throughout: everything in `pypeit14b`).
**PR A is complete pending the new Q2 answer below.**

Profiling infrastructure:
- Added `--ncpu N` to `pypeitdev/speed_up/scripts/profile_kast_blue.py` and
  `profile_deimos.py`: passes `--ncpu N` through to `run_pypeit` and suffixes
  all artifact stems with `_ncpu{N}` (June baselines untouched); `ncpu`
  recorded in the runmeta JSON. `analyze_profile.py --stem` handles the
  suffixed stems unchanged.
- Runs chained strictly sequentially (each pair shares a redux dir and
  cold-cleans it); after each run the QA PNG count was recorded and 3 sample
  PNGs preserved (`Reports/task6_qa_samples/`, `task6_png_counts.txt`) before
  the next cold-clean.

**Crash found on the first `--ncpu 4` DEIMOS run — matplotlib rendering is not
thread-safe.** The run died minutes in, on the main thread inside `calib_one`
QA layout: `ParseFatalException: Unknown symbol: \mathdefault` — matplotlib's
mathtext parser is process-global; worker-thread renders race main-thread
`tight_layout` text metrics and corrupt the parser (log-axis tick labels like
`$10^5$` trigger the parse). Kast had passed at `ncpu=4` only by low volume
(20 vs 971 figures). **Fix (commit `c39818589`, "Render QA on the main thread;
hand only the PNG encode to workers")**: `qa.save_figure`'s deferred path now
renders on the calling thread (`fig.set_dpi(dpi)`; `fig.canvas.draw()`; copy
`buffer_rgba`) and submits only a PIL PNG encode (`_encode_png`) to the pool —
the thread-safe, GIL-releasing "encode half" that the accepted planning Q2
predicted as the recoverable fraction. Non-PNG outputs / `close=False` /
kwargs beyond `dpi` fall back to synchronous `savefig`; `_QA_MAX_PENDING`
lowered 16 → 4 (pending items now hold full RGBA rasters, >100 MB for the
large tilt figures). Docs/release-notes wording updated to match. Verified:
`test_qa.py` 5 passed; the four converted writers pixel-identical at
`ncpu=1` vs `ncpu=4`; a 40-figure log-axis mathtext stress at `ncpu=8` (the
crash trigger) runs clean. Raised as **Q2 in Q&A** (keep vs revert to fully
serial; recommend keep).

Results (all cold, single runs, `pypeit14b`, branch `speed_up_qa` at
`c39818589`; full write-up in **`Reports/speed_up_results.md`**, new file with
baseline / after-A sections):
- **Kast**: `--ncpu 1` **88.5 s**, `--ncpu 4` **82.6 s** (−5.9 s, −6.7%);
  20 QA PNGs in both ✓.
- **DEIMOS**: `--ncpu 1` **10 056.9 s**, `--ncpu 4` **9 834.0 s** (−222.9 s,
  −2.2%); 971 QA PNGs in both ✓; rc 0.
- Mechanism confirmed in the profiles: main-thread `ImagingEncoder.encode`
  self-time is **275.6 s** on DEIMOS at `ncpu=1` and moves to workers at
  `ncpu=4` — the −222.9 s wall gain matches the Q2 "encode half" expectation
  (residue = the deliberately-unconverted synchronous writers).
- vs the June baselines (121.2 s / 12 393.8 s): fresh `ncpu=4` numbers are
  −31.8% (Kast) and −20.7% (DEIMOS), but that bundles the `Agg` fix with the
  env change (June's `pypeit` env no longer runs) and 3 months of branch
  drift — recorded honestly in the results doc; the clean PR-A-isolated
  measurement is the ncpu4-vs-ncpu1 delta.
- Sample QA PNGs from all four runs decode correctly (PIL), sizes/modes sane.
- Per accepted planning-Q2: measurement recorded, machinery kept, no
  escalation to a process pool.

Branch state: `speed_up_qa` = `5120b01ea` + 6 commits (`79d3a55c0`,
`cfeb4007b`, `31f822c9f`, `92feaf418`, `8b6eb3c6f`, `c39818589`); PypeIt tree
clean. Dev-suite artifacts (updated profile scripts, `speed_up_results.md`,
`*_ncpu{1,4}` profiles/logs/analyses, QA samples) left for your commit, per
your pattern of committing the dev-suite side.

**PR A work is done.** Remaining: your answer to Q&A **Q2** (keep the
encode-in-workers design vs revert to fully serial), then open the PR
(`speed_up_qa` → `speed_up`) and start B.md task 1 (branch `speed_up_detpar`
from `speed_up_qa`).

### 2026-09-22 (Task 7 — resolved the develop-merge conflicts; merged, tested, pushed)

Performed the **7th task under Tasks**: you had opened the PR, found conflicts
with `develop`, and started `git merge origin/develop` on `speed_up_qa`
(MERGE_HEAD `f3a1f1d27`, ~288 incoming commits; three `UU` paths). Conflict
resolution was performed by an **Opus** agent per the task instruction; I
briefed it, reviewed the result, and pushed.

Resolutions:
- `doc/releases/2.1.0dev.rst` — two conflict hunks, both merged **additively**
  (develop's bullets — hamspec/echelle parameters, `core.meta.airmass`,
  pinhole-decker tracing, datacube tests — followed by ours: the Agg-backend
  bullet, the `[rdx] ncpu`/`--ncpu` bullet, and the QA-tests bullet).
- `doc/pypeit_par.rst` — **regenerated** from the merged code
  (`doc/scripts/build_par_rst.py`, pypeit14b) rather than hand-merged;
  verified it carries both our `ncpu` row and develop's new parameters
  (`ech_direct_cc`, `force_center_obj`).
- `doc/help/run_pypeit.rst` — **regenerated** via
  `doc/scripts/write_script_help.py`; that script rewrites all help files, so
  the other 7 it touched were restored from the index to keep the merge
  minimal. Verified `--ncpu` present plus develop's new spectrograph names.

Auto-merge sanity checks — **no manual code fixes needed**; all PR-A machinery
survived develop's 288 commits: `qa.py` (`_encode_png`/`init_qa_pool`/
`save_figure`/`flush_qa`; develop lifted the old "cannot import log"
restriction and its `from pypeit import log` coexists with our code),
`pypeit.py` (ncpu override + pool init + both flushes), `exposure.py` /
`pypeit_steps.py` flushes, `run_pypeit.py` `--ncpu`, `ReduxPar.ncpu`, and all
seven converted call-site modules. Conflict-marker sweep clean.

Validation (pypeit14b): targeted `test_qa/test_pypeitpar/test_scriptbase` →
41 passed; **full `pytest pypeit/tests` → 741 passed, 4 skipped, 0 failed**
(the suite grew from 696 with develop's new tests); `run_pypeit -h` shows
`--ncpu`.

Merge commit **`0eca4ee67`** (parents `c39818589` + `f3a1f1d27`), tree clean,
**pushed to `origin/speed_up_qa`** — the PR should now be conflict-free.

Note: **Q&A Q2 (encode-in-workers design vs fully-serial fallback) is still
unanswered** — the merged branch carries the implemented encode-in-workers
design.

### 2026-09-22 (Task 8 — Q2 answered: encode-in-workers design accepted; PR A closed out)

Performed the **8th task under Tasks**: read the **Q2 answer — "Let's follow
your recommendation"** — i.e., **keep the encode-in-workers design** (render on
the main thread, PIL PNG encode in worker threads).

No code change was needed: the accepted design is exactly what is on the
branch (commit `c39818589`, carried through the develop merge `0eca4ee67` and
already pushed to `origin/speed_up_qa`). Verified the tree state: on
`speed_up_qa`, clean, in sync with origin.

Housekeeping: noted the Q2 acceptance in `Reports/speed_up_results.md` (the
"Design deviation" section now records the decision and date).

**All A.md Q&A questions are now answered and PR A is complete**: 7 commits on
top of `speed_up` (6 feature/fix commits + the develop merge), full unit suite
green (741 passed / 4 skipped post-merge), docs/changelog in, re-profiling
recorded in `Reports/speed_up_results.md` (Kast −6.7%, DEIMOS −222.9 s at
`--ncpu 4` vs the fresh serial baselines; PNG counts and content verified).
The PR on GitHub is conflict-free and ready for review/merge. Next in the
stack: **B.md task 1** — branch `speed_up_detpar` from `speed_up_qa` for
detector parallelism v1.
