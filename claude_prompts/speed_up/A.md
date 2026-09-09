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
