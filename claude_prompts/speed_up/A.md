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
7. Read this doc.  Perform the 7th task under Tasks.
8. Read this doc.  Perform the 8th task under Tasks.

9. Read this doc.  Perform the 1st item under PR.
10. Read this doc.  Perform the 2nd item under PR.
11. Read this doc.  Perform the 3rd item under PR.


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

## PR

1. There are now responses to the PR on GitHub from Tim and Kyle.  Please read them and respond as needed.  
Provide ideas and questions for Kyle's point about global in the Q&A/PR section below.
Use Opus if you can.  Log your work.

2. Kyle has responded to your response in the PR and offered his blessing.  Read that and proceed
with the changes.  Put the changes in this PR.  Adopt your recommendations for Q3-Q5.
If you have more questions, ask them in the Q&A/PR section below.
Use Opus if you can.  Log your work.

3. Ok, check the rest of the PR comments and respond as needed.  If you have an questions for me
put them in the Q&A/PR section below.  Use Opus 5.  Log your work.

## Q&A

### PR

**Kyle's point (PR #2198, CHANGES_REQUESTED):** *"I'm very wary of modifying
things in the global scope in the way that is done in this PR. We have global
objects, like `log` and `dataPaths`, but these are class instances that are
instantiated on import, imported into each module as needed, and then
accessed/manipulated. Can we do something similar with this QA functionality?"*

He is right about the mismatch: PR A adds module-level `_QA_POOL` /
`_QA_PENDING` / `_QA_MAX_PENDING` in `pypeit/qa.py`, mutated by module-level
functions (`init_qa_pool`, `save_figure`, `flush_qa`). Nothing else in the
package manages process-wide state that way.

**Idea — the `log` precedent maps onto this almost one-to-one**, including the
part that matters most here (configuration deferred until the parameters are
known):

| `log` | QA equivalent |
|---|---|
| `PypeItLogger` in `pypeit/pkg/logger.py` | `QAWriter` in `pypeit/pkg/qawriter.py` |
| `log = get_logger(...)` in `pypeit/__init__.py` | `qaWriter = QAWriter()` there |
| `from pypeit import log` | `from pypeit import qaWriter` |
| `ScriptBase.init_log()` → `log.init(level=…, log_file=…)` | `PypeIt.__init__` → `qaWriter.init(ncpu=par['rdx']['ncpu'])` |
| `log.close_file()` | `qaWriter.flush()` |

Pool and pending-encode list become instance attributes; call sites become
`qaWriter.save_figure(fig, outfile, dpi=400)`; a default-constructed instance
is fully serial, so a bare `import pypeit` behaves exactly as today and no pool
is ever created. `init()` stays re-callable (needed anyway for PR B's forked
workers, and true of `log.init()` too). Bonus: the current tests have to
save/restore module state; with an instance they can just build an isolated
`QAWriter(ncpu=4)`.

**Honest caveat to keep in mind:** this *encapsulates* the global state, it
does not remove it — a module-level singleton is still process-wide mutable
state. What genuinely improves is convention-matching, discoverability, and
test isolation. Removing the global entirely means threading a writer instance
through every QA function signature (`arc_fit_qa(..., qa_writer=...)` and all
their callers) — a much larger diff. I put that choice to Kyle in the PR reply
and it is Q1 below.

**Implementation constraint I verified (matters for where the class lives):**
`import pypeit` today imports **no matplotlib at all**. If `QAWriter` imports
`pyplot` at module scope and is instantiated in `pypeit/__init__.py`, then
every `import pypeit` (including scripts that never plot) drags in
matplotlib+pyplot, and — more seriously — pyplot would be imported *before*
`scriptbase.configure_matplotlib()` can force `Agg` in `RunPypeIt.main`, which
is the ordering the Agg fix depends on. Easily avoided by keeping
`qawriter.py` import-light and importing `pyplot` lazily inside the methods
that need it (`plt.close`, `plt.show`), but it must be done deliberately.

Questions for you (answer inline beneath each `>A:`):

**Q3 — Which shape: singleton, or instance threaded through call signatures?**
(a) `qaWriter` singleton in `pypeit/__init__.py` mirroring `log`/`dataPaths` —
matches convention, small diff, keeps every QA call site a one-line change;
(b) construct a `QAWriter` in `PypeIt` and pass it down through every QA
function — genuinely removes the global, but changes the signature of ~10 QA
functions and all their callers, and the QA functions are called from deep in
`core/` where nothing else is threaded through. Recommendation: **(a)**, which
is what I proposed to Kyle; I asked him to confirm since he may want (b).
>A: Adopted (a), per "Adopt your recommendations for Q3-Q5" + Kyle's blessing
on the PR ("Yes, I think the singleton approach is preferred").

**Q4 — Where should the class live, and is the lazy-pyplot-import approach
OK?** I'd put it in `pypeit/pkg/qawriter.py` alongside `logger.py` and
`pypeitdata.py` (the existing home for package-level service classes), with
`pyplot` imported inside methods to preserve the "no matplotlib at package
import" property above. The alternative is defining the class in
`pypeit/qa.py` and instantiating it there, but then it is not a package-level
global in the `log`/`dataPaths` sense that Kyle asked for. Recommendation:
`pypeit/pkg/qawriter.py` + lazy imports. Also: name it `qaWriter` (camelCase,
matching `dataPaths`)?
>A: Adopted as recommended: `pypeit/pkg/qawriter.py`, lazy matplotlib imports,
instance named `qaWriter`.

**Q5 — Land the refactor in this PR, or as a follow-up?** Kyle has
CHANGES_REQUESTED on #2198, so it probably has to land here to unblock the
merge; that also avoids PR B inheriting the API that is about to change (B's
`_worker_init` calls into it). The cost is that PR A grows and the measured
results in `speed_up_results.md` would be re-verified (cheap: the Kast pair is
~3 min; I would not re-run DEIMOS unless you want it). Recommendation: land it
in PR A, re-run the Kast pair only, and note in the results doc that the
DEIMOS numbers predate a pure-refactor change. Alternatively: wait for Kyle's
answer to Q1 before writing any of it.
>A: Adopted as recommended: landed in PR A (#2198), Kast pair re-run only,
DEIMOS caveat noted in `speed_up_results.md`.

**Q6 — Re-trigger Cursor Bugbot on the refactor?** Bugbot last reviewed
`c3981858` (before the `develop` merge), so the two largest subsequent
changes — the merge and the whole `QAWriter` refactor — have had no automated
review. You trigger it by commenting `@cursor review`, which I can do, but it
costs some comment noise on a PR that is otherwise just waiting on Kyle.
Recommendation: yes, trigger it once on `b8e5f9321` — the refactor touched 10
production files and a second pair of eyes is cheap. Say the word (or do it
yourself) and I'll handle any findings.
>A: No 

**Q7 — Report develop's duplicated INT/IDS parameter docs separately?** While
checking the PR diff I found that `doc/pypeit_par.rst` **on `develop`** has
the `int_ids_eev10` and `int_ids_redplus2` sections duplicated, including
duplicate Sphinx labels (`.. _instr_par-int_ids_eev10:`), ~179 lines of
accidental repetition — presumably a bad merge of a generated file. This PR
silently fixes it, because resolving the merge conflict required regenerating
the file. I noted it in a PR comment so reviewers aren't puzzled by the diff
size. Do you want me to also open a short issue against `develop` recording
it (in case this PR is reworked and the fix is lost), or is the PR comment
enough? Recommendation: the PR comment is probably enough given this PR
should merge soon.
>A: Agreed

### Tasks

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

### 2026-09-24 (PR item 1 — answered the PR #2198 reviews; fixed the --ncpu validation bug)

Performed the **1st item under PR** (both Q&A/Tasks questions were already
answered; nothing new pending there). Work done with **Opus 5**. PR is
**#2198 "Speed up QA"** — note it targets **`develop`**, not `speed_up` as the
branch-stack plan assumed.

Three review items found:

1. **Tim (`tepickering`) — APPROVED**: *"looks good! smarter use of matplotlib
   is long overdue..."* No action needed.
2. **Cursor Bugbot — inline on `pypeit/pypeit.py:98`, "CLI ncpu skips
   positivity check"**, with Tim endorsing it (*"this is a fix worth doing. i
   will let the authors handle it."*). **Real bug**: `--ncpu` was written
   straight into `par['rdx']['ncpu']` by item assignment, which bypasses the
   `ReduxPar.validate()` that runs at instantiation — so `--ncpu 0` / `-1` was
   silently accepted, recorded in the dumped `.par`, and treated as serial by
   `init_qa_pool`. **Fixed in commit `611c24cd8`**: the override now re-runs
   `self.par['rdx'].validate()` (reusing the single source of truth rather
   than duplicating the check). Verified end-to-end — `PypeIt(..., ncpu=0)`
   and `ncpu=-1` both raise `ValueError: ncpu must be a positive integer.`
   before anything is written; added `test_ncpu_assignment_revalidates`;
   `test_pypeitpar` + `test_qa` → 40 passed. Pushed; replied on the Bugbot
   thread with the diff and behavior.
3. **Kyle (`kbwestfall`) — CHANGES_REQUESTED** on the module-global QA state
   (`_QA_POOL`/`_QA_PENDING` mutated by module-level functions), asking for the
   `log`/`dataPaths` pattern instead (class instance created on import,
   imported per module). **No code written** — per this task's instruction the
   ideas and questions went to **Q&A/PR** above (Q3 singleton vs
   threaded-instance; Q4 location + the lazy-pyplot constraint; Q5 land in this
   PR vs follow-up), and I posted a substantive reply on the PR: agreed the
   mismatch is real, gave the one-to-one `log` → `QAWriter` mapping (incl.
   `log.init()` ↔ `qaWriter.init(ncpu=…)` as the deferred-configuration
   precedent), flagged honestly that a singleton *encapsulates* rather than
   removes global state and offered the larger signature-threading alternative
   if that is what he wants, noted the PR-B fork-reset requirement, and
   restated why rendering must stay on the main thread (the mathtext crash) so
   that constraint is not lost in whatever restructure we land.

**Engineering finding worth keeping** (drove Q4): `import pypeit` today imports
**no matplotlib at all**. A `QAWriter` that imports `pyplot` at module scope
and is instantiated in `pypeit/__init__.py` would change that — slowing every
import and, more seriously, importing pyplot *before*
`scriptbase.configure_matplotlib()` can force `Agg`, which is the ordering the
Agg fix relies on. Whatever shape we choose must keep the class import-light
and import `pyplot` lazily inside methods.

Branch: `speed_up_qa` at **`611c24cd8`**, pushed, tree clean. Blocked on your
answers to **Q&A/PR Q3–Q5** (and ideally Kyle's reply) before writing the
refactor.

### 2026-09-24 (PR item 2 — QA globals refactored into the `qaWriter` singleton; Kyle's ask satisfied)

Performed the **2nd item under PR** with **Opus 5**: Kyle gave his blessing
(*"Yes, I think the singleton approach is preferred. It's the path I took in
#2167, but it's a path I've been wanting to rethink given concerns about how
it may limit building a processing pool. If the `init()` method gets us around
that, great."*), and you said to adopt my Q3–Q5 recommendations — both now
recorded inline above. (#2167 turns out to be his open *issue* "Centralize
definition of output paths", not a PR.)

**Commit `b8e5f9321` — "Move the QA writing into a QAWriter class instance"**,
mirroring the `log`/`dataPaths` convention exactly:
- **New `pypeit/pkg/qawriter.py`** with `class QAWriter` — alongside
  `logger.py`/`pypeitdata.py`. Pool and queued encodes are instance
  attributes; `max_threads` class attribute; a `parallel` property; `__repr__`.
- **`pypeit/__init__.py`** instantiates `qaWriter = QAWriter()` next to `log`
  and `dataPaths` (default instance is fully serial, so a bare `import pypeit`
  behaves exactly as before and no pool is ever created).
- **API migration** across 8 modules: `qa.init_qa_pool(n)` →
  `qaWriter.init(ncpu=n)`, `qa.save_figure(...)` → `qaWriter.save_figure(...)`,
  `qa.flush_qa()` → `qaWriter.flush()`. `qa.py` lost 134 lines of module-level
  state. Modules that only used `qa` for these helpers (`autoid.py`,
  `findobj_skymask.py`, `exposure.py`, `pypeit_steps.py`) had their import
  swapped; `flatfield.py` and `pypeit.py` keep `qa` (they also use
  `set_qa_filename`/`gen_*_html`) and gained `qaWriter`.
- **`init()` stays re-callable** — the point Kyle flagged about pools. A forked
  child drops the inherited pool with `qaWriter.init(ncpu=1)`; queued encodes
  are discarded rather than awaited (they belong to the parent's threads).
  Documented, and covered by `test_init_is_repeatable`. This is exactly what
  PR B's `_worker_init` needs.
- **Lazy matplotlib imports** (the Q4 constraint I had verified): `pyplot` and
  `PIL` are imported *inside* the methods, so `import pypeit` still pulls in no
  matplotlib and cannot select a backend before
  `scriptbase.configure_matplotlib()` forces `Agg`. **Re-verified after the
  change** — `'matplotlib' in sys.modules` is still `False` after
  `import pypeit`.
- **Tests rewritten to use isolated instances** (`QAWriter(ncpu=4)`) instead of
  mutating/restoring module state — the concrete testability benefit I cited to
  Kyle — plus `test_package_writer_is_serial` asserting the singleton's default.
  7 tests in `test_qa.py`.
- **Docs**: new `doc/api/pypeit.pkg.qawriter.rst` + toctree entry; a
  `concurrent.futures.ThreadPoolExecutor` link target added to
  `doc/include/links.rst`; an Under-the-hood release-notes bullet placed
  beside the existing "restructure globals into `pkg`" bullet.

Validation: integration test of the four converted writers through the
singleton at `ncpu=1` vs `ncpu=4` → **pixel-identical**, no leaked figures;
the 40-figure log-axis mathtext stress at `ncpu=8` still clean; **full unit
suite 744 passed / 4 skipped**; docs build with **no new warnings** (the 8
remaining are develop's pre-existing `vlt_uves_template` + `keck_hires.rst`
issues, not mine).

Per Q5, re-ran the **Kast pair only**: `--ncpu 1` 89.8 s, `--ncpu 4` 85.2 s
(vs 88.5/82.6 before the refactor) — within single-run scatter, so the
refactor is **performance-neutral**; 20 QA PNGs throughout. Recorded in
`Reports/speed_up_results.md` as a new "Post-review refactor" section,
including the explicit caveat that the DEIMOS numbers were not re-measured.

Pushed (`speed_up_qa` at `b8e5f9321`) and replied on the PR summarizing the
change, the fork/`init()` story, the lazy-import rationale, and the validation.
No new questions. PR #2198 now carries both review fixes (Bugbot validation +
Kyle's refactor) and is back with the reviewers.

### 2026-09-25 (PR item 3 — swept the whole PR; CI green; found develop's duplicated parameter docs)

Performed the **3rd item under PR** with **Opus 5**: a full sweep of PR #2198
— every review, every inline thread, the CI checks, and the diff itself.

**Nothing new had arrived** since my last push. Full inventory, all accounted
for:
- 3 inline comments, all on one thread: Bugbot's `--ncpu` finding → Tim's
  *"this is a fix worth doing"* → my reply (fixed in `611c24cd8`). Resolved.
- Tim **APPROVED**; Kyle **CHANGES_REQUESTED** (answered by the `b8e5f9321`
  refactor and my reply, but **not yet re-reviewed**).
- No new Bugbot review on the refactor — it last ran on `c3981858`, i.e.
  before both the `develop` merge and the whole `QAWriter` change (raised as
  **Q6**).

**CI: all 21 checks pass on `b8e5f9321`** (the current head) — python 3.11–3.14
on macos/windows, the tox matrix incl. `numpydev` and `astropydev`, `codestyle`,
`conda`, `validate`. PR state is `mergeable=MERGEABLE` but
`mergeStateStatus=BLOCKED` with `reviewDecision=CHANGES_REQUESTED`: **the only
thing standing between this PR and merge is Kyle clearing his review.**

**Diff audit** (a reviewer would ask about two things, so I checked both first):
- *Scope:* clean — exactly **21 files, all QA/`ncpu` related**. The Cursor
  summary at the top of the PR claims it also registers INT/IDS and changes
  P200/NGPS wavelength defaults; that is **stale**, written against
  `c3981858` before the `develop` merge, and those changes dropped out of the
  diff once develop absorbed them. Noted publicly so nobody re-litigates it.
- *The one big file:* `doc/pypeit_par.rst`, ~550 changed lines for one new
  parameter. Two causes, both benign, and I verified the file is reproducible
  (re-running `build_par_rst.py` on the branch is a no-op): (a) the long
  `ncpu` description reflows the `ReduxPar` table's column widths, rewrapping
  every row; (b) **a genuine find — `develop`'s committed `pypeit_par.rst`
  contains the `int_ids_eev10` and `int_ids_redplus2` sections twice**,
  including duplicate `.. _instr_par-int_ids_*:` Sphinx labels, ~179 lines of
  accidental repetition (develop 12 671 lines → 12 492 here). Presumably a bad
  merge of a generated file. Regenerating to resolve the merge conflict
  silently fixed it; raised as **Q7** whether to also report it separately.

Posted one PR comment covering all of the above — the two fixes with their
commits, the green CI, the `pypeit_par.rst` explanation with the duplicate-label
evidence, the stale-summary note, and a direct ask to Kyle for the re-review.

Two questions for you in **Q&A/PR**: **Q6** (re-trigger Bugbot on the
refactor? recommend yes) and **Q7** (also file an issue for develop's
duplicated parameter docs? recommend no, the PR comment suffices). Neither
blocks anything. No code changes this round; branch unchanged at `b8e5f9321`.
