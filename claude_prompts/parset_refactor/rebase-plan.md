# Rebase `parset_refactor` onto `develop`

## Context

`parset_refactor` diverged from `develop` at commit `a838e303e` and has spent 78
commits rewriting `pypeit/par/parset.py` (the `ParSet` base class) and, as a
consequence, every subclass in `pypeit/par/pypeitpar.py`. Meanwhile `develop`
has advanced **851 commits** past that same point, including real parameter
additions to `pypeitpar.py` that were never carried into the refactor. The
branch needs to be rebased onto current `develop` so it doesn't fall further
behind, the refactored `pypeitpar.py` needs to regain every parameter
`develop` added or changed, and the file should end up stylistically uniform.

The `ParSet` rewrite replaced the old pattern (an `__init__` per subclass that
hand-built four parallel dicts — `defaults`, `options`, `dtypes`, `descr` —
matched by list position) with a single declarative class attribute:
`parameters = {'key': parset.set_parameter_definition(dtype=..., default=...,
options=..., descr=...), ...}`. The base `ParSet.__init__(self, **kwargs)`,
`from_dict`, `to_dict`, `to_config`, and the new `to_header`/`from_header`
(FITS round-trip) methods are now fully generic, so subclasses need no
`__init__` or `from_dict` override at all (only `PypeItPar` still overrides
`from_dict`, for its `baseprocess` merge logic). This pattern is *already*
100% consistently applied across all 41 classes in `pypeit/par/pypeitpar.py`
on `parset_refactor` today — the risk is only in what gets bolted back on
during the rebase.

Because 27 of the 78 commits touch `parset.py`/`pypeitpar.py`, and 8 more are
past `Merge branch 'develop' into parset_refactor` merges, a literal `git
rebase develop` would force re-resolving the same conflicts dozens of times
over. **Decision (confirmed with user): squash the 78 commits into a handful
of logical commits first, then rebase that short history onto `develop`.**
This work happens on a **new branch, `parset_refactor-rebase`**, checked out
from the current `parset_refactor` tip — `parset_refactor` itself is left
completely untouched as the fallback (no tag needed; once the rebase is
verified, the user will retire the old branch themselves).

Several other local branches (`funcpar_update`, `wavemask`, `dmod_update`,
`speclist`, `latinhyp`, `standard_update`, `spec_dm`, `sensfunc_merge`, ...)
are stacked on top of `parset_refactor` at various points in its history.
Rebasing those onto the new, squashed `parset_refactor-rebase` is
**out of scope for this pass** but is addressed as future-work
considerations below, since squashing changes their common ancestor's commit
SHAs.

## Investigation findings

- Merge-base: `a838e303e`. `develop` is 851 commits ahead; `parset_refactor`
  is 78 commits ahead.
- 28 files are touched by both branches (conflict candidates); 63 more files
  are touched only by `parset_refactor` (mostly new `doc/include/parset_*.rst`
  files and small call-site updates for the new `ParSet` API — no conflict
  risk).
- Of the 28 overlapping files, **only `pypeit/par/pypeitpar.py` is a true
  content conflict** (both sides rewrote it substantially: `develop`
  +392/-36, `parset_refactor` +5037/-5285). ~24 of the others have trivial
  1-5 line diffs on `parset_refactor`'s side (pure old→new `ParSet`
  attribute-access syntax fixes, e.g. `par.descr['x']` →
  `par.parameters['x']['descr']`) against much larger, unrelated `develop`
  changes — resolution is "take `develop`'s version, reapply the small
  syntax tweak."
- Three files need real (but mechanical) reconciliation:
  - `pypeit/core/qa.py` — `develop` deleted/moved it to `pypeit/qa.py`;
    reapply `parset_refactor`'s one-line type-hint tweak at the new location.
  - `pypeit/scripts/collate_1d.py` — `develop` moved its logic into a new
    `pypeit/collate.py`; take `develop`'s structure, reapply
    `parset_refactor`'s `ParSet`-API-usage edits.
  - `pypeit/telescopes.py` — `develop` only added two new telescope classes;
    `parset_refactor` reformatted the whole file for the new
    `pypeitpar.TelescopePar` import style and bare `super().__init__()`.
    Take `develop`'s new classes, reformat them into `parset_refactor`'s
    style.
- `doc/pypeit_par.rst` and all `doc/include/parset_*.rst` files (plus
  possibly `doc/scripts/base_par.rst`) are **auto-generated** by
  `doc/scripts/build_par_rst.py` — per user direction, **any conflict in an
  `.rst` file gets flagged to the user before attempting a resolution**,
  rather than hand-merged, since regenerating from the finished
  `pypeitpar.py` is almost certainly correct and cheaper. Conflicts in the
  generator scripts themselves (`doc/scripts/*.py`, `pypeit/scripts/*.py`)
  are fixed normally like any other code conflict. The user will run
  `cd doc; make clean; make html` themselves once the whole rebase is done —
  that step is not part of this plan's validation.
- `develop`'s semantic changes to `pypeitpar.py` that must be ported into the
  new declarative style (checklist):
  - `ProcessImagesPar`: new `cr_median_width` (int, default 0)
  - `FlatFieldPar`: new `fiber_pixelflat` (bool, default False)
  - `AlignPar`: new `grow_slit_edge` (default 0.0)
  - `ScatteredLightPar`: new valid method option `'gaps'`
  - `CubePar`: `align` (bool) removed, replaced by `alignment_method` (str,
    default `'phase'`, options `none/user/phase/cc/fit`); new `save_native`,
    `save_separate`, `extraction` (nested, new `CubeExtractionPar` class),
    `weights_init_obj_pos`, `sn_smooth_npix`
  - New class `CubeExtractionPar` (output_filename, whitelight_range, fwhm,
    snr_thresh, manual, boxcar_radius, opt_prof_method, skysub_resid)
  - `SensFuncPar`: new `star_arxiv` (default `'default'`)
  - `WavelengthSolutionPar`: new `lamps_wvrng`, `ech_angle_fits_file`,
    `ech_composite_arc_file`, `ech_direct_cc`
  - `FindObjPar`: new `force_center_obj` (bool, default False)
  - `SkySubPar`: new `joint_fit_use_sci` (bool, default True),
    `sci_exclude_radius`
  - `TelescopePar`: `valid_telescopes` gains `'APO'`, `'INT'`
  - `ReducePar.from_dict`: fix so `trim_edge` is actually parsed from cfg
  - `EdgeTracePar`: description-text-only clarification (non-semantic)
- **Systemic risk beyond git conflicts**: `develop`'s 851 commits were
  written against the *old* `ParSet` API. Any code they added in files
  `parset_refactor` never touched could still use the removed
  `.descr[key]`/`.options[key]`/`.dtype[key]`/`can_call` attribute style,
  which will fail silently (no merge conflict, but breaks at runtime). Must
  be swept for after the rebase, not caught by git.
- The user already has purpose-built tooling for exactly this
  verification, at the repo root / in `pypeit/tests/`:
  - `pypeit_par_accounting.py` — dumps every discoverable `ParSet` subclass's
    keys/dtype/default/options/descr to a text file, for both the old-style
    (`old_parsets`) and new-style (`new_parsets`) API, enabling a line-level
    `diff` between them.
  - `pypeit/tests/test_pypeitpar.py::diff_pars` — recursively diffs two
    instantiated `ParSet` objects value-by-value.
  These should be rerun/reused as the completeness check, not reinvented.

## Plan

1. **Safety net**: `git checkout -b parset_refactor-rebase parset_refactor`
   and do all work on that new branch. `parset_refactor` itself is never
   touched, so it remains the fallback if anything goes wrong; no tag is
   needed. The user will delete `parset_refactor` once satisfied with the
   rebase.

2. **Squash to logical commits** on the new branch via
   `git reset --soft a838e303e` + re-commit in a small number of chunks
   (exact boundaries to be finalized when we do it, but roughly):
   - `ParSet` base class refactor (`pypeit/par/parset.py`, `par/__init__.py`)
   - `pypeitpar.py` restructuring + the small call-site updates it forces
     across the codebase (the ~24 trivial-overlap files plus the 63
     refactor-only files)
   - Auto-generated parameter docs (`doc/pypeit_par.rst`,
     `doc/include/parset_*.rst`, `doc/scripts/build_par_rst.py` and friends)

3. **Rebase onto develop**: `git rebase develop` from the squashed branch.
   Conflicts will concentrate in the pypeitpar.py commit (real conflict) and
   the docs commit. Any `.rst` conflict in that docs commit is surfaced to
   the user rather than resolved unilaterally (see above) — the expectation
   is we drop the conflicting `.rst` content and regenerate it later, but
   that's the user's call each time it comes up.

4. **Resolve `pypeitpar.py`**: don't trust the textual merge — rebuild each
   touched class starting from `develop`'s content, expressed in the new
   `parameters = {...}` / `parset.set_parameter_definition(...)` pattern,
   porting in every item from the checklist above (including the new
   `CubeExtractionPar` class and the `ReducePar.from_dict` fix).

5. **Resolve the three special files** (`qa.py`, `collate_1d.py`,
   `telescopes.py`) per the verdicts above, and mechanically fix the
   remaining ~24 trivial-overlap files (take `develop`, reapply the small
   `ParSet`-API syntax tweak).

6. **Docs**: fix any conflicts in `doc/scripts/*.py` generator code as
   ordinary code conflicts. For the generated `.rst` outputs themselves
   (`doc/pypeit_par.rst`, `doc/include/parset_*.rst`), flag conflicts to the
   user rather than hand-resolving — regeneration via
   `doc/scripts/build_par_rst.py` against the finished `pypeitpar.py` is the
   likely fix, but we confirm before discarding either side's content.

7. **Sweep for the systemic risk**: grep the full post-rebase tree (not just
   the known-overlap files) for old-style `ParSet` access patterns
   (`\.descr\[`, `\.options\[`, `\.dtype\[`, `can_call`, stray `from_dict`
   overrides) that `develop` may have introduced in files `parset_refactor`
   never touched.

8. **Validate**:
   - Rerun `pypeit_par_accounting.py` (`new_parsets`) against the rebased
     `pypeitpar.py`; diff against a fresh `old_parsets` dump taken from
     `develop`'s pre-rebase `pypeitpar.py`, to confirm full parameter parity
     (only *intentional* differences, like the `FrameGroupPar` subclass key
     differences already seen in `pypeit/tests/comparison.txt`, should
     remain).
   - Run the full test suite, especially `pypeit/tests/test_pypeitpar.py`,
     `test_inputfiles.py`, `test_spectrographs.py`.
   - Doc build (`make clean; make html`) is the user's own follow-up step
     once the whole rebase is complete — not part of this plan's validation.

9. **Final style-coherence pass**: re-read the whole of `pypeitpar.py` once
   the port is complete to make sure everything ported in from `develop`
   (new classes/parameters) matches the conventions the refactor already
   established uniformly elsewhere (kwarg ordering in
   `set_parameter_definition`, docstring `.. include::` directive, the
   `FrameGroupPar`-subclass `parameters = FrameGroupPar.parameters | {...}`
   merge idiom, etc.) — not a redesign, just closing any inconsistency the
   manual port introduces.

Once this plan is approved, we'll execute it step by step interactively in
the terminal rather than as one large unattended change, since the
`pypeitpar.py` reconciliation in particular needs a human check at each
class.

## Future work: rebasing the branches stacked on `parset_refactor`

Not part of this pass, but confirmed while investigating (via
`git merge-base parset_refactor <branch>` for every local branch) and worth
planning around, since squashing `parset_refactor`'s history changes the
commit SHAs every dependent branch currently shares with it:

- One lineage shares merge-base `d4dd3aa55` (a past
  "Merge branch 'develop' into parset_refactor" commit) with
  `parset_refactor`: `wavemask` (10 commits ahead) and `funcpar_update` (22)
  branch directly from it; `dmod_update` (42) merges `funcpar_update` +
  `wavemask`; `speclist` (53) builds on `dmod_update`; `latinhyp` (62) builds
  on `speclist`; `standard_update` (61) also branches from the same point.
  This matches the chain already sketched in the user's own `notes` scratch
  file (`parset_refactor → funcpar_update/wavemask → dmod_update → speclist
  → latinhyp`, plus `standard_update`).
- A second lineage shares a later merge-base (`04fbd998a`) with
  `parset_refactor`: `spec_dm` (42 commits, merges `funcpar_update`) and
  `sensfunc_merge` (211, merges `spec_dm`).
- Because these branches' shared history with `parset_refactor` will no
  longer exist verbatim on `parset_refactor-rebase` (it was squashed), a
  plain `git merge`/`git rebase parset_refactor-rebase` from any of them
  would try to replay all of `parset_refactor`'s original 78 commits again.
  The correct tool once we get there is
  `git rebase --onto parset_refactor-rebase <old-parset_refactor-tip> <branch>`
  (rebasing only the commits unique to each branch onto the new base),
  applied in dependency order (e.g. `wavemask` and `funcpar_update` first,
  then `dmod_update`, then `speclist`, then `latinhyp`/`standard_update`; and
  separately `spec_dm` then `sensfunc_merge`).
- Expect renewed `pypeitpar.py`/`parset.py`-adjacent conflicts in several of
  these: the notes file records `wavemask` touching `pypeit/par/pypeitpar.py`
  and `latinhyp` touching `pypeit/par/parset.py`, for example.
- This cascade should be scoped as its own follow-up plan once
  `parset_refactor-rebase` is verified and adopted as the new
  `parset_refactor`.

## Appendix: `PypeIt-development-suite` repo — merge `develop` into `parset_refactor`

This is a **separate, much smaller task** in the `PypeIt-development-suite`
repo (not the main `PypeIt` repo the plan above concerns), prompted by the
fact that `claude_prompts/` (where this file lives) only exists on that
repo's `develop` branch, not yet on its `parset_refactor` branch. Per user
direction, this repo gets a **merge**, not a rebase — there's no equivalent
refactor-collision risk here, so there's no need for the squash/rebase
machinery used above.

Assessment:

- Merge-base of `develop`/`parset_refactor` in this repo: `85a9d0ddc`.
  `develop` is 255 commits ahead; `parset_refactor` is 15 commits ahead —
  both far smaller than the main repo's divergence.
- Only 2 files are touched by both branches, and both diffs are tiny:
  - `unit_tests/test_setup_gui.py`: `develop` +3/-4, `parset_refactor` -1.
  - `vet_tests/test_wavelengths.py`: `develop` +30/-0, `parset_refactor`
    +1/-1.
- A read-only `git merge-tree develop parset_refactor` preview (no working
  tree changes) completed cleanly with no conflict output at all — the two
  small overlapping diffs land in different parts of their respective files.
- **Conclusion: a plain `git merge develop` on the `parset_refactor` branch
  of this repo should go through with zero conflicts.** No special handling
  needed; just run the merge and confirm the resulting `claude_prompts/`
  directory (and anything else new from `develop`) shows up as expected.
