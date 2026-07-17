# PypeIt Output Paths: Implementation Plan (Design A)

## 0. Scope

This document turns **Design A** — the centralized `PypeItOutputPaths`
singleton — into a concrete, execution-ready implementation plan. The
choice between Design A and Design B is settled (see `pypeit_output_paths.md`
§4.1 for the recommendation and rationale); this document does not revisit
that decision or repeat the violation catalog, the current-state analysis,
or Design B. Refer back to `pypeit_output_paths.md` for:

- The original task prompt (§1).
- The full current-state analysis and 9 catalogued path violations (§2).
- Design A's original class sketch, lifecycle, and violation-by-violation
  mapping (§3.1).
- The mutability-guards discussion (§4.2-§4.3): the team agreed to
  implement guards **1** (freeze-on-first-use), **2** (idempotent no-op
  reconfigure), and **5** (logging); to explicitly **skip guard 4**
  (filesystem sentinel file); and to **hold guard 3** (`derive()` for
  non-mutating derived paths) **in reserve** — designed for, but not
  implemented or wired into any caller in this pass.

All file:line citations below were verified directly against the current
source tree (package root `pypeit/pypeit/`) at plan-writing time.

## 1. Exception handling

No new exception class is needed. `PypeItPathError` already exists at
`pypeit/pkg/exceptions.py:21` (`class PypeItPathError(PypeItError)`), is
already in that module's `__all__`, and is already re-exported at package
level via `pypeit/__init__.py:29` (`from .pkg.exceptions import *`), so it's
importable as `from pypeit import PypeItPathError`. `pypeit/outputpaths.py`
will import it directly from `pypeit.pkg.exceptions` (not via the star
re-export) to avoid any import-ordering sensitivity during package init.

## 2. Concrete class design — `pypeit/outputpaths.py`

New module, sibling of `pypeit/outputfiles.py` (whose naming helpers are
reused, not duplicated):

```python
"""
Centralized, singleton-managed output-path resolution for PypeIt reductions.

A single instance, ``pypeit.outputPaths``, is created at import time with
cwd-based defaults (mirroring ``pypeit.log`` and ``pypeit.dataPaths``).
Orchestration-layer code (PypeIt, Calibrations, CoAdd2D, script main()s)
reconfigures it once from a PypeItPar via ``configure()`` and creates
directories via ``make()``. pypeit/core/ code MUST NOT import or query this
object; it only ever receives already-resolved Path/str arguments.
"""
from pathlib import Path

from pypeit import log
from pypeit.pkg.exceptions import PypeItPathError


class PypeItOutputPaths:

    # Attribute names that map to a resolved directory Path, usable with make().
    _MAKE_TARGETS = ('science', 'qa', 'qa_pngs', 'calibrations',
                      'coadd_science', 'coadd_qa', 'coadd_qa_pngs', 'collate')

    def __init__(self, redux_path=None, scidir='Science', qadir='QA',
                 calib_dir='Calibrations', *, coadd_suffix='_coadd',
                 collate_outdir=None):
        self._frozen = False
        self._apply(redux_path, scidir, qadir, calib_dir, coadd_suffix,
                    collate_outdir)

    # ---- internal state management -------------------------------------
    def _apply(self, redux_path, scidir, qadir, calib_dir, coadd_suffix,
               collate_outdir):
        """Set all resolved state in place (no freeze/idempotency checks)."""
        self.redux_path = Path(redux_path).absolute() if redux_path else Path.cwd()
        self._scidir = scidir
        self._qadir = qadir
        self._calib_dir = calib_dir
        self._coadd_suffix = coadd_suffix
        self._collate_outdir = Path(collate_outdir).absolute() if collate_outdir \
                                else None

    def _state(self):
        """Hashable tuple fully describing the resolved state (for guard 2)."""
        return (self.redux_path, self._scidir, self._qadir, self._calib_dir,
                self._coadd_suffix, self._collate_outdir)

    # ---- resolved-path properties --------------------------------------
    @property
    def science(self) -> Path: return self.redux_path / self._scidir
    @property
    def qa(self) -> Path: return self.redux_path / self._qadir
    @property
    def qa_pngs(self) -> Path: return self.qa / 'PNGs'
    @property
    def calibrations(self) -> Path: return self.redux_path / self._calib_dir
    @property
    def coadd_science(self) -> Path:
        return self.redux_path / f'{self._scidir}{self._coadd_suffix}'
    @property
    def coadd_qa(self) -> Path:
        return self.redux_path / f'{self._qadir}{self._coadd_suffix}'
    @property
    def coadd_qa_pngs(self) -> Path: return self.coadd_qa / 'PNGs'
    @property
    def collate(self) -> Path:
        return self._collate_outdir if self._collate_outdir is not None \
                    else self.redux_path
    @property
    def frozen(self) -> bool: return self._frozen

    # ---- reconfiguration: guards 1 (freeze), 2 (idempotent), 5 (log) ---
    def configure(self, par=None, *, redux_path=None, scidir=None, qadir=None,
                  calib_dir=None, coadd_suffix=None, collate_outdir=None,
                  allow_reconfigure=False, caller=None):
        """
        Reconfigure the singleton in place from a PypeItPar and/or explicit
        overrides. Explicit keyword overrides take precedence over ``par``.
        """
        if par is not None:
            rdx, cal = par['rdx'], par['calibrations']
            redux_path = redux_path if redux_path is not None else rdx['redux_path']
            scidir = scidir if scidir is not None else rdx['scidir']
            qadir = qadir if qadir is not None else rdx['qadir']
            calib_dir = calib_dir if calib_dir is not None else cal['calib_dir']
            if collate_outdir is None and 'collate1d' in par.keys() \
                    and par['collate1d']['outdir'] is not None:
                collate_outdir = par['collate1d']['outdir']

        new_redux = Path(redux_path).absolute() if redux_path else self.redux_path
        new = (new_redux,
               self._scidir if scidir is None else scidir,
               self._qadir if qadir is None else qadir,
               self._calib_dir if calib_dir is None else calib_dir,
               self._coadd_suffix if coadd_suffix is None else coadd_suffix,
               (Path(collate_outdir).absolute() if collate_outdir is not None
                    else self._collate_outdir))

        # Guard 2: unchanged resolved state -> silent no-op.
        if new == self._state():
            return

        # Guard 1: frozen (a directory has already been created) + would
        # actually change -> reject unless explicitly allowed.
        if self._frozen and not allow_reconfigure:
            raise PypeItPathError(
                'Output paths are frozen (directories already created) and '
                f'cannot be reconfigured. Current redux_path={self.redux_path}; '
                f'requested redux_path={new[0]}. Pass allow_reconfigure=True to '
                'override (only safe before any output has been written).')

        changed = [name for name, o, n in
                   zip(('redux_path', 'scidir', 'qadir', 'calib_dir',
                        'coadd_suffix', 'collate_outdir'), self._state(), new)
                   if o != n]
        self._apply(*new)

        # Guard 5: log every non-no-op reconfigure.
        ctx = f' [caller={caller}]' if caller else ''
        log.info(f'Output paths reconfigured{ctx}: redux_path={self.redux_path}; '
                 f'changed={changed}')

    # ---- directory creation: freezes the object; guard 5 logging -------
    def make(self, *which, caller=None):
        """
        Create the named output directories and freeze the object against
        further silent reconfiguration. Defaults to the three primary
        directories if none are named.
        """
        if not which:
            which = ('science', 'qa_pngs', 'calibrations')
        made = []
        for name in which:
            if name not in self._MAKE_TARGETS:
                raise PypeItPathError(
                    f'Unknown output-path target: {name!r}. '
                    f'Valid targets: {self._MAKE_TARGETS}.')
            p = getattr(self, name)
            p.mkdir(parents=True, exist_ok=True)
            made.append((name, str(p)))
        self._frozen = True
        ctx = f' [caller={caller}]' if caller else ''
        log.info(f'Created output directories{ctx}: '
                 + ', '.join(f'{n}={p}' for n, p in made))
        return [Path(p) for _, p in made]

    # ---- file-naming convenience (delegates to outputfiles.py) ---------
    def spec_output_file(self, fitstbl, frame, twod=False, ext='.fits',
                          coadd=False):
        """Delegates to outputfiles.spec_output_file using the resolved
        science path; imported locally to avoid a package-init import cycle."""
        from pypeit import outputfiles
        sci_path = self.coadd_science if coadd else self.science
        return outputfiles.spec_output_file(fitstbl, None, frame, twod=twod,
                                             ext=ext, sci_path=sci_path)

    # ---- guard 3: HELD IN RESERVE — stub only, not wired to any caller --
    def derive(self, **overrides):
        """
        Return an independent PypeItOutputPaths seeded from current state
        with the given overrides applied, so derived-path computations
        (coadd, collate, sensfunc) never mutate the shared singleton.

        Reserved for future work. Not called anywhere in this implementation
        pass; Phase 4 instead reads/writes the singleton's own
        coadd_*/collate properties directly (see §5, Phase 4).
        """
        raise NotImplementedError(
            'PypeItOutputPaths.derive() is reserved for future work and is '
            'not yet implemented.')

    def __repr__(self):
        return (f'<{self.__class__.__name__}: redux_path={self.redux_path}, '
                f'frozen={self._frozen}>')
```

Notes:
- `spec_output_file` always supplies `sci_path` explicitly, so it never
  touches `par` — `outputfiles.spec_output_file` only falls back to `par`
  when `sci_path is None`.
- Guard 2's idempotency check compares the fully-resolved `_state()` tuple
  (absolute `Path`s + strings), so defensive repeated `configure()` calls
  from multiple orchestrators (`PypeIt.__init__`, `run_pypeit`, `build_qa`)
  with the same effective values are silent and never trip the freeze guard.
- Guard 1 freezes only in `make()`, never in `configure()`/`__init__` —
  `configure()` must be able to run freely before any directory exists.
- Guard 4 (filesystem sentinel file) is intentionally absent — not stubbed,
  not mentioned as a TODO — per the explicit decision to skip it.

## 3. `pypeit/__init__.py` wiring

Insert immediately after the existing `dataPaths` instantiation
(`pypeit/__init__.py:24-25`):

```python
# Import and instantiate the data path parser
# NOTE: This *MUST* come after log and __version__ are defined above
from .pkg.pypeitdata import PypeItDataPaths
dataPaths = PypeItDataPaths()

# Import and instantiate the output path manager
from .outputpaths import PypeItOutputPaths
outputPaths = PypeItOutputPaths()          # cwd-based defaults at import time

# Import all the exceptions so that they can be directly imported (e.g., `from
# pypeit import PypeItError`) in all package imports.
from .pkg.exceptions import *
```

No import cycle: `outputpaths.py` does `from pypeit import log` (re-entrant,
but `log` is already bound at this point in `__init__.py`, exactly as
`outputfiles.py` already does elsewhere) and `from pypeit.pkg.exceptions
import PypeItPathError` (a standalone module with no dependency back on the
partially-initialized `pypeit` package). The `outputfiles` import used by
`spec_output_file()` is deliberately deferred to function-call time to
avoid `outputfiles.py`'s own top-level `from pypeit import log` happening
before `__init__.py` finishes.

## 4. Corrections/confirmations from source verification

A few facts worth flagging before implementation, since they change the
shape of specific phases below:

- **`qa.py`'s `set_qa_filename` bug is confirmed and slightly worse than a
  simple hardcode**: several `case` branches (e.g. lines 65, 76, 100) still
  actively return `'QA/PNGs/...'` while sibling branches return only
  `'PNGs/...'` (with the alternate form commented out immediately
  above/below — evidence of a partial, abandoned prior fix). Since callers
  join these onto an already-QA-rooted `out_dir`, the `'QA/PNGs'` branches
  currently produce a doubled `<redux>/QA/QA/PNGs/...` path on disk. Phase 3
  must normalize **every** branch to the `'PNGs/...'` form.
- **`find_objects.py` is only half-broken**: `get_findobj_qa_filename`
  (~line 702) already correctly builds
  `os.path.join(par['rdx']['redux_path'], par['rdx']['qadir'])`. The
  hardcoded `os.path.join(..., 'QA')` (~line 1252) lives inside a dead
  `if False:` block (confirmed at ~line 1250, guarded by a `# TODO ::` note
  that flexure QA was never finished). This is a dead-code cleanup, not a
  live-bug fix, unlike `extraction.py:584`, which is live and broken.
- **`Calibrations.__init__`** (`calibrations.py:144-181`) signature is
  `(self, fitstbl, par, spectrograph, caldir, calib_ID, frame, det,
  qadir=None, reuse_calibs=False, show=False, user_slits=None,
  chk_version=True)`. It computes `self.calib_dir`/`self.qa_path` and their
  `mkdir`s directly from the `caldir`/`qadir` arguments it's given — the
  `# TODO: This should only be defined in one place!` comment (line 177)
  sits directly above the `qa_png_path = self.qa_path / 'PNGs'` derivation.
  The signature itself does not need to change; only where its `caldir`/
  `qadir` arguments come from, and whether `Calibrations` itself still does
  the `mkdir`s once the orchestrator (`pypeit_steps.calib_one`) starts
  calling `outputPaths.make()` first.
- **`pypeit_steps.calib_one`** (`pypeit_steps.py:124-179`) is the actual
  call site that resolves `qa_path` (lines ~162-163:
  `os.path.join(par['rdx']['redux_path'], par['rdx']['qadir'])`) when none
  is passed in, then forwards it to `Calibrations.get_instance(...,
  qadir=qa_path, ...)`.
- **`collate.py`**: `par['collate1d']['outdir']` is read directly at
  multiple call sites (lines 219, 229, 693, 802, and via
  `copy_spec1d_to_outdir` at 585-601, which does its own `os.makedirs` at
  line 597); the spec1d→spec2d string-replace is at line 668
  (`filename.replace('spec1d', 'spec2d', 1)`). Note `collate.py` also has a
  *second*, independent `par['collate1d']['spec1d_outdir']` key (lines 767,
  777, 780) used only for copying spec1d files — this is a distinct concept
  from `outdir` and is out of scope for this refactor (it isn't one of the
  9 catalogued violations); it should simply be left alone.

## 5. Phased implementation plan

Guards 1/2/5 land in **Phase 0** with the class itself; they get exercised
against real orchestrator behavior in **Phase 1**.

### Phase 0 — Introduce the module (no behavior change elsewhere)
**Files:** new `pypeit/outputpaths.py`; `pypeit/__init__.py`.
- Add the class exactly as in §2. Wire the singleton per §3.
- Nothing else imports it yet.
- **Verify:** new `pypeit/tests/test_outputpaths.py` (§6). Run
  `pytest pypeit/tests/test_outputpaths.py`. Confirm
  `python -c "import pypeit; print(pypeit.outputPaths)"` succeeds with no
  import-cycle error.

### Phase 1 — Adopt in `PypeIt` and `Calibrations`
**Files:** `pypeit.py`, `calibrations.py`, `pypeit_steps.py`.
- `pypeit.py` `PypeIt.__init__`: after the existing `redux_path` override
  (lines 87-88, which writes `self.par['rdx']['redux_path'] = redux_path`),
  call `outputPaths.configure(self.par, redux_path=redux_path,
  caller='PypeIt.__init__')`. Replace line 119's
  `self.calibrations_path = Path(...) / ...` with
  `self.calibrations_path = outputPaths.calibrations`.
- Change the `science_path` property (lines 144-147) to
  `return outputPaths.science` and the `qa_path` property (lines 149-152)
  to `return outputPaths.qa` — **this fixes the existing `str`-vs-`Path`
  type inconsistency** (both now uniformly return `Path`). Downstream
  consumers (`qa.gen_qa_dir`, `qa.gen_mf_html`, `Calibrations(qadir=...)`)
  already coerce their input through `pathlib.Path(...)`, so this is
  backward-compatible.
- `build_qa` (lines 154-163): call `outputPaths.make('qa_pngs',
  caller='PypeIt.build_qa')` instead of relying on `qa.gen_qa_dir`'s own
  `mkdir`.
- `calibrations.py` `Calibrations.__init__` (144-181): **signature
  unchanged**. Remove the inline `mkdir` calls (lines 171-172, 180-181) and
  the `# TODO` comment (line 177) — directory creation for `calibrations`/
  `qa_pngs` now happens upstream, in the orchestrator, before
  `Calibrations` is constructed.
- `pypeit_steps.py` `calib_one` (124-179): replace the `qa_path` default
  (lines 162-163, `os.path.join`) with `outputPaths.qa`; immediately before
  instantiating `Calibrations.get_instance(...)`, call
  `outputPaths.make('calibrations', 'qa_pngs', caller='calib_one')`.
  Continue passing the already-resolved `calibrations_path`/`qadir=qa_path`
  down as plain arguments (preserves the `core/`-isolation constraint —
  `Calibrations` itself isn't `core/`, but this keeps the pattern
  consistent all the way down).
- **Verify:** `pytest pypeit/tests/test_calibrations.py
  pypeit/tests/test_runpypeit.py`; one MultiSlit dev-suite setup (e.g.
  `shane_kast_blue`) end-to-end, confirming `Science/`, `QA/PNGs/`,
  `Calibrations/` land in the right place and that repeated internal
  `configure()` calls (once from `PypeIt.__init__`, implicitly again if
  `build_qa`/`calib_one` re-derive the same par) stay silent per guard 2.

### Phase 2 — Consolidate `outputfiles.py` as the naming backend
**Files:** `outputfiles.py`.
- No signature changes to `science_path(par)` or `spec_output_file(...)`
  (both have existing callers). This phase is a consolidation checkpoint,
  not a behavior change: confirm the new class's `spec_output_file()`
  convenience method (§2) delegates here rather than reimplementing naming
  logic, and add a regression assertion that the two never diverge.
- **Verify:** new unit test asserting `outputPaths.science ==
  outputfiles.science_path(par)` for a `configure()`-d par.

### Phase 3 — Fix `qa.py`, `extraction.py`, `find_objects.py`, `core/findobj_skymask.py`
**Files:** `qa.py`, `extraction.py`, `find_objects.py`,
`core/findobj_skymask.py`.
- `qa.set_qa_filename` (27-139): normalize **every** branch to
  `'PNGs/...'` (strip the residual active `'QA/PNGs/...'` occurrences,
  e.g. lines 65, 76, 100 — see §4). Delete the `# TODO` (lines 25-26).
  Signature/argument contract (`out_dir`) is unchanged — this remains a
  pure naming helper callable from any layer, no `outputPaths` import
  needed here.
- `qa.gen_exp_html` (494-520): add a required `qa_path` argument; replace
  the hardcoded `Path("QA")` (~line 500) and `f"QA/{uni_name}.html"`
  (~line 506) with paths built under the passed `qa_path`. Update its sole
  caller, `PypeIt.build_qa`, to pass `self.qa_path`.
- `extraction.py:584`: replace `os.path.join(self.par['rdx']['redux_path'],
  'QA')` with `os.path.join(self.par['rdx']['redux_path'],
  self.par['rdx']['qadir'])` — matching the pattern already correct in
  `find_objects.py:702`.
- `find_objects.py`: fix or delete the dead `if False:` block containing
  the hardcoded `'QA'` (~lines 1250-1254) — cosmetic/dead-code cleanup,
  not a live-bug fix (see §4).
- `core/findobj_skymask.py` (lines 1453-1455): remove the
  `qafile.parent.mkdir(parents=True, exist_ok=True)` side effect from
  `core/` entirely. The caller (`find_objects.py`, orchestration layer)
  already ensures the QA directory exists via Phase 1's
  `outputPaths.make('qa_pngs')`, so this becomes redundant, and its removal
  satisfies the "core/ never manages directories" principle exactly.
- **Verify:** extend `tests/test_qa.py` with a regression test asserting
  no `set_qa_filename` branch yields a doubled `QA/QA` segment and every
  returned path contains exactly one `PNGs` component; dev-suite run of a
  QA-plot-generating setup (flexure + object-finding QA) confirming no
  crash from the removed `core/` mkdir and correct `QA/PNGs/` placement.

### Phase 4 — Fix `coadd2d.py` and its callers
**Files:** `coadd2d.py`, `scripts/coadd_2dspec.py`, `scripts/ql.py`.
- `CoAdd2D.output_paths` (475-516): delete the `par['rdx']['qadir'] +=
  '_coadd'` mutation (line 512) entirely — no replacement mutation of
  any kind. Replace the hardcoded `'Science'` (505), the
  `f"{scidir}_coadd"` construction (508), and the hardcoded `'PNGs'` (513)
  with reads of `outputPaths.coadd_science` / `outputPaths.coadd_qa_pngs`
  after the caller has done `outputPaths.configure(par,
  redux_path=coadd_dir, caller='CoAdd2D.output_paths')`. The staticmethod
  keeps its existing call signature for API stability but its body becomes
  a thin wrapper over the singleton's properties plus
  `outputPaths.make('coadd_science', 'coadd_qa_pngs')`.
- `scripts/coadd_2dspec.py`: update the `output_paths(...)` call (84-85) to
  the non-mutating form; replace the hand-rolled `spec1d_{basename}.fits`/
  `.txt` (229, 233) and `spec2d_{basename}.fits` (238) with a shared
  coadd-naming helper in `outputfiles.py` (new, small function reused by
  both this script and `ql.py`, so the naming convention has exactly one
  implementation).
- `scripts/ql.py`: its `output_paths(spec2d_files, par)` call (~1022,
  called *without* `coadd_dir` and therefore also mutating `par` today)
  switches to the same non-mutating form; its hand-rolled `spec2d_
  {basename}.fits` (~1024) is replaced by the same shared helper used in
  `coadd_2dspec.py` (its sibling call at ~1028 already correctly uses
  `outputfiles.spec_output_file`, confirming the target pattern).
- **Verify:** unit test asserting `par['rdx']['qadir']` is bit-for-bit
  unchanged after computing coadd paths (regression guard for the removed
  mutation); dev-suite runs of a `coadd_2dspec` setup and a quicklook (`ql`)
  setup confirming `Science_coadd/` and `QA_coadd/PNGs/` are created and
  spec1d/spec2d names match the shared helper in both scripts.

### Phase 5 — Remaining scripts and `collate`
**Files:** `scripts/trace_edges.py`, `collate.py`, `scripts/collate_1d.py`,
`scripts/sensfunc.py`, `scripts/flux_calib.py`, `scripts/coadd_1dspec.py`,
`scripts/flux_setup.py`.
- `trace_edges.py`: remove the duplicated argparse literal
  `default='Calibrations'` (line 42); source the default from
  `CalibrationsPar`/`outputPaths` instead. Unify the two QA-path branches
  (line 127, correct; line 154, hardcoded `redux_path / 'QA'`) so both
  route through `outputPaths.configure(redux_path=redux_path)` then
  `outputPaths.qa`. Update `calib_dir` construction (line 163) to
  `outputPaths.calibrations`.
- `collate.py`/`collate_1d.py`: call `outputPaths.configure(par,
  collate_outdir=par['collate1d']['outdir'], caller='collate_1d')`, then
  use `outputPaths.collate` at the existing `outdir`-consuming call sites
  (lines 219, 229, 693, 802). Leave `par['collate1d']['spec1d_outdir']`
  (lines 767, 777, 780) untouched — it's a distinct, already-working
  concept, out of scope here. Replace the string-replace spec1d→spec2d
  translation (line 668) with the same shared `outputfiles.py` naming
  helper introduced in Phase 4. `collate_1d.py`'s `os.makedirs` (196-197)
  becomes `outputPaths.make('collate')`.
- `sensfunc.py`/`flux_calib.py`/`coadd_1dspec.py`/`flux_setup.py`: route
  `--par_outfile` defaults and product files (`sens_*`, `.flux`/
  `.coadd1d`/`.tell`, `coadd1d_*`) to be `redux_path`-relative via
  `outputPaths` wherever a redux context exists, defaulting to cwd
  (`outputPaths`'s own default) for genuinely standalone invocations — no
  behavior change for today's cwd-relative usage, just an explicit,
  overridable path instead of an implicit one. Fix the `flux_calib.py:62`
  `action="store_true"` combined with a string `default='fluxing.par'`
  (a pre-existing, unrelated argparse bug — the flag can never actually
  take a filename value) as a drive-by while touching this line.
- **Verify:** dev-suite regression for each of `collate_1d`, `sensfunc`,
  `flux_calib`, `coadd_1dspec`, `flux_setup`, `trace_edges` post-processing
  setups.

### Phase 6 — Cleanup, pathlib migration, full regression
**Files:** files touched in Phases 1-5, plus `metadata.py`, `collate.py`
(remaining `os.path` calls beyond those already touched), `pypeit_steps.py`
stragglers.
- Convert any remaining `os.path.join`/`os.makedirs`/`os.getcwd` calls in
  files already edited by Phases 1-5 to `pathlib`, per this project's
  stated preference for `pathlib` in new/touched code. **Do not** perform a
  blanket, repo-wide `os.path`→`pathlib` sweep beyond files this refactor
  already touches — that's separate cleanup work, not part of this plan.
- Remove now-dead helpers if fully superseded (e.g. `qa.gen_qa_dir` if
  `outputPaths.make('qa_pngs')` has fully replaced its role).
- **Verify:** full `pytest pypeit/tests/`, plus a representative
  multi-pypeline dev-suite subset (see §7).

## 6. Test plan

Convention (per `pypeit/tests/`): plain module-level `def test_*()`
functions, no `TestX` classes, fixed RNG seeds, unit-only (no large data
files). Relevant existing files: `tests/test_pypeitpar.py`,
`tests/test_calibrations.py`, `tests/test_runpypeit.py`, `tests/test_qa.py`,
`tests/test_coadd.py`.

**New file `pypeit/tests/test_outputpaths.py`:**
- `test_default_resolution()` — with `redux_path=None`, confirm
  `science == Path.cwd()/'Science'`, `qa_pngs == .../ 'QA'/'PNGs'`,
  `calibrations == .../'Calibrations'`, `collate == redux_path`,
  `coadd_science == .../'Science_coadd'`, `coadd_qa_pngs ==
  .../'QA_coadd'/'PNGs'`.
- `test_explicit_dirs(tmp_path)` — non-default `scidir`/`qadir`/
  `calib_dir`/`coadd_suffix`/`collate_outdir` all resolve correctly;
  `collate` returns the explicit `collate_outdir` when set.
- `test_configure_idempotent_noop(tmp_path, caplog)` — configure twice
  with identical resolved values; second call logs nothing and does not
  raise (guard 2).
- `test_configure_changes_logs(caplog)` — configure with a changed
  `redux_path`; assert state updates and exactly one `log.info` fires
  listing `redux_path` in `changed` (guard 5).
- `test_freeze_then_reject_reconfigure(tmp_path)` — `make('science')` then
  `configure(redux_path=other)` raises `PypeItPathError` (guard 1); assert
  `.frozen is True`.
- `test_freeze_then_allow_with_override(tmp_path)` — after freeze,
  `configure(redux_path=other, allow_reconfigure=True)` succeeds.
- `test_freeze_noop_when_unchanged(tmp_path)` — after freeze,
  `configure()` with identical values is silently accepted (does not
  raise, despite being frozen) because it's a no-op per guard 2.
- `test_make_creates_and_freezes(tmp_path)` — `make('science', 'qa_pngs',
  'calibrations')` creates all three directories on disk and sets
  `.frozen`.
- `test_make_unknown_target_raises()` — `make('bogus')` raises
  `PypeItPathError`.
- `test_configure_from_par()` — build a default `PypeItPar`, `configure`
  from it, assert properties match `outputfiles.science_path(par)`.
- `test_derive_not_implemented()` — `derive()` raises `NotImplementedError`
  (documents guard 3's held-in-reserve status).
- `@pytest.mark.skip(reason="derive() held in reserve")` `def
  test_derive_independent_copy(): ...` — a skipped placeholder recording
  the intended coverage once guard 3 is implemented (a derived object with
  overridden `redux_path`/`coadd_suffix` must not mutate the parent
  singleton), so the test intent is captured now without being exercised.

**Extensions to existing files:**
- `tests/test_qa.py`: `test_set_qa_filename_no_double_qa()` (every branch
  returns exactly one `PNGs` segment, never `QA/QA`).
- `tests/test_coadd.py`: `test_coadd2d_output_paths_no_par_mutation()` —
  snapshot `par['rdx']['qadir']` before/after computing coadd paths;
  assert unchanged.
- `tests/test_pypeitpar.py`: `test_output_path_defaults_unchanged()` —
  assert `ReduxPar.scidir=='Science'`, `.qadir=='QA'`,
  `.redux_path==os.getcwd()`, `CalibrationsPar.calib_dir=='Calibrations'`,
  `Collate1DPar.outdir==os.getcwd()` — guards the backward-compatibility
  invariant that no ParSet key is renamed/redefaulted by this refactor.
- `tests/test_runpypeit.py`: assert `PypeIt(...).science_path` and
  `.qa_path` are both `pathlib.Path` (regression test for the fixed
  str/Path inconsistency).

## 7. Rollout / PR sequencing

Six PRs, aligned to the six phases — **not** one big-bang PR:

- **PR 1 = Phase 0.** Zero behavior change anywhere else in the codebase;
  fully verifiable with `pypeit/tests` alone. Lowest risk; unblocks
  everything downstream. Merge first.
- **PR 2 = Phases 1-2.** The behavior-sensitive core (`PypeIt`,
  `Calibrations`, `pypeit_steps`, `outputfiles` consolidation). Requires
  dev-suite regression on at least one MultiSlit + one Echelle + one
  SlicerIFU setup (e.g. `shane_kast_blue`, a `keck_hires`/`keck_nires`
  setup, `keck_kcwi`) confirming `Science/`, `QA/PNGs/`, `Calibrations/`
  placement and that `-r/--redux_path` still overrides correctly. Also
  verify against an existing `.pypeit` file with custom `[rdx]
  scidir/qadir` and `[calibrations] calib_dir` values, to confirm the
  backward-compatibility invariant.
- **PR 3 = Phase 3.** QA-path normalization + `core/` mkdir removal.
  Unit-testable for `set_qa_filename`; dev-suite for the QA-generating
  flexure/objfind flow.
- **PR 4 = Phase 4.** `coadd2d`/`coadd_2dspec`/`ql` — the par-mutation
  removal is the headline correctness fix here. Dev-suite: coadd_2dspec
  and quicklook setups.
- **PR 5 = Phase 5.** `trace_edges`, `collate`, `sensfunc`, `flux_calib`,
  `coadd_1dspec`, `flux_setup`. Dev-suite: each post-processing tool's
  existing test setup.
- **PR 6 = Phase 6.** Scoped os.path→pathlib cleanup, dead-helper removal,
  broadest dev-suite regression subset across pypelines.

**Unit tests alone are sufficient for:** all of PR 1; the ParSet-defaults
backward-compatibility guard; `set_qa_filename` normalization; the
par-non-mutation assertion; `configure`/`make`/freeze/idempotency logic;
the `science_path`/`qa_path` typing fix.

**Dev-suite (real instrument data) verification is required starting at
PR 2** for every phase that changes where a file or directory is physically
written — the failure mode (files landing in the wrong absolute directory,
or a frozen-path exception surfacing under real repeated orchestrator
calls) is invisible to unit tests that don't run a full reduction.

Each PR should update the changelog / in-development release notes per
this project's standard changelog convention, noting explicitly that no
user-facing `PypeItPar` key is renamed, removed, or redefaulted by this
work.
