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

### 0.1 Guidelines for where paths may be overridden

These are **guidelines, not hard requirements** (except the first, which
already was one) — deviations are acceptable case-by-case, but only where a
deviation is the least-bad option, and any such deviation should be called
out explicitly, not silently introduced. Every phase below is written to
satisfy them; where a phase doesn't, that is flagged and justified.

1. **(Already in place.)** Paths and output files required by the
   low-level functions in `pypeit/core/` *must* be passed to those
   functions as positional and/or keyword arguments. `core/` code never
   imports or queries `outputPaths`.
2. **Paths required by high-level classes (e.g. `Calibrations`) should use
   `PypeItOutputPaths` directly, and any positional or keyword argument
   that would let a caller override those paths should be removed.** A
   class that both accepts a path argument *and* has access to
   `outputPaths` has two disparate sources of truth that can silently
   diverge — this is exactly the failure mode the whole effort is meant to
   eliminate, and it applies just as much to an argument like
   `Calibrations.__init__`'s `caldir` as it does to `qadir` (see §4).
3. **Configuring `PypeItOutputPaths` should happen once, as soon as the
   `PypeItPar` object is defined — preferably inside the relevant script
   class in `scripts/`, before any high-level class that consumes paths is
   constructed — and then remain unchanged for the rest of the workflow.**
   In practice this means: script `main()` methods call
   `outputPaths.configure(par, ...)` exactly once, immediately after `par`
   is built/parsed; every class or function downstream of that point reads
   `outputPaths` but never calls `configure()` again. This is enforced,
   not just a convention: a second `configure()` call raises (§2) unless
   the instance is in dry-run mode (used only for inspection/testing,
   never in production code paths).

## 1. Exception handling

No new exception class is needed. `PypeItPathError` already exists at
`pypeit/pkg/exceptions.py:21` (`class PypeItPathError(PypeItError)`), is
already in that module's `__all__`, and is already re-exported at package
level via `pypeit/__init__.py:29` (`from .pkg.exceptions import *`), so it's
importable as `from pypeit import PypeItPathError`. `pypeit/pkg/outputpaths.py`
will import it directly via a relative sibling import (`from .exceptions
import PypeItPathError`), matching the existing precedent in
`pypeit/pkg/pypeitdata.py:66`.

## 2. Concrete class design — `pypeit/pkg/outputpaths.py`

New module inside `pypeit/pkg/` — alongside `pypeitdata.py`, `logger.py`,
`exceptions.py`, and `cache.py` — rather than at the top level of the
`pypeit` package, matching where its two closest architectural precedents
(`PypeItDataPaths`, `PypeItLogger`) already live. All imports of other
`pypeit` modules use relative imports, consistent with how `pypeitdata.py`
and `cache.py` already import their `pkg`-sibling modules (e.g.
`pypeitdata.py:66-67`: `from .exceptions import ...` / `from . import cache`).

```python
"""
Centralized, singleton-managed output-path resolution for PypeIt reductions.

A single instance, :obj:`pypeit.outputPaths`, is created at import time with
cwd-based defaults (mirroring :obj:`pypeit.log` and :obj:`pypeit.dataPaths`),
but is not yet *configured*. Orchestration-layer code (a ``scripts/`` script
class, or :class:`~pypeit.pypeit.PypeIt`'s ``__init__`` for entry points
where a script can't build the parameter set itself, e.g. ``run_pypeit``)
calls :func:`PypeItOutputPaths.configure` exactly once per CLI execution to
lock in the real, parameter-derived values. Code in ``pypeit/core/`` must
NOT import or query this object; it only ever receives already-resolved
``Path``/``str`` arguments.

The once-only rule has exactly one narrow, documented exception (the
``force`` parameter of :func:`configure`), used solely by
:class:`~pypeit.pypeit.PypeIt`'s ``__init__`` so that more than one
``PypeIt`` object can be constructed in the same process. See that
parameter's docstring before assuming this is a general pattern -- it is
not; every other caller must rely on the default once-only behavior.

Each managed directory (``redux``, ``science``, ``qa``, ``qa_pngs``,
``calibrations``, ``coadd_science``, ``coadd_qa``, ``coadd_qa_pngs``,
``collate``) is tracked internally as a small :class:`_ManagedPath` record
holding ``parent``, ``name``, ``full``, and ``ready``. :func:`configure`
(and ``__init__``) set only ``parent``/``name`` for every record -- pure
path arithmetic, no I/O, no change to ``ready``. The corresponding
``@property`` (e.g. :attr:`~PypeItOutputPaths.science`) is what resolves and
caches ``full`` (on first access, if not already cached) and creates the
directory on disk exactly once, flipping that record's ``ready`` flag --
but only if the instance has been configured and is not in dry-run mode. In
dry-run mode (or before :func:`configure` has ever been called), a property
still resolves and returns ``full``; it just never touches the filesystem
or sets ``ready``.
"""
import dataclasses
import logging
from pathlib import Path

from .exceptions import PypeItPathError

# Retrieves the already-configured "pypeit" logger by name from Python's
# logging registry -- does NOT recreate or reconfigure it (that already
# happened once, in pypeit/__init__.py, before this module is imported).
# A relative/absolute `from pypeit import log` is deliberately avoided here:
# pkg/ modules do not import back from the partially-initialized top-level
# pypeit package (see the circular-import warning atop pkg/logger.py).
log = logging.getLogger('pypeit')


@dataclasses.dataclass
class _ManagedPath:
    """Internal record for one managed output directory."""
    parent: Path
    """Parent directory."""
    name: str
    """Directory name, relative to `parent`."""
    full: Path = None
    """Full, absolute path (`parent`/`name`); resolved and cached on first
    access."""
    ready: bool = False
    """Whether the directory has been created on disk (or was already
    found to exist) and is available to receive files."""


class PypeItOutputPaths:
    """
    Single source of truth for PypeIt's output directory structure.

    See the module docstring for the configure-once / lazy-create-on-access
    lifecycle. Parameters mirror the ``PypeItPar`` keys that define the same
    concepts (:class:`~pypeit.par.pypeitpar.ReduxPar`,
    :class:`~pypeit.par.pypeitpar.CalibrationsPar`).

    Parameters
    ----------
    redux_path : :obj:`str`, :obj:`~pathlib.Path`, optional
        Root reduction directory. Defaults to the current working
        directory.
    scidir : :obj:`str`, optional
        Name of the science-output subdirectory.
    qadir : :obj:`str`, optional
        Name of the quality-assessment subdirectory.
    calib_dir : :obj:`str`, optional
        Name of the processed-calibrations subdirectory.
    coadd_suffix : :obj:`str`, optional
        Suffix appended to ``scidir``/``qadir`` for the 2D-coadd output
        directories.
    collate_outdir : :obj:`str`, :obj:`~pathlib.Path`, optional
        Output directory for :ref:`pypeit_collate_1d`. Defaults to
        ``redux_path`` if not given.
    dryrun : :obj:`bool`, optional
        If True, properties resolve ``full`` but never create directories
        or set ``ready``, and :func:`configure` may be called more than
        once. Intended for inspection/testing only; always ``False`` in
        production code paths.
    """

    def __init__(self, redux_path=None, scidir='Science', qadir='QA',
                 calib_dir='Calibrations', coadd_suffix='_coadd',
                 collate_outdir=None, dryrun=False):
        self._configured = False
        self._dryrun = dryrun
        self._paths = {}
        self._apply(redux_path, scidir, qadir, calib_dir, coadd_suffix,
                    collate_outdir)

    # ---- internal state management -------------------------------------
    def _apply(self, redux_path, scidir, qadir, calib_dir, coadd_suffix,
               collate_outdir):
        """
        Compute and store ``parent``/``name`` for every managed directory,
        replacing any existing records (so ``full``/``ready`` are always
        reset to unresolved/not-ready here). Pure path arithmetic -- no I/O,
        no change to ``_configured``.
        """
        redux_full = Path(redux_path).absolute() if redux_path is not None \
                        else Path.cwd()
        qa_full = redux_full / qadir
        coadd_qa_name = f'{qadir}{coadd_suffix}'
        coadd_qa_full = redux_full / coadd_qa_name
        collate_full = Path(collate_outdir).absolute() if collate_outdir is not None \
                            else redux_full

        self._paths = {
            'redux':         _ManagedPath(redux_full.parent, redux_full.name),
            'science':       _ManagedPath(redux_full, scidir),
            'qa':            _ManagedPath(redux_full, qadir),
            'qa_pngs':       _ManagedPath(qa_full, 'PNGs'),
            'calibrations':  _ManagedPath(redux_full, calib_dir),
            'coadd_science': _ManagedPath(redux_full, f'{scidir}{coadd_suffix}'),
            'coadd_qa':      _ManagedPath(redux_full, coadd_qa_name),
            'coadd_qa_pngs': _ManagedPath(coadd_qa_full, 'PNGs'),
            'collate':       _ManagedPath(collate_full.parent, collate_full.name),
        }

    def _get(self, key: str) -> Path:
        """Resolve (and, if eligible, create) one managed directory."""
        rec = self._paths[key]
        if rec.full is None:
            rec.full = rec.parent / rec.name
        if not rec.ready and self._configured and not self._dryrun:
            rec.full.mkdir(parents=True, exist_ok=True)
            rec.ready = True
            log.info(f'Output directory ready ({key}): {rec.full}')
        return rec.full

    # ---- resolved-path properties ---------------------------------------
    @property
    def redux(self) -> Path:
        """Root reduction directory."""
        return self._get('redux')

    @property
    def science(self) -> Path:
        """Science-output directory."""
        return self._get('science')

    @property
    def qa(self) -> Path:
        """Quality-assessment directory."""
        return self._get('qa')

    @property
    def qa_pngs(self) -> Path:
        """QA PNG-plot directory."""
        return self._get('qa_pngs')

    @property
    def calibrations(self) -> Path:
        """Processed-calibrations directory."""
        return self._get('calibrations')

    @property
    def coadd_science(self) -> Path:
        """2D-coadd science-output directory."""
        return self._get('coadd_science')

    @property
    def coadd_qa(self) -> Path:
        """2D-coadd quality-assessment directory."""
        return self._get('coadd_qa')

    @property
    def coadd_qa_pngs(self) -> Path:
        """2D-coadd QA PNG-plot directory."""
        return self._get('coadd_qa_pngs')

    @property
    def collate(self) -> Path:
        """Output directory for :ref:`pypeit_collate_1d`."""
        return self._get('collate')

    # ---- object-level state (read-only) ---------------------------------
    @property
    def configured(self) -> bool:
        """Whether :func:`configure` has been called on this instance."""
        return self._configured

    @property
    def dryrun(self) -> bool:
        """Whether this instance is in dry-run (inspection-only) mode."""
        return self._dryrun

    # ---- one-time (re)configuration --------------------------------------
    def configure(self, par=None, redux_path=None, scidir=None, qadir=None,
                  calib_dir=None, coadd_suffix=None, collate_outdir=None,
                  dryrun=None, caller=None, force=False):
        """
        Configure the instance from a :class:`~pypeit.par.pypeitpar.PypeItPar`
        and/or explicit overrides, exactly as if it had just been
        constructed with these arguments -- every managed directory's
        ``full``/``ready`` state is discarded and recomputed from scratch.
        May only be called once per instance unless the instance is in
        dry-run mode (or ``force`` is used, see below).

        Parameters
        ----------
        par : :class:`~pypeit.par.pypeitpar.PypeItPar`, optional
            Full parameter set to draw ``redux_path``/``scidir``/
            ``qadir``/``calib_dir``/``collate_outdir`` from. Explicit
            keyword arguments take precedence over the values in ``par``.
        redux_path : :obj:`str`, :obj:`~pathlib.Path`, optional
            Explicit override for the reduction root.
        scidir : :obj:`str`, optional
            Explicit override for the science-output subdirectory name.
        qadir : :obj:`str`, optional
            Explicit override for the QA subdirectory name.
        calib_dir : :obj:`str`, optional
            Explicit override for the calibrations subdirectory name.
        coadd_suffix : :obj:`str`, optional
            Explicit override for the 2D-coadd directory suffix.
        collate_outdir : :obj:`str`, :obj:`~pathlib.Path`, optional
            Explicit override for the collate1d output directory.
        dryrun : :obj:`bool`, optional
            If not None, set the instance's dry-run mode.
        caller : :obj:`str`, optional
            Free-form string identifying the calling code, included in the
            log message for easier debugging.
        force : :obj:`bool`, optional
            Bypass the once-only restriction and reconfigure even if this
            instance has already been configured (and is not in dry-run
            mode). **This is not a general-purpose escape hatch** -- it
            exists solely so that :class:`~pypeit.pypeit.PypeIt` can be
            instantiated more than once in the same process (e.g.,
            :mod:`~pypeit.scripts.ql` builds a calibration ``PypeIt`` and a
            science ``PypeIt`` in one run; a notebook or driver script may
            loop over several reductions). Constructing a new ``PypeIt``
            object is a deliberate signal that a new, independent reduction
            context is starting, which is categorically different from
            some other function accidentally calling :func:`configure`
            mid-run. Do not pass ``force=True`` from anywhere other than
            :class:`~pypeit.pypeit.PypeIt`'s ``__init__``; every other
            caller should rely on the default once-only behavior.

        Raises
        ------
        PypeItPathError
            If the instance has already been configured and is not in
            dry-run mode, and ``force`` is not used.
        """
        # See the `force` docstring above: this bypass is a narrow carve-out
        # for PypeIt.__init__ only, not a general reconfiguration mechanism.
        if self._configured and not self._dryrun and not force:
            raise PypeItPathError(
                'PypeItOutputPaths has already been configured for this '
                'execution and cannot be reconfigured (only permitted when '
                'the instance is in dry-run mode, for inspection/testing, '
                'or via the narrow force=True carve-out used by '
                'PypeIt.__init__).')

        if par is not None:
            rdx, cal = par['rdx'], par['calibrations']
            redux_path = redux_path if redux_path is not None else rdx['redux_path']
            scidir = scidir if scidir is not None else rdx['scidir']
            qadir = qadir if qadir is not None else rdx['qadir']
            calib_dir = calib_dir if calib_dir is not None else cal['calib_dir']
            if collate_outdir is None and 'collate1d' in par.keys() \
                    and par['collate1d']['outdir'] is not None:
                collate_outdir = par['collate1d']['outdir']

        if dryrun is not None:
            self._dryrun = dryrun

        self._apply(redux_path,
                    scidir if scidir is not None else 'Science',
                    qadir if qadir is not None else 'QA',
                    calib_dir if calib_dir is not None else 'Calibrations',
                    coadd_suffix if coadd_suffix is not None else '_coadd',
                    collate_outdir)
        self._configured = True

        ctx = f' [caller={caller}]' if caller else ''
        dry = ' (dry run)' if self._dryrun else ''
        log.info(f'Output paths configured{ctx}{dry}: '
                 f'redux_path={self._paths["redux"].parent / self._paths["redux"].name}')

    # ---- held in reserve — stub only, not wired to any caller ------------
    def derive(self, **overrides):
        """
        Return an independent :class:`PypeItOutputPaths` seeded from current
        state with the given overrides applied.

        Reserved for future work; not called anywhere in the current
        implementation. Its original motivation (letting coadd/collate
        compute derived paths without touching the shared singleton) is
        largely superseded by treating `coadd_science`/`coadd_qa`/`collate`
        as first-class properties on this same object -- kept as a stub
        rather than removed, pending a decision on whether it's still
        needed for some other future use.

        Raises
        ------
        NotImplementedError
            Always -- this method is not yet implemented.
        """
        raise NotImplementedError(
            'PypeItOutputPaths.derive() is reserved for future work and is '
            'not yet implemented.')

    def __repr__(self):
        redux = self._paths['redux']
        return (f'<{self.__class__.__name__}: '
                f'redux_path={redux.parent / redux.name}, '
                f'configured={self._configured}, dryrun={self._dryrun}>')
```

Notes:
- No file-naming convenience method is included on the class itself.
  Call sites use `outputPaths.science`/`outputPaths.coadd_science`/etc.
  directly as the `sci_path` argument to the existing
  `pypeit.outputfiles.spec_output_file(...)` (and similar helpers) — see
  Phase 2 below for how this plays out concretely at each call site.
- `parent`/`name` are always plain path arithmetic from already-known
  inputs (`redux_path`, `scidir`, `qadir`, `calib_dir`, `coadd_suffix`,
  `collate_outdir`) — none of them depend on another record's resolved
  `full`, so there's no recursive "resolve the parent property first"
  logic needed; `_apply()` computes all nine records' `parent`/`name`
  directly and independently.
- `redux` and `collate` are structurally identical to every other managed
  directory (both decompose a full-path input into `parent`/`name` via
  `Path(...).parent`/`.name`) — `collate` is fully independent of `redux`;
  when `collate_outdir` isn't given, `_apply()` simply reuses the resolved
  `redux_path` value as `collate_outdir`'s input, as if the user had passed
  the same string for both, rather than aliasing the two records together.
- `ready` means *only* "this directory has been created on disk (or
  confirmed to already exist) and is available to receive files" — it says
  nothing about configuration safety. That question is instead answered at
  the object level, by `configured` (see next bullet). This is a
  deliberate change from an earlier draft of this design, which conflated
  the two meanings in a single per-path `frozen` flag.
- `configured` (object-level) is what guards against silent reconfiguration
  now: `configure()` raises if called a second time, unless `dryrun` is
  `True`. There is no more idempotent-no-op comparison logic (an earlier
  draft's "guard 2") — it's unnecessary now that no downstream orchestrator
  is expected to call `configure()` more than once at all (§0.1 guideline
  3), and `configure()` always fully resets internal state rather than
  comparing against what was there before.
- A `_get()` call's `mkdir` (and the `ready` flip) only happens if
  `self._configured and not self._dryrun`; before `configure()` has ever
  been called, or whenever `dryrun` is `True`, a property still resolves
  and returns `full` (useful for inspection), it just never touches disk or
  becomes `ready`.
- (Formerly) guard 4 (filesystem sentinel file) remains intentionally
  absent — not stubbed, not mentioned as a TODO — per the explicit decision
  to skip it.

## 3. `pypeit/__init__.py` wiring

Insert immediately after the existing `dataPaths` instantiation
(`pypeit/__init__.py:24-25`):

```python
# Import and instantiate the data path parser
# NOTE: This *MUST* come after log and __version__ are defined above
from .pkg.pypeitdata import PypeItDataPaths
dataPaths = PypeItDataPaths()

# Import and instantiate the output path manager
from .pkg.outputpaths import PypeItOutputPaths
outputPaths = PypeItOutputPaths()          # cwd-based defaults, unconfigured

# Import all the exceptions so that they can be directly imported (e.g., `from
# pypeit import PypeItError`) in all package imports.
from .pkg.exceptions import *
```

No import cycle, and no absolute `pypeit`-package imports: `pkg/outputpaths.py`
only imports the standard-library `logging` module (retrieving the
already-configured `"pypeit"` logger by name, per §2) and its `pkg`-sibling
`exceptions` module via a relative import (`from .exceptions import
PypeItPathError`) — it never does `from pypeit import ...`, so it carries
none of the partially-initialized-package risk that `pkg/logger.py`'s
docstring warns about.

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
  must normalize **every** *live* branch to the `'PNGs/...'` form — but see
  §4.1: most of the affected branches turn out to be unreachable in
  practice and are removed outright rather than normalized.
- **`find_objects.py` is only half-broken**: `get_findobj_qa_filename`
  (~line 702) already correctly builds
  `os.path.join(par['rdx']['redux_path'], par['rdx']['qadir'])`. The
  hardcoded `os.path.join(..., 'QA')` (~line 1252) lives inside a dead
  `if False:` block (confirmed at ~line 1250, guarded by a `# TODO ::` note
  that flexure QA was never finished). **Decision: leave this dead code in
  place** rather than cleaning it up now — it isn't one of the 9 catalogued
  live violations, since `if False:` already makes it unreachable. Phase 3
  (§5) is scoped only to *confirm* it stays inaccessible (the `if False:`
  guard must not be incidentally touched or re-enabled by any nearby edit
  in this file), not to fix or delete it.
- **Design principle (see §0.1 for the full statement)**: minimize the
  number of places that can override any given output path, with overrides
  limited to the script level wherever practical, and high-level classes
  reading `outputPaths` directly instead of accepting their own path
  arguments. This is squarely within what the original design already
  permits (`pypeit_output_paths.md` §1: the hard constraint that paths
  must be passed as plain arguments applies only to `pypeit/core/`
  functions, never to higher-level classes like `Calibrations`).
- **`Calibrations.__init__`** (`calibrations.py:144-181`), concretely: per
  §0.1 guideline 2, **remove both the `caldir` and `qadir` arguments** from
  the signature, not just `qadir` as first proposed — the project lead
  flagged that `caldir` is exactly the same category of risk as `qadir`
  (a path argument that could be given a value inconsistent with
  `outputPaths`). Inside `__init__`, import `outputPaths` directly (`from
  pypeit import outputPaths`) and derive `self.calib_dir =
  outputPaths.calibrations` and `self.qa_path = outputPaths.qa` from it
  rather than from passed-in arguments. The `# TODO: This should only be
  defined in one place!` comment (line 177) is deleted, since `outputPaths`
  is now that one place — for both paths, not just QA.
  **Behavioral nuance found while checking the `qadir` removal**: today's
  `qadir=None` default is not just "no override" — it's actively used at
  the second `Calibrations.get_instance(...)` call site,
  `pypeit_steps.load_calibrations_for_frame` (~line 519), specifically to
  **skip** writing QA when reloading already-calibrated results during
  science reduction (QA was already written once, during the earlier
  calibration step). Simply deleting `qadir` and unconditionally deriving
  `self.qa_path` would silently start re-writing QA at that reload call
  site too. Recommended fix: replace the removed `qadir` *path* argument
  with a plain `write_qa: bool = True` *flag* argument — still not a path
  override (consistent with §0.1 guideline 2, since it's a boolean, not a
  path), but preserving the existing on/off behavior: `calib_one` uses the
  default `True`; `load_calibrations_for_frame` passes `write_qa=False`
  explicitly. `caldir`'s removal has no equivalent nuance — every call
  site always needs the calibrations directory, so there's no analogous
  "skip" flag to preserve for it.
- **Correction to `pypeit_output_paths.md`'s current-state analysis**: that
  document (§2.2) refers to `pypeit/core/exposure.py`. This module does
  **not** live under `pypeit/core/` — it's the top-level `pypeit/exposure.py`
  (confirmed: `pypeit/core/exposure.py` does not exist; `grep -n
  "^import\|^from" pypeit/pypeit.py` shows `from pypeit import exposure`).
  This matters here: because `exposure.py` is *not* `core/`, §0.1
  guideline 1's hard constraint doesn't apply to it — it's exactly the
  kind of high-level, non-core module guideline 2 is about. Its
  `reduce_exposure(...)` and `save_exposure(...)` functions
  (`exposure.py:350-...`, `:468-...`) both currently take a
  `calibrations_path: str` argument for the sole purpose of locating
  calibration files — per guideline 2, this should be removed in favor of
  `exposure.py` importing `outputPaths` and reading
  `outputPaths.calibrations` directly (see Phase 1, §5, for the knock-on
  simplification this enables in `pypeit_steps.calib_one`/`reduce_calibID`).
- **`pypeit_steps.calib_one`** (`pypeit_steps.py:124-179`) currently
  resolves `qa_path` itself (lines ~162-163:
  `os.path.join(par['rdx']['redux_path'], par['rdx']['qadir'])`) and
  forwards it via `qadir=qa_path`, and also threads a `calibrations_path`
  argument through to `Calibrations.get_instance(...)` as `caldir`. Once
  `Calibrations` and `exposure.py` both read `outputPaths` directly
  (above), **no function in this call chain needs to resolve or forward
  either path as an argument anymore** — `calib_one`, `reduce_calibID`
  (`pypeit.py:250-...`), and `PypeIt.calibrate_all`/`reduce_all`'s calls
  into them all drop their `calibrations_path` parameter/argument
  entirely. `calib_one` doesn't need to do anything explicit about
  directory creation either: accessing `outputPaths.calibrations`/
  `outputPaths.qa` as properties inside `Calibrations.__init__` (§2) both
  resolves *and* creates each directory on first access, so simply
  constructing `Calibrations` is sufficient (see Phase 1, §5).
- **`collate.py`**: `par['collate1d']['outdir']` is read directly at
  multiple call sites (lines 219, 229, 693, 802, and via
  `copy_spec1d_to_outdir` at 585-601, which does its own `os.makedirs` at
  line 597); the spec1d→spec2d string-replace is at line 668
  (`filename.replace('spec1d', 'spec2d', 1)`). Note `collate.py` also has a
  *second*, independent `par['collate1d']['spec1d_outdir']` key (lines 767,
  777, 780) used only for copying spec1d files.
  **Previously flagged as an open question, now resolved by §0.1**: per
  guideline 3, the `scripts/collate_1d.py` script class (`Collate1D.main`)
  is the right place to call `outputPaths.configure(par,
  collate_outdir=par['collate1d']['outdir'])` exactly once, immediately
  after `par` is built (~line 194, before its existing `outdir =
  par['collate1d']['outdir']` / `os.makedirs(outdir, ...)` at lines
  196-197). Per guideline 2, `collate.py`'s functions should then read
  `outputPaths.collate` directly at each of their `outdir`-consuming call
  sites, rather than re-reading `par['collate1d']['outdir']` independently
  in each one. `spec1d_outdir` remains a distinct concept, left alone, as
  before.

### 4.1 Dead `set_qa_filename` methods (checked against the full codebase)

Of the `case` branches in `qa.set_qa_filename` that actively hardcode the
double-`'QA/PNGs/...'` prefix (i.e., the ones flagged above as needing
normalization), a repo-wide search (`grep -rn` across all of `pypeit/`,
excluding `qa.py`'s own definition, plus the PypeIt-development-suite repo)
for each `method` string found that **five of the six are never invoked
anywhere in the codebase**:

| `method` value | Called elsewhere in `pypeit/`? | Called in dev-suite? |
|---|---|---|
| `'slit_profile_qa'` | No | No |
| `'plot_orderfits_Arc'` | No | No |
| `'pca_plot'` | No | No |
| `'pca_arctilt'` | No | No |
| `'plot_orderfits_Blaze'` | No | No |
| `'spat_flexure_qa_corr'` | **Yes** — `pypeit/images/rawimage.py:800` | — |

Only `'spat_flexure_qa_corr'` is live; its double-`QA/QA` bug must still be
fixed (normalized to `'PNGs/...'`) per Phase 3. The other five —
`'slit_profile_qa'`, `'plot_orderfits_Arc'`, `'pca_plot'`, `'pca_arctilt'`,
and `'plot_orderfits_Blaze'` — have no live call site anywhere in `pypeit/`
or the development suite. Rather than spend effort normalizing dead
branches, **these five `case` options will be removed outright from
`set_qa_filename` as part of Phase 3**, not merely fixed. This shrinks the
function and removes speculative/unreachable code along with the path bug,
rather than perpetuating it. (Every other `method` value in the function —
`slit_trace_qa`, `arc_fit_qa`, `arc_fwhm_qa`, `arc_fit2d_global_qa`,
`arc_fit2d_orders_qa`, `arc_tilts_spec_qa`, `arc_tilts_spat_qa`,
`arc_tilts_2d_qa`, `obj_trace_qa`, `obj_profile_qa`, `spec_flexure_qa_corr`,
`spec_flexure_qa_sky`, `spatillum_finecorr`, `detector_structure` — already
returns the correct `'PNGs/...'` form and is unaffected by this change.)

Phase 3 (§5) should be read as: normalize `spat_flexure_qa_corr`, delete
the five dead branches, leave the already-correct branches untouched.

## 5. Phased implementation plan

The `configured`/`dryrun`/`ready` lifecycle (§2) lands in **Phase 0** with
the class itself; it gets exercised against real orchestrator behavior in
**Phase 1**, where `PypeIt.__init__`'s single `configure()` call becomes
the one configuration point for the whole reduction workflow.

**Note on verification scope, per project direction**: this entire project
(Phases 0-7) is implemented and submitted as a **single PR** (see §7) —
the phase breakdown below is a development/review sequence, not a
PR-per-phase plan. Each phase's "Verify" step below covers only unit-level
testing (`pypeit/tests/`), performed as that phase is completed.
**PypeIt-development-suite (real instrument data) regression testing is
not run per phase — it happens exactly once, after Phase 7 is complete**,
as the final gate before the PR is submitted (§7 describes the
consolidated dev-suite matrix). Any per-phase "Verify" bullet below that
mentions a dev-suite setup should be read as "this is what the final,
single dev-suite pass in §7 needs to cover on behalf of this phase," not
as a phase-by-phase dev-suite run.

### Phase 0 — Introduce the module (no behavior change elsewhere) — ✅ COMPLETE
**Files:** new `pypeit/pkg/outputpaths.py`; `pypeit/__init__.py`.
- Add the class exactly as in §2. Wire the singleton per §3.
- Nothing else imports it yet.
- **Verify:** new `pypeit/tests/test_outputpaths.py` (§6). Run
  `pytest pypeit/tests/test_outputpaths.py`. Confirm
  `python -c "import pypeit; print(pypeit.outputPaths)"` succeeds with no
  import-cycle error.

**Completed at commit `90f3ae2376eb4ce5bee2af42a8cfc26f8578ec30`** ("Phase 0
implementation", 2026-07-17), touching exactly the three files above (plus
the docstring/signature polish from the same-day follow-up, folded into
that commit's content). All 9 tests in `test_outputpaths.py` pass; the
import-cycle check succeeds.

### Phase 1 — Adopt in `PypeIt`, `Calibrations`, and `exposure.py` — ✅ COMPLETE
**Files:** `pypeit.py`, `calibrations.py`, `pypeit_steps.py`, `exposure.py`.

**On guideline 3 (§0.1) for this phase specifically**: the ideal per
guideline 3 would be for `scripts/run_pypeit.py`'s `RunPypeIt.main()` to
call `outputPaths.configure(par, ...)` once, before constructing `PypeIt`.
That isn't possible here: for a `run_pypeit` invocation, the `PypeItPar`
object doesn't exist until `PypeIt.__init__` parses the `.pypeit` file
(via `inputfiles.PypeItFile.get_pypeitpar()`) — the script itself never
sees a `par` object at all. `PypeIt.__init__` is therefore the *earliest*
point at which guideline 3 can be satisfied for this entry point, not a
deviation from it. (This is different from `PypeIt`'s notebook/library
use, where a caller could in principle build `par` and call `configure()`
before constructing `PypeIt` — but requiring that of every direct caller
would be a worse ergonomic tradeoff than having `PypeIt.__init__` do it
once itself, so `PypeIt.__init__` remains the single configuration point
for both cases.)

- `pypeit.py` `PypeIt.__init__`: after the existing `redux_path` override
  (lines 87-88, which writes `self.par['rdx']['redux_path'] = redux_path`),
  call `outputPaths.configure(self.par, redux_path=redux_path,
  caller='PypeIt.__init__')` — this is the one and only `configure()` call
  in the entire reduction workflow for a given run; everything else in
  this phase (and Phases 3-5) only reads `outputPaths`. Replace line 119's
  `self.calibrations_path = Path(...) / ...` with
  `self.calibrations_path = outputPaths.calibrations` (kept as a
  convenience attribute for the `.calib` association-summary file and
  logging at lines 128-136 — not an override, just a cached read).
- Change the `science_path` property (lines 144-147) to
  `return outputPaths.science` and the `qa_path` property (lines 149-152)
  to `return outputPaths.qa` — **this fixes the existing `str`-vs-`Path`
  type inconsistency** (both now uniformly return `Path`). Downstream
  consumers (`qa.gen_qa_dir`, `qa.gen_mf_html`) already coerce their input
  through `pathlib.Path(...)`, so this is backward-compatible.
- `build_qa` (lines 154-163): simply access `outputPaths.qa_pngs` (the
  property resolves and creates it on first access, per §2) in place of
  relying on `qa.gen_qa_dir`'s own `mkdir` — no explicit creation call
  needed.
- `calibrations.py` `Calibrations.__init__` (144-181): per §0.1 guideline 2
  and §4, **remove both the `caldir` and `qadir` keyword/positional
  arguments**, adding `write_qa: bool = True` in `qadir`'s place. Add
  `from pypeit import outputPaths` at module level. Replace the
  `caldir`/`qadir`-derived block (lines 169-181) with:
  ```python
  self.calib_dir = outputPaths.calibrations   # property access creates it
  self.write_qa = write_qa
  self.qa_path = outputPaths.qa if write_qa else None   # ditto
  ```
  (deleting the `# TODO` at line 177 along with it, and the two inline
  `mkdir` calls at lines 171-172/180-181 — the property accesses above
  already create both directories on first use, per §2's `_get()`).
  `Calibrations` is now fully self-sufficient for both of its output
  directories — no caller needs to create them or pass them in, and no
  explicit creation call is needed anywhere.
- `pypeit_steps.py` `calib_one` (124-179): delete the `qa_path` resolution
  (lines ~162-163) and the `calibrations_path` parameter entirely — neither
  is needed anymore, since `Calibrations` now resolves both of its own
  paths. Drop the `qadir=qa_path` keyword and the positional
  `calibrations_path` argument from the `Calibrations.get_instance(...)`
  call (the class no longer accepts either).
- `pypeit_steps.py` `load_calibrations_for_frame` (~line 481-524, the
  reload-only call site identified in §4): drop its `calibrations_path`
  parameter for the same reason, and add `write_qa=False` explicitly to
  its `Calibrations.get_instance(...)` call, preserving today's behavior
  of not re-writing QA plots when reloading already-calibrated results
  (today this happens implicitly via the omitted `qadir` argument
  defaulting to `None`).
- `pypeit.py` `reduce_calibID` (~line 250) and its two call sites in
  `PypeIt.calibrate_all`/`reduce_all` (~lines 190, 216, 228): drop the
  `calibrations_path` parameter/argument throughout this chain — per §4's
  correction, `exposure.py` is not `core/`, so it too can read
  `outputPaths` directly rather than receiving this path as an argument
  (next bullet).
- `exposure.py` `reduce_exposure` (~line 350) and `save_exposure`
  (~line 468): remove the `calibrations_path: str` parameter from both
  signatures; add `from pypeit import outputPaths` at module level and
  read `outputPaths.calibrations` directly wherever the parameter was
  used. Update `reduce_calibID`'s calls to both functions (~lines 361,
  372) to drop the now-removed argument.
- **Verify:** `pytest pypeit/tests/test_calibrations.py
  pypeit/tests/test_runpypeit.py`; one MultiSlit dev-suite setup (e.g.
  `shane_kast_blue`) end-to-end, confirming `Science/`, `QA/PNGs/`,
  `Calibrations/` land in the right place; that the single `configure()`
  call in `PypeIt.__init__` is the only one in this phase (nothing
  downstream calls `configure()` again, so `outputPaths.configured` stays
  locked for the rest of the run); and that QA plots are still skipped
  (not regenerated) on the reload-only path through
  `load_calibrations_for_frame`.

**Completed at commit `29bb4b4458b0db1dfecd40fd594c58b2fd498077`** ("Phase 1
complete", 2026-07-17). Actual scope was broader than the bullets above:
`calibrations_path` also had to be dropped from `pypeit_steps.py`'s
`process_one_det`/`findobj_on_det`/`extract_det` and `exposure.py`'s
`process_exposure`/`findobj_on_exposure`/`extract_exposure`, since all of
them forward the value further down this same call chain; two scripts
outside the original file list, `scripts/run_to_calibstep.py` and
`scripts/reduce_by_step.py`, needed their direct calls into
`pypeit_steps`/`exposure` updated to match; and `test_calibrations.py`
needed a new `output_paths` fixture (monkeypatches `calibrations.outputPaths`
with a freshly configured, real, non-dry-run instance) since its tests
construct `Calibrations`/`MultiSlitCalibrations` directly, bypassing
`PypeIt.__init__`. Implementing this phase also surfaced a real conflict
with the once-only `configure()` rule: `PypeIt` is legitimately constructed
more than once per process (`scripts/ql.py` builds a calibration `PypeIt`
and a science `PypeIt` in one run; `test_run_pypeit` calls `RunPypeIt.main()`
twice), which would have raised `PypeItPathError` on the second
construction. Fixed by adding a narrowly-scoped `force=True` parameter to
`configure()` (§2), used only by `PypeIt.__init__` and documented/commented
at both places as a specific carve-out, not a general pattern. All 51
targeted tests and the full `pypeit/tests/` suite (413 tests, excluding the
heavy network-dependent `test_run_pypeit`) pass.

### Phase 2 — Wire call sites through `outputfiles.py` directly — ✅ COMPLETE
**Files:** `outputfiles.py`; call sites identified in Phases 1, 4, 5 below.
- No signature changes to `science_path(par)` or `spec_output_file(...)`
  (both have existing callers, and neither is being wrapped by a method on
  `PypeItOutputPaths` — see §2). This phase is a consolidation checkpoint,
  not a behavior change: at each call site being touched by this refactor,
  pass `sci_path=outputPaths.science` (or `outputPaths.coadd_science` for
  coadd contexts) directly into `outputfiles.spec_output_file(...)`, rather
  than letting it recompute the path from `par`. Add a regression assertion
  that the two never diverge.
- **Verify:** new unit test asserting `outputPaths.science ==
  outputfiles.science_path(par)` for a `configure()`-d par.

**Completed at commit `9d0f21fd5de33482ebbbb6251eae907bdbd79756`** ("Phase 2
complete", 2026-07-17). Scope was limited to the call sites already present
in the files Phase 1 touched: `exposure.py`'s `save_exposure` (all three
`spec_output_file` calls, plus replacing the `science_path(par)` +
manual `.mkdir()` check with a direct `outputPaths.science` access) and
`pypeit.py`'s `reduce_calibID`. `scripts/ql.py`'s two `spec_output_file`
call sites were left untouched, as planned, since that file belongs to
Phases 4/5. The regression assertion called for above already existed from
Phase 0's `test_configure_from_par` (`outputPaths.science ==
outputfiles.science_path(par)`); no new test was needed. Full
`pypeit/tests/` suite (413 tests, excluding the heavy `test_run_pypeit`)
passes.

### Phase 3 — Fix `qa.py`, `extraction.py`, `find_objects.py`, `core/findobj_skymask.py` — ✅ COMPLETE
**Files:** `qa.py`, `extraction.py`, `find_objects.py`,
`core/findobj_skymask.py`.
- `qa.set_qa_filename` (27-139): per §4.1, **delete** the five dead
  `case` branches (`slit_profile_qa`, `plot_orderfits_Arc`, `pca_plot`,
  `pca_arctilt`, `plot_orderfits_Blaze`) rather than fixing them, and
  normalize the one live offender (`spat_flexure_qa_corr`) to `'PNGs/...'`.
  Delete the `# TODO` (lines 25-26). Signature/argument contract
  (`out_dir`) is unchanged — this remains a pure naming helper callable
  from any layer, no `outputPaths` import needed here.

  **Superseded by a same-day follow-up request** (see the completion note
  below): `out_dir` was removed from the signature entirely, and the
  function now returns only a bare file name. Callers join it with
  `outputPaths.qa_pngs` (or their own already-available QA directory)
  themselves.
- `qa.gen_exp_html` (494-520): add a required `qa_path` argument; replace
  the hardcoded `Path("QA")` (~line 500) and `f"QA/{uni_name}.html"`
  (~line 506) with paths built under the passed `qa_path`. Update its sole
  caller, `PypeIt.build_qa`, to pass `self.qa_path`.
- `extraction.py:584`: replace `os.path.join(self.par['rdx']['redux_path'],
  'QA')` with `os.path.join(self.par['rdx']['redux_path'],
  self.par['rdx']['qadir'])` — matching the pattern already correct in
  `find_objects.py:702`.

  **Superseded by a same-day follow-up request**: replaced entirely with
  `outputPaths.qa` (imported directly), rather than reading
  `self.par['rdx']` at all.
- `find_objects.py`: **no change** — per §4, the dead `if False:` block
  containing the hardcoded `'QA'` (~lines 1250-1254) is left in place as
  unreachable code. The only action here is a manual check (not a code
  change) that this file's other edits in this phase don't alter the
  `if False:` guard.
- `core/findobj_skymask.py` (lines 1453-1455): remove the
  `qafile.parent.mkdir(parents=True, exist_ok=True)` side effect from
  `core/` entirely. The caller (`find_objects.py`, orchestration layer)
  already ensures the QA directory exists simply by having accessed
  `outputPaths.qa_pngs` earlier in Phase 1 (property access creates it on
  first use, per §2), so this becomes redundant, and its removal satisfies
  the "core/ never manages directories" principle exactly.
- **Verify:** extend `tests/test_qa.py` with a regression test asserting
  no `set_qa_filename` branch yields a doubled `QA/QA` segment and every
  returned path contains exactly one `PNGs` component; dev-suite run of a
  QA-plot-generating setup (flexure + object-finding QA) confirming no
  crash from the removed `core/` mkdir and correct `QA/PNGs/` placement.

**Completed at commit `64231f23428177332f2ded3e4ccfcce5d7d590f1`** ("Phase 3
complete", 2026-07-17). This single commit includes both the phase as
originally planned above *and* two same-day follow-up requests that
changed the shape of the `qa.py` work substantially:

- **`extraction.py`**: rather than just substituting `qadir` for the
  hardcoded `'QA'` literal, the `self.par['rdx']['redux_path']`/`qadir`
  computation was replaced outright with `outputPaths.qa` (imported
  directly), consistent with §0.1 guideline 2.
- **`qa.set_qa_filename` no longer takes `out_dir` or returns a full path**
  — it now returns only a bare file name, with every `case` branch
  converted to an f-string and a direct `return` (no more
  `outfile = ...` + a single `return str(pathlib.Path(out_dir) / outfile)`
  at the end). This is a bigger contract change than originally planned
  (which assumed `out_dir` stayed as-is), and it rippled through every
  caller:
  - `qa.py`'s own internal helpers (`html_mf_pngs`, `html_exp_pngs`,
    `arc_tilts_2d_qa`, `arc_tilts_spec_qa`, `arc_tilts_spat_qa`,
    `spec_flexure_qa`'s two call sites) keep their own pre-existing
    `out_dir` parameter (still the top-level QA path, unchanged contract
    for *their* callers, e.g. `wavetilts.py` was untouched) and now do
    the `'PNGs'` join locally via pathlib instead of relying on
    `set_qa_filename` to do it.
  - External callers — `find_objects.py`, `flatfield.py` (both call
    sites), `wavecalib.py` (all five call sites) — now import
    `outputPaths` directly and build
    `outputPaths.qa_pngs / qa.set_qa_filename(...)`, dropping their
    `out_dir=self.qa_path` arguments.
  - **One more live call site was found during this follow-up that the
    original Phase 3 research had missed** (grep scope had excluded
    `pypeit/images/`): `pypeit/images/rawimage.py:800`, the actual live
    caller of `spat_flexure_qa_corr`. It derived its `out_dir` via
    `Path(slits.calib_dir).parent` — which resolves to the *redux root*,
    not the QA directory — a likely pre-existing latent bug, since that
    QA PNG was probably never actually landing under `QA/PNGs/` at all.
    Replaced with `outputPaths.qa_pngs` directly, and removed the
    now-unused `from pathlib import Path` import from that file.
  - `tests/test_qa.py` rewritten to match: `test_set_qa_filename_no_double_qa`
    replaced by `test_set_qa_filename_bare_name_only` (asserts every
    surviving method returns a single-component name with no `QA`/`PNGs`
    embedded), and `test_set_qa_filename_dead_methods_removed`'s calls
    drop the now-nonexistent `out_dir` keyword.

Full `pypeit/tests/` suite (415 tests, excluding the heavy
`test_run_pypeit`) passes.

### Phase 4 — Fix `coadd2d.py` and its callers
**Files:** `coadd2d.py`, `scripts/coadd_2dspec.py`, `scripts/ql.py`.

**On guideline 3 (§0.1) for this phase**: unlike `run_pypeit.py` (Phase 1),
the script class here already has `par` in hand well before it needs any
coadd path — `CoAdd2DSpec.main()` (`scripts/coadd_2dspec.py:36-...`) builds
`spectrograph, par, _ = coadd2dFile.get_pypeitpar(...)` at ~line 62, and
only calls `coadd2d.CoAdd2D.output_paths(...)` at ~line 85. So this phase
*can* and does follow guideline 3 exactly: the **script** calls
`outputPaths.configure(...)` once, and `CoAdd2D`/`output_paths()` itself
never calls `configure()` — it only reads already-resolved properties.
(The version of this plan written before §0.1 existed had
`output_paths()` itself call `configure()`; that's corrected here.)

- `CoAdd2D.output_paths` (475-516): delete the `par['rdx']['qadir'] +=
  '_coadd'` mutation (line 512) entirely — no replacement mutation of
  any kind, and no `configure()` call inside this method (see above).
  Replace the hardcoded `'Science'` (505), the `f"{scidir}_coadd"`
  construction (508), and the hardcoded `'PNGs'` (513) with reads of
  `outputPaths.coadd_science` / `outputPaths.coadd_qa_pngs` — accessing
  each property both resolves and creates the corresponding directory on
  first use (§2), so no separate creation call is needed. The staticmethod
  keeps its existing call signature for API stability, but its body now
  assumes `outputPaths` has already been configured by its caller.
- `scripts/coadd_2dspec.py` `CoAdd2DSpec.main()`: immediately after `par`
  is built (~line 62), add `outputPaths.configure(par,
  redux_path=par['rdx']['redux_path'], caller='CoAdd2DSpec.main')` — the
  single configuration point for this script, satisfying guideline 3.
  Update the `output_paths(...)` call (84-85) to the resulting
  non-mutating, non-configuring form; replace the hand-rolled
  `spec1d_{basename}.fits`/`.txt` (229, 233) and `spec2d_{basename}.fits`
  (238) with a shared coadd-naming helper in `outputfiles.py` (new, small
  function reused by both this script and `ql.py`, so the naming
  convention has exactly one implementation).
- `scripts/ql.py`: its `output_paths(spec2d_files, par)` call (~1022,
  called *without* `coadd_dir` and therefore also mutating `par` today)
  needs its own script class to call `outputPaths.configure(par,
  redux_path=par['rdx']['redux_path'])` once beforehand, analogous to
  `CoAdd2DSpec.main()` above, before switching to the non-mutating form;
  its hand-rolled `spec2d_{basename}.fits` (~1024) is replaced by the same
  shared helper used in `coadd_2dspec.py` (its sibling call at ~1028
  already correctly uses `outputfiles.spec_output_file`, confirming the
  target pattern).
- **Verify:** unit test asserting `par['rdx']['qadir']` is bit-for-bit
  unchanged after computing coadd paths (regression guard for the removed
  mutation); a second unit test asserting `CoAdd2D.output_paths` itself
  never calls `outputPaths.configure` (e.g. via a mock/spy), enforcing
  guideline 3 for this call site; dev-suite runs of a `coadd_2dspec` setup
  and a quicklook (`ql`) setup confirming `Science_coadd/` and
  `QA_coadd/PNGs/` are created and spec1d/spec2d names match the shared
  helper in both scripts.

### Phase 5 — Remaining scripts and `collate`
**Files:** `scripts/trace_edges.py`, `collate.py`, `scripts/collate_1d.py`,
`scripts/sensfunc.py`, `scripts/flux_calib.py`, `scripts/coadd_1dspec.py`,
`scripts/flux_setup.py`.

Per guideline 3 (§0.1), every script class in this phase calls
`outputPaths.configure(...)` exactly once, as early as possible after its
`par` is built, and no function it calls does so again — the same pattern
established in Phases 1 and 4.

- `trace_edges.py`: remove the duplicated argparse literal
  `default='Calibrations'` (line 42); source the default from
  `CalibrationsPar`/`outputPaths` instead. Unify the two QA-path branches
  (line 127, correct; line 154, hardcoded `redux_path / 'QA'`) so both
  route through a single `outputPaths.configure(redux_path=redux_path)`
  call in the script class, then read `outputPaths.qa`. Update
  `calib_dir` construction (line 163) to `outputPaths.calibrations`.
- `collate.py`/`scripts/collate_1d.py`: per §4's resolution of the
  previously-open question, `Collate1D.main()`
  (`scripts/collate_1d.py:189-...`) calls `outputPaths.configure(par,
  collate_outdir=par['collate1d']['outdir'], caller='Collate1D.main')`
  once, immediately after `par` is built (~line 194) and before its
  existing `outdir = par['collate1d']['outdir']` / `os.makedirs(outdir,
  ...)` (~lines 196-197) — the latter is simply deleted, since the first
  access to `outputPaths.collate` anywhere below creates the directory
  automatically (§2). All of `collate.py`'s own functions (lines 219, 229,
  693, 802, and
  `copy_spec1d_to_outdir` at 585-601) then read `outputPaths.collate`
  directly instead of independently re-reading
  `par['collate1d']['outdir']` at each call site — this is the concrete
  fix that resolves last turn's open question via guideline 2. Leave
  `par['collate1d']['spec1d_outdir']` (lines 767, 777, 780) untouched —
  it's a distinct, already-working concept, out of scope here. Replace
  the string-replace spec1d→spec2d translation (line 668) with the same
  shared `outputfiles.py` naming helper introduced in Phase 4.
- `sensfunc.py`/`flux_calib.py`/`coadd_1dspec.py`/`flux_setup.py`: each
  script class calls `outputPaths.configure(par, redux_path=...)` once
  (defaulting to cwd, `outputPaths`'s own default, for genuinely
  standalone invocations — no behavior change for today's cwd-relative
  usage), then builds `--par_outfile` defaults and product files
  (`sens_*`, `.flux`/`.coadd1d`/`.tell`, `coadd1d_*`) from
  `outputPaths.redux` rather than the implicit cwd join used today.
  No function below the script class re-derives or overrides these paths.
  Fix the `flux_calib.py:62` `action="store_true"` combined with a string
  `default='fluxing.par'` (a pre-existing, unrelated argparse bug — the
  flag can never actually take a filename value) as a drive-by while
  touching this line.
- **Verify:** dev-suite regression for each of `collate_1d`, `sensfunc`,
  `flux_calib`, `coadd_1dspec`, `flux_setup`, `trace_edges` post-processing
  setups; unit test confirming `collate.py`'s functions no longer read
  `par['collate1d']['outdir']` directly (only `outputPaths.collate`).

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
  accessing `outputPaths.qa_pngs` has fully replaced its role).
- **Verify:** full `pytest pypeit/tests/` (dev-suite regression deferred
  to the single final pass, §7).

### Phase 7 — Update `doc/` to reflect the new output-paths system
**Files:** `doc/outputs.rst` (primary); `doc/qa.rst`, `doc/collate1d.rst`,
`doc/coadd2d.rst`, `doc/quicklook.rst` (consistency review); optionally
`pypeit/par/pypeitpar.py` docstrings + the auto-generated
`doc/pypeit_par.rst`, if a docstring is revised (see below).

- `doc/outputs.rst`'s "Directory Structure" section (`.. _outputs-dir:`
  anchor, currently lines ~17-46):
  - Remove or rewrite the `.. warning::` block stating that overriding
    `scidir`/`qadir`/`calib_dir`/`redux_path` "is not well tested" — after
    Phases 0-6, this is no longer true: `PypeItOutputPaths` (§2) is
    exactly what makes these overrides consistent and exercised by tests.
    Replace with a short, accurate description of the resolution model
    (`redux_path` + `scidir`/`qadir`/`calib_dir`, all existing `PypeItPar`
    keys with unchanged names/defaults — nothing here is new to the user,
    only more reliable).
  - Resolve the `.. TODO: INCLUDE COADD DIRECTORY STRUCTURE` comment
    (~line 45): document the coadd-specific directories produced by
    `pypeit_coadd_2dspec`/quicklook (`${scidir}${coadd_suffix}`/
    `${qadir}${coadd_suffix}`, i.e. `Science_coadd`/`QA_coadd` by default),
    referencing Phase 4's fix and noting the suffix is no longer a silent
    hardcode.
  - Consider a brief cross-reference to `doc/collate1d.rst` for
    `pypeit_collate_1d`'s separately-scoped `outdir` parameter, so the
    page gives a complete picture of every output-path concept covered by
    this project (§2 of `pypeit_output_paths.md`), not just the three
    primary ones.
- Review (not necessarily rewrite) `doc/qa.rst`, `doc/collate1d.rst`,
  `doc/coadd2d.rst`, and `doc/quicklook.rst` for continued accuracy. None
  are expected to need substantive rewrites — no user-facing parameter
  name, default, or semantic changes anywhere in Phases 0-6 (the
  backward-compatibility invariant holds throughout) — but each should be
  read against the final implementation to catch any wording that
  describes the old hardcoded/duplicated behavior as a quirk to work
  around, now that it's fixed.
- If Phase 1's removal of the `# TODO` comments prompts also revising the
  `CalibrationsPar.calib_dir` docstring's "Beware that success when
  changing the default value is not well tested!" caveat
  (`pypeit/par/pypeitpar.py`) to reflect the new, tested behavior, that is
  a source-code docstring change, not a `doc/` edit — regenerate
  `doc/pypeit_par.rst` afterward via this project's standard
  parameter-documentation pipeline (the `build-docs` skill /
  `cd doc; make htmlonly`) rather than hand-editing the auto-generated
  table.
- **Verify:** rebuild the documentation (`cd doc; make htmlonly` for a
  build that doesn't require `PYPEIT_DEV`/internet access, or the full
  `cd doc; make clean; make html` if those are available) and confirm it
  builds cleanly with no new warnings in the touched pages; visually
  check the rendered `outputs.html` directory-structure section —
  including the new coadd-directory description — against the actual
  on-disk output of a representative run from the final dev-suite pass
  (§7), to confirm the documentation matches real behavior, not just the
  code.

## 6. Test plan

Convention (per `pypeit/tests/`): plain module-level `def test_*()`
functions, no `TestX` classes, fixed RNG seeds, unit-only (no large data
files). Relevant existing files: `tests/test_pypeitpar.py`,
`tests/test_calibrations.py`, `tests/test_runpypeit.py`, `tests/test_qa.py`,
`tests/test_coadd.py`.

**New file `pypeit/tests/test_outputpaths.py`:**
- `test_default_resolution(tmp_path)` — a freshly constructed, unconfigured
  instance (`redux_path=tmp_path`) still resolves every property to the
  expected `Path` (`science == tmp_path/'Science'`, `qa_pngs ==
  tmp_path/'QA'/'PNGs'`, `calibrations == tmp_path/'Calibrations'`,
  `collate == tmp_path` (mirrors `redux_path` when `collate_outdir` isn't
  given), `coadd_science == tmp_path/'Science_coadd'`, `coadd_qa_pngs ==
  tmp_path/'QA_coadd'/'PNGs'`), but none of them exist on disk and
  `._paths[key].ready` is `False` for all keys, since the instance was
  never `configure()`-d.
- `test_explicit_dirs(tmp_path)` — non-default `scidir`/`qadir`/
  `calib_dir`/`coadd_suffix`/`collate_outdir` (a different directory from
  `redux_path`) all resolve correctly, and `collate`'s record is verified
  to be an independent `_ManagedPath` from `redux`'s (different `parent`),
  per the "independent, not aliased" design decision.
- `test_configure_creates_on_access(tmp_path)` — after `configure()`,
  accessing `.science` for the first time creates the directory on disk
  and flips `._paths['science'].ready` to `True`; a second access doesn't
  re-log or re-`mkdir` (assert via `caplog`/a monkeypatched `Path.mkdir`
  call count of 1).
- `test_nested_pngs_creates_parent(tmp_path)` — accessing `.qa_pngs` before
  ever accessing `.qa` directly still creates *and marks ready* the parent
  `QA/` directory as well as `QA/PNGs/`, confirming `parent`/`name` being
  set eagerly by `configure()` (not lazily chained through another
  property) still yields the correct nested result.
- `test_configure_twice_raises(tmp_path)` — a second `configure()` call
  (not in dry-run) raises `PypeItPathError`, regardless of whether the new
  values are identical to the current ones (no idempotent no-op — see §2).
- `test_dryrun_configure_repeatable(tmp_path)` — with `dryrun=True`,
  `configure()` can be called more than once; each call fully resets
  `_paths` ("as if reinstantiated"); no property access ever sets `ready`
  or touches disk in dry-run mode, even after `configure()`.
- `test_configure_logs(caplog)` — a successful `configure()` call emits
  exactly one `log.info` recording the resolved `redux_path`.
- `test_configure_from_par()` — build a default `PypeItPar`, `configure()`
  from it, assert properties match `outputfiles.science_path(par)`.
- `test_derive_not_implemented()` — `derive()` raises `NotImplementedError`
  (documents its held-in-reserve status; see §2's note on its reduced but
  not eliminated rationale).

**Extensions to existing files:**
- `tests/test_calibrations.py`: `test_calibrations_no_path_overrides()` —
  confirm `Calibrations.__init__`/`get_instance` no longer accept `caldir`
  or `qadir` (e.g. via `inspect.signature`), and that `self.calib_dir`/
  `self.qa_path` match `outputPaths.calibrations`/`.qa` after
  construction; `test_calibrations_write_qa_flag()` — `write_qa=False`
  yields `self.qa_path is None` and never triggers `outputPaths.qa_pngs`
  to become `ready`.
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

**Single PR**, per project direction — Phases 0-7 (§5) are a development
and review sequence within that one PR, not a series of separate PRs. Each
phase is still completed and unit-tested in order (§5's per-phase
"Verify" steps), so the work remains incrementally reviewable commit-by-
commit even though it lands as one PR.

**Unit testing happens continuously, phase by phase**, exactly as
described in each phase's "Verify" step in §5 and in the test plan in §6 —
there is no change to that part of the plan.

**PypeIt-development-suite (real instrument data) regression testing is
deferred until every phase, including Phase 7's documentation updates, is
complete.** It is then run exactly once, as a single consolidated pass,
covering the union of what each phase's "Verify" step called for:

- One MultiSlit + one Echelle + one SlicerIFU setup (e.g. `shane_kast_blue`,
  a `keck_hires`/`keck_nires` setup, `keck_kcwi`), confirming `Science/`,
  `QA/PNGs/`, `Calibrations/` placement; that `-r/--redux_path` still
  overrides correctly; and, run against an existing `.pypeit` file with
  custom `[rdx] scidir/qadir` and `[calibrations] calib_dir` values, that
  the backward-compatibility invariant holds (Phase 1).
- A QA-plot-generating setup exercising the flexure/object-finding QA path
  (Phase 3).
- A `coadd_2dspec` setup and a quicklook (`ql`) setup, confirming
  `Science_coadd/`/`QA_coadd/PNGs/` placement and that neither script
  mutates `par` (Phase 4).
- One setup each for `collate_1d`, `sensfunc`, `flux_calib`,
  `coadd_1dspec`, `flux_setup`, and `trace_edges` (Phase 5).
- A broader multi-pypeline regression subset for final confidence
  (Phase 6), and a check that the Phase 7 documentation (in particular the
  rebuilt `doc/outputs.rst`) matches the real directory structure produced
  by this same dev-suite pass.

Running this as a single consolidated pass (rather than once per phase)
is reasonable here specifically because all intermediate phases are unit-
tested and merged into the same branch before any dev-suite run — there's
no intermediate state that ships or is reviewed independently, so there's
no benefit to fragmenting the dev-suite verification the way there would
be across separate PRs.

The changelog / in-development release notes are updated once, for the PR
as a whole, per this project's standard changelog convention, noting
explicitly that no user-facing `PypeItPar` key is renamed, removed, or
redefaulted by this work.

## 8. Change log

This section tracks substantive edits to this implementation plan itself,
and, once Phase 0 landed, keeps §2 in sync with the actual
`pypeit/pkg/outputpaths.py` source. Newest entries first.

### 2026-07-17 — Phase 3 complete (plus two same-day follow-up requests)

Phase 3 (§5) is implemented and committed (`64231f23428177332f2ded3e4ccfcce5d7d590f1`,
"Phase 3 complete"). See the completion note at the end of the Phase 3
section in §5 for full detail. Two follow-up requests, made and
implemented in the same commit, changed the shape of this phase
substantially beyond what was originally planned:

- **`extraction.py`**: rather than the originally planned
  `qadir`-substitution fix, its `out_dir` computation was replaced outright
  with a direct `outputPaths.qa` import/read.
- **`qa.set_qa_filename`**: reworked to drop `out_dir` entirely and return
  only a bare file name (every `case` branch converted to an f-string with
  a direct `return`), with callers now responsible for joining
  `outputPaths.qa_pngs` (or their own QA directory) themselves. This
  rippled through every direct caller across `qa.py`, `find_objects.py`,
  `flatfield.py`, and `wavecalib.py`, and surfaced one live call site the
  original Phase 3 research had missed (`pypeit/images/rawimage.py:800`),
  which also turned out to carry a likely pre-existing latent bug (QA PNG
  path derived from the wrong parent directory), now fixed.

Full `pypeit/tests/` suite (415 tests, excluding the heavy
`test_run_pypeit`) passes.

### 2026-07-17 — Phase 2 complete

Phase 2 (§5) is implemented and committed (`9d0f21fd5de33482ebbbb6251eae907bdbd79756`,
"Phase 2 complete"). See the completion note at the end of the Phase 2
section in §5 — scope was limited to the `exposure.py`/`pypeit.py` call
sites already touched by Phase 1; `scripts/ql.py` was correctly left for
Phases 4/5.

### 2026-07-17 — Phase 1 complete

Phase 1 (§5) is implemented and committed (`29bb4b4458b0db1dfecd40fd594c58b2fd498077`,
"Phase 1 complete"). See the completion note at the end of the Phase 1
section in §5 for what was actually touched, including the broader-than-
planned `calibrations_path` removal and the `force=True` carve-out added to
`PypeItOutputPaths.configure()` (§2) to allow constructing more than one
`PypeIt` object per process.

### 2026-07-17 — Phase 0 annotated with its completing commit

Marked Phase 0 (§5) as ✅ COMPLETE and recorded the commit hash it landed
in, `90f3ae2376eb4ce5bee2af42a8cfc26f8578ec30` ("Phase 0 implementation"),
for traceability as later phases build on it.

### 2026-07-17 — Phase 0 implemented; §2 synced to match the merged code

Phase 0 (§5) was implemented: `pypeit/pkg/outputpaths.py`,
`pypeit/__init__.py` wiring, and `pypeit/tests/test_outputpaths.py` (9
tests, all passing; `import pypeit` confirmed cycle-free). Two small
polish requests were then applied directly to the merged code and mirrored
back into §2 so the reference design doesn't drift from what's actually in
the repo:

- All docstrings converted from Google (`Args:`/`Returns:`/`Raises:`) to
  Numpy style (`Parameters\n----------`/`Raises\n------`), matching this
  project's stated docstring preference and the convention already used
  elsewhere (e.g. `alignframe.py`, `inputfiles.py`).
- `_ManagedPath`'s four fields each gained a short attribute docstring.
- `_PATH_KEYS` was removed — nothing in the class itself referenced it
  (only the test suite did, to iterate all managed-path keys); the test
  was updated to iterate `self._paths.values()`/`.keys()` directly instead
  of relying on a separate, independently-maintained constant.
- The `*` keyword-only marker was removed from both `__init__` and
  `configure`'s argument lists, so all parameters are now
  positional-or-keyword.

### 2026-07-17 — Class redesign: per-path lazy `parent`/`name`/`full`/`ready` records; object-level `configured`/`dryrun` replace the old freeze/idempotency guards

Reworked §2 end-to-end around a different lifecycle model proposed by the
project lead: rather than one global `_frozen` flag and an explicit,
separately-called `make()` method, every managed directory is now its own
small internal record (`parent`, `name`, `full`, `ready`), and creation
happens implicitly, the first time each directory's `@property` is
accessed — "make directories when they're needed" rather than requiring a
developer to remember a separate creation call at the right point in the
code. This followed two rounds of back-and-forth surfacing (and then
resolving) several real implications:

- **Shared-root consistency and nested directories**: resolved by
  promoting `redux_path` to its own first-class managed property
  (`redux`), and treating `collate` identically — both decompose a
  full-path input into `parent`/`name` via `Path(...).parent`/`.name`,
  exactly like every other managed directory. Nested directories
  (`qa_pngs` under `qa`, `coadd_qa_pngs` under `coadd_qa`) turned out not
  to need any recursive "resolve the parent property first" logic at all:
  `parent`/`name` for every one of the nine managed paths is set directly
  from already-known inputs at `configure()` time (pure path arithmetic,
  no I/O), so there's nothing to chain lazily.
- **The dual meaning of the old per-path `frozen` flag** ("safe to
  reconfigure?" vs. "ready to use?") is fully separated now: object-level
  `configured` (with a companion `dryrun` mode) answers the first
  question — `configure()` may only be called once per instance unless
  `dryrun=True`, full stop, with no idempotent-no-op comparison logic
  needed (an earlier draft's "guard 2" is gone entirely, since nothing
  downstream is expected to call `configure()` more than once anyway, per
  §0.1 guideline 3). Per-path `ready` (renamed from `frozen`) now means
  *only* "this directory has been created on disk (or already existed)
  and is available" — set solely by the property getter's creation logic,
  never by `configure()`.
- **`dryrun`** (new `__init__`/`configure()` parameter) lets a
  `PypeItOutputPaths` be constructed for inspection/testing without any
  filesystem side effects: `full` still resolves normally, but no `mkdir`
  ever happens and `ready` never becomes `True`. It also exempts an
  instance from the "configure once" restriction, so a dry-run instance
  can be reconfigured repeatedly (each `configure()` call fully resets all
  path records, "as if reinstantiated") — useful for iterative
  inspection/testing without needing a fresh object each time.
- **`make()` and `_MAKE_TARGETS` are removed entirely** — there's no
  longer a batch/explicit creation API to call. This rippled through every
  phase in §5 that previously called `outputPaths.make(...)`: those call
  sites are simply deleted, since the property access that was already
  happening nearby now does the creation implicitly. §6's test plan was
  rewritten to match (dropped the idempotency/freeze-reject tests; added
  tests for lazy creation-on-access, nested-parent creation, dry-run
  repeatability, and the strict "configure twice raises" behavior).
- `derive()` (guard 3) is left as an untouched, still-reserved
  `NotImplementedError` stub — flagged in §2 as having a weaker rationale
  now that `coadd_science`/`coadd_qa`/`collate` are ordinary properties on
  the same object, but not removed, since nothing said today contradicts
  keeping it in reserve.

### 2026-07-17 — Single-PR rollout; dev-suite testing deferred to one final pass; added Phase 7 (documentation)

Per project direction on how this work will actually be executed and
submitted:

- **§7 rewritten for a single PR** covering all phases (0-7), rather than
  one PR per phase. Phases remain a development/review sequence within
  that PR.
- **Dev-suite (real instrument data) regression testing is deferred until
  every phase is complete**, then run once as a consolidated pass, rather
  than once per phase as the previous version of §7 described. Unit
  testing is unaffected — it still happens phase by phase, per §5/§6.
  Added a note at the top of §5 clarifying that per-phase "Verify" bullets
  mentioning dev-suite setups describe what the final single pass needs
  to cover on that phase's behalf, not a per-phase dev-suite run.
- **Added Phase 7** (§5): update `doc/` to reflect the new output-paths
  system, primarily `doc/outputs.rst`'s "Directory Structure" section
  (removing the now-inaccurate "not well tested" override warning and
  resolving its long-standing `TODO: INCLUDE COADD DIRECTORY STRUCTURE`
  comment), plus a consistency review of `doc/qa.rst`, `doc/collate1d.rst`,
  `doc/coadd2d.rst`, and `doc/quicklook.rst`, and a note on regenerating
  `doc/pypeit_par.rst` via the standard doc pipeline if any `PypeItPar`
  docstring text is revised along the way.

### 2026-07-17 — New §0.1 guidelines; `caldir` removed alongside `qadir`; `exposure.py` mislabel corrected; collate open question resolved

Prompted by the project lead pointing out that `caldir` carries the same
override risk as `qadir` did (next entry below) — a sign the requirements
themselves needed supplementing, not just this one call site. Added a new
§0.1 with three guidelines governing where output paths may be overridden
going forward: (1) `core/` functions take paths as arguments (already in
place); (2) high-level classes should read `outputPaths` directly and lose
any path-override arguments; (3) `outputPaths.configure()` should be
called once, as soon as `PypeItPar` exists — preferably in a `scripts/`
script class — and never again afterward. Framed as guidelines with
case-by-case exceptions, not hard requirements, per explicit direction.

Reconsidered the whole document against these guidelines and made the
following adjustments:

- **`Calibrations.__init__`** now loses `caldir` as well as `qadir` (§4,
  Phase 1). Both `self.calib_dir` and `self.qa_path` are derived from
  `outputPaths` directly; `Calibrations` now creates both of its own
  output directories itself, rather than relying on a caller to create
  them first.
- **Corrected a factual error inherited from `pypeit_output_paths.md`**:
  `exposure.py` is `pypeit/exposure.py`, not `pypeit/core/exposure.py` — it
  is not a `core/` module, so guideline 1's hard constraint never applied
  to it. Per guideline 2, its `reduce_exposure`/`save_exposure` functions
  lose their `calibrations_path` argument in favor of reading `outputPaths`
  directly, which in turn lets `pypeit_steps.calib_one`/
  `load_calibrations_for_frame` and `pypeit.reduce_calibID` drop
  `calibrations_path` from their own signatures entirely — nothing in this
  call chain needs to thread that value as an argument anymore.
- **Phase 4 (`coadd2d.py`) corrected against guideline 3**: the version of
  this plan written before §0.1 existed had `CoAdd2D.output_paths()` itself
  call `outputPaths.configure()`. That's now moved to the calling script
  class (`CoAdd2DSpec.main()`, and `ql.py`'s script class), which already
  has `par` in hand before calling `output_paths()` — `output_paths()`
  itself now only reads already-configured properties.
- **Resolved last turn's open question about `collate.py`**: guideline 3
  puts `outputPaths.configure(par, collate_outdir=...)` in
  `Collate1D.main()`, once; guideline 2 has `collate.py`'s functions read
  `outputPaths.collate` directly instead of independently re-reading
  `par['collate1d']['outdir']` at each of their several call sites.
  `spec1d_outdir` remains untouched, as already agreed.
- **Phase 1 gained an explicit discussion of why `PypeIt.__init__` is the
  right (not an exceptional) place to call `configure()` for a
  `run_pypeit` invocation**: the script never sees a `PypeItPar` object at
  all for that entry point (`PypeIt.__init__` is what parses the
  `.pypeit` file), so `PypeIt.__init__` is the earliest point guideline 3
  can apply, not a deviation from it.
- Added `tests/test_calibrations.py` test bullets (§6) and tightened
  PR 2's description (§7) to mention `exposure.py` and the
  single-configure-point verification.

### 2026-07-17 — Leave `find_objects.py` dead code alone; remove `Calibrations.__init__`'s `qadir` argument

Three changes to §4, based on direction from the project lead:

- **`find_objects.py`'s dead `if False:` block** (hardcoded `'QA'`, ~lines
  1250-1254) will be left in place, not cleaned up, as part of this effort.
  Phase 3 (§5) now only checks that it stays inaccessible rather than
  fixing or deleting it.
- **New design principle** added to §4: minimize the number of places that
  can override any given output path, with a first-order goal of limiting
  overrides to the script level; PypeIt's high-level classes should read
  `outputPaths` directly rather than accept their own path-override
  arguments. Concretely, `Calibrations.__init__` **loses its `qadir`
  keyword argument** and instead imports `outputPaths` directly for
  `self.qa_path`. While implementing this, found that today's `qadir=None`
  default is load-bearing at a second call site
  (`pypeit_steps.load_calibrations_for_frame`), where it's relied on to
  *skip* QA regeneration when reloading already-calibrated results.
  Proposed fix: replace `qadir` with a `write_qa: bool = True` flag rather
  than dropping the on/off capability outright. Phase 1 (§5) was rewritten
  to match, including a new `load_calibrations_for_frame` bullet and an
  added dev-suite verification point for the reload-skip behavior.
  `caldir` is explicitly unchanged in this pass.
- **Flagged as unresolved**: whether/how the `collate.py`
  `par['collate1d']['outdir']` handling (last bullet of §4) should change
  in light of the new "minimize override points" principle — noted as an
  open question rather than answered. `par['collate1d']['spec1d_outdir']`
  remains explicitly out of scope and untouched, as before.

### 2026-07-17 — New §4.1: five of six `set_qa_filename` double-QA branches are dead code

Checked whether the `case` branches in `qa.set_qa_filename` that hardcode
the double-`'QA/PNGs/...'` prefix (flagged in §4's first bullet) are
actually reachable, by grepping the full `pypeit/` tree and the
PypeIt-development-suite repo for each `method` string. Result: only
`'spat_flexure_qa_corr'` is called anywhere (`pypeit/images/rawimage.py:800`);
`'slit_profile_qa'`, `'plot_orderfits_Arc'`, `'pca_plot'`, `'pca_arctilt'`,
and `'plot_orderfits_Blaze'` have no live call site at all. Added §4.1
documenting this, and revised Phase 3 (§5) to **delete** those five dead
branches outright rather than normalize them, while still normalizing the
one live offender. §4's first bullet was given a one-line pointer to §4.1.

### 2026-07-17 — Relocate to `pkg/`, relative imports, drop `spec_output_file`

Three changes to §2 (and the sections that reference it), made before any
code was written, based on direction from the project lead:

- **Module location**: `PypeItOutputPaths` moves from a proposed top-level
  `pypeit/outputpaths.py` to `pypeit/pkg/outputpaths.py`, alongside its two
  closest architectural precedents, `pypeitdata.py` and `logger.py` (both
  already in `pkg/`). §2's intro, §3's wiring code (`from .pkg.outputpaths
  import PypeItOutputPaths`), and the Phase 0 file list were updated
  accordingly.
- **Relative imports for `pypeit` modules**: the class no longer imports
  anything from the top-level `pypeit` package by absolute path.
  `PypeItPathError` is now imported as a relative sibling
  (`from .exceptions import PypeItPathError`), matching the existing
  pattern in `pypeitdata.py:66`. The `log` object is obtained via
  `logging.getLogger('pypeit')` (stdlib) rather than `from pypeit import
  log`, since `pkg/` modules deliberately avoid reaching back into the
  partially-initialized top-level package — a risk `pkg/logger.py`'s own
  module docstring explicitly warns about. This also means `outputpaths.py`
  no longer needs the deferred-import trick it previously used to avoid a
  cycle when calling into `outputfiles.py`.
- **Removed the `spec_output_file` convenience method**: the class will
  not wrap `pypeit.outputfiles.spec_output_file(...)`. Call sites will
  instead pass `outputPaths.science`/`outputPaths.coadd_science` directly
  as the `sci_path` argument to the existing `outputfiles.spec_output_file`
  function. §2's notes and Phase 2 (retitled "Wire call sites through
  `outputfiles.py` directly") were updated to describe this instead of a
  delegating method; this also further reduces how much of `outputfiles.py`
  the new class needs to know about, since it no longer imports it at all.

No other section's substance changed — the six-phase plan, test plan, and
PR sequencing in §5-§7 still apply as written, just against the corrected
module path.
