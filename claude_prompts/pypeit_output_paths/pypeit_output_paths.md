# PypeIt Output Paths: Analysis and Proposed Strategies

## 1. Original prompt

> The path definitions in PypeIt are not well maintained.  The paths for science outputs, calibrations, and quality assessment plots are user-defined in some places, only to be hard coded in others.  Additionally, processings scripts nominally executed after the basic processing performed by `run_pypeit` (for example, 2D coadding using `pypeit_coadd_2dspec`) make new directories that are never user defined (like `Science_coadd`).  Finally, we have been migrating away from use of `os.path` toward use of the `pathlib` library; however, this transition has happened slowly.
>
> The goal of this `paths` branch is to propose a new paradigm for how pypeit handles its output paths.
>
> First, create a new markdown file called `pypeit_output_paths.md` in a new directory in the dev-suite called `PypeIt-development-suite/claude_prompts/pypeit_output_paths`.  This will be our primary reference document for the project.  Start the document with a copy of this prompt.  The results from you effort related to the next 2 requests should be logged in this markdown file.
>
> Second, analyze the pypeit codebase to identify where paths are defined, where these definitions are adhered to, and where they are ignored.  In this process, consolidate a list of all output paths and the script or module that defines each one.
>
> Third, propose a new strategy for handling pypeit's output paths.  The strategy should try to limit the amount of book-keeping that is required of developers and provide a clear approach for maintaining paths across modules.  For example, an output paths class could be defined similar to the current `PypeItDataPaths` or `PypeItLogger` objects, which are instantiated in the package `__init__.py` and imported and modified as needed by subsequent functions.  A different approach could be to have classes that need access to output paths inherit from a common base class (similar to the `CalibFrame` class).  Importantly, output paths or file names *must* be passed directly as positional or keyword arguments for functions in the `core/` module; however, this is not a requirement of pypeit's higher-level classes (like `Calibrations`).

---

## 2. Current-state analysis

### 2.1 Where output paths are defined today

All output-directory *names* (as opposed to full paths) are declared as string
parameters in `pypeit/par/pypeitpar.py`. There is no separate "full path"
parameter for Science/QA — only a subdirectory name that gets joined onto a
single root, `redux_path`.

| Concept | Parameter | Default | Location |
|---|---|---|---|
| Reduction root | `ReduxPar.redux_path` | `os.getcwd()` | `pypeit/par/pypeitpar.py:2859-2862` |
| Science subdir name | `ReduxPar.scidir` | `'Science'` | `pypeit/par/pypeitpar.py:2850-2852` |
| QA subdir name | `ReduxPar.qadir` | `'QA'` | `pypeit/par/pypeitpar.py:2854-2857` |
| Calibrations subdir name | `CalibrationsPar.calib_dir` | `'Calibrations'` | `pypeit/par/pypeitpar.py:4803-4808` |
| Collate1D output dir | `Collate1DPar.outdir` | `os.getcwd()` | `pypeit/par/pypeitpar.py:~5660` |

`CalibrationsPar.calib_dir`'s docstring is explicit about the intended
contract and about its own fragility: *"The name of the directory for the
processed calibration frames. The host path for the directory is set by the
redux_path... Beware that success when changing the default value is not
well tested!"*

These parameters are parsed from a `.pypeit` file's `[rdx]`/`[calibrations]`/
`[collate1d]` blocks via each ParSet's `from_dict` (`parkeys` lists at
`pypeitpar.py:2889-2890`, `:4968`). The only command-line override is
`run_pypeit.py`'s `-r/--redux_path` flag (`pypeit/scripts/run_pypeit.py:48-49`),
which overwrites `self.par['rdx']['redux_path']` in `PypeIt.__init__`
(`pypeit/pypeit.py:87-88`) before any downstream path is computed. There is
**no CLI override** for `scidir`, `qadir`, or `calib_dir`.

### 2.2 How the "official" paths are resolved and threaded through

- `PypeIt.__init__` eagerly builds `self.calibrations_path` as a plain
  attribute: `Path(redux_path) / calib_dir` (`pypeit/pypeit.py:119`). This is
  logged and used immediately for calibration-association bookkeeping
  (`pypeit.py:122-136`).
- `science_path` and `qa_path` are instead **properties**, computed on
  demand rather than stored — and inconsistently typed:
  - `science_path` (`pypeit.py:144-147`) delegates to the one canonical
    helper function, `outputfiles.science_path(par)`
    (`pypeit/outputfiles.py:83-95`: `Path(par['rdx']['redux_path']) / par['rdx']['scidir']`),
    and returns a `pathlib.Path`.
  - `qa_path` (`pypeit.py:149-152`) instead inlines
    `os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['qadir'])`
    and returns a plain `str` — the same logic, duplicated, with a different
    return type.
- Directory creation is not centralized: the QA directory is created via
  `PypeIt.build_qa()` → `qa.gen_qa_dir(qa_path)`, which does
  `pathlib.Path(qa_path).mkdir(parents=True, exist_ok=True)`
  (`pypeit/qa.py:426-436`). The Science directory is **not** created in
  `pypeit.py` at all — it's created lazily right before a file is written, in
  `pypeit/core/exposure.py:502-505` (`save_exposure`) and
  `pypeit/pypeit_steps.py:306-309` (intermediate science images).
- `pypeit/calibrations.py`'s `Calibrations` class receives `caldir`/`qadir` as
  constructor arguments (`calibrations.py:144-146`) rather than computing them
  itself, but then **independently recomputes/duplicates** the same
  path-joining logic that `pypeit.py`/`outputfiles.py` already do, storing
  `self.calib_dir = Path(caldir).absolute()` (creating it if missing,
  `calibrations.py:170-172`) and `self.qa_path = Path(qadir).absolute()` plus
  a derived `PNGs` subdirectory (`calibrations.py:175-181`). This duplication
  is **self-documented as a known problem**:

  ```python
  # calibrations.py:177
  # TODO: This should only be defined in one place!  Where?...
  ```

  The caller (`pypeit/pypeit_steps.py:124-180`, `calib_one`) independently
  rebuilds the QA path string from raw pars
  (`pypeit_steps.py:162-163`: `os.path.join(par['rdx']['redux_path'], par['rdx']['qadir'])`)
  if one isn't explicitly passed — a third independent computation of the
  identical value.

**Summary of the core problem:** the same three concepts (redux root +
scidir/qadir/calib_dir) are joined into full paths in at least four different
places (`pypeit.py:119`, `pypeit.py:150-152`, `outputfiles.py:83-95`,
`calibrations.py:170-181`/`pypeit_steps.py:162-163`), with inconsistent
typing (`str` vs `Path`) and inconsistent eagerness (computed once vs.
recomputed as a property vs. recomputed by a downstream caller).

### 2.3 Where the official path system is hardcoded around or ignored

Nine categories of violation were identified:

**(1) `pypeit/coadd2d.py` — `CoAdd2D.output_paths()` (`coadd2d.py:475-513`),
used by `pypeit_coadd_2dspec`.** This is the clearest instance of the
"`Science_coadd` is never user-defined" problem named in the prompt:

```python
# coadd2d.py:502-513
if coadd_dir is not None:
    pypeit_scidir = Path(coadd_dir).absolute() / 'Science'   # (1a) hardcoded literal, and...
else:
    pypeit_scidir = Path(spec2d_files[0]).parent
coadd_scidir = pypeit_scidir.parent / f"{par['rdx']['scidir']}_coadd"  # (1b) ...immediately discarded via .parent
if not coadd_scidir.exists():
    coadd_scidir.mkdir(parents=True)                          # (1c) mkdir buried in coadd2d.py, not centralized
par['rdx']['qadir'] += '_coadd'                                # (1d) MUTATES THE CALLER'S par IN PLACE
qa_path = pypeit_scidir.parent / par['rdx']['qadir'] / 'PNGs'  # (1e) hardcoded 'PNGs' literal
if not qa_path.exists():
    qa_path.mkdir(parents=True)
return str(coadd_scidir), str(qa_path)
```

The `'_coadd'` suffix is completely non-configurable, and line 512's `+=`
permanently mutates the shared `par` object — any code that inspects
`par['rdx']['qadir']` after this call sees the mutated value, and calling
`output_paths()` twice on the same `par` would compound the suffix. Called
from `pypeit/scripts/coadd_2dspec.py:84-85`, which then hand-rolls the
`spec1d_`/`spec2d_` filenames itself
(`coadd_2dspec.py:229,233,238`, e.g. `coadd_scidir / f'spec1d_{basename}.fits'`)
rather than calling the canonical `outputfiles.spec_output_file()`
(`outputfiles.py:97-122`).

**(2) `pypeit/qa.py` — pervasive hardcoded `'QA'`/`'QA/PNGs'` literals.**
The worst offender, and self-documented in the code:

```python
# qa.py:25-26
# TODO: Move these names to the appropriate class.  This always writes
# to QA directory, even if the user sets something else...
```

- `set_qa_filename()` (`qa.py:27-139`) returns a literal string containing
  `'QA/PNGs/...'` or `'PNGs/...'` in every `match`/`case` branch, never
  referencing `par['rdx']['qadir']`. Its `out_dir` argument is only ever
  populated with `pathlib.Path.cwd()` or the literal `'QA'`
  (`qa.py:324`, `qa.py:402`), never with the user's configured `qadir`.
- `gen_exp_html()` (`qa.py:494-520`) takes **no path argument at all**:
  `pathlib.Path("QA") / "PNGs"` (line ~500) and `f"QA/{uni_name}.html"`
  (line ~506) are both hardcoded relative to the process's current working
  directory, independent of `redux_path`/`qadir` entirely.

**(3) `pypeit/extraction.py:584` and `pypeit/find_objects.py:1252`** —
both hardcode `os.path.join(self.par['rdx']['redux_path'], 'QA')`: they
correctly use `redux_path` but hardcode the `'QA'` subdirectory name instead
of `self.par['rdx']['qadir']`.

**(4) `pypeit/scripts/trace_edges.py`** — duplicates the
`'Calibrations'` default as a second, independent argparse literal
(`--calib_dir` default, `trace_edges.py:42`) rather than referencing the
`CalibrationsPar` default — the two could silently diverge. The same script
has **both** a correct path (`qa_path = rdx.qa_path`, `trace_edges.py:127`,
when a `--pypeit_file` is given) and a hardcoded one
(`qa_path = redux_path / 'QA'`, `trace_edges.py:154`, in the no-pypeit-file
branch, which also never loads a full `PypeItPar` and so has no `qadir` to
reference at all).

**(5) `pypeit/scripts/ql.py`** — the *same script* has one branch that
hand-rolls the spec2d filename (`ql.py:1024`:
`spec2d_file = str(coadd_scidir / f'spec2d_{basename}.fits')`, reusing
`coadd2d.py`'s hardcoded logic from (1)) and another branch, a few lines
below, that correctly calls the canonical helper
(`ql.py:1028`: `outputfiles.spec_output_file(pypeIt.fitstbl, pypeIt.par, frame, twod=True)`).

**(6) Several post-processing scripts write outputs relative to the process
cwd, entirely independent of `redux_path`:**
- `pypeit/scripts/sensfunc.py:85,216` — `--par_outfile` defaults to
  `'sensfunc.par'`; output filename is `'sens_' + spec1dname`.
- `pypeit/scripts/flux_calib.py:62` — `--par_outfile` defaults to
  `'fluxing.par'`.
- `pypeit/scripts/coadd_1dspec.py:144` and `:54-55`
  (`build_coadd_file_name`) — `--par_outfile` defaults to `'coadd1d.par'`;
  the coadd1d output file is placed "alongside the first input spec1d file"
  (`os.path.dirname(os.path.abspath(spec1dfiles[0]))`), independent of
  `redux_path`/`scidir` — and one of the few remaining `os.path`-based
  (non-pathlib) functions in this group.
- `pypeit/scripts/flux_setup.py:140-141,210-211,231-232` — `.flux`,
  `.coadd1d`, `.tell` files written with no directory join at all.

**(7) `pypeit/collate.py`/`pypeit/scripts/collate_1d.py`** — routes through
a wholly separate parameter namespace, `par['collate1d']['outdir']`
(own default `os.getcwd()`, `pypeitpar.py:~5660`), never integrated with
`par['rdx']['redux_path']`/`science_path`/`qa_path`. Also duplicates the
spec1d→spec2d filename translation via a hardcoded string-replace
(`collate.py:668`: `filename.replace('spec1d', 'spec2d', 1)`) rather than
using the shared `outputfiles.py` naming convention. Directory creation is
again buried inline (`collate.py:597`, `os.makedirs(outdir, exist_ok=True)`).

**(8) `pypeit/core/findobj_skymask.py:1453-1454`** — a function inside
`pypeit/core/` creates a directory as a side effect:

```python
qafile = Path(objfindQA_filename).absolute()
qafile.parent.mkdir(parents=True, exist_ok=True)
```

This is a mild violation of the "core/ functions should just receive paths,
not manipulate the filesystem" principle — the *path itself* is correctly
passed in as an argument (satisfying the hard constraint on `core/`
function signatures), but the `mkdir` call means `core/` is still taking
directory-management action rather than leaving it entirely to the caller.
By contrast, `pypeit/core/coadd.py` was checked and has no path-construction
or directory-management logic at all — callers pass fully-resolved
strings/paths in, matching the intended `core/` contract.

**(9) `os.path` → `pathlib` migration snapshot.** Of 247 non-test `.py`
files under `pypeit/`, 89 already use `pathlib`; 24 still use
`os.path.join`/`os.makedirs`/`os.getcwd` for path construction, including
(with call-site counts) `metadata.py` (6), `collate.py` (6),
`scripts/coadd_1dspec.py` (4), `io.py` (3), `gui/setup_gui/model.py` (3),
`archive.py` (3), plus single/low-count usages in `coadd3d.py`,
`extraction.py`, `find_objects.py`, `multislit_flexure.py`,
`par/parset.py`, `par/pypeitpar.py`, `pypeit_steps.py`,
`scripts/chk_for_calibs.py`, `scripts/chk_noise_1dspec.py`,
`scripts/chk_noise_2dspec.py`, `scripts/collate_1d.py`,
`scripts/install_ql_calibs.py`, `scripts/obslog.py`, `utils.py`. Notably,
the hardcoding problems in (1)-(7) are independent of which path API is
used — several of the worst offenders (`coadd2d.py`, `sensfunc.py`,
`flux_setup.py`) have already migrated to `pathlib.Path` but still hardcode
directory-name *literals*.

### 2.4 Existing architectural precedents in the codebase

Three patterns already exist in PypeIt that are directly relevant models for
a new output-paths system:

**`PypeItDataPaths`** (`pypeit/pkg/pypeitdata.py`) — a registry of *input*
data-file locations (package-shipped or cache-downloaded), instantiated once
at package-import time in `pypeit/__init__.py:22-25`:
```python
from .pkg.pypeitdata import PypeItDataPaths
dataPaths = PypeItDataPaths()
```
`PypeItDataPaths.__init__` loops over a class-level `defined_paths` dict and
`setattr`s one named `PypeItDataPath` object per registered path (e.g.
`dataPaths.nist`, `dataPaths.arclines`); each has its own `.path`
(`pathlib.Path`), `.get_file_path()`, `.glob()`. ~20 modules import it via
`from pypeit import dataPaths` and call e.g.
`dataPaths.nist.get_file_path('ThAr_vacuum.ascii')`.

**`PypeItLogger`/`log`** (`pypeit/pkg/logger.py`) — instantiated once at
package-import time (`pypeit/__init__.py:17-20`,
`log = get_logger(level=logging.DEBUG)`) with defaults, then **reconfigured
in place** once user CLI args are known:
`ScriptBase.init_log()`/`run_pypeit.py` call `log.init(level=..., log_file=...)`
on the *same* object (it's backed by Python's `logging.getLogger("pypeit")`
registry, so `get_logger()`/`.init()` always mutate the one named logger
rather than creating a new one). ~25+ modules import it via
`from pypeit import log`.

**`CalibFrame`** (`pypeit/calibframe.py`) — the base class that calibration
*data objects* (`WaveCalib`, `FlatImages`, `SlitTraceSet`, `EdgeTraceSet`,
`Alignments`, `ScatteredLight`, `WaveTilts`, and `PypeItCalibrationImage`,
via multiple inheritance — 8 subclasses total) inherit from to get a
standardized output-path/naming contract:
- Class attributes subclasses override: `calib_type` (filename prefix, e.g.
  `'WaveCalib'`), `calib_file_format` (extension, default `'fits'`).
- `set_paths(self, odir, setup, calib_id, detname)` — called **post-hoc**
  (after `__init__`, not passed to the constructor) by the orchestrating
  `Calibrations` object; sets `self.calib_dir` (absolute `Path`, mkdir'd if
  missing), `self.calib_id`, `self.calib_key` (deterministic string).
- `get_path(self)` — read accessor, `construct_file_name(self.calib_key, calib_dir=self.calib_dir)`.
- `construct_file_name(cls, calib_key, calib_dir=None)` (classmethod) —
  `f'{cls.calib_type}_{calib_key}.{cls.calib_file_format}'`, optionally
  joined onto `calib_dir`.
- `to_file(self, file_path=None, ...)` — defaults `file_path` to
  `self.get_path()` if not given, so the deterministic naming convention is
  the default write target.
- `glob(cls, calib_dir, setup, calib_id, detname=None)` (classmethod) —
  discovers existing files on disk matching the naming convention.

The orchestrating `Calibrations` class holds its own `self.calib_dir`
(computed from a constructor `caldir` argument) and calls
`some_calib_object.set_paths(self.calib_dir, setup, calib_id, detname)`
at many call sites (e.g. `calibrations.py:763,995,1896`).

None of these three precedents currently manage *output* directories for
Science/QA/coadd/collate/etc. — that gap is the subject of the two
candidate designs below.

---

## 3. Proposed strategies

Both designs below satisfy the hard constraint stated in the prompt:
**functions in `pypeit/core/` must receive paths/filenames as plain
positional or keyword arguments — never by querying a paths object.** Both
designs place the paths object/mixin strictly at the orchestration layer
(`PypeIt`, `Calibrations`, `CoAdd2D`, script `main()`s), which resolve paths
and pass plain `Path`/`str` primitives down into `core/`. Both designs leave
`CalibFrame`'s per-calibration-file naming convention conceptually
untouched — the real gap is in *directory-level* policy (Science/
Calibrations/QA roots and derived directories like the 2D-coadd output),
not in how individual calibration files are named within those directories.

No single strategy is recommended here — per project discussion, both are
presented in full for the team to weigh and decide between.

### 3.1 Design A — centralized `PypeItOutputPaths` singleton

Mirrors the existing `dataPaths`/`log` precedent directly: one object,
instantiated once, resolving every output-path concept, reconfigured in
place once the user's actual parameters are known.

#### Class structure

New module `pypeit/outputpaths.py` (sibling of `pypeit/outputfiles.py`,
whose pure name-construction helpers stay and are *called by* the new class
rather than duplicated):

```python
class PypeItOutputPaths:
    def __init__(self, redux_path=None, scidir='Science', qadir='QA',
                 calib_dir='Calibrations', *, coadd_suffix='_coadd',
                 collate_outdir=None):
        self.redux_path = Path(redux_path).absolute() if redux_path else Path.cwd()
        self._scidir, self._qadir, self._calib_dir = scidir, qadir, calib_dir
        self._coadd_suffix = coadd_suffix
        self._collate_outdir = Path(collate_outdir).absolute() if collate_outdir else None

    @property
    def science(self)      -> Path: return self.redux_path / self._scidir
    @property
    def qa(self)           -> Path: return self.redux_path / self._qadir
    @property
    def qa_pngs(self)      -> Path: return self.qa / 'PNGs'
    @property
    def calibrations(self) -> Path: return self.redux_path / self._calib_dir
    # derived post-processing outputs -- NO par mutation, suffix configurable
    @property
    def coadd_science(self)  -> Path: return self.redux_path / f'{self._scidir}{self._coadd_suffix}'
    @property
    def coadd_qa_pngs(self)  -> Path: return self.redux_path / f'{self._qadir}{self._coadd_suffix}' / 'PNGs'
    @property
    def collate(self)        -> Path: return self._collate_outdir or self.redux_path
```

Plus:
- **`configure(par, redux_path_override=None, collate_outdir=None)`** — the
  `log`-style in-place reconfigure. Reads `par['rdx']['redux_path'/'scidir'/'qadir']`,
  `par['calibrations']['calib_dir']`, `par['collate1d']['outdir']` and
  mutates the *same* instance, so an early `from pypeit import outputPaths`
  reference stays valid/live.
- **`make(*which)`** — centralized `mkdir(parents=True, exist_ok=True)` for
  named attributes (e.g. `outputPaths.make('science', 'qa_pngs', 'calibrations')`),
  replacing the scattered inline `mkdir` calls currently in `coadd2d.py`,
  `calibrations.py`, `core/exposure.py`, `pypeit_steps.py`, and
  `core/findobj_skymask.py`.
- **File-name convenience methods** delegate to `outputfiles.py` (e.g.
  `spec2d_file()` wraps `outputfiles.spec_output_file(..., sci_path=self.science)`)
  rather than duplicating naming logic.

The object is pathlib-native and does zero I/O at construction time — as
import-safe as `dataPaths`.

#### Lifecycle

**Hybrid: import-time singleton with cwd defaults, reconfigured in place at
each entry point.** This matches the `log` pattern and is the right fit
specifically because multiple independent CLI entry points
(`run_pypeit`, `pypeit_coadd_2dspec`, `pypeit_collate_1d`, `pypeit_sensfunc`, ...)
never share a live `PypeIt` object.

In `pypeit/__init__.py`, right after `log`/`dataPaths` are created:
```python
from .outputpaths import PypeItOutputPaths
outputPaths = PypeItOutputPaths()      # cwd-based defaults
```

- **Live reduction:** `PypeIt.__init__` calls
  `outputPaths.configure(self.par, redux_path_override=redux_path)`,
  replacing the eager `self.calibrations_path = ...` line. The
  `science_path`/`qa_path` properties become thin delegators returning
  `Path` consistently (fixing today's `str`-vs-`Path` inconsistency).
  `Calibrations(...)` is fed `outputPaths.calibrations`/`.qa` directly, so
  `Calibrations.__init__` stops recomputing QA/PNG logic — removing the
  `# TODO` at `calibrations.py:177`.
- **Derived/after-the-fact tools:** each script's `main()` builds its `par`,
  then calls `outputPaths.configure(par, redux_path_override=args.redux_path)`.
  `pypeit_coadd_2dspec` reads `outputPaths.coadd_science`/`.coadd_qa_pngs`;
  its `coadd_dir` CLI argument maps to `redux_path_override`.
- **Threading into `core/`:** orchestration classes pass resolved
  primitives (never the object itself) down into `core/` functions.

#### How this eliminates each violation

1. **`coadd2d.py output_paths()`** — entire body replaced by reads of
   `outputPaths.coadd_science`/`.coadd_qa_pngs` + `outputPaths.make(...)`.
   The `_coadd` suffix becomes the configurable `coadd_suffix` constructor
   arg. Par mutation is structurally impossible — the object copies scalars
   out of `par` at `configure()` time and never writes back to it.
2. **`qa.py`** — `set_qa_filename` stays pure but every branch is
   normalized to emit only `PNGs/<name>` (dropping the leaked `'QA/'`
   prefixes); callers pass `out_dir=outputPaths.qa` explicitly.
   `gen_exp_html` gains a `qa_path` parameter. Removes the `# TODO` at
   `qa.py:25`.
3. **`extraction.py`/`find_objects.py`** — replace the hardcoded
   `os.path.join(redux_path, 'QA')` with `outputPaths.qa`; the resulting
   filename is still passed down into `core/` as a plain string
   (constraint-safe).
4. **`trace_edges.py`** — drop the duplicate `'Calibrations'` argparse
   literal; both CLI branches read `outputPaths.calibrations`/`.qa` after
   `configure()` (the no-pypeit-file branch simply gets cwd defaults,
   which is the current behavior anyway).
5. **`ql.py`** — both branches call the same `outputPaths.spec2d_file(...)`
   helper (which wraps `outputfiles.spec_output_file`); the hand-rolled
   branch is deleted.
6. **`sensfunc.py`/`flux_calib.py`/`coadd_1dspec.py`/`flux_setup.py`** —
   `configure()` then build outputs under `outputPaths.redux_path` via a
   generic `aux_file(name)` helper; defaulting `redux_path` to cwd preserves
   today's behavior exactly while making it overridable.
7. **`collate.py`/`collate_1d.py`** — the existing `par['collate1d']['outdir']`
   folds into `outputPaths.collate` (falling back to `redux_path` when
   unset, so files with an explicit `outdir` keep working unchanged); the
   spec1d→spec2d string-replace is replaced by the shared naming helper.
8. **`core/findobj_skymask.py`** — the function already receives the
   filename as an argument; the `mkdir` call moves up to the
   `find_objects.py` caller (`outputPaths.make('qa_pngs')`), leaving
   `core/` with a purely passive path argument.
9. **os.path → pathlib stragglers** — the object is pathlib-native
   end-to-end; every call site touched by (1)-(8) sheds its
   `os.path.join`/`os.makedirs`/`os.getcwd` calls as a side effect of
   adoption.

#### Migration plan (phased, backward-compatible)

Invariant held throughout: existing `.pypeit` files with
`[rdx] redux_path/scidir/qadir`, `[calibrations] calib_dir`,
`[collate1d] outdir` keep working unchanged — none of those ParSet keys are
renamed, removed, or redefaulted.

- **Phase 0** — add the class + singleton in `__init__.py`; nothing
  consumes it yet. Unit-test its properties against today's
  literal-derived paths for both default and custom `scidir`/`qadir`.
- **Phase 1** — adopt in `pypeit.py` (replace `calibrations_path`
  assembly and the two path properties; thread into `Calibrations`;
  delete the duplicated QA/PNG derivation there).
- **Phase 2** — refactor `outputfiles.py`'s `science_path`/`spec_output_file`
  to accept the resolved path directly (temporary shims for any external
  callers), enabling script-level adoption.
- **Phase 3** — fix `qa.py` + `extraction.py`/`find_objects.py`/
  `core/findobj_skymask.py` (violations 2, 3, 8).
- **Phase 4** — fix `coadd2d.py` + `pypeit_coadd_2dspec` (violation 1);
  optionally add an additive `[rdx] coadd_suffix` ParSet key
  (default `'_coadd'`, via the `add-parameter` skill) so the suffix is
  user-documented and overridable, not just de-hardcoded internally.
- **Phase 5** — fix the remaining post-processing scripts (violations
  4-7); add a shared `-r/--redux_path` CLI option to scripts that lack one
  today, mapping to `configure(redux_path_override=...)`.
- **Phase 6** — remove Phase-2 shims; run the dev-suite 2D-coadd/collate/
  sensfunc/flux setups to confirm an identical output-directory tree to
  today's before/after.

#### Tradeoffs

Centralization is this design's strength: one object, sourced once, hands
every disconnected CLI entry point (the real pain point, since coadd/
collate/sensfunc run outside any live `PypeIt` reduction) a consistent,
mutation-safe `Path`. It maps directly onto the already-adopted `log`/
`dataPaths` idioms, so it is low-surprise to PypeIt developers. Its
weaknesses: it is global mutable state (tests must `configure()`/reset it
between cases; genuinely concurrent in-process reductions would contend,
though this is mitigated by allowing local, non-singleton instances where
needed), and it separates directory *policy* from the *datamodel* of the
product being written there — per-file naming stays a function of
`outputfiles.py`/`CalibFrame`, not of the paths object itself.

---

### 3.2 Design B — generalized `CalibFrame`-style mixin (`OutputPathsMixin`)

Generalizes `CalibFrame`'s existing `set_paths`/`get_path`/
`construct_file_name` contract into a standalone mixin that any orchestrator
class — not just calibration-frame `DataContainer`s — can adopt.

#### Class structure

`CalibFrame` today conflates two things: (1) owning a directory / naming
files deterministically / creating directories on demand, and (2) being a
`DataContainer` with calibration-specific key semantics
(`calib_type`, `calib_key`, `calib_id`, `setup`/`detname`). Design B splits
out the generic half:

```python
# new file: pypeit/outputpaths.py (no dependency on datamodel or core/)
class OutputPathsMixin:
    output_file_prefix = None      # e.g. 'spec2d', 'sens', 'coadd1d'; None => no prefix
    output_file_format = 'fits'    # extension without dot

    def set_output_dir(self, odir, mkdir=True):
        """Resolve, store, and (optionally) create the owned output dir."""
        self.output_dir = None if odir is None else Path(odir).absolute()
        if mkdir and self.output_dir is not None and not self.output_dir.exists():
            self.output_dir.mkdir(parents=True)
        return self.output_dir

    @classmethod
    def construct_file_name(cls, stem, output_dir=None, ext=None):
        _ext = cls.output_file_format if ext is None else ext
        name = stem if cls.output_file_prefix is None else f'{cls.output_file_prefix}_{stem}'
        name = f'{name}.{_ext}' if _ext else name
        return name if output_dir is None else Path(output_dir).absolute() / name

    def get_path(self, stem):
        return self.construct_file_name(stem, output_dir=self.output_dir)
```

`CalibFrame` is then re-expressed as a specialization, gaining the mixin via
multiple inheritance without changing its public API:

```python
class CalibFrame(OutputPathsMixin, datamodel.DataContainer):
    def set_paths(self, odir, setup, calib_id, detname):
        self.set_output_dir(odir)
        self.calib_dir = str(self.output_dir)   # preserve existing str type + header value
        self.calib_id  = self.ingest_calib_id(calib_id)
        self.calib_key = self.construct_calib_key(setup, self.calib_id, detname)

    @classmethod
    def construct_file_name(cls, calib_key, calib_dir=None):
        return super().construct_file_name(calib_key, output_dir=calib_dir, ext=cls.calib_file_format)
        # (output_file_prefix derives from calib_type via a small class-level adapter)
```

No existing `CalibFrame` subclass (`WaveCalib`, `FlatImages`, `SlitTraceSet`,
`EdgeTraceSet`, `Alignments`, `ScatteredLight`, `WaveTilts`,
`PypeItCalibrationImage`) needs to change — this is the design's
"must not regress" guarantee.

Non-`DataContainer` orchestrators adopt the mixin directly:
```python
class CoAdd2D(OutputPathsMixin):
    output_file_prefix = 'spec2d'

class QAPaths(OutputPathsMixin):
    output_file_prefix = None
    output_file_format = 'png'
```
Because the mixin has no metaclass and no `DataContainer` requirement,
multiple inheritance with plain `object`-rooted classes is clean.

#### Lifecycle

Follows the established `CalibFrame` convention: construct the object, then
have the *orchestrator* call `set_output_dir()` (or, for `CalibFrame`,
`set_paths()`) post-hoc rather than passing the directory to `__init__`
(matters because a reconstructed `WaveCalib` gets its directory from a FITS
header via `calib_keys_from_header`, not from a fresh construction call).

A small free function still centralizes "which directory does each owner
get" — the mixin manages *one* directory per instance, but something has to
decide which directory that is:
```python
def resolve_redux_paths(par):
    root = Path(par['rdx']['redux_path']).absolute()
    return dict(science = root / par['rdx']['scidir'],
                qa       = root / par['rdx']['qadir'],
                calib    = root / par['calibrations']['calib_dir'])
```
`PypeIt.__init__` calls this once and hands each owner its slice
(`self.set_output_dir(paths['science'])`, `Calibrations(..., paths['calib'])`,
`QAPaths().set_output_dir(paths['qa'])`).

**`CoAdd2D`'s standalone coadd-directory derivation** (the concrete fix for
the `Science_coadd`-never-user-defined complaint), with no par mutation:
```python
def set_coadd_paths(self, par, spec2d_files, coadd_dir=None, suffix='coadd'):
    base = Path(coadd_dir).absolute() if coadd_dir is not None \
           else Path(spec2d_files[0]).parent.parent
    self.set_output_dir(base / f"{par['rdx']['scidir']}_{suffix}")
    self.qa = QAPaths()
    self.qa.set_output_dir(base / f"{par['rdx']['qadir']}_{suffix}")   # no par mutation
    return self.output_dir, self.qa.output_dir
```
`suffix` defaults to today's hardcoded string (preserving behavior) but is
now an explicit, overridable argument — ideally promoted to a real
parameter for documentation/discoverability.

#### How this eliminates each violation

1. **`coadd2d.py`** — `CoAdd2D` gains the mixin; `output_paths()` is
   replaced by `set_coadd_paths()` above. `'Science'` literal →
   `par['rdx']['scidir']`; hardcoded `_coadd` → configurable `suffix`;
   `'PNGs'` → owned by `QAPaths`; the `par['rdx']['qadir'] += '_coadd'`
   mutation is deleted outright (the suffix is applied to a local `Path`,
   never written back). `mkdir` moves into `set_output_dir`.
2. **`qa.py`** — a new `QAPaths(OutputPathsMixin)`, owned by `PypeIt` (and
   by `Calibrations`, and by `CoAdd2D`), holds the resolved QA dir and
   derives the `PNGs` subdirectory exactly once (replacing the duplicated
   derivation in both `calibrations.py` and `coadd2d.py`). `set_qa_filename`
   stops baking in `'QA/PNGs'`; each branch emits only the relative
   stem+leaf, and `QAPaths.get_path(stem)` does the joining.
   `gen_exp_html`/`html_mf_pngs` receive a `QAPaths` instance instead of a
   hardcoded `'QA'` string.
3. **`extraction.py`/`find_objects.py`** — these are orchestration-layer
   classes already; they receive an injected `QAPaths` owner and call
   `.get_path(stem)` instead of hardcoding the join.
4. **`trace_edges.py`** — the argparse default should reference the
   `CalibrationsPar` default directly rather than restating `'Calibrations'`;
   both CLI branches route through `resolve_redux_paths` +
   `Calibrations`/`QAPaths.set_output_dir`, collapsing the
   correct-vs-hardcoded fork into one path.
5. **`ql.py`** — both branches call `construct_file_name`
   (the mixin's classmethod, or the equivalent `outputfiles` wrapper once
   migrated) as the single naming authority; the hand-rolled branch is
   deleted.
6. **`sensfunc.py`/`flux_calib.py`/`coadd_1dspec.py`/`flux_setup.py`** —
   each script's driver gains the mixin with an appropriate
   `output_file_prefix` (`'sens'`, `'coadd1d'`, etc.) and calls
   `set_output_dir(redux_path)` so outputs land under the reduction root
   instead of cwd (with `Path.cwd()` as an explicit, backward-compatible
   default for genuinely standalone use).
7. **`collate.py`/`collate_1d.py`** — the collate driver gains the mixin;
   the existing `outdir` parameter continues to be the *input* to
   `set_output_dir` (preserving backward compatibility), but joining/mkdir/
   naming now goes through the mixin, and the spec1d→spec2d string-replace
   is replaced by the shared naming helper.
8. **`core/findobj_skymask.py`** — per the hard constraint, `core/` must
   not manage directories. The `mkdir` moves out to the caller (the
   `find_objects` class, itself a mixin owner), which calls
   `set_output_dir(..., mkdir=True)` before invoking the `core/` function;
   the `core/` function keeps receiving the already-existing path as a
   plain argument, neither creating nor querying anything.
9. **os.path → pathlib stragglers** — the mixin returns `Path` uniformly
   and accepts `os.PathLike`/`str`; every site touched by adoption is
   converted to pathlib as a side effect.

#### Migration plan (phased, backward-compatible)

Same backward-compat anchors as Design A (existing `.pypeit`-file keys keep
working unchanged; every existing `CalibFrame` subclass keeps working with
zero edits to the subclasses themselves).

- **Phase 0** — introduce `pypeit/outputpaths.py` with `OutputPathsMixin` +
  `resolve_redux_paths`; no caller changes yet. Unit-test naming/joining
  parity against `CalibFrame.construct_file_name` for a couple of real
  calibration types.
- **Phase 1 (highest-risk, do early, behind tests)** — refactor
  `CalibFrame(OutputPathsMixin, DataContainer)`; reimplement
  `set_paths`/`construct_file_name`/`get_path` as thin adapters. Success
  criterion: byte-identical output filenames and FITS
  `CALIBDIR`/`CALIBKEY` header values for existing dev-suite calibration
  outputs (regression-test against real dev-suite reductions).
- **Phase 2** — route `PypeIt.__init__` and `Calibrations.__init__` through
  `resolve_redux_paths` + `set_output_dir`; introduce `QAPaths`; migrate
  `qa.py`'s `set_qa_filename`/`gen_exp_html`/`html_mf_pngs` to take a
  `QAPaths` (violations 2, 3).
- **Phase 3** — replace `coadd2d.py`'s `output_paths()` with
  `set_coadd_paths()` (violation 1); ideally promote the `suffix` to a real,
  documented `Coadd2DPar` parameter (default `'coadd'`, preserving today's
  behavior) via the `add-parameter` skill.
- **Phase 4** — migrate the remaining scripts and `collate.py`
  (violations 4-7) — mechanical, lowest risk.
- **Phase 5** — `core/` hygiene: move the `mkdir` out of
  `core/findobj_skymask.py` (violation 8); sweep remaining pathlib
  stragglers opportunistically (violation 9).

#### Tradeoffs

The mixin approach is *distributed*: each class owns exactly the one
directory it writes to, matching the mental model developers already have
from `CalibFrame`, and it requires no new global object threaded through
every call site — which keeps the `core/`-takes-paths-as-args constraint
easy to enforce, since the owner (not a shared singleton) computes and
passes the path down. Its weakness is that "which directory does each owner
get" is still decided across several `__init__`s; `resolve_redux_paths`
mitigates this but is a plain helper function, not an enforced single
object, so a global "list every output directory this run will touch" query
is not as trivial as it would be with Design A. It is also a comparatively
larger, riskier refactor for equivalent benefit, since it touches
`CalibFrame` itself — a class with 8 existing subclasses currently working
correctly in production.

---

### 3.3 Side-by-side summary

| | Design A: centralized singleton | Design B: generalized mixin |
|---|---|---|
| Precedent followed | `dataPaths` / `log` | `CalibFrame` |
| Ownership model | One object holds all directories | Each class owns its own directory |
| Risk to existing code | Low — `CalibFrame` untouched | Higher — `CalibFrame` itself refactored (though subclasses unaffected) |
| Global path introspection | Trivial (one object to query) | Harder (must aggregate across owners) |
| Fits multi-entry-point CLI scripts | Directly — `configure()` at each entry point | Indirectly — `resolve_redux_paths()` called at each entry point |
| State model | Global mutable singleton (tests must reset/`configure`) | Distributed instance state (no shared mutable global) |
| `core/` constraint | Satisfied — object lives at orchestration layer only | Satisfied — mixin lives at orchestration layer only |
| `par` mutation bug (coadd2d) | Eliminated (object copies scalars out of `par`) | Eliminated (suffix applied to local `Path`, never written back) |

Both designs fully resolve all 9 catalogued violations and fully satisfy the
hard constraint that `core/` functions must receive paths as arguments.
The decision between them is primarily about state-management philosophy
(one global registry vs. distributed per-class ownership) and appetite for
refactoring `CalibFrame`, not about capability — left open for team
discussion.
