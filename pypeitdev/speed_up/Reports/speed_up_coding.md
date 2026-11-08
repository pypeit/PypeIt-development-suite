# Speeding up PypeIt — Implementation / Coding Document

**Version:** 0.1
**Date:** 2026-09-08
**Author:** JXP and Claude
**Implements:** `speed_up_design.md` **v0.5 (locked)**, §8 Implementation roadmap
**Target branch stack:** `speed_up` (PypeIt repo)

**Changelog**
- 0.1 (2026-09-08): First draft. Concrete implementation guide for the three
  stacked PRs A → B → C, grounded in the `speed_up` branch as of
  `e9ed85c1a` ("Merge pull request #2153 from pypeit/core_refactor").

---

## 0. Scope, conventions, and how to read this document

This document turns the locked design into code. It does **not** revisit design
decisions; where an implementation choice was not fixed by the design, this
document makes a recommendation and records the residual decision in
[§8 Open questions](#8-open-questions).

**Everything below refers to the `speed_up` branch.** The working tree may be on
another branch; inspect `speed_up` non-destructively with:

```bash
cd /mnt/tank/Astronomy/PypeIt/PypeIt
git show speed_up:pypeit/exposure.py | less
git show speed_up:pypeit/par/pypeitpar.py | sed -n '2758,2900p'
git ls-tree --name-only speed_up:pypeit
```

Every line number quoted here is from `speed_up` at `e9ed85c1a`. Re-verify them
before editing (`git show speed_up:<path> | sed -n 'A,Bp'`); they will drift as
the PRs land on top of each other.

### 0.1 Branch stack

| PR | Branch | Branches from | PR target | Contents |
|----|--------|---------------|-----------|----------|
| **A** | `speed_up_qa` | `speed_up` | `speed_up` | `Agg` backend; `par['rdx']['ncpu']` + `run_pypeit --ncpu`; parallel QA PNG writes |
| **B** | `speed_up_detpar` | `speed_up_qa` | `speed_up_qa` (retarget to `speed_up` once A merges) | `pypeit/parallel.py`; stages 1/2/4 of `reduce_exposure` parallelized; failed-detector bug fix |
| **C** | `speed_up_arcfit` | `speed_up_detpar` | `speed_up_detpar` (retarget as above) | Vectorized `fit_arcspec` |

`speed_up` itself eventually merges to `develop` (the B-spline work lands on a
separate branch and is explicitly out of scope here).

Each PR is independently revertible; each has its own kill switch (a parameter
default, a module constant, or a one-line backend call) — see
[§7 Risks and revert plan](#7-risks-and-revert-plan).

### 0.2 Code map (grounding)

The reduction driver on this branch is already factored into module-level
per-detector functions, which is what makes the parallelism cheap:

| File | Symbol | Lines | Role |
|---|---|---:|---|
| `pypeit/pypeit.py` | `PypeIt.__init__` | 60–142 | builds `self.par`, dumps `*_UTC_*.par` at 92–94 |
| | `PypeIt.build_qa` | 154–163 | QA HTML wrappers |
| | `PypeIt.calib_all` | 165–194 | calibration-only loop over detectors (186–190) |
| | `PypeIt.reduce_all` | 196–235 | standard loop (213–220), science loop (226–232) |
| | `reduce_calibID` | 250–384 | per-comb_id loop; calls `reduce_exposure` at 360–365, `save_exposure` at 370–374 |
| `pypeit/exposure.py` | `adjust_for_slitmask` | 21–80 | **cross-detector barrier** (stage 3) |
| | `process_exposure` | 84–138 | **stage 2**; flat detector loop at 125–135 |
| | `findobj_on_exposure` | 140–264 | **stage 3**; fan-out 202–222, barrier 226–231, fan-out 235–261 |
| | `extract_exposure` | 266–348 | **stage 4**; flat detector loop at 313–345 |
| | `reduce_exposure` | 350–466 | driver; detectors at 415–417, **stage 1** calib loop at 421–433 |
| | `save_exposure` | 468–565 | writes `spec1d`/`spec2d` |
| `pypeit/pypeit_steps.py` | `calib_one` | 124–199 | stage-1 body; returns a `Calibrations` object |
| | `process_one_det` | 202–319 | stage-2 body; returns `(sciImg, bkg_redux_sciimg)` |
| | `findobj_on_det` | 321–401 | stage-3 body |
| | `finalize_sky_det` | 402–481 | stage-3 body (post-barrier) |
| | `load_calibrations_for_frame` | 482–536 | reloads calibrations from disk (used by stages 2 & 4) |
| | `extract_det` | 629–801 | stage-4 body; `spec2DObj.gen_qa()` at 798 |
| `pypeit/par/pypeitpar.py` | `ReduxPar` | 2758–2900 | `__init__` 2770–2870, `from_dict` 2872–2888, `validate` 2890–2899 |
| `pypeit/scripts/run_pypeit.py` | `RunPypeIt.get_parser` | 38–65 | CLI |
| | `RunPypeIt.main` | 67–103 | `init_log` 85, `PypeIt(...)` 88–91, `build_qa()` 101 |
| `pypeit/qa.py` | module imports | 6–19 | `from matplotlib import pyplot as plt` at 13 |
| | `set_qa_filename` | 27–141 | QA filename convention |
| | `arc_tilts_2d_qa` / `_spec_` / `_spat_` | 544 / 610 / 693 | 3 PNGs **per slit** (savefig at 600, 684, 750) |
| | `spec_flexure_qa` | 762–918 | savefig 843, 914 |
| | `spat_flexure_qa` | 967–1099 | savefig 1098 |
| `pypeit/core/wavecal/autoid.py` | `arc_fit_qa` | 43–205 | 1 PNG per slit; `plt.close('all')` at 67 **and 200** |
| | `arc_fwhm_qa` | 207–~245 | 1 PNG per slit; `plt.close('all')` at 242 |
| `pypeit/flatfield.py` | fine-corr QA call | 1544–1554 | 1 PNG per slit |
| | `spatillum_finecorr_qa` | ~2000–2062 | savefig 2058 |
| `pypeit/core/findobj_skymask.py` | obj trace/profile QA | 1455, 1620 | 1 PNG per object |
| `pypeit/core/arc.py` | `detect_lines` | 857–1074 | calls `fit_arcspec` at 1028 |
| | `fit_arcspec` | 1077–1143 | **the per-line `curve_fit` loop (1121–1142)** |
| `pypeit/core/fitting.py` | `fit_gauss` | 874–918 | wraps `scipy.optimize.curve_fit` |
| | `gauss_3deg` / `gauss_4deg` / `guess_gauss` | 921 / 936 / 952 | model + initial guess |
| `pypeit/pkg/logger.py` | `PypeItLogger`, `get_logger` | 149–411 | plain `logging` subclass; `FileFormatter` at 138 |
| `pypeit/__init__.py` | `log = get_logger(...)` | 18–20 | the single package-level logger, `logging.getLogger("pypeit")` |

Two facts that shape the whole design:

1. **PypeIt uses stock `logging`.** `pypeit.log` is a `logging.Logger` subclass
   (`pypeit/pkg/logger.py:149`) obtained via `logging.getLogger("pypeit")`
   (`pypeit/pkg/logger.py:397`). Buffered per-worker logging is therefore a
   handler swap, not a bespoke mechanism.
2. **There is no multiprocessing anywhere in `pypeit/` today.** `git grep -n
   "multiprocessing\|ProcessPool\|ThreadPool\|ncpu\|nproc" speed_up -- 'pypeit/'`
   returns nothing. PR A introduces the first concurrency in the package.

### 0.3 Profile numbers used as targets

From `Reports/keck_deimos_600zd_m_6500_profile_report.md` (total wall
**12 390 s**, 4 mosaics) and `Reports/shane_kast_blue_A_profile_report.md`
(total wall **119 s**, 1 detector):

| Target | DEIMOS | Kast |
|---|---:|---:|
| Calibration build (stage 1: flat 18% + edges 15% + tilts 6% + wave 6% + combine 2%) | ~5 800 s (47%) | ~30 s (25%) |
| Local sky + extract (stage 4, `local_skysub_extract` cum) | 2 952 s (24%) | 40.5 s (34%) |
| Image proc / mosaic (stage 2, `build_mosaic` cum 678 s + `rawimage.process`) | ~700–900 s (6–7%) | 7.2 s (6%) |
| Object find + global sky (**stage 3, stays serial**) | ~2 400 s (19%) | ~25 s (21%) |
| QA (matplotlib text metrics cum 894 s + PNG encode 276 s, 971 PNGs) | ~0.9–1.2 ks | ~12 s, 20 PNGs |
| Arc-line Gaussian fits (`fit_gauss` cum 840 s, 1.11 M `curve_fit`, 93 M `gauss_3deg`) | ~840 s (7%) | ~2 s |

---

## 1. PR A — QA cheap wins (+ the `ncpu` plumbing)

**Branch:** `speed_up_qa` → **target:** `speed_up`

Three independent changes, deliberately shipped together because A is where the
`ncpu` control lands for the rest of the stack.

### A.1 New parameter `par['rdx']['ncpu']`

**File:** `pypeit/par/pypeitpar.py`, class `ReduxPar` (line 2758).

Three edits, following the existing conventions in that class exactly.

**(1)** Add the keyword to the signature (line 2770–2772):

```python
    def __init__(self, spectrograph=None, detnum=None, sortroot=None, calwin=None, scidir=None,
                 qadir=None, redux_path=None, ignore_bad_headers=None, slitspatnum=None,
                 maskIDs=None, quicklook=None, chk_version=None, ncpu=None):
```

**(2)** Add the specification block, e.g. immediately after the `quicklook`
block (which ends at line 2802):

```python
        defaults['ncpu'] = 1
        dtypes['ncpu'] = int
        descr['ncpu'] = 'Number of CPUs (worker processes) PypeIt may use to reduce ' \
                        'detectors/mosaics concurrently, and the number of threads used ' \
                        'to write QA figures.  The default, 1, runs the code fully ' \
                        'serially, exactly as in previous versions.  Values greater than ' \
                        '1 are capped at the number of detectors being reduced and at ' \
                        '``os.cpu_count()-1``.  Beware that peak memory usage scales ' \
                        'roughly linearly with the number of detectors reduced at the ' \
                        'same time.  Can be overridden on the command line with ' \
                        '``run_pypeit --ncpu``.'
```

**(3)** Add `'ncpu'` to `parkeys` in `from_dict` (line 2877–2878):

```python
        parkeys = [ 'spectrograph', 'quicklook', 'detnum', 'sortroot', 'calwin', 'scidir', 'qadir',
                    'redux_path', 'ignore_bad_headers', 'slitspatnum', 'maskIDs', 'chk_version',
                    'ncpu']
```

**(4)** Add a bound check to `validate` (line 2890):

```python
        if self.data['ncpu'] is not None and self.data['ncpu'] < 1:
            raise ValueError('ncpu must be a positive integer.')
```

> Follow the repo's `add-parameter` skill: after this change the parameter tables
> must be regenerated (`doc/scripts/build_par_rst.py`, driven by
> `cd doc; make html` or `./update_docs`), which rewrites `doc/pypeit_par.rst`.

### A.2 CLI flag `run_pypeit --ncpu N` (CLI overrides the par value)

**File:** `pypeit/scripts/run_pypeit.py`.

In `get_parser` (after the `-c/--calib_only` argument at line 62–63):

```python
        parser.add_argument('--ncpu', type=int, default=None,
                            help='Number of CPUs to use.  Overrides the [rdx] ncpu '
                                 'parameter.  The default (None) uses the parameter value, '
                                 'which itself defaults to 1 (fully serial).  Values >1 '
                                 'reduce detectors/mosaics concurrently and increase peak '
                                 'memory usage roughly in proportion.')
```

In `main` (line 88–91), pass it through:

```python
        pypeIt = pypeit.PypeIt(
            args.pypeit_file, reuse_calibs=args.reuse_calibs, overwrite=args.overwrite,
            redux_path=args.redux_path, calib_only=args.calib_only, show=args.show,
            ncpu=args.ncpu
        )
```

**File:** `pypeit/pypeit.py`, `PypeIt.__init__` (line 60–63). Add `ncpu=None` to
the signature and apply the override **next to the existing `redux_path`
override at lines 87–88**, i.e. *before* the `.par` file is dumped at 92–94, so
the recorded parameter file reflects what actually ran:

```python
        # Check the output paths are ready
        if redux_path is not None:
            self.par['rdx']['redux_path'] = redux_path
        # The command-line --ncpu overrides the parameter file
        if ncpu is not None:
            self.par['rdx']['ncpu'] = ncpu
```

Also document `ncpu` in the `PypeIt` class docstring `Args:` block (lines 36–51).

### A.3 Force the `Agg` backend

Both profiles record the interactive `qtagg` backend being used to write
non-interactive PNGs (Kast report §4; DEIMOS report §4). `pypeit/qa.py:13`
imports `pyplot` at module scope, and so do ~15 other reduction modules, so the
backend must be selected **before `pypeit.pypeit` is imported** — i.e. at the
top of `RunPypeIt.main`, not in `pypeit/__init__.py` (which would also break the
interactive GUIs in `pypeit/gui/`).

**File:** `pypeit/scripts/run_pypeit.py`, `main` (insert before
`from pypeit import pypeit` at line 73):

```python
        # All QA output from a reduction is non-interactive PNGs.  Force the
        # headless Agg backend, which avoids the GUI toolkit initialisation and
        # the text-metric overhead of the interactive default (qtagg).  This has
        # to happen before any pypeit module imports pyplot.  The one exception
        # is `-s/--show`, which deliberately raises blocking plots.
        import matplotlib
        if not args.show:
            matplotlib.use('Agg', force=True)
```

> The `not args.show` guard is a small deviation from the design's
> "unconditionally"; see **Q1** in §8. It introduces no new parameter.

### A.4 Deferred, optionally threaded QA figure writes

`plt.savefig` is where the expensive work happens (font/text metrics + Agg
render + PNG encode). The plan is therefore: **build figures on the main thread
as today, hand only `Figure.savefig` to a worker thread, and touch `pyplot`
global state (`plt.close`) only on the main thread.**

**File:** `pypeit/qa.py` — add at module scope, after the imports (line 19):

```python
from concurrent.futures import ThreadPoolExecutor

# --------------------------------------------------------------------------
# Deferred QA figure writing
#
# Rendering and encoding a QA PNG (matplotlib text metrics + Agg draw + PIL
# encode) dominates the QA cost and releases the GIL for most of its duration.
# When ncpu>1 the savefig call is handed to a small thread pool; the Figure
# object is only ever *created* and *closed* on the main thread, so pyplot's
# global state is never mutated concurrently.
# --------------------------------------------------------------------------

_QA_POOL = None
"""ThreadPoolExecutor used to write QA figures, or None for serial writes."""

_QA_PENDING = []
"""List of (future, figure, close) tuples that have not yet been reaped."""

_QA_MAX_PENDING = 16
"""Maximum number of un-reaped figures; bounds the memory held by open figures."""


def init_qa_pool(ncpu:int=1):
    """
    (Re)initialise the QA figure-writing thread pool.

    Call once per process, after the parameters are final.  Calling with
    ``ncpu<=1`` restores fully serial, in-line figure writing.  Also used to
    *reset* the pool inside a forked worker process, where the parent's threads
    do not exist.

    Parameters
    ----------
    ncpu
        Number of QA writer threads.  <=1 disables the pool.
    """
    global _QA_POOL, _QA_PENDING
    old = _QA_POOL
    _QA_POOL = None
    _QA_PENDING = []
    if old is not None:
        old.shutdown(wait=False)
    if ncpu is not None and ncpu > 1:
        _QA_POOL = ThreadPoolExecutor(max_workers=min(int(ncpu), 8),
                                      thread_name_prefix='pypeit-qa')


def save_figure(fig, outfile, show:bool=False, close:bool=True, **kwargs):
    """
    Write a matplotlib figure to disk, deferring the write to a background
    thread when the QA pool is active.

    Parameters
    ----------
    fig : `matplotlib.figure.Figure`_
        Figure to write.  Must not be modified after this call.
    outfile : :obj:`str`, `Path`
        Output file.  If None, nothing is written.
    show
        Show the figure interactively.  Forces the synchronous path.
    close
        Close the figure once it has been written.
    **kwargs
        Passed to `matplotlib.figure.Figure.savefig`_ (e.g. ``dpi``).
    """
    if show or _QA_POOL is None:
        if outfile is not None:
            fig.savefig(outfile, **kwargs)
        if show:
            plt.show()
        if close:
            plt.close(fig)
        return
    if outfile is None:
        if close:
            plt.close(fig)
        return
    _QA_PENDING.append((_QA_POOL.submit(fig.savefig, outfile, **kwargs), fig, close))
    if len(_QA_PENDING) >= _QA_MAX_PENDING:
        flush_qa()


def flush_qa():
    """
    Block until every deferred QA figure has been written, then close them.

    Exceptions raised in the writer threads are re-raised here, on the main
    thread.  Safe to call when the pool is inactive (it is then a no-op).
    """
    global _QA_PENDING
    pending, _QA_PENDING = _QA_PENDING, []
    for future, fig, close in pending:
        future.result()
        if close:
            plt.close(fig)
```

**Where the pool is created and drained**

- Create it in `PypeIt.__init__` (`pypeit/pypeit.py`), immediately after the
  `ncpu` override in §A.2: `qa.init_qa_pool(self.par['rdx']['ncpu'])`
  (`qa` is already imported at `pypeit/pypeit.py:18`).
- Drain at the end of `PypeIt.calib_all` (before `self.print_end_time()` at line
  194) and at the end of `PypeIt.reduce_all` (before line 235):
  `qa.flush_qa()`.
- Drain in `RunPypeIt.main` before `pypeIt.build_qa()` (line 101), so the HTML
  wrapper never races the PNGs: `from pypeit import qa; qa.flush_qa()`.
- Drain at the end of `exposure.reduce_exposure` (before the return at line 466)
  and at the end of `pypeit_steps.calib_one` (before the return at line 199), so
  no more than one detector's worth of figures is ever in flight.

### A.5 Converting the call sites

Convert the **high-volume, per-slit** writers only. Everything else (coadd,
sensfunc, telluric, the `PdfPages` writers) keeps `plt.savefig` and is untouched.

Priority order — these account for essentially all 971 DEIMOS PNGs:

| # | File / function | savefig line(s) | PNGs |
|---|---|---:|---|
| 1 | `pypeit/qa.py:arc_tilts_2d_qa` / `arc_tilts_spec_qa` / `arc_tilts_spat_qa` | 600, 684, 750 | **3 per slit** |
| 2 | `pypeit/core/wavecal/autoid.py:arc_fit_qa`, `arc_fwhm_qa` | 152, 199, ~242 | 2 per slit |
| 3 | `pypeit/flatfield.py:spatillum_finecorr_qa` (+ the second writer at 2126) | 2058, 2126 | 1 per slit |
| 4 | `pypeit/core/findobj_skymask.py` object trace/profile QA | 1455, 1620 | 1 per object |
| 5 | `pypeit/qa.py:spec_flexure_qa`, `spat_flexure_qa` | 843, 914, 1098 | per exposure |

The mechanical transformation, using `arc_tilts_2d_qa` (lines 596–606) as the
worked example:

```python
    # BEFORE (pypeit/qa.py:596-606)
    if outfile is not None:
        plt.savefig(outfile, dpi=400)

    if show_QA:
        plt.show()

    plt.close()
    plt.rcdefaults()

    # AFTER
    save_figure(fig, outfile, show=show_QA, dpi=400)
    plt.rcdefaults()
```

Three rules that MUST be observed at every converted site:

1. **The `Figure` object must be captured.** Sites that use
   `fig, ax = plt.subplots(...)` already have it (e.g. `qa.py:578`); sites that
   use `plt.figure()` + `plt.subplot(gs[...])` (e.g. `autoid.arc_fit_qa`,
   `flatfield.spatillum_finecorr_qa`) must be changed to bind
   `fig = plt.figure(...)` and pass that object.
2. **`plt.close('all')` must become `plt.close(fig)`.** `autoid.py:67`,
   `autoid.py:200` and `autoid.py:242` currently call `plt.close('all')`, which
   would destroy figures still queued for writing. This is the single most
   likely source of a silent corruption bug in this PR.
3. **`plt.tight_layout()` / `plt.subplots_adjust()` / `plt.rcdefaults()` must
   run before `save_figure`.** Layout is applied to the current figure; do not
   leave layout calls after the hand-off. `plt.rcdefaults()` after the hand-off
   is acceptable (artists capture rcParams at creation time) but must not be
   relied on for anything savefig-time; pass `dpi` explicitly, as the existing
   code already does.

### A.6 Tests (PR A)

**`pypeit/tests/test_pypeitpar.py`** — extend:

```python
def test_ncpu_default_and_override():
    par = pypeitpar.ReduxPar()
    assert par['ncpu'] == 1, 'ncpu must default to 1 (fully serial)'
    par = pypeitpar.ReduxPar.from_dict({'spectrograph': 'shane_kast_blue', 'ncpu': 4})
    assert par['ncpu'] == 4
    with pytest.raises(ValueError):
        pypeitpar.ReduxPar(ncpu=0)
```

**`pypeit/tests/test_qa.py`** — new tests, no large data needed:

```python
def test_save_figure_serial(tmp_path):
    qa.init_qa_pool(1)
    fig, ax = plt.subplots(); ax.plot([0, 1], [0, 1])
    out = tmp_path / 'serial.png'
    qa.save_figure(fig, out, dpi=50)
    assert out.exists()

def test_save_figure_threaded_matches_serial(tmp_path):
    """Deferred writes must be byte-for-byte identical to in-line writes."""
    import numpy as np
    from PIL import Image
    rng = np.random.default_rng(2718)          # fixed seed: deterministic
    data = rng.normal(size=(64, 64))

    def _draw():
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.imshow(data); ax.set_title('QA')
        return fig

    qa.init_qa_pool(1)
    ref = tmp_path / 'ref.png'
    qa.save_figure(_draw(), ref, dpi=80); qa.flush_qa()

    qa.init_qa_pool(4)
    outs = []
    for i in range(8):
        o = tmp_path / f'par{i}.png'
        qa.save_figure(_draw(), o, dpi=80)
        outs.append(o)
    qa.flush_qa()
    qa.init_qa_pool(1)

    ref_px = np.asarray(Image.open(ref))
    for o in outs:
        assert np.array_equal(np.asarray(Image.open(o)), ref_px)

def test_flush_qa_reraises(tmp_path):
    """A failed background write must surface on the main thread."""
    qa.init_qa_pool(2)
    fig, _ = plt.subplots()
    qa.save_figure(fig, tmp_path / 'nonexistent_dir' / 'x.png', dpi=50)
    with pytest.raises(Exception):
        qa.flush_qa()
    qa.init_qa_pool(1)
```

Compare **decoded pixel arrays**, not raw bytes (matplotlib embeds a `Software`
PNG chunk carrying the matplotlib version). `Pillow` is already available
transitively through matplotlib.

Every test must restore `qa.init_qa_pool(1)` before returning — the pool is
process-global state.

### A.7 Validation and re-profiling (PR A)

```bash
cd $PYPEIT_DEV/pypeitdev/speed_up/scripts
python profile_kast_blue.py    # ~2 min
python analyze_profile.py
python profile_deimos.py       # ~3.5 h
python analyze_profile.py --stem keck_deimos_600zd_m_6500
```

Add an `--ncpu N` argument to both profile scripts (they build the argv at
`profile_deimos.py:85` / the Kast equivalent) and suffix the artifact stems with
`_ncpu{N}` so the baseline is not overwritten.

Expected, from the profile numbers:

- **Agg alone:** removes the GUI-toolkit initialisation (`_tkinter` etc., 1.5 s
  on Kast) and part of the text-metric cost. Expect **1–4 s off Kast (1–3%)**
  and **100–400 s off DEIMOS (1–3%)**.
- **Agg + threaded PNG writes at `--ncpu 4`:** target a further **2–3× on the
  PNG encode/render budget** — i.e. up to another **~300–500 s on DEIMOS**,
  **~4–6 s on Kast**.
- **Combined PR A target: DEIMOS 12 390 s → ~11 500–11 900 s (4–7%); Kast 119 s
  → ~108–114 s (4–9%).**

**Measurement caveat:** matplotlib's FreeType text-metric path may not release
the GIL, in which case the thread pool gives little. Measure before claiming the
win; see **Q2** in §8 for the fallback.

### A.8 Documentation (PR A)

- Regenerate `doc/pypeit_par.rst` (`./update_docs`, or `cd doc; make html`).
- Regenerate `doc/help/run_pypeit.rst` (`doc/scripts/write_script_help.py`, run
  by the doc build).
- `doc/running.rst`: new subsection **"Running on multiple CPUs"** introducing
  `--ncpu` / `[rdx] ncpu`. In PR A state honestly that it currently only affects
  QA figure writing; PR B extends it.
- `doc/qa.rst`: note that reductions force the `Agg` backend and that QA PNGs
  are written concurrently when `ncpu>1`.
- `doc/releases/2.1.0dev.rst`, under **Functionality/Performance Improvements
  and Additions**:

  ```rst
  - ``run_pypeit`` now forces the headless matplotlib ``Agg`` backend (unless
    ``-s/--show`` is used), removing GUI-toolkit overhead from QA figure
    generation.
  - Added the ``[rdx] ncpu`` parameter and the matching ``run_pypeit --ncpu``
    command-line option (the latter takes precedence).  The default, 1, runs
    PypeIt exactly as before.  Values greater than 1 currently write QA figures
    concurrently.
  ```

  and under **Testing**: a bullet for the new QA tests.
- `CHANGES.rst` is marked deprecated by its own header ("USE OF THIS FILE IS NOW
  DEPRECATED"); the `update-changelog` skill still lists it. Add a one-line
  summary there only if the maintainers want it — otherwise the release doc
  alone.

---

## 2. PR B — Detector parallelism v1

**Branch:** `speed_up_detpar` → **target:** `speed_up_qa`

Parallelize stages **1 (calibrations)**, **2 (`process_exposure`)** and
**4 (`extract_exposure`)**. Stage 3 (`findobj_on_exposure`, with the
`adjust_for_slitmask` barrier at `exposure.py:228`) **stays serial**.

### B.1 New module `pypeit/parallel.py`

A single, small, self-contained helper. Nothing else in the package knows about
`concurrent.futures`.

```python
"""
Helpers for distributing per-detector work over multiple processes.

The reduction is already factored into module-level functions that operate on a
single detector/mosaic and return a per-detector result
(:func:`~pypeit.pypeit_steps.calib_one`,
:func:`~pypeit.pypeit_steps.process_one_det`,
:func:`~pypeit.pypeit_steps.extract_det`).  :func:`map_over_detectors` maps such
a function over a list of detectors, either serially (``ncpu<=1``; the literal
pre-existing code path) or over a ``fork``-based process pool (``ncpu>1``).

.. include:: ../include/links.rst
"""
import logging
import multiprocessing as mp
import os
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

from pypeit import log
from pypeit import PypeItError

__all__ = ['nworkers', 'map_over_detectors', 'surviving_detectors']


THREAD_ENV_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                   'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')
"""
Environment variables honoured by the BLAS/OpenMP runtimes numpy and scipy may
be linked against.  Pinned to 1 in each worker so that (a) N workers do not each
spawn N threads, and (b) the floating-point reduction order inside a worker
matches the serial run bit-for-bit.
"""

_FORK_PAYLOAD = None
"""
Module-global set by :func:`map_over_detectors` immediately *before* the pool is
created.  Because the workers are forked, they inherit this copy-on-write and
only the (tiny) detector key has to be pickled on submit.  Never read outside a
worker.
"""


def nworkers(ncpu, nitem):
    """
    Number of worker processes to use.

    Parameters
    ----------
    ncpu : :obj:`int`
        Requested number of CPUs (``par['rdx']['ncpu']``).
    nitem : :obj:`int`
        Number of independent work items (detectors).

    Returns
    -------
    :obj:`int`
        ``min(ncpu, nitem, os.cpu_count()-1)``, at least 1.
    """
    if ncpu is None or ncpu <= 1:
        return 1
    return max(1, min(int(ncpu), int(nitem), max(1, (os.cpu_count() or 1) - 1)))


class _RecordBuffer(logging.Handler):
    """
    Logging handler that accumulates picklable copies of the emitted records.

    The record message is rendered eagerly and ``args``/``exc_info`` are dropped
    (the approach used by `logging.handlers.QueueHandler`_), so the buffer can be
    shipped back to the parent process.  ``record.created`` is preserved, so the
    replayed lines carry the worker's real timestamps.
    """
    def __init__(self):
        super().__init__(level=logging.DEBUG)
        self.records = []

    def emit(self, record):
        rec = logging.makeLogRecord(record.__dict__)
        rec.msg = record.getMessage()
        rec.args = None
        rec.exc_info = None
        rec.exc_text = None
        self.records.append(rec)


def _worker_init():
    """
    Run once in each forked worker process.

    Pins the BLAS/OpenMP thread count, detaches the log handlers inherited from
    the parent (they point at the parent's log file and stderr), and disables the
    inherited QA thread pool (its threads do not survive the fork).
    """
    from pypeit import qa
    for var in THREAD_ENV_VARS:
        os.environ[var] = '1'
    for handler in log.handlers[:]:
        log.removeHandler(handler)
    warn_log = logging.getLogger('py.warnings')
    for handler in warn_log.handlers[:]:
        warn_log.removeHandler(handler)
    # QA figures are written in-line inside a worker; see qa.init_qa_pool.
    qa.init_qa_pool(1)


def _run_one(key):
    """
    Execute one work item in a worker process.

    Returns
    -------
    :obj:`tuple`
        ``(key, result, traceback_or_None, log_records)``
    """
    from threadpoolctl import threadpool_limits
    from pypeit import qa

    func, args, kwargs, payload = _FORK_PAYLOAD

    buf = _RecordBuffer()
    log.addHandler(buf)
    try:
        # threadpool_limits handles runtimes that were already initialised in
        # the parent, which the environment variables above cannot.
        with threadpool_limits(limits=1):
            result = func(key, *args, **kwargs, **payload.get(key, {}))
            qa.flush_qa()
        tb = None
    except BaseException as e:          # reported to the parent, not swallowed
        result = None
        tb = ''.join(traceback.format_exception(type(e), e, e.__traceback__))
    finally:
        log.removeHandler(buf)
    return key, result, tb, buf.records


def map_over_detectors(func, detectors, ncpu=1, args=(), kwargs=None, payload=None,
                       label='detector'):
    """
    Map a per-detector function over a list of detectors.

    ``func`` is called as ``func(det, *args, **kwargs, **payload[det])`` and must
    be a module-level function (a requirement for the future ``spawn`` support;
    with ``fork`` a closure happens to work, but do not rely on it).

    When ``ncpu <= 1`` this is a plain ``for`` loop executed in the calling
    process: no pool, no pickling, no log capture.  That path is byte-identical
    to the pre-parallel code and is the default.

    When ``ncpu > 1`` the pool is entered *even for a single detector*, so that
    the parallel code path is exercised by the single-detector CI reduction.

    Parameters
    ----------
    func : callable
        Per-detector function.  First argument is the detector.
    detectors : :obj:`list`
        Detectors/mosaics.  Ints for single detectors, tuples for mosaics.
        Defines the deterministic result-assembly order.
    ncpu : :obj:`int`, optional
        ``par['rdx']['ncpu']``.
    args : :obj:`tuple`, optional
        Positional arguments common to every detector.
    kwargs : :obj:`dict`, optional
        Keyword arguments common to every detector.
    payload : :obj:`dict`, optional
        ``{det: {kwarg: value}}`` of *per-detector* keyword arguments.  Keeping
        the large per-detector objects here (instead of in ``args``) means each
        worker only ever touches its own slice.
    label : :obj:`str`, optional
        Noun used in the progress log lines.

    Returns
    -------
    :obj:`dict`
        ``{det: result}``, iterated in the order of ``detectors``.

    Raises
    ------
    PypeItError
        Raised, in detector order, if any worker raised.
    """
    global _FORK_PAYLOAD
    _kwargs = {} if kwargs is None else kwargs
    _payload = {} if payload is None else payload

    # ---------------- serial: the literal pre-existing loop ----------------
    if ncpu is None or ncpu <= 1:
        return {det: func(det, *args, **_kwargs, **_payload.get(det, {}))
                for det in detectors}

    # ---------------------------- parallel ---------------------------------
    nw = nworkers(ncpu, len(detectors))
    log.info(f'Processing {len(detectors)} {label}(s) using {nw} worker process(es).')
    _FORK_PAYLOAD = (func, tuple(args), _kwargs, _payload)
    results, errors, records = {}, {}, {}
    try:
        # NOTE: a fresh executor per call; the workers are forked here and
        # inherit _FORK_PAYLOAD copy-on-write.
        with ProcessPoolExecutor(max_workers=nw, mp_context=mp.get_context('fork'),
                                 initializer=_worker_init) as ex:
            futures = {}
            for det in detectors:
                futures[ex.submit(_run_one, det)] = det
                log.info(f'Started {label} {det}')
            for fut in as_completed(futures):
                det, result, tb, recs = fut.result()
                results[det], errors[det], records[det] = result, tb, recs
                log.info(f'Finished {label} {det}')
    finally:
        _FORK_PAYLOAD = None

    # Replay the worker logs grouped in detector order, so each detector's block
    # stays contiguous in the run log.  Formatter timestamps come from
    # record.created, i.e. the worker's real clock.
    for det in detectors:
        recs = records.get(det, [])
        if len(recs) == 0:
            continue
        log.info(f'----- begin log for {label} {det} -----')
        for rec in recs:
            log.handle(rec)
        log.info(f'----- end log for {label} {det} -----')

    # Re-raise the first failure, in detector order (deterministic).
    for det in detectors:
        if errors.get(det) is not None:
            raise PypeItError(f'{label} {det} failed in a worker process:\n{errors[det]}')

    return {det: results[det] for det in detectors}


def surviving_detectors(detectors, success):
    """
    Rebuild the list of detectors that should continue to be reduced.

    Parameters
    ----------
    detectors : :obj:`list`
        The detectors that were attempted, in order.
    success : :obj:`dict`
        ``{det: bool}``.

    Returns
    -------
    :obj:`list`
        The subset of ``detectors``, in the original order, for which
        ``success[det]`` is True.
    """
    return [det for det in detectors if success.get(det, False)]
```

**Notes on the design of the helper**

- **Result determinism.** Futures are reaped with `as_completed` (so live
  progress is timely), but *every* consumer iterates `detectors`, never the
  completion order. Assembly is therefore order-independent by construction.
- **Thread pinning.** `os.environ` is set in `_worker_init`, which only helps
  runtimes that have not yet initialised. `threadpoolctl.threadpool_limits(1)`
  inside `_run_one` also constrains runtimes already loaded at fork time.
  `threadpoolctl` is already installed in every PypeIt environment as a
  transitive dependency of `scikit-learn` (a hard dependency, `pyproject.toml`
  line 45) — **add it explicitly to `[project] dependencies`** rather than
  relying on the transitive path.
- **No large object is pickled on the way in.** `_FORK_PAYLOAD` is inherited via
  `fork`; only the detector key crosses the pipe on submit. Return values still
  pickle (see **Q4**).
- **Logging.** Worker handlers are detached in `_worker_init`, so a forked child
  never writes to the parent's log file descriptor. The buffer records only what
  the logger emitted; the parent's handler levels then apply on replay, exactly
  reproducing the serial filtering.

### B.2 Stage 1 — the calibration loop, and the failed-detector bug

**File:** `pypeit/exposure.py`, `reduce_exposure`, lines **421–433**.

**Current code (buggy):**

```python
    # #####################################
    # Calibrations
    for det in detectors:
        log.info(f'Calibrating detector {det}')
        # run/load calibration
        caliBrate =  pypeit_steps.calib_one(spectrograph, fitstbl, par, det, calib_ID, calibrations_path,
              show=show, run_state=run_state, reuse_calibs=reuse_calibs)
        if not caliBrate.success:
            log.warning(
                f'Calibrations for detector {det} were unsuccessful!  The step that failed was '
                f'{caliBrate.failed_step}.  Continuing by skipping this detector.'
            )
            # Remove from list of detectors
            detectors.remove(det)
            continue
```

`detectors.remove(det)` mutates the list **being iterated**. With
`detectors = [(1,5), (2,6), (3,7), (4,8)]` and mosaic `(2,6)` failing, the list
becomes `[(1,5), (3,7), (4,8)]` while the loop index advances to 2 — so `(3,7)`
is **never calibrated**, yet stays in `detectors` and is passed to stages 2–4,
which will then try to load calibrations that were never built.

**New code:**

```python
    # #####################################
    # Calibrations
    #  NOTE: Workers return only (success, failed_step); the Calibrations object
    #  itself is never returned.  The calibrations are written to
    #  Calibrations/*.fits and the downstream stages reload them from disk (see
    #  pypeit_steps.load_calibrations_for_frame), so nothing heavy has to cross
    #  the process boundary.
    calib_status = parallel.map_over_detectors(
        pypeit_steps.calib_status_one, detectors,
        ncpu=par['rdx']['ncpu'],
        args=(spectrograph, fitstbl, par, calib_ID, calibrations_path),
        kwargs=dict(show=show, run_state=run_state, reuse_calibs=reuse_calibs),
        label='detector')

    # Rebuild the surviving-detector list *after* the loop.  The previous
    # implementation called detectors.remove(det) while iterating over
    # detectors, which silently skipped the following detector.
    for det in detectors:
        success, failed_step = calib_status[det]
        if not success:
            log.warning(
                f'Calibrations for detector {det} were unsuccessful!  The step that failed was '
                f'{failed_step}.  Continuing by skipping this detector.'
            )
    detectors = parallel.surviving_detectors(
        detectors, {det: calib_status[det][0] for det in detectors})

    if len(detectors) == 0:
        log.warning('Calibrations were unsuccessful for every detector; nothing to reduce.')
        return spec2dobj.AllSpec2DObj(), specobjs.SpecObjs()
```

The early return is safe: `reduce_calibID` (`pypeit/pypeit.py:369`) already
guards on `len(this_spec2d.detectors) > 0`, and `AllSpec2DObj.detectors`
(`pypeit/spec2dobj.py:484`) returns `[]` for an empty container.

Add the import at the top of `pypeit/exposure.py` (near line 16):

```python
from pypeit import parallel
```

**New adapter** in `pypeit/pypeit_steps.py`, next to `calib_one` (after line
199). It exists because (a) the mapped function must take the detector first,
and (b) it is where the design's `(det, success, failed_step)` contract is
implemented — the helper supplies the `det` key, so the payload is
`(success, failed_step)`.

```python
def calib_status_one(det, spectrograph, fitstbl, par, calib_ID, calibrations_path:str,
                     reuse_calibs:bool=True, qa_path:str=None, show:bool=False,
                     run_state:dict=None, stop_at_step:str=None):
    """
    Build/load the calibrations for one detector and return only the status.

    Thin wrapper around :func:`calib_one` used by
    :func:`~pypeit.parallel.map_over_detectors`.  The
    :class:`~pypeit.calibrations.Calibrations` object is deliberately *not*
    returned: it is large, and every downstream consumer reloads the
    calibrations from ``Calibrations/*.fits`` anyway.

    Args:
        det (:obj:`int`, :obj:`tuple`):
            Detector or mosaic.  First argument, as required by
            :func:`~pypeit.parallel.map_over_detectors`.
        (remaining arguments as for :func:`calib_one`)

    Returns:
        tuple: ``(success, failed_step)``; a :obj:`bool` and a :obj:`str` (or
        None).
    """
    log.info(f'Calibrating detector {det}')
    caliBrate = calib_one(spectrograph, fitstbl, par, det, calib_ID, calibrations_path,
                          reuse_calibs=reuse_calibs, qa_path=qa_path, show=show,
                          run_state=run_state, stop_at_step=stop_at_step)
    return caliBrate.success, caliBrate.failed_step
```

### B.3 Stage 2 — `process_exposure`

**File:** `pypeit/exposure.py`, lines **119–138**.

New adapter in `pypeit/pypeit_steps.py` (after `process_one_det`, line 319):

```python
def process_one_det_bydet(det, spectrograph, fitstbl, par, frames:list, calib_ID:str,
                          calibrations_path:str, bg_frames:list=None,
                          sci_outfile:str=None, bkg_outfile:str=None):
    """
    Detector-first adapter for :func:`process_one_det`, for use with
    :func:`~pypeit.parallel.map_over_detectors`.

    Returns:
        tuple: ``(sciImg, bkg_redux_sciimg)``; see :func:`process_one_det`.
    """
    log.info(f'Reducing detector {det}')
    return process_one_det(spectrograph, fitstbl, par, frames, det, calib_ID,
                           calibrations_path, bg_frames=bg_frames,
                           sci_outfile=sci_outfile, bkg_outfile=bkg_outfile)
```

Loop body replacement:

```python
    # BEFORE (exposure.py:119-138)
    sciImg_dict = {}
    bkg_redux_sciimg_dict = {}
    for det in detectors:
        log.info(f'Reducing detector {det}')
        sciImg, bkg_redux_sciimg = pypeit_steps.process_one_det(...)
        sciImg_dict[det] = sciImg
        bkg_redux_sciimg_dict[det] = bkg_redux_sciimg
    return sciImg_dict, bkg_redux_sciimg_dict

    # AFTER
    out = parallel.map_over_detectors(
        pypeit_steps.process_one_det_bydet, detectors,
        ncpu=par['rdx']['ncpu'],
        args=(spectrograph, fitstbl, par, frames, calib_ID, calibrations_path),
        kwargs=dict(bg_frames=bg_frames),
        label='detector')

    # Assemble in the fixed `detectors` order (order-independent by key).
    sciImg_dict = {det: out[det][0] for det in detectors}
    bkg_redux_sciimg_dict = {det: out[det][1] for det in detectors}
    return sciImg_dict, bkg_redux_sciimg_dict
```

### B.4 Stage 4 — `extract_exposure`

**File:** `pypeit/exposure.py`, lines **305–348**.

Adapter in `pypeit/pypeit_steps.py` (after `extract_det`, line 801):

```python
def extract_det_bydet(det, spectrograph, fitstbl, par, frames, calib_ID:str,
                      calibrations_path:str, sciImg=None, final_sky=None,
                      sobjs_obj=None, calib_slits=None, bkg_redux_final_sky=None,
                      bkg_redux:bool=False, find_negative:bool=False, show:bool=False):
    """
    Detector-first adapter for :func:`extract_det`, for use with
    :func:`~pypeit.parallel.map_over_detectors`.  The per-detector objects
    (``sciImg``, ``final_sky``, ``sobjs_obj``, ``calib_slits``,
    ``bkg_redux_final_sky``) are supplied through the mapper's ``payload``.

    Returns:
        tuple: ``(spec2DObj, sobjs)``; see :func:`extract_det`.
    """
    return extract_det(spectrograph, fitstbl, par, frames, det, calib_ID,
                       calibrations_path, sciImg, final_sky, sobjs_obj, calib_slits,
                       bkg_redux_final_sky=bkg_redux_final_sky, bkg_redux=bkg_redux,
                       find_negative=find_negative, show=show)
```

Loop body replacement:

```python
    # AFTER (replaces exposure.py:312-345)
    # Slice the per-detector inputs *in the parent*, so each worker receives
    # only its own data.
    payload = {}
    detnames = {}
    for i, det in enumerate(detectors):
        detname = sciImg_dict[det].detector.name
        detnames[det] = detname
        if all_specobjs_objfind.nobj > 0:
            sobjs_on_det = all_specobjs_objfind[all_specobjs_objfind.DET == detname]
        else:
            sobjs_on_det = all_specobjs_objfind
        payload[det] = dict(sciImg=sciImg_dict[det],
                            final_sky=final_sky_dict[det],
                            sobjs_obj=sobjs_on_det,
                            calib_slits=calib_slits[i],
                            bkg_redux_final_sky=bkg_redux_final_sky_dict[det])

    out = parallel.map_over_detectors(
        pypeit_steps.extract_det_bydet, detectors,
        ncpu=par['rdx']['ncpu'],
        args=(spectrograph, fitstbl, par, frames, calib_ID, calibrations_path),
        kwargs=dict(bkg_redux=bkg_redux, find_negative=find_negative),
        payload=payload, label='detector')

    # Assemble in the fixed `detectors` order.
    for det in detectors:
        all_spec2d[detnames[det]], tmp_sobjs = out[det]
        if tmp_sobjs.nobj > 0:
            all_specobjs_extract.add_sobj(tmp_sobjs)
        # NOTE: preserved verbatim from the serial loop, including the fact that
        # this overwrites `calibs` on every iteration so the container ends up
        # holding the *last* detector's association.  Changing that is out of
        # scope for this PR.
        all_specobjs_extract.calibs = calibrations.Calibrations.get_association(
                                fitstbl, spectrograph, calibrations_path,
                                fitstbl[frames[0]]['setup'],
                                fitstbl.find_frame_calib_groups(frames[0])[0], det,
                                must_exist=True, proc_only=True)
```

### B.5 Stage 3 stays serial

`findobj_on_exposure` (`pypeit/exposure.py:140–264`) is untouched. Add a comment
at line 202 recording why:

```python
    # NOTE: This loop is deliberately *not* parallelized.  It fans out per
    # detector, hits a cross-detector barrier (adjust_for_slitmask, line 228,
    # which aggregates objects and slits over all detectors to compute the
    # slitmask offset), and then fans out again.  Parallelizing the
    # fan-out/barrier/fan-out structure is deferred to v2 of the detector
    # parallelism; see pypeitdev/speed_up/Reports/speed_up_design.md §8.
```

`PypeIt.calib_all` (`pypeit/pypeit.py:186–190`) is a second flat detector loop
and is an easy, optional bonus conversion using the same helper. Do it in this
PR only if it costs nothing (it uses `calib_status_one` unchanged); otherwise
leave a TODO.

### B.6 Tests (PR B)

**New file `pypeit/tests/test_parallel.py`** — plain module-level functions per
the repo convention, no large data, deterministic:

```python
"""
Tests of the per-detector parallel mapping helper.
"""
import os
import numpy as np
import pytest

from pypeit import log
from pypeit import PypeItError
from pypeit import parallel


# Module-level so it is importable in a forked worker.
def _work(det, scale=1.0):
    """Deterministic float work with a reduction whose order BLAS can change."""
    rng = np.random.default_rng(1000 + (det if isinstance(det, int) else sum(det)))
    a = rng.normal(size=(200, 200))
    u, s, vt = np.linalg.svd(a @ a.T)
    log.info(f'worked on {det}')
    return float(scale) * s.sum()


def _fail(det):
    if det == 2:
        raise RuntimeError('boom')
    return det


def _env(det):
    from threadpoolctl import threadpool_info
    return os.environ.get('OMP_NUM_THREADS'), [i['num_threads'] for i in threadpool_info()]


DETS = [1, 2, 3, 4]


def test_nworkers():
    assert parallel.nworkers(1, 4) == 1
    assert parallel.nworkers(None, 4) == 1
    assert parallel.nworkers(8, 3) == 3                 # capped by nitem
    assert parallel.nworkers(1000, 100) == max(1, (os.cpu_count() or 1) - 1)


def test_serial_and_parallel_are_bit_identical():
    serial = parallel.map_over_detectors(_work, DETS, ncpu=1, kwargs=dict(scale=2.0))
    par = parallel.map_over_detectors(_work, DETS, ncpu=3, kwargs=dict(scale=2.0))
    assert list(par.keys()) == DETS, 'results must be assembled in detector order'
    for det in DETS:
        assert serial[det] == par[det], f'det {det} differs bit-for-bit'


def test_parallel_is_deterministic():
    a = parallel.map_over_detectors(_work, DETS, ncpu=3)
    b = parallel.map_over_detectors(_work, DETS, ncpu=3)
    assert all(a[d] == b[d] for d in DETS)


def test_payload_is_per_detector():
    payload = {d: dict(scale=float(d)) for d in DETS}
    out = parallel.map_over_detectors(_work, DETS, ncpu=2, payload=payload)
    ref = {d: _work(d, scale=float(d)) for d in DETS}
    assert all(out[d] == ref[d] for d in DETS)


def test_worker_threads_are_pinned():
    out = parallel.map_over_detectors(_env, [1, 2], ncpu=2)
    for det, (omp, nthreads) in out.items():
        assert omp == '1'
        assert all(n == 1 for n in nthreads), 'BLAS/OpenMP threads not pinned'


def test_worker_logs_replayed_in_detector_order(caplog):
    with caplog.at_level('INFO', logger='pypeit'):
        parallel.map_over_detectors(_work, DETS, ncpu=3)
    order = [int(m.split()[-1]) for m in caplog.messages if m.startswith('worked on ')]
    assert order == DETS


def test_worker_exception_names_the_detector():
    with pytest.raises(PypeItError, match='detector 2'):
        parallel.map_over_detectors(_fail, DETS, ncpu=3)


def test_surviving_detectors_rebuild():
    """Regression for the detectors.remove()-during-iteration bug."""
    status = {1: True, 2: False, 3: True, 4: True}
    assert parallel.surviving_detectors([1, 2, 3, 4], status) == [1, 3, 4]
    # The old in-loop remove() would have skipped detector 3 entirely.
```

**Modified `pypeit/tests/test_runpypeit.py`** — the end-to-end identical-output
check. `shane_kast_blue` has a single detector, so `nworkers()` collapses to 1;
`map_over_detectors` still enters the pool because `ncpu>1`, which is exactly why
the helper is written that way — the fork, the thread pinning, the buffered
logging and the result marshalling all get CI coverage on a 2-minute reduction.

```python
def test_run_pypeit_ncpu():
    """
    Reduce shane_kast_blue twice, serially and with --ncpu 2, into separate
    output directories and require the products to match to machine precision.
    """
    if os.environ.get('PYPEIT_SKIP_SLOW_TESTS'):
        pytest.skip('slow test')
    # ... same Setup.main(...) preamble as test_run_pypeit, into two dirs ...
    RunPypeIt.main(RunPypeIt.parse_args([str(pyp_file), '-o', '-r', str(serial_dir)]))
    RunPypeIt.main(RunPypeIt.parse_args([str(pyp_file), '-o', '-r', str(par_dir),
                                         '--ncpu', '2']))
    s = specobjs.SpecObjs.from_fitsfile(serial_spec1d)
    p = specobjs.SpecObjs.from_fitsfile(par_spec1d)
    assert s.nobj == p.nobj
    for key in ('OPT_COUNTS', 'OPT_COUNTS_IVAR', 'OPT_WAVE', 'BOX_COUNTS', 'TRACE_SPAT'):
        np.testing.assert_array_equal(getattr(s[0], key), getattr(p[0], key))
    assert s[0].WAVE_RMS == p[0].WAVE_RMS
```

Use `assert_array_equal` (exact), not `allclose` — the design requires machine
precision, and any difference means the thread pinning is not doing its job.

**Dev suite** — the genuine multi-detector check. Add a `reduce` variant that
runs an existing 2-detector setup (e.g. `keck_lris_blue/multi_600_4000_d560`)
with `--ncpu 2` alongside the serial run, plus a `vet_tests` test that compares
the two `spec1d`/`spec2d` files array-by-array. Use a 2-detector setup rather
than DEIMOS: a 3.5 h duplicate run is too expensive for the suite (see **Q6**).

### B.7 Validation and re-profiling (PR B)

```bash
cd $PYPEIT_DEV/pypeitdev/speed_up/scripts
python profile_deimos.py --ncpu 1     # baseline, wall + cProfile
python profile_deimos.py --ncpu 4     # wall-clock only
python profile_kast_blue.py --ncpu 1
python profile_kast_blue.py --ncpu 2  # must not regress
```

**Important measurement caveat:** the profile scripts run `cProfile` in-process
(`profile_deimos.py:88–99`). Child processes are **not** profiled, so the
`--ncpu>1` `.prof` file is blind to ~70% of the work. For PR B the primary
metric is `wall_clock_s` in `Reports/*.runmeta.json`; keep cProfile only for the
`ncpu=1` baseline. Also record peak RSS (`/usr/bin/time -v`, or a `psutil`
sampler in the driver) — memory is the real constraint here.

**Expected speedup (Amdahl, DEIMOS, 4 mosaics, `--ncpu 4`).** Parallel fraction
p ≈ (5 800 + 2 952 + 800) / 12 390 ≈ **0.77**; serial remainder (stage 3, I/O,
save) ≈ 0.23. Ideal speedup = 1/(0.23 + 0.77/4) = **2.4×**. With load imbalance
(per-mosaic calibration wall spans 1 338–1 584 s, ±10%) and the fork/pickle
overhead, expect a realistic **1.9–2.3×**:

> **DEIMOS: ~11 700 s (post-PR-A) → ~5 200–6 200 s.**

`shane_kast_blue` is single-detector and gains nothing; the requirement is only
that `--ncpu 2` does not *regress* it by more than the fork overhead (a second
or two).

**Memory:** expect peak RSS to scale roughly with `nworkers`. DEIMOS already
uses tens of GB serially; document that `--ncpu 4` on a 32 GB machine may swap
and that the practical cap is memory, not cores.

### B.8 Documentation (PR B)

- `doc/running.rst`: extend the "Running on multiple CPUs" section written in PR
  A — what is parallelized (calibration build, image processing, extraction),
  what is not (object finding + global sky, because of the slitmask barrier),
  the RAM-per-detector tradeoff, and that `ncpu>1` currently requires a
  `fork`-capable platform (Linux; macOS/Windows fall back to serial — see
  **Q3** note in §7).
- `doc/api/pypeit.parallel.rst` is auto-generated by `sphinx-automodapi`; run
  `./update_docs` and `git add doc/api`.
- `doc/releases/2.1.0dev.rst`:

  ```rst
  Functionality/Performance Improvements and Additions
  ----------------------------------------------------
  - PypeIt can now reduce detectors/mosaics concurrently.  With ``[rdx] ncpu``
    (or ``run_pypeit --ncpu``) greater than 1, the calibration build, the
    science-image processing, and the extraction are distributed over worker
    processes, one detector each.  Object finding and global sky subtraction
    remain serial because they require a cross-detector slitmask barrier.
    Outputs are identical to the serial run to machine precision.  Peak memory
    usage scales with the number of detectors reduced simultaneously.

  Bug Fixes
  ---------
  - Fixed a bug in :func:`~pypeit.exposure.reduce_exposure` where a detector
    whose calibrations failed was removed from the detector list *while that
    list was being iterated*, causing the following detector to be skipped
    silently and then reduced with missing calibrations.

  Under-the-hood Improvements
  ---------------------------
  - Added :mod:`~pypeit.parallel`, a small helper that maps a per-detector
    function over the detector list either serially or over a
    ``ProcessPoolExecutor``.  Workers pin their BLAS/OpenMP thread count to 1
    and buffer their log records, which the parent replays grouped in detector
    order.

  Testing
  -------
  - Added :mod:`pypeit.tests.test_parallel` and an ``--ncpu`` end-to-end
    regression to ``test_runpypeit``.
  ```

- Add `threadpoolctl` to `[project] dependencies` in `pyproject.toml` and add a
  bullet under **Dependency Changes**.

---

## 3. PR C — Vectorized arc-line Gaussian fitting

**Branch:** `speed_up_arcfit` → **target:** `speed_up_detpar`

**File:** `pypeit/core/arc.py`, `fit_arcspec` (lines **1077–1143**).

### C.1 What is being replaced

The loop at lines 1121–1142 calls `fitting.fit_gauss` (`pypeit/core/fitting.py:874`)
once per detected line, which calls `scipy.optimize.curve_fit` with
`gauss_3deg` (`fitting.py:921`) and no `sigma`. On DEIMOS this is **1.11 M
`curve_fit` calls, 93 M `gauss_3deg` evaluations, `fit_gauss` cum 840 s**.

Windows are tiny: `nfitpix = round(1.25*fwhm)` (`arc.py:1026`), so for the
typical `fwhm=4` the window is `fitp_even = 6`, `fit_interval = 3`, **7 pixels**.
`yarray` is the *continuum-subtracted* arc (`arc.py:1002`), so wing pixels are
routinely ≤ 0.

`fit_arcspec` and `fit_gauss` are only used together here; `fit_gauss` has one
other, cold caller (`pypeit/multislit_flexure.py:77,84`), which is untouched.

### C.2 The math

For `y = A exp(-(x-c)²/(2σ²))` with `y > 0`, take logs and fit a parabola in the
window-centred coordinate `u = x - x[pixt]`:

```
ln y = a0 + a1 u + a2 u²
```

with the closed-form back-substitution

```
a2 = -1/(2σ²)                 →  σ  = sqrt(-1/(2 a2))      (requires a2 < 0)
a1 = (c-x0)/σ² = -2 a2 u0     →  u0 = -a1 / (2 a2)   ,  c = x0 + u0
a0 = ln A + a2 u0²            →  A  = exp(a0 - a1²/(4 a2))
```

An **unweighted** log-space fit is badly biased: it weights the faint wings as
heavily as the peak. `curve_fit` here minimises `Σ (y_i - g_i)²` in *linear*
space. Writing `y = g·e^δ` with `δ = ln y - p(u)`, the linear-space residual is
`≈ g δ`, so the log-space fit reproduces the same objective to first order when
each row is weighted by `w_i = y_i` (i.e. minimise `Σ y_i² δ_i²`). This is the
standard Guo weighting. The normal equations are then 3×3 per line with moments

```
S_k = Σ_i w_i² u_i^k        (k = 0..4)
T_k = Σ_i w_i² u_i^k ln y_i (k = 0..2)

[S0 S1 S2] [a0]   [T0]
[S1 S2 S3] [a1] = [T1]
[S2 S3 S4] [a2]   [T2]
```

All of which are `np.sum(..., axis=1)` over an `(nline, nwin)` array — one
`np.linalg.solve` on a stacked `(nline, 3, 3)` array replaces 1.11 M
`curve_fit` calls.

Centring on `u` rather than using absolute pixel `x` (up to 4096) is essential:
the absolute-`x` Vandermonde has condition number ~10¹⁴ for a 7-pixel window at
`x ≈ 3000`.

**`centerr`.** The existing code returns `fitcov[1,1]`, i.e. the *variance* of
the centre (despite the name), and it is consumed only as `all_ecent` in
`wvutils.arc_lines_from_spec` (`wvutils.py:353`) → `HolyGrail` pattern weights
(`autoid.py:2116–2136`). Reproduce it by propagating the parabola covariance,
`Cov(a) = s² (XᵀWX)⁻¹` with `s² = Σ w² δ² / (n_good - 3)`:

```
∂c/∂a1 = -1/(2 a2)
∂c/∂a2 =  a1/(2 a2²)
var(c) = J1² C11 + 2 J1 J2 C12 + J2² C22
```

### C.3 Implementation sketch

```python
USE_VECTORIZED_ARCFIT = True
"""
Module switch for the vectorized arc-line fit.  Set to False to restore the
per-line ``scipy.optimize.curve_fit`` implementation (used as the fallback for
individual lines regardless).
"""

ARCFIT_GN_ITER = 3
"""
Number of vectorized Gauss-Newton refinement iterations applied to the analytic
log-parabola solution.  Set to 0 to use the pure analytic fit.
"""

ARCFIT_MIN_GOODPIX = 4
"""Minimum number of positive, in-bounds pixels required for the analytic fit."""


def fit_arcspec(xarray, yarray, pixt, fitp):
    """
    Fit a series of pre-identified arc spectrum lines.

    The implementation is a simple 3-parameter Gaussian (amplitude, centroid,
    width), fit for all lines simultaneously with a weighted analytic
    log-parabola solution (optionally polished with a few vectorized
    Gauss-Newton iterations).  Lines for which the analytic solution is
    unusable fall back to the per-line `scipy.optimize.curve_fit`_ fit.

    (Parameters and Returns unchanged from the previous implementation.)
    """
    fitp_even = fitp if fitp % 2 == 0 else fitp + 1
    fit_interval = fitp_even // 2

    sz_p = pixt.size
    sz_a = yarray.size
    ampl    = np.full(sz_p, -999.0, dtype=float)
    cent    = np.full(sz_p, -999.0, dtype=float)
    widt    = np.full(sz_p, -999.0, dtype=float)
    centerr = np.full(sz_p, -999.0, dtype=float)
    if sz_p == 0:
        return ampl, cent, widt, centerr
    if not USE_VECTORIZED_ARCFIT:
        return _fit_arcspec_curvefit(xarray, yarray, pixt, fit_interval, sz_a,
                                     np.arange(sz_p), ampl, cent, widt, centerr)

    # ---- windows -------------------------------------------------------
    # Identical to the old loop: symmetric about the peak, truncated at the
    # spectrum edges (out-of-range positions are masked, not shifted).
    off = np.arange(-fit_interval, fit_interval + 1)
    idx = pixt[:, None] + off[None, :]                       # (nline, nwin)
    inb = (idx >= 0) & (idx < sz_a)
    idxc = np.clip(idx, 0, sz_a - 1)
    y = yarray[idxc]
    x = xarray[idxc]
    x0 = xarray[pixt]
    u = x - x0[:, None]

    # ---- pixel mask ----------------------------------------------------
    # The log fit needs strictly positive flux; the continuum-subtracted arc
    # goes negative in the wings, and NaNs can arrive from masked pixels.
    gpm = inb & np.isfinite(y) & (y > 0.)
    nwin = inb.sum(axis=1)          # == (pmax - pmin) in the old loop
    ngood = gpm.sum(axis=1)
    # The two 'continue' guards of the old loop, plus a minimum for the parabola
    usable = (nwin > 0) & (nwin >= fit_interval) & (ngood >= ARCFIT_MIN_GOODPIX)

    # ---- batched weighted log-parabola ---------------------------------
    w2 = np.where(gpm, y, 0.)**2
    lny = np.log(np.where(gpm, y, 1.))
    S = [np.sum(w2 * u**k, axis=1) for k in range(5)]
    T = [np.sum(w2 * u**k * lny, axis=1) for k in range(3)]
    A = np.empty((sz_p, 3, 3), dtype=float)
    for j in range(3):
        for k in range(3):
            A[:, j, k] = S[j + k]
    b = np.stack(T, axis=1)

    det3 = np.linalg.det(A)
    solvable = usable & np.isfinite(det3) & (np.abs(det3) > 1e-12 * np.maximum(S[0], 1e-300)**3)
    coef = np.full((sz_p, 3), np.nan)
    if np.any(solvable):
        coef[solvable] = np.linalg.solve(A[solvable], b[solvable])
    a0, a1, a2 = coef.T

    with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
        u0 = -0.5 * a1 / a2
        sig = np.sqrt(-0.5 / a2)
        amp = np.exp(a0 - 0.25 * a1**2 / a2)
        ok = (solvable & np.isfinite(a2) & (a2 < 0.)
              & np.isfinite(u0) & np.isfinite(sig) & np.isfinite(amp)
              & (np.abs(u0) <= fit_interval))      # do not extrapolate out of the window

    # ---- optional vectorized Gauss-Newton polish -----------------------
    # Refines against the *exact* curve_fit objective (linear space, all
    # in-bounds pixels, including negative wings), so the result tracks
    # curve_fit to ~1e-6 instead of only to the log-space approximation.
    if ARCFIT_GN_ITER > 0 and np.any(ok):
        amp, u0, sig, ok = _gauss_newton_refine(u, y, inb, amp, u0, sig, ok,
                                                niter=ARCFIT_GN_ITER)

    # ---- centre variance ------------------------------------------------
    var_c = _center_variance(A, coef, u, w2, lny, ngood, ok)

    ampl[ok]    = amp[ok]
    cent[ok]    = x0[ok] + u0[ok]
    widt[ok]    = sig[ok]
    centerr[ok] = var_c[ok]

    # ---- per-line fallback ----------------------------------------------
    bad = usable & np.logical_not(ok)
    if np.any(bad):
        nbad = int(np.sum(bad))
        if nbad > 0.25 * sz_p:
            log.warning(f'Analytic Gaussian fit failed for {nbad}/{sz_p} arc lines; '
                        'falling back to curve_fit for those lines.')
        _fit_arcspec_curvefit(xarray, yarray, pixt, fit_interval, sz_a,
                              np.where(bad)[0], ampl, cent, widt, centerr)

    return ampl, cent, widt, centerr
```

with the legacy loop preserved verbatim as the fallback (and as the reference
implementation for the tests):

```python
def _fit_arcspec_curvefit(xarray, yarray, pixt, fit_interval, sz_a, which,
                          ampl, cent, widt, centerr):
    """
    The original per-line `scipy.optimize.curve_fit`_ implementation, retained as
    a fallback for lines the analytic fit cannot handle.  Fills ``ampl``,
    ``cent``, ``widt`` and ``centerr`` in place for the indices in ``which``.
    """
    for p in which:
        pmin = max(int(pixt[p]) - fit_interval, 0)
        pmax = min(int(pixt[p]) + fit_interval + 1, sz_a)
        if pmin == pmax or (pmax - pmin) < fit_interval:
            continue
        try:
            fitc, fitcov = fitting.fit_gauss(xarray[pmin:pmax], yarray[pmin:pmax])
            ampl[p], cent[p], widt[p] = fitc
            centerr[p] = fitcov[1, 1]
        except RuntimeError:
            pass
    return ampl, cent, widt, centerr
```

The Gauss-Newton polish, also fully batched:

```python
def _gauss_newton_refine(u, y, gpm, amp, u0, sig, ok, niter=3, lam=1e-8):
    """
    A few damped Gauss-Newton iterations on the exact `curve_fit`_ objective,
    ``sum_i (y_i - A exp(-(u_i-u0)^2 / 2 sigma^2))^2``, for all lines at once.
    Unlike the log-parabola seed, this uses *every* in-bounds pixel, including
    the negative wings of the continuum-subtracted arc.
    """
    m = gpm.astype(float)
    a, c, s = amp.copy(), u0.copy(), sig.copy()
    for _ in range(niter):
        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            d = u - c[:, None]
            g = a[:, None] * np.exp(-0.5 * (d / s[:, None])**2)
            r = (y - g) * m
            J = np.stack([g / a[:, None],
                          g * d / s[:, None]**2,
                          g * d**2 / s[:, None]**3], axis=-1) * m[..., None]
            H = np.einsum('nij,nik->njk', J, J)
            H[:, np.arange(3), np.arange(3)] *= (1. + lam)       # Levenberg damping
            rhs = np.einsum('nij,ni->nj', J, r)
            step = np.zeros_like(rhs)
            good = ok & np.isfinite(H).all(axis=(1, 2)) & (np.abs(np.linalg.det(H)) > 0.)
            if np.any(good):
                step[good] = np.linalg.solve(H[good], rhs[good])
            a = np.where(good, a + step[:, 0], a)
            c = np.where(good, c + step[:, 1], c)
            s = np.where(good, s + step[:, 2], s)
        ok = ok & np.isfinite(a) & np.isfinite(c) & np.isfinite(s) & (s > 0.)
    return a, c, s, ok
```

### C.4 Edge cases (each must be explicitly handled and tested)

| Case | Handling |
|---|---|
| **Negative / zero flux** (continuum-subtracted wings) | masked out of the log fit by `gpm`; the Gauss-Newton polish then re-includes them, so the final answer uses the same pixels as `curve_fit` |
| **NaN pixels** | masked by `np.isfinite(y)` |
| **Peaks at the array edges** | window truncated by `inb`; `nwin >= fit_interval` reproduces the old `continue` |
| **Too few positive pixels** | `ngood >= ARCFIT_MIN_GOODPIX` → falls back to `curve_fit`, matching the old behaviour of always attempting a fit |
| **Upward parabola** (`a2 >= 0`, no peak) | excluded by `ok`; falls back |
| **Singular / ill-conditioned normal matrix** | `det3` guard; falls back |
| **Amplitude overflow** in `exp(a0 - a1²/(4a2))` | `np.errstate(over='ignore')` + `np.isfinite(amp)` → falls back |
| **Saturated / flat-topped lines** | fit is biased, exactly as `curve_fit` is; already rejected downstream by the `tampl_true < nonlinear_counts` test at `arc.py:1043` |
| **Centre outside the fit window** | `abs(u0) <= fit_interval`; downstream `abs(tcent-pixt) < fwhm*0.75` (`arc.py:1044`) still applies |
| **`sz_p == 0`** | early return of the empty sentinel arrays |

The `-999.0` sentinels and the `good` mask at `arc.py:1042–1044` are unchanged,
so a line that the analytic fit cannot handle and that `curve_fit` also fails on
is rejected exactly as before.

### C.5 Tests (PR C)

Extend **`pypeit/tests/test_arc.py`** (keep the existing
`test_detect_lines` assertion `len(arx_w) > 3275` — it is a real-data
regression on the new code):

```python
def _synthetic_arc(nspec=2048, nline=200, fwhm=4.0, seed=42):
    rng = np.random.default_rng(seed)
    x = np.arange(nspec, dtype=float)
    cen = np.sort(rng.uniform(20, nspec - 20, nline))
    cen = cen[np.concatenate(([True], np.diff(cen) > 4 * fwhm))]
    amp = 10.**rng.uniform(1.5, 3.5, cen.size)
    sig = fwhm / 2.355 * rng.uniform(0.9, 1.1, cen.size)
    y = np.zeros(nspec)
    for a, c, s in zip(amp, cen, sig):
        y += a * np.exp(-0.5 * ((x - c) / s)**2)
    y += rng.normal(scale=1.0, size=nspec)
    return x, y, amp, cen, sig


def test_fit_arcspec_recovers_truth():
    x, y, amp, cen, sig = _synthetic_arc()
    pixt = np.round(cen).astype(int)
    a, c, w, ce = arc.fit_arcspec(x, y, pixt, 5)
    good = c > 0
    assert good.sum() > 0.95 * cen.size
    assert np.max(np.abs(c[good] - cen[good])) < 0.05      # pixels
    assert np.max(np.abs(w[good] / sig[good] - 1.)) < 0.05


def test_fit_arcspec_matches_curve_fit():
    """The vectorized fit must agree with the legacy per-line curve_fit."""
    x, y, _, cen, _ = _synthetic_arc()
    pixt = np.round(cen).astype(int)
    new = arc.fit_arcspec(x, y, pixt, 5)
    arc.USE_VECTORIZED_ARCFIT = False
    try:
        old = arc.fit_arcspec(x, y, pixt, 5)
    finally:
        arc.USE_VECTORIZED_ARCFIT = True
    good = (new[1] > 0) & (old[1] > 0)
    assert np.max(np.abs(new[1][good] - old[1][good])) < 0.02   # centroid, pixels
    assert np.max(np.abs(new[2][good] / old[2][good] - 1.)) < 0.02


def test_fit_arcspec_edge_cases():
    x = np.arange(64, dtype=float)
    y = np.zeros(64)
    y[32] = 100.
    y[:8] = -5.                                    # all-negative window
    pixt = np.array([0, 2, 32, 63])
    a, c, w, ce = arc.fit_arcspec(x, y, pixt, 5)
    assert np.all(np.isfinite([a, c, w, ce]))      # no NaN leaks to the caller
    assert c[2] == pytest.approx(32., abs=0.5)


def test_fit_arcspec_no_curve_fit_on_clean_data(monkeypatch):
    """On clean data the vectorized path must not call curve_fit at all."""
    x, y, _, cen, _ = _synthetic_arc()
    def _boom(*a, **k):
        raise AssertionError('fit_gauss should not be called')
    monkeypatch.setattr('pypeit.core.fitting.fit_gauss', _boom)
    arc.fit_arcspec(x, y, np.round(cen).astype(int), 5)
```

Also add a micro-benchmark assertion only if the CI machine is stable enough;
otherwise measure by hand.

### C.6 Numerical vetting via the dev suite

The design accepts noise-level changes vetted by dev-suite RMS checks. Run, at
minimum:

```bash
cd $PYPEIT_DEV
./pypeit_test reduce -i shane_kast_blue shane_kast_red keck_deimos keck_lris_blue
pytest vet_tests/test_wavelengths.py vet_tests/test_wavetilts.py \
       vet_tests/test_extraction.py vet_tests/test_skysub.py \
       --redux_out $PYPEIT_DEV/REDUX_OUT
```

`vet_tests/test_wavelengths.py` asserts per-setup wavelength RMS
(`test_shane_kast_red`, `test_deimos`, `test_keck_lris_blue`,
`test_keck_lris_red`, `test_gmos`, `test_keck_hires`, …) and
`vet_tests/test_wavetilts.py::test_run` asserts the tilt RMS — these are exactly
the two products this change touches. Record before/after RMS values in the PR
description; the acceptance criterion is "no setup's RMS degrades beyond its
existing tolerance, and no setup loses arc lines".

### C.7 Validation and re-profiling (PR C)

Re-profile both setups at `--ncpu 1` (so cProfile sees everything) and confirm:

- `pypeit/core/fitting.py:fit_gauss` disappears from the top of the
  `.pstats.txt` table (DEIMOS baseline: cum 840 s, 1.11 M calls).
- `gauss_3deg` self-time (332 s, 93 M calls) and `scipy` `leastsq` (586 s) drop
  to noise.
- **Expected: ~700–800 s off DEIMOS (~6% of the original wall), ~1–3 s off
  Kast.** The remaining `curve_fit` calls should be the fallback only —
  typically <1% of lines.

### C.8 Documentation (PR C)

- `doc/releases/2.1.0dev.rst`, under **Under-the-hood Improvements**:

  ```rst
  - Replaced the per-line ``scipy.optimize.curve_fit`` Gaussian fit in
    :func:`~pypeit.core.arc.fit_arcspec` with a vectorized weighted analytic
    log-parabola fit (optionally polished with a few batched Gauss-Newton
    iterations), fit for all detected arc lines at once.  This removes ~1.1
    million ``curve_fit`` calls from a Keck/DEIMOS reduction and speeds up both
    wavelength calibration and tilt tracing.  Line centroids change at the
    noise level.
  ```

- No user-facing parameter is added, so `doc/pypeit_par.rst` is unchanged.
  `doc/api/pypeit.core.arc.rst` regenerates automatically.

---

## 4. Cross-cutting validation checklist

Run at the tip of each PR branch, before opening the PR:

1. `pytest pypeit/tests` (the CI-safe unit tests). Must be green, including the
   new `test_parallel.py`, `test_qa.py` and `test_arc.py` additions.
2. `$PYPEIT_DEV/pypeit_test reduce -i shane_kast_blue keck_deimos` — a
   smoke reduction on the two profiled setups.
3. `pytest $PYPEIT_DEV/vet_tests --redux_out $PYPEIT_DEV/REDUX_OUT` for the
   relevant modules (all of them for PR C).
4. Re-profile with `pypeitdev/speed_up/scripts/profile_{kast_blue,deimos}.py`
   and `analyze_profile.py`, and write the numbers into a short
   `Reports/speed_up_results.md` (baseline / after-A / after-B / after-C).
5. **Determinism:** run the same reduction twice at `--ncpu 4` and diff the
   `spec1d`/`spec2d` arrays; they must be bit-identical.
6. **Identical output:** diff `--ncpu 1` vs `--ncpu 4` products; they must be
   bit-identical (this is the locked requirement, and the reason for the BLAS
   pinning).
7. Peak-RSS measurement at `--ncpu 1` and `--ncpu 4` on DEIMOS; record it in the
   docs.

Full dev suite (`./pypeit_test all -t 2`) once, at the tip of `speed_up`, before
merging to `develop`.

---

## 5. Documentation checklist (all PRs)

| Item | A | B | C |
|---|:-:|:-:|:-:|
| `doc/pypeit_par.rst` regenerated (`./update_docs`) | ✔ | | |
| `doc/help/run_pypeit.rst` regenerated | ✔ | | |
| `doc/running.rst` — "Running on multiple CPUs" | ✔ (create) | ✔ (extend) | |
| `doc/qa.rst` — `Agg`, concurrent PNG writes | ✔ | | |
| `doc/api/` regenerated + `git add doc/api` | | ✔ | ✔ |
| `doc/releases/2.1.0dev.rst` bullets | ✔ | ✔ | ✔ |
| `CHANGES.rst` (deprecated per its own header; optional) | ○ | ○ | ○ |
| `pyproject.toml` — explicit `threadpoolctl` | | ✔ | |
| Docstrings on all new functions (Numpy style) | ✔ | ✔ | ✔ |

Docs build check: `cd doc; make htmlonly` (no `$PYPEIT_DEV` needed) or the full
`./update_docs`; `doc/sphinx_warnings.out` must be empty.

---

## 6. Suggested commit sequence

**PR A** (5 commits)
1. `Add [rdx] ncpu parameter` (pypeitpar + test)
2. `Add run_pypeit --ncpu; plumb through PypeIt.__init__`
3. `Force the Agg backend for reductions`
4. `Defer QA figure writes to an optional thread pool` (qa.py machinery + tests)
5. `Convert the per-slit QA writers to qa.save_figure` (+ docs, release notes)

**PR B** (6 commits)
1. `Add pypeit.parallel: fork-based per-detector mapper` (+ tests)
2. `Add detector-first adapters in pypeit_steps`
3. `Parallelize the calibration loop; fix detectors.remove() during iteration`
4. `Parallelize process_exposure`
5. `Parallelize extract_exposure`
6. `Add --ncpu end-to-end regression; docs and release notes`

**PR C** (3 commits)
1. `Vectorize fit_arcspec with an analytic weighted log-parabola fit`
2. `Add Gauss-Newton polish and curve_fit fallback` (+ tests)
3. `Docs, release notes, dev-suite RMS vetting results`

---

## 7. Risks and revert plan

| Risk | Detection | Revert / mitigation |
|---|---|---|
| **BLAS pinning backfires** (a worker's linear algebra is now single-threaded and slower than expected, or `threadpool_limits` deadlocks) | wall-clock at `--ncpu 4` barely improves, or workers hang | Drop the `threadpool_limits` context from `_run_one` (keep the env vars). The identical-output test then becomes the tripwire for FP-order changes. The design explicitly says *"be prepared to revert if it backfires."* |
| **`fork` + OpenMP is undefined behaviour** — a forked child inheriting an active OpenMP thread pool can deadlock | worker hangs, run never completes | Try `OMP_NUM_THREADS=1` set in the *parent's* environment before launching `run_pypeit`; longer term, move to `spawn` (needs module-level `func` and picklable payload — the helper is already written for this, only `_FORK_PAYLOAD` would have to become an explicit argument). |
| **Non-Linux platforms** — `mp.get_context('fork')` is unavailable/unsafe on Windows and unsafe on macOS | `ValueError` from `get_context` | Guard in `map_over_detectors`: if `'fork' not in mp.get_all_start_methods()`, log a warning and fall through to the serial loop. Document that `ncpu>1` is Linux-only in v1. |
| **Memory blow-up** — N detectors in flight ≈ N× peak RSS; DEIMOS already uses tens of GB | OOM kill or swap thrash | The `min(ndet, cpu_count-1)` cap plus documentation. Recommend `--ncpu 2` as the default advice for DEIMOS-class data on ≤32 GB machines. |
| **Stage-4 return pickling** — a `Spec2DObj` carries ~8 full-frame float arrays (~250 MB/detector) through a pipe | measurable gap between the sum of worker times and the wall-clock | See **Q4**: have the worker write the per-detector `spec2d` to `Intermediate/` and return the path. |
| **matplotlib thread-safety** — concurrent `Figure.savefig` on distinct figures is supported but not formally guaranteed | corrupted or truncated PNGs; the `test_save_figure_threaded_matches_serial` pixel comparison catches it | `qa.init_qa_pool(1)` — a one-line revert that restores in-line writing while keeping the deferred-save API. |
| **`plt.close('all')` destroys pending figures** (`autoid.py:67,200,242`) | missing or blank QA PNGs | Covered by the conversion rule in §A.5(2); grep for `close('all')` before merging PR A. |
| **Arc-fit accuracy regression** | `vet_tests/test_wavelengths.py`, `test_wavetilts.py` RMS assertions; `test_detect_lines` line count | `arc.USE_VECTORIZED_ARCFIT = False` restores the exact previous behaviour; the legacy loop is retained verbatim as `_fit_arcspec_curvefit`. |
| **cProfile blindness in the parallel path** | the `--ncpu>1` `.prof` looks implausibly cheap | Documented in §B.7: use `wall_clock_s` from `runmeta.json` as the metric; keep cProfile for `ncpu=1` only. |
| **Nested parallelism** — QA threads inside detector worker processes | oversubscription | `_worker_init` calls `qa.init_qa_pool(1)`; see **Q5**. |

---

## 8. Open questions

These are the residual decisions the design does not settle. None of them blocks
starting the work; each has a recommended default that is what the sketches above
implement.

**Q1 — Does forcing `Agg` have to respect `run_pypeit -s/--show`?**
The design says force `Agg` "unconditionally for headless reductions (no new
param)". But `run_pypeit -s/--show` deliberately raises blocking matplotlib
windows in several debug paths, which `Agg` silently turns into no-ops. The
sketch in §A.3 guards with `if not args.show:` — no new parameter, but not
literally unconditional.
*Recommendation:* keep the `--show` guard.

**Q2 — If the QA thread pool turns out to be GIL-bound, what do we do?**
The dominant QA cost is matplotlib's FreeType text-metrics path
(`_get_text_metrics_with_cache_impl`, cum 894 s on DEIMOS), which may not release
the GIL; PIL's PNG encode (276 s) does. So the thread pool may recover only the
encode half. Options: (a) accept the partial win and keep only `Agg` + deferred
saves; (b) escalate QA to the PR-B process pool; (c) drop QA parallelism.
*Recommendation:* (a) — measure in PR A, keep the machinery (it costs nothing at
`ncpu=1`), and revisit only if QA is still a top-5 hotspot after PR B.

**Q3 — Pure analytic log-parabola, or log-parabola + a few vectorized
Gauss-Newton iterations?**
The design specifies the analytic weighted log-parabola. Adding 2–3 batched
Gauss-Newton steps (§C.3, `ARCFIT_GN_ITER`) is still fully vectorized and calls
no scipy, but it refines against the *exact* `curve_fit` objective using all
pixels including the negative wings, so the results track `curve_fit` to ~1e-6
instead of only to the log-space approximation. Cost: ~3 extra passes over an
`(nline, 7)` array — negligible. Benefit: the "noise-level change" claim becomes
almost trivially true and the dev-suite RMS vet becomes a formality.
*Recommendation:* ship with `ARCFIT_GN_ITER = 3`, and report both variants' RMS
in the PR.

**Q4 — Stage-4 worker returns: pickle the `Spec2DObj`, or write it to disk?**
A `Spec2DObj` carries roughly eight full-frame float arrays (~250 MB for a
DEIMOS mosaic), all of which must be pickled through a pipe back to the parent.
That is probably 1–3 s per detector against a ~700 s stage — acceptable — but it
is the one place the design's "no large objects cross the process boundary"
principle is violated (the design does say stages 2/4 "return as today").
*Recommendation:* measure it in PR B; if the marshalling is >5% of the stage,
switch the worker to `spec2DObj.to_file(Intermediate/...)` and return the path.

**Q5 — Nested `ncpu`: how should QA threads and detector processes compose?**
With `--ncpu 4`, PR A wants 4 QA threads and PR B wants 4 worker processes; done
naively that is 4 processes × 4 QA threads. The sketch resolves this by having
`parallel._worker_init` call `qa.init_qa_pool(1)`, so QA is written in-line
inside a detector worker and threaded only in the (single-process) serial
stages.
*Recommendation:* confirm this — QA serial inside detector workers.

**Q6 — Where should the real multi-detector identical-output regression live?**
`shane_kast_blue` (the only full reduction in `pypeit/tests`) is
single-detector, so CI can exercise the pool machinery but not multi-detector
result assembly. The dev suite is the natural home, but a duplicate DEIMOS
reduction costs ~3.5 h (or ~1.5 h parallel) of suite time.
*Recommendation:* add the check on a cheap 2-detector setup (e.g.
`keck_lris_blue/multi_600_4000_d560`) as a `reduce` variant plus a `vet_tests`
comparison, and run DEIMOS at `--ncpu 4` only manually before the merge to
`develop`.
