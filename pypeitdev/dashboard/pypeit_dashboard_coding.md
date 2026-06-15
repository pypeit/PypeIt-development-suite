# PypeIt Dashboard — Coding Document

**Version:** 1.3.0
**Date:** 2026-06-15
**Author:** JXP and Claude

**Changelog**
- 1.3.0 (2026-06-15): **Stage 6 implemented** — the Science view
  (`view/science_view.py`), model science accessors + science-derive-on-launch,
  `inspect` spec2d/per-object-spec1d/`reduce_by_step` builders, the
  `reduce_by_step` state-write fix + `ScienceFrameState.raw_files`, and the
  Status-view science summary. The six-stage build order is complete.
- 1.2.6 (2026-06-15): **Stage 6 discussion round 2** — design resolved: per-object
  1D via reconstructed `--obj`; science (Re)Build = re-run-with-confirmation +
  prerequisite gating + the `reduce_by_step` state-write fix.
- 1.2.5 (2026-06-15): **Stage 6 discussion** — Science build-order row updated
  for the third **Science tab**, per-object `pypeit_show_1dspec`, and the user-
  requested **(Re)Build** via `pypeit_reduce_by_step` (state-write fix likely
  needed, à la `run_to_calibstep`).
- 1.2.4 (2026-06-15): **Stage 6 prep** — Science build-order row updated to the
  implemented state layer + the `pypeit/state/` package
  (`pypeit.state.science_status.derive_science_from_disk`, `get_science_status()`).
- 1.2.3 (2026-06-15): **Post-Stage-5 consistency pass** — clarified that the
  per-step `*_state.json` writes happen only during a **real run** (a status-only
  derive is read-only); synced `state.rst`/`dashboard.rst`/design X4/X5.
- 1.2.2 (2026-06-15): **Stage 5 Round 2** — a status-only derive no longer
  writes `*_state.json` (`run_the_steps` writes gated by `not status_only`), so
  a Dashboard launch leaves no stale `running` state.
- 1.2.1 (2026-06-15): **Stage 5 Round 1** — `pypeit_run_to_calibstep` now writes
  `*_state.json` (run_state → `calib_one` + final status pass); reload/status
  derive marks missing calibs `undone` not `fail`.
- 1.2.0 (2026-06-15): **Stage 5 implemented** — live monitoring (the
  build-order row marked implemented): `RunLock.stateChanged`,
  `MainWindow._refresh_from_state`, two-channel `ActivityBar`.
- 1.1.4 (2026-06-15): **Stage 5 planning** — §Monitoring updated to the decided
  mechanism (poll `*_state.json` mtime on the `RunLock` timer + a `stateChanged`
  signal; the `.log` tail view deferred).
- 1.1.3 (2026-06-14): **Stage 4 Round 3** — the (Re)Build button turns orange +
  "⏳ Run in progress" + disabled while any run is active (lock visual cue, all
  steps), via a shared `_style_rebuild_button`.
- 1.1.2 (2026-06-14): **Stage 4 Round 2** — processed-image inspect → `ginga`
  directly; (Re)Build button orange-while-running + keeps the rebuilt step
  selected; output filename shown in the detail panel.
- 1.1.1 (2026-06-14): **Stage 4 Round 1** — crash-safe clobber (move-aside +
  restore-on-failure), the `pypeit_run_to_calibstep --calib_group` `IndexError`
  fix, and a bright teal "Inspect output" button beside the blue (Re)Build.
- 1.1.0 (2026-06-14): **Stage 4 implemented.** Marked the *Execution, locking &
  (Re)Build* build-order row implemented: the blue **(Re)Build** control
  (`pypeit_run_to_calibstep`), the `RunLock` single-run lock (`.log`-mtime +
  launched-run lifecycle), the `QMessageBox` clobber confirmation, the
  delete-then-run clobber (in code), and refresh-on-completion. New modules
  `pypeit/dashboard/runlock.py` + `inspect.run_command`.
- 1.0.1 (2026-06-14): **Stage 4 consistency pass.** Corrected the processed-FITS
  viewer to `pypeit_view_fits --inter`/`--proc` (was "`ginga` (processed)"), and
  the Stage 3 build-order row to the requirements actually delivered
  (C1–C9, C11–C14, C16, X4, X5) with C10/C15 listed under Stage 4.
- 1.0.0 (2026-06-14): First **versioned release** of the Dashboard docs
  (semantic versioning starts here), coinciding with Stages 0–3 implemented.
- 0.6 (2026-06-13): **Pulled the Dashboard status/activity area (X4/X5) forward
  from Stage 4 into Stage 3** in the *Developing → Build order* table (per Stage 2
  sign-off feedback): the user wants visible "GUI is busy / waiting on a job"
  feedback as soon as the dashboard launches external tools, and so the
  Calibrations view (Stage 3) can show what Refresh/derive is doing.
- 0.5 (2026-06-13): Clarified the Stage 0 **launch mechanism** in the
  *Developing → Build order* table — the dashboard launches **like every other
  PypeIt script** (e.g. `run_pypeit`) via the `pypeit_dashboard` console script
  (a `ScriptBase` entry point), with the `.pypeit` file as a **required
  positional** argument (`python -m pypeit.dashboard` optional).
- 0.1 (2026-06-13): Initial draft. Added the **GUI Package** section recording
  the decision to build on PyQt6 via `qtpy`, the import convention, the reuse
  policy for `setup_gui/` infrastructure, the external-viewer (subprocess)
  approach, and the resulting coding conventions and dependency posture.
- 0.2 (2026-06-13): Added the **Developing** section — the MVC architecture
  (headless data layer vs. Qt views), the six-stage build order, what to build
  first, the active-run monitoring approach (`*_state.json` polling + `.log`
  mtime/tail; broadcast hook removed from `state.py`), the v1 launch scope
  (`pypeit_run_to_calibstep` only), and the test strategy.
- 0.3 (2026-06-13): Added the **Checking and modifying the layout** section —
  the offscreen render→PNG capability, the three-tier checking workflow
  (render-and-view self-check → user confirmation at milestones → structural
  `pytest-qt` tests written throughout), the light/dark × default/resized render
  matrix, the layout-check harness (real widgets, `shane_kast_blue` 600/4310
  data, in `pypeitdev/dashboard/py/`), and the sketches-vs-renders doc policy.
- 0.4 (2026-06-13): Added the **Debugging** section — automated-first headless
  debugging by bug class (headless `pytest`, `qtbot` interaction sim,
  render-to-PNG, a loud `sys.excepthook`/slot-wrap, deterministic `*_state.json`
  fixtures, subprocess capture, activity-log trace), how correctness is confirmed
  (automated-weighted, user only at milestones), CI-safe test/fixture placement
  (main repo `pypeit/tests/` by default; dev suite only when data demands), and
  the per-incident debug-handoff files in `pypeitdev/dashboard/debugging/`.

---

## Purpose

This is the **coding document** for the PypeIt Dashboard. Where
`pypeit_dashboard_design.md` records *what* the dashboard must do (requirements,
views, behavior), this document records *how* we will build it: the GUI package
and conventions, code organization, reuse policy, and the development/testing
practices the implementation must follow. It is a living document, grown section
by section alongside the design document.

It is grounded in the same context: `pypeit_dashboard_design.md`,
`pypeit_workflow.md`, the prototype `PypeIt/pypeit/dashboard/dashboard.py` (a
reference only — not the basis for the implementation), and the existing
`PypeIt/pypeit/setup_gui/` GUI.

---

## GUI Package

### Decision

The dashboard is built on **PyQt6, accessed through `qtpy`** — the same GUI
stack as the existing `pypeit/setup_gui/`. This was selected against the user's
requirements (Python; few/no new dependencies; PyQt6 already in PypeIt; favor
ease of development and maintenance over performance) and confirmed with the
user before drafting this document.

`qtpy` is a thin abstraction layer that exposes a single, Qt5-style API on top
of whichever Qt binding is installed (PyQt6, PySide6, etc.). PypeIt resolves
that binding to **PyQt6** via its dependencies.

### Why this choice

| Requirement | How PyQt6-via-`qtpy` satisfies it |
|-------------|-----------------------------------|
| Must use Python | PyQt6 is a Python binding for Qt; all dashboard code is Python. |
| Add few/no new dependencies | **Zero** new dependencies — PypeIt's `setup.cfg` already pins `qtpy>=2.2.0`, `pyqt6`, and `pytest-qt`. |
| PyQt6 already included | The design doc and `setup_gui/` both commit to PyQt6; reusing it avoids fragmenting the GUI stack. |
| Ease of development/maintenance over performance | Qt is the most complete, best-documented Python desktop toolkit, and there is an existing, working PyQt6 codebase (`setup_gui/`) to learn from and reuse. |

Alternatives were considered and rejected: **Tkinter** (stdlib, no new dep) is
too weak for the rich, color-coded, tab/grid/scrollable views the design calls
for and shares nothing with `setup_gui`; **web stacks** (Dash/Streamlit/Flask)
add heavy dependencies, a browser/runtime split, and a second language;
**Dear PyGui / Kivy / wxPython** each add a brand-new toolkit with no existing
PypeIt code to reuse. PyQt6 is the only option that scores well on every
requirement.

### Import convention (`qtpy`, not direct `PyQt6`)

All Qt imports go through **`qtpy`**, never directly through `PyQt6`. This
matches `setup_gui/` and keeps the Qt backend swappable. Concretely:

```python
# Correct — backend-agnostic, matches setup_gui
from qtpy.QtWidgets import QWidget, QVBoxLayout, QTableWidget
from qtpy.QtCore import Qt, QObject, Signal, Slot
from qtpy.QtGui import QColor, QIcon
```

```python
# Avoid — locks us to one binding and mismatches setup_gui
from PyQt6.QtCore import pyqtSignal   # use qtpy's `Signal` instead
```

Consequence for signals/slots: use the **Qt5-style names** `Signal` / `Slot`
(as `qtpy` exposes them and as `setup_gui` uses), **not** PyQt6's
`pyqtSignal` / `pyqtSlot`. Note the reference prototype `dashboard.py` mixes the
two styles (`from PyQt6.QtCore import pyqtSignal` alongside `from qtpy...`); the
dashboard implementation will **not** carry that inconsistency forward.

### Reuse of `setup_gui/` infrastructure

Per the user's direction, the dashboard reuses **as much of the `setup_gui`
infrastructure as possible**, refactoring shared pieces into a common location
when that is cleaner than duplicating them. Components identified as directly
reusable:

- **`setup_gui/text_viewer.py`** — `LogWindow` and `TextViewerWindow`, directly
  applicable to the dashboard's log view and any text/QA-text display.
- **`setup_gui/dialog_helpers.py`** — file-open/save dialogs, the
  `DialogResponses` enum, and `PersistentStringListModel`, reusable for the
  dashboard's file selection and confirmation dialogs (including the clobber
  confirmation in the design doc's *Execution, Locking & Status* section).
- **Application bootstrap** — the window-icon, default-font, and theming setup
  used to start the setup GUI; the dashboard should share the same look and
  startup conventions.
- **MVC organization** — `setup_gui` is structured as
  `model.py` / `view.py` / `controller.py`. The dashboard follows the same
  separation so the two GUIs are structurally familiar to maintainers.

Where reuse would force awkward coupling, the dashboard provides its own view
code instead; shared utilities are factored out only when both GUIs genuinely
need them. When a shared component is refactored out of `setup_gui`, the setup
GUI must continue to work unchanged.

### Inspection via external tools (no embedded plotting)

The dashboard **launches the existing PypeIt CLI inspection tools as
subprocesses** rather than embedding live plots:

- Calibration/spectra inspection: `pypeit_chk_*` / `pypeit_show_*` scripts and
  `ginga`.
- Raw input frames: `pypeit_view_fits --proc` (oriented/processed view).
  Processed calibration outputs (`Bias_*`/`Arc_*`/…): opened **directly in
  `ginga`** (Stage 4 Round-2 #1; `pypeit_view_fits` does not display them).

Consequently, **v1 needs no embedded `matplotlib`-Qt canvas or in-window
plotting backend**, and the GUI-package choice does not have to account for one.
QA **PNG** files are still shown *inline*, but those are static images
(displayed with Qt's own image widgets), not live plots. Subprocess launching
uses the standard-library `subprocess` module and must respect the single-run
lock and observability rules in the design document's *Execution, Locking &
Status* section.

### Coding conventions (GUI code)

These apply to all dashboard GUI code (consistent with the prep document's
coding guidelines and the main PypeIt conventions):

- **Language:** Python only; Qt accessed exclusively through `qtpy`.
- **Imports at the top** of each module; no inline imports.
- **Line length ≤ 80 characters.**
- **Docstrings** on every class and method describing purpose and
  inputs/outputs (Numpy style preferred, matching PypeIt).
- Methods authored for the dashboard include **"Generated by JXP and Claude"**
  in their docstring.
- **Inline comments** explain non-obvious GUI wiring (signal/slot connections,
  threading, subprocess launches).
- **Reuse before reinvention:** prefer `setup_gui` components and existing
  PypeIt machinery over new implementations.

### Code location

- **Production dashboard code:** `PypeIt/pypeit/dashboard/`.
- **Design-phase scratch/prototype code:**
  `PypeIt-development-suite/pypeitdev/dashboard/py/`.

### Dependency posture

No additions to PypeIt's dependency list are required for the GUI package: the
dashboard relies entirely on the already-present `qtpy`, `pyqt6`, and (for
tests) `pytest-qt`. Any future need for an embedded plotting backend or other
new dependency must be raised and justified separately, as it would change this
posture.

---

## Developing

This section records *how* the dashboard is built: the code architecture, the
staged build order, what is implemented first, how an active reduction is
monitored, the launch scope for v1, and the test strategy. It reflects the
decisions reached in the *Developing* discussion (see the prep document).

### Architecture (MVC, mirroring `setup_gui/`)

The dashboard follows the same **Model–View–Controller** organization as
`pypeit/setup_gui/`:

- **Model — a headless data/logic layer.** All PypeIt-facing logic lives here
  and is kept **independent of Qt**: loading/deriving the reduction state,
  reducing it to a status table, mapping `(required, status)` to the
  Dashboard-wide palette + glyph, and computing the path-aware step order. Being
  Qt-free, it is unit-testable with plain `pytest` (no display required), and it
  keeps the GUI thin.
- **View — Qt widgets.** The `MainWindow`, the tab bar, the header banner/logo,
  the Status view, and the Calibrations view. Views render what the model
  provides and emit user-intent signals; they hold no PypeIt logic.
- **Controller — wiring + subprocess launching.** Connects view signals to model
  updates, owns the refresh logic, and launches external tools and
  `pypeit_run_to_calibstep` as subprocesses.

**Reuse:** per the GUI Package section, reuse `setup_gui/` infrastructure where
low-friction — the application/font/icon bootstrap, `dialog_helpers.py` (file
dialogs, `DialogResponses`, confirmations), and `text_viewer.py`
(`LogWindow`/`TextViewerWindow` for the log view).

### Build order (six stages)

Development proceeds in stages; earlier stages are prerequisites for later ones.
Each design reference points back to `pypeit_dashboard_design.md`.

| Stage | What is built | Design refs |
|-------|---------------|-------------|
| **0. Walking skeleton** | Launch **like every other PypeIt script** (e.g. `run_pypeit`): the `pypeit_dashboard` console script — a `ScriptBase` subclass (`RunDashboard`) registered as a `pyproject.toml` entry point (optional `python -m pypeit.dashboard` convenience); argparse for the **required positional `.pypeit` argument**; `MainWindow` with the tab bar (Status \| Calibrations), the shared header banner + logo, reusing `setup_gui`'s app/font/icon bootstrap. Views are empty placeholders. | R1, R2, R6 |
| **1. State data layer (headless)** | Load `*_state.json` via `RunPypeItState.load()`; otherwise derive the state via `PypeIt(calib_only=True).calib_all(status_only=True, reload_only=True)`; expose the `get_status()` DataFrame; handle empty/edge states. **No Qt.** | R4, R5, R11 |
| **2. Initialization / Status view** | The landing, *state-first* view: the color + glyph state table, required-vs-optional distinction, global summary strip, configuration-overview navigator grid, scope drop-downs, manual refresh, stale-state indicator, and a **placeholder** Science section (see below). | R3, R7–R13, R15–R18 |
| **3. Calibrations view** | The path-aware step-button row (from `Calibrations.default_steps()`), the detail panel, and **launching external viewers as subprocesses** (`pypeit_chk_*` / `ginga` / `pypeit_view_fits`); per-slit/order drill-down. Builds the reusable subprocess-launch infrastructure. **Also delivers the Dashboard's own status/activity area** (X4/X5) — pulled forward from Stage 4 — so the user can see when the GUI is busy/waiting on a job (e.g. launching a viewer, deriving/reloading state on Refresh). | C1–C9, C11–C14, C16, X4, X5 (C10/C15 deferred to Stage 4) |
| **4. Execution, locking & (Re)Build** *(implemented 2026-06-14)* | The blue **(Re)Build** control (`pypeit_run_to_calibstep`) in the detail panel; a `RunLock` single-run lock via **`.log` mtime watching** (+ the launched-run lifecycle); a `QMessageBox` clobber confirmation naming the file(s) overwritten; clobber = delete-the-step's-output-then-run (in code); state/buttons refresh on completion. (The Dashboard status/activity area, X4/X5, shipped in Stage 3.) | X1–X3, C10, C15 |
| **5. Monitoring (live updates)** *(implemented 2026-06-15)* | Live status refresh while a reduction is active: `RunLock` gained a `stateChanged` signal (poll the `*_state.json` mtime on the one ~2 s timer while active); `MainWindow._refresh_from_state` auto-refreshes both views from the state file (R12), preserving scope/selection; the status bar split into **Build** + **Inspection** channels. (`.log` tail view deferred.) | R14 |
| **6. Science frames** | Populate the Science section. **State layer done (2026-06-15):** `RunPypeItState.science` (`ScienceFrameState`: process/findobj/skysub/extract per `(frame, det)`, products, per-slit/per-object detail) is recorded live via `exposure.reduce_exposure` hooks and derived from disk by `pypeit.state.science_status.derive_science_from_disk` (Science `spec2d`/`spec1d` first, `Intermediate/` fallback); `get_science_status()` returns the per-frame table and `pypeit_status` prints it. The state layer now lives in the **`pypeit/state/` package** (`run_state.py` + `science_status.py`). **Remaining (Stage 6 view):** a third **Science tab** (flat `(frame, det)` rows, `objtype` column, status-as-cells) + per-frame drill-down (per-slit/per-object); a **per-object** `pypeit_show_1dspec` + `pypeit_show_2dspec` launch; and — per the user — a **(Re)Build** control mirroring Calibrations via `pypeit_reduce_by_step` (which, like `run_to_calibstep` pre-Stage-4, needs a state-write fix to be live). | R15, R18 |

### What is implemented first

Stages **0–1** together — a launchable skeleton plus the headless state layer —
followed by Stage **2**, the Status view. This honors the design's *state-first*
principle and produces a runnable artifact early, which the layout-review process
(next discussion) can exercise against real widgets.

### Monitoring an active reduction

For v1 the dashboard monitors a running reduction through the **state file**, not
through an in-process log hook:

- During a **real run** `pypeit/calibrations.py` writes `*_state.json` on
  **every step transition** — it sets a step's `status='running'` before running
  it and `'success'` / `'fail'` afterward, calling `state.write()` each time. So
  `*_state.json` is a **live, per-step status feed** that requires no new
  instrumentation. (A **status-only** pass — the Dashboard/`pypeit_status`
  *deriving* status on launch — is a read and does **not** write, Stage 5 R2.)
- **Plan (Stage 5, decided):** detect that a run is *active* by watching the
  **`.log` mtime** (PypeIt writes to it continuously while running; this is also
  the single-run-lock signal, X1, already in `RunLock`); while active, **poll the
  `*_state.json` mtime** on the **same `RunLock` ~2 s timer** (no
  `QFileSystemWatcher`), and on a change re-read the state file and re-render the
  Status table + calibration buttons, preserving the user's scope/selection.
  `RunLock` gains a `stateChanged` signal for this. A separate **`.log` tail**
  view (reusing `setup_gui/text_viewer.py`'s `LogWindow`) is **deferred** (it may
  not be wanted; see `dashboard_dev_stage_5.md` S5-Q6).
- The previously-prototyped **`log.step` broadcast hook** (a STEP-level log
  record) is **not** used. The commented-out hook has been **removed from
  `pypeit/state.py`** so the source carries no dead path; finer sub-step
  broadcasting can be revisited later if the state-file feed proves too coarse.

### Launch scope (v1)

The dashboard both *observes* and *launches* PypeIt work, but **v1 launching is
limited to `pypeit_run_to_calibstep`** (calibration (re)generation). Launching a
**full `run_pypeit`** reduction from the dashboard is **deferred**; full runs are
only *observed* (via the `.log` mtime / state-file mechanisms above). All
launching is governed by the single-run lock and clobber protection (Stage 4 /
design §*Execution, Locking & Status*).

### Test strategy

Tests are split to match the architecture:

- **Headless tests (plain `pytest`).** Cover the Qt-free data/logic layer — state
  loading/derivation, the status DataFrame, palette/glyph mapping, and path-aware
  step ordering. These live in **`PypeIt/pypeit/tests/`** alongside the rest of
  the unit tests and require no display.
- **Widget-level tests (`pytest-qt`).** Cover the views/controller, **mirroring
  the approach used for the `setup_gui` tests**.

Tests must be deterministic (fixed RNG seeds), consistent with PypeIt's testing
conventions.

---

## Checking and modifying the layout

Building a GUI is iterative: views must be checked against the design document
and adjusted. This section records *how* we verify and refine the dashboard's
layout. The guiding fact is that the dashboard's **real Qt widgets can be
rendered headlessly to a PNG that is then inspected** — so layout checking is
done against the *actual rendered widget*, not a mock-up.

### The enabling capability: offscreen render → image

With Qt's offscreen platform plugin
(`QT_QPA_PLATFORM=offscreen`, set by the dev suite's `source_headless_test.sh`),
a widget can be drawn and captured without a display:

```python
# Build the real widget, then capture it to a PNG (no display needed).
widget = StatusView(example_state)      # a production pypeit.dashboard widget
widget.resize(1650, 900)                # the design's default window size
widget.grab().save("status_view.png")   # QPixmap -> PNG on disk
```

The resulting PNG shows exactly what a user would see — real fonts, colors, and
status glyphs (`✓`, `⏳`, `✗`, `○`, `–`). This is distinct from the
**matplotlib sketches** in `pypeitdev/dashboard/images/` (generated by
`pypeitdev/dashboard/py/dashboard_status_view.py` and
`dashboard_calibrations.py`), which are non-functional mock-ups of the *intended*
look. The sketches remain the design document's reference figures; the renders
described here exercise the *built* widget.

### Three-tier checking workflow

1. **Render-and-view self-check (primary development loop).** A reusable harness
   in `pypeitdev/dashboard/py/` imports the **real** `pypeit.dashboard` widgets,
   builds each view with **representative example data** — the
   `shane_kast_blue` 600/4310 dataset (the same one behind the sketches and the
   workflow worked example) — grabs each to a PNG offscreen, and exits. The
   render is then inspected and compared against the design-doc requirements
   (the R / C / X items) and the sketches, and the widget code is iterated until
   it matches.

2. **User confirmation (final authority on look and intent).** The screenshots
   are presented to the user for sign-off. Confirming the layout therefore uses
   **both** mechanisms: render-and-inspect the real output, *then* present it to
   the user. Screenshots are shown at **view-completion milestones** (when a view
   reaches a reviewable state), not on every intermediate iteration.

3. **Structural `pytest-qt` tests (regression guard).** Tests that assert the
   view's *structure* — e.g. the tab bar is `Status | Calibrations`, the
   calibration step-button row contains the spectrograph's `default_steps()` in
   dependency order, the state table has the `Step | Required | Status | Output`
   columns. These guard against regressions but do **not** judge visual layout
   (tiers 1–2 do). Such tests are written **throughout** development, alongside
   each view, mirroring the `setup_gui` test approach (and per the *Developing*
   test strategy).

### Render matrix

To verify the design's accessibility and responsiveness requirements, each view
is rendered across a small matrix:

- **Themes:** both **light** and **dark** (the design specifies a
  dark-theme-aware palette that must preserve contrast and the color+glyph
  encoding, R10).
- **Window sizes:** the design's **default** size (≈1650×900) **and** a
  **smaller/resized** window, to confirm the layout reflows acceptably (sticky
  headers, scrollable lists, drop-down scoping — R8, R16, R18).

### The layout-check harness

- **Location:** `pypeitdev/dashboard/py/` (beside the existing sketch scripts;
  this is design-phase tooling, not shipped production code).
- **Imports the real widgets** from `pypeit.dashboard` (it does not re-implement
  or mock them), so what it captures is the production layout.
- **Example data:** `shane_kast_blue` 600/4310, for direct comparability with
  the sketches.
- **Output:** PNGs (per view × theme × size) written under
  `pypeitdev/dashboard/py/` output (or a scratch path); these are working
  artifacts for review, not committed reference images.

### Relationship to documentation

The matplotlib **sketches stay in the design document** as the intended-look
reference. Real-widget **renders** are reserved for the *future user
documentation* (where screenshots of the actual application are appropriate),
rather than replacing the design sketches.

---

## Debugging

GUI debugging is normally slow because problems are hard to isolate and the
human is in the loop for every click. Because the bulk of the development is done
by Claude, the strategy here is **automated-first, headless debugging**: most
bugs are reproduced and fixed without a display or any human interaction, and the
user is involved only when genuinely necessary. The thin-view / Qt-free-model
architecture (see *Developing*) is what makes this possible — the heavy logic
lives where it can be exercised by plain `pytest`.

### How bugs are debugged, by class

1. **Logic bugs (the majority) — headless `pytest` on the model layer.** State
   loading/derivation, the status DataFrame, palette/glyph mapping, and
   path-aware step ordering are Qt-free, so they reproduce in plain `pytest` with
   full tracebacks and fast iteration — no display.
2. **Interaction bugs — `pytest-qt` + `qtbot`, offscreen.** User interactions are
   *simulated* programmatically (`qtbot.mouseClick`, `qtbot.keyClicks`,
   `qtbot.waitSignal`) and the resulting state asserted. Questions like "does
   clicking step X open panel Y and launch Z?" become automated tests — no human
   clicking. This mirrors the existing `setup_gui` test harness, which runs
   `pytest-qt` with the **offscreen** platform.
3. **Visual / layout bugs — render-to-PNG self-inspection.** The real widget is
   rendered offscreen to a PNG (the layout-check harness) and inspected directly
   (see *Checking and modifying the layout*).
4. **Crashes — make failures loud (do not let the event loop swallow them).** An
   exception raised inside a Qt slot is, by default, swallowed by the C++ event
   loop: the app silently misbehaves and the traceback is lost. The dashboard
   therefore installs a **global `sys.excepthook`** plus an optional
   **slot-wrapping decorator** so exceptions in slots are **logged with full
   tracebacks** (and surfaced in an error dialog when running in dev mode).
   (`setup_gui` installs no such hook; this is an addition.)
5. **Deterministic reproduction via state fixtures.** A small set of saved
   `*_state.json` fixtures — a **healthy** `shane_kast_blue` 600/4310 state plus
   crafted edge cases (**not-started**, **partial**, **one failed step**,
   **malformed**) — lets state-driven and rendering bugs reproduce **instantly,
   without running `run_pypeit` and without RAW_DATA**.
6. **Subprocess launches — capture and log.** Launches of external tools
   (`pypeit_chk_*`, `ginga`, `pypeit_view_fits`, `pypeit_run_to_calibstep`)
   capture `stdout` / `stderr` / return code, so a failed launch is diagnosable
   from the logs instead of being a silent no-op.
7. **Activity log as a debug trace.** The dashboard's own logging is routed
   through PypeIt's `log`, producing a timestamped record (which signal fired,
   which handler ran, what state changed) that can be read after a headless run —
   more useful than interactive breakpoints inside a running event loop.

### Confirming the GUI works as expected

**Both** mechanisms, weighted to the automated one: correctness is verified
primarily by **running headless** (`pytest-qt`-driven interactions +
render-to-PNG inspection + reading logs/tracebacks); the user is brought in only
at **milestones**, for **interactive "feel"** that is hard to automate, or for
**platform/display-specific** rendering on their machine.

### Minimizing the user's debugging time

- **Self-diagnosing failures** — the loud `sys.excepthook` and structured
  logging mean a failure arrives with a traceback, not a shrug.
- **Deterministic fixtures** — bugs reproduce on any machine without live data,
  so the user rarely needs to reproduce anything.
- **Automated `qtbot` repros** replace "can you try clicking around?".
- **Structured handoff when the user *is* needed** (below).

### Test & fixture placement (CI-safe by default)

- **Prefer the main repo.** Headless `pytest` tests, `pytest-qt` widget tests,
  and the committed `*_state.json` fixtures live in **`PypeIt/pypeit/tests/`**
  (fixtures under `pypeit/tests/files/`). These are CI-safe and self-contained
  (no RAW_DATA), so they run anywhere — this is the default.
- **Fall back to the dev suite only when data demands it.** A test that needs a
  *significant* amount of real data runs in the dev suite (e.g.
  `unit_tests/`, like the `setup_gui` widget tests), not in the main repo.

### Debug handoffs

When a bug genuinely requires the user, the handoff is a **single self-contained
file per incident**, written to
`PypeIt-development-suite/pypeitdev/dashboard/debugging/`. Each handoff file
contains:

- a **one-line reproduction command** (a single `pytest` or harness invocation);
- the **screenshot(s)** of the relevant render;
- a **precise, usually yes/no question** for the user;
- the relevant **activity log** excerpt.

The user is not expected to run the full application or debug open-endedly —
only to answer the specific question (and, if asked, run the one provided
command).

---

## Logs

*(Coding-document work is logged in the parent
`claude_prompts/dashboard/dashboard_coding_prep.md` Logs section.)*
