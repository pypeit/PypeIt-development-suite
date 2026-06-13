# PypeIt Dashboard — Coding Document

**Version:** 0.3
**Date:** 2026-06-13
**Author:** JXP and Claude

**Changelog**
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
- Raw/processed FITS display: `pypeit_view_fits` (raw) and `ginga` (processed).

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
| **0. Walking skeleton** | `pypeit_dashboard` console entry point + `python -m pypeit.dashboard`; argparse for the **required `.pypeit` argument**; `MainWindow` with the tab bar (Status \| Calibrations), the shared header banner + logo, reusing `setup_gui`'s app/font/icon bootstrap. Views are empty placeholders. | R1, R2, R6 |
| **1. State data layer (headless)** | Load `*_state.json` via `RunPypeItState.load()`; otherwise derive the state via `PypeIt(calib_only=True).calib_all(status_only=True, reload_only=True)`; expose the `get_status()` DataFrame; handle empty/edge states. **No Qt.** | R4, R5, R11 |
| **2. Initialization / Status view** | The landing, *state-first* view: the color + glyph state table, required-vs-optional distinction, global summary strip, configuration-overview navigator grid, scope drop-downs, manual refresh, stale-state indicator, and a **placeholder** Science section (see below). | R3, R7–R13, R15–R18 |
| **3. Calibrations view** | The path-aware step-button row (from `Calibrations.default_steps()`), the detail panel, and **launching external viewers as subprocesses** (`pypeit_chk_*` / `ginga` / `pypeit_view_fits`); per-slit/order drill-down. Builds the reusable subprocess-launch infrastructure. | C1–C16 |
| **4. Execution, locking & status** | Single-run lock via **`.log` mtime watching**; clobber confirmation; the Dashboard's own status/activity area; (re)generation via `pypeit_run_to_calibstep`. | X1–X5, C10 |
| **5. Monitoring (live updates)** | Live status refresh while a reduction is active (see *Monitoring*, below). | R14 |
| **6. Science frames** | Populate the Science section once `RunPypeItState` is extended to track per-science-frame status. | R15, R18 |

### What is implemented first

Stages **0–1** together — a launchable skeleton plus the headless state layer —
followed by Stage **2**, the Status view. This honors the design's *state-first*
principle and produces a runnable artifact early, which the layout-review process
(next discussion) can exercise against real widgets.

### Monitoring an active reduction

For v1 the dashboard monitors a running reduction through the **state file**, not
through an in-process log hook:

- `pypeit/calibrations.py` writes `*_state.json` on **every step transition** —
  it sets a step's `status='running'` before running it and `'success'` /
  `'failed'` afterward, calling `state.write()` each time. So `*_state.json` is a
  **live, per-step status feed** that requires no new instrumentation.
- **Plan:** detect that a run is *active* by watching the **`.log` mtime** (PypeIt
  writes to it continuously while running; this is also the single-run-lock
  signal, X1); while active, **poll `*_state.json`** (e.g. a `QTimer` and/or
  `QFileSystemWatcher` with debounce), reload it, and re-render the status table
  and calibration buttons. Separately, **tail the `.log`** into a log view that
  reuses `setup_gui/text_viewer.py`'s `LogWindow`.
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

## Logs

*(Coding-document work is logged in the parent
`claude_prompts/dashboard/dashboard_coding_prep.md` Logs section.)*
