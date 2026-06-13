# PypeIt Dashboard — Design Document

**Version:** 0.7
**Date:** 2026-06-12
**Author:** JXP and Claude

**Changelog**
- 0.1 (2026-06-12): Preamble; first draft of the Initialization section.
- 0.2 (2026-06-12): Revised Initialization to reflect the user's answers in
  Clarifications — science-frame status is forthcoming (design accommodates it),
  metrics live in other windows (not the state table), the `.pypeit` file is a
  required launch argument, blocking the UI on launch is acceptable, and the
  refresh strategy depends on whether PypeIt is running.
- 0.3 (2026-06-12): Added the Calibrations section — the per-(calibration group,
  detector/mosaic) detail view, its color-coded step buttons (grounded in
  `default_steps()` ordering and `state.py` status), the detail panel, and the
  `pypeit_run_to_calibstep` launch path.
- 0.4 (2026-06-12): Revised Calibrations per the user's answers — **unified the
  Dashboard-wide status palette** (running = orange, required-undone = white;
  Initialization key updated to match), `bpm` omitted from the button row, the
  input/output viewers split (raw → `pypeit_view_fits`, processed → `ginga`),
  always launch with `--calib_group`/`--det`, and noted the metric set will grow
  as `state.py` does.
- 0.5 (2026-06-12): Embedded the Status View layout sketch and added the
  anti-clutter scoping design — calib-group/detector **scope drop-downs** (R16),
  a **configuration-overview navigator** grid (R17), and an equally-prominent,
  scalable **Science frames** section (R18).
- 0.6 (2026-06-12): Embedded the Calibrations View sketch; added the **PypeIt
  logo** to the shared header banner (R6); recorded the **magenta selected-step**
  highlight; and refined Calibrations requirements — QA entries are **clickable**
  to a full view (C9), input files and per-slit/order rows are **scrollable
  lists** (C7, C11), and the (re)generate control uses a distinct **action
  color** (C10).
- 0.7 (2026-06-12): Moved the logo to the **top-right corner** (it was clobbering
  the banner) and pointed image references at the new
  `pypeitdev/dashboard/images/` folder (downscaled `pypeit_logo.png` + the two
  sketch PNGs).

---

## Preamble

### What this document is

This is the design document for the **PypeIt Dashboard**: a desktop graphical
application that gives a PypeIt user a single, coherent place to launch, monitor,
and inspect a spectroscopic data reduction. It records the *requirements* and the
*design* of the dashboard — what it must do, how it should be organized, and how
its pieces fit together. It deliberately stays at the level of architecture,
views, and behavior; it does **not** prescribe specific code or APIs (those will
be worked out during development). It is a living document and will grow phase by
phase.

### What the dashboard is for

PypeIt is a powerful but verbose, command-line-driven pipeline. A typical
reduction involves editing a `.pypeit` file, launching a long-running
`run_pypeit` process that emits a torrent of log messages, and then inspecting a
scattered collection of calibration files, 2D/1D spectra, and fixed-format QA
figures using a family of separate `pypeit_chk_*` / `pypeit_show_*` scripts. The
state of a reduction — what has run, what succeeded, what to look at next — is
not surfaced in any one place.

The dashboard's job is to make that workflow **legible and navigable**: to show
at a glance where a reduction stands, what it has produced, where it may have
gone wrong, and to put the existing inspection tools one click away.

### Scope

Per the project goals, the dashboard concentrates on the **core run** of PypeIt —
**Phase 2 (Reduction)** and **Phase 3 (Inspection / QA)** of the workflow
described in `pypeit_workflow.md`. Concretely, it will at least:

- Display the **state** of a PypeIt reduction (a color-coded, formatted status
  view of frames / calibration groups / steps).
- Let the user **run individual steps** of the reduction.
- Let the user **examine standard PypeIt outputs** (calibrations, `spec2d`,
  `spec1d`, QA figures) easily.
- Let the user **launch existing PypeIt scripts** (e.g. `pypeit_chk_edges`,
  `pypeit_show_2dspec`) to inspect files.
- **Monitor the progress** of a running reduction.

**Out of scope (for now):** the Setup phase (frame typing / `.pypeit` file
generation — already served by an existing setup GUI) and the further-processing
phase (fluxing, 1D/2D/3D coaddition, telluric, collation).

### Relationship to other documents

- **`pypeit_workflow.md`** is the companion context document. It captures our
  working understanding of the PypeIt reduction workflow (phases, scripts,
  products, QA touch-points) and includes a concrete worked example
  (Shane Kast blue 600/4310). **All design decisions in this document should be
  grounded in, and cross-referenced to, that workflow.**
- The PypeIt source (`PypeIt/`) and its documentation are the ultimate
  authority on behavior; the dashboard should **reuse existing PypeIt code**
  (e.g. status/state logic, `DataContainer` IO, the `pypeit_chk_*` machinery)
  rather than reimplement it.

### Technology and conventions

- **Language:** Python.
- **GUI toolkit:** PyQt6.
- **Code location:** dashboard code lives in `pypeit/dashboard/`; design-phase
  scratch/prototype code lives in
  `PypeIt-development-suite/pypeitdev/dashboard/py/`.
- **Image assets:** design images — the layout-sketch PNGs and the (downscaled)
  PypeIt logo (`pypeit_logo.png`) — live in
  `PypeIt-development-suite/pypeitdev/dashboard/images/`.
- This document avoids specific code recommendations but assumes the Python +
  PyQt6 stack throughout.

### Design principles (proposed)

These are starting principles to guide the detailed design that follows; they
can be revised as the design matures.

1. **State-first.** The landing view answers "where is this reduction?" before
   anything else.
2. **Reuse, don't reinvent.** Wrap existing PypeIt scripts/logic; treat PypeIt
   as the source of truth for state and products.
3. **Workflow-faithful.** The dashboard's structure mirrors the real PypeIt
   workflow (Reduction → Inspection/QA) so users' mental model transfers.
4. **Path-aware.** MultiSlit / Echelle / IFU differ in products and QA; views
   should adapt to the spectrograph's reduction path.
5. **Non-destructive and observable.** Long-running actions report progress and
   never silently overwrite; the user always knows what the dashboard is doing.

### Document structure (planned)

This document will be built up section by section, following the project's
design phases:

- **Preamble** (this section).
- **Initialization** — how the dashboard launches and what it shows on startup.
- **Calibrations** — the per-(group, detector) calibration detail/inspect view.
- *(subsequent phases: science-frame reduction, Inspection/QA, Monitoring, ...
  — to be added.)*

---

## Initialization

This section sets the requirements and design for the **Initialization phase** —
how the dashboard is launched and what it presents to the user on startup. The
goal of startup is to answer, immediately and unambiguously, the question
*"where does this reduction stand?"* (design principle #1, *state-first*).

### How launch works

- The dashboard is started by running a Python script from within a reduction
  folder — i.e. the per-configuration directory created by `pypeit_setup`
  (e.g. `shane_kast_blue_A/`) that contains the `.pypeit` file. A console entry
  point (working name `pypeit_dashboard`) and/or `python -m pypeit.dashboard` is
  assumed.
- The user provides the name of the `.pypeit` file to use as a **required
  launch argument**. The dashboard does not guess, so the presence of more than
  one `.pypeit` file in the folder is not a problem.
- On startup the dashboard derives the reduction **state**, then renders it as
  the default view: a color-coded, formatted **State Table**. Deriving the state
  may briefly block the UI on launch, which is acceptable (see below).

### State: where it comes from

Grounded in `PypeIt/pypeit/state.py` and `pypeit/scripts/pypeit_status.py`
(see also the worked example in `pypeit_workflow.md` §6):

- The state is modeled by the pydantic class `RunPypeItState`. It records, per
  **calibration group** (`calib_id`) and **detector / mosaic** (`det`), the
  status of each calibration **step**. The steps currently tracked are: `bias`, `dark`,
  `arc`, `tiltimg`, `slits`, `wv_calib`, `tilts`, `scattlight`, `flats`,
  `align`.
- Each step entry carries: a `required` flag, `input_files`, an `output_file`,
  optional `qa_files`, a `status` ∈ {`undone`, `running`, `success`,
  `complete`, `fail`}, and **step-specific metrics**:
  - `bias`: `mean`, `std`
  - `slits`: `nslits`, and per-slit `center`
  - `wv_calib`: per-slit `rms`
  - `tilts`: per-slit `rms`
  - `flats`: `types`

  **Design decision:** these metrics are valuable but are **not** shown in the
  state table. They are surfaced in the Dashboard's dedicated inspection windows
  (Inspection phase). The state table reports only the *status* of each step —
  keeping it a clean, scannable health overview rather than a metrics dump.
- Two sources of state, in priority order:
  1. **State file present** — load `<pypeit_root>_state.json` directly via
     `RunPypeItState.load()`. This is fast and is what a finished or in-progress
     `run_pypeit` leaves behind (confirmed in §6 of the workflow doc).
  2. **No state file** — compute the state the way `pypeit_status.py` does:
     instantiate `PypeIt(..., calib_only=True)` and call
     `calib_all(status_only=True, reload_only=True)`, then read
     `pypeIt.run_state`. No processing is performed. Because the user names the
     `.pypeit` file explicitly, there is no file ambiguity to resolve. This
     computation can take a moment for large multi-slit reductions; it is
     acceptable for it to briefly **block the UI on launch** (no background
     thread is required for initialization).
- **Refresh strategy** (see R12): once running, prefer whichever source is
  authoritative *at that moment* — if a `run_pypeit` reduction is **active**,
  re-read `*_state.json` (the running process owns and updates it); if **no**
  reduction is running, re-derive the state via the `pypeit_status` path above.
- A convenient tabular reduction of the state already exists
  (`RunPypeItState.get_status()` → a pandas DataFrame with columns
  `calibration_group`, `detector`, `steps`, `required`, `status`,
  `output_file`); the dashboard should build on this rather than re-deriving it.

**Important — science-frame status is coming.** Today the state model tracks
**calibration steps only**; it does not yet capture per-science-frame object
finding / extraction status. This is, however, planned: `RunPypeItState` *will*
be extended to record science-frame status. The Initialization design must
therefore **anticipate** science frames from the start — the state view should
extend naturally to a *science* section presented alongside the *calibration*
section (see R8 and R15), rather than being hard-wired to calibrations only.

### Requirements

Requirements are tracked in three buckets so the document doubles as a checklist
as development proceeds.

#### Pending

*Basic requirements (from the project brief):*

- **R1.** The dashboard launches by running a Python script.
- **R2.** It is launched from a folder that contains a `.pypeit` file, and the
  user specifies **which** `.pypeit` file to use as a required launch argument
  (so multiple `.pypeit` files in the folder are unambiguous).
- **R3.** The default view shows the reduction **state** as a formatted table
  with **color-coded** information.
- **R4.** If a state file (`*_state.json`) is present in the folder, load it to
  get the state.
- **R5.** Otherwise, compute the state as `pypeit_status.py` does
  (`PypeIt(calib_only=True).calib_all(status_only=True, reload_only=True)`). It
  is acceptable for this to briefly block the UI on launch.
- **R6. Context header.** Display a header banner with the `.pypeit` file name,
  the spectrograph (`PYP_SPEC`), the setup/configuration ID, the reduction path
  (MultiSlit / Echelle / IFU), and the reduction directory. The **PypeIt logo**
  sits in the **top-right corner** (above the banner, clear of the text). This
  banner+logo is shared by the Status and Calibrations views.
- **R7. At-a-glance summary.** Show an overall progress/health summary derived
  from the state (e.g. "Calibrations: 6/6 required steps success" plus counts of
  fail / running / undone), so health is readable without scanning every row.
  The summary should generalize to also cover science frames once their status
  is tracked (see R15).
- **R8. Grouped, ordered rows.** Show steps in a fixed, logical processing
  order. Rather than stacking every calibration group × detector in one long
  scroll (which does not scale — see R16/R17), the detailed table is **scoped**
  to one `(calibration group, detector/mosaic)` at a time. The grouping scheme
  must extend to a parallel *science-frame* section (see R15/R18) without
  redesign.
- **R9. Required vs. optional clarity.** Visually distinguish *required* steps
  from *optional* ones (e.g. `dark`, `scattlight` were `undone` & not required
  in the Kast run) so an `undone` optional step never reads as a problem.
- **R10. Accessible encoding.** Never rely on color alone: pair every status
  color with a glyph and/or text label (e.g. ✓ success, ✗ fail, ⏳ running,
  ○ undone, – n/a) for colorblind readability.
- **R11. Graceful empty/edge states.** Handle clearly: the named `.pypeit` file
  not found; a valid `.pypeit` file but no state and no outputs ("not started");
  a malformed/partial state file.
- **R12. Manual refresh.** Provide a way to re-read / re-derive the state on
  demand. The refresh **source** depends on whether a reduction is active: if
  `run_pypeit` is **running**, re-read `*_state.json` (the running process owns
  it); if **not running**, re-derive via the `pypeit_status` path. (Automatic,
  continuous refresh while running is R14 / the Monitoring phase.)
- **R13. Stale-state indication.** If `*_state.json` is older than the
  `.pypeit` file or the calibration outputs, indicate the state may be stale.
- **R14. Live status updates (planned).** While a reduction is running, the
  state view should reflect progress in (near) real time. This is planned but
  its **implementation/channel is open** — the commented-out `log.step`
  "broadcast" hook in `state.py` is one candidate (a log-stream channel),
  polling `*_state.json` is another. The detailed design is deferred to the
  **Monitoring** phase; Initialization only needs to not preclude it.
- **R15. Science-frame readiness.** The state view must be designed so that, once
  `RunPypeItState` is extended to track per-science-frame status (object finding
  / extraction), it can be shown as a *science* section alongside the
  *calibration* section — reusing the same status encoding (R10), grouping (R8),
  and summary (R7) — without a redesign.
- **R16. Scope drop-downs.** Provide calibration-group and detector/mosaic
  **drop-downs** that scope the detailed Calibrations and Science tables to a
  single `(group, detector)` — the primary defense against clutter on
  instruments with many calibration groups and/or detectors. (These mirror the
  Calibrations view's selectors, C2.) The always-visible summary strip (R7)
  preserves whole-run health while one scope is shown.
- **R17. Configuration-overview navigator.** Above the scoped tables, show a
  compact **(group × detector) grid**, each cell colored by the *worst* step
  status in that cell (fail > running > to-do > success), per the unified
  palette. It gives whole-run visibility at a glance and acts as a
  **click-to-scope navigator** (clicking a cell sets the drop-downs). For a
  single-group/single-detector run (e.g. Kast) it is a single cell; it scales to
  a heat-map for MOS/mosaic runs.
- **R18. Equally-prominent, scalable Science section.** The Science-frames
  section is given **equal visual weight** to the Calibrations section (same
  heading treatment). For runs with many science frames it is **filterable and
  scrollable** with a visible frame count, so long lists do not overwhelm the
  view. Columns include the per-frame products (`spec2d`, `spec1d`) and will gain
  object-find/extraction status as the state model grows (R15).

#### Implemented

*(None yet — this is a design-phase document; items move here as they are built
and verified.)*

#### Deferred

- **D1.** Editing the `.pypeit` file from within the dashboard — Setup-phase
  concern, out of scope (an existing setup GUI covers it).
- **D2.** Multi-configuration / multi-setup overview (several `.pypeit` runs side
  by side) — possible future top-level view.

### Default State View — visual design

Readability and clarity are the priorities (per the brief). Proposed treatment:

The layout sketch below (generated by
`pypeitdev/dashboard/py/dashboard_status_view.py`, example data:
`shane_kast_blue` 600/4310) illustrates the target:

![PypeIt Dashboard Status View — layout sketch](images/dashboard_status_view.png)

- **Layout.** Top to bottom: tab bar (Status | Calibrations) → header banner
  (R6) → global summary strip (R7) → **scope toolbar** with the calib-group and
  detector drop-downs + filter/refresh (R16) → **configuration-overview
  navigator** grid (R17) → the **Calibrations** section (steps for the selected
  scope) → an equally-prominent **Science frames** section (R18). The two
  sections share the same heading treatment, and the Science section is present
  from the start (populating as the state model grows, R15).
- **Columns.** `Step` | `Required` | `Status` | `Output`. Step names bold; file
  names in a monospace font. Per-step metrics and per-slit detail are
  intentionally **not** columns here — they live in the inspection windows.
- **Status color key** (color + glyph, R10). This is the **Dashboard-wide status
  palette**, reused verbatim by the Calibrations view:

  | Status                  | Color (proposed)            | Glyph |
  |-------------------------|-----------------------------|-------|
  | success / complete      | green (#2E7D32)             | ✓ |
  | running                 | orange (#EF6C00), subtle pulse | ⏳ |
  | fail                    | red (#C62828)               | ✗ |
  | required & not done     | white (#FFFFFF, outlined)   | ○ |
  | optional / not required | grey (#9E9E9E)              | – |
  | not used / n/a          | dimmed light grey (#BDBDBD) | – |

  White needs a subtle border to read on a light background. A dark-theme-aware
  palette should keep equivalent contrast ratios.
- **Scanning aids.** Sticky group/column headers, gentle zebra striping, and a
  filter control (e.g. *required only*, *failures only*) to cut visual load on
  large MOS/echelle reductions.
- **Empty/edge states (R11).** Render as a clear centered message with the next
  action (e.g. *"No reduction outputs found — this reduction has not been run.
  [Run reduction]"*), rather than an empty grid.

### Resolved decisions and remaining open items

The questions raised during the first draft were answered by the user in the
`#### Clarifications` section of the parent `dashboard_design.md`. Resolutions
now folded into this section:

- **Science-frame status** — coming; the design accommodates it (R8, R15).
- **State computation cost** — blocking the UI on launch is acceptable (R5).
- **Which `.pypeit` file** — supplied explicitly by the user at launch (R2).
- **Refresh source** — depends on whether a reduction is running (R12).
- **Metrics** — shown in other (inspection) windows, not the state table.

Remaining open item, carried to a later phase:

- **Live-update channel (R14).** Whether live monitoring rides the `log.step`
  broadcast hook, `*_state.json` polling, or another mechanism is open and will
  be settled in the **Monitoring** phase.

---

## Calibrations

The **Initialization** state table (above) is the at-a-glance overview. The
**Calibrations** view is the drill-down: a separate tab that shows the *detailed*
status of the calibrations for **one** calibration group and **one**
detector / mosaic at a time, and that lets the user inspect each calibration and
(re)generate it. It is where a user answers "is this calibration good, and if
not, what do I do about it?" and allows them to remake it.

This section first lays out the view and the backend it drives, then gives
per-calibration-type specifics, then the formal requirements.

The layout sketch below (generated by
`pypeitdev/dashboard/py/dashboard_calibrations.py`, example data:
`shane_kast_blue` 600/4310 with `wv_calib` selected) illustrates the target:

![PypeIt Dashboard Calibrations View — layout sketch](images/dashboard_calibrations.png)

### Selecting the context

- The view operates on a single `(calibration group, detector/mosaic)` pair,
  chosen by the user via two **drop-down lists** populated from the state
  (`RunPypeItState` enumerates the `(calib_id, det)` pairs; see Initialization).
- Changing either drop-down re-renders the step buttons and clears/refreshes the
  detail panel for the new context.

### The calibration step buttons

The heart of the view is a row of **clickable step buttons**, one per
calibration step, laid out **left → right in dependency order** so the visual
order matches the order in which PypeIt builds them (and the order in which
`pypeit_run_to_calibstep` would generate them). Connectors/arrows between
buttons reinforce "this step precedes that one".

**Path-awareness (which steps, in which order).** The set and order of steps come
from the spectrograph's pipeline, i.e. `Calibrations.default_steps()`:

- **MultiSlit / Echelle:**
  `bias → dark → bpm → slits → arc → tiltimg → wv_calib → tilts → scattlight →
  flats`
- **IFU:**
  `bias → dark → bpm → arc → tiltimg → slits → wv_calib → tilts → align →
  scattlight → flats`

IFU adds an `align` step and orders `slits` after `arc`/`tiltimg`; the button row
must reflect the active spectrograph's pipeline rather than a fixed list.

**`bpm` is not shown.** The bad-pixel-mask step runs internally as part of the
pipeline but produces no standalone output file or QA (its `step_frame_map` entry
is `None`); per the design decision it is **omitted from the button row** (it is
still executed by `pypeit_run_to_calibstep` as a preceding step when needed).

**Button color coding.** Each button's color encodes the step's state, derived
from the `(required, status)` fields in `state.py` together with whether the step
is part of the spectrograph's `default_steps()`. This uses the **Dashboard-wide
status palette** defined in Initialization (the two views are now unified):

| Condition                            | State source                         | Color (glyph) |
|--------------------------------------|--------------------------------------|---------------|
| Step not used by this spectrograph   | step ∉ `default_steps()`             | dimmed light grey #BDBDBD (–), disabled |
| Not required                         | `required = False`                   | grey #9E9E9E (–) |
| Required but not yet generated       | `required = True`, `status = undone` | white #FFFFFF, outlined (○) |
| Currently being generated            | `status = running`                   | orange #EF6C00, subtle pulse (⏳) |
| Generated successfully               | `status ∈ {success, complete}`       | green #2E7D32 (✓) |
| Failed to generate                   | `status = fail`                      | red #C62828 (✗) |

As in the state table (R10), color is **paired with a glyph/label** so status
never depends on color alone.

**Selected step.** The currently selected step (whose detail panel is shown) is
indicated with a thick, high-contrast **magenta ring** (`#D81B60`) plus a
pointer toward the detail panel — chosen because it reads clearly against the
green/grey/white/orange status fills (a blue ring was hard to see).

### The detail panel (on selecting a step)

Clicking a step button selects it and opens a **detail panel** for that step. Per
the user's vision, the panel provides:

- **Metrics.** The step's quality metrics (this is one of the "other windows"
  where metrics belong, per the Initialization decision). From `state.py` today:
  `bias` → `mean`, `std`; `slits` → `nslits` + per-slit `center`; `wv_calib` →
  per-slit `rms`; `tilts` → per-slit `rms`; `flats` → `types`. The panel surfaces
  whatever metrics the state exposes, so this set will **grow as `state.py`
  gains metrics**. Metrics may be threshold-colored using the workflow doc's
  rules of thumb (e.g. wavelength / tilt `rms` < 0.1 px = good).
- **Inspect the output.** A control to launch the appropriate viewer for that
  calibration's processed output file (see the per-type table below): a dedicated
  `pypeit_chk_*` script where one exists (e.g. `slits`, `wv_calib`, `tilts`,
  `flats`), otherwise plain **`ginga`** for the processed image (e.g. `bias`,
  `dark`, `arc`, `tiltimg`). The user can switch between steps and view each.
- **Input files.** The list of `input_files` (raw or processed) used to build the
  calibration (available in `state.py`); the user can **select an input file and
  view it** — **raw** frames via **`pypeit_view_fits`**, **processed**
  calibration frames via **`ginga`**.
- **QA files.** Any related QA PNGs for the step, viewable inline.  For calibrations with many PNG files, a scrollable list of the PNG files should be provided.
- **(Re)generate.** When PypeIt is **not** running, a control to launch
  `pypeit_run_to_calibstep` for the selected step (see below). Disabled while a
  reduction is running.
- **Per-slit / per-order drill-down.** For `slits`, `wv_calib`, and `tilts`, a
  sub-table of per-slit/order rows (each with its own `status` and metric),
  essential for MOS and echelle where there are many slits/orders.

### Backend: `pypeit_run_to_calibstep`

The (re)generate control wraps the existing
`PypeIt/pypeit/scripts/run_to_calibstep.py`:

- **Invocation:**
  `pypeit_run_to_calibstep <pypeit_file> <step> --calib_group <id> --det <det>`.
  Since this view already fixes the group and detector via the drop-downs, the
  dashboard **always** supplies `--calib_group` and `--det` and never uses the
  script's alternative `--science_frame` form.
- **Behavior:** it instantiates `PypeIt(..., calib_only=True)` and calls
  `pypeit_steps.calib_one(..., stop_at_step=<step>)`, which runs
  `Calibrations.run_the_steps(stop_at_step=<step>)` — i.e. it generates the
  requested step **and every step that precedes it** in dependency order, then
  rebuilds the QA HTML. This matches the user's "generate it (and any
  calibrations that precede it)".
- **Valid steps:** `align, arc, bias, bpm, dark, flats, scattlight, slits,
  tiltimg, tilts, wv_calib`.
- **Reuse, don't reinvent** (principle #2): the dashboard should call this
  existing machinery (and `step_frame_map`, `calib_one`, `run_the_steps`) rather
  than reimplement calibration logic. After a run completes, the view refreshes
  state and re-colors the buttons.

### Per-calibration-type specifics

The `step_frame_map` in `calibrations.py` ties each step to its input frametype
and output `DataContainer`. The table below combines that with the inspection
tools from `pypeit_workflow.md`:

"Inspect output" is the viewer for the **processed** output file: a dedicated
`pypeit_chk_*` script where one exists, otherwise plain `ginga`. (Raw input
frames are viewed separately via `pypeit_view_fits`; see the detail panel.)

| Step        | Input frametype | Output (example)            | Inspect output with | Key metric(s) |
|-------------|-----------------|-----------------------------|---------------------|---------------|
| `bias`      | bias            | `Bias_*.fits`               | `ginga` | `mean`, `std` |
| `dark`      | dark            | `Dark_*.fits`               | `ginga` | — |
| `arc`       | arc             | `Arc_*.fits`                | `ginga` | — |
| `tiltimg`   | tilt            | `Tiltimg_*.fits`            | `ginga` | — |
| `slits`     | trace           | `Edges_*.fits.gz`, `Slits_*.fits.gz` | `pypeit_chk_edges` | `nslits`, per-slit `center` |
| `wv_calib`  | arc             | `WaveCalib_*.fits`          | `pypeit_chk_wavecalib`, `pypeit_show_wvcalib` | per-slit `rms` |
| `tilts`     | tilt            | `Tilts_*.fits.gz`           | `pypeit_chk_tilts` | per-slit `rms` |
| `scattlight`| scattlight      | `ScatteredLight_*.fits`     | `pypeit_chk_scattlight` | — |
| `flats`     | illumflat/pixelflat | `Flat_*.fits`           | `pypeit_chk_flats`, `pypeit_show_pixflat` | `types` |
| `align`     | align (IFU)     | `Alignment_*.fits`          | `pypeit_chk_alignments` | — |

(`bpm` is not listed: it produces no output file and is omitted from the button
row, as noted above.)

### Requirements

#### Pending

*General (view-wide):*

- **C1.** Provide a **Calibrations** tab/view, distinct from the Initialization
  state table, scoped to one `(calibration group, detector/mosaic)` at a time.
- **C2.** Two **drop-down** selectors (calibration group; detector / mosaic),
  populated from the state's `(calib_id, det)` pairs; changing either re-renders
  the view.  
- **C3.** A row of **step buttons**, one per calibration step, ordered
  **left→right by dependency** per the spectrograph's `default_steps()`
  (path-aware: MultiSlit/Echelle vs IFU), with connectors conveying precedence.
- **C4.** **Color-code** each button per the **Dashboard-wide status palette**
  (dimmed = not used; grey = not required; white = required-not-generated;
  orange = running; green = success; red = fail), **paired with a glyph/label**
  for accessibility (consistent with R10; same palette as Initialization).
- **C5.** Clicking a button **selects** the step and opens its **detail panel**.
- **C6.** The detail panel shows the step's **metrics** (whatever `state.py`
  exposes — the set grows as the state model does), with optional threshold
  coloring.
- **C7.** The detail panel shows the input files as a **scrollable list** (one
  row per input frame; `bias`/`flats` have many) and lets the user select one
  and **view it** — raw frames via `pypeit_view_fits`, processed frames via
  `ginga`.
- **C8.** The detail panel exposes the **output file** and a control to launch
  the **appropriate viewer** for that calibration type — a `pypeit_chk_*` script
  where one exists, otherwise plain `ginga` (per-type table above).
- **C9.** The detail panel shows related **QA files** as **clickable entries**;
  clicking one opens a **full view** of the PNG.
- **C10.** When PypeIt is **not running**, provide a control to launch
  `pypeit_run_to_calibstep` for the selected step (passing `--calib_group` and
  `--det`), which generates that step **and all preceding steps**; disable it
  while a reduction is running. The control uses a distinct **action color**
  (not a status color, so it is not mistaken for "success").
- **C11.** Support **per-slit / per-order** drill-down for `slits`, `wv_calib`,
  and `tilts`, presented as a **scrollable list** (one row per slit/order; MOS /
  echelle may have many).
- **C12.** Be **path-aware**: show only the steps used by the active
  spectrograph, in its pipeline's order (including IFU's `align`).
- **C13.** **Omit `bpm`** from the button row — it runs internally (and as a
  preceding step for `pypeit_run_to_calibstep`) but has no standalone output
  file/QA to show.
- **C14.** **Reuse existing PypeIt machinery** (`run_to_calibstep`, `calib_one`,
  `run_the_steps`, `step_frame_map`, the `pypeit_chk_*` scripts) rather than
  reimplement calibration logic (principle #2).
- **C15.** Launching a (re)generation must be **observable and non-destructive**
  (principle #5): surface progress/log, and **refresh the state and button
  colors** when it completes.
- **C16.** Prioritize **readability**: clearly labeled, adequately sized buttons;
  dependency order obvious at a glance; using the unified Dashboard-wide status
  palette.

#### Implemented

*(None yet — design phase.)*

#### Deferred

- **CD1.** Launching/queuing calibration (re)generation **while a reduction is
  already running** — out of scope; the (re)generate control is enabled only when
  PypeIt is idle (C10).
- **CD2.** Editing reduction **parameters** from the panel before re-running a
  step — these must be done in the PypeIt file, not the Dashboard.
- **CD3.** **Side-by-side comparison** of a calibration across groups/detectors
  (e.g. diffing two `WaveCalib`) — possible future enhancement.

---

## Logs

*(Design-document work is logged in the parent `dashboard_design.md` Logs
section.)*
