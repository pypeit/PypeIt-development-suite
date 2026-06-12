# PypeIt Dashboard — Design Document

**Version:** 0.1
**Date:** 2026-06-12
**Author:** JXP and Claude

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
- *(subsequent phases: Reduction control, Inspection/QA, ... — to be added.)*

---

## Logs

*(Design-document work is logged in the parent `dashboard_design.md` Logs
section.)*
