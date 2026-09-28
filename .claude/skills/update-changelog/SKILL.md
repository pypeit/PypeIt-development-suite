---
name: update-changelog
description: Record a PypeIt change in the release notes / changelog consistently for a pull request — the in-development release doc and CHANGES.rst. Use when finishing a feature, fix, or datamodel/parameter change that users should be told about before opening a PR.
---

# Update the PypeIt changelog / release notes

Lives in `PypeIt-development-suite/`; operates on the sibling `PypeIt/` repo.

## Reference material

- Per-release notes: `PypeIt/doc/releases/` — one file per version. The current
  in-development file ends in `dev` (e.g. `2.1.0dev.rst`); **add your entry
  here**.
- Top-level summary: `PypeIt/CHANGES.rst`.
- `PypeIt/doc/whatsnew.rst` (user-facing highlights, if relevant).

## Procedure

Work in `PypeIt/` on your feature branch.

1. **Find the in-development release doc** in `doc/releases/` (the `*dev.rst`
   file). Match the existing section structure — typical sections include
   general improvements, bug fixes, and **Instrument-specific** changes.

2. **Add a concise bullet** describing the user-visible change. Mirror the
   wording/voice of existing bullets. For:
   - a new spectrograph/mode → put it under "Instrument-specific";
   - a datamodel version bump → call out the version change and any
     compatibility impact;
   - a new/changed parameter → name the parameter and its effect;
   - a new script → name the `pypeit_*` executable.

3. **Update `CHANGES.rst`** with a matching short entry if the change is
   significant enough to summarize there.

4. **Cross-check related docs** were updated (parameters, datamodels, scripts —
   see the `build-docs`, `add-parameter`, `new-datacontainer` skills). The
   `PypeIt/CLAUDE.md` requires docs to be updated with each PR.

## Verify

- The release doc renders (rebuild via the `build-docs` skill, at least
  `make htmlonly`).
- The bullet is accurate, concise, and in the correct section.
- No version string was edited by hand unless that is the explicit task
  (versioning is managed via `setuptools_scm`).
