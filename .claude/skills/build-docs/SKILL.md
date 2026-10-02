---
name: build-docs
description: Rebuild the PypeIt Sphinx documentation, including the auto-generated parameter/datamodel/script tables. Use when the user changes parameters, datamodels, scripts, or docstrings and needs the docs rebuilt, or wants to check the docs build cleanly.
---

# Build the PypeIt documentation

Lives in `PypeIt-development-suite/`; operates on the sibling `PypeIt/` repo
(`PypeIt/doc/`).

## Reference material

- `PypeIt/doc/Makefile` — defines the `apirst`, `html`, `htmlonly`, `htmlnoex`,
  and `clean` targets.
- `PypeIt/doc/scripts/build_*.py` — the generators that produce auto-built
  tables (parameters, datamodels, detectors, bitmasks, dependencies, script
  help, etc.). `make apirst` runs the ones that do not need the dev suite.

## Choosing a target

From `PypeIt/doc`:

- **Full build** (requires internet and `PYPEIT_DEV` pointing at this dev-suite
  repo with `RAW_DATA`):
  ```console
  cd PypeIt/doc
  make clean
  make html
  ```
  `make html` runs `make apirst` (API docs + auto tables) and `make examples`,
  then builds HTML.

- **Limited/offline build** (no dev-suite data, no examples):
  ```console
  cd PypeIt/doc
  make htmlonly      # skip apirst and examples; just rebuild HTML
  # or
  make htmlnoex      # run apirst but skip examples
  ```

## When to rebuild

- After changing parameters → also covered by the `add-parameter` skill
  (`build_par_rst.py`).
- After changing a datamodel → `build_datacontainer_datamodels.py` /
  `build_afterburn_datamodels_rst.py` (run via `make apirst`).
- After adding/changing a script → `write_script_help.py` (run via `make apirst`).

## Verify

- The build finishes without Sphinx errors (warnings about missing dev-suite
  examples are expected with `htmlonly`/`htmlnoex`).
- Open `PypeIt/doc/_build/html/index.html` and check the changed page.
- Per `CLAUDE.md` in `PypeIt/`, docs should be updated with each pull request.
