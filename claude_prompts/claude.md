# Improving the Claude experinece for PypeIt development

## Init

Generate a CLAUDE.md file to guide development of PypeIt.  The main Repository is named PypeIt/ and it has its own CLAUDE.md file.  The development suite is named PypeIt-development-suite/ 

1. Generate a CLAUDE.md file to be placed in PypeIt-development-suite/.  As useful, have it refer to both repositories.

## Skills

Suggested [Claude Code skills](https://docs.claude.com/en/docs/claude-code/skills)
to add for PypeIt development. **All skills live in the
`PypeIt-development-suite/` repository** under `.claude/skills/<name>/SKILL.md`,
regardless of which repository the skill operates on (the dev suite assumes the
main `PypeIt/` repo is a sibling clone). The "lives in" URL is the location in
this repo; the "based on" URL is the existing PypeIt doc/code the skill codifies.

### Skills that operate on the main `PypeIt/` repository

1. **new-spectrograph** — Scaffold, implement, and register a new
   `Spectrograph` subclass (config keys, detector/`meta` definitions, frame
   typing, default `PypeItPar`, processing path), then wire it into
   `pypeit/spectrographs/__init__.py` and add a matching dev-suite setup.
   - Lives in: https://github.com/pypeit/PypeIt-development-suite/tree/develop/.claude/skills/new-spectrograph
   - Based on: https://github.com/pypeit/PypeIt/blob/develop/doc/dev/new_spectrograph.rst
     and https://github.com/pypeit/PypeIt/tree/develop/pypeit/spectrographs

2. **new-script** — Add a new command-line tool under `pypeit/scripts/`
   following the `scriptbase.ScriptBase` pattern and register its entry point.
   - Lives in: https://github.com/pypeit/PypeIt-development-suite/tree/develop/.claude/skills/new-script
   - Based on: https://github.com/pypeit/PypeIt/blob/develop/doc/dev/new_script.rst
     and https://github.com/pypeit/PypeIt/blob/develop/pypeit/scripts/scriptbase.py

3. **add-parameter** — Add or modify a user-level parameter in a `PypeItPar`
   set, with validation/defaults/docstrings, and regenerate the parameter
   documentation tables.
   - Lives in: https://github.com/pypeit/PypeIt-development-suite/tree/develop/.claude/skills/add-parameter
   - Based on: https://github.com/pypeit/PypeIt/blob/develop/pypeit/par/pypeitpar.py
     and https://github.com/pypeit/PypeIt/blob/develop/doc/scripts/build_par_rst.py

4. **new-datacontainer** — Create or extend a `DataContainer` subclass with a
   strict datamodel and correct FITS I/O, then regenerate the datamodel docs.
   - Lives in: https://github.com/pypeit/PypeIt-development-suite/tree/develop/.claude/skills/new-datacontainer
   - Based on: https://github.com/pypeit/PypeIt/blob/develop/pypeit/datamodel.py
     and https://github.com/pypeit/PypeIt/blob/develop/doc/scripts/build_datacontainer_datamodels.py

5. **build-docs** — Rebuild the Sphinx documentation, including the
   auto-generated tables (`doc/scripts/build_*.py`), via `make html`
   (full, needs `PYPEIT_DEV`) or `make htmlonly` (limited/offline).
   - Lives in: https://github.com/pypeit/PypeIt-development-suite/tree/develop/.claude/skills/build-docs
   - Based on: https://github.com/pypeit/PypeIt/tree/develop/doc/scripts

6. **diagnose-reduction** — Triage a failed or suspect `run_pypeit` reduction
   using the log, QA HTML, and the `pypeit_chk_*` inspection scripts
   (edges, flats, tilts, wavecalib, scattlight, noise).
   - Lives in: https://github.com/pypeit/PypeIt-development-suite/tree/develop/.claude/skills/diagnose-reduction
   - Based on: https://github.com/pypeit/PypeIt/tree/develop/pypeit/scripts
     (the `chk_*.py` scripts) and https://pypeit.readthedocs.io/

7. **wavelength-calibration** — Build, validate, and archive wavelength
   solutions: `pypeit_identify`, line lists, and reusable arc templates.
   - Lives in: https://github.com/pypeit/PypeIt-development-suite/tree/develop/.claude/skills/wavelength-calibration
   - Based on: https://github.com/pypeit/PypeIt/blob/develop/pypeit/scripts/identify.py
     and https://github.com/pypeit/PypeIt/tree/develop/pypeit/core/wavecal

8. **update-changelog** — Update `CHANGES.rst` / `doc/releases` and version
   metadata consistently for a pull request.
   - Lives in: https://github.com/pypeit/PypeIt-development-suite/tree/develop/.claude/skills/update-changelog
   - Based on: https://github.com/pypeit/PypeIt/tree/develop/doc/releases

### Skills that operate on the `PypeIt-development-suite/` repository

9. **run-dev-suite** — Run the development suite correctly: ensure
   `run_pypeit` is on `PATH` and `PYPEIT_DEV` is set, select test types and
   `-i`/`-s` instrument/setup subsets, and produce a report.
   - Lives in: https://github.com/pypeit/PypeIt-development-suite/tree/develop/.claude/skills/run-dev-suite
   - Based on: https://github.com/pypeit/PypeIt-development-suite/blob/develop/README.rst
     and https://github.com/pypeit/PypeIt-development-suite/blob/develop/pypeit_test

10. **add-devsuite-setup** — Add a new instrument/setup or test to the dev
    suite: stage `RAW_DATA/<instrument>/<setup>/` plus the right input file
    (`pypeit_files/`, `coadd1d_files/`, `fluxing_files/`, …) and register it in
    `test_scripts/test_setups.py` (and `pypeit_tests.py` for a new test type).
    - Lives in: https://github.com/pypeit/PypeIt-development-suite/tree/develop/.claude/skills/add-devsuite-setup
    - Based on: https://github.com/pypeit/PypeIt-development-suite/blob/develop/test_scripts/test_setups.py

## Prompts

1. Perform the 1st step under Init.
2. Given your understanding of the code base, provide a list of suggested skills to add for Claude.  Provide the list in the Skills section above and include URLs to their locations on GitHub.
3. Modify your Skill suggestions to all be located in the PypeIt-development-suite/ repository.  Then proceed to generate each of these.