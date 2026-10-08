---
name: new-script
description: Add a new command-line script to PypeIt under pypeit/scripts/ using the ScriptBase pattern and register its entry point. Use when the user wants a new pypeit_* CLI tool or executable.
---

# Add a new PypeIt command-line script

Lives in `PypeIt-development-suite/`; operates on the sibling `PypeIt/` repo.

## Reference material

- Guide: `PypeIt/doc/dev/new_script.rst`
- Base class: `PypeIt/pypeit/scripts/scriptbase.py`
- Examples: any file in `PypeIt/pypeit/scripts/` (e.g. `chk_for_calibs.py`,
  `flux_calib.py`, `run_pypeit.py` for a custom name).

## Procedure

Work in `PypeIt/` on a branch off `develop`.

1. **Create one class per file** in `pypeit/scripts/{module}.py`, subclassing
   `scriptbase.ScriptBase`. The class has no `__init__`. Implement:

   ```python
   from pypeit.scripts import scriptbase

   class NewScript(scriptbase.ScriptBase):

       @classmethod
       def get_parser(cls, width=None):
           parser = super().get_parser(description='A new PypeIt script', width=width)
           # parser.add_argument(...)
           return parser

       @staticmethod
       def main(args):
           # primary work; return value optional
           ...
   ```

2. **Naming**: the default executable name is `pypeit_{module}` (from the file
   name). Override the classmethod `name()` only if you need a different name
   (e.g. `run_pypeit`; see `run_pypeit.py`).

3. **Register the module** in `pypeit/scripts/__init__.py` `__all__`
   (alphabetical).

4. **Register the entry point** in `PypeIt/pyproject.toml` under
   `[options.entry_points]` (the `[project.scripts]`/console-scripts group),
   matching the existing format:

   ```ini
   pypeit_new_script = "pypeit.scripts.new_script:NewScript.entry_point"
   ```

   The base class provides `entry_point` used by `pip` at install time.

## Verify

- Reinstall if needed (`pip install -e .` in `PypeIt/`) so the entry point is
  created.
- `pypeit_new_script -h` prints the parser help.
- If the script help is auto-documented, rebuild docs (see the `build-docs`
  skill — `doc/scripts/write_script_help.py` collects script help).
