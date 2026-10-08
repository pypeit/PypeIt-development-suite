---
name: add-parameter
description: Add or modify a user-level PypeIt parameter in a PypeItPar parameter set (pypeitpar.py), with correct defaults/options/dtypes/docstrings, and regenerate the parameter documentation. Use when the user wants a new tunable reduction option or to change an existing parameter.
---

# Add or modify a PypeIt parameter

Lives in `PypeIt-development-suite/`; operates on the sibling `PypeIt/` repo.

## Reference material

- Parameter sets: `PypeIt/pypeit/par/pypeitpar.py` — a hierarchy of
  `ParSet` subclasses (e.g. `ReduxPar`, `CalibrationsPar`, `WavelengthSolutionPar`,
  `ProcessImagesPar`, `Coadd1DPar`, …) packaged by `PypeItPar`.
- Base machinery: `PypeIt/pypeit/par/parset.py`.
- Doc generator: `PypeIt/doc/scripts/build_par_rst.py` (writes the parameter
  tables in `doc/pypeit_par.rst`).

## Procedure

Work in `PypeIt/` on a branch off `develop`.

1. **Find the right parameter set.** Parameters are grouped by function. Add the
   parameter to the `ParSet` whose scope matches (e.g. extraction options go in
   the relevant extraction par class). Search `pypeitpar.py` for a sibling
   parameter to match the pattern exactly.

2. **Edit the class `__init__`.** Each `ParSet` subclass defines its parameters
   via parallel lists/dicts. For the new key, add entries to all of:
   - `pars[...]` default value (and add the argument to `__init__` with that
     default),
   - `dtypes[...]` allowed type(s),
   - `descr[...]` a clear one-to-several sentence description (this text is what
     appears in the auto-generated docs — write it well),
   - `options[...]` if the value is restricted to a fixed set,
   - any validation in the class body.
   Then add the key to the `pars` dict construction and the `cls.from_dict`
   `parkeys` list so it round-trips through config files.

3. **Keep `from_dict`/`to_config` consistent** — the parameter must be listed in
   the class's `parkeys` so it parses from a `.pypeit`/`.coadd1d`/etc. file.

4. **Wire it into the algorithm** that consumes the parameter (the par object is
   threaded through the reduction; grep for where the sibling parameter is read).

## Regenerate docs

From `PypeIt/doc`:

```console
python scripts/build_par_rst.py    # or rebuild via the build-docs skill
```

Then rebuild HTML (see the `build-docs` skill).

## Verify

- `python -c "from pypeit.par.pypeitpar import PypeItPar; p=PypeItPar(); print(p.to_config())"`
  shows the new key with its default.
- A `.pypeit` file setting the parameter parses without error.
- Add/extend a unit test in `PypeIt/pypeit/tests/` if behavior changed.
