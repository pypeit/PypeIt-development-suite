---
name: new-datacontainer
description: Create or extend a PypeIt DataContainer subclass with a strict datamodel and correct FITS I/O, then regenerate the datamodel documentation. Use when the user wants a new on-disk data product (FITS output) or to change the datamodel of an existing one.
---

# Create or modify a PypeIt DataContainer

Lives in `PypeIt-development-suite/`; operates on the sibling `PypeIt/` repo.

## Reference material

- Base class & docs: `PypeIt/pypeit/datamodel.py` — read the module docstring;
  it documents the `datamodel` dict, versioning, and I/O contract.
- Examples: `pypeit/spec2dobj.py`, `pypeit/wavetilts.py`, `pypeit/flatfield.py`,
  `pypeit/edgetrace.py`, `pypeit/specobj.py` — all `DataContainer` subclasses.
- Doc generator: `PypeIt/doc/scripts/build_datacontainer_datamodels.py`
  (and `build_afterburn_datamodels_rst.py`, `build_specobj_rst.py`).

## Procedure

Work in `PypeIt/` on a branch off `develop`.

1. **Define the datamodel.** Subclass `DataContainer` and set the class-level
   `datamodel` dict: each entry maps an attribute name to `{'otype': ...,
   'atype': ..., 'descr': ...}` (and `'tab_overwrite'`/`'pixel'` etc. as
   appropriate). The `descr` strings are what appear in the generated docs.

2. **Set `version`.** Bump `version` whenever you change the datamodel of an
   existing container — old files are read against the recorded version. Note
   the change for users (see the `update-changelog` skill).

3. **Implement I/O hooks as needed:**
   - `_init_internals()` — declare non-datamodel attributes.
   - `_bundle()` — group attributes into HDUs for writing.
   - `_parse()` (classmethod) — reconstruct from HDUs on read.
   Only override these when the default behavior is insufficient; copy the
   pattern from the closest existing container.

4. **`to_file` / `from_file`** come from the base class; use the standard
   `hdu_prefix`/`output_to_disk` mechanisms rather than custom FITS code.

## Regenerate docs

From `PypeIt/doc`:

```console
python scripts/build_datacontainer_datamodels.py
# afterburn products: python scripts/build_afterburn_datamodels_rst.py
```

Then rebuild HTML (see the `build-docs` skill).

## Verify

- Round-trip: instantiate, `to_file(...)`, `from_file(...)`, and confirm
  equality of all datamodel attributes.
- For a changed existing container, confirm an old file still reads (version
  handling) or that a clear error/upgrade path exists.
- Add a unit test under `PypeIt/pypeit/tests/`.
