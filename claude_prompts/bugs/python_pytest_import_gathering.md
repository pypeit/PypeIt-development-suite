# Bug: `test_flux_calib` passes/fails depending on which test files are collected

- **Date**: 2026-07-23
- **Authors**: Kyle Westfall; Claude (Sonnet 5)
- **PypeIt version**: 2.0.2.dev782+g9abd9367a
- **Branch / commit**: `develop` @ `01bbb4de7775954814942d475b4de93f4e93d610`

## Symptom

From within `pypeit/tests`:

- `pytest . -W ignore` — all tests pass.
- `pytest . -W ignore --ignore=test_runpypeit.py` — `test_fluxspec.py::test_flux_calib`
  fails with:

  ```
  AttributeError: module 'pypeit.scripts' has no attribute 'flux_calib'
  ```

The test's pass/fail outcome depends on whether an unrelated test file
(`test_runpypeit.py`) is included in the same `pytest` invocation, not on
anything the test itself does.

## Initial hypothesis (ruled out)

The wildcard import `from pypeit.scripts import *` in `pypeit/scripts/util.py`
was suspected as the source of namespace pollution. This was **ruled out**:
`pypeit/scripts/util.py` is never imported anywhere during a test run. Its
only consumer in the repo is `doc/scripts/write_script_help.py`, used at
documentation-build time. Confirmed by grepping the codebase for any
`scripts.util` / `scripts import util` reference outside `util.py` itself.

## Root cause

1. `pypeit/scripts/__init__.py` defines `__all__` as a list of submodule
   *name strings* (e.g. `'flux_calib'`, `'run_pypeit'`, ...), but never
   actually imports those submodules. Plain `import pypeit.scripts` does
   **not** make `scripts.flux_calib` available as an attribute:

   ```python
   import pypeit.scripts as scripts
   hasattr(scripts, 'flux_calib')   # False
   ```

2. `pypeit/tests/test_fluxspec.py` does `from pypeit import scripts` and
   then calls `scripts.flux_calib.FluxCalib.parse_args(...)` — it never
   imports `flux_calib` itself. It relies entirely on something *else*
   having already imported `pypeit.scripts.flux_calib` earlier in the
   process.

3. Python's import system has a well-known side effect: whenever a
   submodule is imported by *any* code (`import pkg.sub`,
   `from pkg.sub import X`, `from pkg import sub`), that submodule is
   bound as an attribute on the parent package object and cached in
   `sys.modules` for the remaining lifetime of the process. This binding
   is global — it doesn't matter which piece of code triggered it.

4. Searching the full test suite, the **only** place that directly
   imports `pypeit.scripts.flux_calib` is
   `pypeit/tests/test_runpypeit.py:16`:
   `from pypeit.scripts.flux_calib import FluxCalib`.

5. `pytest` imports every matched test module during its collection
   phase, *before running any test*, regardless of execution order.
   So merely including `test_runpypeit.py` in a run — even though
   `test_flux_calib` executes before it alphabetically — is enough:
   its module-level import fires at collection time and binds
   `pypeit.scripts.flux_calib` as a side effect. Every other test in that
   same session, including `test_fluxspec.py::test_flux_calib`, then
   benefits from that incidental binding.

6. When `test_runpypeit.py` is excluded (`--ignore=test_runpypeit.py`),
   nothing else in the session ever imports that submodule, the
   attribute is never bound, and `scripts.flux_calib` raises
   `AttributeError`.

Verified directly:

```python
import pypeit.scripts as scripts
hasattr(scripts, 'flux_calib')                     # False
from pypeit.scripts.flux_calib import FluxCalib
hasattr(scripts, 'flux_calib')                     # True — global, process-wide side effect
```

and by reproducing both the isolated failure
(`pytest test_fluxspec.py::test_flux_calib` alone) and the isolated pass
(`pytest test_fluxspec.py::test_flux_calib test_runpypeit.py`).

Grepping the rest of the test suite and package code confirmed this
`from pypeit import scripts` + bare-attribute-access pattern exists
**only** in `test_fluxspec.py`; no other test file or production code
shares the same latent hazard.

## Fix

`test_fluxspec.py` should import what it uses directly, e.g.

```python
from pypeit.scripts import flux_calib
```

(or `from pypeit.scripts.flux_calib import FluxCalib`), instead of
depending on `from pypeit import scripts` plus an attribute that only
exists because some other test module happened to import it first. This
makes the test self-sufficient regardless of what else is collected in
the same `pytest` session.

## Can `pytest` itself catch this class of bug?

Not via a simple flag. The pollution lives in `sys.modules` (a
process-global cache) populated at **collection** time, not at
execution time, so:

- Order-randomization plugins (e.g. `pytest-randomly`) only reorder the
  execution of already-collected items; they don't change which files
  get imported during collection, so they would not have caught this.
- `pytest-xdist` workers and `pytest-forked` (`os.fork()`-based) don't
  help either — xdist workers each perform their own full collection
  of the whole suite by default, and forked children inherit the
  parent's already-populated `sys.modules` snapshot taken *after*
  collection.

The only reliable guard is varying **which files are collected at all**,
i.e. running test files/subsets as separate process invocations — e.g. a
CI step or local loop that runs each test module standalone
(`for f in pypeit/tests/test_*.py; do pytest "$f"; done`), or sharding
the suite across CI jobs by file. PypeIt's current CI
(`.github/workflows/ci_tests.yml`) runs `pytest` once against the whole
`pypeit` tree (per `testpaths` in `pyproject.toml`), so nothing in the
existing setup would expose this kind of bug today.

The durable fix remains removing the implicit coupling (explicit
imports) rather than relying on test-runner configuration to catch it
after the fact.
