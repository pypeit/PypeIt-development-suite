---
name: run-dev-suite
description: Run the PypeIt development suite correctly — set up the environment, pick test types and instrument/setup subsets, run in parallel, and produce a report. Use when the user wants to validate PypeIt against real data, run regression tests for a setup, or reproduce a dev-suite failure.
---

# Run the PypeIt development suite

This skill operates on **this** repo (`PypeIt-development-suite/`).

## Reference material

- `README.rst` (this repo) — full documentation of `pypeit_test`.
- `pypeit_test` — the runner (wraps `test_scripts/test_main.py`).
- `test_scripts/test_setups.py` — the catalog of instruments/setups/tests.

## Preconditions

1. `run_pypeit` must be on `PATH` (i.e. PypeIt installed, e.g.
   `pip install -e ".[dev]"` in the sibling `PypeIt/` repo). The runner errors
   out if it is missing.
2. `PYPEIT_DEV` must point to the top level of this repo.
3. `RAW_DATA/` must be present (it is git-ignored; usually a symlink to the
   Google-Drive-synced data). `unit`/`reduce`/`afterburn`/`ql`/`vet` need it;
   `pypeit_tests` does not.
4. For GUI/headless runs: `source source_headless_test.sh` (sets
   `QT_QPA_PLATFORM=offscreen`).

## Running

```console
cd $PYPEIT_DEV

# Everything, in order (pypeit_tests, unit, reduce, afterburn, ql, vet)
./pypeit_test all

# Specific test types
./pypeit_test reduce afterburn ql

# Restrict to instruments / setups (prefer a narrow subset while iterating)
./pypeit_test reduce -i shane_kast_blue
./pypeit_test reduce -s shane_kast_blue/600_4310_d55

# List all supported instruments/setups
./pypeit_test list

# Parallelism (watch memory — see the table in README.rst), report, perf, coverage
./pypeit_test -t 2 all
./pypeit_test all -r test_report.txt
./pypeit_test all --csv performance.csv
./pypeit_test all --coverage coverage_report.txt
```

The pytest phases (`pypeit_tests`, `unit`, `vet`) can also be run directly:

```console
pytest unit_tests
pytest vet_tests --redux_out $PYPEIT_DEV/REDUX_OUT
```

## Guidance

- **Iterate narrow, then go wide.** While debugging, target a single
  `instrument/setup`; the full suite takes 10+ hours and tens of GiB of RAM.
- `reduce` produces the outputs that `afterburn`/`vet` consume — run them in
  order (or use `all`).
- Outputs land in `REDUX_OUT/` (git-ignored).
- Use `-p` (prep only) to validate a setup without running the full reduction,
  and `-m` to ignore cached calibrations.

## Verify / interpret

- Read the test summary (and `-r` report) for PASSED/FAILED per setup.
- For a failure, open the per-setup `*.test.log` in `REDUX_OUT/...` and switch
  to the `diagnose-reduction` skill.
