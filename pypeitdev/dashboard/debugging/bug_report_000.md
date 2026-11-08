# Bug Report 000

## Bug Description

When running the Dashboard on a `keck_mosfire` dataset, I was trying to create the Flats using the (Re)Build button.  The process crashed but
the Dashboard did not indicate this and the step remained orange and the Activity Bar indicated the process was running.

## Desired outcome

When the process crashes, the Dashboard should indicate this and the step should turn red. The Activity Bar should indicate the process has failed.

## Reproduce

You can reproduce this by running the Dashboard on the data in /home/xavier/Projects/PypeIt/Redux/Keck/MOSFIRE/Fynbo/test_A/ using the (Re)Build button on the Flats step.  I have artificially crashed the process by adding a raise PypeItError("This is only a test") to the flatfield.py file in PypeIt/pypeit/flatfield.py.

## Report

**Status:** Fixed (Dashboard v1.4.0).

### Reproduce

Confirmed on the reporter's own data
(`/home/xavier/Projects/PypeIt/Redux/Keck/MOSFIRE/Fynbo/test_A/`).  The crashed
(Re)Build left `keck_mosfire_A_state.json` with:

```
current_step: flats
flats -> status = 'running'      # orange, never updated
```

i.e. the flats step was stuck on `running` (orange) even though the process had
already died.

### Diagnose (root cause)

There were **two** distinct defects, at two layers.

1. **Pipeline (the step stays orange).**  In
   `pypeit/calibrations.py::Calibrations.run_the_steps`, each step is recorded as
   `running` (and written to `*_state.json`) *before* `get_<step>()` is called:

   ```python
   self.state.safe_update_calib(step, ..., 'status', 'running')
   self.state.safe_write()
   getattr(self, f'get_{step}')(force=force)   # <- not wrapped
   # ... the 'fail' write below is only reached via self.success
   ```

   The only path that wrote `fail` was the `if not self.success:` branch **after**
   the call returns.  When `get_flats()` *raised* (the artificial
   `PypeItError` in `flatfield.py`), the exception propagated straight out of
   `run_the_steps` → `calib_one` → the `pypeit_run_to_calibstep` script, which
   exited non-zero.  The `fail` write was skipped, so `*_state.json` was left with
   `flats = running`.  The Dashboard re-reads the state on completion, sees
   `running`, and colors the step orange — forever.  (This affects any client of
   the state: `pypeit_status`, too.)

2. **Dashboard (the Build channel says "running" / "Idle", never "failed").**
   The Build channel of the activity bar is driven by *two* competing sources:
   the `Launcher` (which reports "exited with code N" on finish) and the
   live monitor in `view/main_window.py` (`_on_state_changed` /
   `_on_lock_changed`, which write "Monitoring run — updating live…" / "Idle").
   After a crash the reduction `.log` mtime is still "recent"
   (`RunLock.ACTIVE_WINDOW_S = 10 s`), so `is_locked()` stays `True` for several
   more polls: the next `_on_state_changed` tick overwrites the launcher's
   message back to "Monitoring… (busy)", and when the lock finally releases,
   `_on_lock_changed(False)` sets it to "Idle".  The failure is never surfaced.

### Fix

- **Root cause, pipeline** (`pypeit/calibrations.py`): wrapped the
  `get_<step>()` call in `run_the_steps` in `try/except` that marks the step
  `fail` (and `safe_write`s, except in status-only mode) before **re-raising**,
  so the run still aborts with a non-zero exit but the state reflects the
  failure.  The step now turns **red** on the next refresh for every client.
- **Dashboard** (`pypeit/dashboard/view/main_window.py`): added a sticky
  `_last_run_failed` flag, set in `_on_run_finished(code)` when `code != 0`,
  which (a) writes "(Re)Build failed (exit code N) — see the log." on the Build
  channel, (b) suppresses the live monitor's "Monitoring…"/"Idle" overwrites
  while it is set, and (c) is cleared when a new run starts (`lockChanged →
  True`).

### Verify

- Regression tests (both fail before the fix, pass after):
  `pypeit/tests/test_state.py::test_run_the_steps_marks_fail_on_exception`
  (pipeline) and
  `pypeit/tests/test_dashboard.py::test_main_window_failed_run_marks_build_channel`
  (Dashboard).  `pytest test_dashboard.py test_state.py` → **100 passed**
  (offscreen).
- End-to-end on the reporter's `keck_mosfire` data: re-ran
  `pypeit_run_to_calibstep keck_mosfire_A.pypeit flats --calib_group 0 --det 1`.
  The run now exits **1** and `keck_mosfire_A_state.json` records
  `flats → 'fail'` (which `palette.step_style` renders red, `#C62828`, glyph ✗).
- Docs rebuilt the ReadTheDocs way (`sphinx-build … -W --keep-going`):
  **exit 0, 0 warnings**.

> Note: `pypeit/flatfield.py` still carries the reporter's artificial
> `raise PypeItError("This is only a test")`; revert it before normal use.