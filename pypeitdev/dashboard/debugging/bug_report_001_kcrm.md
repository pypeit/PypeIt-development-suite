# Bug Report KCRM

## Preliminary bugs and quick fixes
Initial attempts for a calibration-only run on a `keck_kcrm` failed at the step of generating the `Alignment` frame,
and then later at creating the flat frames. 

1. **Minor fix to get_align() function in calibration.py** When trying to find the calibrations via `self.find_calibrations()`, `self.raw_files` is returned, but the next immediate conditional checks the length of `raw_files`. The first fix is just to make these consistent with each other.

2. **Adding an align_state() function in calibration.py** Consistent with other functions like `tiltimg_state()` and `arc_state()`, I have added an `align_state()`:

   ```python
   def align_state(self):
        self.base_state('align', self.alignments)
   ```
This enables the status for the alignment frames to update in the dashboard (i.e., green checkmark). 
Prior to this, the frame would be created, but the Dashboard would only show that the frame was being created.

3. **Removing the raise PypeItError("This is only a test") line in flatfield.py** This allowed the pipeline to the end. That is,
the reduction of the science frame(s).

## Current Status
I am able to complete a reduction run of a science frame, and I am also able to re-run individual calibration steps.