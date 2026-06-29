#!/usr/bin/env python
"""
Profile a clean COLD run of ``run_pypeit keck_deimos_600zd_m_6500.pypeit -o`` with
cProfile.

Part of the PypeIt speed-up effort (see
``claude_prompts/speed_up/speed_up_design_prompts.md``). This is the second,
heavier profiling target after shane_kast_blue: keck_deimos 600ZD_M_6500 is a
MultiSlit slitmask reduction that defaults to **4 detector mosaics**
``[(1,5),(2,6),(3,7),(4,8)]`` (the 8 chips paired), so it exercises the full
multi-detector / multi-slit workload.

What it does
------------
1. Cleans the redux directory so the run is cold: removes any prior calibration
   and science products (``Calibrations/``, ``Science/``, ``QA/``,
   ``Intermediate/``), the stale ``*_state.json``, the previous run log, the
   ``*.calib`` association file, and the old ``*_UTC_*.par`` dumps. The
   ``.pypeit`` input and the raw data are left untouched. (For this setup the dir
   is essentially clean already, so this is largely a no-op.)
2. Runs the reduction in-process under ``cProfile`` (single process; PypeIt core
   uses no multiprocessing, so the profile is complete).
3. Saves the binary profile to ``Reports/keck_deimos_600zd_m_6500.prof``, a copy
   of the timestamped run log, and a small ``*.runmeta.json`` with the overall
   wall-clock time.

Run with the pypeit env python::

    /home/xavier/miniconda3/envs/pypeit/bin/python profile_deimos.py

Execution (the reduction itself) is the slow part; expect tens of minutes for the
full 4-mosaic run, plus cProfile overhead.
"""
import os
import sys
import time
import json
import glob
import shutil
import cProfile
from pathlib import Path

REDUX = Path("/home/xavier/Projects/PypeIt/PypeIt-development-suite/"
             "REDUX_OUT/keck_deimos/600ZD_M_6500")
PYPEIT_FILE = "keck_deimos_600zd_m_6500.pypeit"
STEM = "keck_deimos_600zd_m_6500"
REPORTS = Path("/home/xavier/Projects/PypeIt/PypeIt-development-suite/"
               "pypeitdev/speed_up/Reports")
PROF_OUT = REPORTS / f"{STEM}.prof"
META_OUT = REPORTS / f"{STEM}.runmeta.json"
LOG_COPY = REPORTS / f"{STEM}.run.log"


def clean_cold():
    """Remove prior products so the profiled run rebuilds everything."""
    removed = []
    for sub in ("Calibrations", "Science", "QA", "Intermediate"):
        p = REDUX / sub
        if p.is_dir():
            shutil.rmtree(p)
            removed.append(str(p.relative_to(REDUX)) + "/")
    patterns = ["*_state.json", "*.calib", f"{STEM}.log",
                f"{STEM}_UTC_*.par"]
    for pat in patterns:
        for f in glob.glob(str(REDUX / pat)):
            os.remove(f)
            removed.append(os.path.basename(f))
    return removed


def main():
    REPORTS.mkdir(parents=True, exist_ok=True)
    os.chdir(REDUX)

    removed = clean_cold()
    print("Cold-clean removed:")
    for r in removed:
        print(f"  - {r}")
    if not removed:
        print("  (nothing to remove)")

    # Import after chdir so PypeIt picks up the right cwd defaults.
    from pypeit.scripts.run_pypeit import RunPypeIt

    args = RunPypeIt.parse_args([PYPEIT_FILE, "-o"])

    print(f"\nRunning (cold) under cProfile: run_pypeit {PYPEIT_FILE} -o")
    pr = cProfile.Profile()
    t0 = time.perf_counter()
    rc = 0
    pr.enable()
    try:
        RunPypeIt.main(args)
    except SystemExit as e:
        rc = int(e.code) if e.code is not None else 0
    finally:
        pr.disable()
        wall = time.perf_counter() - t0
        pr.dump_stats(str(PROF_OUT))

    print(f"\nWall-clock (cProfile-instrumented): {wall:.2f} s")
    print(f"Profile written to: {PROF_OUT}")

    # Copy the run log (timestamped on this branch) for timeline analysis.
    runlog = REDUX / f"{STEM}.log"
    if runlog.exists():
        shutil.copy(str(runlog), str(LOG_COPY))
        print(f"Run log copied to: {LOG_COPY}")

    meta = {
        "redux_dir": str(REDUX),
        "pypeit_file": PYPEIT_FILE,
        "command": f"run_pypeit {PYPEIT_FILE} -o",
        "cold_run": True,
        "removed_before_run": removed,
        "wall_clock_s": round(wall, 3),
        "return_code": rc,
        "prof": str(PROF_OUT),
        "run_log": str(LOG_COPY),
    }
    with open(META_OUT, "w") as fp:
        json.dump(meta, fp, indent=2)
    print(f"Run metadata written to: {META_OUT}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
