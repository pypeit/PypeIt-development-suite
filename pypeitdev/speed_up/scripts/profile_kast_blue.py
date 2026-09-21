#!/usr/bin/env python
"""
Profile a clean COLD run of ``run_pypeit shane_kast_blue_A.pypeit -o`` with
cProfile.

Part of the PypeIt speed-up effort (see
``claude_prompts/speed_up/speed_up_design_prompts.md``).

What it does
------------
1. Cleans the redux directory so the run is cold: removes the calibration and
   science products (``Calibrations/``, ``Science/``, ``QA/``,
   ``Intermediate/``), the stale ``*_state.json``, the previous run log, the
   ``*.calib`` association file, and the old ``*_UTC_*.par`` dumps. The
   ``.pypeit`` input and the raw data are left untouched.
2. Runs the reduction in-process under ``cProfile`` (single process; the
   shane_kast_blue path uses no multiprocessing, so the profile is complete).
3. Saves the binary profile to ``Reports/shane_kast_blue_A.prof`` and a small
   ``Reports/shane_kast_blue_A.runmeta.json`` with the overall wall-clock time.

Run with the pypeit env python::

    /home/xavier/miniconda3/envs/pypeit/bin/python profile_kast_blue.py

Execution (the reduction itself) is the slow part; expect ~1-2 min.
"""
import os
import sys
import time
import json
import glob
import shutil
import argparse
import cProfile
from pathlib import Path

REDUX = Path("/home/xavier/Projects/PypeIt/PypeIt-development-suite/"
             "REDUX_OUT/shane_kast_blue/600_4310_d55/shane_kast_blue_A")
PYPEIT_FILE = "shane_kast_blue_A.pypeit"
REPORTS = Path("/home/xavier/Projects/PypeIt/PypeIt-development-suite/"
               "pypeitdev/speed_up/Reports")
STEM = "shane_kast_blue_A"


def clean_cold():
    """Remove prior products so the profiled run rebuilds everything."""
    removed = []
    for sub in ("Calibrations", "Science", "QA", "Intermediate"):
        p = REDUX / sub
        if p.is_dir():
            shutil.rmtree(p)
            removed.append(str(p.relative_to(REDUX)) + "/")
    patterns = ["*_state.json", "*.calib", "shane_kast_blue_A.log",
                "shane_kast_blue_A_UTC_*.par"]
    for pat in patterns:
        for f in glob.glob(str(REDUX / pat)):
            os.remove(f)
            removed.append(os.path.basename(f))
    return removed


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ncpu", type=int, default=None,
                    help="Pass --ncpu N through to run_pypeit and suffix the "
                         "output artifact stems with _ncpu{N} so existing "
                         "baselines are not overwritten.")
    cli = ap.parse_args()

    stem = STEM + (f"_ncpu{cli.ncpu}" if cli.ncpu is not None else "")
    prof_out = REPORTS / f"{stem}.prof"
    meta_out = REPORTS / f"{stem}.runmeta.json"
    log_copy = REPORTS / f"{stem}.run.log"

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

    argv = [PYPEIT_FILE, "-o"]
    if cli.ncpu is not None:
        argv += ["--ncpu", str(cli.ncpu)]
    args = RunPypeIt.parse_args(argv)

    print(f"\nRunning (cold) under cProfile: run_pypeit {' '.join(argv)}")
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
        pr.dump_stats(str(prof_out))

    print(f"\nWall-clock (cProfile-instrumented): {wall:.2f} s")
    print(f"Profile written to: {prof_out}")

    # Copy the run log (timestamped on this branch) for timeline analysis.
    runlog = REDUX / "shane_kast_blue_A.log"
    if runlog.exists():
        shutil.copy(str(runlog), str(log_copy))
        print(f"Run log copied to: {log_copy}")

    meta = {
        "redux_dir": str(REDUX),
        "pypeit_file": PYPEIT_FILE,
        "command": f"run_pypeit {' '.join(argv)}",
        "ncpu": cli.ncpu,
        "cold_run": True,
        "removed_before_run": removed,
        "wall_clock_s": round(wall, 3),
        "return_code": rc,
        "prof": str(prof_out),
        "run_log": str(log_copy),
    }
    with open(meta_out, "w") as fp:
        json.dump(meta, fp, indent=2)
    print(f"Run metadata written to: {meta_out}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
