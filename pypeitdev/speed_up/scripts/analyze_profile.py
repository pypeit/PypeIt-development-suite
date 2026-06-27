#!/usr/bin/env python
"""
Analyze the cProfile output + run log produced by ``profile_kast_blue.py``.

Produces two machine-readable text artifacts in ``Reports/``:

- ``shane_kast_blue_A.pstats.txt`` — top functions by cumulative time and by
  internal (tottime) time, plus a few targeted callers.
- ``shane_kast_blue_A.timeline.txt`` — per-phase wall-clock timeline parsed from
  the timestamped run log, mapped onto the ``pypeit_workflow.md`` phases.

Run with the pypeit env python::

    /home/xavier/miniconda3/envs/pypeit/bin/python analyze_profile.py
"""
import io
import re
import json
import pstats
import datetime
from pathlib import Path

REPORTS = Path("/home/xavier/Projects/PypeIt/PypeIt-development-suite/"
               "pypeitdev/speed_up/Reports")
PROF = REPORTS / "shane_kast_blue_A.prof"
LOG = REPORTS / "shane_kast_blue_A.run.log"
META = REPORTS / "shane_kast_blue_A.runmeta.json"
PSTATS_OUT = REPORTS / "shane_kast_blue_A.pstats.txt"
TIMELINE_OUT = REPORTS / "shane_kast_blue_A.timeline.txt"

# ---------------------------------------------------------------------------
# 1. pstats tables
# ---------------------------------------------------------------------------


def write_pstats():
    buf = io.StringIO()
    st = pstats.Stats(str(PROF), stream=buf)

    buf.write("=" * 78 + "\n")
    buf.write("TOP 40 BY CUMULATIVE TIME\n")
    buf.write("=" * 78 + "\n")
    st.sort_stats("cumulative").print_stats(40)

    buf.write("\n" + "=" * 78 + "\n")
    buf.write("TOP 40 BY INTERNAL (tottime) TIME\n")
    buf.write("=" * 78 + "\n")
    st.sort_stats("tottime").print_stats(40)

    PSTATS_OUT.write_text(buf.getvalue())
    print(f"Wrote {PSTATS_OUT}")


# ---------------------------------------------------------------------------
# 2. Phase timeline from the timestamped run log
# ---------------------------------------------------------------------------

# FileFormatter on this branch:
#   "%(levelname)8s | %(asctime)s | %(filename)s:%(funcName)s:%(lineno)s | msg"
# asctime datefmt = "%Y-%m-%d %H:%M:%S"
LOGLINE = re.compile(
    r"^\s*(?P<level>[A-Z]+)\s*\|\s*"
    r"(?P<ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\s*\|\s*"
    r"(?P<file>[^:]+):(?P<func>[^:]+):(?P<line>\d+)\s*\|\s*"
    r"(?P<msg>.*)$"
)

# Phase boundary markers: (regex on "file:func", or on message) -> phase label.
# Ordered; first match wins for labeling a transition point.
MARKERS = [
    # ---- Initialization ----
    (r"inputfiles\.py:from_file", "Initialization: load .pypeit"),
    (r"metadata\.py:_build", "Initialization: build metadata"),
    (r"pypeit\.py:__init__|x_pypeit\.py:__init__", "Initialization: PypeIt setup"),
    # ---- Reduction driver ----
    (r"reduce_all", "Reduction: reduce_all begins"),
    # ---- Calibrations (per calib_one) ----
    (r"calib_one", "Calibrations: build/load"),
    (r"buildimage|combine", "Calibrations: image combine"),
    (r"edgetrace|TraceImage|trace", "Calibrations: slit edges"),
    (r"wavecalib|autoid|WaveCalib", "Calibrations: wavelength"),
    (r"wavetilts|BuildWaveTilts|fit2tiltimg", "Calibrations: tilts"),
    (r"flatfield|FlatImages|flat", "Calibrations: flat"),
    (r"run_the_steps", "Calibrations: complete"),
    # ---- Object find / extract ----
    (r"find_objects", "Reduce: object finding"),
    (r"global_skysub|skysub", "Reduce: sky subtraction"),
    (r"extract_one|extraction\.py", "Reduce: extraction"),
    (r"local_skysub_extract", "Reduce: local sky + extract"),
    (r"save|spec2dobj|onespec|all_spec2d|save_exposure", "Reduce: write products"),
    # ---- QA ----
    (r"qa\.py|gen_qa|build_qa|html", "QA"),
]
COMPILED = [(re.compile(p), lab) for p, lab in MARKERS]


def label_for(file_, func):
    key = f"{file_}:{func}"
    for rx, lab in COMPILED:
        if rx.search(key):
            return lab
    return None


def parse_timeline():
    if not LOG.exists():
        print(f"!! No run log at {LOG}; skipping timeline.")
        return
    events = []  # (datetime, label, raw "file:func", msg)
    first_ts = last_ts = None
    for line in LOG.read_text(errors="replace").splitlines():
        m = LOGLINE.match(line)
        if not m:
            continue
        ts = datetime.datetime.strptime(m["ts"], "%Y-%m-%d %H:%M:%S")
        if first_ts is None:
            first_ts = ts
        last_ts = ts
        lab = label_for(m["file"], m["func"])
        if lab is not None:
            events.append((ts, lab, f"{m['file']}:{m['func']}:{m['line']}", m["msg"]))

    buf = io.StringIO()
    buf.write("=" * 78 + "\n")
    buf.write("PHASE TIMELINE (from timestamped run log; 1-s log resolution)\n")
    buf.write("=" * 78 + "\n")
    if first_ts and last_ts:
        buf.write(f"Log span: {first_ts}  ->  {last_ts}  "
                  f"({(last_ts - first_ts).total_seconds():.0f} s total)\n\n")

    # Collapse consecutive identical labels into phase intervals.
    buf.write(f"{'phase':<34}{'first seen':<21}{'Δ to next (s)':>14}\n")
    buf.write("-" * 78 + "\n")
    collapsed = []
    for ts, lab, loc, msg in events:
        if collapsed and collapsed[-1][1] == lab:
            continue
        collapsed.append((ts, lab, loc, msg))
    for i, (ts, lab, loc, msg) in enumerate(collapsed):
        if i + 1 < len(collapsed):
            dt = (collapsed[i + 1][0] - ts).total_seconds()
        else:
            dt = (last_ts - ts).total_seconds() if last_ts else 0
        buf.write(f"{lab:<34}{str(ts):<21}{dt:>14.0f}\n")

    # Aggregate elapsed time per phase label (sum of intervals labeled the same).
    buf.write("\n" + "-" * 78 + "\n")
    buf.write("AGGREGATE wall-clock per phase label (sum of intervals)\n")
    buf.write("-" * 78 + "\n")
    agg = {}
    for i, (ts, lab, loc, msg) in enumerate(collapsed):
        nxt = collapsed[i + 1][0] if i + 1 < len(collapsed) else last_ts
        agg[lab] = agg.get(lab, 0.0) + (nxt - ts).total_seconds()
    for lab, sec in sorted(agg.items(), key=lambda kv: -kv[1]):
        buf.write(f"  {lab:<40}{sec:>8.0f} s\n")

    TIMELINE_OUT.write_text(buf.getvalue())
    print(f"Wrote {TIMELINE_OUT}")


def main():
    if META.exists():
        meta = json.loads(META.read_text())
        print(f"Run: {meta['command']}  wall={meta['wall_clock_s']} s  "
              f"rc={meta['return_code']}")
    write_pstats()
    parse_timeline()


if __name__ == "__main__":
    main()
