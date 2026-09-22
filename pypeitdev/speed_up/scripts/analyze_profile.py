#!/usr/bin/env python
"""
Analyze the cProfile output + run log produced by ``profile_kast_blue.py`` or
``profile_deimos.py``.

Produces two machine-readable text artifacts in ``Reports/``:

- ``<stem>.pstats.txt`` — top functions by cumulative time and by internal
  (tottime) time.
- ``<stem>.timeline.txt`` — per-phase wall-clock timeline parsed from the
  timestamped run log, mapped onto the ``pypeit_workflow.md`` phases. For
  multi-detector setups (e.g. keck_deimos) it also reports a per-detector/mosaic
  wall-clock breakdown.

The dataset is selected with ``--stem`` (default ``shane_kast_blue_A``)::

    /home/xavier/miniconda3/envs/pypeit/bin/python analyze_profile.py \
        --stem keck_deimos_600zd_m_6500
"""
import io
import re
import sys
import json
import pstats
import argparse
import datetime
from pathlib import Path

REPORTS = Path("/home/xavier/Projects/PypeIt/PypeIt-development-suite/"
               "pypeitdev/speed_up/Reports")

# Set by main() from --stem.
PROF = LOG = META = PSTATS_OUT = TIMELINE_OUT = None


def set_paths(stem):
    global PROF, LOG, META, PSTATS_OUT, TIMELINE_OUT
    PROF = REPORTS / f"{stem}.prof"
    LOG = REPORTS / f"{stem}.run.log"
    META = REPORTS / f"{stem}.runmeta.json"
    PSTATS_OUT = REPORTS / f"{stem}.pstats.txt"
    TIMELINE_OUT = REPORTS / f"{stem}.timeline.txt"

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


# ---------------------------------------------------------------------------
# 3. Per-detector / per-mosaic wall-clock breakdown (multi-detector setups)
# ---------------------------------------------------------------------------

# keck_deimos (and other multi-detector setups) emit, inside the science
# reduction loop, one "Calibrating detector (a, b)" line per mosaic; calibrations
# are then built lazily for that detector before it is reduced. These lines are
# clean per-mosaic boundaries.
DET_MARKER = re.compile(
    r"reduce_exposure:\d+ \| Calibrating detector (?P<det>.+)$")


def per_detector_breakdown():
    if not LOG.exists():
        return
    rows = []  # (datetime, det_label)
    first_ts = last_ts = None
    for line in LOG.read_text(errors="replace").splitlines():
        m = LOGLINE.match(line)
        if not m:
            continue
        ts = datetime.datetime.strptime(m["ts"], "%Y-%m-%d %H:%M:%S")
        if first_ts is None:
            first_ts = ts
        last_ts = ts
        dm = DET_MARKER.search(line)
        if dm:
            rows.append((ts, dm["det"].strip()))

    if not rows:
        return  # single-detector setup; nothing to break down

    buf = io.StringIO()
    buf.write("\n" + "=" * 78 + "\n")
    buf.write("PER-DETECTOR / PER-MOSAIC WALL-CLOCK (science reduction loop)\n")
    buf.write("Interval = 'Calibrating detector X' -> next detector (or log end).\n")
    buf.write("Includes that detector's lazily-built calibrations + reduction.\n")
    buf.write("=" * 78 + "\n")
    buf.write(f"{'detector/mosaic':<22}{'start':<21}{'wall (s)':>10}\n")
    buf.write("-" * 78 + "\n")
    total = 0.0
    for i, (ts, det) in enumerate(rows):
        nxt = rows[i + 1][0] if i + 1 < len(rows) else last_ts
        dt = (nxt - ts).total_seconds()
        total += dt
        buf.write(f"{det:<22}{str(ts):<21}{dt:>10.0f}\n")
    buf.write("-" * 78 + "\n")
    buf.write(f"{'sum (reduction loop)':<22}{'':<21}{total:>10.0f}\n")
    if first_ts is not None:
        pre = (rows[0][0] - first_ts).total_seconds()
        buf.write(f"{'(pre-loop: init/setup)':<22}{str(first_ts):<21}{pre:>10.0f}\n")

    with open(TIMELINE_OUT, "a") as fp:
        fp.write(buf.getvalue())
    print(f"Appended per-detector breakdown to {TIMELINE_OUT}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stem", default="shane_kast_blue_A",
                    help="Report file stem (e.g. keck_deimos_600zd_m_6500)")
    args = ap.parse_args()
    set_paths(args.stem)

    if META.exists():
        meta = json.loads(META.read_text())
        print(f"Run: {meta['command']}  wall={meta['wall_clock_s']} s  "
              f"rc={meta['return_code']}")
    write_pstats()
    parse_timeline()
    per_detector_breakdown()


if __name__ == "__main__":
    main()
