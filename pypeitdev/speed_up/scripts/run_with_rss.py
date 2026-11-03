#!/usr/bin/env python
"""
Run a command, time it, and sample the peak memory of its whole process tree.

Written for the PR B (detector parallelism) validation, where the reduction
forks one worker process per detector and the interesting question is how peak
memory scales with ``--ncpu``.

Why not ``/usr/bin/time -v``
----------------------------
``/usr/bin/time -v`` reports ``ru_maxrss``, which for a process with children is
the maximum RSS of *any single* process, not the total across the tree.  For a
4-worker reduction that under-reports the real footprint by roughly a factor of
four, which is precisely the number being measured.

Why PSS rather than RSS
-----------------------
The workers are forked, so they share a large amount of copy-on-write memory
with the parent.  Summing RSS over the tree double-counts every shared page and
over-reports.  Summing **PSS** (proportional set size, which divides each shared
page by the number of processes mapping it) counts each physical page exactly
once, so the tree total is the true footprint.  Both are reported, since PSS is
Linux-only and RSS is the number people are used to seeing.

Usage::

    python run_with_rss.py --json out.json --interval 2.0 -- run_pypeit file.pypeit -o --ncpu 4
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import psutil


def sample_tree(proc):
    """
    Sum the memory of a process and all of its descendants.

    Processes that exit while being walked are skipped, which is expected: the
    pool's workers come and go.

    Args:
        proc (:class:`psutil.Process`): Root of the tree.

    Returns:
        :obj:`tuple`: Total RSS and total PSS in bytes, and the number of
        processes that contributed.  PSS is None if it is unavailable.
    """
    rss = 0
    pss = 0
    pss_ok = True
    n = 0
    try:
        procs = [proc] + proc.children(recursive=True)
    except psutil.Error:
        return 0, None, 0
    for p in procs:
        try:
            rss += p.memory_info().rss
            if pss_ok:
                try:
                    pss += p.memory_full_info().pss
                except (psutil.AccessDenied, AttributeError):
                    pss_ok = False
            n += 1
        except psutil.Error:
            continue
    return rss, (pss if pss_ok else None), n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--json', type=str, default=None,
                    help='Write the measurements to this JSON file.')
    ap.add_argument('--interval', type=float, default=2.0,
                    help='Sampling interval in seconds.')
    ap.add_argument('--log', type=str, default=None,
                    help='Write the command stdout+stderr here instead of inheriting.')
    ap.add_argument('cmd', nargs=argparse.REMAINDER,
                    help='The command to run, after a bare --.')
    cli = ap.parse_args()

    cmd = cli.cmd[1:] if cli.cmd and cli.cmd[0] == '--' else cli.cmd
    if len(cmd) == 0:
        ap.error('no command given')

    out = open(cli.log, 'w') if cli.log else None
    t0 = time.perf_counter()
    child = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT if out else None)
    proc = psutil.Process(child.pid)

    peak_rss = 0
    peak_pss = 0
    peak_nproc = 0
    nsample = 0
    while child.poll() is None:
        rss, pss, n = sample_tree(proc)
        peak_rss = max(peak_rss, rss)
        if pss is not None:
            peak_pss = max(peak_pss, pss)
        peak_nproc = max(peak_nproc, n)
        nsample += 1
        time.sleep(cli.interval)
    wall = time.perf_counter() - t0
    rc = child.returncode
    if out:
        out.close()

    meta = {'command': ' '.join(cmd),
            'wall_clock_s': round(wall, 3),
            'return_code': rc,
            'peak_rss_gb': round(peak_rss / 1024**3, 3),
            'peak_pss_gb': round(peak_pss / 1024**3, 3) if peak_pss else None,
            'peak_process_count': peak_nproc,
            'n_samples': nsample,
            'sample_interval_s': cli.interval}

    print(json.dumps(meta, indent=2))
    if cli.json:
        Path(cli.json).parent.mkdir(parents=True, exist_ok=True)
        with open(cli.json, 'w') as fp:
            json.dump(meta, fp, indent=2)
    return rc


if __name__ == '__main__':
    sys.exit(main())
