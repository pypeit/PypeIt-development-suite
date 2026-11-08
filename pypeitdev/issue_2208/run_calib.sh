#!/bin/bash
# Cold, calibrations-only (-c) run of keck_lris_blue/multi_600_4000_slitmask
# into $REDUX/<name>, with BLAS/OpenMP pinned to one thread.
#
# Usage: run_calib.sh <name> [VAR=value ...]
#   e.g. run_calib.sh hsw_a OPENBLAS_CORETYPE=Haswell
#
# Written for PypeIt Issue #2208.
set -euo pipefail

NAME=$1; shift
for kv in "$@"; do export "$kv"; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

ENVBIN=/home/xavier/miniconda3/envs/pypeit14b/bin
DEV=/mnt/tank/Astronomy/PypeIt/PypeIt-development-suite
REDUX=/mnt/tank/Astronomy/PypeIt/REDUX/issue_2208
SETUP=${SETUP:-keck_lris_blue_multi_600_4000_slitmask}
INST_SETUP=${INST_SETUP:-keck_lris_blue/multi_600_4000_slitmask}

OUT=$REDUX/$NAME
rm -rf "$OUT"
mkdir -p "$OUT"

# Point the data block at the raw data
sed "s#^ *path .*# path $DEV/RAW_DATA/$INST_SETUP/#" \
    "$DEV/pypeit_files/$SETUP.pypeit" > "$OUT/$SETUP.pypeit"

# Record the BLAS runtimes as this process tree will see them
"$ENVBIN/python" -c "
import os, numpy, scipy.linalg, sklearn.decomposition, threadpoolctl
print({k: v for k, v in os.environ.items() if k.endswith('_THREADS') or k.startswith('OPENBLAS') or k.startswith('MKL')})
for i in threadpoolctl.threadpool_info():
    print(i['filepath'].split('/')[-1], i['internal_api'], i.get('architecture'), i['version'], i['num_threads'])
" > "$OUT/blas_info.txt"

start=$(date +%s.%N)
"$ENVBIN/run_pypeit" "$OUT/$SETUP.pypeit" -c -o -r "$OUT" > "$OUT/run.log" 2>&1
status=$?
end=$(date +%s.%N)
echo "$NAME status=$status wall=$(echo "$end - $start" | bc) s" | tee -a "$OUT/timing.txt"
