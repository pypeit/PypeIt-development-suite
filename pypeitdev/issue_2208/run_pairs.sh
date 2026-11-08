#!/bin/bash
# Two cold, full, serial run_pypeit reductions of each setup, compared with
# compare_calibs.py (Calibrations and Science).  BLAS/OpenMP pinned to one
# thread; wall time and peak memory from /usr/bin/time -v.
#
# Usage: run_pairs.sh <tag> [instrument/setup ...]
#   default setups: keck_lris_blue/multi_600_4000_slitmask
#                   keck_deimos/830G_M_8500
#                   mmt_binospec/Longslit_G1000
# Outputs: $REDUX/pairs_<tag>/<instr>_<setup>_{a,b}/, and
#          $REDUX/pairs_<tag>/summary.txt
#
# Written for PypeIt Issue #2208.

TAG=$1; shift
SETUPS=("$@")
if [ ${#SETUPS[@]} -eq 0 ]; then
    SETUPS=(keck_lris_blue/multi_600_4000_slitmask keck_deimos/830G_M_8500
            mmt_binospec/Longslit_G1000)
fi

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH I2208_PCA_FIX I2208_PCA_CAPTURE

ENVBIN=/home/xavier/miniconda3/envs/pypeit14b/bin
DEV=/mnt/tank/Astronomy/PypeIt/PypeIt-development-suite
HERE=$DEV/pypeitdev/issue_2208
TOP=/mnt/tank/Astronomy/PypeIt/REDUX/issue_2208/pairs_$TAG
mkdir -p "$TOP"
SUMMARY=$TOP/summary.txt
echo "run_pairs $TAG started $(date)  pypeit=$("$ENVBIN/python" -c 'import pypeit; print(pypeit.__file__)')" >> "$SUMMARY"

for is in "${SETUPS[@]}"; do
    instr=${is%%/*}
    setup=${is#*/}
    lower=$(echo "${instr}_${setup}" | tr '[:upper:]' '[:lower:]')
    for run in a b; do
        OUT=$TOP/${instr}_${setup}_$run
        rm -rf "$OUT"
        mkdir -p "$OUT"
        # Point the data block at the raw data (as fix_pypeit_file_directory does)
        sed "s#^ *path .*# path $DEV/RAW_DATA/$instr/$setup/#" \
            "$DEV/pypeit_files/$lower.pypeit" > "$OUT/$lower.pypeit"
        /usr/bin/time -v -o "$OUT/time.txt" \
            "$ENVBIN/run_pypeit" "$OUT/$lower.pypeit" -o -r "$OUT" > "$OUT/run.log" 2>&1
        status=$?
        wall=$(grep "Elapsed (wall clock)" "$OUT/time.txt" | awk '{print $NF}')
        rss=$(grep "Maximum resident set size" "$OUT/time.txt" | awk '{printf "%.2f GB", $NF/1048576}')
        echo "$instr/$setup $run status=$status wall=$wall maxrss=$rss" | tee -a "$SUMMARY"
    done
    "$ENVBIN/python" "$HERE/compare_calibs.py" "$TOP/${instr}_${setup}_a" "$TOP/${instr}_${setup}_b" \
        > "$TOP/cmp_${instr}_${setup}.txt" 2>&1
    ndiff=$(grep -cE "^\S.*DIFFER|only in" "$TOP/cmp_${instr}_${setup}.txt")
    nfiles=$(grep -cE "arrays" "$TOP/cmp_${instr}_${setup}.txt")
    first=$(grep -m1 -E "^\S.*DIFFER" "$TOP/cmp_${instr}_${setup}.txt")
    echo "$instr/$setup compare: $ndiff of $nfiles files differ; first: ${first:-none}" | tee -a "$SUMMARY"
done
echo "run_pairs $TAG ALLDONE $(date)" | tee -a "$SUMMARY"
