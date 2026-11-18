#!/bin/bash
# PR C, "PR" task 1: full dev suite (./pypeit_test all) on
#   develop : origin/develop (git worktree, detached)
#   arcfit  : origin/develop + speed_up_arcfit squash-merged, uncommitted
#             (git worktree), i.e. exactly what PR C would merge.
# The worktree is selected with PYTHONPATH (the pypeit14b editable install
# points at the main checkout).  The two suites run one after the other.
# Run as a systemd user service with a memory cap.
set -u
DEV=/home/xavier/Projects/PypeIt/PypeIt-development-suite
ENV=/home/xavier/miniconda3/envs/pypeit14b
WT=/mnt/tank/Astronomy/PypeIt/worktrees
OUT=$DEV/REDUX_OUT_PR
THREADS=${THREADS:-4}
VARIANTS=${VARIANTS:-"develop arcfit"}
export PATH=$ENV/bin:/usr/bin:/bin
export PYPEIT_DEV=$DEV/
export QT_QPA_PLATFORM=offscreen
unset PYPEIT_ARCFIT_VARIANT
mkdir -p $OUT
cd $DEV
for v in $VARIANTS; do
    export PYTHONPATH=$WT/$v
    mkdir -p $OUT/$v
    echo "$(date -Is) $v start; worktree $(git -C $WT/$v rev-parse --short HEAD) + $(git -C $WT/$v status --porcelain -- pypeit | wc -l) changed pypeit files; imports $(cd /tmp && python -c 'import pypeit; print(pypeit.__file__)' 2>/dev/null | tail -1)" >> $OUT/progress.txt
    python ./pypeit_test all -o $OUT/$v -t $THREADS -r $OUT/$v/report.txt -c $OUT/$v/perf.csv \
        > $OUT/$v/stdout.txt 2>&1
    echo "$(date -Is) $v rc=$? $(grep -E 'PASSED|FAILED' $OUT/$v/report.txt | grep -i 'suite' | tail -1)" >> $OUT/progress.txt
done
echo "$(date -Is) done" >> $OUT/progress.txt
