#!/bin/bash
# PR C, "PR" task 1 after Q16 = (c) "Bxa": full dev suite (./pypeit_test all,
# which includes the vet tests) on the arcfit worktree with the Bxa change
# (get_xcorr_arc always uses the curve_fit-matching line fit), into
# REDUX_OUT_PR/arcfitBxa.  The earlier B run stays in REDUX_OUT_PR/arcfit for
# reference.  Then the controls on three setups: develop2 (develop again, run-
# to-run noise) and arcfitA (Bxa code with fast_line_fit forced off via the
# variant hook).  Safe recipe: one test at a time (-t 1) under MemoryMax=80G
# (C.md, Running).
set -u
DEV=/home/xavier/Projects/PypeIt/PypeIt-development-suite
ENV=/home/xavier/miniconda3/envs/pypeit14b
WT=/mnt/tank/Astronomy/PypeIt/worktrees
HOOK=$DEV/pypeitdev/speed_up/scripts/arcfit/variant_hook
OUT=$DEV/REDUX_OUT_PR
export PATH=$ENV/bin:/usr/bin:/bin
export PYPEIT_DEV=$DEV/
export QT_QPA_PLATFORM=offscreen
unset PYPEIT_ARCFIT_VARIANT
cd $DEV

export PYTHONPATH=$WT/arcfit
mkdir -p $OUT/arcfitBxa
echo "$(date -Is) arcfitBxa start (-t 1); $(cd /tmp && python -c 'import pypeit; print(pypeit.__file__)' 2>/dev/null | tail -1)" >> $OUT/progress.txt
python ./pypeit_test all -o $OUT/arcfitBxa -t 1 -r $OUT/arcfitBxa/report.txt -c $OUT/arcfitBxa/perf.csv \
    > $OUT/arcfitBxa/stdout.txt 2>&1
echo "$(date -Is) arcfitBxa rc=$? $(grep -E 'SUITE' $OUT/arcfitBxa/report.txt | tail -1)" >> $OUT/progress.txt

CTRL="gemini_gmos/GN_E2V_MULTI_R400_600 gemini_gmos/GS_HAM_MULTI_R400_700 keck_mosfire/J_multi"
for v in develop2 arcfitA; do
    unset PYPEIT_ARCFIT_VARIANT
    case $v in
        develop2) export PYTHONPATH=$WT/develop ;;
        arcfitA)  export PYTHONPATH=$WT/arcfit:$HOOK; export PYPEIT_ARCFIT_VARIANT=A ;;
    esac
    rm -rf $OUT/$v; mkdir -p $OUT/$v
    echo "$(date -Is) $v start; $CTRL; $(cd /tmp && python -c 'import pypeit, pypeit.core.arc as a; print(pypeit.__file__, getattr(a, "ARCFIT_MATCH_CURVEFIT", "n/a"))' 2>/dev/null | tail -1)" >> $OUT/progress.txt
    python ./pypeit_test reduce -s $CTRL -o $OUT/$v -t 1 -r $OUT/$v/report.txt > $OUT/$v/stdout.txt 2>&1
    echo "$(date -Is) $v rc=$?" >> $OUT/progress.txt
done
echo "$(date -Is) all done" >> $OUT/progress.txt
