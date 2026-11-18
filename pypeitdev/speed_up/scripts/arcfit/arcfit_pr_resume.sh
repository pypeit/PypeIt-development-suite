#!/bin/bash
# PR C, "PR" task 1: resume the develop + PR C full suite after systemd-oomd
# killed it on 2026-11-04 21:38 (32 tests in 23 setups unfinished, vet tests
# not run), then run the controls (develop2, arcfitA).  Safe recipe: one test
# at a time (-t 1) under MemoryMax=80G (C.md, Running).
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
SETUPS=$(cat $OUT/arcfit_resume_setups.txt | tr '\n' ' ')
echo "$(date -Is) arcfit resume start (-t 1): $SETUPS" >> $OUT/progress.txt
python ./pypeit_test reduce afterburn ql -s $SETUPS -o $OUT/arcfit -t 1 -r $OUT/arcfit/report_resume.txt \
    -c $OUT/arcfit/perf_resume.csv > $OUT/arcfit/stdout_resume.txt 2>&1
echo "$(date -Is) arcfit resume rc=$? $(grep -E 'SUITE' $OUT/arcfit/report_resume.txt | tail -1)" >> $OUT/progress.txt
python ./pypeit_test vet -o $OUT/arcfit -t 1 -r $OUT/arcfit/report_vet.txt > $OUT/arcfit/stdout_vet.txt 2>&1
echo "$(date -Is) arcfit vet rc=$? $(grep -E 'VET' $OUT/arcfit/report_vet.txt | tail -1)" >> $OUT/progress.txt
CTRL="gemini_gmos/GN_E2V_MULTI_R400_600 gemini_gmos/GS_HAM_MULTI_R400_700 keck_mosfire/J_multi"
for v in develop2 arcfitA; do
    unset PYPEIT_ARCFIT_VARIANT
    case $v in
        develop2) export PYTHONPATH=$WT/develop ;;
        arcfitA)  export PYTHONPATH=$WT/arcfit:$HOOK; export PYPEIT_ARCFIT_VARIANT=A ;;
    esac
    mkdir -p $OUT/$v
    echo "$(date -Is) $v start; $CTRL; $(cd /tmp && python -c 'import pypeit, pypeit.core.arc as a; print(pypeit.__file__, getattr(a, "ARCFIT_MATCH_CURVEFIT", "n/a"))' 2>/dev/null | tail -1)" >> $OUT/progress.txt
    python ./pypeit_test reduce -s $CTRL -o $OUT/$v -t 1 -r $OUT/$v/report.txt > $OUT/$v/stdout.txt 2>&1
    echo "$(date -Is) $v rc=$?" >> $OUT/progress.txt
done
echo "$(date -Is) all done" >> $OUT/progress.txt
