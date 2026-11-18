#!/bin/bash
# PR C, "PR" task 1, controls for the setups whose products differ between
# develop and develop + PR C:
#   develop2 : develop again (run-to-run reproducibility)
#   arcfitA  : develop + PR C with option A forced for every caller (variant
#              hook: ARCFIT_MATCH_CURVEFIT = True, hence no width floor,
#              arcfit_options neutralised) -> should reproduce develop
# Waits for the main comparison (arcfit_pr_compare.sh) to finish first.
set -u
DEV=/home/xavier/Projects/PypeIt/PypeIt-development-suite
ENV=/home/xavier/miniconda3/envs/pypeit14b
WT=/mnt/tank/Astronomy/PypeIt/worktrees
HOOK=$DEV/pypeitdev/speed_up/scripts/arcfit/variant_hook
OUT=$DEV/REDUX_OUT_PR
SETUPS=${SETUPS:-"gemini_gmos/GN_E2V_MULTI_R400_600 gemini_gmos/GS_HAM_MULTI_R400_700 keck_mosfire/J_multi"}
export PATH=$ENV/bin:/usr/bin:/bin
export PYPEIT_DEV=$DEV/
export QT_QPA_PLATFORM=offscreen
until grep -q " done$" $OUT/progress.txt 2>/dev/null; do sleep 120; done
cd $DEV
for v in develop2 arcfitA; do
    unset PYPEIT_ARCFIT_VARIANT
    case $v in
        develop2) export PYTHONPATH=$WT/develop ;;
        arcfitA)  export PYTHONPATH=$WT/arcfit:$HOOK; export PYPEIT_ARCFIT_VARIANT=A ;;
    esac
    mkdir -p $OUT/$v
    echo "$(date -Is) $v start; setups $SETUPS; $(cd /tmp && python -c 'import pypeit, pypeit.core.arc as a; print(pypeit.__file__, getattr(a, "ARCFIT_MATCH_CURVEFIT", "n/a"))' 2>/dev/null | tail -1)" >> $OUT/progress.txt
    python ./pypeit_test reduce -s $SETUPS -o $OUT/$v -t 3 -r $OUT/$v/report.txt > $OUT/$v/stdout.txt 2>&1
    echo "$(date -Is) $v rc=$?" >> $OUT/progress.txt
done
echo "$(date -Is) controls done" >> $OUT/progress.txt
