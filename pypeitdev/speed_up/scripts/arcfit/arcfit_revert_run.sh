#!/bin/bash
# "Revert to CurveFit?" item 2 (Q20 = R): re-run the 13 Keck/DEIMOS setups
# whose arc line lists differ from develop under R (curve_fit + median-width
# floor + archive-sky FWHM cache) on the revert worktree (develop 90c41e1a9 +
# PR C reverted to curve_fit) into REDUX_OUT_PR/revertR.  Then the DEIMOS vet
# tests and the arcfit vet tests on that output only.  Safe recipe: one test
# at a time (-t 1); run under systemd-run with MemoryMax (C.md, Running).
set -u
DEV=/home/xavier/Projects/PypeIt/PypeIt-development-suite
ENV=/home/xavier/miniconda3/envs/pypeit14b
WT=/mnt/tank/Astronomy/PypeIt/worktrees
OUT=$DEV/REDUX_OUT_PR
export PATH=$ENV/bin:/usr/bin:/bin
export PYPEIT_DEV=$DEV/
export QT_QPA_PLATFORM=offscreen
export PYTHONPATH=$WT/revert
unset PYPEIT_ARCFIT_VARIANT
cd $DEV
SETUPS="keck_deimos/600ZD_M_6500 keck_deimos/1200G_M_7750 keck_deimos/830G_LVM_8400 keck_deimos/1200B_M_5200 keck_deimos/900ZD_M_6000 keck_deimos/830G_M_8500 keck_deimos/1200G_M_5500 keck_deimos/600ZD_tilted keck_deimos/830G_M_9000_dither keck_deimos/1200B_LVM_5200 keck_deimos/900ZD_LVM_5500 keck_deimos/830G_L_8100 keck_deimos/830G_M_8100_26"
F=$OUT/revertR
rm -rf $F; mkdir -p $F
echo "$(date -Is) revertR start; $(cd /tmp && python -c 'import pypeit, pypeit.core.arc as a; print(pypeit.__file__, hasattr(a, "USE_VECTORIZED_ARCFIT"))' 2>/dev/null | tail -1); $SETUPS" >> $OUT/progress.txt
python ./pypeit_test reduce afterburn -s $SETUPS -o $F -t 1 -r $F/report.txt -c $F/perf.csv > $F/stdout.txt 2>&1
echo "$(date -Is) revertR rc=$? $(grep -E 'SUITE' $F/report.txt | tail -1)" >> $OUT/progress.txt
python -m pytest vet_tests --redux_out $F -k "deimos or arcfit" -p no:cacheprovider -q > $F/vet.txt 2>&1
echo "$(date -Is) revertR vet: $(tail -1 $F/vet.txt)" >> $OUT/progress.txt
echo "$(date -Is) revertR all done" >> $OUT/progress.txt
