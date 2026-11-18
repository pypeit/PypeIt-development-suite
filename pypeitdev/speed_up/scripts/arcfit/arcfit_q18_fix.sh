#!/bin/bash
# "PR" task 3 (Q18 = (a)): re-run the setups touched by the three fixes
# (one-copy MMIRS templates, GNIRS 32_SB_SXD order-8 tilt spec_order = 2)
# on the arcfit worktree (develop 90c41e1a9 + PR C Bxa + fixes) into
# REDUX_OUT_PR/arcfitBxaFix.  Then build REDUX_OUT_PR/bxafix_merged, which
# links every setup to arcfitBxa except the re-run ones, and run the vet tests
# and the dev-suite unit test test_maskdesign_matching against it.
set -u
DEV=/home/xavier/Projects/PypeIt/PypeIt-development-suite
ENV=/home/xavier/miniconda3/envs/pypeit14b
WT=/mnt/tank/Astronomy/PypeIt/worktrees
OUT=$DEV/REDUX_OUT_PR
export PATH=$ENV/bin:/usr/bin:/bin
export PYPEIT_DEV=$DEV/
export QT_QPA_PLATFORM=offscreen
export PYTHONPATH=$WT/arcfit
unset PYPEIT_ARCFIT_VARIANT
cd $DEV
SETUPS="mmt_mmirs/H3000_H mmt_mmirs/HK_HK3 mmt_mmirs/HK_zJ mmt_mmirs/J_zJ mmt_mmirs/K_K mmt_mmirs/nep_as1_MOS gemini_gnirs_echelle/32_SB_SXD"
F=$OUT/arcfitBxaFix
rm -rf $F; mkdir -p $F
echo "$(date -Is) arcfitBxaFix start; $SETUPS" >> $OUT/progress.txt
python ./pypeit_test reduce afterburn -s $SETUPS -o $F -t 1 -r $F/report.txt > $F/stdout.txt 2>&1
echo "$(date -Is) arcfitBxaFix rc=$? $(grep -E 'SUITE' $F/report.txt | tail -1)" >> $OUT/progress.txt

M=$OUT/bxafix_merged
rm -rf $M; mkdir -p $M
for inst in $OUT/arcfitBxa/*/; do
    i=$(basename $inst); mkdir -p $M/$i
    for s in $inst*; do ln -s $s $M/$i/$(basename $s); done
done
for s in $SETUPS; do rm -f $M/$s; ln -s $F/$s $M/$s; done
python -m pytest vet_tests --redux_out $M -p no:cacheprovider -q > $OUT/arcfitBxaFix/vet.txt 2>&1
echo "$(date -Is) vet: $(tail -1 $F/vet.txt)" >> $OUT/progress.txt
python -m pytest unit_tests/test_maskdesign_matching.py -p no:cacheprovider -q > $F/unit_maskdesign.txt 2>&1
echo "$(date -Is) unit maskdesign: $(tail -1 $F/unit_maskdesign.txt)" >> $OUT/progress.txt
echo "$(date -Is) q18 all done" >> $OUT/progress.txt
