#!/bin/bash
# PR C (Task 5) dev-suite vetting: reduce the four instruments once per
# fit_arcspec variant, each into its own output directory, then run the
# wavelength/tilt/extraction/skysub vet tests against it.  The variant is
# selected with variant_hook/sitecustomize.py, so pypeit/core/arc.py is not
# edited.  Run as a systemd user service (see C.md, Task 5 log).
set -u
DEV=/home/xavier/Projects/PypeIt/PypeIt-development-suite
ENV=/home/xavier/miniconda3/envs/pypeit14b
OUT=$DEV/REDUX_OUT_C
THREADS=${THREADS:-4}
VARIANTS=${VARIANTS:-"legacy gn10 gn0"}
export PATH=$ENV/bin:/usr/bin:/bin
export PYPEIT_DEV=$DEV/
export PYTHONPATH=$DEV/pypeitdev/speed_up/scripts/arcfit/variant_hook
export QT_QPA_PLATFORM=offscreen
mkdir -p $OUT
cd $DEV
echo "$(date -Is) start; PypeIt HEAD $(git -C /mnt/tank/Astronomy/PypeIt/PypeIt rev-parse --short HEAD); dirty: $(git -C /mnt/tank/Astronomy/PypeIt/PypeIt status --porcelain -- pypeit | wc -l) files" >> $OUT/progress.txt
for v in $VARIANTS; do
    export PYPEIT_ARCFIT_VARIANT=$v
    mkdir -p $OUT/$v
    echo "$(date -Is) $v: hook -> $(python -c 'import pypeit.core.arc as a; print(a.USE_VECTORIZED_ARCFIT, a.ARCFIT_GN_ITER, a.ARCFIT_RESCUE_ITER)' 2>/dev/null | tail -1)" >> $OUT/progress.txt
    echo "$(date -Is) $v: reduce start" >> $OUT/progress.txt
    python ./pypeit_test reduce -i shane_kast_blue shane_kast_red keck_deimos keck_lris_blue \
        -o $OUT/$v -t $THREADS -r $OUT/$v/reduce_report.txt -c $OUT/$v/reduce_perf.csv \
        > $OUT/$v/reduce_stdout.txt 2>&1
    echo "$(date -Is) $v: reduce rc=$?" >> $OUT/progress.txt
    python -m pytest -p no:cacheprovider -q -rA vet_tests/test_wavelengths.py vet_tests/test_wavetilts.py \
        vet_tests/test_extraction.py vet_tests/test_skysub.py --redux_out $OUT/$v \
        > $OUT/$v/vet_stdout.txt 2>&1
    echo "$(date -Is) $v: vet rc=$? $(tail -1 $OUT/$v/vet_stdout.txt)" >> $OUT/progress.txt
done
echo "$(date -Is) done" >> $OUT/progress.txt
