#!/bin/bash
# PR C (Task 11) dev-suite vetting of the new default: option B (no curve_fit
# hand-off) plus the median-width floor in detect_lines.  No variant hook: the
# code defaults are what is tested.  Reduces the four instruments into
# REDUX_OUT_C/Bfloor, then runs the wavelength/tilt/extraction/skysub vet tests
# and the new vet_tests/test_arcfit.py.  Run as a systemd user service.
set -u
DEV=/home/xavier/Projects/PypeIt/PypeIt-development-suite
ENV=/home/xavier/miniconda3/envs/pypeit14b
OUT=$DEV/REDUX_OUT_C
V=${V:-Bfloor}
THREADS=${THREADS:-4}
export PATH=$ENV/bin:/usr/bin:/bin
export PYPEIT_DEV=$DEV/
unset PYTHONPATH PYPEIT_ARCFIT_VARIANT
export QT_QPA_PLATFORM=offscreen
mkdir -p $OUT/$V
cd $DEV
echo "$(date -Is) $V start; PypeIt HEAD $(git -C /mnt/tank/Astronomy/PypeIt/PypeIt rev-parse --short HEAD); dirty: $(git -C /mnt/tank/Astronomy/PypeIt/PypeIt status --porcelain -- pypeit | wc -l) files; defaults -> $(python -c 'import pypeit.core.arc as a; print(a.USE_VECTORIZED_ARCFIT, a.ARCFIT_MATCH_CURVEFIT, a.ARCFIT_MIN_FRAC_MED_WIDTH)' 2>/dev/null | tail -1)" >> $OUT/floor_progress.txt
python ./pypeit_test reduce -i shane_kast_blue shane_kast_red keck_deimos keck_lris_blue \
    -o $OUT/$V -t $THREADS -r $OUT/$V/reduce_report.txt -c $OUT/$V/reduce_perf.csv \
    > $OUT/$V/reduce_stdout.txt 2>&1
echo "$(date -Is) $V: reduce rc=$?" >> $OUT/floor_progress.txt
python -m pytest -p no:cacheprovider -q -rA vet_tests/test_wavelengths.py vet_tests/test_wavetilts.py \
    vet_tests/test_extraction.py vet_tests/test_skysub.py vet_tests/test_arcfit.py --redux_out $OUT/$V \
    > $OUT/$V/vet_stdout.txt 2>&1
echo "$(date -Is) $V: vet rc=$? $(tail -1 $OUT/$V/vet_stdout.txt)" >> $OUT/floor_progress.txt
echo "$(date -Is) $V done" >> $OUT/floor_progress.txt
