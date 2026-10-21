#!/bin/bash
# DevSuite reduce + afterburn on an edge-tracing subset, then the vet tests
# whose inputs that subset fully covers.  BLAS/OpenMP pinned to one thread
# per reduction.
#
# Usage: run_devsuite_2208.sh <tag> [threads]
# Outputs: $REDUX/devsuite_<tag>/ (REDUX_OUT-style tree), plus
#          pypeit_test_report.txt, vet_report.txt and done markers there.
#
# Written for PypeIt Issue #2208.

TAG=$1
THREADS=${2:-4}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH I2208_PCA_FIX I2208_PCA_CAPTURE
export PYPEIT_DEV=/mnt/tank/Astronomy/PypeIt/PypeIt-development-suite
export PATH=/home/xavier/miniconda3/envs/pypeit14b/bin:$PATH
export QT_QPA_PLATFORM=offscreen

OUT=/mnt/tank/Astronomy/PypeIt/REDUX/issue_2208/devsuite_$TAG
mkdir -p "$OUT"

SETUPS=(
    # Multi-slit and slit-mask (edge-tracing PCA)
    keck_lris_blue/multi_600_4000_slitmask
    keck_lris_blue/long_600_4000_d560
    keck_lris_red/multi_400_8500_d560
    keck_lris_red_mark4/multi_600_10000_slitmask
    keck_deimos/600ZD_M_6500
    keck_deimos/830G_M_8500
    keck_deimos/830G_M_9000_dither
    keck_mosfire/J_multi
    gemini_gmos/GS_HAM_B600_MOS
    mmt_binospec/Multislit_G270
    # Wide long slit (tilt-tracing PCA)
    mmt_binospec/Longslit_G1000
    # Echelle (edge PCA and echelle object-finding PCA)
    keck_esi/Ech_1x1
    magellan_mage/1x1
    vlt_xshooter/VIS_1x1
    keck_nires/ABBA_wstandard
    # Long-slit control
    shane_kast_blue/600_4310_d55
)

# Vet tests whose inputs are all produced by the setups above
VET=(
    vet_tests/test_slitmask.py
    # The test_chk_*_slits tests read $PYPEIT_DEV/Cooked, which is not present
    vet_tests/test_edgetrace.py::test_addrm_slit
    vet_tests/test_edgetrace.py::test_sobel_enhance
    vet_tests/test_findobj.py
    vet_tests/test_wavetilts.py
    vet_tests/test_skysub.py
    vet_tests/test_fluxspec.py
    vet_tests/test_collate1d.py
    vet_tests/test_echelles.py::test_magellan_mage_orders
)

cd "$PYPEIT_DEV"
echo "devsuite $TAG started $(date); threads=$THREADS; pypeit=$(python -c 'import pypeit; print(pypeit.__file__)')"
start=$(date +%s)
./pypeit_test reduce afterburn -t "$THREADS" -o "$OUT" -s "${SETUPS[@]}" \
    -r "$OUT/pypeit_test_report.txt" -c "$OUT/performance.csv"
echo "pypeit_test exit=$? after $(( ($(date +%s) - start) / 60 )) min" | tee "$OUT/pypeit_test.done"

python -m pytest "${VET[@]}" --redux_out "$OUT" -ra -p no:cacheprovider \
    > "$OUT/vet_report.txt" 2>&1
echo "vet exit=$? ; $(tail -1 "$OUT/vet_report.txt")" | tee "$OUT/vet.done"
echo "devsuite $TAG ALLDONE $(date); total $(( ($(date +%s) - start) / 60 )) min"
