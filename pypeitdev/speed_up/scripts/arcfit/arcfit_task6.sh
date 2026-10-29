#!/bin/bash
# PR C Tasks 5 (shipped option A) and 6 (re-profiling): one job at a time.
#  1. dev-suite reduce + vet for variant A (into REDUX_OUT_C/A)
#  2. cProfile, --ncpu 1: shane_kast_blue under legacy, A and B
#  3. cProfile, --ncpu 1: keck_deimos 600ZD_M_6500 under A, then legacy
#  4. uninstrumented --ncpu 4 DEIMOS (run_with_rss + instr_run_pypeit, as in
#     B's port validation) under A, then legacy
#  5. cProfile, --ncpu 1: DEIMOS under B
# The variant is selected with variant_hook/sitecustomize.py.
set -u
DEV=/home/xavier/Projects/PypeIt/PypeIt-development-suite
ENV=/home/xavier/miniconda3/envs/pypeit14b
PY=$ENV/bin/python
SCR=$DEV/pypeitdev/speed_up/scripts
OUT=$DEV/REDUX_OUT_C
P=$OUT/task6_progress.txt
export PATH=$ENV/bin:/usr/bin:/bin
export PYPEIT_DEV=$DEV/
export PYTHONPATH=$SCR/arcfit/variant_hook
export QT_QPA_PLATFORM=offscreen
log() { echo "$(date -Is) $*" >> $P; }
log "start; PypeIt HEAD $(git -C /mnt/tank/Astronomy/PypeIt/PypeIt rev-parse --short HEAD); dirty pypeit/ files: $(git -C /mnt/tank/Astronomy/PypeIt/PypeIt status --porcelain -- pypeit | wc -l)"

# 1. dev suite, variant A
VARIANTS=A THREADS=4 $SCR/arcfit/arcfit_vet.sh
log "devsuite A done: $(grep ' A: ' $OUT/progress.txt | tail -2 | tr '\n' ' ')"

prof() {  # script variant extra-args...
    local script=$1 v=$2; shift 2
    log "START $script $v $*"
    PYPEIT_ARCFIT_VARIANT=$v $PY $SCR/$script --tag C_$v "$@" > $OUT/_profiles/$script.$v.stdout.txt 2>&1
    log "END $script $v rc=$?"
}
mkdir -p $OUT/_profiles
# 2. Kast
for v in legacy A B; do
    prof profile_kast_blue.py $v
    PYPEIT_ARCFIT_VARIANT=$v $PY $SCR/analyze_profile.py --stem shane_kast_blue_A_C_$v >> $OUT/_profiles/analyze.txt 2>&1
done
# 3. DEIMOS serial, A then legacy
for v in A legacy; do
    prof profile_deimos.py $v
    $PY $SCR/analyze_profile.py --stem keck_deimos_600zd_m_6500_C_$v >> $OUT/_profiles/analyze.txt 2>&1
done
# 4. DEIMOS --ncpu 4, uninstrumented
DEIM=$DEV/REDUX_OUT/keck_deimos/600ZD_M_6500/keck_deimos_600zd_m_6500.pypeit
INSTR=$DEV/REDUX_OUT/_port_validation/instr_run_pypeit.py
for v in A legacy; do
    d=$OUT/_ncpu4/$v; rm -rf $d; mkdir -p $d
    log "START deimos ncpu4 $v"
    PYPEIT_ARCFIT_VARIANT=$v $PY $SCR/run_with_rss.py --json $d/rss.json --interval 2.0 -- \
        $PY $INSTR $DEIM -o -r $d --ncpu 4 --log_file $d/run.log > $d/stdout.txt 2> $d/stderr.txt
    log "END deimos ncpu4 $v rc=$?"
done
# 5. DEIMOS serial, B
prof profile_deimos.py B
$PY $SCR/analyze_profile.py --stem keck_deimos_600zd_m_6500_C_B >> $OUT/_profiles/analyze.txt 2>&1
log "done"
