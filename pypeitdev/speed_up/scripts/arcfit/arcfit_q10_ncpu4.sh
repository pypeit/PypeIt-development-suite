#!/bin/bash
# PR C, Q10: uninstrumented DEIMOS 600ZD --ncpu 4 run of the shipped variant A
# with the archive-sky FWHM cache (B's run_with_rss + instr_run_pypeit recipe),
# to compare with REDUX_OUT_C/_ncpu4/{A,legacy} from arcfit_task6.sh.
set -u
DEV=/home/xavier/Projects/PypeIt/PypeIt-development-suite
ENV=/home/xavier/miniconda3/envs/pypeit14b
PY=$ENV/bin/python
SCR=$DEV/pypeitdev/speed_up/scripts
OUT=$DEV/REDUX_OUT_C
P=$OUT/q10_progress.txt
export PATH=$ENV/bin:/usr/bin:/bin
export PYPEIT_DEV=$DEV/
export PYTHONPATH=$SCR/arcfit/variant_hook
export QT_QPA_PLATFORM=offscreen
log() { echo "$(date -Is) $*" >> $P; }
DEIM=$DEV/REDUX_OUT/keck_deimos/600ZD_M_6500/keck_deimos_600zd_m_6500.pypeit
INSTR=$DEV/REDUX_OUT/_port_validation/instr_run_pypeit.py
d=$OUT/_ncpu4/A_cache; rm -rf $d; mkdir -p $d
log "START deimos ncpu4 A_cache"
PYPEIT_ARCFIT_VARIANT=A $PY $SCR/run_with_rss.py --json $d/rss.json --interval 2.0 -- \
    $PY $INSTR $DEIM -o -r $d --ncpu 4 --log_file $d/run.log > $d/stdout.txt 2> $d/stderr.txt
log "END deimos ncpu4 A_cache rc=$?"
log "ncpu4 done"
