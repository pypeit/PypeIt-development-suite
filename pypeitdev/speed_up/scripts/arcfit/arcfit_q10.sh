#!/bin/bash
# PR C, Q10: re-profile DEIMOS 600ZD (cProfile, --ncpu 1) after caching the
# archive-sky FWHM in pypeit.core.flexure.get_archive_spectrum.
#  1. variant A (shipped default), tag C_A_cache, then analyze
#  2. variant B (fast_line_fit for all callers), tag C_B_cache, then analyze
# Compare with the C_A / C_B profiles from arcfit_task6.sh.
# The variant is selected with variant_hook/sitecustomize.py.
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
log "start Q10 re-profile; flexure cache present: $(grep -c '_archive_sky_fwhm' /mnt/tank/Astronomy/PypeIt/PypeIt/pypeit/core/flexure.py) refs"

mkdir -p $OUT/_profiles
for v in A B; do
    log "START profile_deimos.py $v tag C_${v}_cache"
    PYPEIT_ARCFIT_VARIANT=$v $PY $SCR/profile_deimos.py --tag C_${v}_cache \
        > $OUT/_profiles/profile_deimos.py.${v}_cache.stdout.txt 2>&1
    log "END profile_deimos.py $v rc=$?"
    $PY $SCR/analyze_profile.py --stem keck_deimos_600zd_m_6500_C_${v}_cache >> $OUT/_profiles/analyze_q10.txt 2>&1
    log "analyze $v rc=$?"
done
log "done"
