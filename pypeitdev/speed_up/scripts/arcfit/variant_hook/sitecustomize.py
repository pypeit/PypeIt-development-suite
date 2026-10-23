"""
Select a fit_arcspec variant for a whole dev-suite run, including the
run_pypeit subprocesses and their forkserver/spawn workers, without editing
pypeit/core/arc.py.  Put this directory on PYTHONPATH and set
PYPEIT_ARCFIT_VARIANT to one of:

    legacy  -> USE_VECTORIZED_ARCFIT = False   (the pre-PR-C code path)
    gn0     -> ARCFIT_GN_ITER = 0               (analytic seed + rescue)
    gnN     -> ARCFIT_GN_ITER = N               (e.g. gn10, the shipped value)

Unset or empty: no change.
"""
import os

_v = os.environ.get('PYPEIT_ARCFIT_VARIANT', '').strip()
if _v:
    import pypeit.core.arc as _arc
    if _v == 'legacy':
        _arc.USE_VECTORIZED_ARCFIT = False
    elif _v.startswith('gn') and _v[2:].isdigit():
        _arc.USE_VECTORIZED_ARCFIT = True
        _arc.ARCFIT_GN_ITER = int(_v[2:])
    else:
        raise ValueError(f'Unknown PYPEIT_ARCFIT_VARIANT={_v!r}')
