"""
Select a fit_arcspec variant for a whole dev-suite run, including the
run_pypeit subprocesses and their forkserver/spawn workers, without editing
pypeit/core/arc.py.  Put this directory on PYTHONPATH and set
PYPEIT_ARCFIT_VARIANT to one of:

    legacy  -> USE_VECTORIZED_ARCFIT = False   (the pre-PR-C code path)
    A       -> ARCFIT_MATCH_CURVEFIT = True     (curve_fit hand-offs, no width floor)
    B       -> ARCFIT_MATCH_CURVEFIT = False    (no hand-off + width floor; the default since Task 11)
    gnN     -> ARCFIT_GN_ITER = N               (e.g. gn0: analytic seed + rescue)

Note: the 2026-10-22 Task 5 run labelled ``gn10`` predates
ARCFIT_MATCH_CURVEFIT, and is equivalent to ``B``.

Unset or empty: no change.
"""
import os

_v = os.environ.get('PYPEIT_ARCFIT_VARIANT', '').strip()
if _v:
    import pypeit.core.arc as _arc
    if _v == 'legacy':
        _arc.USE_VECTORIZED_ARCFIT = False
    elif _v == 'A':
        _arc.USE_VECTORIZED_ARCFIT = True
        _arc.ARCFIT_MATCH_CURVEFIT = True
    elif _v == 'B':
        _arc.USE_VECTORIZED_ARCFIT = True
        _arc.ARCFIT_MATCH_CURVEFIT = False
    elif _v.startswith('gn') and _v[2:].isdigit():
        _arc.USE_VECTORIZED_ARCFIT = True
        _arc.ARCFIT_GN_ITER = int(_v[2:])
    else:
        raise ValueError(f'Unknown PYPEIT_ARCFIT_VARIANT={_v!r}')

    # Since Task 11 (2026-11-03), BuildWaveCalib and BuildWaveTilts set
    # ARCFIT_MATCH_CURVEFIT from [calibrations][wavelengths] fast_line_fit via
    # arc.arcfit_options, which would override the variant there.  Make that
    # context manager a no-op so the variant applies to every caller.
    if _v in ('A', 'B'):
        from contextlib import contextmanager as _cm

        @_cm
        def _fixed_options(match_curvefit=None):
            yield
        _arc.arcfit_options = _fixed_options
