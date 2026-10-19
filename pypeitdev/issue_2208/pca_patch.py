"""
Runtime patch of :mod:`pypeit.core.pca` for diagnosing PypeIt Issue #2208.

PypeIt itself is not modified.  :func:`install` replaces the module-level
``PCA`` name used by :func:`pypeit.core.pca.pca_decomposition` with a factory
that can force the SVD solver or the random seed, and wraps
``pca_decomposition`` so that every call's inputs, the solver sklearn actually
used, and the caller are recorded (and optionally saved).

Modes (``fix``):
    - ``None``: unchanged behaviour (just records)
    - ``'full'``: ``svd_solver='full'``
    - ``'seed'``: ``random_state=0`` (solver left to sklearn's ``'auto'``)

To apply it to ``run_pypeit`` without code changes, put ``site/`` (which holds a
``sitecustomize.py``) on ``PYTHONPATH`` and set ``I2208_PCA_FIX`` and/or
``I2208_PCA_CAPTURE``.
"""
import functools
import inspect
import os
from pathlib import Path

import numpy as np

# One entry per pca_decomposition call: dict(shape, npca, solvers, caller)
CALLS = []


def install(fix=None, capture_dir=None):
    from sklearn.decomposition import PCA as SklearnPCA
    from pypeit.core import pca as pcamod

    if getattr(pcamod, '_i2208_installed', False):
        return
    created = []

    def pca_factory(*args, **kwargs):
        if fix == 'full':
            kwargs.setdefault('svd_solver', 'full')
        elif fix == 'seed':
            kwargs.setdefault('random_state', 0)
        elif fix is not None:
            raise ValueError(f'Unknown fix: {fix}')
        obj = SklearnPCA(*args, **kwargs)
        created.append(obj)
        return obj

    orig = pcamod.pca_decomposition

    @functools.wraps(orig)
    def pca_decomposition(vectors, npca=None, pca_explained_var=99.0, mean=None):
        created.clear()
        result = orig(vectors, npca=npca, pca_explained_var=pca_explained_var, mean=mean)
        stack = [f.function for f in inspect.stack()[1:6]]
        entry = dict(i=len(CALLS), shape=list(vectors.shape), npca=npca,
                     npca_used=int(result[0].shape[1]),
                     solvers=[getattr(p, '_fit_svd_solver', None) for p in created],
                     caller='<'.join(stack))
        CALLS.append(entry)
        if capture_dir is not None:
            d = Path(capture_dir)
            d.mkdir(parents=True, exist_ok=True)
            np.savez(d / f'pca_call_{os.getpid()}_{entry["i"]:03d}.npz',
                     vectors=vectors, npca=-1 if npca is None else npca,
                     pca_explained_var=np.nan if pca_explained_var is None else pca_explained_var,
                     mean=np.array([]) if mean is None else mean,
                     coeffs=result[0], components=result[1], pca_mean=result[2],
                     solvers=np.array(entry['solvers']), caller=entry['caller'])
        return result

    pcamod.PCA = pca_factory
    pcamod.pca_decomposition = pca_decomposition
    pcamod._i2208_installed = True


def install_from_env():
    fix = os.environ.get('I2208_PCA_FIX') or None
    capture = os.environ.get('I2208_PCA_CAPTURE') or None
    if fix is None and capture is None:
        return
    install(fix=fix, capture_dir=capture)
