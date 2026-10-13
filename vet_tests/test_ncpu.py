"""
Module to verify that reducing with multiple CPUs changes nothing that matters.

PypeIt distributes the calibration build, the science-image processing and the
extraction over one worker process per detector when ``[rdx] ncpu > 1``.  This
module checks that doing so does not change the reduction.

.. note::

    These tests deliberately do **not** require bit-identical products, even
    though detector parallelism is designed to be exact.  PypeIt is currently
    not reproducible run-to-run on multi-slit data at all: edge centroiding
    returns a value differing by one float32 ULP from bit-identical input, and
    the sigma-clipped iterative fits downstream amplify that into visibly
    different sky models.  See `Issue 2208`_.  Two *serial* runs of this setup
    differ by as much as a serial/parallel pair does, so an exact comparison
    here would be a test of that bug rather than of the parallelism.

    What *is* asserted here is chosen by measurement rather than by guesswork.
    The structural properties are asserted exactly, because they are what the
    parallelism could plausibly break -- a dropped detector, a reordered object
    list, a missing extension -- and they do not depend on floating-point noise.
    The numerical comparisons are restricted to the quantities measured to be
    robust against `Issue 2208`_, for which the serial/parallel difference is
    the same size as the serial/serial difference:

    ========================  ==================  =================
    quantity                  serial vs. serial   serial vs. ncpu 2
    ========================  ==================  =================
    ``spec1d`` ``OPT_WAVE``   2.103e-05           2.103e-05
    ``spec2d`` ``WAVEIMG``    1.114e-10           1.318e-10
    ``spec2d`` ``TILTS``      2.129e-05           1.702e-05
    ``spec2d`` ``SCIIMG``     3.687e-07           1.192e-07
    ========================  ==================  =================

    These validate that every detector received the right calibrations and the
    right slice of the payload, which is what detector parallelism could get
    wrong; misrouting would change them by order unity, not by 1e-5.

    The flux-like arrays (``OPT_COUNTS``, ``SKYMODEL``, ``OBJMODEL``) cannot be
    compared numerically at all while `Issue 2208`_ stands: two *serial* runs of
    this setup differ by up to **350 sigma** in ``OPT_COUNTS``, so any threshold
    loose enough to pass would be far too loose to catch a real bug.  They are
    covered by :func:`test_ncpu_flux_arrays_identical`, which is marked
    ``xfail`` and will start passing -- and should then be un-marked -- once
    `Issue 2208`_ is fixed.

.. _Issue 2208: https://github.com/pypeit/PypeIt/issues/2208

The two reductions are produced by the dev suite's ``reduce`` phase: the
ordinary reduce test writes the serial products into the setup's redux
directory, and :class:`~test_scripts.pypeit_tests.PypeItReduceNcpuTest` writes
the ``--ncpu 2`` products into the ``ncpu`` subdirectory beneath it.
"""
from pathlib import Path

from IPython import embed

import numpy as np
import pytest

from astropy.io import fits

from pypeit import specobjs


#: Fallback registration, used only if the dev-suite one cannot be loaded.
#: Kept in step by :func:`test_ncpu_registration_matches_dev_suite`.
_FALLBACK_SETUPS = [('keck_deimos', '830G_M_8500'),
                    ('keck_lris_blue', 'multi_600_4000_slitmask'),
                    ('mmt_binospec', 'Longslit_G1000')]


def _load_ncpu_reduce():
    """
    Load ``_ncpu_reduce`` from the dev suite's ``test_scripts/test_setups.py``.

    ``test_scripts`` cannot simply be imported from here: it has no
    ``__init__.py``, so it is a *namespace* package, and ``vet_tests`` contains
    a regular module also called ``test_scripts`` -- which wins the import
    regardless of ``sys.path`` order.  The file is therefore loaded directly,
    under a private name, with a synthetic parent package so that its relative
    imports still resolve.

    Returns
    -------
    :obj:`dict` or None
        The ``_ncpu_reduce`` mapping, or None if it could not be loaded.
    """
    import importlib.util
    import sys
    import types

    pkg_dir = Path(__file__).resolve().parent.parent / 'test_scripts'
    name = '_devsuite_test_scripts'
    try:
        if name not in sys.modules:
            pkg = types.ModuleType(name)
            pkg.__path__ = [str(pkg_dir)]
            sys.modules[name] = pkg
        spec = importlib.util.spec_from_file_location(f'{name}.test_setups',
                                                      pkg_dir / 'test_setups.py')
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        return mod._ncpu_reduce
    except Exception:
        return None


def _registered_setups():
    """
    The ``(instrument, setup)`` pairs reduced with ``--ncpu``.

    Taken from the dev-suite registration where possible, so this module cannot
    drift out of step with what the suite actually runs; see
    :func:`_load_ncpu_reduce`.

    Returns
    -------
    :obj:`list`
        The ``(instrument, setup)`` tuples, sorted.
    """
    reg = _load_ncpu_reduce()
    if reg is None:
        return sorted(_FALLBACK_SETUPS)
    return sorted((instr, setup) for instr, setups in reg.items() for setup in setups)


#: The setups compared by this module; see :func:`_registered_setups`.
NCPU_SETUPS = _registered_setups()

#: Subdirectory holding the parallel reduction; see
#: ``test_scripts.pypeit_tests.NCPU_SUBDIR``.
NCPU_SUBDIR = 'ncpu'

#: spec1d columns that are robust against Issue 2208 (see the module note).
SPEC1D_ROBUST_KEYS = ('OPT_WAVE', 'BOX_WAVE')

#: spec2d extensions that are robust against Issue 2208.
SPEC2D_ROBUST_EXT = ('WAVEIMG', 'TILTS', 'SCIIMG', 'IVARRAW')

#: spec1d/spec2d arrays that Issue 2208 makes untestable numerically.
FLUX_KEYS = ('OPT_COUNTS', 'BOX_COUNTS')
FLUX_EXT = ('SKYMODEL', 'OBJMODEL')

#: Quantile of the relative-difference distribution that is compared; see
#: :func:`_rel_quantile` for why this is not a maximum.
QUANTILE = 0.999

#: Tolerance on that quantile.  The largest value measured across all three
#: setups is 8.2e-05 (``DET07-IVARRAW`` on ``keck_deimos/830G_M_8500``), and a
#: genuine misrouting bug would show up at order unity, so 1e-3 sits an order of
#: magnitude above the noise and four below a real failure.
RTOL = 1e-3

#: Permitted change in an array's non-zero count; see :func:`_nonzero_delta`.
#: A handful of slit-edge pixels may move (Issue 2208); a zeroed slit would not
#: be anywhere near this small.
MAX_NONZERO_DELTA = 64


def _science_dirs(redux_out, instrument, setup):
    """
    Return the serial and parallel ``Science/`` directories for a setup.

    Args:
        redux_out (:obj:`str`): The dev-suite ``REDUX_OUT`` directory.
        instrument (:obj:`str`): Instrument name.
        setup (:obj:`str`): Setup name.

    Returns:
        :obj:`tuple`: The serial and parallel `Path`_ objects.
    """
    root = Path(redux_out).absolute() / instrument / setup
    return root / 'Science', root / NCPU_SUBDIR / 'Science'


def _rel_quantile(a, b, q=QUANTILE):
    """
    A high quantile of the relative difference between two arrays.

    A *maximum* is the wrong statistic for these comparisons.  Issue 2208 can
    shift a slit edge by a fraction of a pixel, which moves a single pixel in or
    out of the illuminated region; that pixel is then zero in one reduction and
    finite in the other, giving a relative difference of exactly 1.0 and
    dominating the maximum no matter how well the other eight million pixels
    agree.  Measured on ``keck_deimos/830G_M_8500``, DET07 has exactly one such
    pixel in ``WAVEIMG`` and ``TILTS``, while DET03 has none.

    A high quantile ignores that handful of boundary pixels and still catches
    anything real: a detector reduced with the wrong calibrations, or a payload
    delivered to the wrong worker, disagrees across essentially the whole array,
    so its quantile is of order unity rather than the 1e-4 seen here.
    :func:`_nonzero_delta` covers the complementary case of a region being
    zeroed wholesale.

    Args:
        a, b (`numpy.ndarray`_): Arrays to compare.
        q (:obj:`float`, optional): Quantile to report.

    Returns:
        :obj:`float`: The quantile of the relative difference, or ``inf`` if the
        shapes differ.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return np.inf
    scale = np.maximum(np.abs(a), np.abs(b))
    good = np.isfinite(a) & np.isfinite(b) & (scale > 0)
    if not np.any(good):
        return 0.0
    return float(np.quantile(np.abs(a[good] - b[good]) / scale[good], q))


def _nonzero_delta(a, b):
    """
    Difference in the number of non-zero elements between two arrays.

    Guards the blind spot of :func:`_rel_quantile`: if a whole slit or detector
    were zeroed out, most pixels would agree trivially (zero against zero) and
    the quantile would stay small, but the non-zero count would move sharply.

    Args:
        a, b (`numpy.ndarray`_): Arrays to compare.

    Returns:
        :obj:`int`: The absolute difference in non-zero counts.
    """
    return int(abs(np.count_nonzero(np.asarray(a)) - np.count_nonzero(np.asarray(b))))


@pytest.mark.parametrize('instrument,setup', NCPU_SETUPS)
def test_ncpu_same_detectors(redux_out, instrument, setup):
    """
    Both reductions must cover every detector.

    This is the regression for the bug where a detector whose calibrations
    failed was removed from the detector list *while that list was being
    iterated*, so the following detector was silently skipped; that bug would
    leave one missing here.  Asserted exactly: it does not depend on
    floating-point noise.

    It also guards against the test itself going vacuous.  A setup whose
    ``.pypeit`` file selects a single detector or a single mosaic -- note that
    ``detnum = (3,7)`` is one *mosaic* while ``detnum = 3,7`` is two detectors
    -- gives the reduction nothing to distribute, so ``--ncpu 2`` would quietly
    be a no-op and every other test in this module would pass without testing
    anything.  Requiring at least two therefore keeps the registration honest.
    """
    serial_dir, par_dir = _science_dirs(redux_out, instrument, setup)
    assert serial_dir.is_dir(), f'serial Science directory missing: {serial_dir}'
    assert par_dir.is_dir(), f'--ncpu Science directory missing: {par_dir}'

    dets = {}
    for tag, sdir in (('serial', serial_dir), ('--ncpu', par_dir)):
        spec2d = sorted(sdir.glob('spec2d_*.fits'))
        assert len(spec2d) > 0, f'{tag} reduction produced no spec2d files'
        with fits.open(spec2d[0]) as hdu:
            dets[tag] = {h.name.split('-')[0] for h in hdu if '-' in h.name}
        assert len(dets[tag]) > 1, \
                f'{tag} reduction of {instrument}/{setup} covered only ' \
                f'{sorted(dets[tag])}; a single detector means --ncpu has nothing to ' \
                f'distribute and this comparison is vacuous'
    assert dets['serial'] == dets['--ncpu'], \
            f"{instrument}/{setup}: detectors differ between the reductions -- serial " \
            f"{sorted(dets['serial'])} vs --ncpu {sorted(dets['--ncpu'])}"


@pytest.mark.parametrize('instrument,setup', NCPU_SETUPS)
def test_ncpu_same_structure(redux_out, instrument, setup):
    """
    The two reductions must produce the same files, objects and extensions.

    All asserted exactly.  These are the properties detector parallelism could
    plausibly break: results assembled in completion order rather than detector
    order would reorder the objects, and a payload routed to the wrong worker
    would put objects on the wrong detector.
    """
    serial_dir, par_dir = _science_dirs(redux_out, instrument, setup)

    serial_1d = sorted(serial_dir.glob('spec1d_*.fits'))
    assert len(serial_1d) > 0, f'no serial spec1d files in {serial_dir}'

    for sfile in serial_1d:
        pfile = par_dir / sfile.name
        assert pfile.exists(), f'--ncpu reduction did not produce {pfile.name}'

        s = specobjs.SpecObjs.from_fitsfile(str(sfile), chk_version=False)
        p = specobjs.SpecObjs.from_fitsfile(str(pfile), chk_version=False)

        # Catch stale accumulation before reporting a count mismatch.
        # save_exposure merges detectors into an existing spec1d rather than
        # replacing it, so reducing into a directory that already holds
        # products duplicates every object.  That is leftover state, not a
        # parallelism defect, and saying so saves the next person the hour it
        # cost to work out the first time.
        for tag, obj, path in (('serial', s, sfile), ('--ncpu', p, pfile)):
            if obj.nobj == 0:
                continue
            names = list(obj.NAME)
            assert len(names) == len(set(names)), \
                    f'{path.name}: the {tag} reduction contains {len(names)} objects but only ' \
                    f'{len(set(names))} distinct names, i.e. every object appears more than ' \
                    f'once.  This is stale state, not a parallelism failure: the reduction ' \
                    f'was run into a directory that already held products.  Remove ' \
                    f'{path.parent} and reduce again.'

        assert s.nobj == p.nobj, \
                f'{sfile.name}: serial found {s.nobj} objects, --ncpu found {p.nobj}'
        if s.nobj == 0:
            continue
        # Same objects, on the same detectors, in the same order.
        assert list(s.DET) == list(p.DET), \
                f'{sfile.name}: objects are on different detectors, or in a different order'
        assert list(s.NAME) == list(p.NAME), \
                f'{sfile.name}: object names differ, or are in a different order'

    for sfile in sorted(serial_dir.glob('spec2d_*.fits')):
        pfile = par_dir / sfile.name
        assert pfile.exists(), f'--ncpu reduction did not produce {pfile.name}'
        with fits.open(sfile) as sh, fits.open(pfile) as ph:
            assert [h.name for h in sh] == [h.name for h in ph], \
                    f'{sfile.name}: the two reductions wrote different extensions'


@pytest.mark.parametrize('instrument,setup', NCPU_SETUPS)
def test_ncpu_wavelengths_agree(redux_out, instrument, setup):
    """
    The wavelength solution must agree between the two reductions.

    This is the sharpest available check that each detector was reduced with
    its own calibrations: the wavelengths are robust against Issue 2208 (the
    serial/parallel difference equals the serial/serial difference to three
    significant figures), while a payload routed to the wrong worker would
    change them by order unity.
    """
    serial_dir, par_dir = _science_dirs(redux_out, instrument, setup)

    bad = []
    ncompared = 0
    for sfile in sorted(serial_dir.glob('spec1d_*.fits')):
        s = specobjs.SpecObjs.from_fitsfile(str(sfile), chk_version=False)
        p = specobjs.SpecObjs.from_fitsfile(str(par_dir / sfile.name), chk_version=False)
        if s.nobj == 0:
            continue
        for i in range(s.nobj):
            for key in SPEC1D_ROBUST_KEYS:
                a = getattr(s[i], key, None)
                b = getattr(p[i], key, None)
                if a is None or b is None:
                    assert a is None and b is None, \
                            f'{sfile.name} object {i}: {key} present in only one reduction'
                    continue
                ncompared += 1
                rel = _rel_quantile(a, b)
                if rel > RTOL:
                    bad.append(f'{sfile.name} object {i} {key}: {QUANTILE:.1%} quantile of '
                               f'relative difference is {rel:.3e} (> {RTOL:.0e})')

    assert ncompared > 0, 'no spec1d wavelength arrays were actually compared'
    assert len(bad) == 0, 'serial and --ncpu wavelengths disagree:\n  ' + '\n  '.join(bad[:20])


@pytest.mark.parametrize('instrument,setup', NCPU_SETUPS)
def test_ncpu_spec2d_geometry_agrees(redux_out, instrument, setup):
    """
    The spec2d calibration and input images must agree between the reductions.

    ``WAVEIMG`` and ``TILTS`` come from the calibrations and ``SCIIMG`` and
    ``IVARRAW`` from the image processing, i.e., from stages 1 and 2 of the
    parallelised reduction.  All four are robust against Issue 2208, so this
    checks those two stages directly.
    """
    serial_dir, par_dir = _science_dirs(redux_out, instrument, setup)

    bad = []
    ncompared = 0
    for sfile in sorted(serial_dir.glob('spec2d_*.fits')):
        with fits.open(sfile) as sh, fits.open(par_dir / sfile.name) as ph:
            for sext, pext in zip(sh, ph):
                if sext.data is None or not isinstance(sext.data, np.ndarray) \
                        or sext.data.dtype.kind not in 'fiu':
                    continue
                if not any(sext.name.endswith(e) for e in SPEC2D_ROBUST_EXT):
                    continue
                ncompared += 1
                rel = _rel_quantile(sext.data, pext.data)
                if rel > RTOL:
                    bad.append(f'{sfile.name} {sext.name}: {QUANTILE:.1%} quantile of '
                               f'relative difference is {rel:.3e} (> {RTOL:.0e})')
                nzd = _nonzero_delta(sext.data, pext.data)
                if nzd > MAX_NONZERO_DELTA:
                    bad.append(f'{sfile.name} {sext.name}: non-zero element count differs by '
                               f'{nzd} (> {MAX_NONZERO_DELTA}), suggesting a region was '
                               f'zeroed rather than merely perturbed')

    assert ncompared > 0, 'no spec2d geometry extensions were actually compared'
    assert len(bad) == 0, \
            'serial and --ncpu spec2d calibration/input images disagree:\n  ' \
            + '\n  '.join(bad[:20])


@pytest.mark.xfail(reason='PypeIt is not reproducible run-to-run on multi-slit data; '
                          'see https://github.com/pypeit/PypeIt/issues/2208',
                   strict=False)
@pytest.mark.parametrize('instrument,setup', NCPU_SETUPS)
def test_ncpu_flux_arrays_identical(redux_out, instrument, setup):
    """
    The extracted fluxes and sky/object models must be identical.

    This is what detector parallelism is *designed* to deliver, and it holds on
    deterministic data (``shane_kast_blue`` is bit-identical between
    ``--ncpu 1`` and ``--ncpu 2``).  It cannot hold here while Issue 2208
    stands, so this is marked ``xfail``; when that issue is fixed this test
    should start passing and the marker should be removed.
    """
    serial_dir, par_dir = _science_dirs(redux_out, instrument, setup)

    for sfile in sorted(serial_dir.glob('spec1d_*.fits')):
        s = specobjs.SpecObjs.from_fitsfile(str(sfile), chk_version=False)
        p = specobjs.SpecObjs.from_fitsfile(str(par_dir / sfile.name), chk_version=False)
        for i in range(s.nobj):
            for key in FLUX_KEYS:
                a = getattr(s[i], key, None)
                b = getattr(p[i], key, None)
                if a is None or b is None:
                    continue
                np.testing.assert_array_equal(
                        np.asarray(a), np.asarray(b),
                        err_msg=f'{sfile.name} object {i}: {key} differs')

    for sfile in sorted(serial_dir.glob('spec2d_*.fits')):
        with fits.open(sfile) as sh, fits.open(par_dir / sfile.name) as ph:
            for sext, pext in zip(sh, ph):
                if sext.data is None or not isinstance(sext.data, np.ndarray) \
                        or sext.data.dtype.kind not in 'fiu':
                    continue
                if not any(sext.name.endswith(e) for e in FLUX_EXT):
                    continue
                np.testing.assert_array_equal(
                        sext.data, pext.data,
                        err_msg=f'{sfile.name} {sext.name} differs')


def test_ncpu_registration_matches_dev_suite():
    """
    The fallback list must match the dev suite's actual registration.

    If someone adds an ``--ncpu`` setup to ``_ncpu_reduce`` the parametrisation
    picks it up automatically, but the fallback would silently go stale; this
    catches that.
    """
    reg = _load_ncpu_reduce()
    if reg is None:
        pytest.skip('dev-suite registration could not be loaded')
    registered = sorted((instr, setup) for instr, setups in reg.items() for setup in setups)
    assert registered == sorted(_FALLBACK_SETUPS), \
            f'vet_tests/test_ncpu.py fallback list is stale: dev suite registers ' \
            f'{registered}, fallback has {sorted(_FALLBACK_SETUPS)}'
