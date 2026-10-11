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


#: Setup used for the comparison: cheap, and genuinely two-detector.
#: NOTE: most dev-suite LRIS-blue setups pin ``detnum`` to a single detector in
#: their .pypeit file, which would defeat the point of this test; this one does
#: not, so it reduces both DET01 and DET02.
INSTRUMENT = 'keck_lris_blue'
SETUP = 'multi_600_4000_slitmask'

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

#: Relative tolerance for the robust quantities.  The measured run-to-run
#: scatter in these is ~2e-5 (see the module note), and a genuine misrouting
#: bug would show up at order unity, so 1e-3 is comfortably clear of the noise
#: while still catching anything real.
RTOL = 1e-3


def _science_dirs(redux_out):
    """
    Return the serial and parallel ``Science/`` directories for the setup.

    Args:
        redux_out (:obj:`str`): The dev-suite ``REDUX_OUT`` directory.

    Returns:
        :obj:`tuple`: The serial and parallel `Path`_ objects.
    """
    root = Path(redux_out).absolute() / INSTRUMENT / SETUP
    return root / 'Science', root / NCPU_SUBDIR / 'Science'


def _max_rel(a, b):
    """
    Maximum relative difference between two arrays, ignoring non-finite values.

    Args:
        a, b (`numpy.ndarray`_): Arrays to compare.

    Returns:
        :obj:`float`: The maximum relative difference, or ``inf`` if the shapes
        differ.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return np.inf
    scale = np.maximum(np.abs(a), np.abs(b))
    good = np.isfinite(a) & np.isfinite(b) & (scale > 0)
    if not np.any(good):
        return 0.0
    return float(np.max(np.abs(a[good] - b[good]) / scale[good]))


def test_ncpu_same_detectors(redux_out):
    """
    Both reductions must cover every detector.

    This is the regression for the bug where a detector whose calibrations
    failed was removed from the detector list *while that list was being
    iterated*, so the following detector was silently skipped.  With two
    detectors, that bug would leave one of them missing here.  Asserted
    exactly: it does not depend on floating-point noise.
    """
    serial_dir, par_dir = _science_dirs(redux_out)
    assert serial_dir.is_dir(), f'serial Science directory missing: {serial_dir}'
    assert par_dir.is_dir(), f'--ncpu Science directory missing: {par_dir}'

    dets = {}
    for tag, sdir in (('serial', serial_dir), ('--ncpu', par_dir)):
        spec2d = sorted(sdir.glob('spec2d_*.fits'))
        assert len(spec2d) > 0, f'{tag} reduction produced no spec2d files'
        with fits.open(spec2d[0]) as hdu:
            dets[tag] = {h.name.split('-')[0] for h in hdu if '-' in h.name}
        assert len(dets[tag]) == 2, \
                f'{tag} reduction covered detectors {sorted(dets[tag])}; expected both of them'
    assert dets['serial'] == dets['--ncpu'], \
            f"detectors differ: serial {sorted(dets['serial'])} vs " \
            f"--ncpu {sorted(dets['--ncpu'])}"


def test_ncpu_same_structure(redux_out):
    """
    The two reductions must produce the same files, objects and extensions.

    All asserted exactly.  These are the properties detector parallelism could
    plausibly break: results assembled in completion order rather than detector
    order would reorder the objects, and a payload routed to the wrong worker
    would put objects on the wrong detector.
    """
    serial_dir, par_dir = _science_dirs(redux_out)

    serial_1d = sorted(serial_dir.glob('spec1d_*.fits'))
    assert len(serial_1d) > 0, f'no serial spec1d files in {serial_dir}'

    for sfile in serial_1d:
        pfile = par_dir / sfile.name
        assert pfile.exists(), f'--ncpu reduction did not produce {pfile.name}'

        s = specobjs.SpecObjs.from_fitsfile(str(sfile), chk_version=False)
        p = specobjs.SpecObjs.from_fitsfile(str(pfile), chk_version=False)
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


def test_ncpu_wavelengths_agree(redux_out):
    """
    The wavelength solution must agree between the two reductions.

    This is the sharpest available check that each detector was reduced with
    its own calibrations: the wavelengths are robust against Issue 2208 (the
    serial/parallel difference equals the serial/serial difference to three
    significant figures), while a payload routed to the wrong worker would
    change them by order unity.
    """
    serial_dir, par_dir = _science_dirs(redux_out)

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
                rel = _max_rel(a, b)
                if rel > RTOL:
                    bad.append(f'{sfile.name} object {i} {key}: max relative difference '
                               f'{rel:.3e} (> {RTOL:.0e})')

    assert ncompared > 0, 'no spec1d wavelength arrays were actually compared'
    assert len(bad) == 0, 'serial and --ncpu wavelengths disagree:\n  ' + '\n  '.join(bad[:20])


def test_ncpu_spec2d_geometry_agrees(redux_out):
    """
    The spec2d calibration and input images must agree between the reductions.

    ``WAVEIMG`` and ``TILTS`` come from the calibrations and ``SCIIMG`` and
    ``IVARRAW`` from the image processing, i.e., from stages 1 and 2 of the
    parallelised reduction.  All four are robust against Issue 2208, so this
    checks those two stages directly.
    """
    serial_dir, par_dir = _science_dirs(redux_out)

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
                rel = _max_rel(sext.data, pext.data)
                if rel > RTOL:
                    bad.append(f'{sfile.name} {sext.name}: max relative difference '
                               f'{rel:.3e} (> {RTOL:.0e})')

    assert ncompared > 0, 'no spec2d geometry extensions were actually compared'
    assert len(bad) == 0, \
            'serial and --ncpu spec2d calibration/input images disagree:\n  ' \
            + '\n  '.join(bad[:20])


@pytest.mark.xfail(reason='PypeIt is not reproducible run-to-run on multi-slit data; '
                          'see https://github.com/pypeit/PypeIt/issues/2208',
                   strict=False)
def test_ncpu_flux_arrays_identical(redux_out):
    """
    The extracted fluxes and sky/object models must be identical.

    This is what detector parallelism is *designed* to deliver, and it holds on
    deterministic data (``shane_kast_blue`` is bit-identical between
    ``--ncpu 1`` and ``--ncpu 2``).  It cannot hold here while Issue 2208
    stands, so this is marked ``xfail``; when that issue is fixed this test
    should start passing and the marker should be removed.
    """
    serial_dir, par_dir = _science_dirs(redux_out)

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
