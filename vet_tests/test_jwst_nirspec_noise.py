"""
Vet the noise model of the JWST/NIRSpec reductions.

The errors of the reduced data are vetted with the residuals of the data from
the sky and object models, normalized by the model errors: if the errors are
correct, the normalized residuals have unit standard deviation (chi = 1).  In
the 2D spectra, the residuals are ``sciimg - skymodel - objmodel``, normalized
by ``ivarmodel``; in the 1D spectra, they are the differences of the
extracted counts from their running median, normalized by their errors.  To
measure the noise rather than the errors of the sky and object models or the
spectral features, the residuals are high-pass filtered by subtracting their
running median along the spectral direction, and the robust standard deviation
of the filtered residuals is corrected for the filtering.  In the 2D spectra,
the pixels with significant source signal are also excluded: the spectral
structure of bright sources (including the negative traces of the nod
subtraction) is not removed by the high-pass filter and would be counted as
noise, especially in the 2D coadds, which have a higher S/N.
"""
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage
from scipy import stats

from pypeit import spec2dobj
from pypeit import specobjs

# Length in pixels of the running median along the spectral direction
NMED = 21

# Pixels whose sky-subtracted signal, smoothed over NMED spectral pixels and 3
# spatial pixels, has an absolute S/N above this threshold are excluded from
# the 2D chi
SIGNAL_SNR = 2.

# JWST/NIRSpec setups of the dev suite, and those with a 2D coadd
SETUPS = ['S200A2_G395M_F290LP', 'MSA_PRISM_nomosaic', 'MSA_PRISM_2visits',
          'MSA_G395M_F290LP_mosaic']
COADD_SETUPS = ['MSA_PRISM_nomosaic', 'MSA_PRISM_2visits', 'MSA_G395M_F290LP_mosaic']

# Bounds of the 2D chi of each exposure and 2D coadd.  With the JWST/NIRSpec
# noise model of PypeIt 3.1 (shot noise from the counts and the excess
# readnoise measured on the exposure differences), the 2D chi is 1.00-1.07 for
# the exposures and 1.04-1.12 for the coadds; before, it was 1.15-1.46 for the
# exposures and 1.28-1.50 for the coadds.
CHI_2D_BOUNDS = (0.9, 1.2)

# Bounds of the median 1D chi of the objects of each setup.  The 1D chi is
# noisier (fewer pixels) and, for the brighter sources, is increased by the
# spectral features narrower than the running median; the median 1D chi is
# 1.06-1.19, and it was 1.26-1.63 before.
CHI_1D_BOUNDS = (0.9, 1.25)


def _robust_std(x):
    """Standard deviation from the median absolute deviation."""
    return stats.median_abs_deviation(x, scale='normal')


def _highpass_correction():
    """
    Ratio of the standard deviation of white noise to the robust standard
    deviation of the noise high-pass filtered by subtracting its running
    median, measured on simulated noise.
    """
    rng = np.random.default_rng(seed=1)
    noise = rng.normal(size=(4000, 250))
    filtered = noise - ndimage.median_filter(noise, size=(NMED, 1), mode='nearest')
    return 1 / _robust_std(filtered.ravel())


def _chi_2d(spec2d_file, correction):
    """
    Chi of the 2D residuals of each detector (or mosaic) in a spec2d file,
    excluding the pixels with significant source signal (see SIGNAL_SNR).
    """
    allspec2d = spec2dobj.AllSpec2DObj.from_fits(str(spec2d_file), chk_version=False)
    # For difference (A-B) images, the sky is removed by the subtraction and
    # the sky model only fits residual structure, so the source signal is the
    # image itself; otherwise, it is the image minus the sky model
    bkg_redux = allspec2d['meta']['bkg_redux']
    chi = []
    for det in allspec2d.detectors:
        s = allspec2d[det]
        gpm = (s.bpmmask.mask == 0) & (s.ivarmodel > 0) & (s.slits.slit_img() >= 0)
        # Pixels with significant (positive or negative) source signal
        signal = s.sciimg if bkg_redux else s.sciimg - s.skymodel
        snr = np.where(gpm, signal * np.sqrt(s.ivarmodel), 0.)
        signal = np.absolute(ndimage.median_filter(snr, size=(NMED, 3), mode='nearest')) \
                    > SIGNAL_SNR
        model = s.skymodel if s.objmodel is None else s.skymodel + s.objmodel
        resid = np.where(gpm, (s.sciimg - model) * np.sqrt(s.ivarmodel), 0.)
        resid -= ndimage.median_filter(resid, size=(NMED, 1), mode='nearest')
        chi.append(_robust_std(resid[gpm & ~signal]) * correction)
    return chi


def _chi_1d(spec1d_file, correction):
    """
    Chi of the optimally extracted spectra of the objects in a spec1d file,
    with their median S/N.
    """
    chi = []
    for sobj in specobjs.SpecObjs.from_fitsfile(str(spec1d_file), chk_version=False):
        if sobj.OPT_COUNTS is None:
            continue
        gpm = sobj.OPT_MASK & (sobj.OPT_COUNTS_IVAR > 0)
        if np.sum(gpm) < 3*NMED:
            continue
        counts = sobj.OPT_COUNTS[gpm]
        resid = (counts - ndimage.median_filter(counts, size=NMED, mode='nearest')) \
                    * np.sqrt(sobj.OPT_COUNTS_IVAR[gpm])
        chi.append((sobj.NAME, sobj.S2N, _robust_std(resid) * correction))
    return chi


@pytest.mark.parametrize('setup', SETUPS)
def test_chi_2d(redux_out, setup):
    # The errors of the 2D spectra of each exposure must match the scatter of
    # the data around the sky and object models
    sci_dir = Path(redux_out).absolute() / 'jwst_nirspec' / setup / 'Science'
    spec2d_files = sorted(sci_dir.glob('spec2d*.fits'))
    assert len(spec2d_files) > 0, f'No spec2d files found in {sci_dir}'
    correction = _highpass_correction()
    for f in spec2d_files:
        for chi in _chi_2d(f, correction):
            assert CHI_2D_BOUNDS[0] < chi < CHI_2D_BOUNDS[1], \
                f'The 2D chi of {f.name} is {chi:.3f}: the errors are ' \
                + ('overestimated.' if chi < CHI_2D_BOUNDS[0] else 'underestimated.')


@pytest.mark.parametrize('setup', COADD_SETUPS)
def test_chi_2d_coadd(redux_out, setup):
    # The errors of the 2D coadds must match the scatter of the data around the
    # sky and object models
    coadd_dir = Path(redux_out).absolute() / 'jwst_nirspec' / setup / 'Science_coadd'
    spec2d_files = sorted(coadd_dir.glob('spec2d*.fits'))
    assert len(spec2d_files) > 0, f'No spec2d files found in {coadd_dir}'
    correction = _highpass_correction()
    for f in spec2d_files:
        for chi in _chi_2d(f, correction):
            assert CHI_2D_BOUNDS[0] < chi < CHI_2D_BOUNDS[1], \
                f'The 2D chi of the coadd {f.name} is {chi:.3f}: the errors are ' \
                + ('overestimated.' if chi < CHI_2D_BOUNDS[0] else 'underestimated.')


@pytest.mark.parametrize('setup', SETUPS)
def test_chi_1d(redux_out, setup):
    # The errors of the optimally extracted 1D spectra must match their scatter
    # (median over the objects of the setup)
    sci_dir = Path(redux_out).absolute() / 'jwst_nirspec' / setup / 'Science'
    spec1d_files = sorted(sci_dir.glob('spec1d*.fits'))
    assert len(spec1d_files) > 0, f'No spec1d files found in {sci_dir}'
    correction = _highpass_correction()
    chi = [c for f in spec1d_files for _, _, c in _chi_1d(f, correction)]
    assert len(chi) > 0, f'No extracted spectra found in {sci_dir}'
    median_chi = np.median(chi)
    assert CHI_1D_BOUNDS[0] < median_chi < CHI_1D_BOUNDS[1], \
        f'The median 1D chi of {setup} is {median_chi:.3f} (values: ' \
        + ', '.join(f'{c:.2f}' for c in chi) + '): the errors are ' \
        + ('overestimated.' if median_chi < CHI_1D_BOUNDS[0] else 'underestimated.')
