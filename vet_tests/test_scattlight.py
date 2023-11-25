"""
Module to run tests on arcoadd
"""
import os

import pytest
import numpy as np

from pypeit import scattlight
from pypeit.core import scattlight as core_scattlight
from pypeit import edgetrace
from pypeit.spectrographs.util import load_spectrograph

from IPython import embed

import warnings
warnings.simplefilter("ignore", UserWarning)


def scattlight_resid(droot, par, tol=0.1):
    """ Calculate the residual, and check that it is within the allowed tolerance.

    Args:
        droot (:obj:`str`):
            The root directory that points to the calibrations folder
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            Parameters required by all of PypeIt methods
        tol (:obj:`float`, optional):
            Allowed tolerance (default = 0.1 = 10% tolerance)
    """
    # Load the calibration frames
    scattlightframe = scattlight.ScatteredLight.from_file(os.path.join(droot, 'ScatteredLight_A_0_DET01.fits.gz'))
    slits = edgetrace.EdgeTraceSet.from_file(os.path.join(droot, 'Edges_A_0_DET01.fits.gz')).get_slits()

    # Get the central trace of each slit
    left, right, _ = slits.select_edges(initial=True, flexure=None)
    centrace = 0.5 * (left + right)

    # Find all pixels off the slit
    spatbin = slits.binspat
    pad = par['scienceframe']['process']['scattlight']['finecorr_pad'] // spatbin
    mask_regions = par['scienceframe']['process']['scattlight']['finecorr_mask']
    offslitmask = slits.slit_img(pad=pad, initial=True, flexure=None) == -1
    # Mask the user-specified inter-slit regions
    offslitmask = core_scattlight.mask_slit_regions(offslitmask, centrace, mask_regions=mask_regions)

    # Calculate the residual correction
    ww = np.where(offslitmask)
    resid = (scattlightframe.scattlight_raw - scattlightframe.scattlight_model)/scattlightframe.scattlight_model

    # Calculate some statistics
    medval = np.median(resid[ww])
    stdval = 1.4826 * np.median(np.abs(resid[ww]-medval))

    # Check the scattered light tolerance
    assert(np.abs(medval) < stdval)  # Check that the median value is within 1 standard deviation of 0
    assert(stdval < tol)  # Check that the correction is better than 10%


def test_scattlight_keckesi(redux_out):
    """ Calculate the residuals of the scattered light subtraction for Keck ESI"""
    droot = os.path.join(redux_out,
                         'keck_esi',
                         'Ech_1x1',
                         'Calibrations')

    # Grab the spectrograph and parset
    spec = load_spectrograph("keck_esi")
    par = spec.default_pypeit_par()

    # Check the residuals
    scattlight_resid(droot, par, tol=0.1)


def test_scattlight_keckkcwi(redux_out):
    """ Calculate the residuals of the scattered light subtraction for Keck/KCWI """
    droot = os.path.join(redux_out,
                         'keck_kcwi',
                         'small_bh2_4200',
                         'Calibrations')

    # Grab the spectrograph and parset
    spec = load_spectrograph("keck_kcwi")
    par = spec.default_pypeit_par()

    # Check the residuals
    scattlight_resid(droot, par, tol=0.1)
