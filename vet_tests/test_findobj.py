import numpy as np
import pytest
from pathlib import Path

from pypeit import specobjs, spec2dobj
from pypeit.core import findobj_skymask
from pypeit.spectrographs.util import load_spectrograph


@pytest.fixture
def esi_data(redux_out):
    """
    Load the ESI Ech_1x1 spec1d and spec2d for the BD+28d4211 exposure and
    return the objects most useful for testing:

    - ``sobjs``:       SpecObjs for ECH_OBJID==1 (one object per order)
    - ``order_vec``:   sorted array of echelle order numbers
    - ``slit_left``:   tweaked left slit edges  (nspec, norders)
    - ``slit_right``:   tweaked right slit edges (nspec, norders)
    - ``slit_spat_id``: integer spatial ID per slit
    - ``image``:       sky-subtracted science image from the spec2d
    - ``ivar``:        inverse-variance model from the spec2d
    - ``slitmask``:    integer slit image (-1 off slit) from the spec2d
    - ``spec_min_max``: (2, norders) spectral range array
    """
    sci_dir = Path(redux_out) / 'keck_esi' / 'Ech_1x1' / 'Science'
    spec1d_file = 'spec1d_e221024_0040-BD+28d4211_ESI_20221024T044338.293.fits'
    spec2d_file = 'spec2d_e221024_0040-BD+28d4211_ESI_20221024T044338.293.fits'

    # --- spec1d: one SpecObj per order for ECH_OBJID==1 ---
    sobjs_all = specobjs.SpecObjs.from_fitsfile(sci_dir / spec1d_file)
    sobjs = sobjs_all[sobjs_all.ECH_OBJID == 1]
    order_vec = np.sort(np.unique(sobjs.ECH_ORDER))

    # --- spec2d: sky-subtracted image + inverse variance + slit geometry ---
    s2d = spec2dobj.Spec2DObj.from_file(sci_dir / spec2d_file, 'DET01')
    slits = s2d.slits
    # Prefer tweaked slit edges (used during the actual reduction)
    slit_left = slits.left_tweak  if slits.left_tweak  is not None else slits.left
    slit_right = slits.right_tweak if slits.right_tweak is not None else slits.right
    slit_spat_id = slits.spat_id
    slitmask = slits.slit_img()
    nspec = s2d.sciimg.shape[0]
    norders = slits.nslits
    # specmin/specmax are ±inf in the file → use full spectral range
    spec_min_max = np.array([[0] * norders, [nspec - 1] * norders])

    return dict(
        sobjs=sobjs,
        order_vec=order_vec,
        slit_left=slit_left,
        slit_right=slit_right,
        slit_spat_id=slit_spat_id,
        image=s2d.sciimg,
        ivar=s2d.ivarmodel,
        slitmask=slitmask,
        spec_min_max=spec_min_max,
        norders=norders,
        binning=sobjs_all[0].DETECTOR.binning
    )


# ---------------------------------------------------------------------------
# Tests for ech_fill_in_orders — plate-scale conversion (ech_pscale_issue)
# ---------------------------------------------------------------------------

def test_fill_in_orders_plate_scale_conversion(esi_data):
    """
    Test that ech_fill_in_orders correctly scales FWHM, maskwidth, and
    BOX_R_PIX when filling in a missing order using the plate-scale ratio.

    Two calls are made with the same input (highest order missing):

    1. Uniform plate scale → pscale_conv = 1, values are unchanged.
    2. Missing order has half the plate scale of the source order →
       pscale_conv = plate_scale_ord[source] / plate_scale_ord[missing]
                   = 0.154 / 0.077 = 2.0, so values are doubled.
    """
    sobjs        = esi_data['sobjs']
    order_vec    = esi_data['order_vec']
    slit_left    = esi_data['slit_left']
    slit_right    = esi_data['slit_right']
    slit_spat_id = esi_data['slit_spat_id']
    norders      = esi_data['norders']

    # Simulate one missing order by dropping the highest order from the input
    missing_order = order_vec[-1]
    sobjs_partial = sobjs[sobjs.ECH_ORDER != missing_order]
    # obj_id must be an integer array with one element per specobj
    obj_id_arr = np.ones(len(sobjs_partial), dtype=int)
    # Source = nearest detected order (the one just below the missing order)
    src = sobjs_partial[sobjs_partial.ECH_ORDER == order_vec[-2]][0]

    # --- case 1: uniform plate scale → pscale_conv == 1 ---
    plate_scale_uniform = np.ones(norders) * 0.154
    sobjs_filled = findobj_skymask.ech_fill_in_orders(
        sobjs_partial, slit_left, slit_right, slit_spat_id,
        order_vec, plate_scale_uniform, obj_id_arr)

    filled = sobjs_filled[sobjs_filled.ECH_ORDER == missing_order]
    assert len(filled) == 1
    # With pscale_conv == 1, all three pixel-valued quantities are unchanged
    assert np.isclose(filled[0].FWHM,      src.FWHM)
    assert np.isclose(filled[0].maskwidth,  src.maskwidth)
    assert np.isclose(filled[0].BOX_R_PIX,  src.BOX_R_PIX)

    # --- case 2: missing order has half the plate scale → pscale_conv == 2 ---
    plate_scale_half = np.ones(norders) * 0.154
    plate_scale_half[-1] = 0.077
    pscale_conv = 0.154 / 0.077  # == 2.0
    sobjs_filled2 = findobj_skymask.ech_fill_in_orders(
        sobjs_partial, slit_left, slit_right, slit_spat_id,
        order_vec, plate_scale_half, obj_id_arr)

    filled2 = sobjs_filled2[sobjs_filled2.ECH_ORDER == missing_order]
    assert len(filled2) == 1
    # Pixel-valued quantities must be multiplied by pscale_conv == 2
    assert np.isclose(filled2[0].FWHM,      src.FWHM      * pscale_conv)
    assert np.isclose(filled2[0].maskwidth,  src.maskwidth  * pscale_conv)
    assert np.isclose(filled2[0].BOX_R_PIX,  src.BOX_R_PIX  * pscale_conv)


# ---------------------------------------------------------------------------
# Test for maxshift parameter (ech_pscale_issue)
# ---------------------------------------------------------------------------

def test_ech_findobj_ineach_order_maxshift(esi_data):
    """
    Verify that maxshift actually constrains the iterative trace centroiding
    in ech_findobj_ineach_order.

    Two runs are made on the same real ESI image (BD+28d4211):

    - ``maxshift=0.0``: centroids cannot move at all from the initial peak
      position, so the trace stays flat (constant in the spectral direction).
    - ``maxshift=2.0``: centroids are allowed to shift up to 2 pixels per
      step, so the trace follows the real curvature of the spectrum.

    The mean absolute difference in TRACE_SPAT between the two runs must
    exceed 0.1 pixels on every order, confirming that the parameter is wired
    in and has a measurable effect.
    """
    spec = load_spectrograph('keck_esi')
    plate_scale = spec.order_platescale(esi_data['order_vec'], binning=esi_data['binning'])
    inmask = esi_data['ivar'] > 0
    fof_link = 1.5

    # Run with no allowed centroid shift → traces remain at the initial peak
    sobjs_tight = findobj_skymask.ech_findobj_ineach_order(
        esi_data['image'], esi_data['ivar'], esi_data['slitmask'],
        esi_data['slit_left'], esi_data['slit_right'], esi_data['slit_spat_id'],
        esi_data['order_vec'], esi_data['spec_min_max'], plate_scale,
        maxshift=0.0, inmask=inmask, snr_thresh=10.0)
    obj_id = findobj_skymask.ech_fof_sobjs(
        sobjs_tight, esi_data['slit_left'], esi_data['slit_right'], esi_data['order_vec'],
        plate_scale, fof_link=fof_link
    )
    uniq_ids, cnts = np.unique(obj_id, return_counts=True)
    sobjs_tight = sobjs_tight[obj_id == uniq_ids[np.argmax(cnts)]]

    # Run with a permissive shift → traces follow the real spectral curvature
    sobjs_loose = findobj_skymask.ech_findobj_ineach_order(
        esi_data['image'], esi_data['ivar'], esi_data['slitmask'],
        esi_data['slit_left'], esi_data['slit_right'], esi_data['slit_spat_id'],
        esi_data['order_vec'], esi_data['spec_min_max'], plate_scale,
        maxshift=2.0, inmask=inmask, snr_thresh=10.0)
    obj_id = findobj_skymask.ech_fof_sobjs(
        sobjs_loose, esi_data['slit_left'], esi_data['slit_right'], esi_data['order_vec'],
        plate_scale, fof_link=fof_link
    )
    uniq_ids, cnts = np.unique(obj_id, return_counts=True)
    sobjs_loose = sobjs_loose[obj_id == uniq_ids[np.argmax(cnts)]]

    # Both runs must detect one object per order
    assert len(sobjs_tight) == esi_data['norders'], (
        'Mismatch between number of detected objects and orders when no trace shift allowed'
    )
    assert len(sobjs_loose) == esi_data['norders'], (
        'Mismatch between number of detected objects and orders when permissive shift allowed'
    )

    # The traces must differ on every order: a tighter maxshift prevents the
    # centroiding from correcting the trace position, so TRACE_SPAT diverges
    # from the fully refined (loose) result by at least 0.1 pixels on average.
    for sobj_tight, sobj_loose in zip(sobjs_tight, sobjs_loose):
        mean_diff = np.abs(sobj_loose.TRACE_SPAT - sobj_tight.TRACE_SPAT).mean()
        assert mean_diff > 0.1, (
            f"order {sobj_tight.ECH_ORDER}: expected TRACE_SPAT to differ by "
            f">0.1 pix between maxshift=0 and maxshift=2, got {mean_diff:.4f}")


