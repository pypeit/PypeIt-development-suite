from pathlib import Path

import numpy as np

from pypeit.slittrace import SlitTraceSet
from pypeit.flatfield import FlatImages
from pypeit import spec2dobj
from pypeit import specobjs

def compare_float_arrays(array1, array2, min_tolerance=16):
    # None values cause issues, to fix this we replace them with nans
    if np.any(array1==None):       
        array1[array1==None] = np.nan

    if np.any(array2==None):       
        array2[array2==None] = np.nan

    rtols = [1* (10**-i) for i in range(min_tolerance, -1, -1)]
    for rtol in rtols:
        if np.allclose(array1, array2, rtol=rtol, equal_nan=True):
            return rtol
        
    return np.inf

def compare_calibrations(full_calib_path, by_steps_calib_path):
    # Compare slits between a full reduction and one done step by step

    good_slits = 0
    good_flats = 0
    for by_steps_slits_file in by_steps_calib_path.glob("Slits_*.fits.gz"):
        # Look for its counter part in the full reduction path
        full_slit_file = full_calib_path / by_steps_slits_file.name
        assert full_slit_file.exists(), f"Slits file {by_steps_slits_file} does not exist in the full reduction path: {full_calib_path}"

        full_slits = SlitTraceSet.from_file(full_slit_file)
        by_steps_slits = SlitTraceSet.from_file(by_steps_slits_file)
        assert full_slits.center is not None
        assert by_steps_slits.center is not None
        tol= compare_float_arrays(full_slits.center, by_steps_slits.center)
        assert tol <= 1e-8, f"Slits file {by_steps_slits_file} does not match the one in the full reduction path {full_calib_path}. Tolerance: {tol}"
        good_slits += 1

    # Compare flats between a full reduction and one done step by step

    for by_steps_flat_file in by_steps_calib_path.glob("Flat_*.fits"):
        # Look for its counter part in the full reduction path
        full_flat_file = full_calib_path / by_steps_flat_file.name
        assert full_flat_file.exists(), f"Flat file {by_steps_flat_file} does not exist in the full reduction path: {full_calib_path}"

        full_flats = FlatImages.from_file(full_flat_file)
        by_steps_flats = FlatImages.from_file(by_steps_flat_file)

        tol = compare_float_arrays(full_flats.pixelflat_norm, by_steps_flats.pixelflat_norm) 
        assert tol <=1e-8, f"Flats file {by_steps_flat_file} does not match the one in the full reduction path {full_calib_path}. Tolerance: {tol}"
        good_flats += 1

    assert good_slits >= 1, f"No slits file found in {by_steps_calib_path}"
    assert good_flats >= 1, f"No flats file found in {by_steps_calib_path}"

def compare_science_frames(raw_file, det, full_science_path, by_steps_science_path):

    # Check the 'sciimg' field in the spec2ds to see if the by_steps processing agrees with the full reduction
    full_spec2d_paths = list(full_science_path.glob(f"spec2d_{raw_file}*"))
    assert len(full_spec2d_paths) > 0, f"No spec2d file found for {raw_file} in {full_science_path}"    
    full_spec2d = spec2dobj.Spec2DObj.from_file(full_spec2d_paths[0], detname=det)

    by_steps_spec2d_paths = list(by_steps_science_path.glob(f"spec2d_{raw_file}*"))
    assert len(by_steps_spec2d_paths) > 0, f"No spec2d file found for {raw_file} in {by_steps_science_path}"
    by_steps_spec2d = spec2dobj.Spec2DObj.from_file(by_steps_spec2d_paths[0], detname=det)

    tol = compare_float_arrays(full_spec2d.sciimg, by_steps_spec2d.sciimg)
    print(f"{det} Spec2D sciimg is within {tol}")
    assert tol < 1, f"Spec2D file {by_steps_spec2d_paths[0]} does not match file in {full_science_path}, tolerance: {tol}"

    # Now check the spec1d file
    full_spec1d_paths = list(full_science_path.glob(f"spec1d_{raw_file}*.fits"))
    assert len(full_spec1d_paths) > 0, f"No spec1d file found for {raw_file} in {full_science_path}"    
    full_spec1d = specobjs.SpecObjs.from_fitsfile(full_spec1d_paths[0])

    by_steps_spec1d_paths = list(by_steps_science_path.glob(f"spec1d_{raw_file}*.fits"))
    assert len(by_steps_spec1d_paths) > 0, f"No spec1d file found for {raw_file} in {by_steps_science_path}"
    by_steps_spec1d = specobjs.SpecObjs.from_fitsfile(by_steps_spec1d_paths[0])


    det_idx = full_spec1d.DET == det
    full_spec1d_det = full_spec1d[det_idx]

    det_idx = by_steps_spec1d.DET == det
    by_steps_spec1d_det = by_steps_spec1d[det_idx]

    assert len(full_spec1d_det) == len(by_steps_spec1d_det), f"{full_spec1d_paths[0]} has {len(full_spec1d_det)} objects for detector {det} but {by_steps_spec1d_paths[0]} has {len(by_steps_spec1d_det)}"

    for full_sobj in full_spec1d_det:
        idx = by_steps_spec1d_det.NAME == full_sobj.NAME
        if np.sum(idx) == 0:
            print(f"{det} {full_sobj.NAME} is not in {by_steps_spec1d_paths[0]}")
            continue
        if np.sum(idx) > 1:
            print(f"{det} {full_sobj.NAME} is duplicated in {by_steps_spec1d_paths[0]}")
            continue

        by_steps_sobj = by_steps_spec1d_det[idx][0]

        array_fields = ["BOX_COUNTS", "BOX_WAVE", "OPT_COUNTS", "OPT_WAVE"]

        for field in array_fields:
            tol = compare_float_arrays(full_sobj[field], by_steps_sobj[field])
            print(f"{det} {full_sobj.NAME}.{field} is within {tol}")
            assert tol < 1, f"Object {full_sobj.NAME}.{field} does not match between {full_spec1d_paths[0]} and {by_steps_spec1d_paths[0]}. Tolerance: {tol}"




def test_shane_kast_blue(redux_out):
    base_dir = Path(redux_out, 'shane_kast_blue', '600_4310_d55','shane_kast_blue_A')

    full_science_dir = Path(base_dir, "Science")
    by_steps_science_dir = Path(base_dir, "step_by_step", "Science")

    compare_science_frames("b24", "DET01", full_science_dir, by_steps_science_dir)
    compare_science_frames("b27", "DET01", full_science_dir, by_steps_science_dir)

def test_keck_deimos(redux_out):
    base_dir = Path(redux_out, 'keck_deimos', '600ZD_tilted')
    
    full_science_dir = Path(base_dir, "Science")
    by_steps_science_dir = Path(base_dir, "step_by_step", "Science")

    compare_science_frames("d0225_0054", "MSC01", full_science_dir, by_steps_science_dir)

def test_keck_mosfire(redux_out):
    base_dir = Path(redux_out, 'keck_mosfire', 'Y_long')

    full_science_dir = Path(base_dir, "Science")
    by_steps_science_dir = Path(base_dir, "step_by_step", "Science")

    compare_science_frames("m191118_0064", "DET01", full_science_dir, by_steps_science_dir)
    compare_science_frames("m191120_0043", "DET01", full_science_dir, by_steps_science_dir)
