"""
Tests for the Binospec IFU illumination correction and cube building tools.
"""
from pathlib import Path

import numpy as np
from astropy.io import fits

from pypeit.scripts.binospec_ifu_illumcorr import BinospecIFUIllumCorr
from pypeit.scripts.binospec_ifu_cube import BinospecIFUCube


def test_binospec_ifu_illumcorr(redux_out):
    """Test that pypeit_binospec_ifu_illumcorr produces a corrected spec1d file."""
    sci_dir = Path(redux_out) / 'mmt_binospec_ifu' / 'G270' / 'Science'
    spec1d_files = sorted(sci_dir.glob('spec1d_*.fits'))
    assert len(spec1d_files) > 0, f"No spec1d files found in {sci_dir}"

    spec1d_file = spec1d_files[0]

    # Run illumination correction
    pargs = BinospecIFUIllumCorr.parse_args([str(spec1d_file)])
    BinospecIFUIllumCorr.main(pargs)

    # Check the output file was created
    illumcorr_file = spec1d_file.with_name(
        spec1d_file.stem + '_illumcorr' + spec1d_file.suffix)
    assert illumcorr_file.exists(), \
        f"Illumination-corrected file not created: {illumcorr_file}"

    # Check ILLUMCOR header flag is set
    with fits.open(illumcorr_file) as hdu:
        assert hdu[0].header.get('ILLUMCOR', False), \
            "ILLUMCOR header flag not set in output file"

        # Check that at least some BinTableHDU extensions exist (fiber spectra)
        n_tables = sum(1 for ext in hdu[1:] if isinstance(ext, fits.BinTableHDU))
        assert n_tables > 0, "No BinTableHDU extensions found in output"

    # Clean up
    illumcorr_file.unlink()


def test_binospec_ifu_cube_from_spec1d(redux_out):
    """Test building a datacube from a Binospec IFU spec1d file."""
    sci_dir = Path(redux_out) / 'mmt_binospec_ifu' / 'G270' / 'Science'
    spec1d_files = sorted(sci_dir.glob('spec1d_*.fits'))
    assert len(spec1d_files) > 0, f"No spec1d files found in {sci_dir}"

    spec1d_file = spec1d_files[0]
    output_file = sci_dir / 'test_cube_spec1d.fits'

    # Build cube from spec1d
    pargs = BinospecIFUCube.parse_args([str(spec1d_file), '-o', str(output_file)])
    BinospecIFUCube.main(pargs)

    # Validate the output datacube
    assert output_file.exists(), f"Datacube not created: {output_file}"

    with fits.open(output_file) as hdu:
        # Check expected HDU structure
        assert 'FLUX' in hdu, "FLUX HDU missing from datacube"
        assert 'VAR' in hdu, "VAR HDU missing from datacube"

        flux = hdu['FLUX'].data
        var = hdu['VAR'].data

        # Check 3D shape (n_wave, ny, nx)
        assert flux.ndim == 3, f"FLUX is not 3D: {flux.ndim}D"
        assert var.shape == flux.shape, "VAR shape does not match FLUX shape"

        # Check that the cube contains real data
        assert np.any(flux != 0), "FLUX cube is all zeros"
        assert np.any(var > 0), "VAR cube has no positive values"

        # Check WCS keywords
        hdr = hdu['FLUX'].header
        assert hdr['CTYPE1'] == 'RA---TAN'
        assert hdr['CTYPE2'] == 'DEC--TAN'
        assert hdr['CTYPE3'] == 'WAVE'

    # Clean up
    output_file.unlink()


def test_binospec_ifu_cube_from_spec2d(redux_out):
    """Test building a datacube from a Binospec IFU spec2d file."""
    sci_dir = Path(redux_out) / 'mmt_binospec_ifu' / 'G270' / 'Science'
    spec2d_files = sorted(sci_dir.glob('spec2d_*.fits'))
    assert len(spec2d_files) > 0, f"No spec2d files found in {sci_dir}"

    spec2d_file = spec2d_files[0]
    output_file = sci_dir / 'test_cube_spec2d.fits'

    # Build cube from spec2d
    pargs = BinospecIFUCube.parse_args([str(spec2d_file), '-o', str(output_file)])
    BinospecIFUCube.main(pargs)

    # Validate the output datacube
    assert output_file.exists(), f"Datacube not created: {output_file}"

    with fits.open(output_file) as hdu:
        assert 'FLUX' in hdu, "FLUX HDU missing from datacube"
        assert 'VAR' in hdu, "VAR HDU missing from datacube"

        flux = hdu['FLUX'].data
        var = hdu['VAR'].data

        assert flux.ndim == 3, f"FLUX is not 3D: {flux.ndim}D"
        assert var.shape == flux.shape, "VAR shape does not match FLUX shape"

        assert np.any(flux != 0), "FLUX cube is all zeros"
        assert np.any(var > 0), "VAR cube has no positive values"

        # Check WCS keywords
        hdr = hdu['FLUX'].header
        assert hdr['CTYPE1'] == 'RA---TAN'
        assert hdr['CTYPE2'] == 'DEC--TAN'
        assert hdr['CTYPE3'] == 'WAVE'

    # Clean up
    output_file.unlink()
