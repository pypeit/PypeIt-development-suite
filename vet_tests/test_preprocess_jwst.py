"""
Vet the products of pypeit_preprocess_jwst.

The pypeit_preprocess_jwst dev-suite test (see PypeItPreprocessJWSTTest in
test_scripts/pypeit_tests.py) runs the JWST pipeline (Stage 1 and Stage 2) on
the first NRS1 exposure of the JWST/NIRSpec fixed-slit S200A2 data, with
``--source_type POINT``.  The input files of the S200A2_G395M_F290LP setup were
produced in the same way, so they are used as reference.
"""
from pathlib import Path
import importlib.util
import json
import os

import numpy as np
import pytest

from astropy.io import fits

from pypeit.pypeitsetup import PypeItSetup
from pypeit.scripts.preprocess_jwst import PreprocessJWST
from pypeit.spectrographs import jwst_nirspec

# The preprocessing test only runs if the jwst package is installed
pytestmark = pytest.mark.skipif(importlib.util.find_spec('jwst') is None,
                                reason='jwst package not installed; pypeit_preprocess_jwst '
                                       'was not run')

SETUP = 'S200A2_G395M_F290LP_preprocess'
REF_SETUP = 'S200A2_G395M_F290LP'
ROOT = 'jw01967012001_03102_00001_nrs1'
FS_SLITS = ['S1600A1', 'S200A1', 'S200A2', 'S400A1']


def _output_dir(redux_out):
    out = Path(redux_out).absolute() / 'jwst_nirspec' / SETUP / 'jwst_output'
    assert out.is_dir(), f'Missing pypeit_preprocess_jwst output directory: {out}'
    return out


def _slit_extensions(hdul, extname='SCI'):
    # Extensions of each slit, selected by the EXTVER of the SCI extension of the
    # slit (only the SCI extensions have the SLTNAME keyword)
    return {h.header['SLTNAME']: hdul[extname, h.ver] for h in hdul if h.name == 'SCI'}


def test_preprocess_jwst_products(redux_out):
    out = _output_dir(redux_out)

    # Stage 1 and Stage 2 products used by PypeIt
    for suffix in ['rate', 'rateints', 'assign_wcs', 'extract_2d', 'cal', 'interpolatedflat']:
        assert (out / f'{ROOT}_{suffix}.fits').is_file(), f'Missing {suffix} file'
    # resample_spec and extract_1d are skipped
    assert not any(out.glob('*_s2d.fits')) and not any(out.glob('*_x1d.fits')), \
        'Resampled or 1D-extracted products should not be produced'

    # The parameters used are saved and are the defaults (with source type POINT)
    with open(out / 'det1_config.json') as f:
        assert json.load(f) == PreprocessJWST.default_det1_params(), \
            'Wrong Detector1Pipeline parameters'
    with open(out / 'spec2_config.json') as f:
        assert json.load(f) == PreprocessJWST.default_spec2_params(source_type='POINT'), \
            'Wrong Spec2Pipeline parameters'

    # The parameters reached the JWST pipeline: check the step status keywords
    rate_hdr = fits.getheader(out / f'{ROOT}_rate.fits')
    assert rate_hdr['S_CLNFNS'] == 'COMPLETE', 'clean_flicker_noise was not applied'
    assert rate_hdr['S_RAMP'] == 'COMPLETE', 'ramp_fit was not applied'
    cal_hdr = fits.getheader(out / f'{ROOT}_cal.fits')
    for key, status in [('S_WCS', 'COMPLETE'), ('S_EXTR2D', 'COMPLETE'), ('S_FLAT', 'COMPLETE'),
                        ('S_PHOTOM', 'COMPLETE'), ('S_BKDSUB', 'SKIPPED'),
                        ('S_NSCLEN', 'SKIPPED')]:
        assert cal_hdr[key] == status, f'Wrong status of the {key} step: {cal_hdr[key]}'

    # All the fixed slits are extracted, and the target slit (S200A2) is
    # calibrated as a point source
    with fits.open(out / f'{ROOT}_cal.fits') as hdul:
        assert hdul[0].header['FXD_SLIT'] == 'S200A2', 'Wrong target slit'
        sci = _slit_extensions(hdul)
        assert sorted(sci) == FS_SLITS, 'Wrong extracted slits'
        assert sci['S200A2'].header['SRCTYPE'] == 'POINT', 'Wrong source type for S200A2'
        assert all(sci[s].header['SRCTYPE'] == 'EXTENDED' for s in FS_SLITS if s != 'S200A2'), \
            'Wrong source type for the slits without the target'
        # The slit GWCS is stored in the file (used by PypeIt for the slit
        # platescale and the dither offsets)
        assert 'ASDF' in hdul, 'Missing ASDF extension with the GWCS'


def test_preprocess_jwst_pypeit_setup(redux_out):
    # PypeIt can ingest the products
    out = _output_dir(redux_out)
    ps = PypeItSetup.from_file_root(out, 'jwst_nirspec')
    ps.run(setup_only=True)
    # The other products (e.g., rate) are listed, but not used (frame type None)
    frametypes = {Path(f.lstrip('# ')).name: str(t)
                  for f, t in zip(ps.fitstbl['filename'], ps.fitstbl['frametype'])}
    used = {f: t for f, t in frametypes.items() if t != 'None'}
    assert used == {f'{ROOT}_cal.fits': 'arc,tilt',
                    f'{ROOT}_interpolatedflat.fits': 'pixelflat,trace',
                    f'{ROOT}_assign_wcs.fits': 'science'}, f'Wrong frame types: {frametypes}'
    assert sorted(set(frametypes) - set(used)) == \
           [f'{ROOT}_{s}.fits' for s in ['extract_2d', 'rate', 'rateints']], \
        f'Unexpected files: {frametypes}'

    # The along-slit platescale of the target slit from the GWCS
    if not jwst_nirspec.stdatamodels_installed():
        pytest.skip('stdatamodels not installed; cannot read the slit GWCS')
    platescale = jwst_nirspec.slit_platescale_from_gwcs(out / f'{ROOT}_cal.fits', 'S200A2')
    assert np.isclose(platescale, 0.1053, atol=5e-4), f'Wrong S200A2 platescale: {platescale}'


def test_preprocess_jwst_vs_reference(redux_out):
    # The products must reproduce the input files of the S200A2_G395M_F290LP
    # setup, which were produced in the same way.  The comparison is only
    # meaningful if the same version of the JWST pipeline and the same CRDS
    # context are used.
    out = _output_dir(redux_out)
    ref_dir = Path(os.environ['PYPEIT_DEV']) / 'RAW_DATA' / 'jwst_nirspec' / REF_SETUP
    if not (ref_dir / f'{ROOT}_cal.fits').is_file():
        pytest.skip(f'Reference files not found in {ref_dir}')
    hdr, ref_hdr = fits.getheader(out / f'{ROOT}_cal.fits'), fits.getheader(ref_dir / f'{ROOT}_cal.fits')
    if (hdr['CAL_VER'], hdr['CRDS_CTX']) != (ref_hdr['CAL_VER'], ref_hdr['CRDS_CTX']):
        pytest.skip(f"JWST pipeline version / CRDS context ({hdr['CAL_VER']}, {hdr['CRDS_CTX']}) "
                    f"differ from those of the reference files ({ref_hdr['CAL_VER']}, "
                    f"{ref_hdr['CRDS_CTX']}); regenerate the S200A2_G395M_F290LP input files")

    for suffix in ['assign_wcs', 'cal', 'interpolatedflat']:
        with fits.open(out / f'{ROOT}_{suffix}.fits') as hdul, \
                fits.open(ref_dir / f'{ROOT}_{suffix}.fits') as ref_hdul:
            if suffix == 'assign_wcs':
                # Full-frame image
                pairs = {('SCI', 'full frame'): (hdul['SCI'], ref_hdul['SCI'])}
            else:
                pairs = {}
                for extname in ['SCI', 'WAVELENGTH']:
                    ext, ref_ext = _slit_extensions(hdul, extname), _slit_extensions(ref_hdul, extname)
                    assert sorted(ext) == sorted(ref_ext), f'Different slits in the {suffix} files'
                    pairs.update({(extname, slit): (ext[slit], ref_ext[slit]) for slit in ext})
            for (extname, slit), (h, ref_h) in pairs.items():
                data, ref_data = h.data, ref_h.data
                assert data.shape == ref_data.shape, f'Different shape of {suffix} {extname} ({slit})'
                gpm = np.isfinite(ref_data)
                assert np.array_equal(np.isfinite(data), gpm), \
                    f'Different NaN mask in {suffix} {extname} ({slit})'
                # Differences are normalized by the typical value of the data
                norm = np.maximum(np.abs(ref_data[gpm]), np.median(np.abs(ref_data[gpm])))
                rel_diff = np.max(np.abs(data[gpm] - ref_data[gpm]) / norm)
                assert rel_diff < 1e-3, \
                    f'{suffix} {extname} ({slit}) differs from the reference: {rel_diff:.2e}'
