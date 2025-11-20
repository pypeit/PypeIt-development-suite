# Copied from ../jwst_diff_coadd2d.py 10 Jan 2023

# This only performs the basis spec2 reduction using the JWST pipeline

import os
from pathlib import Path

from IPython import embed

# set environment variables
if os.getenv('CRDS_PATH') is None:
    os.environ['CRDS_PATH'] = str(Path.home().resolve() / 'crds_cache')
if os.getenv('CRDS_SERVER_URL') is None:
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'

## JWST imports
from jwst.pipeline import Spec2Pipeline

raw_dir = Path('/Users/westfall/Work/JWST/data/10Jan2023/2756_nirspec/raw').resolve()
detectors = ['nrs1', 'nrs2']
mode = 'MSA'

jwst_rdx_dir = raw_dir.parent / 'jwst_rdx'
if not jwst_rdx_dir.exists():
    jwst_rdx_dir.mkdir(parents=True)

exp_seq = [#'jw02756001001_02101',      # CRDS errors on these exposures
           'jw02756001001_03101',
           'jw02756001001_03103']

# Image B of source JD1 has source_name = '2756_202', which is index 32 in the
# list of slits in jw02756001001_03101_00001_nrs2_cal.fits

# "Bright" sources are: 2756_410044, 2756_410045, 2756_80065

overwrite_jwst_spec2 = False

exposures = []
for _exp_seq in exp_seq:
    for d in detectors:
        exposures += list(sorted(raw_dir.glob(f'{_exp_seq}*{d}_rate.fits')))

# TODO: Flat-field usage unclear.
# TODO: Skip pathloss? (or is that done by setting extended source type)
jwst_spec2_par = {
    'extract_2d': {'save_results': True},
    #'residual_fringe': {'skip': True},
    'bkg_subtract': {'skip': True},
    'imprint_subtract': {'save_results': True},
    'msa_flagging': {'save_results': True},
    'master_background_mos': {'skip': True},
    'srctype': {'source_type': 'EXTENDED'},
    'resample_spec': {'skip': True},
    'extract_1d': {'skip': True},
    'flat_field': {'save_interpolated_flat': True}
}

# NOTE: In the order they're produced
jwst_rdx_products = ['msa_flagging', 'extract_2d', 'interpolatedflat', 'cal']

# TODO: Parallelize this
for sci in exposures:
    ofiles = [jwst_rdx_dir / sci.name.replace('rate', product) 
                    for product in jwst_rdx_products]
    if not overwrite_jwst_spec2 and all([f.exists() for f in ofiles]):
        continue
    Spec2Pipeline.call(sci, save_results=True, output_dir=str(jwst_rdx_dir),
                        steps=jwst_spec2_par)

