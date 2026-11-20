"""
Write a NOTES.txt file in each RAW_DATA/jwst_nirspec/<setup> folder (and in the
CRDS cache of the preprocessing test) explaining how the JWST products were
produced with pypeit_preprocess_jwst.

Usage: python make_notes.py [jwst_nirspec_dir]

jwst_nirspec_dir defaults to $PYPEIT_DEV/RAW_DATA/jwst_nirspec.  Update DATE,
VERSIONS, CONTEXT and REQUIREMENTS when the products are regenerated for a new
version of the JWST pipeline.
"""
import argparse
import os
from pathlib import Path

SCRIPT = 'PypeIt-development-suite/pypeitdev/jwst/preprocessing/make_notes.py'

DATE = '2026-11-19'
CONTEXT = 'jwst_1584.pmap'
VERSIONS = 'jwst 3.0.0, stdatamodels 6.0.0, crds 14.0.6, PypeIt 3.0.1.dev'
REQUIREMENTS = 'jwst>=3.0.0,<3.1.0 and stdatamodels>=6.0.0,<6.1.0'

UPDATE = f"""\
PypeIt requires {REQUIREMENTS}.  If these requirements
change, regenerate these files (and all the other jwst_nirspec setups) with the
new versions and their default CRDS context, and update this file with
{SCRIPT}."""

SETUPS = {
    'S200A2_G395M_F290LP': (
        'jw01967012001_03102_0000[1-3]_nrs[12]_uncal.fits',
        '--source_type POINT',
        'Also the reference of vet_tests/test_preprocess_jwst.py: regenerate it '
        'together with\nS200A2_G395M_F290LP_preprocess.'),
    'MSA_PRISM_nomosaic': (
        'jw01287001001_1[13]101_0000[1-3]_nrs[12]_uncal.fits, jw01287001001_02_msa.fits',
        '', ''),
    'MSA_PRISM_2visits': (
        'jw0128700100[13]_11101_0000[1-3]_nrs1_uncal.fits, jw0128700100[13]_02_msa.fits',
        '', ''),
    'MSA_G395M_F290LP_mosaic': (
        'jw03543001001_07101_0000[2-4]_nrs[12]_uncal.fits, jw03543001001_03_msa.fits',
        '', ''),
}


def setup_notes(raw, options, extra):
    cmd = ' '.join(['pypeit_preprocess_jwst raw -o jwst_output', options,
                    f'--crds_context {CONTEXT}'])
    cmd = ' '.join(cmd.split())
    text = f"""\
Produced with pypeit_preprocess_jwst on {DATE}
  Software:      {VERSIONS}
  CRDS context:  {CONTEXT} (default of jwst 3.0)
  Raw data:      {raw}
  Command:       {cmd}
Only the cal, interpolatedflat and assign_wcs products are kept.
"""
    if extra:
        text += f'{extra}\n'
    return text + f'\n{UPDATE}\n'


PREPROCESS = f"""\
Input of the pypeit_preprocess_jwst dev-suite test ({DATE}; {VERSIONS}).
  jw01967012001_03102_00001_nrs1_uncal.fits  raw exposure (first NRS1 exposure
                                             of S200A2_G395M_F290LP)
  crds/  CRDS cache for this exposure with context {CONTEXT} (default of
         jwst 3.0); the test runs offline with --crds_path crds
         --crds_context {CONTEXT} --crds_offline (see crds/NOTES.txt)

PypeIt requires {REQUIREMENTS}.  If these requirements
change, the uncal file stays the same; only crds/ must be replaced with the
files of the new default CRDS context (see crds/NOTES.txt).  Then also update
the crds_context of the test in test_scripts/test_setups.py, regenerate
S200A2_G395M_F290LP (the reference of the vet test), and update this file with
{SCRIPT}.
"""

CRDS = f"""\
Minimal CRDS cache with the default CRDS context of jwst 3.0 ({CONTEXT}),
for jw01967012001_03102_00001_nrs1_uncal.fits ({DATE}; {VERSIONS}):
config/jwst, all the mappings of {CONTEXT}, and the reference files that
the pipeline reports in its log ("Prefetch for ..." and "PARS-... parameters
found" lines) when run on this exposure against a complete CRDS cache.

If PypeIt's jwst requirements ({REQUIREMENTS}) change,
download the files of the new default context in the same way and check them by
running the test with --crds_offline.  This file is written by {SCRIPT}.
"""


def main():
    parser = argparse.ArgumentParser(description='Write the NOTES.txt files of the JWST/NIRSpec '
                                                 'dev-suite data.')
    parser.add_argument('jwst_nirspec_dir', nargs='?', type=Path,
                        default=Path(os.environ.get('PYPEIT_DEV', '.')) / 'RAW_DATA' / 'jwst_nirspec',
                        help='Directory with the jwst_nirspec setups (default: '
                             '$PYPEIT_DEV/RAW_DATA/jwst_nirspec)')
    out = parser.parse_args().jwst_nirspec_dir
    for setup, args in SETUPS.items():
        (out / setup / 'NOTES.txt').write_text(setup_notes(*args))
    (out / 'S200A2_G395M_F290LP_preprocess' / 'NOTES.txt').write_text(PREPROCESS)
    (out / 'S200A2_G395M_F290LP_preprocess' / 'crds' / 'NOTES.txt').write_text(CRDS)
    print(f'NOTES.txt files written in {out}')


if __name__ == '__main__':
    main()
