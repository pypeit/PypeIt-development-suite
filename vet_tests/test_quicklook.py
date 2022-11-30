"""
Module to run tests on scripts
"""
import os
import glob
import numpy as np
import pytest

from pypeit.pypmsgs import PypeItError
from pypeit.inputfiles import PypeItFile
from pypeit import specobjs

def test_shane_kast_ql(redux_out):
    instr = 'shane_kast_blue' 
    outroot = os.path.join(redux_out,
                             instr,
                             '600_4310_d55')

    for test in ['std', 'masters', 'calibs', 'multi', 'nostack', 'boxcar']:
        outdir = os.path.join(outroot, f'QL_{test}')
        if test in ['std', 'masters', 'calibs', 'boxcar']:
            rdxfolder = 'b27'
        else:
            rdxfolder = 'b27-b28'
        rdxdir = os.path.join(outdir, rdxfolder)
        scidir = os.path.join(rdxdir, 'Science')

        # PypeIt File
        pypeit_file = os.path.join(rdxdir, f'{instr}_A.pypeit')
        assert os.path.isfile(pypeit_file)
        pypeitFile = PypeItFile.from_file(pypeit_file)
        assert pypeitFile.config['rdx']['quicklook']

        # Outputs
        spec2d_files = glob.glob(os.path.join(scidir, 'spec2d*')) 
        nfiles = 2 if test in ['nostack'] else 1
        assert len(spec2d_files) == nfiles
        spec1d_files = glob.glob(os.path.join(scidir, 'spec1d*.fits')) 
        assert len(spec1d_files) == nfiles

        sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_files[0])
        assert sobjs.nobj == 1

        # Additional
        if test == 'boxcar':
            assert np.isclose(sobjs.BOX_RADIUS[0], 4.651162790697675)
        else:
            assert not np.isclose(sobjs.BOX_RADIUS[0], 4.651162790697675)

    # Check boxcar -- 4.65116279