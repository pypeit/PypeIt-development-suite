

import numpy as np
# Read in a wavelength solution and compute
from pypeit import wavecalib

calibfile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_NIRES/NIRES/MF_keck_nires/MasterWaveCalib_A_01_aa.json'
wv_calib, par = wavecalib.load_wv_calib(calibfile)
steps= wv_calib.pop('steps')
par_dum = wv_calib.pop('par')
nslits = len(wv_calib)
nspec = wv_calib['0']['spec'].size

# def archive_reidentify(spec_new, spec_detns, spec_archve, wv_calib, wvdata, rms_threshold = 0.15)

# assignments
spec = np.zeros((nspec, nslits))
for islit in range(nslits):
    spec[:,islit] = wv_calib[str(islit)]['spec']