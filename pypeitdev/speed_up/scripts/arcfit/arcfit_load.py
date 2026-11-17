import glob
import numpy as np
from pypeit.core import skyspec
from pypeit.wavecalib import WaveCalib
def load(which):
    if which == 'sky': return [skyspec.load_sky_spectrum('paranal_sky.fits').flux]
    out = []
    for f in sorted(glob.glob('/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT/keck_deimos/600ZD_M_6500/Calibrations/WaveCalib_A_0_MSC0*.fits')):
        wc = WaveCalib.from_file(f, chk_version=False); out += [wc.arc_spectra[:, i] for i in range(wc.arc_spectra.shape[1])]
    return out
