

# Read in a wavelength solution and compute
from pypeit import wavecalib

calibfile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_NIRES/NIRES/MF_keck_nires/MasterWaveCalib_A_01_aa.json'
wv_calib, par = wavecalib.load_wv_calib(calibfile)

