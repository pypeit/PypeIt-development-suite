

import numpy as np
# Read in a wavelength solution and compute
from pypeit import wavecalib
from pypeit.core.wavecal import autoid



instrument = 'NIRES'
if instrument == 'NIRES':
    calibfile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_NIRES/NIRES/MF_keck_nires/MasterWaveCalib_A_01_aa.json'
    wv_calib_arxiv, par = wavecalib.load_wv_calib(calibfile)
    steps= wv_calib_arxiv.pop('steps')
    par_dum = wv_calib_arxiv.pop('par')

    datafile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_NIRES/NIRES/MF_keck_nires/MasterWaveCalib_A_01_ac.json'
    wv_calib_data, par = wavecalib.load_wv_calib(datafile)
    steps= wv_calib_data.pop('steps')
    par_dum = wv_calib_data.pop('par')
elif instrument == 'LRIS-R':
    # Use one detector as the arxiv the other as the data
    calibfile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_LRIS_red/multi_400_8500_d560/MF_keck_lris_red/MasterWaveCalib_A_01_aa.json'
    wv_calib_arxiv, par = wavecalib.load_wv_calib(calibfile)
    steps= wv_calib_arxiv.pop('steps')
    par_dum = wv_calib_arxiv.pop('par')

    datafile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_LRIS_red/multi_400_8500_d560/MF_keck_lris_red/MasterWaveCalib_A_02_aa.json'
    wv_calib_data, par = wavecalib.load_wv_calib(datafile)
    steps= wv_calib_data.pop('steps')
    par_dum = wv_calib_data.pop('par')
elif instrument == 'LRIS-B':
    # Use one detector as the arxiv the other as the data
    calibfile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_LRIS_blue/multi_600_4000_d560/MF_keck_lris_blue/MasterWaveCalib_A_02_aa.json'
    wv_calib_arxiv, par = wavecalib.load_wv_calib(calibfile)
    steps= wv_calib_arxiv.pop('steps')
    par_dum = wv_calib_arxiv.pop('par')

    datafile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_LRIS_blue/multi_600_4000_d560/MF_keck_lris_blue/MasterWaveCalib_A_01_aa.json'
    wv_calib_data, par = wavecalib.load_wv_calib(datafile)
    steps= wv_calib_data.pop('steps')
    par_dum = wv_calib_data.pop('par')

nslits = len(wv_calib_data)
# assignments
spec = np.zeros((wv_calib_data['0']['spec'].size, nslits))
for slit in range(nslits):
    spec[:,slit] = wv_calib_data[str(slit)]['spec']

match_toler = par['match_toler']
n_first = par['n_first']
sigrej_first = par['sigrej_first']
n_final = par['n_final']
sigrej_final = par['sigrej_final']
func = par['func']
nonlinear_counts=par['nonlinear_counts']
sigdetect = par['lowest_nsig']
rms_threshold = par['rms_threshold']
lamps = par['lamps']

nreid_min = 1
wv_calib_out, patt_dict, bad_slits = autoid.reidentify(spec, wv_calib_arxiv, lamps, nreid_min, rms_threshold=rms_threshold,
                                                nonlinear_counts=nonlinear_counts,sigdetect=sigdetect,use_unknowns=True,
                                                match_toler=match_toler,func='legendre',n_first=n_first,
                                                sigrej_first=sigrej_first,n_final=n_final, sigrej_final=sigrej_final,
                                                debug_xcorr=True, debug_reid=True)