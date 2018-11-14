

import numpy as np
# Read in a wavelength solution and compute
from pypeit import wavecalib
from pypeit.core.wavecal import autoid, waveio
from pypeit import utils


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

narxiv = len(wv_calib_arxiv)
nspec = wv_calib_arxiv['0']['spec'].size
# assignments
spec_arxiv = np.zeros((nspec, narxiv))
for iarxiv in range(narxiv):
    spec_arxiv[:,iarxiv] = wv_calib_arxiv[str(iarxiv)]['spec']

det_arxiv = {}
wave_soln_arxiv = np.zeros((nspec, narxiv))
xrng = np.arange(nspec)
for iarxiv in range(narxiv):
    spec_arxiv[:, iarxiv] = wv_calib_arxiv[str(iarxiv)]['spec']
    fitc = wv_calib_arxiv[str(iarxiv)]['fitc']
    fitfunc = wv_calib_arxiv[str(iarxiv)]['function']
    fmin, fmax = wv_calib_arxiv[str(iarxiv)]['fmin'], wv_calib_arxiv[str(iarxiv)]['fmax']
    wave_soln_arxiv[:, iarxiv] = utils.func_val(fitc, xrng, fitfunc, minv=fmin, maxv=fmax)
    det_arxiv[str(iarxiv)] = wv_calib_arxiv[str(iarxiv)]['xfit']



match_toler = 2.0 #par['match_toler']
n_first = par['n_first']
sigrej_first = par['sigrej_first']
n_final = par['n_final']
sigrej_final = par['sigrej_final']
func = par['func']
nonlinear_counts=par['nonlinear_counts']
sigdetect = par['lowest_nsig']
rms_threshold = par['rms_threshold']
lamps = par['lamps']

line_list = waveio.load_line_lists(lamps)


cc_thresh =0.8
cc_local_thresh = 0.8
n_local_cc =11


nreid_min = 1

new = False
if new:
    all_patt_dict={}
    all_detections = {}
    for islit in range(nslits):
        all_detections[str(islit)], all_patt_dict[str(islit)] = autoid.reidentify(spec[:,islit], spec_arxiv, wave_soln_arxiv, det_arxiv, line_list, nreid_min,
                                                                                  detections=None, cc_thresh=cc_thresh,cc_local_thresh=cc_local_thresh,
                                                                                  match_toler=match_toler, nlocal_cc=11, nonlinear_counts=nonlinear_counts,sigdetect=sigdetect,
                                                                                  debug_xcorr=True, debug_reid=True)
else:
    wv_calib_out, patt_dict, bad_slits = autoid.reidentify_old(spec, wv_calib_arxiv, lamps, nreid_min, rms_threshold=rms_threshold,
                                                               nonlinear_counts=nonlinear_counts,sigdetect=sigdetect,use_unknowns=True,
                                                               match_toler=match_toler,func='legendre',n_first=n_first,
                                                               sigrej_first=sigrej_first,n_final=n_final, sigrej_final=sigrej_final,
                                                               debug_xcorr=True, debug_reid=True)