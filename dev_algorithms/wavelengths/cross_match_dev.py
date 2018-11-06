

import numpy as np
# Read in a wavelength solution and compute
from pypeit import wavecalib

calibfile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_NIRES/NIRES/MF_keck_nires/MasterWaveCalib_A_01_aa.json'
wv_calib, par = wavecalib.load_wv_calib(calibfile)
steps= wv_calib.pop('steps')
par_dum = wv_calib.pop('par')

datafile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_NIRES/NIRES/MF_keck_nires/MasterWaveCalib_A_01_ac.json'
wv_calib_data, par = wavecalib.load_wv_calib(datafile)
steps= wv_calib_data.pop('steps')
par_dum = wv_calib_data.pop('par')
nslits = len(wv_calib_data)
# assignments
spec = np.zeros((wv_calib_data['0']['spec'].size, nslits))
nonlinear_counts=par['nonlinear_counts']
sigdetect = par['lowest_nsig']
detections = None

# assignments
lamps = par['lamps']
# def archive_reidentify(spec, wv_calib, lamps, detections = None,
# rms_threshold = 0.15, nonlinear_counts =par['nonlinear_counts'], sigdetect = par['lowest_nsig'])



nspec = spec.shape[1]
narxiv = len(wv_calib)
nspec_arxiv = wv_calib['0']['spec'].size
if nspec_arxiv != nspec:
    msgs.error('Different spectral binning is not supported yet but it will be soon')

# If the detections were not passed in find the lines in each spectrum
if detections is None:
    detections = {}
    for slit in range(nslits):
        tcent, ecent, cut_tcent, icut = wvutils.arc_lines_from_spec(spec[:, slit], min_nsig=sigdetect,nonlinear_counts=nonlinear_counts)
        detections[str(slit)] = [tcent[icut].copy(), ecent[icut].copy()]
else:
    if len(detections) != nslits:
        msgs.error('Detections must be a dictionary with nslit elements')

# For convenience pull out all the spectra from the wv_calib archive
spec_arxiv = np.zeros((nspec, narxiv))
for islit in range(narxiv):
    spec_arxiv[:,islit] = wv_calib[str(islit)]['spec']

# JFH Changed this to take the median which is more robust. Could even reject outliers
disp_arxiv = np.zeros(narxiv, dtype=float)
wvc_arxiv = np.zeros(narxiv, dtype=float)
for islit in range(narxiv):
    fitc = wv_calib[str(good_slits[islit])]['fitc']
    xfit = xrng / (self._npix - 1)
    fitfunc = self._all_final_fit[str(good_slits[islit])]['function']
    fmin, fmax = 0.0, 1.0
    wave_soln = utils.func_val(fitc, xfit, fitfunc, minv=fmin, maxv=fmax)
    wvc_good[islit] = wave_soln[self._npix // 2]
    disp_good[islit] = np.median(wave_soln - np.roll(wave_soln, 1))


