

import numpy as np
# Read in a wavelength solution and compute
from pypeit import wavecalib
from pypeit.core.wavecal import waveio
from pypeit.core.wavecal import wvutils
from pypeit import utils
from astropy import table
from pypeit import msgs


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
for slit in range(nslits):
    spec[:,slit] = wv_calib_data[str(slit)]['spec']

nonlinear_counts=par['nonlinear_counts']
sigdetect = par['lowest_nsig']
detections = None

# assignments
lamps = par['lamps']
use_unknowns=True
debug = True
# def archive_reidentify(spec, wv_calib, lamps, detections = None,
# rms_threshold = 0.15, nonlinear_counts =par['nonlinear_counts'], sigdetect = par['lowest_nsig'], use_unknowns=True, debug=True)

#

# Generate the line list
line_lists = waveio.load_line_lists(lamps)
unknwns = waveio.load_unknown_list(lamps)
if use_unknowns:
    tot_list = table.vstack([line_lists, unknwns])
else:
    tot_list = line_lists
# Generate the final linelist and sort
wvdata = np.array(tot_list['wave'].data)  # Removes mask if any
wvdata.sort()

nspec = spec.shape[0]
narxiv = len(wv_calib)
nspec_arxiv = wv_calib['0']['spec'].size
if nspec_arxiv != nspec:
    msgs.error('Different spectral binning is not supported yet but it will be soon')

# If the detections were not passed in find the lines in each spectrum
if detections is None:
    detections = {}
    for islit in range(nslits):
        tcent, ecent, cut_tcent, icut = wvutils.arc_lines_from_spec(spec[:, islit], min_nsig=sigdetect,nonlinear_counts=nonlinear_counts)
        detections[str(islit)] = [tcent[icut].copy(), ecent[icut].copy()]
else:
    if len(detections) != nslits:
        msgs.error('Detections must be a dictionary with nslit elements')

# For convenience pull out all the spectra from the wv_calib archive
spec_arxiv = np.zeros((nspec, narxiv))
wave_soln_arxiv = np.zeros((nspec, narxiv))
wvc_arxiv = np.zeros(narxiv, dtype=float)
disp_arxiv = np.zeros(narxiv, dtype=float)
xrng = np.arange(nspec_arxiv)
for iarxiv in range(narxiv):
    spec_arxiv[:,slit] = wv_calib[str(iarxiv)]['spec']
    fitc = wv_calib[str(iarxiv)]['fitc']
    xfit = xrng/(nspec_arxiv - 1)
    fitfunc = wv_calib[str(iarxiv)]['function']
    fmin, fmax = wv_calib[str(iarxiv)]['fmin'],wv_calib[str(iarxiv)]['fmax']
    wave_soln_arxiv[:,iarxiv] = utils.func_val(fitc, xfit, fitfunc, minv=fmin, maxv=fmax)
    wvc_arxiv[iarxiv] = wave_soln_arxiv[nspec_arxiv//2, slit]
    disp_arxiv[iarxiv] = np.median(wave_soln_arxiv[:,iarxiv] - np.roll(wave_soln_arxiv[:,iarxiv], 1))


# Loop over the slits in the spectrum and cross-correlate each with each arxiv spectrum to identify lines
for islit in range(nslits):
    slit_det = detections[str(islit)][0]
    lindex = np.array([], dtype=np.int)
    dindex = np.array([], dtype=np.int)
    wcen = np.zeros(narxiv)
    disp = np.zeros(narxiv)
    shift_vec = np.zeros(narxiv)
    stretch_vec = np.zeros(narxiv)
    ccorr_vec = np.zeros(narxiv)
    for iarxiv in range(narxiv):
        msgs.info('Cross-correlating slit # {:d}'.format(islit + 1) + ' with arxiv slit # {:d}'.format(iarxiv + 1))
        # Match the peaks between the two spectra.
        success, shift_vec[cntr], stretch_vec[cntr], ccorr_vec[cntr], _, _ = \
            wvutils.xcorr_shift_stretch(spec[:, islit], spec_arxiv[:, iarxiv], debug=debug)