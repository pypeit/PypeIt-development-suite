

import numpy as np
# Read in a wavelength solution and compute
from pypeit import wavecalib
from pypeit.core.wavecal import waveio, wvutils, fitting, patterns, qa
from pypeit import utils
from astropy import table
from pypeit import msgs
from matplotlib import pyplot as plt
import copy


def fit_slit(spec, patt_dict, tcent, line_lists, outroot=None, slittxt="Slit", thar=False,match_toler=3.0,
             func='legendre', n_first=2,sigrej_first=2.0,n_final=4,sigrej_final=3.0,verbose=False):
    """ Perform a fit to the wavelength solution

    Parameters
    ----------
    spec : ndarray
      arc spectrum
    patt_dict : dict
      dictionary of patterns
    tcent: ndarray
      List of the detections in this slit to be fit using the patt_dict
    outroot : str
      root directory to save QA
    slittxt : str
      Label used for QA

    Returns
    -------
    final_fit : dict
      A dictionary containing all of the information about the fit
    """
    # Check that patt_dict and tcent refer to each other
    if patt_dict['mask'].shape != tcent.shape:
        msgs.error('patt_dict and tcent do not refer to each other. Something is very wrong')

    # Perform final fit to the line IDs
    if thar:
        NIST_lines = (line_lists['NIST'] > 0) & (np.char.find(line_lists['Source'].data, 'MURPHY') >= 0)
    else:
        NIST_lines = line_lists['NIST'] > 0
    ifit = np.where(patt_dict['mask'])[0]

    if outroot is not None:
        plot_fil = outroot + slittxt + '_fit.pdf'
    else:
        plot_fil = None
    # Purge UNKNOWNS from ifit
    imsk = np.ones(len(ifit), dtype=np.bool)
    for kk, idwv in enumerate(np.array(patt_dict['IDs'])[ifit]):
        if np.min(np.abs(self._line_lists['wave'][NIST_lines] - idwv)) > 0.01:
            imsk[kk] = False
    ifit = ifit[imsk]
    # Fit
    try:
        final_fit = fitting.iterative_fitting(spec, tcent, ifit,
                                              np.array(patt_dict['IDs'])[ifit], line_lists[NIST_lines],
                                              patt_dict['bdisp'],match_toler=match_toler, func=func, n_first=n_first,
                                              sigrej_first=sigrej_first,n_final=n_final, sigrej_final=sigrej_final,
                                              plot_fil=plot_fil, verbose=verbose)
    except TypeError:
        # A poor fitting result, this can be ignored.
        return None

    if plot_fil is not None:
        print("Wrote: {:s}".format(plot_fil))

    # Return
    return final_fit


calibfile ='/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_NIRES/NIRES/MF_keck_nires/MasterWaveCalib_A_01_aa.json'
wv_calib_arxiv, par = wavecalib.load_wv_calib(calibfile)
steps= wv_calib_arxiv.pop('steps')
par_dum = wv_calib_arxiv.pop('par')

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
cc_thresh = 0.8
# def archive_reidentify(spec, wv_calib_arxiv, lamps, detections = None, cc_thresh = 0.8
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
narxiv = len(wv_calib_arxiv)
nspec_arxiv = wv_calib_arxiv['0']['spec'].size
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

# For convenience pull out all the spectra from the wv_calib_arxiv archive
spec_arxiv = np.zeros((nspec, narxiv))
wave_soln_arxiv = np.zeros((nspec, narxiv))
wvc_arxiv = np.zeros(narxiv, dtype=float)
disp_arxiv = np.zeros(narxiv, dtype=float)
xrng = np.arange(nspec_arxiv)
for iarxiv in range(narxiv):
    spec_arxiv[:,iarxiv] = wv_calib_arxiv[str(iarxiv)]['spec']
    fitc = wv_calib_arxiv[str(iarxiv)]['fitc']
    xfit = xrng
    fitfunc = wv_calib_arxiv[str(iarxiv)]['function']
    fmin, fmax = wv_calib_arxiv[str(iarxiv)]['fmin'],wv_calib_arxiv[str(iarxiv)]['fmax']
    wave_soln_arxiv[:,iarxiv] = utils.func_val(fitc, xfit, fitfunc, minv=fmin, maxv=fmax)
    wvc_arxiv[iarxiv] = wave_soln_arxiv[nspec_arxiv//2, iarxiv]
    disp_arxiv[iarxiv] = np.median(wave_soln_arxiv[:,iarxiv] - np.roll(wave_soln_arxiv[:,iarxiv], 1))

wv_calib = {}
patt_dict = {}
bad_slits = np.array([], dtype=np.int)
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
        # Match the peaks between the two spectra. This code attempts to compute the stretch if cc > cc_thresh
        success, shift_vec[iarxiv], stretch_vec[iarxiv], ccorr_vec[iarxiv], _, _ = \
            wvutils.xcorr_shift_stretch(spec[:, islit], spec_arxiv[:, iarxiv], cc_thresh=cc_thresh, debug=True)
        # If cc < cc_thresh or if this optimization failed, don't reidentify from this arxiv spectrum
        if success != 1:
            continue
        # Estimate wcen and disp for this slit based on its shift/stretch relative to the archive slit
        disp[iarxiv] = disp_arxiv[iarxiv] / stretch_vec[iarxiv]
        wcen[iarxiv] = wvc_arxiv[iarxiv] - shift_vec[iarxiv]*disp[iarxiv]
        # For each peak in the arxiv spectrum, identify the corresponding peaks in the input islit spectrum. Do this by
        # transforming these arxiv slit line pixel locations into the (shifted and stretched) input islit spectrum frame
        arxiv_det = wv_calib_arxiv[str(iarxiv)]['xfit']
        arxiv_det_ss = arxiv_det*stretch_vec[iarxiv] + shift_vec[iarxiv]
        if debug:
            plt.figure(figsize=(14, 6))
            tampl_slit = np.interp(slit_det, xrng, spec[:, islit])
            plt.plot(xrng, spec[:, islit], color='red', drawstyle='steps-mid', label='input arc',linewidth=1.0, zorder=10)
            plt.plot(slit_det, tampl_slit, 'r.', markersize=10.0, label='input arc lines', zorder=10)
            tampl_arxiv = np.interp(arxiv_det, xrng, spec_arxiv[:, iarxiv])
            plt.plot(xrng, spec_arxiv[:, iarxiv], color='black', drawstyle='steps-mid', linestyle=':',
                     label='arxiv arc', linewidth=0.5)
            plt.plot(arxiv_det, tampl_arxiv, 'k+', markersize=8.0, label='arxiv arc lines')
            spec_arxiv_ss = wvutils.shift_and_stretch(spec_arxiv[:, iarxiv], shift_vec[iarxiv], stretch_vec[iarxiv])
            # tampl_ss = np.interp(gsdet_ss, xrng, gdarc_ss)
            for iline in range(arxiv_det_ss.size):
                plt.plot([arxiv_det[iline], arxiv_det_ss[iline]], [tampl_arxiv[iline], tampl_arxiv[iline]],
                         color='cornflowerblue', linewidth=1.0)
            plt.plot(xrng, spec_arxiv_ss, color='black', drawstyle='steps-mid', label='arxiv arc shift/stretch',linewidth=1.0)
            plt.plot(arxiv_det_ss, tampl_arxiv, 'k.', markersize=10.0, label='predicted arxiv arc lines')
            plt.title(
                'Cross-correlation of input slit # {:d}'.format(islit + 1) + ' and arxiv slit # {:d}'.format(iarxiv + 1) +
                ': ccor = {:5.3f}'.format(ccorr_vec[iarxiv]) +
                ', shift = {:6.1f}'.format(shift_vec[iarxiv]) +
                ', stretch = {:5.4f}'.format(stretch_vec[iarxiv]) +
                ', wv_cen = {:7.1f}'.format(wcen[iarxiv]) +
                ', disp = {:5.3f}'.format(disp[iarxiv]))
            plt.ylim(-5.0, 1.5 *spec[:, islit].max())
            plt.legend()
            plt.show()

        # Calculate wavelengths for all of the gsdet detections
        wvval_arxiv= utils.func_val(wv_calib_arxiv[str(iarxiv)]['fitc'], arxiv_det,wv_calib_arxiv[str(iarxiv)]['function'],
                                    minv=wv_calib_arxiv[str(iarxiv)]['fmin'], maxv=wv_calib_arxiv[str(iarxiv)]['fmax'])

        # Loop over the current slit line pixel detections and find the nearest arxiv spectrum line
        for iline in range(slit_det.size):
            pdiff = np.abs(slit_det[iline] - arxiv_det_ss)
            bstpx = np.argmin(pdiff)
            # If a match is found within 2 pixels, consider this a successful match
            if pdiff[bstpx] < 2.0:
                # Using the arxiv arc wavelength solution, search for the nearest line in the line list
                bstwv = np.abs(wvdata - wvval_arxiv[bstpx])
                # This is probably not a good match
                if bstwv[np.argmin(bstwv)] > 2.0 * disp_arxiv[iarxiv]:
                    continue
                lindex = np.append(lindex, np.argmin(bstwv))  # index in the line list self._wvdata
                dindex = np.append(dindex, iline)

    if np.all(wcen == 0.0) or len(np.unique(lindex)) < 3:
        wv_calib[str(islit)] = {}
        bad_slits = np.append(bad_slits,islit)
        continue

    # Finalize the best guess of each line
    # Initialise the patterns dictionary, min_nsig not used anywhere
    patt_dict_slit = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0., min_nsig=sigdetect,
                     mask=np.zeros(slit_det.size, dtype=np.bool))
    patt_dict_slit['sign'] = 1 # This is not used anywhere
    patt_dict_slit['bwv'] = np.median(wcen[wcen != 0.0])
    patt_dict_slit['bdisp'] = np.median(disp[disp != 0.0])
    patterns.solve_triangles(slit_det, wvdata, dindex, lindex, patt_dict=patt_dict_slit)

    if debug:
        tmp_list = table.vstack([line_lists, unknwns])
        qa.match_qa(spec[:, islit], slit_det, tmp_list, patt_dict_slit['IDs'], patt_dict_slit['scores'])


    # Use only the perfect IDs
    iperfect = np.array(patt_dict_slit['scores']) != 'Perfect'
    patt_dict_slit['mask'][iperfect] = False
    patt_dict_slit['nmatch'] = np.sum(patt_dict_slit['mask'])
    if patt_dict_slit['nmatch'] < 3:
        patt_dict_slit['acceptable'] = False

    # Check if a solution was found
    if not patt_dict_slit['acceptable']:
        wv_calib[str(islit)] = {}
        bad_slits = np.append(bad_slits,islit)
        continue
    final_fit = fit_slit(spec[:,islit], patt_dict_slit, slit_det)
    if final_fit is None:
        # This pattern wasn't good enough
        wv_calib[str(islit)] = {}
        bad_slits = np.append(bad_slits, islit)
        continue
    if final_fit['rms'] > rms_threshold:
        msgs.warn('---------------------------------------------------' + msgs.newline() +
                  'Reidentify report for slit {0:d}/{1:d}:'.format(islit + 1, nslits) + msgs.newline() +
                  '  Poor RMS ({0:.3f})! Need to add additional spectra to arxiv to improve fits'.format(final_fit['rms']) + msgs.newline() +
                  '---------------------------------------------------')
        bad_slits = np.append(bad_slits, islit)
        # Note this result in new_bad_slits, but store the solution since this might be the best possible
    patt_dict[str(islit)] = copy.deepcopy(patt_dict)
    wv_calib[str(islit)] = copy.deepcopy(final_fit)
    if debug:
        yplt = utils.func_val(final_fit['fitc'], xrng, 'legendre', minv=0.0, maxv=1.0)
        plt.plot(final_fit['xfit'], final_fit['yfit'], 'bx')
        plt.plot(xrng, yplt, 'r-')
        plt.show()
    # return wv_calib, patt_dict, bad_slits