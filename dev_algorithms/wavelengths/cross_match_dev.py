

import numpy as np
# Read in a wavelength solution and compute
from pypeit import wavecalib
from pypeit.core.wavecal import autoid, waveio
from pypeit import utils




def reidentify_old(spec, wv_calib_arxiv, lamps, nreid_min, detections=None, cc_thresh=0.8,cc_local_thresh = 0.8,
               line_pix_tol=2.0, nlocal_cc=11, rms_threshold=0.15, nonlinear_counts=1e10,sigdetect = 5.0,
               use_unknowns=True,match_toler=3.0,func='legendre',n_first=2,sigrej_first=3.0,n_final=4, sigrej_final=2.0,
               seed=None, debug_xcorr=False, debug_reid=False):

    """ Determine  a wavelength solution for a set of spectra based on archival wavelength solutions

    Parameters
    ----------
    spec :  float ndarray (nspec, nslits)
       Array of arc spectra for which wavelength solutions are desired.

    wv_calib_arxiv: dict
       Dictionary containing archival wavelength solutions for a collection of slits/orders to be used to reidentify
       lines and  determine the wavelength solution for spec. This dict is a standard format for PypeIt wavelength solutions
       as created by pypeit.core.wavecal.fitting.iterative_fitting

    lamps: list of strings
       The are line lamps that are on or the name of the linelist that should be used. For example for Shane Kast blue
       this would ['CdI','HgI','HeI']. For X-shooter NIR which calibrates of a custom OH sky line list,
       it is the name of the line list, i.e. ['OH_XSHOOTER']


    Optional Parameters
    -------------------
    detections: float ndarray, default = None
       An array containing the pixel centroids of the lines in the arc as computed by the pypeit.core.arc.detect_lines
       code. If this is set to None, the line detection will be run inside the code.

    cc_thresh: float, default = 0.8
       Threshold for the *global* cross-correlation coefficient between an input spectrum and member of the archive required to
       attempt reidentification. Spectra from the archive with a lower cross-correlation are not used for reidentification

    cc_local_thresh: float, default = 0.8
       Threshold for the *local* cross-correlation coefficient, evaluated at each reidentified line,  between an input
       spectrum and the shifted and stretched archive spectrum above which a line must be to be considered a good line for
       reidentification. The local cross-correlation is evaluated at each candidate reidentified line
       (using a window of nlocal_cc), and is then used to score the the reidentified lines to arrive at the final set of
       good reidentifications

    line_pix_tol: float, default = 2.0
       Matching tolerance in pixels for a line reidentification. A good line match must match within this tolerance to the
       the shifted and stretched archive spectrum, and the archive wavelength solution at this match must be within
       line_pix_tol dispersion elements from the line in line list.

    n_local_cc: int, defualt = 11
       Size of pixel window used for local cross-correlation computation for each arc line. If not an odd number one will
       be added to it to make it odd.

    rms_threshold: float, default = 0.15
       Minimum rms for considering a wavelength solution to be an acceptable good fit. Slits/orders with a larger RMS
       than this are flagged as bad slits

    nonlinear_counts: float, default = 1e10
       Arc lines above this saturation threshold are not used in wavelength solution fits because they cannot be accurately
       centroided

    sigdetect: float, default 5.0
       Sigma threshold above fluctuations for arc-line detection. Arcs are continuum subtracted and the fluctuations are
       computed after continuum subtraction.

    use_unknowns : bool, default = True
       If True, arc lines that are known to be present in the spectra, but have not been attributed to an element+ion,
       will be included in the fit.

    match_toler: float, default = 3.0
       Matching tolerance when searching for new lines. This is the difference in pixels between the wavlength assigned to
       an arc line by an iteration of the wavelength solution to the wavelength in the line list.

    func: str, default = 'legendre'
       Name of function used for the wavelength solution

    n_first: int, default = 2
       Order of first guess to the wavelength solution.

    sigrej_first: float, default = 2.0
       Number of sigma for rejection for the first guess to the wavelength solution.

    n_final: int, default = 4
       Order of the final wavelength solution fit

    sigrej_final: float, default = 3.0
       Number of sigma for rejection for the final fit to the wavelength solution.

    seed: int or np.random.RandomState, optional, default = None
       Seed for scipy.optimize.differential_evolution optimizer. If not specified, the calculation will be seeded
       in a deterministic way from the input arc spectrum spec.

    debug_xcorr: bool, default = False
       Show plots useful for debugging the cross-correlation used for shift/stretch computation

    debug_reid: bool, default = False
       Show plots useful for debugging the line reidentification

    Returns
    -------
    (wv_calib, patt_dict, bad_slits)

    wv_calib: dict
       Wavelength solution for the input arc spectra spec. These are stored in standard pypeit format, i.e.
       each index of spec[:,slit] corresponds to a key in the wv_calib dictionary wv_calib[str(slit)] which yields
       the final_fit dictionary for this slit

    patt_dict: dict
       Arc lines pattern dictionary with some information about the IDs as well as the cross-correlation values

    bad_slits: ndarray, int
       Numpy array with the indices of the bad slits. These are the indices in the input arc spectrum array spec[:,islit]


    Revision History
    ----------------
    November 2018 by J.F. Hennawi. Based on an initial version of this code written by Ryan Cooke.
    """

    # Determine the seed for scipy.optimize.differential_evolution optimizer
    if seed is None:
        # If no seed is specified just take the sum of all the elements and round that to an integer
        seed = np.fmin(int(np.sum(spec)),2**32-1)

    random_state = np.random.RandomState(seed = seed)


    nlocal_cc_odd = nlocal_cc + 1 if nlocal_cc % 2 == 0 else nlocal_cc
    window = 1.0/nlocal_cc_odd* np.ones(nlocal_cc_odd)

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

    nspec, nslits = spec.shape
    narxiv = len(wv_calib_arxiv)
    nspec_arxiv = wv_calib_arxiv['0']['spec'].size
    if nspec_arxiv != nspec:
        msgs.error('Different spectral binning is not supported yet but it will be soon')

    # If the detections were not passed in find the lines in each spectrum
    if detections is None:
        detections = {}
        for islit in range(nslits):
            tcent, ecent, cut_tcent, icut = wvutils.arc_lines_from_spec(spec[:, islit], sigdetect=sigdetect,nonlinear_counts=nonlinear_counts)
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
        fitfunc = wv_calib_arxiv[str(iarxiv)]['function']
        fmin, fmax = wv_calib_arxiv[str(iarxiv)]['fmin'],wv_calib_arxiv[str(iarxiv)]['fmax']
        wave_soln_arxiv[:,iarxiv] = utils.func_val(fitc, xrng, fitfunc, minv=fmin, maxv=fmax)
        wvc_arxiv[iarxiv] = wave_soln_arxiv[nspec_arxiv//2, iarxiv]
        disp_arxiv[iarxiv] = np.median(wave_soln_arxiv[:,iarxiv] - np.roll(wave_soln_arxiv[:,iarxiv], 1))

    wv_calib = {}
    patt_dict = {}
    bad_slits = np.array([], dtype=np.int)

    marker_tuple = ('o','v','<','>','8','s','p','P','*','X','D','d','x')
    color_tuple = ('black','green','red','cyan','magenta','blue','darkorange','yellow','dodgerblue','purple','lightgreen','cornflowerblue')
    marker = itertools.cycle(marker_tuple)
    colors = itertools.cycle(color_tuple)

    # Loop over the slits in the spectrum and cross-correlate each with each arxiv spectrum to identify lines
    for islit in range(nslits):
        slit_det = detections[str(islit)][0]
        line_indx = np.array([], dtype=np.int)
        det_indx = np.array([], dtype=np.int)
        line_cc = np.array([], dtype=float)
        line_iarxiv = np.array([], dtype=np.int)
        wcen = np.zeros(narxiv)
        disp = np.zeros(narxiv)
        shift_vec = np.zeros(narxiv)
        stretch_vec = np.zeros(narxiv)
        ccorr_vec = np.zeros(narxiv)
        for iarxiv in range(narxiv):
            msgs.info('Cross-correlating slit # {:d}'.format(islit + 1) + ' with arxiv slit # {:d}'.format(iarxiv + 1))
            # Match the peaks between the two spectra. This code attempts to compute the stretch if cc > cc_thresh
            success, shift_vec[iarxiv], stretch_vec[iarxiv], ccorr_vec[iarxiv], _, _ = \
                wvutils.xcorr_shift_stretch(spec[:, islit], spec_arxiv[:, iarxiv], cc_thresh=cc_thresh, seed = random_state,
                                            debug=debug_xcorr)
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
            spec_arxiv_ss = wvutils.shift_and_stretch(spec_arxiv[:, iarxiv], shift_vec[iarxiv], stretch_vec[iarxiv])

            if debug_xcorr:
                plt.figure(figsize=(14, 6))
                tampl_slit = np.interp(slit_det, xrng, spec[:, islit])
                plt.plot(xrng, spec[:, islit], color='red', drawstyle='steps-mid', label='input arc',linewidth=1.0, zorder=10)
                plt.plot(slit_det, tampl_slit, 'r.', markersize=10.0, label='input arc lines', zorder=10)
                tampl_arxiv = np.interp(arxiv_det, xrng, spec_arxiv[:, iarxiv])
                plt.plot(xrng, spec_arxiv[:, iarxiv], color='black', drawstyle='steps-mid', linestyle=':',
                         label='arxiv arc', linewidth=0.5)
                plt.plot(arxiv_det, tampl_arxiv, 'k+', markersize=8.0, label='arxiv arc lines')
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
                plt.ylim(1.2*spec[:, islit].min(), 1.5 *spec[:, islit].max())
                plt.legend()
                plt.show()

            # Calculate wavelengths for all of the gsdet detections
            wvval_arxiv= utils.func_val(wv_calib_arxiv[str(iarxiv)]['fitc'], arxiv_det,wv_calib_arxiv[str(iarxiv)]['function'],
                                        minv=wv_calib_arxiv[str(iarxiv)]['fmin'], maxv=wv_calib_arxiv[str(iarxiv)]['fmax'])
            # Compute a "local" zero lag correlation of the slit spectrum and the shifted and stretch arxiv spectrum over a
            # a nlocal_cc_odd long segment of spectrum. We will then uses spectral similarity as a further criteria to
            # decide which lines are good matches
            prod_smooth = scipy.ndimage.filters.convolve1d(spec[:, islit]*spec_arxiv_ss, window)
            spec2_smooth = scipy.ndimage.filters.convolve1d(spec[:, islit]**2, window)
            arxiv2_smooth = scipy.ndimage.filters.convolve1d(spec_arxiv_ss**2, window)
            denom = np.sqrt(spec2_smooth*arxiv2_smooth)
            corr_local = np.zeros_like(denom)
            corr_local[denom > 0] = prod_smooth[denom > 0]/denom[denom > 0]
            corr_local[denom == 0.0] = -1.0

            # Loop over the current slit line pixel detections and find the nearest arxiv spectrum line
            for iline in range(slit_det.size):
                # match to pixel in shifted/stretch arxiv spectrum
                pdiff = np.abs(slit_det[iline] - arxiv_det_ss)
                bstpx = np.argmin(pdiff)
                # If a match is found within 2 pixels, consider this a successful match
                if pdiff[bstpx] < line_pix_tol:
                    # Using the arxiv arc wavelength solution, search for the nearest line in the line list
                    bstwv = np.abs(wvdata - wvval_arxiv[bstpx])
                    # This is a good wavelength match if it is within line_pix_tol disperion elements
                    if bstwv[np.argmin(bstwv)] < line_pix_tol*disp_arxiv[iarxiv]:
                        line_indx = np.append(line_indx, np.argmin(bstwv))  # index in the line list array wvdata of this match
                        det_indx = np.append(det_indx, iline)             # index of this line in the detected line array slit_det
                        line_cc = np.append(line_cc,np.interp(slit_det[iline],xrng,corr_local)) # local cross-correlation at this match
                        line_iarxiv = np.append(line_iarxiv,iarxiv)

        narxiv_used = np.sum(wcen != 0.0)
        if (narxiv_used == 0) or (len(np.unique(line_indx)) < 3):
            wv_calib[str(islit)] = {}
            patt_dict[str(islit)] = {}
            bad_slits = np.append(bad_slits,islit)
            continue

        if debug_reid:
            plt.figure(figsize=(14, 6))
            # Plot a summary of the local x-correlation values for each line on each slit
            for iarxiv in range(narxiv):
                # Only plot those that we actually tried to reidentify (i.e. above cc_thresh)
                if wcen[iarxiv] != 0.0:
                    this_iarxiv = line_iarxiv == iarxiv
                    plt.plot(wvdata[line_indx[this_iarxiv]],line_cc[this_iarxiv],marker=next(marker),color=next(colors),
                             linestyle='',markersize=5.0,label='arxiv slit={:d}'.format(iarxiv))

            plt.hlines(cc_local_thresh, wvdata[line_indx].min(), wvdata[line_indx].max(), color='red', linestyle='--',label='Local xcorr threshhold')
            plt.title('slit={:d}'.format(islit + 1) + ': Local x-correlation for reidentified lines from narxiv_used={:d}'.format(narxiv_used) +
                      ' arxiv slits. Requirement: nreid_min={:d}'.format(nreid_min) + ' matches > threshold')
            plt.xlabel('wavelength from line list')
            plt.ylabel('Local x-correlation coefficient')
            #plt.ylim((0.0, 1.2))
            plt.legend()
            plt.show()

        # Finalize the best guess of each line
        # Initialise the patterns dictionary, min_nsig not used anywhere
        patt_dict_slit = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0., min_nsig=sigdetect,mask=np.zeros(slit_det.size, dtype=np.bool))
        patt_dict_slit['sign'] = 1 # This is not used anywhere
        patt_dict_slit['bwv'] = np.median(wcen[wcen != 0.0])
        patt_dict_slit['bdisp'] = np.median(disp[disp != 0.0])
        patterns.solve_xcorr(slit_det, wvdata, det_indx, line_indx, line_cc, patt_dict=patt_dict_slit,nreid_min=nreid_min,
                             cc_local_thresh=cc_local_thresh)

        if debug_reid:
            tmp_list = table.vstack([line_lists, unknwns])
            qa.match_qa(spec[:, islit], slit_det, tmp_list, patt_dict_slit['IDs'], patt_dict_slit['scores'])

        # Use only the perfect IDs
        iperfect = np.array(patt_dict_slit['scores']) != 'Perfect'
        patt_dict_slit['mask'][iperfect] = False
        patt_dict_slit['nmatch'] = np.sum(patt_dict_slit['mask'])
        if patt_dict_slit['nmatch'] < 3:
            patt_dict_slit['acceptable'] = False

        # Check if an acceptable reidentification solution was found
        if not patt_dict_slit['acceptable']:
            wv_calib[str(islit)] = {}
            patt_dict[str(islit)] = copy.deepcopy(patt_dict_slit)
            bad_slits = np.append(bad_slits,islit)
            continue
        # Perform the fit
        final_fit = fitting.fit_slit(spec[:,islit], patt_dict_slit, slit_det, line_lists, match_toler=match_toler,
                             func=func, n_first=n_first,sigrej_first=sigrej_first,n_final=n_final,
                             sigrej_final=sigrej_final)

        # Did the fit succeed?
        if final_fit is None:
            # This pattern wasn't good enough
            wv_calib[str(islit)] = {}
            patt_dict[str(islit)] = copy.deepcopy(patt_dict_slit)
            bad_slits = np.append(bad_slits, islit)
            continue
        # Is the RMS below the threshold?
        if final_fit['rms'] > rms_threshold:
            msgs.warn('---------------------------------------------------' + msgs.newline() +
                      'Reidentify report for slit {0:d}/{1:d}:'.format(islit + 1, nslits) + msgs.newline() +
                      '  Poor RMS ({0:.3f})! Need to add additional spectra to arxiv to improve fits'.format(final_fit['rms']) + msgs.newline() +
                      '---------------------------------------------------')
            bad_slits = np.append(bad_slits, islit)
            # Note this result in new_bad_slits, but store the solution since this might be the best possible

        # Add the patt_dict and wv_calib to the output dicts
        patt_dict[str(islit)] = copy.deepcopy(patt_dict_slit)
        wv_calib[str(islit)] = copy.deepcopy(final_fit)
        if debug_reid:
            qa.arc_fit_qa(wv_calib[str(islit)])
            #yplt = utils.func_val(final_fit['fitc'], xrng, final_fit['function'], minv=final_fit['fmin'], maxv=final_fit['fmax'])
            #plt.plot(final_fit['xfit'], final_fit['yfit'], 'bx')
            #plt.plot(xrng, yplt, 'r-')
            #plt.show()

    return wv_calib, patt_dict, bad_slits



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