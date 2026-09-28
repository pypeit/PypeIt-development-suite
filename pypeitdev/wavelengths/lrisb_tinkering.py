# imports
from matplotlib import pyplot as plt
import os
import numpy as np
import scipy
import pdb

from astropy.table import Table

from linetools import utils as ltu

from pypeit import utils
from pypeit.core.wavecal import wvutils
from pypeit.core.wavecal import waveio
from pypeit.core.wavecal import fitting
from pypeit.core.wavecal import autoid
from pypeit.core import pydl

lamps = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']
llist = waveio.load_line_lists(lamps)

def load_300(wv_dict_only=False):
    # Build the 'master'
    wfile = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT/Keck_LRIS_blue/multi_300_5000_d680', 'MF_keck_lris_blue',
                         'MasterWaveCalib_A_1_01.json')
    wv_dict = ltu.loadjson(wfile)
    if wv_dict_only:
        return wv_dict
    # Load template
    template = Table.read('keck_lris_blue_300_d680.fits')
    return wv_dict, template

def load_600(wv_dict_only=False):
    # Build the 'master'
    wfile = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT/Keck_LRIS_blue/multi_600_4000_d560',
                         'MF_keck_lris_blue', 'MasterWaveCalib_A_1_01.json')
    wv_dict = ltu.loadjson(wfile)
    if wv_dict_only:
        return wv_dict
    # Load template
    template = Table.read('keck_lris_blue_600_d560.fits')
    return wv_dict, template

def build_600(outfile='keck_lris_blue_600_d560.fits'):
    wv_dict = load_600(wv_dict_only=True)
    #
    good = [0,7]
    lcut = 4500.
    yvals = []
    lvals = []
    for kk, slit in enumerate(good):
        iwv_calib = wv_dict[str(slit)]
        x = np.arange(iwv_calib['nspec'])
        tmpwv = utils.func_val(iwv_calib['fitc'], x, iwv_calib['function'],
                                           minx=iwv_calib['fmin'], maxx=iwv_calib['fmax'])
        #
        if kk == 0:
            gdi = tmpwv < lcut
        else:
            gdi = tmpwv > lcut
        # Save
        yvals.append(np.array(iwv_calib['spec'])[gdi])
        lvals.append(tmpwv[gdi])
    nwspec = np.concatenate(yvals)
    nwwv = np.concatenate(lvals)
    # Write to disk
    tbl = Table()
    tbl['wave'] = nwwv
    tbl['flux'] = nwspec
    tbl.meta['INSTR'] = 'keck_lris_blue_600_d560'
    tbl.meta['BINNING'] = 2
    tbl.write(outfile, overwrite=True)
    print("Wrote: {}".format(outfile))

def get_shift(wv_dict, slit, nwspec, debug=True, subpix=None):
    ncomb = nwspec.size
    tspec = np.array(wv_dict[str(slit)]['spec'])  # list
    if subpix is not None:
        tspec = tspec[subpix[0]:subpix[1]]
    # Pad
    pspec = np.zeros_like(nwspec)
    nspec = len(tspec)
    npad = ncomb - nspec
    pspec[npad // 2:npad // 2 + len(tspec)] = tspec
    #
    shift_cc, corr_cc = wvutils.xcorr_shift(nwspec, pspec)
    print("Shift = {}; cc = {}".format(shift_cc, corr_cc))
    if debug:
        xvals = np.arange(ncomb)
        plt.clf()
        ax = plt.gca()
        #
        ax.plot(xvals, nwspec)
        ax.plot(xvals, np.roll(pspec, int(shift_cc)), 'k')
        plt.show()
    # Return
    return npad, shift_cc


def target_lines(wv_dict, slit, nwwv, nwspec, trim=0.05, verbose=True, subpix=None, **kwargs):
    tspec = np.array(wv_dict[str(slit)]['spec'])
    if subpix is not None:
        tspec = tspec[subpix[0]:subpix[1]]
    nspec = len(tspec)
    # Shift
    npad, shift_cc = get_shift(wv_dict, slit, nwspec, subpix=subpix)
    i0 = npad // 2 + int(shift_cc)
    # Master
    if i0 < 0: # Pad?
        mspec = np.concatenate([np.zeros(-1*i0), nwspec[0:i0+nspec]])
        mwv = np.concatenate([np.zeros(-1*i0), nwwv[0:i0+nspec]])
    elif (i0+nspec) > nwspec.size: # Pad?
        mspec = np.concatenate([nwspec[i0:], np.zeros(nspec-nwspec.size+i0)])
        mwv = np.concatenate([nwwv[i0:], np.zeros(nspec-nwspec.size+i0)])
    else:
        mspec = nwspec[i0:i0 + nspec]
        mwv = nwwv[i0:i0 + nspec]
    # Peaks
    all_tcent, all_ecent, cut_tcent, icut, arc_cont_sub = wvutils.arc_lines_from_spec(mspec, **kwargs)
    # Wavelengths
    gdwv = mwv[np.round(all_tcent).astype(int)]
    # Match to line list
    gdm = np.where(mwv > 0.)[0][0] # To deal with initial padding
    dwv = mwv[gdm+1] - mwv[gdm]
    lmask = np.zeros(gdwv.size, dtype=bool)
    for kk, iwv in enumerate(gdwv):
        if np.min(np.abs(iwv - llist['wave'])) < dwv:
            lmask[kk] = True
    # Trim off tfrac on each side
    tmask = (all_tcent / mspec.size > trim) & (all_tcent / mspec.size < (1 - trim))
    # Targets
    mask = tmask & lmask
    targ_wv = gdwv[mask]
    targ_x = all_tcent[mask]
    targ_y = mspec[np.round(all_tcent).astype(int)][mask]

    if verbose:
        print("Target lines: {}".format(targ_wv))

    # Return
    return tspec, mspec, mwv, targ_wv, targ_x, targ_y

'''
# GLOBAL ME
tspec, mspec, mwv, targ_wv, targ_y = target_lines(2)
all_tcent, all_ecent, cut_tcent, icut, arc_cont_sub = wvutils.arc_lines_from_spec(tspec)
targ_grid = np.outer(targ_wv, np.ones(all_tcent.size))
'''



def loss_func(theta, all_tcent, targ_grid, targ_y, orig_coeff, nfit, verbose=False, sigclip=True):
    func = 'legendre'
    fmin, fmax = 0., 1.
    # Replace
    coeff = orig_coeff.copy()
    coeff[0:nfit] = theta
    new_wv = utils.func_val(coeff, all_tcent, func, minx=fmin, maxx=fmax)
    #
    new_grid = np.outer(np.ones(targ_grid.shape[0]), new_wv)
    diff = np.abs(targ_grid - new_grid)
    mgrid = np.min(diff, axis=1)
    if verbose:
        print(mgrid)
    # Could sigclip here
    if sigclip:
        iIter = 0
        maxiter = mgrid.size//2
        qdone = False
        thismask = np.ones(mgrid.size, dtype=bool)
        while (not qdone) and (iIter < maxiter):
            thismask, qdone = pydl.djs_reject(np.zeros(mgrid.size), mgrid,
                                              outmask=thismask, lower=3.,upper=3.,
                                              maxrej=1, sticky=True)
            iIter += 1

    # Loss metric
    loss = mgrid ** 2 * targ_y

    return np.sum(loss)


def fitme(bounds, all_tcent, targ_grid, targ_y, orig_coeff):
    args=(all_tcent, targ_grid, targ_y, orig_coeff, len(bounds))
    result = scipy.optimize.differential_evolution(loss_func, bounds, tol=1e-4,
                                                   args=args,
                                                   disp=True, polish=True, seed=None)
    # return
    return result


def do_it_all(slit, instr, plot_fil=None, IDtol=1., subpix=None):
    if instr == '600':
        wv_dict, template = load_600()
        fwhm = 4.
        n_order = 3
        n_first = 2
    elif instr == '300':
        wv_dict, template = load_300()
        fwhm = 8.
        n_order = 3
        n_first = 3
    else:
        pdb.set_trace()
    nwwv = template['wave'].data
    nwspec = template['flux'].data
    # Snippet
    tspec, mspec, mwv, targ_wv, targ_x, targ_y = target_lines(wv_dict, slit, nwwv, nwspec,
                                                              fwhm=fwhm, subpix=subpix)
    # Find new lines
    all_tcent, all_ecent, cut_tcent, icut, arc_cont_sub = wvutils.arc_lines_from_spec(tspec, fwhm=fwhm,
                                                                                      debug=True)
    # Grab initial guess
    func = 'legendre'
    sigrej_first = 3.
    fmin, fmax = 0., 1.
    xfit = np.arange(mwv.size)
    yfit = mwv
    initial_mask = (mwv > targ_wv.min()) & (mwv < targ_wv.max())
    i1 = np.where(initial_mask)[0][0]
    disp = mwv[i1 + 1] - mwv[i1]
    #
    mask, fit = utils.robust_polyfit(xfit, yfit, n_order, function=func, sigma=sigrej_first,
                                     initialmask=~initial_mask,
                                     minx=fmin, maxx=fmax, verbose=True)  # , weights=wfit)
    rms_ang = utils.calc_fit_rms(xfit[mask == 0], yfit[mask == 0], fit, func, minx=fmin, maxx=fmax)
    #                                     weights=wfit[mask == 0])
    rms_pix = rms_ang / disp
    print("Initial fit: {};  RMS = {}".format(fit, rms_pix))
    # Target grid
    targ_grid = np.outer(targ_wv, np.ones(all_tcent.size))
    # Check the loss
    loss = loss_func(fit[0:2], all_tcent, targ_grid, targ_y, fit, 2)
    # Bounds
    bounds = [[fit[0] - 50., fit[0] + 50],
              [fit[1]*0.9, fit[1]*1.1]]
    #for ifit in fit[2:]:
    #    bounds.append([np.abs(ifit)*-2, np.abs(ifit)*2])
    # Differential evolution
    #diff_result = fitme(bounds, all_tcent, targ_grid, targ_y)
    diff_result = fitme(bounds, all_tcent, targ_grid, targ_y, fit.copy())
    new_fit = fit.copy()
    new_fit[:len(bounds)] = diff_result.x
    loss = loss_func(new_fit[0:len(bounds)], all_tcent, targ_grid, targ_y, fit, len(bounds), verbose=True)
    print("Refined fit: {}".format(diff_result.x))
    # Now cut on saturation
    all_tcent, all_ecent, cut_tcent, icut, arc_cont_sub = wvutils.arc_lines_from_spec(tspec, fwhm=fwhm,
                                                                                      nonlinear_counts=55000)
                                                                                      #debug=True)
    final_guess = utils.func_val(new_fit, all_tcent, func, minx=fmin, maxx=fmax)
    # Match to line list
    IDs = []
    dwv = mwv[i1 + 1] - mwv[i1]
    lmask = np.zeros(final_guess.size, dtype=bool)
    for kk, iwv in enumerate(final_guess):
        imin = np.argmin(np.abs(iwv - llist['wave']))
        if np.abs(iwv - llist['wave'][imin]) < dwv*IDtol:
            lmask[kk] = True
            IDs.append(llist['wave'][imin])
    print("IDs: {}".format(IDs))
    # Final fit
    final_fit = fitting.iterative_fitting(tspec, all_tcent, np.where(lmask)[0], np.array(IDs), llist, dwv,
                              verbose=True, plot_fil=plot_fil, n_first=n_first)
    # Return
    return final_fit


def reid_in_steps(slit, instr, plot_fil=None, IDtol=1., subpix=None):
    if instr == '600':
        wv_dict, template = load_600()
        fwhm = 4.
        n_first = 2
        nsnippet = 2
    elif instr == '300':
        wv_dict, template = load_300()
        fwhm = 8.
        n_first = 3
        nsnippet = 2
    else:
        pdb.set_trace()
    # Template
    nwwv = template['wave'].data
    nwspec = template['flux'].data

    # Full snippet (padded as need be)
    tspec, mspec, mwv, targ_wv, targ_x, targ_y = target_lines(wv_dict, slit, nwwv, nwspec,
                                                              fwhm=fwhm, subpix=subpix)
    # Loop on snippets
    nsub = tspec.size // nsnippet
    sv_det, sv_IDs = [], []
    for kk in range(nsnippet):
        # Construct
        i0 = nsub * kk
        i1 = min(nsub*(kk+1), tspec.size)
        tsnippet = tspec[i0:i1]
        msnippet = mspec[i0:i1]
        mwvsnippet = mwv[i0:i1]
        # Run reidentify
        detections, spec_cont_sub, patt_dict = autoid.reidentify(tsnippet, msnippet, mwvsnippet,
                                                                 llist, 1, debug_xcorr=False,
                                                                 nonlinear_counts=55000.,
                                                                 debug_reid=True,  # verbose=True,
                                                                 match_toler=2.5,
                                                                 cc_thresh=0.1, fwhm=fwhm)
        # Deal with IDs
        sv_det.append(i0 + detections)
        sv_IDs.append(patt_dict['IDs'])

    # Collate and proceed
    dets = np.concatenate(sv_det)
    IDs = np.concatenate(sv_IDs)
    #
    i1 = np.where(mwv > 0.)[0][0]
    dwv = mwv[i1+1] - mwv[i1]
    #
    gd_det = np.where(IDs > 0.)[0]
    final_fit = fitting.iterative_fitting(tspec, dets, gd_det, IDs[gd_det], llist, dwv,
                                          verbose=True, plot_fil=plot_fil, n_first=n_first)
    return final_fit

# Command line execution
if __name__ == '__main__':
    build_600()





