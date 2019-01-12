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

lamps = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']
llist = waveio.load_line_lists(lamps)

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
    tbl.write(outfile, overwrite=True)
    print("Wrote: {}".format(outfile))

def get_shift(wv_dict, slit, nwspec):
    ncomb = nwspec.size
    tspec = wv_dict[str(slit)]['spec']  # list
    # Pad
    pspec = np.zeros_like(nwspec)
    nspec = len(tspec)
    npad = ncomb - nspec
    pspec[npad // 2:npad // 2 + len(tspec)] = tspec
    #
    result_out, shift_out, stretch_out, corr_out, shift_cc, corr_cc = wvutils.xcorr_shift_stretch(
        nwspec, pspec)
    # Return
    return npad, shift_cc


def target_lines(wv_dict, slit, nwwv, nwspec, trim=0.05):
    tspec = wv_dict[str(slit)]['spec']
    nspec = len(tspec)
    # Shift
    npad, shift_cc = get_shift(wv_dict, slit, nwspec)
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
    all_tcent, all_ecent, cut_tcent, icut, arc_cont_sub = wvutils.arc_lines_from_spec(mspec)
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
    targ_y = mspec[np.round(all_tcent).astype(int)][mask]

    # Return
    return np.array(tspec), mspec, mwv, targ_wv, targ_y

'''
# GLOBAL ME
tspec, mspec, mwv, targ_wv, targ_y = target_lines(2)
all_tcent, all_ecent, cut_tcent, icut, arc_cont_sub = wvutils.arc_lines_from_spec(tspec)
targ_grid = np.outer(targ_wv, np.ones(all_tcent.size))
'''



def loss_func(theta, all_tcent, targ_grid, targ_y):
    func = 'legendre'
    fmin, fmax = 0., 1.
    new_wv = utils.func_val(theta, all_tcent, func, minx=fmin, maxx=fmax)
    #
    new_grid = np.outer(np.ones(targ_grid.shape[0]), new_wv)
    diff = np.abs(targ_grid - new_grid)
    mgrid = np.min(diff, axis=1)
    # Could sigclip here

    # Loss metric
    loss = mgrid ** 2 * targ_y

    return np.sum(loss)


def fitme(bounds, all_tcent, targ_grid, targ_y):
    args=(all_tcent, targ_grid, targ_y)
    result = scipy.optimize.differential_evolution(loss_func, bounds, tol=1e-4,
                                                   args=args,
                                                   disp=True, polish=True, seed=None)
    # return
    return result

def do_it_all(slit, instr, plot_fil=None):
    if instr == '600':
        wv_dict, template = load_600()
    nwwv = template['wave'].data
    nwspec = template['flux'].data
    # Snippet
    tspec, mspec, mwv, targ_wv, targ_y = target_lines(wv_dict, slit, nwwv, nwspec)
    # Find new lines
    all_tcent, all_ecent, cut_tcent, icut, arc_cont_sub = wvutils.arc_lines_from_spec(tspec)
    # Grab initial guess
    n_order = 3
    func = 'legendre'
    sigrej_first = 3.
    fmin, fmax = 0., 1.
    xfit = np.arange(mwv.size)
    yfit = mwv
    initial_mask = mwv <= 0.
    i1 = np.where(~initial_mask)[0][0]
    disp = mwv[i1 + 1] - mwv[i1]
    #
    mask, fit = utils.robust_polyfit(xfit, yfit, n_order, function=func, sigma=sigrej_first,
                                     initialmask=initial_mask,
                                     minx=fmin, maxx=fmax, verbose=True)  # , weights=wfit)
    # Target grid
    targ_grid = np.outer(targ_wv, np.ones(all_tcent.size))
    # Bounds
    bounds = [[fit[0] - 50., fit[0] + 50],
              [fit[1] * 0.8, fit[1] * 1.2],
              [fit[2] * -2, fit[2] * 2],
              [fit[3] * -2, fit[3] * 2]]
    # Differential evolution
    diff_result = fitme(bounds, all_tcent, targ_grid, targ_y)
    final_guess = utils.func_val(diff_result.x, all_tcent, func, minx=fmin, maxx=fmax)
    # Match to line list
    IDs = []
    dwv = mwv[i1 + 1] - mwv[i1]
    lmask = np.zeros(final_guess.size, dtype=bool)
    for kk, iwv in enumerate(final_guess):
        imin = np.argmin(np.abs(iwv - llist['wave']))
        if np.abs(iwv - llist['wave'][imin]) < dwv:
            lmask[kk] = True
            IDs.append(llist['wave'][imin])
    # Final fit
    final_fit = fitting.iterative_fitting(tspec, all_tcent, np.where(lmask)[0], np.array(IDs), llist, dwv,
                              verbose=True, plot_fil=plot_fil, n_first=2)
    # Return
    return final_fit



# Command line execution
if __name__ == '__main__':
    build_600()





