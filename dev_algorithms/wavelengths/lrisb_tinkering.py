# imports
from matplotlib import pyplot as plt
import os
import numpy as np
import scipy
import pdb

from linetools import utils as ltu

from pypeit import utils
from pypeit.core.wavecal import wvutils
from pypeit.core.wavecal import waveio

wfile = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT/Keck_LRIS_blue/multi_600_4000_d560',
                     'MF_keck_lris_blue', 'MasterWaveCalib_A_1_01.json')
wv_dict = ltu.loadjson(wfile)

lamps = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']
llist = waveio.load_line_lists(lamps)

# Build the 'master'
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
ncomb = nwspec.size
xvals = np.arange(ncomb)

def get_shift(slit):
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


def target_lines(slit, trim=0.05):
    tspec = wv_dict[str(slit)]['spec']
    nspec = len(tspec)
    # Shift
    npad, shift_cc = get_shift(slit)
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



