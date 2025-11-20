""" Make a plot or two to check the gains """
from astropy.io.fits.util import _extract_number
import numpy as np
import os

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from pypeit import flatfield

from IPython import embed

def calc_cuts(master_file:str, namp:int=4):
    if namp != 4:
        raise IOError("Not ready for this value yet!")
    # Load
    flatImages = flatfield.FlatImages.from_file(master_file)

    # Mid points
    spat_mid = flatImages.pixelflat_raw.shape[1]//2
    spec_mid = flatImages.pixelflat_raw.shape[0]//2

    # Cuts
    L1L2_spec = np.median(flatImages.pixelflat_raw[:,spat_mid-3:spat_mid], axis=1)
    U1U2_spec = np.median(flatImages.pixelflat_raw[:,spat_mid:spat_mid+3], axis=1)

    L2U2_spat = np.median(flatImages.pixelflat_raw[spec_mid:spec_mid+3, :], axis=0)
    L1U1_spat = np.median(flatImages.pixelflat_raw[spec_mid-3:spec_mid, :], axis=0)

    return L1L2_spec, U1U2_spec, L1U1_spat, L2U2_spat

def calc_gains(master_file:str, namp:int=4):
    if namp != 4:
        raise IOError("Not ready for this value yet!")
    # Load
    flatImages = flatfield.FlatImages.from_file(master_file)

    # Cut me
    L1L2_spec, U1U2_spec, L1U1_spat, L2U2_spat = calc_cuts(master_file, namp=namp)

    # L2 / L1
    rtio_spat = L2U2_spat / L1U1_spat
    half_spat = rtio_spat.size // 2

    # Assume longslit 
    L2_L1 = np.median(rtio_spat[half_spat-250:half_spat])
    correct_L2L1 = 1./L2_L1

    # U2 / U1
    U2_U1 = np.median(rtio_spat[half_spat:half_spat+250])
    correct_U2U1 = 1./U2_U1

    # U1 / L1 
    rtio_spec = U1U2_spec / L1L2_spec
    half_spec = rtio_spec.size // 2

    U1_L1 = np.median(rtio_spec[half_spec-250:half_spec])
    correct_U1L1 = 1./U1_L1

    # U1 / L1 

    # Time to print
    print(f"Correct L2 by: {correct_L2L1}")
    print(f"Correct U1 by: {correct_U1L1}")
    print(f"Correct U2 by: {correct_U1L1 * correct_U2U1}")

def plot_gains(master_file:str, namp:int=4, outfile:str='gain_plots.png'):
    if namp != 4:
        raise IOError("Not ready for this value yet!")

    # Cut me
    left_spec, right_spec, bot_spat, top_spat = calc_cuts(master_file, namp=namp)

    # Ratios
    rtio_spec = left_spec/right_spec
    rtio_spat = top_spat/bot_spat

    # Plot
    fig = plt.figure(figsize=(9, 5))
    plt.clf()
    gs = gridspec.GridSpec(1,2)

    # Spec
    for kk, rtio, lbl, clr in zip(np.arange(2), [rtio_spec, rtio_spat], 
                                  ['spec', 'spat'], ['k', 'r']):
        ax= plt.subplot(gs[kk])
        # 
        ax.plot(rtio, color=clr)
        ax.set_xlabel(lbl)
        ax.set_ylabel('ratio')
        ax.set_ylim(0.95, 1.05)
        if lbl == 'spat':
            ax.set_xlim(1400, 2600)


    plt.tight_layout(pad=0.5, h_pad=0.5, w_pad=0.5)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))


# Calc gains for corrections using a longslit
#rdx_path = './'
rdx_path = '/scratch/REDUX/Keck/LRIS/new_LRISr/keck_lris_red_mark4_A'
orig_master_file = os.path.join(rdx_path, 'Masters', 'MasterFlat_A_1_01.fits')
calc_gains(orig_master_file)

# Plot Gains for the corrected image
rdx_path = '/scratch/REDUX/Keck/LRIS/new_LRISr/keck_lris_red_mark4_A'
#rdx_path = '/scratch/REDUX/Keck/LRIS/new_LRISr/keck_lris_red_mark4_C'
new_master_file = os.path.join(rdx_path, 'Masters', 'MasterFlat_A_1_01.fits')
plot_gains(new_master_file)
