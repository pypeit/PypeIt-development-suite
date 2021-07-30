
import sys
import argparse
from pathlib import Path

import numpy as np

from pypeit.io import fits_open
from pypeit.sensfunc import SensFunc
from pypeit.spectrographs import util
from pypeit.core import fitting




def parse_args(options=None, return_parser=False):
    parser = argparse.ArgumentParser(description='Combine DEIMOS sensitivity functions to create a general purpose one.')
    parser.add_argument("--output", type=str,
                        help="Full path name of the output file")        
    parser.add_argument("filelist", type=str, nargs="+",
                        help="List of sensfunc files to combine.")
            

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """ Executes sensitivity function computation.
    """
    sflist = []
    for file in args.filelist:
        sf = SensFunc.from_file(file)
        sf.sensfile = file
        sflist.append(sf)


    sflist.sort(key=lambda x: x.sens['WAVE_MIN'][0])

    if len(sflist)!=3:
        print("Failed. Currently require exaclty 3 sensfunc files.")
        sys.exit(1)

    # Take everything in the first sensfunc
    combined_wave = np.concatenate((sflist[0].sens['SENS_WAVE'][0], sflist[0].sens['SENS_WAVE'][1]))
    combined_zp = np.concatenate((sflist[0].sens['SENS_ZEROPOINT'][0], sflist[0].sens['SENS_ZEROPOINT'][1]))
    combined_zp_gpm = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_GPM'][0], sflist[0].sens['SENS_ZEROPOINT_GPM'][1]))
    combined_zp_fit = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT'][0], sflist[0].sens['SENS_ZEROPOINT_FIT'][1]))
    combined_zp_fit_gpm = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][0], sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][1]))

    # For detector 3 in the second, take everything that doesn't overlap detector 7 in the first
    # Take all of detecter 7 in the second.
    non_overlap = sflist[1].sens['SENS_WAVE'][0] > np.max(sflist[0].sens['SENS_WAVE'][1])
    combined_wave = np.concatenate((combined_wave, sflist[1].sens['SENS_WAVE'][0][non_overlap], sflist[1].sens['SENS_WAVE'][1]))
    combined_zp = np.concatenate((combined_zp, sflist[1].sens['SENS_ZEROPOINT'][0][non_overlap], sflist[1].sens['SENS_ZEROPOINT'][1]))
    combined_zp_gpm = np.concatenate((combined_zp_gpm, sflist[1].sens['SENS_ZEROPOINT_GPM'][0][non_overlap], sflist[1].sens['SENS_ZEROPOINT_GPM'][1]))
    combined_zp_fit = np.concatenate((combined_zp_fit, sflist[1].sens['SENS_ZEROPOINT_FIT'][0][non_overlap], sflist[1].sens['SENS_ZEROPOINT_FIT'][1]))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm, sflist[1].sens['SENS_ZEROPOINT_FIT_GPM'][0][non_overlap], sflist[1].sens['SENS_ZEROPOINT_FIT_GPM'][1]))

    # For the third sens func, ignore detector 3 because it overlaps everything.
    # Take the non overlapping parts of detector 7
    non_overlap = sflist[2].sens['SENS_WAVE'][1] > np.max(sflist[1].sens['SENS_WAVE'][1])
    combined_wave = np.concatenate((combined_wave, sflist[2].sens['SENS_WAVE'][1][non_overlap]))
    combined_zp = np.concatenate((combined_zp, sflist[2].sens['SENS_ZEROPOINT'][1][non_overlap]))
    combined_zp_gpm = np.concatenate((combined_zp_gpm, sflist[2].sens['SENS_ZEROPOINT_GPM'][1][non_overlap]))
    combined_zp_fit = np.concatenate((combined_zp_fit, sflist[2].sens['SENS_ZEROPOINT_FIT'][1][non_overlap]))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm, sflist[2].sens['SENS_ZEROPOINT_FIT_GPM'][1][non_overlap]))

    newsens = SensFunc.empty_sensfunc_table(1, len(combined_wave))
    newsens['SENS_WAVE'] = combined_wave
    newsens['SENS_ZEROPOINT'] = combined_zp
    newsens['SENS_ZEROPOINT_GPM'] = combined_zp_gpm
    newsens['SENS_ZEROPOINT_FIT'] = combined_zp_fit
    newsens['SENS_ZEROPOINT_FIT_GPM'] = combined_zp_fit_gpm
    sflist[0].sens = newsens
    sflist[0].to_file(args.output + "_stitched.fits", overwrite=True)
    #for order in [10, 11, 12, 13, 14, 15]:
    order = 8
    zppf = fitting.robust_fit(combined_wave[combined_zp_gpm], combined_zp[combined_zp_gpm], order)
    zpfpf = fitting.robust_fit(combined_wave[combined_zp_fit_gpm], combined_zp[combined_zp_fit_gpm], order)

    zppf.to_file(args.output + f"_zp_fit_{order}.fits", overwrite=True)
    zpfpf.to_file(args.output + f"_zpfit_fit_{order}.fits", overwrite=True)

def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
