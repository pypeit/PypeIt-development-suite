from pypeit.io import fits_open
from pypeit.sensfunc import SensFunc
from pypeit.spectrographs import util
from pypeit.core import fitting
from pathlib import Path



import sys
import argparse

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
    for sf in sflist:
        print(sf.sensfile)
        base_name = Path(sf.sensfile).stem
        for order in [5, 6, 7, 8, 9, 10]:
            pf = fitting.robust_fit(sf.sens['SENS_WAVE'][0], sf.sens['SENS_ZEROPOINT_FIT'][0], order)
            pf.to_file(str(base_name) + f"_fit_{order}.fits")

def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
