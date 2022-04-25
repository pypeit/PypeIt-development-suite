
import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy.table import Table
from matplotlib import pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages

from pypeit.io import fits_open
from pypeit.sensfunc import SensFunc
from pypeit.spectrographs import util
from pypeit import utils
from pypeit.core import fitting


import sys
import argparse
from pathlib import Path





def build_figure():

    utils.pyplot_rcparams()
    fig = plot.figure(figsize=(12,8))
    axis = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    return fig, axis

def create_zp_plot(axis, orig_file_metadata, det, sf):

    dispangle = orig_file_metadata['CENTRAL_WAVE'][0]
    filter = orig_file_metadata['BLOCKING'][0]
    airmass = orig_file_metadata['AIRMASS'][0]
    orig_path = Path(orig_file_metadata['FILENAME'][0])
    converted_filename = str(orig_path.parent.name + "_" + orig_path.stem)
    std_name = orig_file_metadata['STD_NAME'][0]

    i = 0 if det == 3 else 1

    wave = sf.sens[i]['SENS_WAVE']
    zp = sf.sens[i]['SENS_ZEROPOINT_FIT']
    gpm = sf.sens[i]['SENS_ZEROPOINT_FIT_GPM']
    
    label = f'{converted_filename} {std_name} det {det}, dispangle {dispangle:.3f} filter {filter} airmass {airmass:.3f}'
    return create_plot(axis, None, label, wave[gpm], zp[gpm])

def create_plot(axis, color, label, x, y, linewidth=2.5, linestyle='solid'):

    axis.plot(x, y, color=color, linewidth=linewidth, linestyle=linestyle, label=label)
    xmin = (0.98*x.min())
    xmax = (1.02*x.max())
    ymin = 14.0
    ymax = (1.05*y.max())

    return xmin, xmax, ymin, ymax



def setup_axes(axis,xmin, ymin, xmax, ymax, title):
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    axis.legend()
    axis.set_xlabel('Wavelength (Angstroms)')
    axis.set_ylabel('Zeropoint Fit (AB mag)')
    axis.set_title('PypeIt Zeropoint Fit for' + title)

    
def parse_args(options=None, return_parser=False):
    parser = argparse.ArgumentParser(description='Graph a sensitivity function')
    parser.add_argument("orig_root_dir", type=str,
                        help="Original path of fits files processed by Greg Wirths scripts")        
    parser.add_argument("sens_root_dir", type=str,
                        help="Path of sensfunc fits files.")        
    parser.add_argument("pdfFile", type=str)
    parser.add_argument("filelist", type=str, nargs="+",
                        help="List of files to parse.")

            

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """ Executes sensitivity function computation.
    """
    orig_root_dir = Path(args.orig_root_dir)
    sens_root_dir = Path(args.sens_root_dir)
     
    with PdfPages(args.pdfFile) as pdf:
        for file in args.filelist:
            title = str(Path(file).stem).replace("_", " ")
            with open(file, "r") as f:
                fits_files = [Path(x.strip()) for x in f.readlines() if x.strip() != '']

            fig, axis = build_figure()
            xmin = np.zeros(len(fits_files)*2)
            ymin = np.zeros(len(fits_files)*2)
            xmax = np.zeros(len(fits_files)*2)
            ymax = np.zeros(len(fits_files)*2)

            i = 0  
            for fits_file in fits_files:
                hdul = fits_open(orig_root_dir / fits_file)
                orig_file_metadata = Table(hdul[1].data)
                sens_file = sens_root_dir / fits_file.parent / (f"sens_" + fits_file.name)
                print('Reading sensitivity function from file: {:}'.format(sens_file))
                sf = SensFunc.from_file(sens_file)

                for det in [3, 7]:
                    (xmin[i], xmax[i], ymin[i], ymax[i]) = create_zp_plot(axis, orig_file_metadata, det, sf)
                    i += 1

                #(xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, None, f'{sens_file} zeropoint', sf.wave, sf.zeropoint)
                #i += 1


            title = f' keck_deimos {title}'

            setup_axes(axis, np.min(xmin), np.min(ymin), np.max(xmax),np.max(ymax), title)

            pdf.savefig(fig)

            plot.show()

def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
