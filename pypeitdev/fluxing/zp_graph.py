
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

def create_zp_plot(axis, label, det_num, sf):


    wave = sf.sens[det_num]['SENS_WAVE']
    zp = sf.sens[det_num]['SENS_ZEROPOINT_FIT']
    gpm = sf.sens[det_num]['SENS_ZEROPOINT_FIT_GPM']
    
    return create_plot(axis, None, label, wave[gpm], zp[gpm])

def create_plot(axis, color, label, x, y, linewidth=2.5, linestyle='solid'):

    axis.plot(x, y, color=color, linewidth=linewidth, linestyle=linestyle, label=label)
    xmin = (0.98*x.min())
    xmax = (1.02*x.max())
    #ymin = 14.0
    ymin = (0.98*y.min())
    ymax = min((1.05*y.max()),21.5)

    return xmin, xmax, ymin, ymax



def setup_axes(axis,xmin, ymin, xmax, ymax, title):
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    axis.legend()
    axis.set_xlabel('Wavelength (Angstroms)')
    axis.set_ylabel('Zeropoint Fit (AB mag)')
    axis.set_title(title)

    
def parse_args(options=None, return_parser=False):
    parser = argparse.ArgumentParser(description='Graph a sensitivity function')
    parser.add_argument("--title", type=str,
                        help="Optional Title for the plot", default=None)
    parser.add_argument("--labels", default=[], type=str, nargs="+")
    parser.add_argument("--no-labels", default=False, action="store_true")
    parser.add_argument("pdfFile", type=str, help="Name of PDF file to save.")
    parser.add_argument("filelist", type=Path, nargs="+",
                        help="List of sens funcs to plot.")

            

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """ Executes sensitivity function computation.
    """
    if args.title is None:
        title = 'PypeIt Zeropoint Fit for' + args.filelist[0].stem.replace("_", " ")
    else:
        title = args.title

    sflist =  [SensFunc.from_file(str(file)) for file in args.filelist]
    num_plots = sum([len(sf.sens) for sf in sflist])
    with PdfPages(args.pdfFile) as pdf:

        xmin = np.zeros(num_plots)
        ymin = np.zeros(num_plots)
        xmax = np.zeros(num_plots)
        ymax = np.zeros(num_plots)
        fig, axis = build_figure()
        i = 0  
        for n, sf in enumerate(sflist):

            for det_num in range(len(sf.sens)):
                if args.no_labels:
                    label = None
                elif i >= len(args.labels):
                    label = f"{args.filelist[n]} det {i}"
                else:
                    label = args.labels[i]                    
                (xmin[i], xmax[i], ymin[i], ymax[i]) = create_zp_plot(axis, label, det_num, sf)
                i += 1


        setup_axes(axis, np.min(xmin), np.min(ymin), np.max(xmax),np.max(ymax), title)

        pdf.savefig(fig)

        plot.show()

def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
