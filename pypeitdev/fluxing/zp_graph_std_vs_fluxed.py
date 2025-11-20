
import numpy as np
from matplotlib import pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units

from pypeit.io import fits_open
from pypeit.sensfunc import SensFunc
from pypeit import utils
from pypeit.core import fitting, flux_calib
from pypeit.specobjs import SpecObjs

import argparse
import os.path

from pypeit.spectrographs.util import load_spectrograph




def build_figure():

    utils.pyplot_rcparams()
    fig = plot.figure(figsize=(12,8))
    axis = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    return fig, axis


def create_plot(axis, color, label, x, y, linewidth=2.5, linestyle='solid'):



    axis.plot(x, y, color=color, linewidth=linewidth, linestyle=linestyle, label=label)
    xmin = (0.98*x.min())
    xmax = (1.02*x.max())
    ymin = (0.98*y.min())
    ymax = (1.05*y.max())

    return xmin, xmax, ymin, ymax

def setup_axes(axis,xmin, ymin, xmax, ymax, args):
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    axis.legend()
    axis.set_xlabel('Wavelength (Angstroms)')
    axis.set_ylabel(f'Flux (1e-17 erg/s/cm^2/Ang)')
    axis.set_title(f'PypeIt DEIMOS fluxed standard star spectrum vs archived spectrum')


def get_sobj_legend(args, sobj):
    spec = load_spectrograph("keck_deimos")
    headarr = spec.get_headarr(args.raw_data)
    dispangle = spec.get_meta_value(headarr, "dispangle")
    filter = spec.get_meta_value(headarr, "filter1")
    grating = spec.get_meta_value(headarr, "dispname")
    target = spec.get_meta_value(headarr, "target")
    det = sobj['DET'][0]
    return f"Fluxed {os.path.basename(args.fluxed_standard)} x {args.scale} {target} {grating} dispangle {round(dispangle)} filter {filter} det {det}"

def find_nearby_value(x, y, value):
    sorted_idx = np.argsort(x)
    i = np.searchsorted(x, [value], sorter=sorted_idx)[0]
    return y[sorted_idx][i]

def parse_args(options=None, return_parser=False):
    parser = argparse.ArgumentParser(description='Graph a fluxed spectrum of a standard star vs to the archived spectrum of that star.')
    parser.add_argument("fluxed_standard", type=str,
                        help="Flux calibrated spec1d file of a standard star.")
    parser.add_argument("raw_data", type=str)
    parser.add_argument("scale", type=float)
    parser.add_argument("pdfFile", type=str)
    parser.add_argument("names", type=str, nargs="+")
    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """ Executes sensitivity function computation.
    """
    sobjs = SpecObjs.from_fitsfile(args.fluxed_standard)

    ra = sobjs.header['RA']
    dec = sobjs.header['DEC']

    std_dict = flux_calib.get_standard_spectrum(ra=float(ra), dec=float(dec))
    std_wave = std_dict['wave'] / (1. * units.AA)
    std_flux = std_dict['flux'] / (1. * std_dict['flux'].unit)

    nplots = len(args.names)# + 1
    xmin = np.zeros(nplots)
    ymin = np.zeros(nplots)
    xmax = np.zeros(nplots)
    ymax = np.zeros(nplots)
    i = 0
    with PdfPages(args.pdfFile) as pdf:

            fig, axis = build_figure()

            for name in args.names:
                idx = sobjs.NAME == name
                sobj = sobjs[idx]
                mask = sobj['OPT_MASK']

            for name in args.names:
                idx = sobjs.NAME == name
                sobj = sobjs[idx]
                mask = sobj['OPT_MASK']
                x = sobj['OPT_WAVE'][mask]
                y = sobj['OPT_FLAM'][mask] * args.scale
                legend = get_sobj_legend(args, sobj)
                (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, None, legend, x, y)#, linestyle='dashed')
                i+=1 
            
            wave_min = np.min(xmin[0:i])
            wave_max = np.max(xmax[0:i])
            wave_mask = (std_wave >= wave_min) & (std_wave <= wave_max)
            #(xmin[i], xmax[i], ymin[i], ymax[i]) = 
            create_plot(axis, None, f"Archived spectrum for {std_dict['name']} from {std_dict['std_source']}.", std_wave[wave_mask], std_flux[wave_mask])#, linestyle='dashed')
            #i+=1 
           
            setup_axes(axis, np.min(xmin), np.min(ymin), np.max(xmax),np.max(ymax), args)

            pdf.savefig(fig)

            plot.show()

def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
