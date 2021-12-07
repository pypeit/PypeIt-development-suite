
import argparse
import copy
import os
import datetime

import numpy as np
from astropy.io import fits
from astropy.table import Table

from pypeit.sensfunc import SensFunc
from pypeit.core import fitting
from pypeit.core.wavecal import wvutils
from pypeit.spectrographs.util import load_spectrograph
from pypeit.spectrographs.keck_deimos import load_wmko_std_spectrum
import pypeit.io

from pypeit import utils
from matplotlib import pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages
from stitchdeimos import stitch_sensfunc

def write_stitched_sensfunc(sflist, args, combined_wave, combined_zp_fit,combined_zp_fit_gpm):
    
    newsf = copy.deepcopy(sflist[0])
    newsf.sens =  None
    newsf.splice_multi_det = False
    newsf.wave_splice = None
    newsf.zeropoint_splice = None

    newsens = SensFunc.empty_sensfunc_table(1, len(combined_wave))
    newsens['SENS_WAVE'] = combined_wave
    newsens['SENS_ZEROPOINT_FIT'] = combined_zp_fit
    newsens['SENS_ZEROPOINT_FIT_GPM'] = combined_zp_fit_gpm
    newsens['WAVE_MIN'] = np.min(combined_wave)
    newsens['WAVE_MAX'] = np.max(combined_wave)

    newsf.sens = newsens

    newsf.spec1df = None
    newsf.std_name = None
    newsf.std_cal = None
    newsf.std_ra = None
    newsf.std_dec = None
    newsf.airmass = None
    newsf.telluric = None
    newsf.wave = np.empty((combined_wave.size,1))
    newsf.wave[:,0] = combined_wave
    newsf.zeropoint = np.empty((combined_zp_fit.size,1))
    newsf.zeropoint[:,0] = combined_zp_fit

    newsf.spectrograph = load_spectrograph("keck_deimos")    
    newsf.throughput = newsf.compute_throughput()[0]

    file_name = os.path.join(args.dest_path, f"keck_deimos_{args.grating}_sensfunc.fits")

    hdr = pypeit.io.initialize_header(primary=True)
    hdr['HISTORY'] = "This DEIMOS sensfunc was stitched together from the following source files."
    # Trim off spec1d_ from base of each spec1d file
    hdr['HISTORY'] = os.path.basename(sflist[0].spec1df)[7:]
    hdr['HISTORY'] = os.path.basename(sflist[1].spec1df)[7:]
    hdr['HISTORY'] = os.path.basename(sflist[2].spec1df)[7:]
    
    newsf.to_file(file_name, primary_hdr = hdr, overwrite=True)
    return (newsf, file_name)

def build_figure():

    utils.pyplot_rcparams()
    fig = plot.figure(figsize=(12,8))
    axis = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    return fig, axis

def create_plot(axis, color, label, x, y, linewidth=2.5, marker = '', linestyle='solid'):

    axis.plot(x, y, color=color, linewidth=linewidth, marker=marker, linestyle=linestyle, label=label)
    xmin = (0.98*x.min())
    #xmin = 3000
    xmax = (1.02*x.max())
    #ymin = 14.0
    ymin = 16.0
    ymax = (1.05*y.max())
    #ymax = 22

    return xmin, xmax, ymin, ymax


def setup_axes(axis,xmin, ymin, xmax, ymax, title):
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    axis.legend()
    axis.set_xlabel('Wavelength (Angstroms)')
    axis.set_ylabel('Zeropoint Fit (AB mag)')
    axis.set_title('PypeIt stitched SensFunc ' + title)

def plot_stitch_results(sf, polyfit_areas, sflist, filename, showQA):

    fig, axis = build_figure()

    sens_gpm = sf.sens['SENS_ZEROPOINT_FIT_GPM'][0]
    bpm = np.logical_not(sens_gpm)

    num_plots = 1

    if np.sum(bpm) != 0: # Only plot bad pixels if there are some
        num_plots += 1
        
    if polyfit_areas is not None:
        num_plots +=1

    if sflist is not None:
        num_plots += (2*len(sflist))
    xmin = np.zeros(num_plots)
    ymin = np.zeros(num_plots)
    xmax = np.zeros(num_plots)
    ymax = np.zeros(num_plots)

    i = 0

    # Plot the translated/stitched zero points without the "bad pixels" marked when generating sensfuncs
    x = sf.wave[sens_gpm]
    y = sf.zeropoint[sens_gpm]
    (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, (1, .7, .7, .5), f'Translated/stitched SENS_ZEROPOINT_FIT', x,y, marker='.', linestyle='none')
    i+=1

    # Plot the bad pixels
    if np.sum(bpm) != 0: # Only plot bad pixels if there are some
        x = sf.wave[bpm]
        y = sf.zeropoint[bpm]
        (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, "red", f'Inverse SENS_ZEROPOINT_FIT_GPM', x, y, marker='.', linestyle='none')
        i+=1

    # Plot the areas of the zeropoint that came from polynomial fit in blue
    if polyfit_areas is not None:
        x = sf.wave[polyfit_areas]
        y = sf.zeropoint[polyfit_areas]
        (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, "blue", f'Polynomial fit', x, y, marker='.', linestyle='none')
        i+=1

    # Plot the original sensfuncs in light gray in the background
    for bksf in sflist:
        for det in [0, 1]:
            gpm = bksf.sens['SENS_ZEROPOINT_FIT_GPM'][det]
            x = bksf.sens['SENS_WAVE'][det][gpm]
            y = bksf.sens['SENS_ZEROPOINT_FIT'][det][gpm]
            if i == 3:
                label = "Original sensfuncs"
            else:
                label = None
            (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, (.5, .5, .5), label, x,y, linewidth = 1, linestyle='solid')#linestyle='dashed')
            i+=1


    setup_axes(axis, np.min(xmin), np.min(ymin), np.max(xmax),np.max(ymax), filename)


    if showQA:
        plot.show()

    pdf_filename = os.path.splitext(filename)[0] + ".pdf"
    with PdfPages(pdf_filename) as pdf:
        pdf.savefig(fig)

def get_basename(source_file, source_meta):
    obsdate = datetime.datetime.strptime(source_meta['DATE'][0], "%d%b%Y")

    # Strip off path and .fits extension from the source name
    source_basename = os.path.splitext(os.path.basename(source_file))[0]

    # Build a new base name of the format "source-std_instrument_date"
    # For example: 2005aug27_d0827_0046-Feige110_DEIMOS_20050827
    return f"{source_basename}-{source_meta['STD_NAME'][0].replace(' ', '')}_{source_meta['INSTRUMENT'][0]}_{obsdate.strftime('%Y%m%d')}"

def get_source_meta(source_file):
    hdul = fits.open(source_file)
    return Table(hdul[1].data)

def create_spec1d_files(args, source_files):

    source_meta = []
    spec1d_files = []

    for source_file in source_files:
        meta = get_source_meta(source_file)
        spec1d_name = f"spec1d_{get_basename(source_file, meta)}.fits"
        dest_file = os.path.join(args.dest_path, spec1d_name)
        load_wmko_std_spectrum(source_file, str(dest_file), True)
        source_meta.append(meta)
        spec1d_files.append(dest_file)

    return source_meta, spec1d_files

def create_sens_files(args, spec1d_files):

    sflist = []


    spectrograph = load_spectrograph("keck_deimos")
    par = spectrograph.default_pypeit_par()
    par['sensfunc']['multi_spec_det'] = [3,7]
    par['sensfunc']['algorithm'] = "IR"
    par['sensfunc']['extrap_blu'] = 0.
    par['sensfunc']['extrap_red'] = 0.

    par_outfile = os.path.join(args.dest_path, "sensfunc.par")
    print(f'Writing the sensfunc parameters to {par_outfile}')
    par['sensfunc'].to_config(par_outfile, section_name='sensfunc', include_descr=False)


    for spec1d_file in spec1d_files:
        dest_file = os.path.join(args.dest_path, "sens_" + os.path.basename(spec1d_file)[7:])

        # Instantiate the relevant class for the requested algorithm
        sensobj = SensFunc.get_instance(spec1d_file, dest_file, par['sensfunc'])

        # Generate the sensfunc
        sensobj.run()

        # Write it out to a file
        sensobj.to_file(dest_file, overwrite=True)

        sflist.append(sensobj)




    return sflist 

def create_stitched_sensfuncs(args, sflist):

    (combined_wave,  combined_zp_fit, combined_zp_fit_gpm, polyfit_areas) = stitch_sensfunc(args.grating, sflist)

    (newsf, newfile) = write_stitched_sensfunc(sflist, args, combined_wave, combined_zp_fit, combined_zp_fit_gpm)

    plot_stitch_results(newsf, polyfit_areas, sflist, newfile, args.showQA)

    return newfile



def parse_args(options=None, return_parser=False):
    parser = argparse.ArgumentParser(description='Combine DEIMOS sensitivity functions to create a general purpose one.')
    parser.add_argument("grating", type=str, choices=['1200G', '1200B', '600ZD', '830G', '900ZD'])
    parser.add_argument("source_path", type=str, help="Path of the throughput fits files generated by Greg Writh's IDL script")
    parser.add_argument("dest_path", type=str, help = "Path to place generated spec1ds, sensfunc files, and the final stitched sensfunc.")
    parser.add_argument("--showQA", action="store_true", default=False, help="Show the the QA plots before saving the QA pdf.")
            

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """ Executes sensitivity function computation.
    """

    source_map = {"1200G": ['extract/1200G/2005aug27_d0827_0046.fits',
                            'extract/1200G/2005aug27_d0827_0047.fits',
                            'extract/1200G/2005aug27_d0827_0048.fits'],
                  "1200B": ['extract/1200B/2017oct03_d1003_0055.fits',
                            'extract/1200B/2017oct03_d1003_0056.fits',
                            'extract/1200B/2017oct03_d1003_0057.fits'],
                  "900ZD": ['extract/900ZD/2010oct06_d1006_0138.fits',
                            'extract/900ZD/2010oct06_d1006_0139.fits',
                            'extract/900ZD/2010oct06_d1006_0140.fits'],
                  "600ZD": ['extract/600ZD/2010sep24_d0924_0008.fits',
                            'extract/600ZD/2010sep24_d0924_0009.fits',
                            'extract/600ZD/2010sep24_d0924_0010.fits'],
                  "830G":  ['extract/830G/2010oct06_d1006_0135.fits',
                            'extract/830G/2010oct06_d1006_0136.fits',
                            'extract/830G/2010oct06_d1006_0137.fits']}

    source_files = [os.path.join(args.source_path, x) for x in source_map[args.grating]]

    source_meta, spec1d_files = create_spec1d_files(args, source_files)

    sflist = create_sens_files(args, spec1d_files)

    create_stitched_sensfuncs(args, sflist)



def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
