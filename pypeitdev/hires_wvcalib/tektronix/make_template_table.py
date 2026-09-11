"""
The template table currently available does not correctly label the orders. This script
generates a new template table with the correct order labels.
"""
import os
import glob
from IPython import embed

import numpy as np
from astropy.table import Table
import astropy.io.fits as fits

import pypeit
pypeit_path = os.path.dirname(os.path.realpath(pypeit.__file__))

def empty_design_table(nrows):
    """
    Construct an empty arxiv table.

    Args:
        nrows (:obj:`int`):
            Number of rows in the table.

    Returns:
        `astropy.table.Table`_: Instance of the empty arxiv table.
    """
    return Table([np.zeros(nrows, dtype="<U30"),
                        np.zeros(nrows, dtype=int),
                        np.zeros(nrows, dtype="<U6"),
                        np.zeros(nrows, dtype=float),
                        np.zeros(nrows, dtype=float),
                        np.zeros(nrows, dtype="<U3"),
                        np.zeros(nrows, dtype=int),
                        np.zeros(nrows, dtype=int),
                        np.zeros(nrows, dtype=int),
                        np.zeros(nrows, dtype=int)],
                        names=('Name', 'Chip', 'XDISP', 'ECH', 'XDAng', 'Deck', 'Rbin', 'Cbin', 'IOrder', 'EOrder'))


if __name__ == '__main__':
    # Allow for command line arguments to specify an overwrite option
    import argparse
    parser = argparse.ArgumentParser(description='Ingest HIRES (original detector) makee wavelength calibrations and perform fits to the wavelength solution coefficients vs. XDISP and ECH angle')
    parser.add_argument('-o', '--overwrite', action='store_true', help='Whether to overwrite existing files')
    args = parser.parse_args()
    overwrite = args.overwrite

    # Define the input and output files
    intempl = "hires_orig_templ_input.dat"
    outfile = "hires_orig_templ_makee.dat"

    # Grad a list of all files
    templ_dirc = "templates_makee/"
    files = sorted(glob.glob(f"{templ_dirc}*.fits"))
    outtable = empty_design_table(len(files))

    # Read in a template spectrum from the new HIRES to use to align the orders.
    hires_file = os.path.join(pypeit_path, 'data', 'arc_lines', 'reid_arxiv',
                             'keck_hires_composite_arc.fits')
    print("Loading template file: ", hires_file)
    hdu = fits.open(hires_file, chk_version=False)
    orders = np.arange(hdu[1].data['order_min'][0], hdu[1].data['order_max'][0]+1)
    assert(orders.size == hdu[1].data['norders'][0])
    # Grab a wavelength array and mask out the zero values
    mskarr = np.ma.array(hdu[2].data, mask=(hdu[2].data == 0))
    mid_waves = np.ma.mean(mskarr, axis=0).data
    # Make a fitting function to get the order number as a function of wavelength
    order_fit = np.polyfit(mid_waves, orders, 6)
    order_mod = np.polyval(order_fit, mid_waves)
    # plt.plot(mid_waves, orders-order_mod)
    # plt.show()

    # Loop through the template spectra and fill in the table with the information from the header
    for ii, file in enumerate(files):
        # Read in the template spectrum
        hdul = fits.open(file)
        header = hdul[0].header
        data = hdul[0].data

        # Get the XDISP value from the filename
        if "adrd97-" in file:
            xdisp = "RED97"
        elif "adrd-" in file:
            xdisp = "RED"
        elif "aduv-" in file:
            xdisp = "UV"
        else:
            print("File loading error: ", file)
            assert False

        # Load the wavelength solution
        all_wave, mid_wave, mid_ordr = [], [], []
        for oo in range(header['NAXIS2']):
            coeff = np.array(header['WV_0_{0:02d}'.format(oo+1)].split() + \
                     header['WV_4_{0:02d}'.format(oo+1)].split()).astype(float)[::-1]
            pixarr = np.arange(header['NAXIS1'])
            all_wave.append(np.polyval(coeff, pixarr))
            mid_wave.append(np.mean(all_wave[oo]))
            mid_ordr.append(int(np.round(np.polyval(order_fit, mid_wave[oo]))))
        mid_wave = np.array(mid_wave)
        mid_ordr = np.array(mid_ordr)

        # Fill in the table with the information from the header
        outtable[ii]['Name'] = file.split('/')[-1]
        outtable[ii]['Chip'] = 1
        outtable[ii]['XDISP'] = xdisp
        outtable[ii]['ECH'] = header['ECHANGL']
        outtable[ii]['XDAng'] = header['XDANGL']
        outtable[ii]['Deck'] = header['DECKNAME'].strip()
        outtable[ii]['Rbin'] = header['BINNING'].strip().split(",")[0]
        outtable[ii]['Cbin'] = header['BINNING'].strip().split(",")[1]
        outtable[ii]['IOrder'] = np.min(mid_ordr)
        outtable[ii]['EOrder'] = np.max(mid_ordr)

    # Write the table to a file
    print("Writing output table to: ", outfile)
    outtable.write(outfile, format='ascii.fixed_width', overwrite=overwrite)