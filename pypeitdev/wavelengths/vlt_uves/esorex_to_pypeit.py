import numpy as np
from scipy.io import readsav
import glob
from IPython import embed

import astropy.io.fits as fits
from astropy.table import Table
from matplotlib import pyplot as plt

from pypeit.core.wavecal import wvutils
from pypeit.core.arc import fit2darc


def convert_esorex_to_reid(fname, outname, to_cache=False):
    """
    Convert the ESOREX outputs to a FITS file for use in PypeIt.
    """
    # Load the data from the .sav file
    allfils = glob.glob("{0:s}/thar_*_sci_{0:s}_??_?.fits".format(fname))
    if len(allfils) == 0:
        raise IOError("No ThAr files found in {0:s}".format(fname))
    elif len(allfils) > 1:
        raise IOError("Multiple ThAr files found in {0:s}".format(fname))
    # Read the ThAr file
    tharname = allfils[0]
    print("Reading ThAr file: {0:s}".format(tharname))
    specdata = fits.open(tharname)[0].data
    nord, nspec = specdata.shape

    wallfils = glob.glob("{0:s}/wpol_*_sci_{0:s}_??_?.fits".format(fname))
    if len(wallfils) == 0:
        raise IOError("No ThAr files found in {0:s}".format(fname))
    elif len(wallfils) > 1:
        raise IOError("Multiple ThAr files found in {0:s}".format(fname))
    # Read the wavelength fits file
    wpolname = wallfils[0]
    print("Reading WPOL file: {0:s}".format(wpolname))
    hdu = fits.open(wpolname)
    all_pix = hdu[1].data['X']
    all_orders = hdu[1].data['Order']
    all_wv = hdu[1].data['Ident']
    wgd = np.where(np.logical_not(np.isnan(all_wv)))
    orderref = np.min(all_orders)
    if False:
        test_order = 130
        wgd = np.where(np.logical_not(np.isnan(all_wv)) & (all_orders==test_order))
        all_peak = hdu[1].data['Peak']
        wave = np.arange(nspec)
        spec = specdata[orderref-test_order-1,:]
        plt.plot(wave, spec, 'k-', drawstyle='steps-mid', lw=1)
        plt.plot(all_pix[wgd], all_peak[wgd], 'ro', markersize=1)
        plt.show()

    print("Manually update the order list for each setting here")
    if fname == "346":
        nspec_coeff = 6
        norder_coeff = 4
        order_list = (orderref + np.arange(nord))[::-1]
    elif fname == "390":
        nspec_coeff = 6
        norder_coeff = 4
        order_list = (orderref + np.arange(nord))[::-1]
    else:
        print("Unknown order list")
        embed()
        assert False

    # Fit the 2D arc
    arcfit = fit2darc(all_wv[wgd], all_pix[wgd], all_orders[wgd], nspec, nspec_coeff=nspec_coeff, norder_coeff=norder_coeff, debug=True)

    # Use this 2D fit to get the wavelength solution for each order
    pixlist = np.arange(nspec)[np.newaxis, :].repeat(nord, axis=0)/(nspec-1)
    ordlist = order_list[:, np.newaxis].repeat(nspec, axis=1)
    wavelengths = arcfit.eval(pixlist, x2=ordlist)
    wavelengths /= ordlist

    if True:
        for oo in range(nord):
            plt.plot(wavelengths[oo,:], specdata[oo,:])
        plt.xlim(np.median(wavelengths), np.median(wavelengths)+50)
        plt.ylim(-100, 5000)
        plt.show()

    # Check the binning
    # BLUE XDs
    if np.abs(nspec-3000) < 3000/100:
        binspec = 1
        print("npix=", nspec, " --> BLUE XD --> Binning = 1")
    elif np.abs(nspec-1500) < 1500/100:
        binspec = 2
        print("npix=", nspec, " --> BLUE XD --> Binning = 2")
    # RED XDs
    elif np.abs(nspec - 4060) < 1:#4060 / 100:
        binspec = 1
        print("npix=", nspec, " --> RED XD --> Binning = 1")
    elif np.abs(nspec - 2030) < 1:#2030 / 100:
        binspec = 1
        print("npix=", nspec, " --> RED XD --> Binning = 2")
    else:
        print("Error :: Binning unknown", nspec)
        assert(False)

    # Write out the data to a FITS file
    wvutils.write_template(wavelengths, specdata, binspec, './', outname, to_cache=to_cache, order=order_list)


if __name__ == "__main__":
    # List of all files to convert
    fils = ["346", "390"]#, "760"]
    outnames = ["vlt_uves_346_1x1.fits", "vlt_uves_390_1x1.fits"]
    for ff, fil in enumerate(fils):
        # Convert the ESOREX reduction files to a FITS file
        convert_esorex_to_reid(fil, outnames[ff], to_cache=False)
