import numpy as np
from scipy.io import readsav

import astropy.io.fits as fits
from astropy.table import Table

from pypeit.core.wavecal import wvutils


def convert_idl_to_reid(fname, outname, to_cache=False):
    """
    Convert the IDL .sav file to a FITS file for use in PypeIt.
    """

    order = data['guess_ordr']
    nord = order.size
    npix = data['sv_aspec'].shape[1]
    xdat = np.arange(npix)
    nrm = data['all_arcfit']['nrm']

    # Check the binning
    # BLUE XDs
    if np.abs(npix-2940) < 2940/100:
        binspec = 1
        print("npix=", npix, " --> Binning = 1")
    elif np.abs(npix-1470) < 1470/100:
        binspec = 2
        print("npix=", npix, " --> Binning = 2")
    # RED XDs
    elif np.abs(npix - 4060) < 4060 / 100:
        binspec = 1
        print("npix=", npix, " --> Binning = 1")
    elif np.abs(npix - 2030) < 2030 / 100:
        binspec = 1
        print("npix=", npix, " --> Binning = 2")
    else:
        print("Error :: Binning unknown", npix)
        assert(False)

    # Create a binary table HDU containing the spectrum and wavelength solution, one row at a time
    # The first column is the order number, the second column is the wavelength, and the third column is the flux
    orders = np.array([])
    wavelengths, _specdata = None, None
    for oo in range(nord):
        this_nrm = nrm[oo]
        # Check if the normalization is valid
        if this_nrm[1] == 0:
            # Bad fit... skip this order
            continue


if __name__ == "__main__":
    print("This is not the best way to go. Use the UVES pipeline results instead in esorex_to_pypeit.py")
    assert False
    # List of all files to convert
    fname = "arccen_346.npy"
    outname = "vlt_uves_346_1x1.fits"
    # Load the data from the .npy file that was output from within PypeIt.
    fsave = np.load(fname)

    # Load a comparison frame to build a partial solution
    hdu = fits.open("vlt_uves_390_1x1.fits")
    data = hdu[1].data
    binspec = hdu[1].header['BINSPEC']
    ref_wave, ref_flux, ref_order = data['wave'], data['flux'], data['order']

    # Compare the reference spectrum to each of the new spectra
    for oo in range(fsave.shape[1]):
        # Plot the spectrum
        plt.plot(fsave[:,oo])
        plt.plot(ref_flux[0,:], 'k-', lw=3, label='Ref')
        plt.show()


    # Store in the arrays:
    if wavelengths is None:
        wavelengths = wsave.reshape((1, -1))
        _specdata = fsave.reshape((1, -1))
    else:
        wavelengths = np.append(wavelengths, wsave.reshape((1, -1)), axis=0)
        _specdata = np.append(_specdata, fsave.reshape((1, -1)), axis=0)
    orders = np.append(orders, order[oo])

    # Write out the data to a FITS file
    wvutils.write_template(wavelengths, _specdata, binspec, './', outname, to_cache=to_cache, order=orders)
