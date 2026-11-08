import numpy as np
from scipy.io import readsav

import astropy.io.fits as fits
from astropy.table import Table

from pypeit.core.wavecal import wvutils


def convert_idl_to_reid(fname, outname, to_cache=False):
    """
    Convert the IDL .sav file to a FITS file for use in PypeIt.
    """
    # Load the data from the .sav file
    data = readsav(fname)
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
        # Normalize the xdat to the range of the fit
        xnrm = 2. * (xdat - this_nrm[0]) / this_nrm[1]
        coeffs = data['all_arcfit']['ffit'][oo].flatten()[::-1]
        wave = 10.0**np.polyval(coeffs, xnrm)
        if wave[1]-wave[0] < 0:
            # Reverse the direction of the wavelength solution
            wsave = wave.copy()[::-1]
            fsave = data['sv_aspec'][oo, :][::-1]
        else:
            wsave = wave.copy()
            fsave = data['sv_aspec'][oo, :]
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


if __name__ == "__main__":
    # List of all files to convert
    fils = ["uves_arc_blue_390.idl", "uves_arc_red_580.idl"]
    outnames = ["vlt_uves_390_1x1.fits", "vlt_uves_580_1x1.fits"]
    for ff, fil in enumerate(fils):
        # Convert the IDL file to a FITS file
        convert_idl_to_reid(fil, outnames[ff], to_cache=True)
