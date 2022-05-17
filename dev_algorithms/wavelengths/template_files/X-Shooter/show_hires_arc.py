''' Simple script to show a HIRES arc'''

import os
import numpy as np
from scipy.io import readsav

from matplotlib import pyplot as plt

from astropy import units

from pypeit.core.wave import airtovac
from pypeit.core.wavecal import templates
from pypeit.core import arc

from IPython import embed


def xidl_hires(xidl_file, specbin=1):
    """
    Read an XIDL format solution for Keck/HIRES
    Note:  They used air
    Args:
        xidl_file (str):
            Keck/HIRES save file
    Returns:
    """
    xidl_dict = readsav(xidl_file)
    order_vec = xidl_dict['guess_ordr']
    norders = order_vec.size
    nspec = xidl_dict['sv_aspec'].shape[1]

    # Wavelengths
    wave = np.zeros((norders, specbin*nspec))
    spec = np.zeros((norders, specbin*nspec))

    calib = xidl_dict['all_arcfit']
    order_mask = np.ones(norders, dtype=bool)

    # Here we go on the fits
    for kk in range(norders):
        # Generate the wavelengths
        if calib['FUNC'][kk] == b'CHEBY':
            log10_wv_air = cheby_val(calib['FFIT'][kk], 
                               np.arange(nspec),
                        calib['NRM'][kk], calib['NORD'][kk])
        elif calib['FUNC'][kk] == b'POLY':
            log10_wv_air = templates.poly_val(calib['FFIT'][kk], 
                              np.arange(nspec),
                              calib['NRM'][kk])
        else:
            order_mask[kk]=False
            continue

        wv_vac = airtovac(10**log10_wv_air * units.AA).value
        ispec = xidl_dict['sv_aspec'][kk,:]
        # Flip to blue to red?
        if wv_vac[1] < wv_vac[0]:
            wv_vac = wv_vac[::-1]
            ispec = ispec[::-1]
        # Fill
        if specbin != 1:
            wave[kk,:] = arc.resize_spec(wv_vac, nspec*specbin)
            spec[kk,:] = arc.resize_spec(ispec, nspec*specbin)
        else:
            wave[kk,:] = wv_vac
            spec[kk,:] = ispec
    # Return

    return order_vec[order_mask], wave[order_mask,:], spec[order_mask,:]

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to fuss with FRB galaxies [v1.1]')
    parser.add_argument("xidl_file", type=str, help="XIDL filename")
    parser.add_argument("order", type=int, help="Order number")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs

def main(pargs):
    # Load up the spectrum
    full_path = os.path.join(os.getenv('XIDL'), 'Keck', 'HIRES', 'CALIBS',
        'ARCS', pargs.xidl_file)
    orders, wave, spec = xidl_hires(full_path)

    this_order = np.where(orders == pargs.order)[0][0]

    # Plot
    plt.clf()
    ax = plt.gca()
    ax.plot(wave[this_order], spec[this_order], 'k')
    plt.show()

if __name__ == "__main__":

    # get the arguments
    pargs = parser()

    # Run
    main(pargs)


# Examples
# python show_hires_arc.py hires_tmpl2x1B0.idl 107