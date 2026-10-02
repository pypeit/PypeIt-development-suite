import glob

import numpy as np
import astropy.io.fits as fits
from matplotlib import pyplot as plt

from pypeit.core.wavecal import waveio


def check_orders(infile, ordcheck):
    """
    Check that the orders are correct in the files output by esorex_to_pypeit.py
    """
    info, par = waveio.load_reid_arxiv(infile)
    # Find the correct order
    foundit = False
    for ii in range(len(info)):
        if info[str(ii)]['order'] == ordcheck:
            flux = info[str(ii)]['spec']
            wave = info[str(ii)]['wave_soln']
            foundit = True
    if not foundit:
        print(f"WARNING!!  Did not find order {ordcheck} in {infile}")
        return
    # Check one order (and assume the rest are correct)
    plt.plot(wave, flux, 'k-', drawstyle='steps-mid')
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Flux')
    plt.title(f'Order {ordcheck} for {infile}')
    plt.show()


if __name__ == "__main__":
    # List of all files to check
    ordcheck = [120, 100,
                120, 100,
                100, 80,
                80, 70]
    outnames = ["vlt_uves_564l_1x1.fits", "vlt_uves_564u_1x1.fits",
                "vlt_uves_580l_1x1.fits", "vlt_uves_580u_1x1.fits",
                "vlt_uves_760l_1x1.fits", "vlt_uves_760u_1x1.fits",
                "vlt_uves_860l_1x1.fits", "vlt_uves_860u_1x1.fits"]
    for ff, fil in enumerate(outnames):
        check_orders(fil, ordcheck[ff])
