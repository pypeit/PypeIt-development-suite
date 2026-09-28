""" Codes to help port ESI wavelengths from XIDL """
import os
from astropy.table import Table

import numpy as np

from scipy.io import readsav
from matplotlib import pyplot as plt

from pypeit.core.wavecal import templates

from IPython import embed

def ingest_xidl_archive(outfile, debug=False):

    # Read archive file
    #tbl = Table.read('esi_archive_orders.fits')

    # Read SAV file
    xidl_dict = readsav('CALIBS/ECH_arcfit.idl')
    calib = xidl_dict['all_arcfit']
    

    # Orders
    order_vec = np.arange(15, 6 - 1, -1)
    order_vec_raw, wave, arc = templates.xidl_esihires(
        'CALIBS/ECH_arcfit.idl', order_vec=order_vec,
        log10=False)

    # Build a Table
    tbl = Table()
    tbl['order'] = order_vec
    tbl['wave'] = wave
    tbl['flux'] = arc

    # Write
    tbl.write(outfile, overwrite=True)
    print(f'Wrote: {outfile}')

    # Check
    if debug:
        plt.clf()
        plt.plot(wave[0,:], arc[0,:])
        plt.show()

# Command line execution
if __name__ == '__main__':
    ingest_xidl_archive('keck_esi_ECH.fits')