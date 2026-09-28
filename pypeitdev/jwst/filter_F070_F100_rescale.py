import os
from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy import interpolate


filterpath = '/Users/joe/python/PypeIt-development-suite/dev_algorithms/jwst/filters/nirspec/'
F070LP_file = os.path.join(filterpath, 'jwst_nirspec_f070lp_trans.fits')
F100LP_file = os.path.join(filterpath, 'jwst_nirspec_f100lp_trans.fits')

F070LP = Table.read(F070LP_file)
F100LP = Table.read(F100LP_file)

_wave_F070LP = F070LP['WAVELENGTH']
_trans_F070LP = F070LP['THROUGHPUT']
wave_F100LP = F100LP['WAVELENGTH']
trans_F100LP = F100LP['THROUGHPUT']

# Interpolate F070LP transmission curve to F100LP wavelength grid
wave_mask  = wave_F100LP > 0.9600
wave_trans = wave_F100LP[wave_mask]
trans_F070LP = interpolate.interp1d(_wave_F070LP, _trans_F070LP, bounds_error=False, fill_value=np.nan)(wave_trans)
trans_ratio = trans_F070LP/trans_F100LP[wave_mask]

# Save the transmission ratio to a file
outfile = os.path.join(filterpath, 'F070overF100_trans_ratio.fits')
# Create a new astropy table to hold output
out_table = Table()
out_table['wavelength'] = wave_trans
out_table['F070overF100_trans_ratio'] = trans_ratio
out_table.write(outfile, overwrite=True)

show=True
if show:
    plt.plot(wave_trans, trans_ratio)
    plt.show()
