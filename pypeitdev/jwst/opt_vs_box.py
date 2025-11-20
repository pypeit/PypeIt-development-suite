import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

sci_path = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/pypeit/Science'
box_file = os.path.join(sci_path, 'spec1d_diff_jw02736007001_03103_00001-00003_slit64.fits')
opt_file = os.path.join(sci_path, 'spec1d_jw02736007001_03103_00001-00003_slit64.fits')
box_tbl = Table.read(box_file, hdu=1)
opt_tbl = Table.read(opt_file, hdu=1)


fx = plt.figure(1, figsize=(12, 6))
# left, bottom, width, height
rect = [0.16, 0.14, 0.82, 0.73]
ax = fx.add_axes(rect)

box_wv_mask = box_tbl['BOX_WAVE'] > 1.0
opt_wv_mask = opt_tbl['OPT_WAVE'] > 1.0
ax.plot(box_tbl['BOX_WAVE'][box_wv_mask], box_tbl['BOX_COUNTS']*box_tbl['BOX_MASK'][box_wv_mask],
        drawstyle='steps-mid', color='black', label='diff-boxcar')
ax.plot(box_tbl['BOX_WAVE'][box_wv_mask], box_tbl['BOX_COUNTS_SIG']*box_tbl['BOX_MASK'][box_wv_mask],
        drawstyle='steps-mid', color='red', label='diff-boxcar-error')
ax.plot(opt_tbl['OPT_WAVE'][opt_wv_mask], opt_tbl['OPT_COUNTS']*opt_tbl['OPT_MASK'][opt_wv_mask]
        , drawstyle='steps-mid', color='green', label='optimal')
ax.plot(opt_tbl['OPT_WAVE'][opt_wv_mask], opt_tbl['OPT_COUNTS_SIG']*opt_tbl['OPT_MASK'][opt_wv_mask]
        , drawstyle='steps-mid', color='purple', label='optimal-error')
ax.set_xlabel('Wavelength (Angstroms)')
ax.set_ylabel('Counts')
plt.legend()
plt.show()
