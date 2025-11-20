
import os
import numpy as np
from pypeit import specobjs
from pypeit import flatfield
from pypeit.core import extract
from pypeit.utils import inverse
from pypeit.core.moment import moment1d
from pypeit.core import flat
from pypeit.core import coadd
from matplotlib import pyplot as plt

redux_path = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/shane_kast_red/600_7500_d55_ret'
spec1d_file =  os.path.join(redux_path, 'Science/spec1d_r136-G191b2b_KASTr_20150123T024320.750.fits')
flat_file = os.path.join(redux_path, 'Calibrations/Flat_A_0_DET01.fits')
#redux_path = '/Users/joe/python/PypeIt-development-suite//REDUX_OUT/shane_kast_red/600_7500_d57'
#spec1d_file = os.path.join(redux_path, 'Science/spec1d_r121-bd284211_KASTr_20181014T034326.900.fits')
#flat_file = os.path.join(redux_path, 'Calibrations/Flat_A_0_DET01.fits')

sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file, chk_version=False)
sobj_std = sobjs.get_std()[0]
sobj_flat = sobj_std.copy()
sobj_std.BOX_RADIUS = 10.0
flatImages = flatfield.FlatImages.from_file(flat_file)

pixelflat_raw = flatImages.pixelflat_raw
pixelflat_norm = flatImages.pixelflat_norm
pixelflat_proc, flat_bpm = flat.flatfield(pixelflat_raw, pixelflat_norm)


inmask = pixelflat_raw > 0.0
invvar = inverse(pixelflat_raw)

extract.extract_boxcar(pixelflat_proc, invvar, np.logical_not(flat_bpm), 1000*np.ones_like(pixelflat_raw),
                       0.0*pixelflat_raw, sobj_flat, base_var=None,
                       count_scale=None, noise_floor=None)

flux_new, ivar_new, gpm_new, blaze_new = coadd.interp_oned(sobj_std.BOX_WAVE, sobj_std.BOX_WAVE, sobj_std.BOX_COUNTS,
                                                           sobj_std.BOX_COUNTS_IVAR, sobj_std.BOX_MASK, blaze_function=sobj_flat.BOX_COUNTS,
                                                           sensfunc=True)

blaze_norm = blaze_new/blaze_new.max()
std_norm = sobj_std.BOX_COUNTS/sobj_std.BOX_COUNTS.max()
std_by_blaze = sobj_std.BOX_COUNTS*inverse(blaze_norm)
std_by_blaze_norm = std_by_blaze/std_by_blaze.max()
plt.plot(sobj_std.BOX_WAVE,blaze_norm, color='b', label='blaze')
plt.plot(sobj_std.BOX_WAVE,std_norm, color='r', label='std')
plt.plot(sobj_std.BOX_WAVE,std_by_blaze_norm, color='g', label='std/blaze')
#plt.plot(sobj_std.BOX_WAVE,sobj_std.BOX_COUNTS/sobj_std.BOX_COUNTS.max())
plt.legend()
plt.show()