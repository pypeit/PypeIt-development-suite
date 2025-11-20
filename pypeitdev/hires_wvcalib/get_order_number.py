


import os
import sys
import numpy as np
#import jax
#from jax import numpy as jnp
#from jax import jit
#import optax
import itertools
from tqdm.auto import trange
from astropy.table import Table
from astropy.io import fits
from pkg_resources import resource_filename
from matplotlib import pyplot as plt
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.wavecal import templates
from pypeit.core.wavecal import wvutils
from pypeit.core.fitting import robust_fit
from pypeit.core import coadd
from pypeit.core import fitting
from pypeit.core.wavecal import autoid, waveio, wv_fitting
from pypeit.core.wavecal.wvutils import  smooth_ceil_cont, xcorr_shift
from pypeit import utils
from pypeit import msgs
from astropy import table
from scipy import interpolate
from IPython import embed
from astropy import constants as const

c_kms = const.c.to('km/s').value



xidl_arxiv_file = os.path.join(os.getenv('PYPEIT_DEV'), 'dev_algorithms', 'hires_wvcalib', 'hires_wvcalib_xidl.fits')

arxiv_params = Table.read(xidl_arxiv_file, hdu=1)[0]
arxiv = Table.read(xidl_arxiv_file, hdu=2)
order_vec = np.arange(arxiv_params['order_min'], arxiv_params['order_max'] + 1, 1)

# Plot of coefficient vs ech_angle
#iord =35
#UV  = arxiv['populated_and_good'][:, iord] & (arxiv['xdisp'][:, iord] == 'UV')
#RED = arxiv['populated_and_good'][:, iord] & (arxiv['xdisp'][:, iord] == 'RED')
#for icoeff in range(arxiv['coeff'].shape[-1]):
#    plt.plot(arxiv['ech_angle'][UV, iord], arxiv['coeff'][UV, iord, icoeff], 'b.', label='UV')
#    plt.plot(arxiv['ech_angle'][RED, iord], arxiv['coeff'][RED, iord, icoeff], 'r+', label='RED')
#    plt.title('Order = {:d} and Coeff = {:d}'.format(order_vec[iord], icoeff))
#    plt.show()


# Plot of bluest order vs xd_angle
for det in range(1, 4, 1):
    UV  = (arxiv['det_file'] == det) & (arxiv['xdisp_file'] == 'UV')
    RED = (arxiv['det_file'] == det)  & (arxiv['xdisp_file'] == 'RED')
    plt.plot(arxiv['xd_angle_file'][UV], arxiv['reddest_order'][UV], 'b.', label='UV, bluest')
    plt.plot(arxiv['xd_angle_file'][UV], arxiv['reddest_good_order'][UV], 'b+', label='UV, bluest good')
    plt.plot(arxiv['xd_angle_file'][RED], arxiv['bluest_order'][RED], 'r.', label='RED, bluest')
    plt.plot(arxiv['xd_angle_file'][RED], arxiv['bluest_good_order'][RED], 'r+', label='RED, bluest good')
    plt.legend()
    plt.show()
