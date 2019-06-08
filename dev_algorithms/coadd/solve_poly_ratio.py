

import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
from pypeit import utils
from pypeit.core import coadd
from scipy import interpolate
from coadd1d import interp_spec
from scipy import stats
from pypeit import msgs
import IPython
from coadd1d import solve_poly_ratio

#def solve_poly_fn(theta, xvector, polyfunc, nback = None):
#
#    if nback is None:
#        acoeff = theta
#    else:
#        acoeff = theta[0:-nback]
#        bcoeff = theta[-nback:]
#
#    ymult = utils.func_val(acoeff, xvector, polyfunc, minx=wave_min, maxx=wave_max)
#    if nback is not None:
#        ymult = utils.func_val(acoeff, xvector, polyfunc, minx=wave_min, maxx=wave_max)



from astropy.table import Table
from astropy.io import fits
from matplotlib import pyplot as plt
import sys

fitfunc = poly_ratio_fitfunc
stackfile  = '/Users/joe/REDUX/lris_redux/Nov_2004/Final/SDSSJ073522.43+295710.1_N.fits'
hdu = fits.open(stackfile)
flux_ref = hdu[0].data
sig_ref = hdu[1].data
ivar_ref = utils.calc_ivar(sig_ref)
mask_ref = (ivar_ref > 0.0)
wave_ref = hdu[2].data

infiles = ['/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0063.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0064.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0065.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0066.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0219.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0220.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0221.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0222.fits.gz',
    '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0063.fits.gz']
nfiles = len(infiles)
objid = 0
for idx, file in enumerate(infiles):
    obj = Table.read(file, hdu=5)
    flux = obj[objid]['FLUX_OPT']
    ivar = obj[objid]['IVAR_OPT']
    wave = obj[objid]['WAVE_OPT']
    if idx == 0:
        nspec = flux.size
        flux_arr = np.zeros((nfiles, nspec))
        wave_arr = np.zeros((nfiles, nspec))
        ivar_arr = np.zeros((nfiles, nspec))
    flux_arr[idx,:] = flux
    ivar_arr[idx,:] = ivar
    wave_arr[idx,:] = wave

mask_arr = (ivar_arr > 0.0)
#plt.plot(wave_ref, flux_ref, drawstyle='steps-mid')
#for ii in range(nfiles):
#    plt.plot(wave_arr[:,ii], flux_arr[:,ii],drawstyle='steps-mid')

flux_inter, ivar_inter, mask_inter = interp_spec(wave_ref, wave_arr, flux_arr, ivar_arr, mask_arr)

idx = 5
wave = wave_ref
flux = flux_inter[idx, :]
ivar = ivar_inter[idx, :]
mask = ivar > 0
norder = 3

ymult,flux_rescale, ivar_rescale, outmask = solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, norder,
                                                             mask=mask, mask_ref=mask_ref, debug=True)