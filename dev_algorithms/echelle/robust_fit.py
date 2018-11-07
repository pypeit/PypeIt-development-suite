
import numpy as np
from pypeit.core import extract
from pypeit import msgs
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

from pypeit import ginga
from pypeit.core import pydl
from pypeit import utils
from pypeit.core import pixels
from sklearn.decomposition import PCA
from pypeit import specobjs
from pypeit.core import extract
from astropy.stats import SigmaClip
from pydl.pydlutils.spheregroup import spheregroup
import matplotlib.pyplot as plt

### Test polynomial fit
from numpy.polynomial import polynomial as poly
xfit = np.linspace(-5,5,101)
rand = np.random.RandomState(3)
yfit = xfit**3 - xfit + rand.randn(len(xfit))*20.0
yreal = xfit**3 - xfit
nout =50
indrand = rand.choice(len(yfit),nout)
std = np.std(yreal-yfit)
nsigma = rand.uniform(2,5,nout)
sign = rand.choice([-1,1.],nout)

yfit[indrand] =  yreal[indrand] + sign*nsigma*std
#yfit[0] = yfit[0]+3*np.std(yreal-yfit)
#yfit[10] = yfit[10]+4*np.std(yreal-yfit)
#yfit[19] = yfit[19]-3.5*np.std(yreal-yfit)
#yfit[49] = yfit[49]-2.8*np.std(yreal-yfit)
#yfit[69] = yfit[69]+5*np.std(yreal-yfit)
norder =3

xvec = np.linspace(xfit.min(), xfit.max(), num=200)

## Old robust_olyfit
msk, poly_coeff = utils.robust_polyfit(xfit, yfit, norder, sigma=3.0, function='polynomial')

msk_new, poly_coeff_new = utils.robust_polyfit_djs(xfit,yfit,norder, \
                                           function = 'polynomial', minv = None, maxv = None, bspline_par = None,\
                                           guesses = None, maxiter = 10, inmask = None, sigma = None,invvar = None,\
                                           lower = 2, upper = 2,maxdev=None,maxrej=None,groupdim=None,groupsize=None,\
                                           groupbadpix=False, grow=0,sticky=True,use_mad=True)

msk_nosticky, poly_coeff_nosticky = utils.robust_polyfit_djs(xfit,yfit,norder, \
                                           function = 'polynomial', minv = None, maxv = None, bspline_par = None,\
                                           guesses = None, maxiter = 10, inmask = None, sigma = None,invvar = None,\
                                           lower = 2, upper = 2,maxdev=None,maxrej=None,groupdim=None,groupsize=None,\
                                           groupbadpix=False, grow=0,sticky=False,use_mad=True)
robust_mask = msk == 0
robust_mask_new = msk_new == 1
robust_mask_nosticky = msk_nosticky == 1
maskdiff = np.any(robust_mask_nosticky != robust_mask_new)
coeffdiff = np.any(poly_coeff_nosticky != poly_coeff_new)

print(maskdiff)

plt.figure()
plt.plot(xfit,yfit,'ko',mfc='None',label='Good Points')
#plt.plot(xfit[[0,10,19,49,69]],yfit[[0,10,19,49,69]],'bo',label='Outlier')
plt.plot(xfit[indrand],yfit[indrand],'bo',label='Outlier')
plt.plot(xfit,yreal,'k-',lw=3,label='Real')
plt.plot(xfit[~robust_mask], yfit[~robust_mask], 'ms', markersize=10.0,mfc='None', label='robust_polyfit rejected')
plt.plot(xfit[~robust_mask_new],yfit[~robust_mask_new],'r+', markersize = 20.0, label = 'robust_polyfit_djs rejected')
plt.plot(xfit[~robust_mask_nosticky],yfit[~robust_mask_nosticky],'go', mfc='None',markersize = 30.0, label = 'robust_polyfit_djs rejected')

plt.plot(xvec, utils.func_val(poly_coeff, xvec, 'polynomial'),lw=2,ls='-.', color='m', label='robust polyfit')
plt.plot(xvec, utils.func_val(poly_coeff_new, xvec, 'polynomial'),lw=1,ls='--', color='r',label = 'new robust polyfit')
plt.plot(xvec, utils.func_val(poly_coeff_nosticky, xvec, 'polynomial'),lw=1,ls=':', color='g',label = 'new robust polyfit')

plt.legend()
plt.show()

sys.exit(-1)

### Test some results from the PCA
#xfit = np.array([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
#yfit = np.array([ 5205.0605,3524.0981,1974.9368,694.22455,-359.67508,-1217.6045,-1898.94,-2371.3154,-2726.7856,-2823.9968 ])
xfit = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9.])
yfit = np.array([  5.781704 ,  -4.7644916,  -2.8626044,  -2.2049518,  -1.45643  ,-2.9384856,   4.4096513,  22.384567 ,  -6.0114756, -12.337624 ])
norder = 3
xvec = np.linspace(xfit.min(), xfit.max(), num=100)

msk, poly_coeff = utils.robust_polyfit(xfit, yfit, norder, sigma=3.0, function='polynomial')

msk_new, poly_coeff_new = utils.robust_polyfit_djs(xfit,yfit,norder, \
                                           function = 'polynomial', minv = None, maxv = None, bspline_par = None,\
                                           guesses = None, maxiter = 10, inmask = None, sigma = None,invvar = None,\
                                           lower = 3, upper = 3,maxdev=None,maxrej=None,groupdim=None,groupsize=None,\
                                           groupbadpix=False, grow=0,sticky=False)
robust_mask = msk == 0
robust_mask_new = msk_new == 1

plt.plot(xfit, yfit, 'ko', mfc='None',markersize=8.0, label='pca coeff')
plt.plot(xfit[~robust_mask], yfit[~robust_mask], 'ms', mfc='None',markersize=10.0, label='robust_polyfit rejected')
plt.plot(xfit[~robust_mask_new],yfit[~robust_mask_new],'r+', markersize = 20.0, label = 'robust_polyfit_djs rejected')
plt.plot(xvec, utils.func_val(poly_coeff, xvec, 'polynomial'), color='m', label='robust polyfit')
plt.plot(xvec, utils.func_val(poly_coeff_new, xvec, 'polynomial'), color='r',label = 'new robust polyfit')
plt.legend()
plt.show()


