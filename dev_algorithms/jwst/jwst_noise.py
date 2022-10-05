
import os
import numpy as np
import scipy
from matplotlib import pyplot as plt
from astropy.stats import sigma_clipped_stats

from jwst import datamodels
DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']



# PypeIt imports
from jwst_utils import compute_diff, get_cuts
from pypeit.display import display
from pypeit.utils import inverse
from pypeit.core import findobj_skymask
from pypeit.core import skysub, coadd

# G395M data
rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/Raw'
output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/output'

det = 'nrs1'
# NIRSPEC 3-point dither
scifile = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + det + '_rate.fits')
bkgfile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + det + '_rate.fits')
bkgfile2 = os.path.join(rawpath_level2, 'jw02736007001_03103_00003_' + det + '_rate.fits')

sci = datamodels.open(scifile)
bkg1 = datamodels.open(bkgfile1)

sci_rate = sci.data
sci_err = sci.err
sci_var_poisson = sci.var_poisson
sci_var_rnoise = sci.var_rnoise


# Good pixel masks
gpm_sci = np.logical_not(sci.dq & DO_NOT_USE)
gpm_bkg1 = np.logical_not(bkg1.dq & DO_NOT_USE)

bkg1_rate = bkg1.data
bkg1_var_poisson = bkg1.var_poisson
bkg1_var_rnoise = bkg1.var_rnoise

gpm_diff = gpm_sci & gpm_bkg1
diff = sci_rate - bkg1_rate
#var_diff = sci.err**2 + bkg1.err**2
var_diff = sci_var_poisson + sci_var_rnoise + bkg1_var_poisson + bkg1_var_rnoise
sig_diff =np.sqrt(var_diff)
ivar_diff = inverse(var_diff)
chi = diff*np.sqrt(ivar_diff)

xlo = 700
xhi = 2048
ylo = 950
yhi = 971

#nspat, nspec = sci_rate.shape
cuts = get_cuts(sci_rate[ylo:yhi, xlo:xhi])
cuts_var = get_cuts(sci_err[ylo:yhi, xlo:xhi]**2)

# yvals = ylo + np.arange(yhi - ylo)
#slit_left = np.full(nspec, ylo)
#slit_righ = np.full(nspec, yhi)
#spec_val = xlo + np.arange(xhi - xlo)
#viewer_sci, ch_sci = display.show_image(sci_rate.T, cuts=get_cuts(sci_rate), chname='raw', clear=True)
#display.show_slits(viewer_sci, ch_sci, slit_left, slit_righ, spec_vals=spec_val, pstep=1)

display.show_image(sci_rate[ylo:yhi, xlo:xhi], cuts=cuts, chname='science', wcs_match=True)
display.show_image(sci_err[ylo:yhi, xlo:xhi]**2, cuts=cuts_var, chname='sci_var', wcs_match=True)
display.show_image(sci_var_poisson[ylo:yhi, xlo:xhi], cuts=cuts_var, chname='sci_var_poi', wcs_match=True)
display.show_image(sci_var_rnoise[ylo:yhi, xlo:xhi], cuts=cuts_var, chname='sci_var_rn', wcs_match=True)
display.show_image(bkg1_rate[ylo:yhi, xlo:xhi], cuts=cuts, chname='bkg1', wcs_match=True)
display.show_image(diff[ylo:yhi, xlo:xhi], cuts=get_cuts(diff), chname='diff', wcs_match=True)
display.show_image(sig_diff[ylo:yhi, xlo:xhi], cuts=get_cuts(sig_diff), chname='sigma', wcs_match=True)
display.show_image(chi[ylo:yhi, xlo:xhi], cuts=(-5.0,5.0), chname='chi', wcs_match=True)

# Grab only non-flagged pixels
chi_sub = chi[ylo:yhi, xlo:xhi]
gpm_diff_sub = gpm_diff[ylo:yhi, xlo:xhi]
chi_flat = chi_sub[gpm_diff_sub].flatten()
# Compute the mean and standard deviation of the chi distribution
mean, med, sigma = sigma_clipped_stats(chi_flat, sigma_lower=5.0,sigma_upper=5.0)
sig_range = 5.0*sigma

n_bins = 50
binsize = 2.0 * sig_range / n_bins
bins_histo = -sig_range + np.arange(n_bins) * binsize + binsize / 2.0
xvals = np.arange(-10.0, 10, 0.02)
gauss = scipy.stats.norm(loc=0.0, scale=1.0)

plt.figure(figsize=(12, 8))

plt.hist(chi_flat, bins_histo, density=True, histtype='step', align='mid', color='k', linewidth=3,
         label='Chi distribution, sigma={:5.3f}'.format(sigma))
plt.plot(xvals, gauss.pdf(xvals), 'c-', lw=3, label='sigma=1')
plt.ylabel('Residual distribution')
plt.xlabel('chi')
plt.xlim([-sig_range, sig_range])
plt.legend(fontsize=13, loc=2)
plt.show()

rate_test = sci_rate[ylo:yhi, xlo:xhi][gpm_sci[ylo:yhi, xlo:xhi]]
var_test = (sci_err[ylo:yhi, xlo:xhi][gpm_sci[ylo:yhi, xlo:xhi]]**2)
var_poi_test = sci_var_poisson[ylo:yhi, xlo:xhi][gpm_sci[ylo:yhi, xlo:xhi]]

plt.plot(rate_test, var_poi_test, color='black', marker='o', markersize=3, linestyle='None', label='Poisson')
plt.plot(rate_test, 1/2946.0*rate_test, color='blue',linestyle='solid', label='1/4000*rate')
plt.xlim([-0.02,0.15])
plt.ylim([-1e-5,1.5e-4])
plt.show()