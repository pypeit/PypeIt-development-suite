
import numpy as np
from astropy import stats
from numpy.random import randn

# Standard sigma clipping of a masked array
randvar = np.ma.MaskedArray(randn(1000,50),randn(1000,50) < -5.0)
sigclip = stats.SigmaClip(sigma=2.0, maxiters=5)
filtered_data, lower, upper = sigclip(randvar, masked=True,return_bounds=True,axis=0)
print(lower)
# The same operation but using astropy.stats.mad_std as the stdfunc
sigclip_mad = stats.SigmaClip(sigma=2.0, maxiters=5,stdfunc=stats.mad_std)
filtered_data_mad, lower_mad, upper_mad = sigclip_mad(randvar, masked=True,return_bounds=True,axis=0)
print(lower_mad)

# Note that stats.mad_stad is supposed to be usable with masked arrays
rand_test = np.ma.MaskedArray(randn(10), np.arange(10) > 5)
mad_test = stats.mad_std(rand_test)
print(mad_test)

# The problem is that according to the documentation for astropy.stats.SigmaClip:
"""
stdfunc : {'std'} or callable, optional
    The statistic or callable function/object used to compute the
    standard deviation about the center value.  If set to ``'std'``
    then having the optional `bottleneck`_ package installed will
    result in the best performance.  If using a callable
    function/object and the ``axis`` keyword is used, then it must
    be callable that can ignore NaNs (e.g. `numpy.nanstd`) and has
    an ``axis`` keyword to return an array with axis dimension(s)
    removed.  The default is ``'std'``.
"""

# The same operation but using np.ma.std as the stdfunc
sigclip_np = stats.SigmaClip(sigma=2.0, maxiters=5,stdfunc=np.ma.std)
filtered_data_np, lower_np, upper_np = sigclip_np(randvar, masked=True,return_bounds=True,axis=0)
print(lower_np)

# The same operation now using np.nanstd as the docs suggested
sigclip_nan = stats.SigmaClip(sigma=2.0, maxiters=5,stdfunc=np.nanstd)
filtered_data_nan, lower_nan, upper_nan = sigclip_nan(randvar, masked=True,return_bounds=True,axis=0)
print(lower_nan)

def nan_mad_std(data, axis=None, func=None):
    return stats.mad_std(data, axis=axis, func=func, ignore_nan=True)


# The same operation now using np.nanstd as the docs suggested
sigclip_fix = stats.SigmaClip(sigma=2.0, maxiters=5,stdfunc=nan_mad_std)
filtered_data_fix, lower_fix, upper_fix = sigclip_fix(randvar, masked=True,return_bounds=True,axis=0)
print(lower_fix)





