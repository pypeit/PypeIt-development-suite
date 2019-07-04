
from matplotlib import pyplot as plt
import numpy as np
from numpy.lib.stride_tricks import as_strided
from scipy import signal
from IPython import embed
from pypeit import utils

seed=5
rand = np.random.RandomState(seed=seed)

x = rand.randn(1000)
y = np.roll(x,10)

maxlag = 100
nsize = x.size
lags = np.arange(-nsize + 1, nsize)
corr_scipy = signal.correlate(x, y, mode='full')

lags_maxlag, corr = utils.cross_correlate(x, y, maxlag)
#corr2 = crosscorrelation2(x, y, maxlag)

#lags_maxlag = np.arange(-maxlag, maxlag + 1,dtype=float)


plt.plot(lags,corr_scipy)
plt.plot(lags_maxlag,corr)
#plt.plot(lags_maxlag,corr2)
plt.show()
