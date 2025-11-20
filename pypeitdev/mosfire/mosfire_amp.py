

from pypeit import ginga
from matplotlib import pyplot as plt
import glob
import numpy as np
from pypeit import spec2dobj

scipath = '/Users/joe/Dropbox/PypeIt_Redux/MOSFIRE/May20/Yband_redux/Science/'
spec2dfiles = glob.glob(scipath + 'spec2d*.fits')
nfiles = len(spec2dfiles)
for indx, ifile in enumerate(spec2dfiles):
    spec2DObj = spec2dobj.Spec2DObj.from_file(ifile, 1)
    if indx == 0:
        resids = np.zeros((nfiles, spec2DObj.sciimg.shape[0]))
    resids[indx, :] = np.median(spec2DObj.sciimg - spec2DObj.skymodel - spec2DObj.objmodel, axis=0)
    plt.plot(resids[indx,:], alpha=0.5, drawstyle='steps-mid')

plt.show()