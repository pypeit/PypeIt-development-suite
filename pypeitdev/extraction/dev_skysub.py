"""
Port the global sky sub routine from LowRedux
"""

import numpy as np
import sys
import os

from astropy.io import fits

from pypit import msgs
from pypit import ardebug as debugger

debug = debugger.init()
debug['develop'] = True
msgs.reset(debug=debug, verbosity=2)

from pypit import ginga

from pydl.pydlutils.bspline import bspline

sys.path.append(os.path.abspath("./"))
import dev_extract

maskval=-999999.9,

def global_skysub(sciimg, sciivar, piximg, slitmask, edgmask,
                  skymask=None, bsp=0.6, islit=None, sigrej=3.,
                  debug=False):

    # Python indexing
    ny = sciimg.shape[0]
    #
    nslit = np.max(slitmask)
    sky_image = np.zeros_like(sciimg)

    if skymask is None:
        skymask = np.ones_like(slitmask, dtype=int)

    # Mask edges and more
    sky_slitmask = slitmask * (skymask * (sciivar > 0) & (edgmask == 0) & (sciimg != maskval))

    if islit is None:
        nreduce = nslit
        slit_vec = np.arange(nslit) + 1
    else:
        nreduce = 1
        slit_vec = [islit]

    for jj in range(nreduce):
        slitid = slit_vec[jj]
        # Select only the pixels on this slit
        all = (slitmask == slitid) & (sciimg != maskval) & (sciivar > 0.)
        isky = sky_slitmask == slitid
        debugger.set_trace() # The 2 lines above seem a bit wrong
        if (np.sum(isky) < 10):
            msgs.warn('Not enough sky pixels found in slit ', slitid,
                      np.sum(isky), np.sum(all))
            continue

        # Setup (sort)
        isrt = np.argsort(piximg[isky])
        wsky = piximg[isky][isrt]
        sky = sciimg[isky][isrt]
        sky_ivar = sciivar[isky][isrt]

        pos_sky = np.where((sky > 1.0) & (sky_ivar > 0.))[0]
        # Pre-fit!
        if len(pos_sky) > ny:
            lsky = np.log(sky[pos_sky])
            lsky_ivar = lsky * 0. + 0.1

            # Init bspline to get the sky breakpoints (kludgy)
            tmp = bspline(wsky[pos_sky], nord=4, bkspace=bsp)

            #skybkpt = bspline_bkpts(wsky[pos_sky], nord=4, bkspace=bsp $
            #, / silent)
            lskyset, outmask, lsky_fit, red_chi = dev_extract.bspline_longslit(
                wsky[pos_sky], lsky, lsky_ivar, np.ones_like(pos_sky),
                fullbkpt = tmp.breakpoints, upper=sigrej, lower=sigrej,
                kwargs_reject={'groupbadpix':True})
            res = (sky[pos_sky] - np.exp(lsky_fit)) * np.sqrt(sky_ivar[pos_sky])
            lmask = (res < 5.0) & (res > -4.0)
            sky_ivar[pos_sky] = sky_ivar[pos_sky] * lmask

        # Full
        full_bspline = bspline(wsky, nord=4, bkspace=bsp)
        skyset, full_out, yfit, _ = dev_extract.bspline_longslit(
            wsky, sky, sky_ivar, np.ones_like(sky),
            fullbkpt=full_bspline.breakpoints,
            upper=sigrej, lower=sigrej, kwargs_reject={'groupbadpix':True, 'maxrej': 10})
        sky_image[all] = skyset.value(piximg[all])[0] #, skyset)
        #debugger.set_trace()
        if debug:
            from matplotlib import pyplot as plt
            plt.clf()
            ax = plt.gca()
            ax.scatter(wsky[full_out], sky[full_out])
            ax.scatter(wsky[~full_out], sky[~full_out], color='red')
            ax.plot(wsky, yfit, color='green')
            plt.show()
    # Return
    return sky_image

# Command line execution
if __name__ == '__main__':
    # Load the test image
    hdul = fits.open('data/LRIS/Sky/lrisr_sky_test.fits')
    scifrcp = hdul[0].data
    ivar = hdul[1].data
    ordpix = hdul[2].data
    tilts = hdul[3].data*(scifrcp.shape[0])
    #
    edgemask = np.zeros_like(scifrcp)
    # Run me
    sky_image = global_skysub(scifrcp, ivar, tilts, ordpix, edgemask)
    # Show
    ginga.show_image(scifrcp-sky_image)
    debugger.set_trace()


