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

from pydl.pydlutils.bspline import bspline

sys.path.append(os.path.abspath("./"))
import dev_extract


def global_skysub(sciimg, sciivar, piximg, slitmask, edgmask,
                  skymask=None, bsp=0.6, islit=None, sigrej=3.):

    # Python indexing
    ny = sciimg.shape[0]
    #
    nslit = np.max(slitmask)
    sky_image = np.zeros_like(sciimg)

    if skymask is None:
        skymask = np.ones_like(slitmask, dtype=int)

    # Mask edges and more
    sky_slitmask = slitmask * (skymask * (sciivar > 0) & (edgmask == 0))

    if islit is None:
        nreduce = nslit
        slit_vec = np.arange(nslit) + 1
    else:
        nreduce = 1
        slit_vec = [islit]

    for jj in range(nreduce):
        slitid = slit_vec[jj]
        # Select only the pixels on this slit
        all = slitmask == slitid
        isky = sky_slitmask == slitid
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
            debugger.set_trace()

        # Full
        full_bspline = bspline(wsky, nord=4, bkspace=bsp)
        #fullbkpt = bspline_bkpts(wsky, nord=4, bkspace=bsp, / silent)
        skyset, full_out, yfit, _ = dev_extract.bspline_longslit(
            wsky, sky, sky_ivar, np.ones_like(isky),
            fullbkpt=full_bspline.breakbpoints,
            upper=sigrej, lower=sigrej, kwargs_reject={'groupbadpix':True, 'maxrej': 10})
        sky_image[all] = full_bspline.value(piximg[all])#, skyset)
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
    global_skysub(scifrcp, ivar, tilts, ordpix, edgemask)
    debugger.set_trace()


