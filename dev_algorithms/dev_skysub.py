"""
Port the global sky sub routine from LowRedux
"""

import numpy as np

from pypit import msgs
debug = debugger.init()
debug['develop'] = True
msgs.reset(debug=debug, verbosity=2)

def global_skysub(sciimg, sciivar, piximg, slitmask, edgmask,
                  skymask=None, subsample=None, npoly=None, nbkpts=None,
                   bsp=0.6, islit=None, sigrej=3.):

    # Python indexing
    nx = sciimg.shape[1]
    ny = sciimg.shape[0]
    #
    nslit = max(slitmask)
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
        isky = np.where(sky_slitmask == slitid)[0]
        if (len(isky) < 10):
            msgs.warn('Not enough sky pixels found in slit ', slitid, nsky, nall)
            continue

        # Setup (sort)
        isky = isky[np.argsort(piximg[isky])]
        wsky = piximg[isky]
        sky = sciimg[isky]
        sky_ivar = sciivar[isky]

        pos_sky = np.where((sky > 1.0) & (sky_ivar > 0.))[0]
        if len(pos_sky) > ny:
            lsky = np.log(sky[pos_sky])
            lsky_ivar = lsky * 0. + 0.1

            skybkpt = bspline_bkpts(wsky[pos_sky], nord=4, bkspace=bsp $
            , / silent)
            lskyset = bspline_longslit(wsky[pos_sky], lsky, lsky_ivar, $
            pos_sky * 0. + 1, fullbkpt = skybkpt $
            , upper = sigrej, lower = sigrej $
            , / silent, yfit = lsky_fit $
            , / groupbadpix)
            res = (sky[pos_sky] - exp(lsky_fit)) * sqrt(sky_ivar[pos_sky])
            lmask = (res LT 5.0 AND res GT -4.0)
            sky_ivar[pos_sky] = sky_ivar[pos_sky] * lmask

    fullbkpt = bspline_bkpts(wsky, nord=4, bkspace=bsp, / silent)
    skyset = bspline_longslit(wsky, sky, sky_ivar, isky * 0. + 1. $
    , / groupbadpix, maxrej = 10 $
    , fullbkpt = fullbkpt, upper = sigrej $
    , lower = sigrej, / silent, yfit = yfit)
    ;;;;;;;;;;;;;;;;;;;
    ;; JXP - - Have
    had
    to
    kludge
    this
    when
    using
    a
    kludged
    Arc
    frame
    ;      skyset = bspline_iterfit(wsky, sky $;nvvar = sky_ivar  $
    ;, / groupbadpix, maxrej = 10 $
    ;, everyn = 255L, upper = sigrej $
    ;, lower = sigrej, / silent, yfit = yfit)
    ;;;;;;;;;;;;;;;;;;;
    sky_image[all] = bspline_valu(piximg[all], skyset)
    IF
    KEYWORD_SET(CHK)
    THEN $
    x_splot, wsky, sky, psym1 = 3, xtwo = wsky, ytwo = yfit, / block
    endfor
    ;   stop
    return, sky_image
    end
