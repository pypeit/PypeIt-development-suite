



from pypeit import msgs
import numpy as np
from IPython import embed
from pypeit import flatfield
from astropy.io import fits
import os
from pypeit.spectrographs.util import load_spectrograph
from pypeit import ginga
from pypeit import traceslits
from pypeit.core import parse
from pypeit.core import pixels
from pypeit import scienceimage
from matplotlib import pyplot as plt
from pypeit.core import flat

# Function compute_flats(flat, mstilts, slit_left, sligt_righ, thismask, inmask = None
# spec_samp_fine = 0.8, spec_samp_coarse = 50, spat_samp =  5.0, spat_illum_thresh = 0.03, trim_edg = (3.0,3.0), debug = True

#spec_samp_fine = 0.8
#spec_samp_coarse = 50.0
#spat_samp = 5.0
#trim_edg = (3.0, 3.0)
#spat_illum_thresh = 0.02
#npoly = None
#spatbkpt = None
#debug = True

# Imports for this fast running median routine
from collections import deque
from itertools import islice
from bisect import insort, bisect_left


'''def fit_flat(flat, mstilts, slit_left, slit_righ, thismask, inmask = None,spec_samp_fine = 0.8,  spec_samp_coarse = 50,
             spat_samp = 5.0, spat_illum_thresh = 0.02, npoly = None, trim_edg = (3.0,3.0), spat_bkpt = None, debug = False):


    nspec = flat.shape[0]
    nspat = flat.shape[1]
    piximg = mstilts * (nspec-1)
    # Compute the approximate number of pixels sampling each spatial pixel for this slit
    npercol = np.fmax(np.floor(np.sum(thismask)/nspec),1.0)
    # Demand at least 10 pixels per row (on average) per degree of the polynomial
    if npoly is None:
        npoly_in = 7
        npoly = np.fmax(np.fmin(npoly_in, (np.ceil(npercol/10.)).astype(int)),1)


    ximg, edgmask = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg=trim_edg)
    if inmask is None:
        inmask = np.copy(thismask)

    log_flat = np.log(np.fmax(flat, 1.0))
    inmask_log = ((flat > 1.0) & inmask)
    log_ivar = inmask_log.astype(float) # set errors to just be 1.0

    # Flat field pixels for fitting spectral direction
    fit_spec = thismask & inmask & (edgmask == False)
    isrt_spec = np.argsort(piximg[fit_spec])
    pix_fit = piximg[fit_spec][isrt_spec]
    log_flat_fit = log_flat[fit_spec][isrt_spec]
    log_ivar_fit = log_ivar[fit_spec][isrt_spec]
    inmask_log_fit = inmask_log[fit_spec][isrt_spec]
    nfit_spec = np.sum(fit_spec)
    logrej = 0.5 # rejectino threshold for spectral fit in log(image)
    msgs.info('Spectral fit of flatfield for {:}'.format(nfit_spec) + ' pixels')

    # Fit the spectral direction of the blaze
    # ToDo Figure out how to deal with the fits going crazy at the edges of the chip in spec direction
    spec_set_fine, outmask_spec, specfit, _ = utils.bspline_profile(pix_fit, log_flat_fit, log_ivar_fit,
                                                                    np.ones_like(pix_fit), inmask = inmask_log_fit,
                                                                    nord = 4, upper=logrej, lower=logrej,
                                                                    kwargs_bspline = {'bkspace':spec_samp_fine},
                                                                    kwargs_reject={'groupbadpix':True, 'maxrej': 5})


    # Debugging/checking spectral fit
    if debug:
        goodbk = spec_set_fine.mask
        specfit_bkpt, _ = spec_set_fine.value(spec_set_fine.breakpoints[goodbk])
        was_fit_and_masked = (outmask_spec == False)
        plt.clf()
        ax = plt.gca()
        ax.plot(pix_fit,log_flat_fit, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full',
                linestyle='None', label = 'all pixels')
        ax.plot(pix_fit[was_fit_and_masked],log_flat_fit[was_fit_and_masked], color='red', marker='+',
                markersize=1.5, mfc='red', fillstyle='full', linestyle='None', label='masked')
        ax.plot(pix_fit, specfit, color='cornflowerblue', label = 'fit to blaze')
        ax.plot(spec_set_fine.breakpoints[goodbk], specfit_bkpt, color='lawngreen', marker='o', markersize=2.0,
                mfc='lawngreen', fillstyle='full', linestyle='None', label='bspline breakpoints')
        ax.set_ylim((0.99*specfit.min(),1.01*specfit.max()))
        plt.legend()
        plt.xlabel('Spectral Pixel')
        plt.ylabel('log(flat counts)')
        plt.show()

    # Evaluate and save
    spec_model = np.ones_like(flat)
    spec_model[thismask], _ = np.exp(spec_set_fine.value(piximg[thismask]))
    norm_spec = np.ones_like(flat)
    norm_spec[thismask] = flat[thismask]/np.fmax(spec_model[thismask], 1.0)

    # Flat field pixels for fitting spatial direction
    slitwidth = np.median(slit_righ - slit_left) # How many pixels wide is the slit at each Y?

    fit_spat = thismask & inmask & (spec_model > 1.0)
    isrt_spat = np.argsort(ximg[fit_spat])
    ximg_fit = ximg[fit_spat][isrt_spat]
    norm_spec_fit = norm_spec[fit_spat][isrt_spat]
    norm_spec_ivar = np.ones_like(norm_spec_fit)/(spat_illum_thresh**2)
    nfit_spat = np.sum(fit_spat)

    ximg_resln = spat_samp/slitwidth
    isamp = (np.arange(nfit_spat//10)*10.0).astype(int)
    samp_width = (np.ceil(isamp.size*ximg_resln)).astype(int)

    illumquick1 = utils.fast_running_median(norm_spec_fit[isamp], samp_width)
    #illumquick1 = scipy.ndimage.filters.median_filter(norm_spec_fit[isamp], size=samp_width, mode = 'reflect')
    statinds = (ximg_fit[isamp] > 0.1) & (ximg_fit[isamp] < 0.9)
    mean = np.mean(illumquick1[statinds])
    illum_max_quick = (np.abs(illumquick1[statinds]/mean-1.0)).max()
    npad = 10000

    imed = None
    no_illum=False
    if(illum_max_quick <= spat_illum_thresh/3.0):
        ximg_in = np.concatenate((-0.2 + 0.2*np.arange(npad)/(npad - 1), ximg_fit, 1.0 + 0.2*np.arange(npad)/(npad - 1)))
        normin = np.ones(2*npad + nfit_spat)
        #msgs.info('illum_max={:7.3f}'.format(illum_max_quick))
        msgs.info('Subsampled illum fluctuations = {:7.3f}'.format(illum_max_quick) +
                  '% < spat_illum_thresh/3={:4.2f}'.format(100.0*spat_illum_thresh/3.0) +'%')
        msgs.info('Slit illumination function set to unity for this slit')
        no_illum=True
    else:
        illumquick = np.interp(ximg_fit,ximg_fit[isamp],illumquick1)
        chi_illum = (norm_spec_fit - illumquick)*np.sqrt(norm_spec_ivar)
        imed = np.abs(chi_illum) < 10.0 # 10*spat_illum_thresh ouliters, i.e. 30%
        nmed = np.sum(imed)
        med_width = (np.ceil(nmed*ximg_resln)).astype(int)
        normimg_raw = utils.fast_running_median(norm_spec_fit[imed],med_width)
        #normimg_raw = scipy.ndimage.filters.median_filter(norm_spec_fit[imed], size=med_width, mode='reflect')
        sig_res = np.fmax(med_width/15.0,0.5)
        normimg = scipy.ndimage.filters.gaussian_filter1d(normimg_raw,sig_res, mode='nearest')
        statinds = (ximg_fit[imed] > 0.1) & (ximg_fit[imed] < 0.9)
        mean = np.mean(normimg[statinds])
        normimg = normimg/mean
        # compute median value of normimg edge pixels
        if(normimg.size > 12):
            lmed = np.median(normimg[0:10])
            rmed = np.median(normimg[-1:-11:-1])
        else:
            lmed = normimg[0]
            rmed = normimg[-1]
        # Bmask regions where illumination function takes on extreme values
        if np.any(~np.isfinite(normimg)):
            msgs.error('Inifinities in slit illumination function computation normimg')
        illum_max = (np.abs(normimg[statinds]/mean - 1.0)).max()
        if (illum_max <= spat_illum_thresh):
            ximg_in = np.concatenate((-0.2 + 0.2*np.arange(npad)/(npad - 1), ximg_fit,1.0 + 0.2*np.arange(npad)/(npad-1)))
            normin = np.ones(2*npad + nfit_spat)
            #msgs.info('illum_max={:7.3f}'.format(illum_max))
            msgs.info('Illum fluctuations = {:7.3f}'.format(illum_max*100) +
                      ' < spat_illum_thresh={:4.2f}'.format(100.0*spat_illum_thresh)+'%')
            msgs.info('Slit illumination function set to unity for this slit%')
            no_illum = True
        else:
            ximg_in = np.concatenate((-0.2 + 0.2*np.arange(npad)/(npad-1), ximg_fit[imed], 1.0 + 0.2*np.arange(npad)/(npad-1)))
            normin =  np.concatenate((np.full(npad, lmed), normimg, np.full(npad,rmed)))

    if spat_bkpt is not None:
        fullbkpt = spatbkpt
    else:
        ximg_samp = np.median(ximg_fit - np.roll(ximg_fit,1))
        bsp_set = pydl.bspline(ximg_in,nord=4, bkspace=ximg_samp*100.0)
        fullbkpt = bsp_set.breakpoints
    spat_set, outmask_spat, spatfit, _ = utils.bspline_profile(ximg_in, normin, np.ones_like(normin),np.ones_like(normin),
                                                               nord=4,upper=5.0, lower=5.0,fullbkpt = fullbkpt)

    if debug:
        plt.clf()
        ax = plt.gca()
        ax.plot(ximg_fit, norm_spec_fit, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full',linestyle='None',
                label = 'all pixels')
        if imed is not None: # If we computed a median show the pixels we used
            ax.plot(ximg_fit[~imed], norm_spec_fit[~imed], color='red', marker='+',markersize=1.5, mfc='red',
            fillstyle='full', linestyle='None', label = 'masked')
            if no_illum is True:
                label =  'median spatial profile, NOT USED'
            else:
                label = 'median spatial profile'
            ax.plot(ximg_fit[imed], normimg, color='lawngreen', label = label)
        ax.plot(ximg_in, spatfit, color='cornflowerblue', label = 'final slit illumination function')
        ax.set_ylim((np.fmax(0.8 * spatfit.min(), 0.5), 1.2 * spatfit.max()))
        ax.set_xlim(-0.02, 1.02)
        plt.legend()
        plt.xlabel('Normalized Slit Position')
        plt.ylabel('Normflat Spatial Profile')
        plt.show()

    # Evaluate and save
    illumflat = np.ones_like(flat)
    illumflat[thismask], _ = spat_set.value(ximg[thismask])
    norm_spec_spat = np.ones_like(flat)
    norm_spec_spat[thismask] = flat[thismask]/np.fmax(spec_model[thismask], 1.0)/np.fmax(illumflat[thismask],0.01)
    msgs.info('Performing illumination +scattered light flat field fit')

    # Flat field pixels for fitting spectral direction
    isrt_spec = np.argsort(piximg[thismask])
    pix_twod = piximg[thismask][isrt_spec]
    ximg_twod = ximg[thismask][isrt_spec]
    norm_twod = norm_spec_spat[thismask][isrt_spec]

    fitmask = inmask[thismask][isrt_spec] & (np.abs(norm_twod - 1.0) < 0.30)
    # Here we ignore the formal photon counting errors and simply assume that a typical error per pixel.
    # This guess is somewhat aribtrary. We then set the rejection threshold with sigrej_illum
    var_value = 0.01
    norm_twod_ivar = fitmask.astype(float)/(var_value**2)
    sigrej_illum = 3.0

    poly_basis = pydl.fpoly(2.0*ximg_twod - 1.0, npoly).T

    # Perform the fill 2d fit now
    twod_set, outmask_twod, twodfit, _ = utils.bspline_profile(pix_twod, norm_twod, norm_twod_ivar,poly_basis,
                                                               inmask = fitmask, nord = 4,
                                                               upper=sigrej_illum, lower=sigrej_illum,
                                                               kwargs_bspline = {'bkspace':spec_samp_coarse},
                                                               kwargs_reject={'groupbadpix':True, 'maxrej': 10})

    if debug:
        resid = (norm_twod  - twodfit)
        badpix = (outmask_twod == False) & fitmask
        goodpix = outmask_twod & fitmask
        plt.clf()
        ax = plt.gca()
        ax.plot(pix_twod[goodpix], resid[goodpix], color='k', marker='o', markersize=0.2, mfc='k', fillstyle='full',linestyle='None',
                label = 'good points')
        ax.plot(pix_twod[badpix],resid[badpix], color='red', marker='+',markersize=0.5, mfc='red', fillstyle='full', linestyle='None', label='masked')
        plt.hlines(sigrej_illum*var_value,pix_twod.min(),pix_twod.max(), color='lawngreen',linestyle='--',
                   label='rejection thresholds',zorder=10,linewidth=2.0)
        plt.hlines(-sigrej_illum*var_value,pix_twod.min(),pix_twod.max(), color='lawngreen',linestyle='--',
                   zorder=10,linewidth=2.0)
        ax.set_ylim((-0.05,0.05))
        ax.set_xlim((pix_twod.min(), pix_twod.max()))
        plt.legend()
        plt.xlabel('Spectral Pixel')
        plt.ylabel('Residuals from pixelflat 2-d fit')
        plt.show()

        plt.clf()
        ax = plt.gca()
        ax.plot(ximg_twod[goodpix], resid[goodpix], color='k', marker='o', markersize=0.2, mfc='k', fillstyle='full',
                linestyle='None',
                label='good points')
        ax.plot(ximg_twod[badpix], resid[badpix], color='red', marker='+', markersize=0.5, mfc='red', fillstyle='full',
                linestyle='None', label='masked')
        plt.hlines(sigrej_illum * var_value, ximg_twod.min(), ximg_twod.max(), color='lawngreen', linestyle='--',
                   label='rejection thresholds', zorder=10,linewidth=2.0)
        plt.hlines(-sigrej_illum * var_value, ximg_twod.min(), ximg_twod.max(), color='lawngreen', linestyle='--',
                   zorder=10,linewidth=2.0)
        ax.set_ylim((-0.05, 0.05))
        ax.set_xlim(-0.02, 1.02)
        plt.legend()
        plt.xlabel('Normalized Slit Position')
        plt.ylabel('Residuals from pixelflat 2-d fit')
        plt.show()

    # Evaluate and save
    twod_model = np.ones_like(flat)
    twod_this = np.zeros_like(twodfit)
    twod_this[isrt_spec] = twodfit
    twod_model[thismask] = twod_this

    # Compute all the final output images output
    pixelflat = np.ones_like(flat)
    flat_model = np.ones_like(flat)
    flat_model[thismask] = twod_model[thismask]*np.fmax(illumflat[thismask],0.05)*np.fmax(spec_model[thismask],1.0)
    pixelflat[thismask] = flat[thismask]/flat_model[thismask]

    # ToDo Add some code here to treat the edges and places where fits go bad?

    return (pixelflat[thismask], illumflat[thismask], flat_model[thismask])

'''

type = 'ESI'
devpath = os.getenv('PYPEIT_DEV')

if type == 'LRIS_red':
    det = 1
    sdet = parse.get_dnum(det, prefix=False)
    rawpath = devpath + '/RAW_DATA/Keck_LRIS_red/multi_400_8500_d560/'
    masterpath = devpath + '/REDUX_OUT/Keck_LRIS_red/multi_400_8500_d560/MF_keck_lris_red/'

    # Read in the msbias for bias subtraction
    biasfile = masterpath + 'MasterBias_A_' + sdet + '_aa.fits'
    msbias = fits.getdata(biasfile)
    # Read in and process flat field images
    pixflat_image_files = np.core.defchararray.add(rawpath, ['r170320_2057.fits','r170320_2058.fits','r170320_2059.fits']).tolist()
    spectro_name = 'keck_lris_red'
    spectrograph = load_spectrograph(spectrograph=spectro_name)
    par = spectrograph.default_pypeit_par()
    flatField = flatfield.FlatField(spectrograph, file_list=pixflat_image_files,det=det, par=par['calibrations']['pixelflatframe']
                                    , msbias = msbias)
    flatimg = flatField.build_pixflat()
    # Read in the tilts
    tiltsfile = masterpath + 'MasterTilts_A_' + sdet + '_aa.fits'
    mstilts = fits.getdata(tiltsfile)
    # Read in the tslits_dict
    traceslitsroot = masterpath + 'MasterTrace_A_' + sdet + '_aa'
    Tslits = traceslits.TraceSlits.from_master_files(traceslitsroot)
    tslits_dict = {}
    tslits_dict['lcen']=Tslits.lcen
    tslits_dict['rcen']=Tslits.rcen
    tslits_dict['slitpix'] = pixels.slit_pixels(tslits_dict['lcen'],tslits_dict['rcen'], flatimg.shape, Tslits.par['pad'])
elif type == 'ESI':
    flatfile = '/Users/joe/REDUX/esi_redux/Mar_2008/Flats/FlatECH10_1x1_D.fits.gz'
    piximgfile = '/Users/joe/REDUX/esi_redux/Mar_2008/Final/f_esi1044.fits.gz'
    waveimgfile = '/Users/joe/REDUX/esi_redux/Mar_2008/Arcs/ArcECH10_1x1IMG.fits.gz'
    sedg_file = '/Users/joe/REDUX/esi_redux/Mar_2008/Flats/SEdgECH10_1x1.fits.gz'
    flatimg = fits.getdata(flatfile)
    (nspec, nspat) = flatimg.shape
    piximg = fits.getdata(piximgfile, 3)
    mstilts = piximg/nspec
    slit_edges = fits.getdata(sedg_file,0)
    tslits_dict = {}
    tslits_dict['lcen']=slit_edges[0,:,:].T
    tslits_dict['rcen']=slit_edges[1,:,:].T
    tslits_dict['slitpix'] = pixels.slit_pixels(tslits_dict['lcen'],tslits_dict['rcen'], flatimg.shape, 0)

# View the slits?
#slit_ids = [trace_slits.get_slitid(flat.shape, tslits_dict['lcen'], tslits_dict['rcen'], ii)[0]
#            for ii in range(tslits_dict['lcen'].shape[1])]
#viewer, ch = ginga.show_image(flat, cuts = (1000.0,30000.0),clear=True)
#ginga.show_slits(viewer, ch, tslits_dict['lcen'], tslits_dict['rcen'], slit_ids)  # , args.det)

maskslits = np.zeros(tslits_dict['lcen'].shape[1], dtype=bool)
gdslits = np.where(~maskslits)[0]

pixelflat = np.ones_like(flatimg)
illumflat = np.ones_like(flatimg)
flat_model = np.zeros_like(flatimg)

debug = True
# Loop on slits
for slit in gdslits:
    msgs.info("Computing flat field image for slit: {:d}".format(slit + 1))
    slit_left = tslits_dict['lcen'][:, slit]
    slit_righ = tslits_dict['rcen'][:, slit]
    thismask = (tslits_dict['slitpix'] == slit + 1)
    inmask = None # in the future set this to the bpm
    sys.exit(-1)
    pixelflat[thismask], illumflat[thismask], flat_model[thismask] = flat.fit_flat(flatimg, mstilts, thismask,
                                                                                   slit_left, slit_righ,inmask=inmask,
                                                                                   debug = debug)


ginga.show_image(pixelflat,cuts = (0.9,1.1),chname='pixeflat', wcs_match=True, clear=True)
ginga.show_image(illumflat,cuts = (0.9,1.1), chname='illumflat', wcs_match=True)
ginga.show_image(flatimg,chname='flat', wcs_match=True)
ginga.show_image(flat_model,chname='flat_model',wcs_match = True)

'''
    bkspace = 1.0/nsamp # This is the spatial sampling interval in units of fractional slit width

    fit_spat = thismask & inmask
    isrt_spat = np.argsort(ximg[fit_spat])
    ximg_fit = ximg[fit_spat][isrt_spat]
    norm_spec_fit = norm_spec[fit_spat][isrt_spat]
    norm_spec_ivar = np.ones_like(norm_spec_fit)/(0.03**2)
    sigrej_illum = 3.0
    nfit_spat = np.sum(fit_spat)
    msgs.info('Fit to flatfield slit illumination function {:}'.format(nfit_spat) + ' pixels')
    # Fit the slit illumination function now
    spat_set, outmask_spat, spatfit, _ = utils.bspline_profile(ximg_fit, norm_spec_fit, norm_spec_ivar,
                                                               np.ones_like(ximg_fit),nord = 4,
                                                               upper=sigrej_illum,
                                                               lower=sigrej_illum,
                                                               kwargs_bspline = {'bkspace':bkspace},
                                                               kwargs_reject={'groupbadpix':True, 'maxrej': 5})

    # Debugging/checking
    if debug:
        goodbk = spat_set.mask
        spatfit_bkpt, _ = spat_set.value(spat_set.breakpoints[goodbk])
        plt.clf()
        ax = plt.gca()
        was_fit_and_masked = (outmask_spat == False)
        ax.plot(ximg_fit, norm_spec_fit, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full',
                linestyle='None')
        ax.plot(ximg_fit[was_fit_and_masked],norm_spec_fit[was_fit_and_masked], color='red', marker='+',
                markersize=1.5, mfc='red', fillstyle='full', linestyle='None')
        ax.plot(ximg_fit, spatfit, color='cornflowerblue')
        ax.plot(spat_set.breakpoints[goodbk], spatfit_bkpt, color='lawngreen', marker='o', markersize=2.0, mfc='lawngreen', fillstyle='full', linestyle='None')
        ax.set_ylim((np.fmax(0.8*spatfit.min(),0.5),1.2*spatfit.max()))
        ax.set_xlim(-0.2, 1.2)
        plt.show()

'''
