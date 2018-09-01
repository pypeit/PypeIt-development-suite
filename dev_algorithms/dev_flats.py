



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
from pypeit import utils
from matplotlib import pyplot as plt
imp

type = 'LRIS_red'
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
    flat = flatField.build_pixflat()
    # Read in the tilts
    tiltsfile = masterpath + 'MasterTilts_A_' + sdet + '_aa.fits'
    mstilts = fits.getdata(tiltsfile)
    # Read in the tslits_dict
    traceslitsroot = masterpath + 'MasterTrace_A_' + sdet + '_aa'
    Tslits = traceslits.TraceSlits.from_master_files(traceslitsroot)
    tslits_dict = {}
    tslits_dict['lcen']=Tslits.lcen
    tslits_dict['rcen']=Tslits.rcen
    tslits_dict['slitpix'] = pixels.slit_pixels(tslits_dict['lcen'],tslits_dict['rcen'], flat.shape, Tslits.par['pad'])
elif type == 'ESI':
    flatfile = '/Users/joe/REDUX/esi_redux/Mar_2008/Flats/FlatECH10_1x1_D.fits.gz'
    piximgfile = '/Users/joe/REDUX/esi_redux/Mar_2008/Final/f_esi1044.fits.gz'
    waveimgfile = '/Users/joe/REDUX/esi_redux/Mar_2008/Arcs/ArcECH10_1x1IMG.fits.gz'
    sedg_file = '/Users/joe/REDUX/esi_redux/Mar_2008/Flats/SEdgECH10_1x1.fits.gz'
    flat = fits.getdata(flatfile)
    (nspec, nspat) = flat.shape
    piximg = fits.getdata(piximgfile, 3)
    mstilts = piximg/nspec
    slit_edges = fits.getdata(sedg_file,0)
    tslits_dict = {}
    tslits_dict['lcen']=slit_edges[0,:,:].T
    tslits_dict['rcen']=slit_edges[0,:,:].T
    tslits_dict['slitpix'] = pixels.slit_pixels(tslits_dict['lcen'],tslits_dict['rcen'], flat.shape, 0)


maskslits = np.zeros(tslits_dict['lcen'].shape[1], dtype=bool)
gdslits = np.where(~maskslits)[0]

# Loop on slits
for slit in gdslits:
    msgs.info("Computing flat field image for slit: {:d}".format(slit + 1))
    slit_left = tslits_dict['lcen'][:, slit]
    slit_righ = tslits_dict['rcen'][:, slit]
    thismask = (tslits_dict['slitpix'] == slit + 1)
    inmask = None # in the future set this to the bpm

    # Function compute_flats(flat, mstilts, slit_left, sligt_righ, thismask, inmask = None
    # spec_samp_fine = 0.8, spec_samp_coarse = 50, spat_samp =  5.0, spat_illum_thresh = 0.03, trim_edg = (3.0,3.0), debug = True
    spec_samp_fine = 0.8
    trim_edg = (3.0,3.0)
    debug = True


    nspec = flat.shape[0]
    nspat = flat.shape[1]
    piximg = mstilts * (nspec-1)


    ximg, edgmask = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg=trim_edg)
    if inmask == None:
        inmask = np.copy(thismask)

    log_flat = np.log(np.fmax(flat, 1.0))
    inmask_log = ((flat > 1.0) & inmask)
    log_ivar = inmask_log.astype(float) # set errors to just be 1.0

    # Flat field pixels for fitting spectral direction
    fit_spec = thismask & inmask & (edgmask == False)
    isrt_spec = np.argsort(piximg[fit_spec])
    pix_fit_spec = piximg[fit_spec][isrt_spec]
    log_flat_fit = log_flat[fit_spec][isrt_spec]
    log_ivar_fit = log_ivar[fit_spec][isrt_spec]
    nfit_spec = np.sum(fit_spec)
    logrej = 0.5 # rejectino threshold for spectral fit in log(image)
    msgs.info('Spectral fit of flatfield for {:}'.format(nfit_spec) + ' pixels')

    # Fit the Full fit now
    spec_set_fine, outmask_spec, specfit, _ = utils.bspline_profile(pix_fit_spec, log_flat_fit, log_ivar_fit, np.ones_like(pix_fit),
                                                     nord = 4, upper=logrej, lower=logrej,
                                                     kwargs_bspline = {'bkspace':spec_samp_fine},
                                                     kwargs_reject={'groupbadpix':True, 'maxrej': 5})


    # Debugging/checking
    if debug:
        goodbk = spec_set_fine.mask
        yfit_bkpt = spec_set_fine.value(spec_set_fine.breakpoints[goodbk])
        plt.clf()
        ax = plt.gca()
        was_fit_and_masked = (outmask_spec == False)
        ax.plot(pix_fit_spec,log_flat_fit, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full', linestyle='None')
        ax.plot(pix_fit_spec[was_fit_and_masked],log_flat_fit[was_fit_and_masked], color='red', marker='+', markersize=1.5, mfc='red', fillstyle='full', linestyle='None')
        ax.plot(pix_fit_spec, specfit, color='cornflowerblue')
        ax.plot(spec_set_fine.breakpoints[goodbk], yfit_bkpt, color='lawngreen', marker='o', markersize=4.0, mfc='lawngreen', fillstyle='full', linestyle='None')
        ax.set_ylim((0.99*specfit.min(),1.01*specfit.max()))
        plt.show()


    # Evaluate and save
    spec_model = np.ones_like(flat)
    spec_model[thismask], _ = np.exp(spec_set_fine.value(piximg[thismask]))
    norm_spec = flat/(np.fmax(spec_model, 1.0))

    # Flat field pixels for fitting spatial direction
    fit_spat = thismask & inmask
    isrt = np.argsort(piximg[fit_spat])
    pix_fit_spat = piximg[fit_spat][isrt]
    norm_spec_fit = norm_spec[fit_spat][isrt]
    norm_spec_ivar = 1.0/(0.03**2)
    sigrej_illum = 3.0
    nfit_spat = np.sum(fit_spat)

    msgs.info('Spectral fit of flatfield for {:}'.format(nfit) + ' pixels')
    logrej = 0.5

    # Fit the Full fit now
    spec_set_fine, outmask_spec, specfit, _ = utils.bspline_profile(pix_fit, log_flat_fit, log_ivar_fit, np.ones_like(pix_fit),
                                                     nord = 4, upper=logrej, lower=logrej,
                                                     kwargs_bspline = {'bkspace':spec_samp_fine},
                                                     kwargs_reject={'groupbadpix':True, 'maxrej': 5})



    # Now we divide out this model and perform a fit to the slit illumination function
    msgs.info('Spectral fit of flatfield for {:}'.format(nfit) + ' pixels')


    sys.exit(-1)