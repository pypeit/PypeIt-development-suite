


import numpy as np
from pypeit.core import extract
from pypeit import msgs
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from matplotlib import pyplot as plt

from pypeit import ginga
from pypeit.core import pydl
from pypeit import utils
from pypeit.core import pixels
from sklearn.decomposition import PCA
from pypeit import specobjs
from pypeit.core import extract
from astropy.stats import SigmaClip
from pydl.pydlutils.spheregroup import spheregroup


def pca_trace(xcen, usepca = None, npca = None, pca_explained_var=99.0,coeff_npoly = None, cen_npoly = 3, debug=True):
    """
    Using sklearn PCA tracing
    xcen: flux weighted centroid, numpy array
    usepca: bool array with the size equal to the number of orders. True are orders that need to be PCAed.
                    default is None and hence all orders will be PCAed.
    npca: number of PCA components you want to keep. default is None and it will be assigned automatically by
                    calculating the number of components contains approximately 99% of the variance
    pca_explained_var: explained variance cut used for auto choosing npca.
    coeff_npoly: order of polynomial used for PCA coefficients fitting
    cen_npoly: order of polynomail for center of traces
    debug:
    :return: pca fitting result
    """

    nspec = xcen.shape[0]
    norders = xcen.shape[1]

    if usepca is None:
        usepca = np.zeros(norders,dtype=bool)

    # use_order = True orders used to predict the usepca = True bad orders
    use_order = np.invert(usepca)
    ngood = np.sum(use_order)

    if npca is None:
        pca_full = PCA()
        xcen_use = (xcen[:, use_order] - np.mean(xcen[:, use_order], 0)).T
        pca_full.fit(xcen_use)
        var = np.cumsum(np.round(pca_full.explained_variance_ratio_, decimals=3) * 100)
        if var[0]>=pca_explained_var:
            npca = 1
            msgs.info('The first PCA component contains more than {:5.3f} of the information'.format(pca_explained_var))
        else:
            npca = int(np.ceil(np.interp(pca_explained_var, var,np.arange(norders)+1)))
            msgs.info('Truncated PCA to contain {:5.3f}'.format(pca_explained_var) + '% of the total variance. ' +
                      'Number of components to keep is npca = {:d}'.format(npca))
    else:
        npca = int(npca)

    if ngood < npca:
        msgs.warn('Not enough good traces for a PCA fit: ngood = {:d}'.format(ngood) + ' is < npca = {:d}'.format(npca))
        msgs.warn('Using the input trace f or now')
        return xcen

    if coeff_npoly is None:
        coeff_npoly = int(np.fmin(np.fmax(np.floor(3.3*ngood/norders),1.0),3.0))

    # Polynomail coefficient for center of traces
    if cen_npoly is None:
        cen_npoly = int(np.fmin(np.fmax(np.floor(3.3 * ngood / norders), 1.0), 3.0))

    # Polynomial coefficient for PCA coefficients
    npoly_vec =np.zeros(npca, dtype=int)
    # Fit first pca dimension (with largest variance) with a higher order npoly depending on number of good orders.
    # Fit all higher dimensions (with lower variance) with a line
    # Cascade down and use lower order polynomial for PCA directions that contain less variance
    for ipoly in range(npca):
        npoly_vec[ipoly] = np.fmax(coeff_npoly - ipoly,1)

    pca = PCA(n_components=npca)
    xcen_use = (xcen[:,use_order] - np.mean(xcen[:,use_order],0)).T
    pca_coeffs_use = pca.fit_transform(xcen_use)
    pca_vectors = pca.components_

    order_vec = np.arange(norders,dtype=float)
    # pca_coeffs = np.zeros((norders, npca))
    pca_coeffs_new = np.zeros((norders, npca))
    # Now loop over the dimensionality of the compression and perform a polynomial fit to
    for idim in range(npca):
        # Only fit the use_order orders, then use this to predict the others
        xfit = order_vec[use_order]
        yfit = pca_coeffs_use[:,idim]
        norder = npoly_vec[idim]

        ## No longer used, replaced by the robust_polyfit_djs
        # msk, poly_coeff = utils.robust_polyfit(xfit, yfit, norder, sigma = 3.0, function='polynomial')
        # pca_coeffs[:,idim] = utils.func_val(poly_coeff, order_vec, 'polynomial')

        # TESTING traceset fitting
        #xtemp = xfit.reshape(1, xfit.size)
        #ytemp = yfit.reshape(1, yfit.size)
        #tset = pydl.xy2traceset(xtemp, ytemp, ncoeff=norder,func='polynomial')
        #tset_yfit = tset.yfit.reshape(tset.yfit.shape[1])

        ## Test new robust fitting with djs_reject
        msk_new, poly_coeff_new = utils.robust_polyfit_djs(xfit, yfit, norder, \
                                                   function='polynomial', minv=None, maxv=None, bspline_par=None, \
                                                   guesses=None, maxiter=10, inmask=None, sigma=None, invvar=None, \
                                                   lower=5, upper=5, maxdev=None, maxrej=None, groupdim=None,
                                                   groupsize=None, \
                                                   groupbadpix=False, grow=0, sticky=False)
        pca_coeffs_new[:,idim] = utils.func_val(poly_coeff_new, order_vec, 'polynomial')

        if debug:
            # Evaluate the fit
            xvec = np.linspace(order_vec.min(),order_vec.max(),num=100)
            robust_mask_new = msk_new == 1
            plt.plot(xfit, yfit, 'ko', mfc='None', markersize=8.0, label='pca coeff')
            plt.plot(xfit[~robust_mask_new], yfit[~robust_mask_new], 'r+', markersize=20.0,label='robust_polyfit_djs rejected')
            plt.plot(xvec, utils.func_val(poly_coeff_new, xvec, 'polynomial'),ls='-.', color='steelblue',
                     label='robust_polyfit_djs norder=%s'%str(norder))
            #(_,tset_fit) = tset.xy(xpos=xvec.reshape(1,xvec.size))
            #plt.plot(xfit[~robust_mask], yfit[~robust_mask], 'ms', mfc='None', markersize=10.0,label='robust_polyfit rejected')
            #plt.plot(xfit[~robust_mask_new], yfit[~robust_mask_new], 'r+', markersize=20.0,label='robust_polyfit_djs rejected')
            #plt.plot(xfit[~tset_mask],yfit[~tset_mask], 'bo', markersize = 10.0, label = 'traceset rejected')
            #plt.plot(xvec, utils.func_val(poly_coeff, xvec, 'polynomial'),ls='--', color='m', label='robust polyfit')
            #yfit_tset = tset_fit[0,:]
            #robust_mask = msk == 0
            #tset_mask = tset.outmask[0,:]
            #plt.plot(xvec, yfit_tset,ls=':', color='b',label='traceset')
            plt.xlabel('Order Vector', fontsize=14)
            plt.ylabel('PCA Fitting', fontsize=14)
            plt.title('PCA Fitting on Good Orders')
            plt.legend()
            plt.show()

    #ToDo should we be masking the bad orders here and interpolating/extrapolating?
    spat_mean = np.mean(xcen,0)

    #msk_spat, poly_coeff_spat = utils.robust_polyfit(order_vec, spat_mean, ncen, sigma = 3.0, function = 'polynomial')
    msk_spat, poly_coeff_spat = utils.robust_polyfit_djs(order_vec, spat_mean, cen_npoly, \
                                                       function='polynomial', minv=None, maxv=None, bspline_par=None, \
                                                       guesses=None, maxiter=10, inmask=None, sigma=None, invvar=None, \
                                                       lower=3, upper=3, maxdev=None, maxrej=None, groupdim=None,
                                                       groupsize=None, \
                                                       groupbadpix=False, grow=0, sticky=False)
    if debug:
        robust_mask_spat = msk_spat == 1
        plt.plot(order_vec, spat_mean, 'ko', mfc='None', markersize=8.0, label='order center')
        plt.plot(order_vec[~robust_mask_new], spat_mean[~robust_mask_new], 'r+', markersize=20.0,
                 label='robust_polyfit_djs rejected')
        plt.plot(order_vec, utils.func_val(poly_coeff_spat, order_vec, 'polynomial'), ls='-.', color='steelblue',
                 label='robust_polyfit_djs norder=%s'%str(cen_npoly))
        plt.xlabel('Order Vector',fontsize=14)
        plt.ylabel('Order Center Position',fontsize=14)
        plt.title('Order Center Fitting')
        plt.legend()
        plt.show()

    ibad = np.where(msk_spat == 1)
    spat_mean[ibad] = utils.func_val(poly_coeff_spat,order_vec[ibad],'polynomial')

    #pca_fit = np.outer(np.ones(nspec), spat_mean) + np.outer(pca.mean_,np.ones(norders)) + (np.dot(pca_coeffs, pca_vectors)).T
    pca_fit = np.outer(np.ones(nspec), spat_mean) + np.outer(pca.mean_,np.ones(norders)) + (np.dot(pca_coeffs_new, pca_vectors)).T

    return pca_fit


def ech_objfind(image, ivar, ordermask, slit_left, slit_righ,inmask=None,plate_scale=0.2,ncoeff = 5,
                npca=None,coeff_npoly=None, cen_npoly=3,
                min_snr=0.0,nabove_min_snr=0,pca_explained_var=99.0,pca_percentile=20.0,snr_pca=3.0,
                box_radius=2.0,sig_thresh=5.,show_peaks=False,show_fits=False,show_trace=False,debug=True):
    """
    Object finding for Echelle spectragraph
    image :  float ndarray
        Image to search for objects from. This image has shape (nspec, nspat) image.shape where the first dimension (nspec)
        is spectral, and second dimension (nspat) is spatial. Note this image can either have the sky background in it, or have already been sky subtracted.
        Object finding works best on sky-subtracted images, but often one runs on the frame with sky first to identify the brightest
        objects which are then masked (see skymask below) in sky subtraction.
    ivar: ivar for your science image, same shape with image
    ordermask: ordermask, 2D array
    slit_left:  float ndarray
        Left boundary of slit/order to be extracted (given as floating pt pixels). This a 1-d array with shape (nspec, 1)
        or (nspec)
    slit_righ:  float ndarray
        Left boundary of slit/order to be extracted (given as floating pt pixels). This a 1-d array with shape (nspec, 1)
        or (nspec)
    inmask: mask for each order, parsed into robust poly fit
    plate_scale: plate scale of your detector, in unit of arcsec/pix
    ncoeff: order number for trace
    npca: number of PCA components you want to keep. default is None and it will be assigned automatically by
                    calculating the number of components contains approximately 99% of the variance
    npoly_cen: order of polynomial used for PCA coefficients fitting
    npca: Number of PCA components you want to keep
    min_snr: minimum SNR for object finding
    nabove_min_snr: how many orders with SNR>min_snr you require for a TRUE object
    pca_percentile: percentile used for determining which order is a bad order
    snr_pca: SNR used for determining which order is a bad order
                    if an order with ((SNR_now < np.percentile(SNR_now, pca_percentile)) &
                    (SNR_now < snr_pca)) | sobjs_trim[indx].ech_usepca, then this order will be PCAed
    box_radius: box_car extraction radius
    sig_thresh: threshord for finding objects
    show_peaks: whether plotting the QA of peak finding of your object in each order
    show_fits: Plot trace fitting
    show_trace: whether display the resulting traces on top of the image
    debug:
    :return: all objects found
    """


    if inmask is None:
        inmask = (ordermask > 0)

    frameshape = image.shape
    nspec = frameshape[0]
    norders = slit_left.shape[1]

    if isinstance(plate_scale,(float, int)):
        plate_scale_ord = np.full(norders, plate_scale)  # 0.12 binned by 3 spatially for HIRES
    elif isinstance(plate_scale,(np.ndarray, list, tuple)):
        if len(plate_scale) == norders:
            plate_scale_ord = plate_scale
        elif len(plate_scale) == 1:
            plate_scale_ord = np.full(norders, plate_scale[0])
        else:
            msgs.error('Invalid size for plate_scale. It must either have one element or norders elements')
    else:
        msgs.error('Invalid type for plate scale')

    specmid = nspec // 2
    slit_width = slit_righ - slit_left
    spec_vec = np.arange(nspec)
    slit_spec_pos = nspec/2.0
    slit_spat_pos = np.zeros((norders, 2))
    for iord in range(norders):
        slit_spat_pos[iord, :] = (np.interp(slit_spec_pos, spec_vec, slit_left[:,iord]), np.interp(slit_spec_pos, spec_vec, slit_righ[:,iord]))

    # Loop over orders and find objects
    sobjs = specobjs.SpecObjs()
    # ToDo replace orderindx with the true order number here? Maybe not. Clean up slitid and orderindx!
    for iord  in range(norders):
        msgs.info('Finding objects on slit # {:d}'.format(iord + 1))
        thismask = ordermask == (iord + 1)
        inmask_iord = inmask & thismask
        specobj_dict = {'setup': 'HIRES', 'slitid': iord + 1, 'scidx': 0,'det': 1, 'objtype': 'science'}
        # ToDO Some of the objfind parameters probably need to be passed in here, and also be argumentus to this routine
        sobjs_slit, skymask[thismask], objmask[thismask], proc_list = \
            extract.objfind(image, thismask, slit_left[:,iord], slit_righ[:,iord], inmask=inmask_iord,show_peaks=show_peaks,
                            show_fits=show_fits, show_trace=show_trace, specobj_dict = specobj_dict, sig_thresh = sig_thresh)
        # ToDO make the specobjs _set_item_ work with expressions like this spec[:].orderindx = iord
        for spec in sobjs_slit:
            spec.ech_orderindx = iord
        sobjs.add_sobj(sobjs_slit)


    nfound = len(sobjs)

    # Compute the FOF linking length based on the instrument place scale and matching length FOFSEP = 1.0"
    FOFSEP = 1.0 # separation of FOF algorithm in arcseconds
    FOF_frac = FOFSEP/(np.median(slit_width)*np.median(plate_scale_ord))

    # Feige: made the code also works for only one object found in one order
    # Run the FOF. We use fake coordinaes
    fracpos = sobjs.spat_fracpos
    ra_fake = fracpos/1000.0 # Divide all angles by 1000 to make geometry euclidian
    dec_fake = 0.0*fracpos
    if nfound>1:
        (ingroup, multgroup, firstgroup, nextgroup) = spheregroup(ra_fake, dec_fake, FOF_frac/1000.0)
        group = ingroup.copy()
        uni_group, uni_ind = np.unique(group, return_index=True)
        nobj = len(uni_group)
        msgs.info('FOF matching found {:d}'.format(nobj) + ' unique objects')
    elif nfound==1:
        group = np.zeros(1,dtype='int')
        uni_group, uni_ind = np.unique(group, return_index=True)
        nobj = len(group)
        msgs.warn('Only find one object no FOF matching is needed')

    gfrac = np.zeros(nfound)
    for jj in range(nobj):
        this_group = group == uni_group[jj]
        gfrac[this_group] = np.median(fracpos[this_group])

    uni_frac = gfrac[uni_ind]

    sobjs_align = sobjs.copy()
    # Now fill in the missing objects and their traces
    for iobj in range(nobj):
        for iord in range(norders):
            # Is there an object on this order that grouped into the current group in question?
            on_slit = (group == uni_group[iobj]) & (sobjs_align.ech_orderindx == iord)
            if not np.any(on_slit):
                # Add this to the sobjs_align, and assign required tags
                thisobj = specobjs.SpecObj(frameshape, slit_spat_pos[iord,:], slit_spec_pos, det = sobjs_align[0].det,
                                           setup = sobjs_align[0].setup, slitid = (iord + 1),
                                           scidx = sobjs_align[0].scidx, objtype=sobjs_align[0].objtype)
                thisobj.ech_orderindx = iord
                thisobj.spat_fracpos = uni_frac[iobj]
                thisobj.trace_spat = slit_left[:,iord] + slit_width[:,iord]*uni_frac[iobj] # new trace
                thisobj.trace_spec = spec_vec
                thisobj.spat_pixpos = thisobj.trace_spat[specmid]
                thisobj.set_idx()
                # Use the real detections of this objects for the FWHM
                this_group = group == uni_group[iobj]
                # Assign to the fwhm of the nearest detected order
                imin = np.argmin(np.abs(sobjs_align[this_group].ech_orderindx - iord))
                thisobj.fwhm = sobjs_align[imin].fwhm
                thisobj.maskwidth = sobjs_align[imin].maskwidth
                thisobj.ech_fracpos = uni_frac[iobj]
                thisobj.ech_group = uni_group[iobj]
                thisobj.ech_usepca = True
                sobjs_align.add_sobj(thisobj)
                group = np.append(group, uni_group[iobj])
                gfrac = np.append(gfrac, uni_frac[iobj])
            else:
                # ToDo fix specobjs to get rid of these crappy loops!
                for spec in sobjs_align[on_slit]:
                    spec.ech_fracpos = uni_frac[iobj]
                    spec.ech_group = uni_group[iobj]
                    spec.ech_usepca = False

    # Some code to ensure that the objects are sorted in the sobjs_align by fractional position on the order and by order
    # respectively
    sobjs_sort = specobjs.SpecObjs()
    for iobj in range(nobj):
        this_group = group == uni_group[iobj]
        this_sobj = sobjs_align[this_group]
        sobjs_sort.add_sobj(this_sobj[np.argsort(this_sobj.ech_orderindx)])

    # Loop over the objects and perform a quick and dirty extraction to assess S/N.
    varimg = utils.calc_ivar(ivar)
    flux_box = np.zeros((nspec, norders, nobj))
    ivar_box = np.zeros((nspec, norders, nobj))
    mask_box = np.zeros((nspec, norders, nobj))
    SNR_arr = np.zeros((norders, nobj))
    for iobj in range(nobj):
        for iord in range(norders):
            indx = (sobjs_sort.ech_group == uni_group[iobj]) & (sobjs_sort.ech_orderindx == iord)
            spec = sobjs_sort[indx]
            thismask = ordermask == (iord + 1)
            inmask_iord = inmask & thismask
            box_rad_pix = box_radius/plate_scale_ord[iord]
            flux_tmp  = extract.extract_boxcar(image*inmask_iord, spec.trace_spat,box_rad_pix, ycen = spec.trace_spec)
            var_tmp  = extract.extract_boxcar(varimg*inmask_iord, spec.trace_spat,box_rad_pix, ycen = spec.trace_spec)
            ivar_tmp = utils.calc_ivar(var_tmp)
            pixtot  = extract.extract_boxcar(ivar*0 + 1.0, spec.trace_spat,box_rad_pix, ycen = spec.trace_spec)
            mask_tmp = (extract.extract_boxcar(ivar*inmask_iord == 0.0, spec.trace_spat,box_rad_pix, ycen = spec.trace_spec) != pixtot)
            flux_box[:,iord,iobj] = flux_tmp*mask_tmp
            ivar_box[:,iord,iobj] = np.fmax(ivar_tmp*mask_tmp,0.0)
            mask_box[:,iord,iobj] = mask_tmp
            (mean, med_sn, stddev) = sigma_clipped_stats(flux_box[mask_tmp,iord,iobj]*np.sqrt(ivar_box[mask_tmp,iord,iobj]),
                                                         sigma_lower=5.0,sigma_upper=5.0)
            SNR_arr[iord,iobj] = med_sn



    # Purge objects with low SNR and that don't show up in enough orders
    keep_obj = np.zeros(nobj,dtype=bool)
    sobjs_trim = specobjs.SpecObjs()
    uni_group_trim = np.array([],dtype=int)
    uni_frac_trim =  np.array([],dtype=float)
    for iobj in range(nobj):
        if (np.sum(SNR_arr[:,iobj] > min_snr) >= nabove_min_snr):
            keep_obj[iobj] = True
            ikeep = sobjs_sort.ech_group == uni_group[iobj]
            sobjs_trim.add_sobj(sobjs_sort[ikeep])
            uni_group_trim = np.append(uni_group_trim, uni_group[iobj])
            uni_frac_trim = np.append(uni_frac_trim, uni_frac[iobj])
        else:
            msgs.info('Purging object #{:d}'.format(iobj) + ' which does not satisfy min_snr > {:5.2f}'.format(min_snr) +
                      ' on at least nabove_min_snr >= {:d}'.format(nabove_min_snr) + ' orders')

    nobj_trim = np.sum(keep_obj)
    if nobj_trim == 0:
        return specobjs.SpecObjs()

    SNR_arr_trim = SNR_arr[:,keep_obj]

    # Do a final loop over objects and make the final decision about which orders will be interpolated/extrapolated by the PCA
    for iobj in range(nobj_trim):
        SNR_now = SNR_arr_trim[:,iobj]
        indx = (sobjs_trim.ech_group == uni_group_trim[iobj])
        # PCA interp/extrap if:
        #      (SNR is below pca_percentile of the total SNRs) AND (SNR < snr_pca)
        #                                 OR
        #      (if this order was not originally traced by the object finding, see above)
        usepca = ((SNR_now < np.percentile(SNR_now, pca_percentile)) & (SNR_now < snr_pca)) | sobjs_trim[indx].ech_usepca
        # ToDo fix specobjs to get rid of these crappy loops!
        for iord, spec in enumerate(sobjs_trim[indx]):
            spec.ech_usepca = usepca[iord]
            if usepca[iord]:
                msgs.info('Using PCA to predict trace for object #{:d}'.format(iobj) + ' on order #{:d}'.format(iord))

    sobjs_final = sobjs_trim.copy()
    # Loop over the objects one by one and adjust/predict the traces
    npoly_cen = 3
    pca_fits = np.zeros((nspec, norders, nobj_trim))
    for iobj in range(nobj_trim):
        igroup = sobjs_final.ech_group == uni_group_trim[iobj]
        # PCA predict the masked orders which were not traced
        pca_fits[:, :, iobj] = pca_trace((sobjs_final[igroup].trace_spat).T,usepca = None,
                                         npca=npca, pca_explained_var=pca_explained_var, coeff_npoly=coeff_npoly,
                                         cen_npoly=cen_npoly, debug=debug)
        #pca_fits[:,:,iobj] = pca_trace((sobjs_final[igroup].trace_spat).T,usepca = usepca,
        #                               npca = npca, pca_explained_var=pca_explained_var,coeff_npoly = coeff_npoly,
        #                               cen_npoly = cen_npoly, debug=debug)
        # usepca = sobjs_final[igroup].ech_usepca,
        # Perform iterative flux weighted centroiding using new PCA predictions
        xinit_fweight = pca_fits[:,:,iobj].copy()
        inmask_now = inmask & (ordermask > 0)
        xfit_fweight = extract.iter_tracefit(image, xinit_fweight, ncoeff, inmask = inmask_now, show_fits=show_fits)
        # Perform iterative Gaussian weighted centroiding
        xinit_gweight = xfit_fweight.copy()
        xfit_gweight = extract.iter_tracefit(image, xinit_gweight, ncoeff, inmask = inmask_now, gweight=True,show_fits=show_fits)
        # Assign the new traces
        for iord, spec in enumerate(sobjs_final[igroup]):
            spec.trace_spat = xfit_gweight[:,iord]
            spec.spat_pixpos = spec.trace_spat[specmid]


    # Set the IDs
    sobjs_final.set_idx()
    if show_trace:
        viewer, ch = ginga.show_image(objminsky*(ordermask > 0))
        for iobj in range(nobj_trim):
            for iord in range(norders):
                ginga.show_trace(viewer, ch, pca_fits[:,iord, iobj], str(uni_frac[iobj]), color='yellow')

        for spec in sobjs_trim:
            color = 'green' if spec.ech_usepca else 'magenta'
            ginga.show_trace(viewer, ch, spec.trace_spat, spec.idx, color=color)

        #for spec in sobjs_final:
        #    color = 'red' if spec.ech_usepca else 'green'
        #    ginga.show_trace(viewer, ch, spec.trace_spat, spec.idx, color=color)

    return sobjs_final

# HIRES
spectro = 'MAGE'
if spectro == 'HIRES':
    hdu = fits.open('/Users/feige//Dropbox/hires_fndobj/f_hires0184G.fits.gz')
    objminsky =hdu[2].data
    ivar  = hdu[1].data
    mask = (ivar > 0.0)
    order_str = Table.read('/Users/feige//Dropbox/hires_fndobj/OStr_G_04.fits')
    slit_left = (order_str['LHEDG']).T
    slit_righ = (order_str['RHEDG']).T
    plate_scale = 0.36
elif spectro == 'MAGE':
    from scipy.io import readsav
    hdu = fits.open('/Users/feige/Dropbox/hires_fndobj/MAGE/f_mage2013.fits.gz')
    objminsky = hdu[0].data - hdu[2].data
    ivar  = hdu[4].data
    #ivar = utils.calc_ivar(var)
    mask = (ivar > 0.0)
    slit_left = readsav('/Users/feige/Dropbox/hires_fndobj/MAGE/left_edge.sav',python_dict=False)['left_edge'].T
    slit_righ = readsav('/Users/feige/Dropbox/hires_fndobj/MAGE/right_edge.sav',python_dict=False)['right_edge'].T
    plate_scale = 0.3
elif spectro == 'ESI':
    # ESI
    hdu = fits.open('/Users/feige//Dropbox/Cowie_2002-02-17/Final/fringe_ES.20020217.35453.fits.gz')
    sciimg = hdu[0].data
    var  = hdu[1].data
    ivar = utils.calc_ivar(var)
    mask = (var > 0.0)
    skyimg = hdu[2].data
    objminsky = sciimg - skyimg
    hdu_sedg = fits.open('/Users/feige//Dropbox/Cowie_2002-02-17/Flats/SEdgECH75_1x1.fits')
    data = hdu_sedg[0].data
    slit_left = data[0,:,:].T
    slit_righ = data[1,:,:].T
    plate_scale = 0.149
elif spectro == 'NIRES':
    from linetools import utils as ltu
    jdict = ltu.loadjson('/Users/feige//Dropbox/hires_fndobj/tilt_nires.json')
    slit_left = np.array(jdict['lcen'])
    slit_righ = np.array(jdict['rcen'])
    hdu = fits.open('/Users/feige/Dropbox/hires_fndobj/spec2d_J1724+1901_NIRES_2018Jun04T130207.856.fits')
    objminsky = hdu[1].data - hdu[3].data
    ivar = hdu[2].data
    mask = (ivar>0)
    plate_scale = 0.123
elif spectro == 'GNIRS':
    from scipy.io import readsav
    #hdu = fits.open('/Users/feige//Dropbox/hires_fndobj/sci-N20170331S0216-219.fits')
    hdu = fits.open('/Users/feige//Dropbox/hires_fndobj/GNIRS/J021514.76+004223.8/Science/J021514.76+004223.8_1/sci-N20170927S0294-297.fits')
    hdu = fits.open('/Users/feige/Dropbox/hires_fndobj/GNIRS/J005424.45+004750.2/Science/J005424.45+004750.2_7/sci-N20171021S0264-267.fits')
    hdu = fits.open('/Users/feige/Dropbox/hires_fndobj/GNIRS/J002407.02-001237.2/Science/J002407.02-001237.2_5/sci-N20171006S0236-239.fits')
    obj = hdu[0].data
    #objminsky = obj - hdu[1].data
    objminsky = hdu[1].data - obj # test negative trace
    ivar  = hdu[2].data
    #ivar = utils.calc_ivar(var)
    mask = (ivar > 0.0)
    #slit_left = readsav('/Users/feige//Dropbox/hires_fndobj/left_edge.sav', python_dict=False)['left_edge'].T
    #slit_righ = readsav('/Users/feige//Dropbox/hires_fndobj/right_edge.sav', python_dict=False)['right_edge'].T
    #slit_left = readsav('/Users/feige//Dropbox/hires_fndobj/GNIRS/J021514.76+004223.8/left_edge_J0215.sav',python_dict=False)['left_edge'].T
    #slit_righ = readsav('/Users/feige//Dropbox/hires_fndobj/GNIRS/J021514.76+004223.8/right_edge_J0215.sav',python_dict=False)['right_edge'].T
    #slit_left = readsav('/Users/feige//Dropbox/hires_fndobj/GNIRS/J005424.45+004750.2/left_edge_J0054.sav',python_dict=False)['left_edge'].T
    #slit_righ = readsav('/Users/feige//Dropbox/hires_fndobj/GNIRS/J005424.45+004750.2/right_edge_J0054.sav',python_dict=False)['right_edge'].T
    slit_left = readsav('/Users/feige/Dropbox/hires_fndobj/GNIRS/J002407.02-001237.2/left_edge_J0024.sav',python_dict=False)['left_edge'].T
    slit_righ = readsav('/Users/feige/Dropbox/hires_fndobj/GNIRS/J002407.02-001237.2/right_edge_J0024.sav',python_dict=False)['right_edge'].T
    plate_scale = 0.15


ordermask = pixels.slit_pixels(slit_left, slit_righ, objminsky.shape, 0)

#viewer, ch = ginga.show_image(objminsky)
#ginga.show_slits(viewer,ch, slit_left, slit_righ)

#yvec = np.outer(np.ones(norders), spec_vec)

#tset_left = pydl.xy2traceset(yvec, slit_left.T, ncoeff=ncoeff)
#slit_left_fit = tset_left.yfit.T

#tset_righ = pydl.xy2traceset(yvec, slit_righ.T, ncoeff=ncoeff)
#slit_righ_fit = tset_righ.yfit.T

#ginga.show_slits(viewer,ch, slit_left_fit, slit_righ_fit)

slit_mid = (slit_left + slit_righ)/2.0
#usepca = np.zeros(slit_mid.shape[1],dtype=bool)
#usepca[0:8] = True
#pca_out = pca_trace(slit_mid, npca = 4, npoly_cen = 3)
#sys.exit(-1)
#viewer, ch = ginga.show_image(ordermask > 0)
#ginga.show_slits(viewer,ch, slit_mid, pca_out)

# create the ouptut images skymask and objmask
skymask = np.zeros_like(objminsky, dtype=bool)
objmask = np.zeros_like(objminsky, dtype=bool)
image = objminsky.copy()
inmask = mask.copy()
#Routine starts here.
#------
ncoeff = 5
box_radius = 2.0  # arcseconds
min_snr = 0.2
nabove_min_snr = 2
pca_percentile = 20.0
snr_pca = 3.0
npca = 3
pca_explained_var = 99.8
sig_thresh = 50
box_radius = 2.0
ncoeff = 5
cen_npoly = 3


sobjs_final = ech_objfind(image, ivar, ordermask, slit_left, slit_righ,inmask=inmask,plate_scale=plate_scale,ncoeff = ncoeff,
                npca=npca,pca_explained_var=pca_explained_var, coeff_npoly=None, cen_npoly=cen_npoly,
                min_snr=min_snr,nabove_min_snr=nabove_min_snr,pca_percentile=pca_percentile,
                snr_pca=snr_pca,box_radius=box_radius,sig_thresh=sig_thresh,
                show_peaks=False,show_fits=False,show_trace=True,debug=True)