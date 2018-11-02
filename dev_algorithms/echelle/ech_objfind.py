


import numpy as np
from pypeit.core import extract
from pypeit import msgs
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

from pypeit import ginga
from pypeit.core import pydl
from pypeit import utils
from pypeit.core import pixels
from sklearn.decomposition import PCA
from pypeit import specobjs
from pypeit.core import extract
from astropy.stats import SigmaClip
from pydl.pydlutils.spheregroup import spheregroup


from matplotlib import pyplot as plt


#  xcen = flux weighted centroid, xcen_fit = fit to flux weighted centroid, inmask = mask of same dimension as xcen,
# usepca = None, npca =2, nspat_poly = 3
#coeff = tset_left.coeff # (norder, ncoeff)
#xcen_fit = slit_left_fit # (nspec, norder)
#ncoeff = coeff.shape[1]

# usepca = False, good order used to predict bad order
# usepca = True, bad order predicted by the good orders



def pca_trace(xcen, usepca = None, npca = 2, npoly_cen = 3, debug=True):

    nspec = xcen.shape[0]
    norders = xcen.shape[1]
    if usepca is None:
        usepca = np.zeros(norders,dtype=bool)

    # use_order = True orders used to predict the usepca = True bad orders
    use_order = np.invert(usepca)
    ngood = np.sum(use_order)
    if ngood < npca:
        msgs.warn('Not enough good traces for a PCA fit: ngood = {:d}'.format(ngood) + ' is < npca = {:d}'.format(npca))
        msgs.warn('Using the input trace for now')
        return xcen

    pca = PCA(n_components=npca)
    xcen_use = (xcen[:,use_order] - np.mean(xcen[:,use_order],0)).T
    pca_coeffs_use = pca.fit_transform(xcen_use)
    pca_vectors = pca.components_

    # Fit first pca dimension (with largest variance) with a higher order npoly depending on number of good orders.
    # Fit all higher dimensions (with lower variance) with a line
    #npoly_vec = np.ones(npca,dtype=int)
    npoly = int(np.fmin(np.fmax(np.floor(3.3*ngood/norders),1.0),3.0))
    #npoly_vec[0]=npoly
    npoly_vec = np.full(npca, npoly)

    order_vec = np.arange(norders,dtype=float)
    pca_coeffs = np.zeros((norders, npca))
    # Now loop over the dimensionality of the compression and perform a polynomial fit to
    for idim in range(npca):
        # ToDO robust_polyfit is garbage remove it entirely from PypeIT!
        xfit = order_vec[use_order]
        yfit = pca_coeffs_use[:,idim]
        norder = npoly_vec[idim]
        msk, poly_coeff = utils.robust_polyfit(xfit, yfit, norder, sigma = 3.0, function='polynomial')
        # TESTING
        xtemp = xfit.reshape(1, xfit.size)
        ytemp = yfit.reshape(1, yfit.size)
        tset = pydl.xy2traceset(xtemp, ytemp, ncoeff=norder,function='poly')
        tset_yfit = tset.yfit.reshape(tset.yfit.shape[1])
        #pca_coeffs_tset = tset.yfit
        pca_coeffs[:,idim] = utils.func_val(poly_coeff, order_vec, 'polynomial')

        if debug:
            # Evaluate the fit
            xvec = np.linspace(order_vec.min(),order_vec.max(),num=100)
            (_,tset_fit) = tset.xy(xpos=xvec.reshape(1,xvec.size))
            yfit_tset = tset_fit[0,:]
            robust_mask = msk == 0
            tset_mask = tset.outmask[0,:]
            plt.plot(xfit[robust_mask],yfit[robust_mask],'ko', markersize = 8.0, label = 'pca coeff fit')
            plt.plot(xfit[~robust_mask]*1.05,yfit[~robust_mask],'k+', markersize = 10.0, label = 'pca coeff rejected')
            plt.plot(xfit[tset_mask],yfit[tset_mask],'ro', markersize = 8.0, label = 'pca coeff fit')
            plt.plot(xfit[~tset_mask]*1.05,yfit[~tset_mask], 'r+', markersize = 10.0, label = 'pca coeff rejected')
            plt.plot(xvec, utils.func_val(poly_coeff, xvec, 'polynomial'), color='black',label = 'robust polyfit')
            plt.plot(xvec, yfit_tset, color='red',label='pydl')
            plt.legend()
            from IPython import embed
            embed()
            plt.show()

    #ToDo should we be masking the bad orders here and interpolating/extrapolating?
    spat_mean = np.mean(xcen,0)
    msk_spat, poly_coeff_spat = utils.robust_polyfit(order_vec, spat_mean, npoly_cen, sigma = 3.0, function = 'polynomial')
    ibad = np.where(msk_spat == 1)
    spat_mean[ibad] = utils.func_val(poly_coeff_spat,order_vec[ibad],'polynomial')

    pca_fit = np.outer(np.ones(nspec), spat_mean) + np.outer(pca.mean_,np.ones(norders)) + (np.dot(pca_coeffs, pca_vectors)).T

    return pca_fit


def ech_objfind(image, ivar, ordermask, slit_left, slit_righ,inmask=None,plate_scale=0.2,npca=2,ncoeff = 5,min_snr=0.0,nabove_min_snr=0,
                pca_percentile=20.0,snr_pca=3.0,box_radius=2.0,show_peaks=False,show_fits=False,show_trace=False):


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
    show_peaks=True
    show_fits=True
    # ToDo replace orderindx with the true order number here? Maybe not. Clean up slitid and orderindx!
    for iord  in range(norders):
        msgs.info('Finding objects on slit # {:d}'.format(iord + 1))
        thismask = ordermask == (iord + 1)
        inmask_iord = inmask & thismask
        specobj_dict = {'setup': 'HIRES', 'slitid': iord + 1, 'scidx': 0,'det': 1, 'objtype': 'science'}
        sobjs_slit, skymask[thismask], objmask[thismask], proc_list = \
            extract.objfind(image, thismask, slit_left[:,iord], slit_righ[:,iord], inmask=inmask_iord,show_peaks=show_peaks,
                            show_fits=show_fits, show_trace=False, specobj_dict = specobj_dict)
        # ToDO make the specobjs _set_item_ work with expressions like this spec[:].orderindx = iord
        for spec in sobjs_slit:
            spec.ech_orderindx = iord
        sobjs.add_sobj(sobjs_slit)


    nfound = len(sobjs)

    # Compute the FOF linking length based on the instrument place scale and matching length FOFSEP = 1.0"
    FOFSEP = 1.0 # separation of FOF algorithm in arcseconds
    FOF_frac = FOFSEP/(np.median(slit_width)*np.median(plate_scale_ord))

    # Run the FOF. We use fake coordinaes
    fracpos = sobjs.spat_fracpos
    ra_fake = fracpos/1000.0 # Divide all angles by 1000 to make geometry euclidian
    dec_fake = 0.0*fracpos
    (ingroup, multgroup, firstgroup, nextgroup) = spheregroup(ra_fake, dec_fake, FOF_frac/1000.0)
    group = ingroup.copy()
    uni_group, uni_ind = np.unique(group, return_index=True)
    nobj = len(uni_group)
    msgs.info('FOF matching found {:d}'.format(nobj) + ' unique objects')
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
        pca_fits[:,:,iobj] = pca_trace((sobjs_final[igroup].trace_spat).T, usepca = None, npca = npca, npoly_cen = npoly_cen)
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
    #for spec in sobjs:
    #    ginga.show_trace(viewer, ch, spec.trace_spat, spec.idx, color='orange')






#xcen = slit_left        # (nspec, norder)
#inmask = np.ones_like(xcen,dtype=bool) # (nspec, norder) boolean


# HIRES
spectro = 'ESI'
if spectro == 'HIRES':
    hdu = fits.open('/Users/joe/Dropbox/hires_fndobj/f_hires0181G.fits.gz')
    objminsky =hdu[2].data
    ivar  = hdu[1].data
    mask = (ivar > 0.0)
    order_str = Table.read('/Users/joe/Dropbox/hires_fndobj/OStr_G_02.fits')
    slit_left = (order_str['LHEDG']).T
    slit_righ = (order_str['RHEDG']).T
    plate_scale = 0.36
elif spectro == 'ESI':
    # ESI
    hdu = fits.open('/Users/joe/Dropbox/Cowie_2002-02-17/Final/fringe_ES.20020217.35453.fits.gz')
    sciimg = hdu[0].data
    var  = hdu[1].data
    ivar = utils.calc_ivar(var)
    mask = (var > 0.0)
    skyimg = hdu[2].data
    objminsky = sciimg - skyimg
    hdu_sedg = fits.open('/Users/joe/Dropbox/Cowie_2002-02-17/Flats/SEdgECH75_1x1.fits')
    data = hdu_sedg[0].data
    slit_left = data[0,:,:].T
    slit_righ = data[1,:,:].T
    plate_scale = 0.149


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
pca_out = pca_trace(slit_mid, npca = 4, npoly_cen = 3)
sys.exit(-1)
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
min_snr = 0.3
nabove_min_snr = 2
pca_percentile = 20.0
snr_pca = 3.0
npca = 4

sobjs_final = ech_objfind(image, ivar, ordermask, slit_left, slit_righ,inmask=inmask,plate_scale=plate_scale, npca = npca,
                          ncoeff = 5,min_snr=min_snr,nabove_min_snr=nabove_min_snr,pca_percentile=pca_percentile,snr_pca=snr_pca,
                          box_radius=box_radius,show_peaks=True,show_fits=False,show_trace=True)


