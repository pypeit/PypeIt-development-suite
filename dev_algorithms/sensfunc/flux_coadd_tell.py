import os
import numpy as np

from astropy.io import fits

from pypeit.core import telluric, save, load
from pypeit.core.flux_calib import apply_sensfunc
from pypeit.core import coadd1d
from pypeit import msgs

basedir = os.getenv('HOME')
#basedir = '/d2/Feige'

def get_sens_from_file(std1dfile=None, instrument='GNIRS', star_type=None, star_mag=None,star_ra=None,
                       star_dec=None, sens_polyorder=8, mask_abs_lines=True, disp=True, debug=False):


    # sensfunction output file name
    if '.sens.fits' in std1dfile:
        sensfile = std1dfile
    else:
        sensfile = std1dfile.replace('.fits','.sens.fits')

    # get the pca pickle file and atmosphere model grid
    pca_file = os.path.join(basedir, 'Dropbox/PypeIt_Redux/qso_pca_1200_3100.pckl')

    if (instrument=='GNIRS') or (instrument=='NIRES'):
        telgridfile = os.path.join(basedir, 'Dropbox/PypeIt_Redux/TelFit_MaunaKea_3100_26100_R20000.fits')
    elif (instrument == 'XSHOOTER_VIS') or (instrument == 'GMOS-S'):
        telgridfile = os.path.join(basedir,
                                   'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_VIS_4900_11100_R25000.fits')
    elif instrument == 'XSHOOTER_NIR':
        telgridfile = os.path.join(basedir,
                                   'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
    else:
        telgridfile = os.path.join(basedir, 'Dropbox/PypeIt_Redux/TelFit_MaunaKea_3100_26100_R20000.fits')
        msgs.warn('No telluric grid is found. Using MaunaKea!')
    msgs.info('Using {:}'.format(telgridfile))

    # run telluric.sensfunc_telluric to get the sensfile
    TelSens = telluric.sensfunc_telluric(std1dfile, telgridfile, sensfile, star_type=star_type, star_mag=star_mag,
                                         star_ra=star_ra, star_dec=star_dec, mask_abs_lines=mask_abs_lines,
                                         polyorder=sens_polyorder, disp=disp, debug=debug)
    return sensfile, telgridfile

def apply_tell_from_file(z_obj, stackfilename, tell_method='qso', instrument='NIRES', telgridfile=None,
                         polyorder=3, fit_region_min=[9200.0], fit_region_max=[9700.0], mask_lyman_a=True,
                         show=True, debug=False):

    # output names
    if tell_method == 'poly':
        telloutfile = stackfilename.replace('.fits', '_polytellmodel.fits')
        outfile = stackfilename.replace('.fits', '_polytellcorr.fits')
    else:
        telloutfile = stackfilename.replace('.fits','_tellmodel.fits')
        outfile = stackfilename.replace('.fits','_tellcorr.fits')

    if tell_method=='qso':
        pca_file = os.path.join(basedir, 'Dropbox/PypeIt_Redux/qso_pca_1200_3100.pckl')
        # run telluric.qso_telluric to get the final results
        # TODO: add other modes here
        TelQSO = telluric.qso_telluric(stackfilename, telgridfile, pca_file, z_obj, telloutfile, outfile,
                                       create_bal_mask=None, disp=show, debug=debug, show=show)
    elif tell_method=='poly':
        TelPoly = telluric.poly_telluric(stackfilename, telgridfile, telloutfile, outfile, z_obj=z_obj, polyorder=polyorder,
                                         fit_region_min=fit_region_min, fit_region_max=fit_region_max, func='legendre',
                                         model='exp', mask_lyman_a=mask_lyman_a, debug_init=debug, debug=debug, show=show)

def flux_tell(sci_path, stdfile, spec1dfiles=None, std_path=None, fileroot=None, outroot=None, z_qso=None, tell_method='qso',
              instrument=None, star_type=None, star_mag=None, star_ra=None, star_dec=None, mask_abs_lines=True,
              sens_polyorder=8, objids=None, ex_value='OPT', polyorder=3, fit_region_min=[9200.0], fit_region_max=[9700.0],
              scale_method=None, hand_scale=None, const_weights=False, wave_grid_min=None, wave_grid_max=None,
              mask_lyman_a=True, do_sens=False, do_flux=False, do_stack=False, do_tell=False, use_exist_sens=True,
              disp=False, debug=False):

    if std_path is None:
        std_path = sci_path
    if fileroot is None:
        fileroot = 'spec1d_flux_tell.fits'
    if outroot is None:
        outroot = fileroot

    std1dfile = os.path.join(std_path, stdfile)
    header = fits.getheader(std1dfile, 0)
    pypeline = header['PYPELINE']

    if spec1dfiles is None:
        scifiles_all = os.listdir(sci_path)
        spec1dfiles = []
        for i in range(len(scifiles_all)):
            if ('spec1d' in scifiles_all[i]) and (fileroot in scifiles_all[i]):
                spec1dfiles.append(scifiles_all[i])

    nfiles = len(spec1dfiles)
    fnames = []
    for ifile in range(nfiles):
        fnames.append(os.path.join(sci_path, spec1dfiles[ifile]))
    if objids is None:
        objids = ['OBJ0001'] * nfiles
    elif len(objids) != nfiles:
        msgs.error('The length of objids should be exactly the same with the number of spec1d files.')
    #fnames = np.sort(fnames)

    ### get sensfunc
    if (star_ra is None) and (star_dec is None) and (star_mag is None) and (star_type is None):
        star_ra, star_dec = header['RA'], header['DEC']

    if '.sens.fits' not in stdfile:
        sensfile = std1dfile.replace('.fits', '.sens.fits')
    telgridfile = None # value it to None
    if do_sens:
        if os.path.exists(sensfile) and (use_exist_sens):
            msgs.warn('{:} is already exists. Skip doing sensfunc.'.format(sensfile))
        else:
            sensfile, telgridfile = get_sens_from_file(std1dfile=std1dfile, instrument=instrument, star_type=star_type,
                                                       star_mag=star_mag, star_ra=star_ra, star_dec=star_dec,
                                                       sens_polyorder = sens_polyorder,
                                                       mask_abs_lines=mask_abs_lines, disp=disp, debug=debug)
    if telgridfile is None:
        msgs.info('Loading sensfile {:}'.format(sensfile))

        if (instrument=='GNIRS') or (instrument=='NIRES'):
            telgridfile = os.path.join(basedir, 'Dropbox/PypeIt_Redux/TelFit_MaunaKea_3100_26100_R20000.fits')
        elif (instrument == 'XSHOOTER_VIS') or (instrument == 'GMOS-S'):
            telgridfile = os.path.join(basedir,
                                       'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_VIS_4900_11100_R25000.fits')
        elif instrument == 'XSHOOTER_NIR':
            telgridfile = os.path.join(basedir,
                                       'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
        else:
            telgridfile = os.path.join(basedir, 'Dropbox/PypeIt_Redux/TelFit_MaunaKea_3100_26100_R20000.fits')
            msgs.warn('No telluric grid is found. Using MaunaKea!')

    ### Apply the sensfunc to all spectra (only sensfunc but not tellluric)
    if do_flux:
        apply_sensfunc(fnames, sensfile, extinct_correct=False, tell_correct=False, debug=debug, show=disp)
        # fnames_flux = [f.replace('.fits', '_flux.fits') for f in fnames]
    else:
        msgs.warn('You skiped the fluxing step, make sure you have applied sensfunction to your 1D spectra.')

    # The name of the final stacked 1d spectrum
    if len(outroot.split('.')) == 1:
        outroot = outroot+'.fits'
    elif outroot.split('.')[-1] != 'fits':
        outroot = outroot + '.fits'
    stackfile = os.path.join(sci_path, outroot)

    if do_stack:
        ## Let's coadd all the fluxed spectra
        # you should get a coadded spectrum named as '{:}_stack.fits'.format(qsoname)
        #                a straight merge of individual order stacked spectra named as '{:}_merge.fits'.format(qsoname)
        #                a individual order stacked spectra (multi-extension) named as '{:}_order.fits'.format(qsoname)
        if pypeline == 'Echelle':
            wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames, objids, show=disp,
                                                                sensfile=sensfile, ex_value=ex_value, outfile=stackfile,
                                                                scale_method=scale_method, hand_scale=hand_scale,
                                                                wave_grid_min=wave_grid_min,wave_grid_max=wave_grid_max,
                                                                const_weights=const_weights, debug=debug)
        else:
            wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.multi_combspec(fnames, objids, show=disp,
                                                                ex_value=ex_value, outfile=stackfile,
                                                                wave_grid_min=wave_grid_min,wave_grid_max=wave_grid_max,
                                                                scale_method=scale_method, hand_scale=hand_scale,
                                                                const_weights=const_weights, debug=debug, debug_scale=debug)
    elif os.path.exists(stackfile):
        msgs.info('Loading stacked 1d spectrum {:}'.format(stackfile))
    else:
        msgs.warn('No stacked 1d spectrum was found. Please set do_stack=True!')

    ### Telluric correction
    if do_tell:
        apply_tell_from_file(z_qso, stackfile, tell_method=tell_method, instrument=instrument, telgridfile=telgridfile,
                             polyorder=polyorder, fit_region_min=fit_region_min, fit_region_max=fit_region_max,
                             mask_lyman_a=mask_lyman_a, show=disp, debug=debug)

def stack_multinight(sci_path,fileroot, outroot=None, spec1dfiles=None, objids=None, wave_method='log10', ex_value='OPT',
                     wave_grid_min=None, wave_grid_max=None, dwave=None, dv=None, dloglam=None, samp_fact=1.0,
                     scale_method='poly', hand_scale=None, const_weights=False, ref_percentile=80.0, sn_smooth_npix=None,
                     debug=False, show=False):

    if spec1dfiles is None:
        #spec1dfiles = np.genfromtxt(spec1dlist,dtype='str')
        scifiles_all = os.listdir(sci_path)
        spec1dfiles = []
        for i in range(len(scifiles_all)):
            if ('tellcorr' in scifiles_all[i]) and (fileroot in scifiles_all[i]):
                spec1dfiles.append(scifiles_all[i])

    nfiles = len(spec1dfiles)
    fnames = []
    nspec_array = []
    for ifile in range(nfiles):
        fname = os.path.join(sci_path,spec1dfiles[ifile])
        fnames.append(fname)
        nspeci = fits.getheader(fname,1)['NAXIS2']
        nspec_array.append(nspeci)
    nspec = np.max(nspec_array)
    hdulist = fits.open(fnames[0])
    header = hdulist[0].header

    if objids is None:
        objids = ['OBJ0001-SPEC0001-OPT']*nfiles

    waves = np.zeros((nspec, nfiles))
    fluxes, ivars, masks = np.zeros_like(waves), np.zeros_like(waves), np.zeros_like(waves, dtype=bool)
    for ifile in range(nfiles):
        hdulist = fits.open(fnames[ifile])
        ext_id = objids[ifile]
        wave, flux, ivar, mask = load.load_ext_to_array(hdulist, ext_id, ex_value=ex_value,
                                                        flux_value=True, nmaskedge=None)
        waves[:len(wave),ifile], fluxes[:len(wave),ifile] = wave, flux
        ivars[:len(wave), ifile], masks[:len(wave), ifile] = ivar, mask

    # Decide how much to smooth the spectra by if this number was not passed in
    if sn_smooth_npix is None:
        nspec, nexp = waves.shape
        # This is the effective good number of spectral pixels in the stack
        nspec_eff = np.sum(waves > 1.0)/nexp
        sn_smooth_npix = int(np.round(0.1*nspec_eff))
        msgs.info('Using a sn_smooth_npix={:d} to decide how to scale and weight your spectra'.format(sn_smooth_npix))

    wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.combspec(waves, fluxes, ivars, masks, sn_smooth_npix,
             wave_method=wave_method, scale_method=scale_method, const_weights=const_weights,ref_percentile=ref_percentile,
             wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max, dwave=dwave, dv=dv, dloglam=dloglam, samp_fact=samp_fact,
             hand_scale=hand_scale, debug=debug, show=show, debug_scale=debug, show_scale=show)

    if outroot is None:
        outroot = fileroot.copy()
    if len(outroot.split('.')) == 1:
        outroot = outroot + '.fits'
    elif outroot.split('.')[-1] != 'fits':
        outroot = outroot + '.fits'
    outfile = os.path.join(sci_path, outroot)

    save.save_coadd1d_to_fits(outfile, wave_stack, flux_stack, ivar_stack, mask_stack, header=header,
                              ex_value=ex_value, overwrite=True)

def merge_vis_nir(outfile, spec1dvis, spec1dnir, sci_path='./',stack_region = [10150.0,10200.0],
                  wave_method = 'log10', dwave=None, dv=None, dloglam=None, samp_fact=1.0, wave_grid_min=None,
                  wave_grid_max=None, const_weights=True, sn_smooth_npix=None, ref_percentile=20.0,
                  maxiter_scale=5, sigrej_scale=3, scale_method='none', hand_scale=None, sn_max_medscale=2.0,
                  sn_min_medscale=0.5, sn_clip=30.0, lower=3.0, upper=3.0, maxrej=None, maxiter_reject=5,
                  qafile=None, title='', debug=False, show=True):

    ## Read VIS
    hdu_vis = fits.open(os.path.join(sci_path,spec1dvis))
    wave_vis, flux_vis, ivar_vis, mask_vis = load.load_ext_to_array(hdu_vis, 'OBJ0001-SPEC0001-OPT',
                                                                    ex_value='OPT', flux_value=True, nmaskedge=None)

    ## Read NIR
    hdu_nir = fits.open(os.path.join(sci_path,spec1dnir))
    wave_nir, flux_nir, ivar_nir, mask_nir = load.load_ext_to_array(hdu_nir, 'OBJ0001-SPEC0001-OPT',
                                                                    ex_value='OPT', flux_value=True, nmaskedge=None)

    ## combine the VIS and NIR in the overlap region
    # data preparation
    mask_stack1 = (wave_vis>stack_region[0]) & (wave_vis<stack_region[1])
    nstack1 = np.sum(mask_stack1)
    mask_stack2 = (wave_nir>stack_region[0]) & (wave_nir<stack_region[1])
    nstack2 = np.sum(mask_stack2)

    waves = np.zeros((np.max([nstack1,nstack2]), 2))
    fluxes, ivars, masks = np.zeros_like(waves), np.zeros_like(waves), np.zeros_like(waves, dtype=bool)

    waves[:nstack1,0], fluxes[:nstack1,0] = wave_vis[mask_stack1], flux_vis[mask_stack1]
    ivars[:nstack1,0], masks[:nstack1,0] = ivar_vis[mask_stack1], mask_vis[mask_stack1]

    waves[:nstack2,1], fluxes[:nstack2,1] = wave_nir[mask_stack2], flux_nir[mask_stack2]
    ivars[:nstack2,1], masks[:nstack2,1] = ivar_nir[mask_stack2], mask_nir[mask_stack2]

    ## Start doing the combine in the overlap region
    # Decide how much to smooth the spectra by if this number was not passed in
    if sn_smooth_npix is None:
        nspec, nexp = waves.shape
        # This is the effective good number of spectral pixels in the stack
        nspec_eff = np.sum(waves > 1.0)/nexp
        sn_smooth_npix = int(np.round(0.1*nspec_eff))
        msgs.info('Using a sn_smooth_npix={:d} to decide how to scale and weight your spectra'.format(sn_smooth_npix))

    # merge VIS and NIR
    wave_grid, _, _ = coadd1d.get_wave_grid(waves, masks=masks, wave_method=wave_method, wave_grid_min=wave_grid_min,
                                        wave_grid_max=wave_grid_max,dwave=dwave, dv=dv, dloglam=dloglam, samp_fact=samp_fact)

    # Evaluate the sn_weights. This is done once at the beginning
    rms_sn, weights = coadd1d.sn_weights(waves, fluxes, ivars, masks, sn_smooth_npix, const_weights=const_weights, verbose=True)

    fluxes_scale, ivars_scale, scales, scale_method_used = coadd1d.scale_spec_stack(
            wave_grid, waves, fluxes, ivars, masks, rms_sn, weights, ref_percentile=ref_percentile, maxiter_scale=maxiter_scale,
            sigrej_scale=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale,
            sn_max_medscale=sn_max_medscale, sn_min_medscale=sn_min_medscale, debug=debug, show=show)
    scale_blue = np.median(scales[:,0])
    scale_red = np.median(scales[:,1])

    # Rejecting and coadding
    wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused = coadd1d.spec_reject_comb(
            wave_grid, waves, fluxes_scale, ivars_scale, masks, weights, sn_clip=sn_clip, lower=lower, upper=upper,
            maxrej=maxrej, maxiter_reject=maxiter_reject, debug=debug, title=title)

    ## Merge blue, overlap stacked spectrum and the red
    wave = np.hstack([wave_vis[wave_vis<stack_region[0]],wave_stack,wave_nir[wave_nir>stack_region[1]]])
    flux = np.hstack([scale_blue*flux_vis[wave_vis<stack_region[0]], flux_stack,
                      scale_red*flux_nir[wave_nir>stack_region[1]]])
    ivar = np.hstack([1.0/scale_blue**2*ivar_vis[wave_vis<stack_region[0]], ivar_stack,
                      1.0/scale_red**2*ivar_nir[wave_nir>stack_region[1]]])
    mask = np.hstack([mask_vis[wave_vis<stack_region[0]]>0, mask_stack,
                      mask_nir[wave_nir>stack_region[1]]>0])
    nused_vis = np.ones_like(wave_vis,dtype=int)
    nused_nir = np.ones_like(wave_nir,dtype=int)
    nused_all = np.hstack([nused_vis[wave_vis<stack_region[0]],nused,nused_nir[wave_nir>stack_region[1]]])
    if show:
        coadd1d.coadd_qa(wave, flux, ivar, nused_all, mask=mask, title='Final merged spectrum', qafile=qafile)

    save.save_coadd1d_to_fits(os.path.join(sci_path,outfile), wave, flux, ivar, mask,
                              header=hdu_vis[0].header, ex_value='OPT', overwrite=True)
