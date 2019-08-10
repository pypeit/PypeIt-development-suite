import os
import numpy as np

from astropy.io import fits

from pypeit.core import telluric, save, load
from pypeit.core.flux_calib import apply_sensfunc
from pypeit.core import coadd1d
from pypeit import msgs


def get_sens_from_file(std1dfile=None, instrument='GNIRS', star_type=None, star_mag=None,star_ra=None,
                       star_dec=None, mask_abs_lines=True, disp=True, debug=False):
    # sensfunction output file name
    sensfile = std1dfile.replace('.fits','.sens.fits')

    # get the pca pickle file and atmosphere model grid
    pca_file = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/qso_pca_1200_3100.pckl')

    if (instrument=='GNIRS') or (instrument=='NIRES'):
        telgridfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/TelFit_MaunaKea_3100_26100_R20000.fits')
    elif instrument == 'XSHOOTER_VIS':
        telgridfile = os.path.join(os.getenv('HOME'),
                                   'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_VIS_4900_11100_R25000.fits')
    elif instrument == 'XSHOOTER_NIR':
        telgridfile = os.path.join(os.getenv('HOME'),
                                   'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
    else:
        msgs.error('Telluric Grid TBD.')

    # run telluric.sensfunc_telluric to get the sensfile
    TelSens = telluric.sensfunc_telluric(std1dfile, telgridfile, sensfile, star_type=star_type, star_mag=star_mag,
                                         star_ra=star_ra, star_dec=star_dec, mask_abs_lines=mask_abs_lines,
                                         disp=disp, debug=debug)
    return sensfile, telgridfile

def apply_tell_from_file(z_obj, stackfilename, tell_method='qso', instrument='NIRES', telgridfile=None,
                         polyorder=3, fit_region_min=[9200.0], fit_region_max=[9700.0], mask_lyman_a=True,
                         show=True, debug=False):

    # output names
    telloutfile = stackfilename.replace('.fits','_tellmodel.fits')
    outfile = stackfilename.replace('.fits','_tellcorr.fits')

    if tell_method=='qso':
        pca_file = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/qso_pca_1200_3100.pckl')
        # run telluric.qso_telluric to get the final results
        # TODO: add other modes here
        TelQSO = telluric.qso_telluric(stackfilename, telgridfile, pca_file, z_obj, telloutfile, outfile,
                                       create_bal_mask=None, disp=show, debug=debug, show=show)
    elif tell_method=='poly':
        TelPoly = telluric.poly_telluric(stackfilename, telgridfile, telloutfile, outfile, z_obj=z_obj, polyorder=polyorder,
                                         fit_region_min=fit_region_min, fit_region_max=fit_region_max, func='legendre',
                                         model='exp', mask_lyman_a=mask_lyman_a, debug_init=debug, debug=debug, show=show)

def stack_multinight(sci_path, spec1dfiles,fileroot, objids=None, wave_method='log10', ex_value='OPT',
                     sn_smooth_npix=None, debug=False, show=False):

    #spec1dfiles = np.genfromtxt(spec1dlist,dtype='str')
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
             wave_method=wave_method, debug=debug, show=show)

    if fileroot is not None:
        if len(fileroot.split('.')) == 1:
            fileroot = fileroot + '.fits'
        elif fileroot.split('.')[-1] != 'fits':
            fileroot = fileroot + '.fits'
        outfile = os.path.join(sci_path, fileroot)

        save.save_coadd1d_to_fits(outfile, wave_stack, flux_stack, ivar_stack, mask_stack, header=header,
                                  ex_value=ex_value, overwrite=True)

def flux_tell(sci_path, stdfile, fileroot=None, z_qso=None, tell_method='qso', instrument=None,
              star_type=None, star_mag=None, star_ra=None, star_dec=None, mask_abs_lines=True,
              objids=None, ex_value='OPT', polyorder=3, fit_region_min=[9200.0], fit_region_max=[9850.0],
              mask_lyman_a=True, do_sens=True, do_flux=True, do_stack=True, do_tell=True, disp=False, debug=False):

    if fileroot is None:
        fileroot = 'spec1d_flux_tell.fits'

    ### get sensfunc
    std1dfile = os.path.join(sci_path, stdfile)
    if (star_ra is None) and (star_dec is None) and (star_mag is None) and (star_type is None):
        header = fits.getheader(std1dfile,0)
        star_ra, star_dec = header['RA'], header['DEC']

    if do_sens:
        sensfile, telgridfile = get_sens_from_file(std1dfile=std1dfile, instrument=instrument, star_type=star_type,
                                                   star_mag=star_mag, star_ra=star_ra, star_dec=star_dec,
                                                   mask_abs_lines=mask_abs_lines, disp=disp, debug=debug)
    else:
        sensfile = std1dfile.replace('.fits', '.sens.fits')
        msgs.info('Loading sensfile {:}'.format(sensfile))

        if (instrument=='GNIRS') or (instrument=='NIRES'):
            telgridfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/TelFit_MaunaKea_3100_26100_R20000.fits')
        elif instrument == 'XSHOOTER_VIS':
            telgridfile = os.path.join(os.getenv('HOME'),
                                       'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_VIS_4900_11100_R25000.fits')
        elif instrument == 'XSHOOTER_NIR':
            telgridfile = os.path.join(os.getenv('HOME'),
                                       'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
        else:
            msgs.error('Telluric Grid TBD.')

    ### Apply the sensfunc to all spectra (only sensfunc but not tellluric)
    if do_flux:
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
            objids = ['OBJ0001']*nfiles
        elif len(objids) != nfiles:
            msgs.error('The length of objids should be exactly the same with the number of spec1d files.')

        apply_sensfunc(fnames, sensfile, extinct_correct=False, tell_correct=False, debug=debug, show=disp)
        # fnames_flux = [f.replace('.fits', '_flux.fits') for f in fnames]
    else:
        msgs.warn('You skiped the fluxing step, make sure you have applied sensfunction to your 1D spectra.')

    # The name of the final stacked 1d spectrum
    if len(fileroot.split('.')) == 1:
        fileroot = fileroot+'.fits'
    elif fileroot.split('.')[-1] != 'fits':
        fileroot = fileroot + '.fits'
    stackfile = os.path.join(sci_path, fileroot)

    if do_stack:
        ## Let's coadd all the fluxed spectra
        # you should get a coadded spectrum named as '{:}_stack.fits'.format(qsoname)
        #                a straight merge of individual order stacked spectra named as '{:}_merge.fits'.format(qsoname)
        #                a individual order stacked spectra (multi-extension) named as '{:}_order.fits'.format(qsoname)

        wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames, objids, show=disp, sensfile=sensfile,
                                                                              ex_value=ex_value, outfile=stackfile, debug=debug)
    elif os.path.exists(stackfile):
        msgs.info('Loading stacked 1d spectrum {:}'.format(stackfile))
    else:
        msgs.warn('No stacked 1d spectrum was found. Please set do_stack=True!')

    ### Telluric correction
    if do_tell:
        apply_tell_from_file(z_qso, stackfile, tell_method=tell_method, instrument=instrument, telgridfile=telgridfile,
                             polyorder=polyorder, fit_region_min=fit_region_min, fit_region_max=fit_region_max,
                             mask_lyman_a=mask_lyman_a, show=disp, debug=debug)

