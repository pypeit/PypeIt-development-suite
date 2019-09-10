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
    elif instrument=='XSHOOTER':
        telgridfile = os.path.join(os.getenv('HOME'),
                                   'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')

    # run telluric.sensfunc_telluric to get the sensfile
    TelSens = telluric.sensfunc_telluric(std1dfile, telgridfile, sensfile, star_type=star_type, star_mag=star_mag,
                                         star_ra=star_ra, star_dec=star_dec, mask_abs_lines=mask_abs_lines,
                                         disp=disp, debug=debug)
    return sensfile

def apply_sens_from_file(spec1dfiles, datapath, sensfile, outfile, objids=None, tell_correct=False,
                         ex_value='OPT', debug=False, show=False):
    #spec1dfiles = np.genfromtxt(spec1dlist,dtype='str')
    nfiles = len(spec1dfiles)
    fnames = []
    for ifile in range(nfiles):
        fnames.append(os.path.join(datapath,spec1dfiles[ifile]))
    if objids is None:
        objids = ['OBJ0001']*nfiles

    ## Apply the sensfunc to all spectra (only sensfunc but not tellluric)
    apply_sensfunc(fnames, sensfile, extinct_correct=False, tell_correct=tell_correct, debug=debug, show=show)

    #fnames_flux = [f.replace('.fits', '_flux.fits') for f in fnames]

    ## Let's coadd all the fluxed spectra
    # you should get a coadded spectrum named as 'spec1d_stack_{:}.fits'.format(qsoname)
    #                a straight merge of individual order stacked spectra named as 'spec1d_merge_{:}.fits'.format(qsoname)
    #                a individual order stacked spectra (multi-extension) named as 'spec1d_order_{:}.fits'.format(qsoname)
    # TODO: change the outfile to work with datapath. It's a hard coding on these names in coadd1d
    wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames, objids, show=show, sensfile=sensfile,
                                                                          ex_value=ex_value, outfile=outfile, debug=debug)

def apply_qsotell_from_file(z_qso, fileroot, instrument='NIRES', show=True, debug=False):

    pca_file = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/qso_pca_1200_3100.pckl')

    if (instrument=='GNIRS') or (instrument=='NIRES'):
        telgridfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/TelFit_MaunaKea_3100_26100_R20000.fits')
    elif instrument=='XSHOOTER':
        telgridfile = os.path.join(os.getenv('HOME'),
                                   'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')

    # run telluric.qso_telluric to get the final results
    spec1dfluxfile = 'spec1d_stack_{:}.fits'.format(fileroot)
    telloutfile = 'spec1d_stack_{:}_tellmodel.fits'.format(fileroot)
    outfile = 'spec1d_stack_{:}_tellcorr.fits'.format(fileroot)

    # TODO: add other modes here
    TelQSO = telluric.qso_telluric(spec1dfluxfile, telgridfile, pca_file, z_qso, telloutfile, outfile,
                                   create_bal_mask=None, debug=debug, show=show)

def stack_multinight(spec1dfiles, datapath, outfile, objids=None, wave_method='log10', ex_value='OPT',
                     sn_smooth_npix=None, debug=False, show=False):

    #spec1dfiles = np.genfromtxt(spec1dlist,dtype='str')
    nfiles = len(spec1dfiles)
    fnames = []
    nspec_array = []
    for ifile in range(nfiles):
        fname = os.path.join(datapath,spec1dfiles[ifile])
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

    if outfile is not None:
        save.save_coadd1d_to_fits(outfile, wave_stack, flux_stack, ivar_stack, mask_stack, header=header,
                                  ex_value=ex_value, overwrite=True)

def flux_tell(sci_path, stdfile, objname=None, outroot=None, z_qso=None, tell_method='qso', instrument=None,
              star_type=None, star_mag=None, star_ra=None, star_dec=None, mask_abs_lines=True,
              objids=None, disp=False, debug=False):

    # get sensfunc
    std1dfile = os.path.join(sci_path, stdfile)
    if (star_ra is None) and (star_dec is None) and (star_mag is None) and (star_type is None):
        header = fits.getheader(std1dfile,0)
        star_ra, star_dec = header['RA'], header['DEC']

    from IPython import embed
    embed()

    sensfile = get_sens_from_file(std1dfile=std1dfile, instrument=instrument, star_type=star_type, star_mag=star_mag,
                                  star_ra=star_ra, star_dec=star_dec, mask_abs_lines=mask_abs_lines,
                                  disp=disp, debug=debug)
    scifiles_all = os.listdir(sci_path)
    spec1dfiles = []
    for i in range(len(scifiles_all)):
        if ('spec1d' in scifiles_all[i]) and (objname in scifiles_all[i]):
            spec1dfiles.append(scifiles_all[i])

    # apply sensfunc
    apply_sens_from_file(spec1dfiles, sci_path, sensfile, outroot, objids=objids, tell_correct=False, debug=debug,
                         show=disp)
    apply_qsotell_from_file(z_qso, outroot, instrument='NIRES', show=disp, debug=debug)


sci_path = '/Users/feige/Dropbox/OBS_DATA/NIRES/ut190218/Science'
stdfile = 'spec1d_s190218_0078-HD101060_NIRES_2019Feb18T105456.083.fits'
objname = 'J1129'
outroot = 'J1129_NIRES'
instrument = 'NIRES'
star_ra, star_dec = '05:05:30.61', '+52:49:51.9'
star_type, star_mag = None, None

#spec1dlist = os.path.join(datapath, 'spec1dlist_J1129')
#sensfile = os.path.join(datapath, 'spec1d_s190218_0030-G191B2B_NIRES_2019Feb18T050628.633.sens.fits')

flux_tell(sci_path, stdfile, objname=objname, outroot=outroot, instrument=instrument, star_type=star_type,
              star_mag=star_mag, star_ra=star_ra, star_dec=star_dec, mask_abs_lines=True, disp=False, debug=False)