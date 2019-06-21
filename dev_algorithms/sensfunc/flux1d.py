import os
import scipy
import numpy as np
from pkg_resources import resource_filename
import matplotlib.pyplot as plt

from astropy import units
from astropy import constants
from astropy import coordinates
from astropy.table import Table, Column

from pypeit import utils
from pypeit import msgs
from pypeit.core import load, save, coadd1d
from pypeit.spectrographs.util import load_spectrograph

def load_extinction_data(longitude, latitude, toler=5. * units.deg):
    """
    Find the best extinction file to use, based on longitude and latitude
    Loads it and returns a Table

    Parameters
    ----------
    toler : Angle, optional
      Tolerance for matching detector to site (5 deg)

    Returns
    -------
    ext_file : Table
      astropy Table containing the 'wavelength', 'extinct' data for AM=1.
    """
    # Mosaic coord
    mosaic_coord = coordinates.SkyCoord(longitude, latitude, frame='gcrs', unit=units.deg)
    # Read list
    extinct_path = resource_filename('pypeit', '/data/extinction/')
    extinct_summ = extinct_path + 'README'
    extinct_files = Table.read(extinct_summ, comment='#', format='ascii')
    # Coords
    ext_coord = coordinates.SkyCoord(extinct_files['Lon'], extinct_files['Lat'], frame='gcrs',
                                     unit=units.deg)
    # Match
    idx, d2d, d3d = coordinates.match_coordinates_sky(mosaic_coord, ext_coord, nthneighbor=1)
    if d2d < toler:
        extinct_file = extinct_files[int(idx)]['File']
        msgs.info("Using {:s} for extinction corrections.".format(extinct_file))
    else:
        msgs.warn("No file found for extinction corrections.  Applying none")
        msgs.warn("You should generate a site-specific file")
        return None
    # Read
    extinct = Table.read(extinct_path + extinct_file, comment='#', format='ascii',
                         names=('iwave', 'mag_ext'))
    wave = Column(np.array(extinct['iwave']) * units.AA, name='wave')
    extinct.add_column(wave)
    # Return
    return extinct[['wave', 'mag_ext']]


def extinction_correction(wave, airmass, extinct):
    """
    Derive extinction correction
    Based on algorithm in LowRedux (long_extinct)

    Parameters
    ----------
    wave : ndarray
      Wavelengths for interpolation. Should be sorted
      Assumes Angstroms
    airmass : float
      Airmass
    extinct : Table
      Table of extinction values

    Returns
    -------
    flux_corr : ndarray
      Flux corrections at the input wavelengths
    """
    # Checks
    if airmass < 1.:
        msgs.error("Bad airmass value in extinction_correction")
    # Interpolate
    f_mag_ext = scipy.interpolate.interp1d(extinct['wave'],extinct['mag_ext'], bounds_error=False, fill_value=0.)
    mag_ext = f_mag_ext(wave)#.to('AA').value)

    # Deal with outside wavelengths
    gdv = np.where(mag_ext > 0.)[0]

    if len(gdv) == 0:
        msgs.warn("No valid extinction data available at this wavelength range. Extinction correction not applied")
    elif gdv[0] != 0:  # Low wavelengths
        mag_ext[0:gdv[0]] = mag_ext[gdv[0]]
        msgs.warn("Extrapolating at low wavelengths using last valid value")
    elif gdv[-1] != (mag_ext.size - 1):  # High wavelengths
        mag_ext[gdv[-1] + 1:] = mag_ext[gdv[-1]]
        msgs.warn("Extrapolating at high wavelengths using last valid value")
    else:
        msgs.info("Extinction data covered the whole spectra. Correct it!")
    # Evaluate
    flux_corr = 10.0 ** (0.4 * mag_ext * airmass)
    # Return
    return flux_corr

def apply_sensfunc_spec(wave, counts, ivar, sensfunc, airmass, exptime, mask=None, extinct_correct=True, telluric=None,
                   longitude=None, latitude=None, debug=False):

    if mask is None:
        mask = ivar > 0.0

    # Did the user request a telluric correction from the same file?
    if telluric is not None:
        # This assumes there is a separate telluric key in this dict.
        msgs.info('Applying telluric correction')
        sensfunc = sensfunc*(telluric > 1e-10)/(telluric + (telluric < 1e-10))

    if extinct_correct:
        if longitude is None or latitude is None:
            msgs.error('You must specify longitude and latitude if we are extinction correcting')
        # Apply Extinction if optical bands
        msgs.info("Applying extinction correction")
        msgs.warn("Extinction correction applyed only if the spectra covers <10000Ang.")
        extinct = load_extinction_data(longitude,latitude)
        ext_corr = extinction_correction(wave* units.AA, airmass, extinct)
        senstot = sensfunc * ext_corr
    else:
        senstot = sensfunc.copy()

    flam = counts * senstot/ exptime
    flam_ivar = ivar / (senstot / exptime) **2

    # Mask bad pixels
    msgs.info(" Masking bad pixels")
    outmask =  mask & (senstot>0.)

    # debug
    if debug:
        wave_mask = wave > 1.0
        fig = plt.figure(figsize=(12, 8))
        ymin, ymax = coadd1d.get_ylim(flam, flam_ivar, outmask)
        plt.plot(wave[wave_mask], flam[wave_mask], color='black', drawstyle='steps-mid', zorder=1, alpha=0.8)
        plt.plot(wave[wave_mask], np.sqrt(utils.calc_ivar(flam_ivar[wave_mask])), zorder=2, color='red', alpha=0.7,
                       drawstyle='steps-mid', linestyle=':')
        plt.ylim([ymin,ymax])
        plt.xlim([wave[wave_mask].min(),wave[wave_mask].max()])
        plt.xlabel('Wavelength (Angstrom)')
        plt.ylabel('Flux')
        plt.show()

    return flam, flam_ivar, outmask

def apply_sensfunc_specobjs(specobjs, sens_table, func, airmass, exptime, extinct_correct=True, tell_correct=False,
                            longitude=None, latitude=None, debug=False, show=False):

    nspec = len(specobjs)

    if show:
        fig = plt.figure(figsize=(12, 8))
        xmin, xmax = [], []
        ymin, ymax = [], []

    for ispec in range(nspec):
        # get the ECH_ORDER, ECH_ORDERINDX, WAVELENGTH from your science
        sobj_ispec = specobjs[ispec]
        try:
            iord, iord_indx, iord_idx = sobj_ispec.ech_order, sobj_ispec.ech_orderindx, sobj_ispec.idx
            msgs.info('Applying sensfunc to Echelle data')
        except:
            iord_idx = 0
            msgs.info('Applying sensfunc to Longslit/Multislit data')

        for extract_type in ['boxcar', 'optimal']:
            extract = getattr(sobj_ispec, extract_type)

            if len(extract) == 0:
                continue
            msgs.info("Fluxing {:s} extraction for:".format(extract_type) + msgs.newline() + "{}".format(sobj_ispec.idx))
            wave = np.copy(np.array(extract['WAVE']))
            wave_mask = wave>1.0
            counts = np.copy(np.array(extract['COUNTS']))
            counts_ivar = np.copy(np.array(extract['COUNTS_IVAR']))
            mask = np.copy(np.array(extract['MASK']))

            # get sensfunc from the sens_table
            sens_table_ispec = sens_table[iord_indx]
            coeff = sens_table_ispec['SENS_COEFF']
            sensfunc = np.zeros_like(wave)
            sensfunc[wave_mask] = np.exp(utils.func_val(coeff, wave[wave_mask], func,
                                             minx=wave[wave_mask].min(), maxx=wave[wave_mask].max()))

            # get telluric from the sens_table
            if tell_correct:
                msgs.work('Evaluate telluric!')
                telluric = None
            else:
                telluric = None

            flam, flam_ivar, outmask = apply_sensfunc_spec(wave, counts, counts_ivar, sensfunc, airmass, exptime,
                                                      mask=mask, extinct_correct=extinct_correct, telluric=telluric,
                                                      longitude=longitude, latitude=latitude, debug=debug)
            flam_sig = np.sqrt(utils.calc_ivar(flam_ivar))
            # The following will be changed directly in the specobjs, so do not need to return anything.
            extract['MASK'] = outmask
            extract['FLAM'] = flam
            extract['FLAM_SIG'] = flam_sig
            extract['FLAM_IVAR'] = flam_ivar

            if show:
                xmin_ispec = wave[wave_mask].min()
                xmax_ispec = wave[wave_mask].max()
                xmin.append(xmin_ispec)
                xmax.append(xmax_ispec)
                ymin_ispec, ymax_ispec = coadd1d.get_ylim(flam, flam_ivar, outmask)
                ymin.append(ymin_ispec)
                ymax.append(ymax_ispec)
                if extract_type == 'boxcar':
                    plt.plot(wave[wave_mask], flam[wave_mask], color='black', drawstyle='steps-mid', zorder=1, alpha=0.8)
                    plt.plot(wave[wave_mask], np.sqrt(utils.calc_ivar(flam_ivar[wave_mask])), zorder=2, color='m',
                             alpha=0.7, drawstyle='steps-mid', linestyle=':')
                else:
                    plt.plot(wave[wave_mask], flam[wave_mask], color='blue', drawstyle='steps-mid', zorder=1, alpha=0.8)
                    plt.plot(wave[wave_mask], np.sqrt(utils.calc_ivar(flam_ivar[wave_mask])), zorder=2, color='red',
                             alpha=0.7, drawstyle='steps-mid', linestyle=':')
    if show:
        xmin_final, xmax_final = np.min(xmin), np.max(xmax)
        ymax_final = 1.3*np.median(ymax)
        ymin_final = -0.15*ymax_final
        plt.xlim([xmin_final, xmax_final])
        plt.ylim([ymin_final, ymax_final])
        plt.title('Blue is Optimal extraction and Black is Boxcar extraction',fontsize=16)
        plt.xlabel('Wavelength (Angstrom)')
        plt.ylabel('Flux')
        plt.show()

def apply_sensfunc(fnames, sensfile, extinct_correct=True, tell_correct=False, debug=False, show=False):

    sens_param = Table.read(sensfile, 1)
    sens_table = Table.read(sensfile, 2)
    func = sens_param['FUNCTION']

    nexp = np.size(fnames)
    for iexp in range(nexp):
        spec1dfile = fnames[iexp]
        outfile = spec1dfile[:-5] + '_flux.fits'
        sobjs, head = load.load_specobjs(spec1dfile)
        instrument = head['INSTRUME']
        spectrograph = load_spectrograph(instrument)

        airmass, exptime = head['AIRMASS'], head['EXPTIME']
        longitude, latitude = head['LON-OBS'], head['LAT-OBS']

        apply_sensfunc_specobjs(sobjs, sens_table, func, airmass, exptime, extinct_correct=extinct_correct,
                                tell_correct=tell_correct, longitude=longitude, latitude=latitude,
                                debug=debug, show=show)
        save.save_1d_spectra_fits(sobjs, head, spectrograph, outfile, helio_dict=None, overwrite=True)


## Standard star
#datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Science/')
#fnames = [datapath + 'spec1d_XSHOO.2018-11-08T00:11:57.074-Feige110_XShooter_NIR_2018Nov08T001157.074.fits',
#          datapath + 'spec1d_XSHOO.2018-11-08T00:16:56.583-Feige110_XShooter_NIR_2018Nov08T001656.583.fits']

## Lensed Quasar
datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0439_old/vlt_xshooter_nir/Science/')
fnames = [datapath + 'J0439_XSHOOTER_NIR_01.fits', datapath + 'J0439_XSHOOTER_NIR_02.fits',
          datapath + 'J0439_XSHOOTER_NIR_03.fits',
          datapath + 'J0439_XSHOOTER_NIR_04.fits', datapath + 'J0439_XSHOOTER_NIR_05.fits',
          datapath + 'J0439_XSHOOTER_NIR_06.fits',
          datapath + 'J0439_XSHOOTER_NIR_07.fits', datapath + 'J0439_XSHOOTER_NIR_08.fits',
          datapath + 'J0439_XSHOOTER_NIR_09.fits']

## DES z~6.9 Quasar
#datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0020m3653/NIR/Science/')
#fnames = [datapath + 'spec1d_VHSJ0020-3653OffsetstarB_XShooter_NIR_2017Dec17T024443.537.fits',
#          datapath + 'spec1d_VHSJ0020-3653OffsetstarB_XShooter_NIR_2017Dec17T030550.032.fits',
#          datapath + 'spec1d_VHSJ0020-3653OffsetstarB_XShooter_NIR_2017Oct26T001535.660.fits',
#          datapath + 'spec1d_VHSJ0020-3653OffsetstarB_XShooter_NIR_2017Oct26T002641.612.fits']

## DES z~6.5 Quasar
#datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0224-4711/pypeit_nir/Science/')
#fnames = [
#    datapath + 'spec1d_XSHOO.2017-11-23T06:52:51.782-VDESJ0224-4711blindoffset_XShooter_NIR_2017Nov23T065251.782.fits',
#    datapath + 'spec1d_XSHOO.2017-11-23T07:13:18.374-VDESJ0224-4711blindoffset_XShooter_NIR_2017Nov23T071318.374.fits',
#    datapath + 'spec1d_XSHOO.2018-01-19T01:57:51.708-VDESJ0224-4711blindoffset_XShooter_NIR_2018Jan19T015751.708.fits',
#    datapath + 'spec1d_XSHOO.2018-01-19T02:18:18.297-VDESJ0224-4711blindoffset_XShooter_NIR_2018Jan19T021818.297.fits']

sensfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Feige110_sens_tell.fits')
apply_sensfunc(fnames, sensfile, extinct_correct=True, tell_correct=False, debug=False, show=True)

