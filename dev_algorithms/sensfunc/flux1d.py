import os, glob
import scipy
import numpy as np
from pkg_resources import resource_filename
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy import units
from astropy import constants
from astropy import coordinates
from astropy.table import Table, Column

from pypeit import utils
from pypeit import msgs
from pypeit.core import load, save, coadd1d
from pypeit.spectrographs.util import load_spectrograph
from IPython import embed

PYPEIT_FLUX_SCALE = 1e-17

def find_standard_file(ra, dec, toler=20.*units.arcmin, check=False):
    """
    Find a match for the input file to one of the archived
    standard star files (hopefully).  Priority is by order of search.

    Args:
        ra (str):
            Object right-ascension in hh:mm:ss string format (e.g.,
            '05:06:36.6').
        dec (str):
            Object declination in dd:mm:ss string format (e.g.,
            52:52:01.0')
        toler (:class:`astropy.units.quantity.Quantity`, optional):
            Tolerance on matching archived standards to input.  Expected
            to be in arcmin.
        check (:obj:`bool`, optional):
            If True, the routine will only check to see if a standard
            star exists within the input ra, dec, and toler range.

    Returns:
        dict, bool: If check is True, return True or False depending on
        if the object is matched to a library standard star.  If check
        is False and no match is found, return None.  Otherwise, return
        a dictionary with the matching standard star with the following
        meta data::
            - 'file': str -- Filename
              table
            - 'name': str -- Star name
            - 'ra': str -- RA(J2000)
            - 'dec': str -- DEC(J2000)
    """
    # Priority
    std_sets = ['xshooter', 'calspec']

    # SkyCoord
    if ':' in ra:
        obj_coord = coordinates.SkyCoord(ra, dec, unit=(units.hourangle, units.deg))
    else:
        obj_coord = coordinates.SkyCoord(ra, dec, unit=(units.deg, units.deg))

    # Loop on standard sets
    closest = dict(sep=999 * units.deg)

    for sset in std_sets:
        path = 'data/standards/{:}/'.format(sset)
        star_file = resource_filename('pypeit', path+'{:}_info.txt'.format(sset, sset))
        star_tbl = Table.read(star_file, comment='#', format='ascii')

        star_coords = coordinates.SkyCoord(star_tbl['RA_2000'], star_tbl['DEC_2000'],
                                           unit=(units.hourangle, units.deg))
        idx, d2d, d3d = coordinates.match_coordinates_sky(obj_coord, star_coords, nthneighbor=1)

        if d2d < toler:
            if check:
                return True
            else:
                # Generate a dict
                _idx = int(idx)
                std_dict = dict(cal_file=os.path.join(path,star_tbl[_idx]['File']),
                                name=star_tbl[_idx]['Name'],
                                std_ra=star_tbl[_idx]['RA_2000'],
                                std_dec=star_tbl[_idx]['DEC_2000'])
                if os.path.exists(star_file):
                    root = resource_filename('pypeit', std_dict['cal_file'] + '*')
                    fil = glob.glob(root)
                    if len(fil) == 0:
                        msgs.error("No standard star file: {:s}".format(fil))
                    else:
                        fil = fil[0]
                        msgs.info("Loading standard star file: {:s}".format(fil))
                    if sset == 'xshooter':
                        std_spec = Table.read(fil, format='ascii')
                        std_dict['wave'] = std_spec['col1'] * units.AA
                        std_dict['flux'] = std_spec['col2'] / PYPEIT_FLUX_SCALE * \
                                           units.erg / units.s / units.cm ** 2 / units.AA
                        msgs.info("Fluxes are flambda, normalized to 1e-17")
                    elif sset == 'calspec':
                        std_spec = fits.open(fil)[1].data
                        std_dict['wave'] = std_spec['WAVELENGTH'] * units.AA
                        std_dict['flux'] = std_spec['FLUX'] / PYPEIT_FLUX_SCALE \
                                           * units.erg / units.s / units.cm ** 2 / units.AA
                        msgs.info("Fluxes are flambda, normalized to 1e-17")
                else:
                    msgs.error("No standard star file found: {:s}".format(star_file))

                return std_dict
        else:
            # Save closest found so far
            imind2d = np.argmin(d2d)
            mind2d = d2d[imind2d]
            if mind2d < closest['sep']:
                closest['sep'] = mind2d
                closest.update(dict(name=star_tbl[int(idx)]['Name'],
                                    ra=star_tbl[int(idx)]['RA_2000'],
                                    dec=star_tbl[int(idx)]['DEC_2000']))
    # Standard star not found
    if check:
        return False

    msgs.error("No standard star was found within a tolerance of {:g}".format(toler) + msgs.newline()
               + "Closest standard was {:s} at separation {:g}".format(closest['name'], closest['sep'].to('arcmin')))

    return None

def stellar_model(V, sptype):
    """Parse Kurucz SED given T and g
    Also convert absolute/apparent magnitudes

    Parameters:
    ----------
    V: float
      Apparent magnitude of the telluric star
    sptype: str
      Spectral type of the telluric star

    Returns:
    ----------
    loglam: ndarray
      log wavelengths
    flux: ndarray
      SED f_lambda (cgs units, I think, probably per Ang)
    """

    # Grab telluric star parameters
    # log(g) of the Sun
    logg_sol = np.log10(6.67259e-8) + np.log10(1.989e33) - 2.0 * np.log10(6.96e10)

    # Load Schmidt-Kaler (1982) table
    sk82_file = resource_filename('pypeit', 'data/standards/kurucz93/schmidt-kaler_table.txt')
    sk82_tab = ascii.read(sk82_file, names=('Sp', 'logTeff', 'Teff', '(B-V)_0', 'M_V', 'B.C.', 'M_bol', 'L/L_sol'))

    # Match input type
    mti = np.where(sptype == sk82_tab['Sp'])[0]
    if len(mti) != 1:
        raise ValueError('Not ready to interpolate yet.')

    # Calculate final quantities
    # Relation between radius, temp, and bolometric luminosity
    logR = 0.2 * (42.26 - sk82_tab['M_bol'][mti[0]] - 10.0 * sk82_tab['logTeff'][mti[0]])

    # Mass-bolometric luminosity relation from schimdt-kaler p28 valid for M_bol < 7.5
    logM = 0.46 - 0.10 * sk82_tab['M_bol'][mti[0]]
    logg = logM - 2.0 * logR + logg_sol
    M_V = sk82_tab['M_V'][mti[0]]
    Teff = sk82_tab['Teff'][mti[0]]

    # Flux factor (absolute/apparent V mag)
    # Define constants
    parsec = constants.pc.cgs  # 3.086e18
    R_sol = constants.R_sun.cgs  # 6.96e10

    # Distance modulus
    logd = 0.2 * (V - M_V) + 1.0
    D = parsec * 10. ** logd
    R = R_sol * 10. ** logR

    # Factor converts the kurucz surface flux densities to flux observed on Earth
    flux_factor = (R / D.value) ** 2

    # Grab closest T in Kurucz SEDs
    T1 = 3000. + np.arange(28) * 250
    T2 = 10000. + np.arange(6) * 500
    T3 = 13000. + np.arange(22) * 1000
    T4 = 35000. + np.arange(7) * 2500
    Tk = np.concatenate([T1, T2, T3, T4])
    indT = np.argmin(np.abs(Tk - Teff))

    # Grab closest g in Kurucz SEDs
    loggk = np.arange(11) * 0.5
    indg = np.argmin(np.abs(loggk - logg))

    # Grab Kurucz filename
    std_file = resource_filename('pypeit', '/data/standards/kurucz93/kp00/kp00_{:d}.fits.gz'.format(int(Tk[indT])))
    std = Table.read(std_file)

    # Grab specific spectrum
    loglam = np.array(np.log10(std['WAVELENGTH']))
    gdict = {0: 'g00', 1: 'g05', 2: 'g10', 3: 'g15', 4: 'g20',
             5: 'g25', 6: 'g30', 7: 'g35', 8: 'g40', 9: 'g45',
             10: 'g50'}
    flux = std[gdict[indg]]

    # scale the model to the V-band magnitude
    star_lam = 10 ** loglam
    star_flux = flux.data * flux_factor
    # Generate a dict matching the output of find_standard_file
    std_dict = dict(cal_file='KuruczTelluricModel', name=sptype, Vmag=V, std_ra=None, std_dec=None)
    std_dict['wave'] = star_lam * units.AA
    std_dict['flux'] = star_flux / PYPEIT_FLUX_SCALE * units.erg / units.s / units.cm ** 2 / units.AA

    return std_dict

def get_standard_spectrum(star_type=None, star_mag=None, ra=None, dec=None):
    '''
    Get the standard spetrum using given information of your standard/telluric star.

    Parameters:
      star_type: str
         Spectral type of your standard/telluric star
      star_mag: float
       Apparent magnitude of the telluric star
      ra: str
        Standard right-ascension in hh:mm:ss string format (e.g.,'05:06:36.6').
      dec: str
        Object declination in dd:mm:ss string format (e.g., 52:52:01.0')
    Return: dict
        Dictionary containing the information you provided and the standard/telluric spectrum.
    '''
    # Create star model
    if (ra is not None) and (dec is not None) and (star_mag is None) and (star_type is None):
        # Pull star spectral model from archive
        msgs.info("Getting archival standard spectrum")
        # Grab closest standard within a tolerance
        std_dict = find_standard_file(ra, dec)

    elif (star_mag is not None) and (star_type is not None):
        ## using vega spectrum
        if 'A0' in star_type:
            msgs.info('Getting vega spectrum')
            std_dict={'stellar_type':star_type , 'Vmag': star_mag}

            ## Vega model from TSPECTOOL
            vega_file = resource_filename('pypeit', '/data/standards/vega_tspectool_vacuum.dat')
            vega_data = Table.read(vega_file, comment='#', format='ascii')
            std_dict = dict(cal_file='vega_tspectool_vacuum', name=star_type, std_ra=None, std_dec=None)
            std_dict['wave'] = vega_data['col1'] * units.AA

            # vega is V=0.03
            std_dict['flux'] = vega_data['col2'] * 10**(0.4*(0.03-star_mag)) / PYPEIT_FLUX_SCALE * \
                               units.erg / units.s / units.cm ** 2 / units.AA
        ## using Kurucz stellar model
        else:
            # Create star spectral model
            msgs.info("Getting kurucz+93 stellar model")
            std_dict = stellar_model(star_mag, star_type)
    else:
        msgs.error('Insufficient information provided for fluxing. '
                   'Either the coordinates of the standard or a stellar type and magnitude are needed.')

    return std_dict


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

def apply_sensfunc_specobjs(specobjs, sens_meta, sens_table, airmass, exptime, extinct_correct=True, tell_correct=False,
                            longitude=None, latitude=None, debug=False, show=False):

    # TODO This function should operate on a single object
    func = sens_meta['FUNC'][0]
    polyorder_vec = sens_meta['POLYORDER_VEC'][0]
    nimgs = len(specobjs)

    if show:
        fig = plt.figure(figsize=(12, 8))
        xmin, xmax = [], []
        ymin, ymax = [], []

    for ispec in range(nimgs):
        # get the ECH_ORDER, ECH_ORDERINDX, WAVELENGTH from your science
        sobj_ispec = specobjs[ispec]
        ## TODO Comment on the logich here. Hard to follow
        try:
            ech_order, ech_orderindx, idx = sobj_ispec.ech_order, sobj_ispec.ech_orderindx, sobj_ispec.idx
            msgs.info('Applying sensfunc to Echelle data')
        except:
            ech_orderindx = 0
            idx = sobj_ispec.idx
            msgs.info('Applying sensfunc to Longslit/Multislit data')

        for extract_type in ['boxcar', 'optimal']:
            extract = getattr(sobj_ispec, extract_type)

            if len(extract) == 0:
                continue
            msgs.info("Fluxing {:s} extraction for:".format(extract_type) + msgs.newline() + "{}".format(idx))
            wave = extract['WAVE'].value.copy()
            wave_mask = wave > 1.0
            counts = extract['COUNTS'].copy()
            counts_ivar = extract['COUNTS_IVAR'].copy()
            mask = extract['MASK'].copy()

            # get sensfunc from the sens_table
            coeff = sens_table[ech_orderindx]['OBJ_THETA'][0:polyorder_vec[ech_orderindx] + 2]
            wave_min = sens_table[ech_orderindx]['WAVE_MIN']
            wave_max = sens_table[ech_orderindx]['WAVE_MAX']
            sensfunc = np.zeros_like(wave)
            sensfunc[wave_mask] = np.exp(utils.func_val(coeff, wave[wave_mask], func,
                                             minx=wave_min, maxx=wave_max))

            # get telluric from the sens_table
            if tell_correct:
                msgs.work('Evaluate telluric!')
                telluric = None
            else:
                telluric = None

            flam, flam_ivar, outmask = apply_sensfunc_spec(wave, counts, counts_ivar, sensfunc, airmass, exptime,
                                                      mask=mask, extinct_correct=extinct_correct, telluric=telluric,
                                                      longitude=longitude, latitude=latitude, debug=debug)
            flam_sig = np.sqrt(utils.inverse(flam_ivar))
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

                med_width = (2.0 * np.ceil(0.1 / 10.0 * np.size(wave[outmask])) + 1).astype(int)
                flam_med, flam_ivar_med = coadd1d.median_filt_spec(flam, flam_ivar, outmask, med_width)
                if extract_type == 'boxcar':
                    plt.plot(wave[wave_mask], flam_med[wave_mask], color='black', drawstyle='steps-mid', zorder=1, alpha=0.8)
                    #plt.plot(wave[wave_mask], np.sqrt(utils.calc_ivar(flam_ivar_med[wave_mask])), zorder=2, color='m',
                    #         alpha=0.7, drawstyle='steps-mid', linestyle=':')
                else:
                    plt.plot(wave[wave_mask], flam_med[wave_mask], color='dodgerblue', drawstyle='steps-mid', zorder=1, alpha=0.8)
                    #plt.plot(wave[wave_mask], np.sqrt(utils.calc_ivar(flam_ivar_med[wave_mask])), zorder=2, color='red',
                    #         alpha=0.7, drawstyle='steps-mid', linestyle=':')
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

    sens_meta = Table.read(sensfile, 1)
    sens_table = Table.read(sensfile, 2)

    nexp = np.size(fnames)
    for iexp in range(nexp):
        spec1dfile = fnames[iexp]
        outfile = spec1dfile[:-5] + '_flux.fits'
        sobjs, head = load.load_specobjs(spec1dfile)
        instrument = head['INSTRUME']
        spectrograph = load_spectrograph(instrument)
        airmass, exptime = head['AIRMASS'], head['EXPTIME']
        longitude, latitude = head['LON-OBS'], head['LAT-OBS']

        apply_sensfunc_specobjs(sobjs, sens_meta, sens_table, airmass, exptime, extinct_correct=extinct_correct,
                                tell_correct=tell_correct, longitude=longitude, latitude=latitude,
                                debug=debug, show=show)
        save.save_1d_spectra_fits(sobjs, head, spectrograph, outfile, helio_dict=None, overwrite=True)


## Standard star
#datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Science/')
#fnames = [datapath + 'spec1d_XSHOO.2018-11-08T00:11:57.074-Feige110_XShooter_NIR_2018Nov08T001157.074.fits',
#          datapath + 'spec1d_XSHOO.2018-11-08T00:16:56.583-Feige110_XShooter_NIR_2018Nov08T001656.583.fits']

## Lensed Quasar
'''
datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0439_old/vlt_xshooter_nir/Science/')
fnames = [datapath + 'J0439_XSHOOTER_NIR_01.fits', datapath + 'J0439_XSHOOTER_NIR_02.fits',
          datapath + 'J0439_XSHOOTER_NIR_03.fits', datapath + 'J0439_XSHOOTER_NIR_04.fits',
          datapath + 'J0439_XSHOOTER_NIR_05.fits', datapath + 'J0439_XSHOOTER_NIR_06.fits',
          datapath + 'J0439_XSHOOTER_NIR_07.fits', datapath + 'J0439_XSHOOTER_NIR_08.fits',
          datapath + 'J0439_XSHOOTER_NIR_09.fits', datapath + 'J0439_XSHOOTER_NIR_10.fits',
          datapath + 'J0439_XSHOOTER_NIR_11.fits', datapath + 'J0439_XSHOOTER_NIR_12.fits',
          datapath + 'J0439_XSHOOTER_NIR_13.fits', datapath + 'J0439_XSHOOTER_NIR_14.fits',
          datapath + 'J0439_XSHOOTER_NIR_15.fits', datapath + 'J0439_XSHOOTER_NIR_16.fits',
          datapath + 'J0439_XSHOOTER_NIR_17.fits', datapath + 'J0439_XSHOOTER_NIR_18.fits',
          datapath + 'J0439_XSHOOTER_NIR_19.fits', datapath + 'J0439_XSHOOTER_NIR_20.fits',
          datapath + 'J0439_XSHOOTER_NIR_21.fits', datapath + 'J0439_XSHOOTER_NIR_22.fits',
          datapath + 'J0439_XSHOOTER_NIR_23.fits', datapath + 'J0439_XSHOOTER_NIR_24.fits',
          datapath + 'J0439_XSHOOTER_NIR_25.fits', datapath + 'J0439_XSHOOTER_NIR_26.fits',
          datapath + 'J0439_XSHOOTER_NIR_27.fits', datapath + 'J0439_XSHOOTER_NIR_28.fits',
          datapath + 'J0439_XSHOOTER_NIR_29.fits', datapath + 'J0439_XSHOOTER_NIR_30.fits',
          datapath + 'J0439_XSHOOTER_NIR_31.fits', datapath + 'J0439_XSHOOTER_NIR_32.fits',
          datapath + 'J0439_XSHOOTER_NIR_33.fits', datapath + 'J0439_XSHOOTER_NIR_34.fits',
          datapath + 'J0439_XSHOOTER_NIR_35.fits', datapath + 'J0439_XSHOOTER_NIR_36.fits',
          datapath + 'J0439_XSHOOTER_NIR_37.fits', datapath + 'J0439_XSHOOTER_NIR_38.fits']
'''

'''
## Pisco
datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/Pypeit_files/PISCO_NIR/Science/')
fnames = [datapath + 'spec1d_XSHOO.2017-06-28T23:51:39.115-PSOJ205p09_1_XShooter_NIR_2017Jun28T235139.115.fits',
          datapath + 'spec1d_XSHOO.2017-06-29T00:12:24.992-PSOJ205p09_1_XShooter_NIR_2017Jun29T001224.992.fits',
          datapath + 'spec1d_XSHOO.2018-02-15T08:07:00.164-PSOJ205p09_2_XShooter_NIR_2018Feb15T080700.164.fits',
          datapath + 'spec1d_XSHOO.2018-02-15T08:27:44.713-PSOJ205p09_2_XShooter_NIR_2018Feb15T082744.713.fits',
          datapath + 'spec1d_XSHOO.2018-03-10T05:08:40.055-PSOJ205p09_3_XShooter_NIR_2018Mar10T050840.055.fits',
          datapath + 'spec1d_XSHOO.2018-03-10T05:29:25.933-PSOJ205p09_3_XShooter_NIR_2018Mar10T052925.933.fits',
          datapath + 'spec1d_XSHOO.2018-03-10T05:58:33.093-PSOJ205p09_4_XShooter_NIR_2018Mar10T055833.093.fits',
          datapath + 'spec1d_XSHOO.2018-03-10T06:19:22.960-PSOJ205p09_4_XShooter_NIR_2018Mar10T061922.960.fits',
          datapath + 'spec1d_XSHOO.2018-03-11T08:41:01.429-PSOJ205p09_5_XShooter_NIR_2018Mar11T084101.429.fits',
          datapath + 'spec1d_XSHOO.2018-03-11T09:01:47.309-PSOJ205p09_5_XShooter_NIR_2018Mar11T090147.309.fits',
          datapath + 'spec1d_XSHOO.2018-03-18T06:26:29.402-PSOJ205p09_6_XShooter_NIR_2018Mar18T062629.402.fits',
          datapath + 'spec1d_XSHOO.2018-03-18T06:47:15.945-PSOJ205p09_6_XShooter_NIR_2018Mar18T064715.945.fits',
          datapath + 'spec1d_XSHOO.2018-03-18T08:32:32.302-PSOJ205p09_7_XShooter_NIR_2018Mar18T083232.302.fits',
          datapath + 'spec1d_XSHOO.2018-03-18T08:53:18.845-PSOJ205p09_7_XShooter_NIR_2018Mar18T085318.845.fits',
          datapath + 'spec1d_XSHOO.2018-03-20T05:22:01.367-PSOJ205p09_8_XShooter_NIR_2018Mar20T052201.367.fits',
          datapath + 'spec1d_XSHOO.2018-03-20T05:42:48.575-PSOJ205p09_8_XShooter_NIR_2018Mar20T054248.575.fits',
          datapath + 'spec1d_XSHOO.2018-03-21T06:06:09.751-PSOJ205p09_9_XShooter_NIR_2018Mar21T060609.751.fits',
          datapath + 'spec1d_XSHOO.2018-03-21T06:26:55.631-PSOJ205p09_9_XShooter_NIR_2018Mar21T062655.631.fits',
          datapath + 'spec1d_XSHOO.2018-03-21T07:13:08.797-PSOJ205p09_10_XShooter_NIR_2018Mar21T071308.797.fits',
          datapath + 'spec1d_XSHOO.2018-03-21T07:33:56.008-PSOJ205p09_10_XShooter_NIR_2018Mar21T073356.008.fits',
          datapath + 'spec1d_XSHOO.2018-03-22T05:00:30.229-PSOJ205p09_11_XShooter_NIR_2018Mar22T050030.229.fits',
          datapath + 'spec1d_XSHOO.2018-03-22T05:21:17.439-PSOJ205p09_11_XShooter_NIR_2018Mar22T052117.439.fits',
          datapath + 'spec1d_XSHOO.2018-03-22T06:01:30.465-PSOJ205p09_12_XShooter_NIR_2018Mar22T060130.465.fits',
          datapath + 'spec1d_XSHOO.2018-03-22T06:22:17.675-PSOJ205p09_12_XShooter_NIR_2018Mar22T062217.675.fits',
          datapath + 'spec1d_XSHOO.2018-03-23T04:45:05.305-PSOJ205p09_13_XShooter_NIR_2018Mar23T044505.305.fits',
          datapath + 'spec1d_XSHOO.2018-03-23T05:05:50.520-PSOJ205p09_13_XShooter_NIR_2018Mar23T050550.520.fits',
          datapath + 'spec1d_XSHOO.2018-03-23T05:43:50.910-PSOJ205p09_14_XShooter_NIR_2018Mar23T054350.910.fits',
          datapath + 'spec1d_XSHOO.2018-03-23T06:04:37.452-PSOJ205p09_14_XShooter_NIR_2018Mar23T060437.452.fits',
          datapath + 'spec1d_XSHOO.2018-03-23T06:37:00.347-PSOJ205p09_15_XShooter_NIR_2018Mar23T063700.347.fits',
          datapath + 'spec1d_XSHOO.2018-03-23T06:57:45.561-PSOJ205p09_16_XShooter_NIR_2018Mar23T065745.561.fits',
          datapath + 'spec1d_XSHOO.2018-03-23T07:36:36.220-PSOJ205p09_16_XShooter_NIR_2018Mar23T073636.220.fits',
          datapath + 'spec1d_XSHOO.2018-03-23T07:57:23.428-PSOJ205p09_16_XShooter_NIR_2018Mar23T075723.428.fits',
          datapath + 'spec1d_XSHOO.2018-03-24T05:20:10.198-PSOJ205p09_17_XShooter_NIR_2018Mar24T052010.198.fits',
          datapath + 'spec1d_XSHOO.2018-03-24T05:40:56.740-PSOJ205p09_17_XShooter_NIR_2018Mar24T054056.740.fits',
          datapath + 'spec1d_XSHOO.2018-03-24T06:17:18.332-PSOJ205p09_18_XShooter_NIR_2018Mar24T061718.332.fits',
          datapath + 'spec1d_XSHOO.2018-03-24T06:38:05.542-PSOJ205p09_18_XShooter_NIR_2018Mar24T063805.542.fits',
          datapath + 'spec1d_XSHOO.2018-03-24T07:13:35.465-PSOJ205p09_19_XShooter_NIR_2018Mar24T071335.465.fits',
          datapath + 'spec1d_XSHOO.2018-03-24T07:34:22.671-PSOJ205p09_19_XShooter_NIR_2018Mar24T073422.671.fits',
          datapath + 'spec1d_XSHOO.2018-03-25T05:18:51.790-PSOJ205p09_20_XShooter_NIR_2018Mar25T051851.790.fits',
          datapath + 'spec1d_XSHOO.2018-03-25T05:39:38.333-PSOJ205p09_20_XShooter_NIR_2018Mar25T053938.333.fits',
          datapath + 'spec1d_XSHOO.2018-03-26T06:36:32.911-PSOJ205p09_21_XShooter_NIR_2018Mar26T063632.911.fits',
          datapath + 'spec1d_XSHOO.2018-03-26T06:57:20.118-PSOJ205p09_21_XShooter_NIR_2018Mar26T065720.118.fits',
          datapath + 'spec1d_XSHOO.2018-03-27T06:57:29.269-PSOJ205p09_22_XShooter_NIR_2018Mar27T065729.269.fits',
          datapath + 'spec1d_XSHOO.2018-03-27T07:18:15.146-PSOJ205p09_22_XShooter_NIR_2018Mar27T071815.146.fits',
          datapath + 'spec1d_XSHOO.2018-04-08T04:59:57.666-PSOJ205p09_23_XShooter_NIR_2018Apr08T045957.666.fits',
          datapath + 'spec1d_XSHOO.2018-04-08T05:20:44.206-PSOJ205p09_23_XShooter_NIR_2018Apr08T052044.206.fits',
          datapath + 'spec1d_XSHOO.2018-04-08T05:57:54.622-PSOJ205p09_24_XShooter_NIR_2018Apr08T055754.622.fits',
          datapath + 'spec1d_XSHOO.2018-04-08T06:18:42.492-PSOJ205p09_24_XShooter_NIR_2018Apr08T061842.492.fits',
          datapath + 'spec1d_XSHOO.2018-04-09T04:59:20.273-PSOJ205p09_25_XShooter_NIR_2018Apr09T045920.273.fits',
          datapath + 'spec1d_XSHOO.2018-04-09T05:20:06.817-PSOJ205p09_25_XShooter_NIR_2018Apr09T052006.817.fits',
          datapath + 'spec1d_XSHOO.2018-04-09T05:58:52.879-PSOJ205p09_26_XShooter_NIR_2018Apr09T055852.879.fits',
          datapath + 'spec1d_XSHOO.2018-04-09T06:19:39.420-PSOJ205p09_26_XShooter_NIR_2018Apr09T061939.420.fits',
          datapath + 'spec1d_XSHOO.2018-04-10T05:07:58.959-PSOJ205p09_27_XShooter_NIR_2018Apr10T050758.959.fits',
          datapath + 'spec1d_XSHOO.2018-04-10T05:28:46.165-PSOJ205p09_27_XShooter_NIR_2018Apr10T052846.165.fits',
          datapath + 'spec1d_XSHOO.2018-04-10T06:05:42.360-PSOJ205p09_28_XShooter_NIR_2018Apr10T060542.360.fits',
          datapath + 'spec1d_XSHOO.2018-04-10T06:26:30.231-PSOJ205p09_28_XShooter_NIR_2018Apr10T062630.231.fits',
          datapath + 'spec1d_XSHOO.2018-04-15T04:58:27.291-PSOJ205p09_29_XShooter_NIR_2018Apr15T045827.291.fits',
          datapath + 'spec1d_XSHOO.2018-04-15T05:19:15.167-PSOJ205p09_29_XShooter_NIR_2018Apr15T051915.167.fits',
          datapath + 'spec1d_XSHOO.2018-04-15T05:59:39.364-PSOJ205p09_30_XShooter_NIR_2018Apr15T055939.364.fits',
          datapath + 'spec1d_XSHOO.2018-04-15T06:20:25.908-PSOJ205p09_30_XShooter_NIR_2018Apr15T062025.908.fits',
          datapath + 'spec1d_XSHOO.2018-04-16T04:23:25.898-PSOJ205p09_31_XShooter_NIR_2018Apr16T042325.898.fits',
          datapath + 'spec1d_XSHOO.2018-04-16T04:44:12.441-PSOJ205p09_31_XShooter_NIR_2018Apr16T044412.441.fits',
          datapath + 'spec1d_XSHOO.2018-04-16T05:26:58.617-PSOJ205p09_32_XShooter_NIR_2018Apr16T052658.617.fits',
          datapath + 'spec1d_XSHOO.2018-04-16T05:47:45.160-PSOJ205p09_32_XShooter_NIR_2018Apr16T054745.160.fits']
'''
## DES z~6.9 Quasar
datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0020-3653/NIR/Science/')
fnames = [datapath + 'spec1d_XSHOO.2017-10-26T00:26:41.612-VHSJ0020-3653_XShooter_NIR_2017Oct26T002641.612.fits',
          datapath + 'spec1d_XSHOO.2017-10-26T00:48:55.512-VHSJ0020-3653_XShooter_NIR_2017Oct26T004855.512.fits',
          datapath + 'spec1d_XSHOO.2017-12-17T02:44:43.537-VHSJ0020-3653_XShooter_NIR_2017Dec17T024443.537.fits',
          datapath + 'spec1d_XSHOO.2017-12-17T03:05:50.032-VHSJ0020-3653_XShooter_NIR_2017Dec17T030550.032.fits']

## DES z~6.5 Quasar J0224
#datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0224-4711/pypeit_nir/Science/')
#fnames = [
#    datapath + 'spec1d_XSHOO.2017-11-23T06:52:51.782-VDESJ0224-4711blindoffset_XShooter_NIR_2017Nov23T065251.782.fits',
#    datapath + 'spec1d_XSHOO.2017-11-23T07:13:18.374-VDESJ0224-4711blindoffset_XShooter_NIR_2017Nov23T071318.374.fits',
#    datapath + 'spec1d_XSHOO.2018-01-19T01:57:51.708-VDESJ0224-4711blindoffset_XShooter_NIR_2018Jan19T015751.708.fits',
#    datapath + 'spec1d_XSHOO.2018-01-19T02:18:18.297-VDESJ0224-4711blindoffset_XShooter_NIR_2018Jan19T021818.297.fits']

##J0226
#datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0226+0302/pypeit_nir/Science/')
#fnames = [datapath + 'spec1d_XSHOO.2017-12-17T03:50:47.125-PSOJ036.5078blindoffset_XShooter_NIR_2017Dec17T035047.125.fits',
#          datapath + 'spec1d_XSHOO.2017-12-17T04:11:13.716-PSOJ036.5078blindoffset_XShooter_NIR_2017Dec17T041113.716.fits',
#          datapath + 'spec1d_XSHOO.2018-01-14T02:12:34.014-PSOJ036.5078blindoffset_XShooter_NIR_2018Jan14T021234.014.fits',
#          datapath + 'spec1d_XSHOO.2018-01-14T02:33:00.603-PSOJ036.5078blindoffset_XShooter_NIR_2018Jan14T023300.603.fits']

## Standard star
#datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0020-3653/NIR/Science/')
#fnames = [datapath + 'spec1d_STD,FLUX_XShooter_NIR_2017Dec17T081653.582.fits',
#          datapath + 'spec1d_STD,FLUX_XShooter_NIR_2017Dec17T082243.751.fits']

## Telluric star
#datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0224-4711/Test_tell/')
#fnames = [datapath + 'spec1d_XSHOO.2017-11-23T07:44:02.747-STD,TELLURIC_XShooter_NIR_2017Nov23T074402.747.fits',
#          datapath + 'spec1d_XSHOO.2017-11-23T07:44:57.633-STD,TELLURIC_XShooter_NIR_2017Nov23T074457.633.fits',
#          datapath + 'spec1d_XSHOO.2017-11-23T07:46:53.532-STD,TELLURIC_XShooter_NIR_2017Nov23T074653.532.fits',
#          datapath + 'spec1d_XSHOO.2017-11-23T07:47:33.917-STD,TELLURIC_XShooter_NIR_2017Nov23T074733.917.fits']

##J1048
#datapath = os.path.join(os.getenv('HOME'), 'Dropbox/OBS_DATA/XSHOOTER/NIR/ut20170202/Science/')
#fnames = [datapath + 'spec1d_XSHOO.2017-02-02T04:19:28.545-VIKJ1048m0109_XShooter_NIR_2017Feb02T041928.545.fits',
#          datapath + 'spec1d_XSHOO.2017-02-02T04:40:14.422-VIKJ1048m0109_XShooter_NIR_2017Feb02T044014.422.fits',
#          datapath + 'spec1d_XSHOO.2017-02-02T05:17:52.162-VIKJ1048m0109_XShooter_NIR_2017Feb02T051752.162.fits']

#sensfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/NIR_Stack/Feige110_sens_tell_wang.fits')
sensfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/LTT3218_sens_tell.fits')
apply_sensfunc(fnames, sensfile, extinct_correct=False, tell_correct=False, debug=False, show=False)

'''
## test find standard
star_mag  = None
star_type = None
spec1dfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Science/spec1d_XSHOO.2018-11-08T00:16:56.583-Feige110_XShooter_NIR_2018Nov08T001656.583.fits')
header = fits.getheader(spec1dfile)
ra=header['RA']
dec=header['DEC']
std_dict = get_standard_spectrum(star_type=star_type, star_mag=star_mag, ra=ra, dec=dec)

### Grid
from astropy.io import fits
telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
#from telluric import read_telluric_grid
#tell_dict = read_telluric_grid(tellgridfile, wave_min=None, wave_max=None, pad = 0)

wave_min = None
wave_max = None
pad = 0

hdul = fits.open(telgridfile)
wave_grid_full = 10.0 * hdul[1].data
model_grid_full = hdul[0].data
nspec_full = wave_grid_full.size

if wave_min is not None:
    ind_lower = np.argmin(np.abs(wave_grid_full - wave_min)) - pad
else:
    ind_lower = 0
if wave_max is not None:
    ind_upper = np.argmin(np.abs(wave_grid_full - wave_max)) + pad
else:
    ind_upper = nspec_full
wave_grid = wave_grid_full[ind_lower:ind_upper]
model_grid = model_grid_full[:, :, :, :, ind_lower:ind_upper]
'''
