import os, glob
import scipy
import numpy as np
from pkg_resources import resource_filename
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.io import ascii
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
    try:
        ra, dec = float(ra), float(dec)
        obj_coord = coordinates.SkyCoord(ra, dec, unit=(units.deg, units.deg))
    except:
        obj_coord = coordinates.SkyCoord(ra, dec, unit=(units.hourangle, units.deg))

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
                        # TODO let's add the star_mag here and get a uniform set of tags in the std_dict
                        std_spec = Table.read(fil, format='ascii')
                        std_dict['std_source'] = sset
                        std_dict['wave'] = std_spec['col1'] * units.AA
                        std_dict['flux'] = std_spec['col2'] / PYPEIT_FLUX_SCALE * \
                                           units.erg / units.s / units.cm ** 2 / units.AA
                        msgs.info("Fluxes are flambda, normalized to 1e-17")
                    elif sset == 'calspec':
                        std_dict['std_source'] = sset
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
    std_dict['std_source'] = 'KuruczTelluricModel'
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
            ## Vega model from TSPECTOOL
            vega_file = resource_filename('pypeit', '/data/standards/vega_tspectool_vacuum.dat')
            vega_data = Table.read(vega_file, comment='#', format='ascii')
            std_dict = dict(cal_file='vega_tspectool_vacuum', name=star_type, Vmag=star_mag, std_ra=ra, std_dec=dec)
            std_dict['std_source'] = 'VEGA'
            std_dict['wave'] = vega_data['col1'] * units.AA

            # vega is V=0.03
            std_dict['flux'] = vega_data['col2'] * 10**(0.4*(0.03-star_mag)) / PYPEIT_FLUX_SCALE * \
                               units.erg / units.s / units.cm ** 2 / units.AA
        ## using Kurucz stellar model
        else:
            # Create star spectral model
            msgs.info("Getting kurucz+93 stellar model")
            std_dict = stellar_model(star_mag, star_type)
            std_dict['std_ra'] = ra
            std_dict['std_dec'] = dec
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

