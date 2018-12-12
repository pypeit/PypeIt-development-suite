# Module for echelle fluxing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import linetools
import os
import json
import matplotlib.pyplot as plt

from astropy import units
from astropy.io import fits
from astropy.table import Table

from pypeit import msgs
from pypeit.core import flux
from pypeit.core import load
from pypeit.core import save
from pypeit import utils
from pypeit import masterframe
from pypeit import specobjs

from pypeit.spectrographs.util import load_spectrograph
from pypeit.par.pypeitpar import TelescopePar

from pypeit import debugger

def ech_load_specobj(fname,order=None):

    """ Load a spec1d file into a list of SpecObjExp objects
    Parameters
    ----------
    fname : str

    Returns
    -------
    specObjs : list of SpecObjExp
    head0
    """
    #if order is None:
    #    msgs.warn('You did not specify an order. Return specObjs with all orders.')
    #    specObjs, head0 = load.load_specobj(fname)
    #    return specObjs, head0

    speckeys = ['WAVE', 'SKY', 'MASK', 'FLAM', 'FLAM_IVAR', 'FLAM_SIG', 'COUNTS_IVAR', 'COUNTS']
    #
    specObjs = []
    hdulist = fits.open(fname)
    head0 = hdulist[0].header
    for hdu in hdulist:
        if hdu.name == 'PRIMARY':
            continue
        #elif hdu.name[8:17] != 'ORDER'+'{0:04}'.format(order):
        #    continue
        # Parse name
        idx = hdu.name
        objp = idx.split('-')
        if objp[-2][0:3] == 'DET':
            det = int(objp[-2][3:])
        else:
            det = int(objp[-2][1:])
        if objp[-3][:5] == 'ORDER':
            iord = int(objp[-3][5:])
        else:
            msgs.warn('Loading longslit data ?')
            iord = int(-1)
        # if order is not None and iord !=order then do not return this extenction
        # if order is None return all extensions
        # if order is not None and iord ==order then only return the specific order you want.
        if (order is not None) and (iord !=order):
            continue
        # Load data
        spec = Table(hdu.data)
        shape = (len(spec), 1024)  # 2nd number is dummy
        # New and wrong
        try:
            specobj = specobjs.SpecObj(shape, None, None, idx = idx)
        except:
            debugger.set_trace()
            msgs.error("BUG ME")
        # Add order number
        specobj.ech_orderindx = iord
        # ToDo: need to changed to the real order number?
        specobj.ech_order = iord
        # Add trace
        try:
            specobj.trace_spat = spec['TRACE']
        except:
            # KLUDGE!
            specobj.trace_spat = np.arange(len(spec['BOX_WAVE']))
        # Add spectrum
        if 'BOX_COUNTS' in spec.keys():
            for skey in speckeys:
                try:
                    specobj.boxcar[skey] = spec['BOX_{:s}'.format(skey)].data
                except KeyError:
                    pass
            # Add units on wave
            specobj.boxcar['WAVE'] = specobj.boxcar['WAVE'] * units.AA

        if 'OPT_COUNTS' in spec.keys():
            for skey in speckeys:
                try:
                    specobj.optimal[skey] = spec['OPT_{:s}'.format(skey)].data
                except KeyError:
                    pass
            # Add units on wave
            specobj.optimal['WAVE'] = specobj.optimal['WAVE'] * units.AA
        # Append
        specObjs.append(specobj)
    # Return
    return specObjs, head0

def ech_save_master(sens_dicts, outfile=None):
    """
    Over-load the ech_save_master() method in MasterFrame to write a JSON file

    Parameters
    ----------
    outfile : str, optional
      Use this input instead of the 'proper' (or unattainable) MasterFrame name

    Returns
    -------

    """
    # Allow one to over-ride output name
    if outfile is None:
        outfile = 'MasterSensFunc.fits'
        #outfile = self.ms_name

    norder = sens_dicts['norder']
    # Add steps
    #self.sens_dict['steps'] = self.steps
    # Do it
    prihdu = fits.PrimaryHDU()
    hdus = [prihdu]

    for iord in range(norder):
        sens_dict = sens_dicts[str(iord)]
        cols = []
        cols += [fits.Column(array=sens_dict['wave'], name=str('WAVE'), format=sens_dict['wave'].dtype)]
        cols += [
            fits.Column(array=sens_dict['sensfunc'], name=str('SENSFUNC'), format=sens_dict['sensfunc'].dtype)]
        # Finish
        coldefs = fits.ColDefs(cols)
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = 'SENSFUNC-ORDER{0:04}'.format(iord)
        # Add critical keys from sens_dict to header
        for key in ['wave_min', 'wave_max', 'exptime', 'airmass', 'std_file', 'std_ra',
                    'std_dec', 'std_name', 'cal_file', 'ech_orderindx']:
            try:
                tbhdu.header[key.upper()] = sens_dict[key].value
            except AttributeError:
                tbhdu.header[key.upper()] = sens_dict[key]
            except KeyError:
                pass  # Will not require all of these
        hdus += [tbhdu]

    # Add critical keys from sens_dict to primary header
    for key in ['exptime', 'airmass', 'std_file', 'std_ra',
                'std_dec', 'std_name', 'cal_file']:
        try:
            prihdu.header[key.upper()] = sens_dict[key].value
        except AttributeError:
            prihdu.header[key.upper()] = sens_dict[key]
        except KeyError:
            pass  # Will not require all of these
    prihdu.header['NORDER'] = norder

        # Finish
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(outfile, overwrite=True)

    # Finish
    msgs.info("Wrote sensfunc to MasterFrame: {:s}".format(outfile))


def ech_load_master(filename, force=False):

    # Does the master file exist?
    if not os.path.isfile(filename):
        #msgs.warn("No Master frame found of type {:s}: {:s}".format(self.frametype, filename))
        msgs.warn("No Master frame found of {:s}".format(filename))
        if force:
            msgs.error("Crashing out because reduce-masters-force=True:" + msgs.newline() + filename)
        return None
    else:
        #msgs.info("Loading a pre-existing master calibration frame of type: {:}".format(self.frametype) + " from filename: {:}".format(filename))
        msgs.info("Loading a pre-existing master calibration frame of SENSFUNC from filename: {:}".format(filename))

        hdu = fits.open(filename)
        norder = hdu[0].header['NORDER']
        sens_dicts = {}
        for iord in range(norder):
            head = hdu[iord+1].header
            tbl = hdu['SENSFUNC-ORDER{0:04}'.format(iord)].data
            sens_dict = {}
            sens_dict['wave'] = tbl['WAVE']
            sens_dict['sensfunc'] = tbl['SENSFUNC']
            for key in ['wave_min','wave_max','exptime','airmass','std_file','std_ra','std_dec',
                        'std_name','cal_file', 'ech_orderindx']:
                try:
                    sens_dict[key] = head[key.upper()]
                except:
                    pass
            sens_dicts[str(iord)] = sens_dict
        sens_dicts['norder'] = norder
        return sens_dicts

#def find_standard(std_specobjs):
#    std_idx = flux.find_standard(std_specobjs)
#    std = std_specobjs[std_idx]
#    return std

def ech_generate_sensfunc(stdframe,spectrograph=None, telluric=True, star_type=None,
                      star_mag=None, ra=None, dec=None, std_file = None, BALM_MASK_WID=5., nresln=None,debug=False):

    if spectrograph is None:
        std_specobjs, std_header = ech_load_specobj(stdframe, order=0)
        spectrograph = std_header['INSTRUME']
        msgs.info('You are working on {:s}'.format(spectrograph))
    ext_final = fits.getheader(stdframe, -1)
    norder = ext_final['ORDER'] + 1

    sens_dicts = {}
    for iord in range(norder):
        std_specobjs, std_header = ech_load_specobj(stdframe, order=iord)
        std_idx = flux.find_standard(std_specobjs)
        std = std_specobjs[std_idx]
        sens_dict = flux.generate_sensfunc(std.boxcar['WAVE'],
                                           std.boxcar['COUNTS'],
                                           std.boxcar['COUNTS_IVAR'],
                                           std_header['AIRMASS'],
                                           std_header['EXPTIME'],
                                           spectrograph,star_type=star_type,star_mag=star_mag,
                                           telluric=telluric,ra=ra,dec=dec,BALM_MASK_WID=BALM_MASK_WID,
                                           nresln=nresln,std_file=std_file,debug=debug)
        sens_dict['ech_orderindx'] = iord
        sens_dicts[str(iord)] = sens_dict
    sens_dicts['norder'] = norder

    return sens_dicts

def ech_flux_science(sci_specobjs,sens_dicts,sci_header,spectrograph=None):
    """
    Flux the internal list of sci_specobjs

    Wrapper to flux.apply_sensfunc()

    Returns
    -------

    """
    norder = sens_dicts['norder']
    if spectrograph is None:
        spectrograph = sci_header['INSTRUME']
        msgs.info('You are working on {:s}'.format(spectrograph))
    for iord in range(norder):
        sens_dict = sens_dicts[str(iord)]
        for sci_obj in sci_specobjs:
            if sci_obj.ech_orderindx == iord:
                flux.apply_sensfunc(sci_obj, sens_dict, sci_header['AIRMASS'],
                                      sci_header['EXPTIME'], spectrograph)


def write_science(sci_specobjs, sci_header, outfile):
    """
    Write the flux-calibrated science spectra

    Parameters
    ----------
    outfile : str

    Returns
    -------

    """
    if len(sci_specobjs) == 0:
        msgs.warn("No science spectra to write to disk!")
    #
    if 'VEL-TYPE' in sci_header.keys():
        helio_dict = dict(refframe=sci_header['VEL-TYPE'],
                          vel_correction=sci_header['VEL'])
    else:
        helio_dict = None
    telescope=None
    if 'LON-OBS' in sci_header.keys():
        telescope = TelescopePar(longitude=sci_header['LON-OBS'],
                                 latitude=sci_header['LAT-OBS'],
                                 elevation=sci_header['ALT-OBS'])
    # KLUDGE ME
    if isinstance(sci_specobjs, list):
        specObjs = specobjs.SpecObjs(sci_specobjs)
    elif isinstance(sci_specobjs, specobjs.SpecObjs):
        specObjs = sci_specobjs
    else:
        msgs.error("BAD INPUT")
    save.save_1d_spectra_fits(specObjs, sci_header, outfile,
                              helio_dict=helio_dict,
                              telescope=telescope, overwrite=True)
    # Step
    #self.steps.append(inspect.stack()[0][3])

