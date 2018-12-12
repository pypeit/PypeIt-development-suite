
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from pypeit.core import coadd
from pypeit.core import load
from pypeit import msgs

from linetools.spectra.utils import collate

def load_single_order(fname,objid=None,order=None,extract='OPT',flux=True):
    """
    Loading single order spectrum from a PypeIt 1D specctrum fits file
    :param file:
    :param objid:
    :param norder:
    :param order:
    :param extract:
    :param flux:
    :return:
    """
    if objid is None:
        objid = 0
    if order is None:
        msgs.error('Please specify which order you want to load')

    # read extension name into a list
    primary_header = fits.getheader(fname, 0)
    nspec = primary_header['NSPEC']
    extnames = [primary_header['EXT0001']] * nspec
    for kk in range(nspec):
        extnames[kk] = primary_header['EXT' + '{0:04}'.format(kk + 1)]
    extnameroot = extnames[0]

    # Figure out which extension is the required data
    ordername = '{0:04}'.format(order)
    extname = extnameroot.replace('OBJ0000', objid)
    extname = extname.replace('ORDER0000', 'ORDER' + ordername)
    try:
        exten = extnames.index(extname) + 1
        msgs.info("Loading extension {:s} of spectrum {:s}".format(extname, fname))
    except:
        msgs.error("Spectrum {:s} does not contain {:s} extension".format(fname, extname))

    spectrum = load.load_1dspec(fname, exten=exten, extract=extract, flux=flux)
    # Polish a bit -- Deal with NAN, inf, and *very* large values that will exceed
    #   the floating point precision of float32 for var which is sig**2 (i.e. 1e38)
    bad_flux = np.any([np.isnan(spectrum.flux), np.isinf(spectrum.flux),
                       np.abs(spectrum.flux) > 1e30,
                       spectrum.sig ** 2 > 1e10,
                       ], axis=0)
    if np.sum(bad_flux):
        msgs.warn("There are some bad flux values in this spectrum.  Will zero them out and mask them (not ideal)")
        spectrum.data['flux'][spectrum.select][bad_flux] = 0.
        spectrum.data['sig'][spectrum.select][bad_flux] = 0.

    return spectrum

def ech_load_spec(files,objid=None,norder=None,order=None,extract='OPT',flux=True):
    """
    files: A list of file names
    objid:
    norder:
    extract:
    flux:
    """

    nfiles = len(files)
    if objid is None:
        objid = ['OBJ0000'] * nfiles
    elif len(objid) == 1:
        objid = objid * nfiles
    elif len(objid) != nfiles:
        msgs.error('The length of objid should be either 1 or equal to the number of spectra files.')

    # Load spectra
    spectra_list = []
    for ii, fname in enumerate(files):
        if norder is None:
            ext_final = fits.getheader(fname,-1)
            norder = ext_final['ORDER']+1
            msgs.info('spectrum {:s} has {:d} orders'.format(fname, norder))
        elif norder <=1:
            msgs.error('The number of orders have to be greater than one for echelle. Longslit data?')

        if order is None:
            msgs.info('Loading all orders into a gaint spectra')
            for iord in range(norder):
                spectrum = load_single_order(fname,objid=objid[ii],order=iord,extract=extract,flux=flux)
                # Append
                spectra_list.append(spectrum)
        elif order >= norder:
            msgs.error('order number cannot greater than the total number of orders')
        else:
            spectrum = load_single_order(fname, objid=objid[ii], order=order, extract=extract, flux=flux)
            # Append
            spectra_list.append(spectrum)
    # Join into one XSpectrum1D object
    spectra = collate(spectra_list)
    # Return
    return spectra

# return spectrum from arrays of wave, flux and sigma
def spec_from_array(wave,flux,sig,**kwargs):

    from linetools.spectra.xspectrum1d import XSpectrum1D
    ituple = (wave, flux, sig)
    spectrum = XSpectrum1D.from_tuple(ituple, **kwargs)
    # Polish a bit -- Deal with NAN, inf, and *very* large values that will exceed
    #   the floating point precision of float32 for var which is sig**2 (i.e. 1e38)
    bad_flux = np.any([np.isnan(spectrum.flux), np.isinf(spectrum.flux),
                       np.abs(spectrum.flux) > 1e30,
                       spectrum.sig ** 2 > 1e10,
                       ], axis=0)
    if np.sum(bad_flux):
        msgs.warn("There are some bad flux values in this spectrum.  Will zero them out and mask them (not ideal)")
        spectrum.data['flux'][spectrum.select][bad_flux] = 0.
        spectrum.data['sig'][spectrum.select][bad_flux] = 0.
    return spectrum

def ech_coadd_spectra(spectra, wave_grid_method='velocity', niter=5,
                  wave_grid_min=None, wave_grid_max=None,v_pix=None,
                  scale_method='auto', do_offset=False, sigrej_final=3.,
                  do_var_corr=False, qafile=None, outfile=None,
                  do_cr=True, **kwargs):


    ech_kwargs = {'echelle':True,'wave_grid_min': wave_grid_min, 'wave_grid_max': wave_grid_max, 'v_pix': v_pix}
    kwargs.update(ech_kwargs)
    spec1d = coadd.coadd_spectra(spectra, wave_grid_method=wave_grid_method, niter=niter,
                        scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                        do_var_corr=do_var_corr, qafile=qafile, outfile=outfile,
                        do_cr=do_cr, debug=False,**kwargs)
    return spec1d

def ech_coadd(files,objids=None,norder=None,extract='OPT',flux=True,giantcoadd=False,
              wave_grid_method='velocity', niter=5,wave_grid_min=None, wave_grid_max=None,v_pix=None,
              scale_method='auto', do_offset=False, sigrej_final=3.,do_var_corr=False,
              qafile=None, outfile=None,do_cr=True, **kwargs):

    if giantcoadd:
        msgs.info('Coadding all orders and exposures at once')
        spectra = ech_load_spec(files, objid=objids, norder=norder, order=None, extract=extract, flux=flux)
        wave_grid = np.zeros((2,spectra.nspec))
        for i in range(spectra.nspec):
            wave_grid[0, i] = spectra[i].wvmin.value
            wave_grid[1, i] = spectra[i].wvmax.value
        ech_kwargs = {'echelle': True, 'wave_grid_min': np.min(wave_grid), 'wave_grid_max': np.max(wave_grid),
                      'v_pix': v_pix}
        kwargs.update(ech_kwargs)
        # Coadding
        spec1d = coadd.coadd_spectra(spectra, wave_grid_method=wave_grid_method, niter=niter,
                                          scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                                          do_var_corr=do_var_corr, qafile=qafile, outfile=outfile,
                                          do_cr=do_cr, **kwargs)
    else:
        msgs.info('Coadding individual orders first and then merge order')
        spectra_list = []
        # Keywords for Table
        rsp_kwargs = {}
        rsp_kwargs['wave_tag'] = '{:s}_WAVE'.format(extract)
        rsp_kwargs['flux_tag'] = '{:s}_FLAM'.format(extract)
        rsp_kwargs['sig_tag'] = '{:s}_FLAM_SIG'.format(extract)
        wave_grid = np.zeros((2,norder))
        for iord in range(norder):
            spectra = ech_load_spec(files, objid=objids, norder=norder, order=iord, extract=extract, flux=flux)
            ech_kwargs = {'echelle': False, 'wave_grid_min': spectra.wvmin.value, 'wave_grid_max': spectra.wvmax.value, 'v_pix': v_pix}
            wave_grid[0,iord] = spectra.wvmin.value
            wave_grid[1,iord] = spectra.wvmax.value
            kwargs.update(ech_kwargs)
            # Coadding the individual orders
            if qafile is not None:
                qafile_iord = qafile+'_%s'%str(iord)
            else:
                qafile_iord =  None
            spec1d_iord = coadd.coadd_spectra(spectra, wave_grid_method=wave_grid_method, niter=niter,
                                       scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                                       do_var_corr=do_var_corr, qafile=qafile_iord, outfile=outfile,
                                       do_cr=do_cr, **kwargs)
            spectrum = spec_from_array(spec1d_iord.wavelength, spec1d_iord.flux, spec1d_iord.sig,**rsp_kwargs)
            spectra_list.append(spectrum)
        # Join into one XSpectrum1D object
        spectra_coadd = collate(spectra_list)
        kwargs['echelle'] = True
        kwargs['wave_grid_min'] = np.min(wave_grid)
        kwargs['wave_grid_max'] = np.max(wave_grid)
        # ToDo: Currently I'm not using the first order due to some problem. Need to add it back after fix the problem.
        spec1d = coadd.coadd_spectra(spectra_coadd[1:], wave_grid_method=wave_grid_method, niter=niter,
                                          scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                                          do_var_corr=do_var_corr, qafile=qafile, outfile=outfile,
                                          do_cr=do_cr, **kwargs)

    return spec1d

def test_nires(giantcoadd=False):
    #scifiles = ['/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T100254.698_FLUX.fits',
    #            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T100949.328_FLUX.fits',
    #            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T101642.428_FLUX.fits',
    #            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T102337.058_FLUX.fits']
    #objids = ['OBJ0001','OBJ0002','OBJ0002','OBJ0001']
    norder =5
    datapath = '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/'
    cat = np.genfromtxt(datapath+'J0252_objinfo.txt',dtype=str)
    filenames = cat[:,0]
    scifiles = []
    for i in range(len(filenames)):
        filename = datapath+filenames[i]
        scifiles += [filename.replace('.fits','_FLUX.fits')]
    objids = cat[:,1]

    # Coadding
    kwargs={}
    spec1d = ech_coadd(scifiles, objids=objids, norder=norder, extract='OPT', flux=True,giantcoadd=giantcoadd,
              wave_grid_method='velocity', niter=5,wave_grid_min=None, wave_grid_max=None, v_pix=None,
              scale_method='median', do_offset=False, sigrej_final=3.,
              do_var_corr=False, qafile='test_nires', outfile=None, do_cr=True,**kwargs)

def ech_load_xidl(files,extn=4,norder=6,order=None,extract='OPT',sensfile=None,AB=True):
    """
    ONLY tested for GNIRS for now
    files: A list of file names
    objid:
    norder:
    extract:
    flux:
    """

    nfiles = len(files)

    # Keywords for Table
    rsp_kwargs = {}
    rsp_kwargs['wave_tag'] = '{:s}_WAVE'.format(extract)
    rsp_kwargs['flux_tag'] = '{:s}_FLAM'.format(extract)
    rsp_kwargs['sig_tag'] = '{:s}_FLAM_SIG'.format(extract)

    # Reading magsens func
    if sensfile is not None:
        magfunc  = fits.getdata(sensfile, 0)
        loglams = fits.getdata(sensfile, 1)
        exptime = np.zeros(nfiles)
        for ii, fname in enumerate(files):
            exptime[ii] = fits.getheader(fname,0)['EXPTIME']

    # Load spectra
    spectra_list = []
    for ii, fname in enumerate(files):
        ext_final = fits.getdata(fname, -1)
        head_final = fits.getheader(fname,-1)
        if norder is None:
            if AB:
                norder = head_final['NAXIS2']/2
            else:
                norder = head_final['NAXIS2']
            msgs.info('spectrum {:s} has {:d} orders'.format(fname, norder))
        elif norder <=1:
            msgs.error('The number of orders have to be greater than one for echelle. Longslit data?')

        if AB:
            if order is None:
                msgs.info('Loading all orders into a gaint spectra')
                for iord in range(norder):
                    for aa in range(2):
                        wave = ext_final['WAVE_{:s}'.format(extract)][2 * iord + aa]
                        flux = ext_final['FLUX_{:s}'.format(extract)][2 * iord + aa]
                        sig = ext_final['SIG_{:s}'.format(extract)][2 * iord + aa]
                        if sensfile is not None:
                            magfunc1 = np.interp(np.log10(wave), loglams, magfunc[iord, :])
                            sensfunc = 10.0 ** (0.4 * magfunc1)
                            scale = sensfunc / exptime[ii]
                            flux = flux * scale
                            sig = sig * scale
                        spectrum = spec_from_array(wave, flux, sig, **rsp_kwargs)
                        # Append
                        spectra_list.append(spectrum)
            elif order >= norder:
                msgs.error('order number cannot greater than the total number of orders')
            else:
                for aa in range(2):
                    wave = ext_final['WAVE_{:s}'.format(extract)][2 * order + aa]
                    flux = ext_final['FLUX_{:s}'.format(extract)][2 * order + aa]
                    sig = ext_final['SIG_{:s}'.format(extract)][2 * order + aa]
                    if sensfile is not None:
                        magfunc1 = np.interp(np.log10(wave), loglams, magfunc[order, :])
                        sensfunc = 10.0 ** (0.4 * magfunc1)
                        scale = sensfunc / exptime[ii]
                        flux = flux * scale
                        sig = sig * scale
                    spectrum = spec_from_array(wave, flux, sig, **rsp_kwargs)
                    # Append
                    spectra_list.append(spectrum)
        else:
            if order is None:
                msgs.info('Loading all orders into a gaint spectra')
                for iord in range(norder):
                    wave = ext_final['WAVE_{:s}'.format(extract)][iord]
                    flux = ext_final['FLUX_{:s}'.format(extract)][iord]
                    sig = ext_final['SIG_{:s}'.format(extract)][iord]
                    if sensfile is not None:
                        magfunc1 = np.interp(np.log10(wave), loglams, magfunc[iord, :])
                        sensfunc = 10.0 ** (0.4 * magfunc1)
                        scale = sensfunc / exptime[ii]
                        flux = flux * scale
                        sig = sig * scale
                    spectrum = spec_from_array(wave, flux, sig, **rsp_kwargs)
                    # Append
                    spectra_list.append(spectrum)
            elif order >= norder:
                msgs.error('order number cannot greater than the total number of orders')
            else:
                wave = ext_final['WAVE_{:s}'.format(extract)][order]
                flux = ext_final['FLUX_{:s}'.format(extract)][order]
                sig = ext_final['SIG_{:s}'.format(extract)][order]
                if sensfile is not None:
                    magfunc1 = np.interp(np.log10(wave), loglams, magfunc[order, :])
                    sensfunc = 10.0 ** (0.4 * magfunc1)
                    scale = sensfunc / exptime[ii]
                    flux = flux * scale
                    sig = sig * scale
                spectrum = spec_from_array(wave, flux, sig, **rsp_kwargs)
                # Append
                spectra_list.append(spectrum)

    # Join into one XSpectrum1D object
    spectra = collate(spectra_list)
    # Return
    return spectra


def test_gnirs(giantcoadd=False,qafile='testgnirs'):
    norder=6
    order = 0

    scifiles = ['/Users/feige/Work/Observations/GN-2015A-Q-28_P338+29/Redux/Science/PSO338+29_0/sci-cN20150707S0189-192.fits',
                '/Users/feige/Work/Observations/GN-2015A-Q-28_P338+29/Redux/Science/PSO338+29_0/sci-cN20150707S0193-196.fits',
                '/Users/feige/Work/Observations/GN-2015A-Q-28_P338+29/Redux/Science/PSO338+29_0/sci-cN20150707S0220-223.fits',
                '/Users/feige/Work/Observations/GN-2015A-Q-28_P338+29/Redux/Science/PSO338+29_0/sci-cN20150707S0224-223.fits']
    sensfile = '/Users/feige/Work/Observations/GN-2015A-Q-28_P338+29/Redux/Combine/HIP111538_0_sens.fits'
    if giantcoadd:
        spectra = ech_load_xidl(scifiles,extn=4,norder=norder,order=None,extract='OPT',sensfile=sensfile,AB=True)
        kwargs={}
        ech_kwargs = {'echelle': True, 'wave_grid_min': None, 'wave_grid_max': None}
        kwargs.update(ech_kwargs)
        # Coadding
        spec1d = coadd.coadd_spectra(spectra, wave_grid_method='velocity', niter=5,
                  scale_method='auto', do_offset=False, sigrej_final=3.,
                  do_var_corr=False, qafile=qafile, outfile=None, do_cr=True,**kwargs)
    else:
        msgs.info('Coadding individual orders first and then merge order')
        spectra_list = []
        # Keywords for Table
        extract='OPT'
        rsp_kwargs = {}
        rsp_kwargs['wave_tag'] = '{:s}_WAVE'.format(extract)
        rsp_kwargs['flux_tag'] = '{:s}_FLAM'.format(extract)
        rsp_kwargs['sig_tag'] = '{:s}_FLAM_SIG'.format(extract)
        wave_grid = np.zeros((2,norder))
        for iord in range(norder):
            spectra = ech_load_xidl(scifiles, extn=4, norder=norder, order=iord, extract='OPT', sensfile=sensfile,
                                    AB=True)
            kwargs = {}
            ech_kwargs = {'echelle': False, 'wave_grid_min': spectra.wvmin.value, 'wave_grid_max': spectra.wvmax.value, 'v_pix': None}
            wave_grid[0,iord] = spectra.wvmin.value
            wave_grid[1,iord] = spectra.wvmax.value
            kwargs.update(ech_kwargs)
            # Coadding the individual orders
            if qafile is not None:
                qafile_iord = qafile+'_%s'%str(iord)
            else:
                qafile_iord =  None
            spec1d_iord = coadd.coadd_spectra(spectra, wave_grid_method='velocity', niter=5,
                  scale_method='auto', do_offset=False, sigrej_final=3.,
                  do_var_corr=False, qafile=qafile_iord, outfile=None, do_cr=True,**kwargs)
            spectrum = spec_from_array(spec1d_iord.wavelength, spec1d_iord.flux, spec1d_iord.sig,**rsp_kwargs)
            spectra_list.append(spectrum)
        # Join into one XSpectrum1D object
        spectra_coadd = collate(spectra_list)
        kwargs['echelle'] = True
        kwargs['wave_grid_min'] = np.min(wave_grid)
        kwargs['wave_grid_max'] = np.max(wave_grid)
        spec1d = coadd.coadd_spectra(spectra_coadd, wave_grid_method='velocity', niter=5,
                  scale_method='auto', do_offset=False, sigrej_final=3.,
                  do_var_corr=False, qafile=qafile, outfile=None, do_cr=True,**kwargs)


#test_nires(giantcoadd=False)
test_gnirs(giantcoadd=False)
