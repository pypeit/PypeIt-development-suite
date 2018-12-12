from ech_fluxspec import *
from ech_coadd import *

## test ech_fluxspec

def flux_example(debug=True):
    """
    test single NIRES frame
    :param debug:
    :return:
    """
    stdframe = '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_HIP13917_V8p6_NIRES_2018Oct01T094225.598.fits'
    sciframe = '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T100254.698.fits'
    sens_dicts = ech_generate_sensfunc(stdframe,telluric=True, star_type='A0',
                          star_mag=8.6, ra=None, dec=None, std_file = None, BALM_MASK_WID=15., nresln=None,debug=debug)
    ech_save_master(sens_dicts, outfile='MasterSensFunc_NIRES.fits')
    sens_dicts = ech_load_master('MasterSensFunc_NIRES.fits')
    sci_specobjs, sci_header = ech_load_specobj(sciframe)
    ech_flux_science(sci_specobjs,sens_dicts,sci_header,spectrograph=None)
    write_science(sci_specobjs, sci_header, sciframe[:-5]+'_FLUX.fits')


def flux_example2(debug=False):
    """
    test NIRES with a list of files
    :return:
    """
    datapath = '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/'
    stdframe = datapath+'spec1d_HIP13917_V8p6_NIRES_2018Oct01T094225.598.fits'
    cat = np.genfromtxt(datapath+'J0252_objinfo.txt',dtype=str)
    filenames = cat[:,0]

    sens_dicts = ech_generate_sensfunc(stdframe,telluric=True, star_type='A0',
                                       star_mag=8.6, ra=None, dec=None, std_file=None, BALM_MASK_WID=10., nresln=None,
                                       debug=debug)
    ech_save_master(sens_dicts, outfile='MasterSensFunc_NIRES.fits')
    for i in range(len(filenames)):
        sciframe = datapath+filenames[i]
        sci_specobjs, sci_header = ech_load_specobj(sciframe)
        ech_flux_science(sci_specobjs,sens_dicts,sci_header,spectrograph=None)
        write_science(sci_specobjs, sci_header, sciframe[:-5]+'_FLUX.fits')


def ech_coadd_spectra(spectra, wave_grid_method='velocity', niter=5,
                  wave_grid_min=None, wave_grid_max=None,v_pix=None,
                  scale_method='auto', do_offset=False, sigrej_final=3.,
                  do_var_corr=False, qafile=None, outfile=None,
                  do_cr=True, **kwargs):
    """
    Deprecated
    """
    ech_kwargs = {'echelle':True,'wave_grid_min': wave_grid_min, 'wave_grid_max': wave_grid_max, 'v_pix': v_pix}
    kwargs.update(ech_kwargs)
    spec1d = coadd.coadd_spectra(spectra, wave_grid_method=wave_grid_method, niter=niter,
                        scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                        do_var_corr=do_var_corr, qafile=qafile, outfile=outfile,
                        do_cr=do_cr, debug=False,**kwargs)
    return spec1d


def coadd_nires(giantcoadd=False,debug=False):
    #scifiles = ['/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T100254.698_FLUX.fits',
    #            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T100949.328_FLUX.fits',
    #            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T101642.428_FLUX.fits',
    #            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T102337.058_FLUX.fits']
    #objids = ['OBJ0001','OBJ0002','OBJ0002','OBJ0001']
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
    spec1d = ech_coadd(scifiles, objids=objids,extract='OPT', flux=True,giantcoadd=giantcoadd,
              wave_grid_method='velocity', niter=5,wave_grid_min=None, wave_grid_max=None, v_pix=None,
              scale_method='median', do_offset=False, sigrej_final=3.,
              do_var_corr=False, qafile='test_nires', outfile=None, do_cr=True,debug=debug,**kwargs)
    return spec1d


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


def coadd_gnirs(giantcoadd=False,qafile='testgnirs',debug=False):
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
                  do_var_corr=False, qafile=qafile, outfile=None, do_cr=True,debug=debug,**kwargs)
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
                  do_var_corr=False, qafile=qafile_iord, outfile=None, do_cr=True,debug=debug,**kwargs)
            spectrum = spec_from_array(spec1d_iord.wavelength, spec1d_iord.flux, spec1d_iord.sig,**rsp_kwargs)
            spectra_list.append(spectrum)
        # Join into one XSpectrum1D object
        spectra_coadd = collate(spectra_list)
        kwargs['echelle'] = True
        kwargs['wave_grid_min'] = np.min(wave_grid)
        kwargs['wave_grid_max'] = np.max(wave_grid)
        spec1d = coadd.coadd_spectra(spectra_coadd, wave_grid_method='velocity', niter=5,
                  scale_method='auto', do_offset=False, sigrej_final=3.,
                  do_var_corr=False, qafile=qafile, outfile=None, do_cr=True,debug=debug,**kwargs)
    return spectra_coadd,spec1d

def try_mergeorder():
    spectra_coadd, spec1d = coadd_gnirs(giantcoadd=False)

    norder = spectra_coadd.nspec

    ## Sort spectra
    wvmin = np.zeros(norder)
    for iord in range(norder):
        wvmin[iord] = spectra_coadd[iord].wvmin.value
    indsort = np.argsort(wvmin)
    spectra_coadd_sort = spectra_coadd[indsort]

    for iord in range(norder):
        wave = spectra_coadd_sort[iord].data['wave'].filled(0.)[0]
        flux = spectra_coadd_sort[iord].data['flux'].filled(0.)[0]
        plt.plot(wave, flux)
    plt.show()
    for ii in range(norder - 1):
        wvmin1, wvmax1 = spectra_coadd_sort[ii].wvmin, spectra_coadd_sort[ii].wvmax
        wvmin2, wvmax2 = spectra_coadd_sort[ii + 1].wvmin, spectra_coadd_sort[ii + 1].wvmax

        wave1 = spectra_coadd_sort[iord].data['wave'].filled(0.)[0]
        flux1 = spectra_coadd_sort[iord].data['flux'].filled(0.)[0]



#flux_example(debug=True)
#flux_example2()

coadd_nires(giantcoadd=False,debug=True)
#spectra_coadd,spec1d = coadd_gnirs(giantcoadd=False,debug=False)

