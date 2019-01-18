from ech_fluxspec import *
#from ech_coadd import *
from pypeit.ech_coadd import *
from pypeit import fluxspec
import matplotlib.pyplot as plt
## test ech_fluxspec


def flux_single(debug=True,datapath='./',stdframe=None,sciframe=None,star_type='A0',star_mag=8.0,
                 ra=None, dec=None, std_file=None, BALM_MASK_WID=15., nresln=None,outfile=None):
    """
    test single frame
    :param debug:
    :return:
    """
    sens_dicts = ech_generate_sensfunc(datapath+stdframe,telluric=True, star_type=star_type,
                                       star_mag=star_mag,ra=ra, dec=dec, std_file = std_file,
                                       BALM_MASK_WID=BALM_MASK_WID, nresln=nresln,debug=debug)
    ech_save_master(sens_dicts, outfile=datapath+outfile)
    #sens_dicts = ech_load_master(datapath+'MasterSensFunc_NIRES.fits')
    sci_specobjs, sci_header = ech_load_specobj(datapath+sciframe)
    ech_flux_science(sci_specobjs,sens_dicts,sci_header,spectrograph=None)
    write_science(sci_specobjs, sci_header, datapath+sciframe[:-5]+'_FLUX.fits')
    return sens_dicts


def flux_example(debug=True,datapath='./'):
    """
    test single NIRES frame
    :param debug:
    :return:
    """
    stdframe = datapath+'/spec1d_HIP13917_V8p6_NIRES_2018Oct01T094225.598.fits'
    sciframe = datapath+'/spec1d_J0252-0503_NIRES_2018Oct01T100254.698.fits'
    sens_dicts = ech_generate_sensfunc(stdframe,telluric=True, star_type='A0',
                          star_mag=8.6, ra=None, dec=None, std_file = None, BALM_MASK_WID=15., nresln=None,debug=debug)
    ech_save_master(sens_dicts, outfile=datapath+'MasterSensFunc_NIRES.fits')
    sens_dicts = ech_load_master(datapath+'MasterSensFunc_NIRES.fits')
    sci_specobjs, sci_header = ech_load_specobj(sciframe)
    ech_flux_science(sci_specobjs,sens_dicts,sci_header,spectrograph=None)
    write_science(sci_specobjs, sci_header, sciframe[:-5]+'_FLUX.fits')


def flux_example2(debug=False,datapath='./'):
    """
    test NIRES with a list of files
    :return:
    """
    stdframe = datapath+'spec1d_HIP13917_V8p6_NIRES_2018Oct01T094225.598.fits'
    cat = np.genfromtxt(datapath+'J0252_objinfo.txt',dtype=str)
    filenames = cat[:,0]

    sens_dicts = ech_generate_sensfunc(stdframe,telluric=True, star_type='A0',
                                       star_mag=8.6, ra=None, dec=None, std_file=None, BALM_MASK_WID=15., nresln=None,
                                       debug=debug)
    ech_save_master(sens_dicts, outfile='MasterSensFunc_NIRES.fits')
    for i in range(len(filenames)):
        sciframe = datapath+filenames[i]
        sci_specobjs, sci_header = ech_load_specobj(sciframe)
        ech_flux_science(sci_specobjs,sens_dicts,sci_header,spectrograph=None)
        write_science(sci_specobjs, sci_header, sciframe[:-5]+'_flux.fits')

def flux_nires(debug=False,datapath='./',objinfo='J1135_info.txt',stdframe='spec1d_HIP53735_NIRES_2018Jun04T055220.216.fits'):
    """
    test NIRES with a list of files
    :return:
    """
    stdframe = datapath+stdframe
    cat = np.genfromtxt(datapath+objinfo,dtype=str)
    filenames = cat[:,0]

    sens_dicts = ech_generate_sensfunc(stdframe,telluric=True, star_type='A0',
                                       star_mag=8.89, ra=None, dec=None, std_file=None, BALM_MASK_WID=15., nresln=None,
                                       debug=debug)
    ech_save_master(sens_dicts, outfile='MasterSensFunc_NIRES.fits')
    for i in range(len(filenames)):
        sciframe = datapath+filenames[i]
        sci_specobjs, sci_header = ech_load_specobj(sciframe)
        ech_flux_science(sci_specobjs,sens_dicts,sci_header,spectrograph=None)
        write_science(sci_specobjs, sci_header, sciframe[:-5]+'_FLUX.fits')

def ech_flux_new(spectragraph='keck_nires',debug=False,datapath='./',star_type='A0',star_mag=8.6,
                   objinfo='J1135_info.txt',stdframe='spec1d_HIP53735_NIRES_2018Jun04T055220.216.fits'):
    """
    test NIRES with a list of files
    :return:
    """
    cat = np.genfromtxt(datapath+objinfo,dtype=str)
    filenames = cat[:,0]

    for i in range(len(filenames)):
        sciframe = datapath+filenames[i]

        FxSpec = fluxspec.EchFluxSpec(std_spec1d_file=datapath+stdframe,
                                      sci_spec1d_file=sciframe,
                                      spectrograph=spectragraph,
                                      telluric=True,
                                      sens_file=datapath+'sens_'+stdframe,
                                      star_type=star_type,
                                      star_mag=star_mag,
                                      debug=debug)
        if i==0:
            _ = FxSpec.generate_sensfunc()
            _ = FxSpec.save_master(FxSpec.sens_dict, outfile=datapath+'sens_'+stdframe)
        FxSpec.flux_science()
        FxSpec.write_science(sciframe[:-5]+'_flux.fits')


def ech_coadd_new(giantcoadd=False,debug=False,datapath='./',objinfo='J0252_objinfo.txt',qafile='ech_coadd',
                  outfile='J0910_GNIRS.fits'):

    cat = np.genfromtxt(datapath+objinfo,dtype=str)
    filenames = cat[:,0]
    scifiles = []
    for i in range(len(filenames)):
        filename = datapath+filenames[i]
        scifiles += [filename.replace('.fits','_flux.fits')]
    objids = cat[:,1]

    # Coadding
    kwargs={}
    spec1d = ech_coadd(scifiles, objids=objids,extract='OPT', flux=True,giantcoadd=giantcoadd,
              wave_grid_method='velocity', niter=5,wave_grid_min=None, wave_grid_max=None, v_pix=None,
              scale_method='median', do_offset=False, sigrej_final=3.,
              do_var_corr=False, qafile=datapath+qafile, outfile=datapath+outfile, do_cr=True,debug=debug,**kwargs)
    return spec1d

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


def coadd_nires(giantcoadd=False,debug=False,datapath='./',objinfo='J0252_objinfo.txt'):
    #scifiles = ['/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T100254.698_FLUX.fits',
    #            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T100949.328_FLUX.fits',
    #            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T101642.428_FLUX.fits',
    #            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T102337.058_FLUX.fits']
    #objids = ['OBJ0001','OBJ0002','OBJ0002','OBJ0001']
    cat = np.genfromtxt(datapath+objinfo,dtype=str)
    filenames = cat[:,0]
    scifiles = []
    for i in range(len(filenames)):
        filename = datapath+filenames[i]
        scifiles += [filename.replace('.fits','_flux.fits')]
    objids = cat[:,1]

    # Coadding
    kwargs={}
    spec1d = ech_coadd(scifiles, objids=objids,extract='OPT', flux=True,giantcoadd=giantcoadd,
              wave_grid_method='velocity', niter=5,wave_grid_min=None, wave_grid_max=None, v_pix=None,
              scale_method='median', do_offset=False, sigrej_final=3.,
              do_var_corr=False, qafile='ech_coad', outfile=None, do_cr=True,debug=debug,**kwargs)
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


def coadd_gnirs(giantcoadd=False,qafile='testgnirs',debug=False,datapath='./'):
    norder=6
    order = 0

    scifiles = [datapath+'/sci-cN20150707S0189-192.fits',
                datapath+'/sci-cN20150707S0193-196.fits',
                datapath+'/sci-cN20150707S0220-223.fits',
                datapath+'/sci-cN20150707S0224-223.fits']
    sensfile = datapath+'/HIP111538_0_sens.fits'
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
    return spec1d,spectra_coadd,kwargs


def readgnirsxidl(filename):
    fitsfile = fits.open(filename)
    header = fitsfile[0].header
    wave, flux, error = fitsfile[2].data, fitsfile[0].data, fitsfile[1].data

    return wave, flux, error, header


def try_mergeorder(spectra_coadd,wave_grid_method='velocity',kwargs=None):

    norder = spectra_coadd.nspec

    ## Sort spectra
    wvmin = np.zeros(norder)
    for iord in range(norder):
        wvmin[iord] = spectra_coadd[iord].wvmin.value
    indsort = np.argsort(wvmin)
    spectra_coadd_sort = spectra_coadd[indsort]

    # Final wavelength array
    new_wave = coadd.new_wave_grid(spectra_coadd_sort.data['wave'], wave_method=wave_grid_method, **kwargs)

    # Rebin
    rspec = spectra_coadd_sort.rebin(new_wave*units.AA, all=True, do_sig=True, grow_bad_sig=True,
                          masking='none')

    # Define mask -- THIS IS THE ONLY ONE TO USE
    rmask = rspec.data['sig'].filled(0.) <= 0.

    # S/N**2, weights
    sn2, weights = coadd.sn_weight(rspec, rmask)

    for iord in range(norder):
        wave = spectra_coadd_sort[iord].data['wave'].filled(0.)[0]
        flux = spectra_coadd_sort[iord].data['flux'].filled(0.)[0]
        sig = spectra_coadd_sort[iord].data['sig'].filled(0.)[0]
        plt.plot(wave, flux/sig)
    plt.show()
    for ii in range(norder - 1):
        wvmin1, wvmax1 = spectra_coadd_sort[ii].wvmin, spectra_coadd_sort[ii].wvmax
        wvmin2, wvmax2 = spectra_coadd_sort[ii + 1].wvmin, spectra_coadd_sort[ii + 1].wvmax

        wave1 = spectra_coadd_sort[iord].data['wave'].filled(0.)[0]
        flux1 = spectra_coadd_sort[iord].data['flux'].filled(0.)[0]


def other_tests():
    # flux XSHOOTER
    datapath = '/Users/feige/Dropbox/PypeIt_DATA/XSHOOTER/J0439/NIR/Science/'
    stdframe = 'spec1d_STD,TELLURIC_XShooter_NIR_2018Oct08T232940.178.fits'
    sciframe = 'spec1d_STD,TELLURIC_XShooter_NIR_2018Oct08T233037.583.fits' # test fluxing with telluric star

    scilist = ['spec1d_J0439_XShooter_NIR_2018Oct09T063213.112.fits','spec1d_J0439_XShooter_NIR_2018Oct09T064021.266.fits',
               'spec1d_J0439_XShooter_NIR_2018Oct09T064829.420.fits','spec1d_J0439_XShooter_NIR_2018Oct09T065654.883.fits',
               'spec1d_J0439_XShooter_NIR_2018Oct09T070503.037.fits','spec1d_J0439_XShooter_NIR_2018Oct09T071311.193.fits']
    objids = ['OBJ0000','OBJ0000','OBJ0000','OBJ0001','OBJ0001','OBJ0001']
    star_mag = 5.644
    star_type = 'B9'
    outfile = 'HIP095619_XSHOOTER_NIR_SENS'

    flux_single(debug=False,datapath=datapath,stdframe=stdframe,sciframe=sciframe,
                star_type=star_type,star_mag=star_mag,BALM_MASK_WID=15., nresln=None,outfile=outfile)
    sci_specobjs, sci_header = ech_load_specobj(datapath + sciframe[:-5] + '_FLUX.fits')
    wavemask = sci_specobjs[13].optimal['WAVE']>1000.0*units.AA
    plt.plot(sci_specobjs[13].optimal['WAVE'][wavemask],sci_specobjs[13].optimal['FLAM'][wavemask])
    plt.show()

    sens_dicts = ech_generate_sensfunc(datapath + stdframe, telluric=True, star_type=star_type,
                                       star_mag=star_mag, ra=None, dec=None, std_file=None,
                                       BALM_MASK_WID=5.0, nresln=None, debug=True)
    ech_save_master(sens_dicts, outfile=datapath + outfile)


    scifiles = []
    for i in range(len(scilist)):
        sciframe = scilist[i]
        sci_specobjs, sci_header = ech_load_specobj(datapath + sciframe)
        ech_flux_science(sci_specobjs, sens_dicts, sci_header, spectrograph=None)
        write_science(sci_specobjs, sci_header, datapath + sciframe[:-5] + '_FLUX.fits')

        filename = datapath+sciframe
        scifiles += [filename.replace('.fits','_FLUX.fits')]

    # Coadding
    kwargs={}
    spec1d = ech_coadd(scifiles, objids=objids,extract='OPT', flux=True,giantcoadd=False,
              wave_grid_method='velocity', niter=5,wave_grid_min=None, wave_grid_max=None, v_pix=None,
              scale_method='median', do_offset=False, sigrej_final=3.,
              do_var_corr=False, qafile='test_xshooter.png', outfile=None, do_cr=True,debug=True,**kwargs)


    #Test GNIRS
    spec1d,spectra_coadd,kwargs = coadd_gnirs(giantcoadd=False,debug=True,datapath='/Users/feige/Dropbox/PypeIt_Redux/GNIRS/')
    wave, flux, error, header = readgnirsxidl('/Users/feige/Dropbox/PypeIt_Redux/GNIRS/PSO338+29_forpypeit.fits')
    cat = np.genfromtxt('/Users/feige/Dropbox/PypeIt_Redux/GNIRS/mods_spectrum_p338.txt')
    wave_opt, flux_opt, error_opt = cat[:, 0], cat[:, 1]/1e-17, cat[:, 2]/1e-17

    plt.figure(figsize=(12,4))
    plt.plot(wave_opt, flux_opt,'k-',label='MODS')
    plt.plot(wave,flux,'b-',label='XIDL')
    plt.plot(wave,error,'b-',alpha=0.5)
    plt.plot(spec1d.wavelength,spec1d.flux,'r-',label='PypeIt')
    plt.plot(spec1d.wavelength,spec1d.sig,'r-',alpha=0.5)
    plt.xlim([9000.,24800.])
    plt.ylim([-0.5,2.0])
    plt.legend(loc=1,fontsize=16)
    plt.show()

#flux_example(debug=True,datapath='/Users/feige/Dropbox/PypeIt_Redux/NIRES/')
#flux_example2(debug=False,datapath='/Users/feige/Dropbox/PypeIt_Redux/NIRES/')
#coadd_nires(giantcoadd=False,debug=True,datapath='/Users/feige/Dropbox/PypeIt_Redux/NIRES/')

#flux_example(debug=True,datapath='/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930_20190103/Science/')
#flux_example2(debug=False,datapath='/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930_20190103/Science/')
#coadd_nires(giantcoadd=False,debug=True,datapath='/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930_20190103/Science/')

### 2018 June data reduction
#flux_nires(debug=False,datapath='/Users/feige/Work/Observations/201805_NIRES/pypeit_redux/20180604/Flux/',\
#           objinfo='J1135_info.txt',stdframe='spec1d_HIP53735_NIRES_2018Jun04T055220.216.fits')
#coadd_nires(giantcoadd=False,debug=True,datapath='/Users/feige/Work/Observations/201805_NIRES/pypeit_redux/20180604/Flux/',\
#            objinfo='J1135_info.txt')
#flux_nires(debug=False,datapath='/Users/feige/Work/Observations/201805_NIRES/pypeit_redux/20180604/Flux/',\
#           objinfo='J2125_info.txt',stdframe='spec1d_HIP107535_NIRES_2018Jun04T150456.336.fits')
#coadd_nires(giantcoadd=False,debug=True,datapath='/Users/feige/Work/Observations/201805_NIRES/pypeit_redux/20180604/Flux/',\
#            objinfo='J2125_info.txt')
#flux_nires(debug=False,datapath='/Users/feige/Work/Observations/201805_NIRES/pypeit_redux/20180604/Flux/',\
#           objinfo='J1316_info.txt',stdframe='spec1d_HIP61471_NIRES_2018Jun04T084258.466.fits')
#coadd_nires(giantcoadd=False,debug=True,datapath='/Users/feige/Work/Observations/201805_NIRES/pypeit_redux/20180604/Flux/',\
#           objinfo='J1316_info.txt')
#flux_nires(debug=False,datapath='/Users/feige/Work/Observations/201805_NIRES/pypeit_redux/20180604/Flux/',\
#           objinfo='J1724_info.txt',stdframe='spec1d_HIP87643_NIRES_2018Jun04T141653.306.fits')
#coadd_nires(giantcoadd=False,debug=True,datapath='/Users/feige/Work/Observations/201805_NIRES/pypeit_redux/20180604/Flux/',\
#           objinfo='J1724_info.txt')


### 2018 Oct NIRES data reduction
ech_flux_new(debug=True,datapath='/Users/feige/Work/Observations/NIRES/pypeit_redux/ut181001/Science/',\
            objinfo='J0252_info.txt',stdframe='spec1d_HIP13917_V8p6_NIRES_2018Oct01T094225.598.fits',
            star_type='A0',star_mag=8.6)
#coadd_nires(giantcoadd=False,debug=True,datapath='/Users/feige/Work/Observations/NIRES/pypeit_redux/ut181001/Science/',\
#           objinfo='J0252_info.txt')

### 2018 May GNIRS reduction
#ech_flux_new(debug=False,datapath='/Users/feige/Dropbox/PypeIt_Redux/GNIRS/ut180517/Science/',\
#            objinfo='J0910_info.txt',stdframe='spec1d_HIP43018_GNIRS_2018May16T225142.936.fits',
#            star_type='A0',star_mag=8.72)
#ech_coadd_new(giantcoadd=False,debug=True,datapath='/Users/feige/Dropbox/PypeIt_Redux/GNIRS/ut180517/Science/',
#              objinfo='J0910_info.txt',qafile='ech_coadd',outfile='J0910_GNIRS.fits')

#ech_flux_new(debug=False,datapath='/Users/feige/Dropbox/PypeIt_Redux/GNIRS/P006p39/Flux/',\
#            objinfo='P006_info.txt',stdframe='spec1d_A2VStar_GNIRS_2016Aug05T002219.142.fits',
#            star_type='A2',star_mag=6.688)
#ech_coadd_new(giantcoadd=False,debug=True,datapath='/Users/feige/Dropbox/PypeIt_Redux/GNIRS/P006p39/Flux/',
#              objinfo='P006_info.txt',qafile='ech_coadd',outfile='P006_GNIRS.fits')
