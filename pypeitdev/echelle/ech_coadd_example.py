from ech_fluxspec import *
from pypeit.ech_coadd import *
from pypeit import fluxspec
import matplotlib.pyplot as plt

"""
The following two functions just small wrappers for the flux and ech_coadd. 
There are three examples below:
    Pisco GNIRS data
    NIRES z~7 quasar data
    XSHOOTER lensed quasar data
"""
def ech_flux_new(spectragraph='keck_nires',debug=False,datapath='./',star_type='A0',star_mag=8.6,BALM_MASK_WID=20.0,
                 norder=5,resolution=3000,polycorrect=True,polysens=False, objinfo='J1135_info.txt',telluric=True,
                 stdframe='spec1d_HIP53735_NIRES_2018Jun04T055220.216.fits'):
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
                                      telluric=telluric,
                                      sens_file=datapath+'sens_'+stdframe,
                                      star_type=star_type,
                                      star_mag=star_mag,
                                      BALM_MASK_WID = BALM_MASK_WID,
                                      resolution = resolution,
                                      polycorrect = polycorrect,
                                      polysens = polysens,
                                      norder = norder,
                                      debug=debug)
        if i==0:
            _ = FxSpec.generate_sensfunc()
            _ = FxSpec.save_master(FxSpec.sens_dict, outfile=datapath+'sens_'+stdframe)
        FxSpec.flux_science()
        FxSpec.write_science(sciframe[:-5]+'_flux.fits')


def ech_coadd_new(giantcoadd=False,debug=False,datapath='./',objinfo='J0252_objinfo.txt',qafile='ech_coadd',
                  outfile='J0910_GNIRS.fits',flux=True):

    cat = np.genfromtxt(datapath+objinfo,dtype=str)
    filenames = cat[:,0]
    scifiles = []
    for i in range(len(filenames)):
        filename = datapath+filenames[i]
        if flux:
            scifiles += [filename.replace('.fits','_flux.fits')]
        else:
            scifiles +=[filename]
    objids = cat[:,1]

    # Coadding
    kwargs={}
    spec1d = ech_coadd(scifiles, objids=objids,extract='OPT', flux=flux,giantcoadd=giantcoadd,
              wave_grid_method='velocity', niter=5,wave_grid_min=None, wave_grid_max=None, v_pix=None,
              scale_method='median', do_offset=False, sigrej_final=3.,
              do_var_corr=False, qafile=datapath+qafile, outfile=datapath+outfile, do_cr=True,debug=debug,**kwargs)
    return spec1d


############################### Begin GNIRS reduction
### Flux and Coadd Pisco, observed on ut170331 and ut170403
#ech_flux_new(debug=True,datapath='/Users/feige/Dropbox/PypeIt_Redux/GNIRS/ut170331/Science/',\
#            objinfo='pisco_info.txt',stdframe='spec1d_HIP68868_GNIRS_2017Mar31T112201.531.fits',
#            star_type='A0',star_mag=8.63,norder=5,resolution=1000,BALM_MASK_WID=100.0,polycorrect=True,polysens=False)
#ech_coadd_new(giantcoadd=False,debug=True,datapath='/Users/feige/Dropbox/PypeIt_Redux/GNIRS/ut170331/Science/',
#              objinfo='pisco_info.txt',qafile='Pisco_GNIRS_ut170331_20190121',outfile='Pisco_GNIRS_ut170331_20190121.fits')
#ech_flux_new(debug=False,datapath='/Users/feige/Dropbox/PypeIt_Redux/GNIRS/ut170403/Science/',\
#            objinfo='pisco_info.txt',stdframe='spec1d_HIP68868_GNIRS_2017Apr03T120855.582.fits',
#            star_type='A0',star_mag=8.63,norder=5,resolution=1000,BALM_MASK_WID=100.0,polycorrect=True,polysens=False)
#ech_coadd_new(giantcoadd=False,debug=False,datapath='/Users/feige/Dropbox/PypeIt_Redux/GNIRS/ut170403/Science/',
#              objinfo='pisco_info.txt',qafile='Pisco_GNIRS_ut170403_20190121',outfile='Pisco_GNIRS_ut170403_20190121.fits')
#ech_coadd_new(giantcoadd=False,debug=True,datapath='/Users/feige/Dropbox/PypeIt_Redux/GNIRS/FSpec/Pisco/',
#              objinfo='pisco_info.txt',qafile='Pisco_GNIRS_20190121',outfile='Pisco_GNIRS_20190121.fits')


############################### Begin NIRES reduction
######### 2018 Oct NIRES data reduction
### J0252
#ech_flux_new(debug=False,datapath='/Users/feige/Dropbox/PypeIt_Redux/NIRES/J0252/ut181001/Science/',\
#            objinfo='J0252_info.txt',stdframe='spec1d_HIP13917_V8p6_NIRES_2018Oct01T125246.538.fits',telluric=True,
#            star_type='A0',star_mag=8.6,norder=5,resolution=3000,BALM_MASK_WID=100.0,polycorrect=True,polysens=False)
#ech_coadd_new(giantcoadd=False,debug=True,datapath='/Users/feige/Dropbox/PypeIt_Redux/NIRES/J0252/ut181001/Science/',\
#           objinfo='J0252_info.txt',qafile='J0252_NIRES_ut181001_20190122',outfile='J0252_NIRES_ut181001_20190122.fits')


############################### Begin XSHOOTER reduction
#ech_flux_new(debug=False,datapath='/Users/feige/Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Science/',\
#            objinfo='J0439_info.txt',stdframe='spec1d_STD,TELLURIC_XShooter_NIR_2018Oct08T232940.178.fits',
#            star_type='B8',star_mag=5.644)
#ech_coadd_new(giantcoadd=False,debug=True,datapath='/Users/feige/Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Science/',
#              objinfo='J0439_info.txt',qafile='ech_coadd',outfile='J0439_XSHOOTER.fits')

