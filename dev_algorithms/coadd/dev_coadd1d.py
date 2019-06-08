import os
import coadd1d
import numpy as np
import os
from pypeit import utils
from astropy.table import Table
from astropy.io import fits
from matplotlib import pyplot as plt
import sys
import IPython

def read_lris_stack():
    stackfile = '/Users/joe/REDUX/lris_redux/Nov_2004/Final/SDSSJ073522.43+295710.1_N.fits'
    hdu = fits.open(stackfile)
    flux_ref = hdu[0].data
    sig_ref = hdu[1].data
    ivar_ref = utils.calc_ivar(sig_ref)
    mask_ref = (ivar_ref > 0.0)
    wave_ref = hdu[2].data

    infiles = ['/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0063.fits.gz',
               '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0064.fits.gz',
               '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0065.fits.gz',
               '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0066.fits.gz',
               '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0219.fits.gz',
               '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0220.fits.gz',
               '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0221.fits.gz',
               '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0222.fits.gz',
               '/Users/joe/REDUX/lris_redux/Nov_2004/0735+2957/Blue/Science/sci-lblue0063.fits.gz']
    nfiles = len(infiles)
    objid = 0
    for idx, file in enumerate(infiles):
        obj = Table.read(file, hdu=5)
        flux = obj[objid]['FLUX_OPT']
        ivar = obj[objid]['IVAR_OPT']
        wave = obj[objid]['WAVE_OPT']
        if idx == 0:
            nspec = flux.size
            flux_arr = np.zeros((nfiles, nspec))
            wave_arr = np.zeros((nfiles, nspec))
            ivar_arr = np.zeros((nfiles, nspec))
        flux_arr[idx, :] = np.flip(flux)
        ivar_arr[idx, :] = np.flip(ivar)
        wave_arr[idx, :] = np.flip(wave)

    mask_arr = (ivar_arr > 0.0)

    return wave_arr, flux_arr, ivar_arr, mask_arr

def read_gmos_stack():
    datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/GMOS/R400_Flux/')
    fnames = [datapath+'spec1d_flux_S20180903S0136-J0252-0503_GMOS-S_1864May27T160716.387.fits',\
              datapath+'spec1d_flux_S20180903S0137-J0252-0503_GMOS-S_1864May27T160719.968.fits',\
              datapath+'spec1d_flux_S20180903S0138-J0252-0503_GMOS-S_1864May27T160723.353.fits',\
              datapath+'spec1d_flux_S20180903S0141-J0252-0503_GMOS-S_1864May27T160727.033.fits',\
              datapath+'spec1d_flux_S20180903S0142-J0252-0503_GMOS-S_1864May27T160730.419.fits',\
              datapath+'spec1d_flux_S20181015S0140-J0252-0503_GMOS-S_1864May27T185252.770.fits']
    gdobj = ['SPAT1073-SLIT0001-DET03','SPAT1167-SLIT0001-DET03','SPAT1071-SLIT0001-DET03','SPAT1072-SLIT0001-DET03',
             'SPAT1166-SLIT0001-DET03','SPAT1073-SLIT0001-DET03']

    # parameters for load_1dspec_to_array
    ex_value = 'OPT'
    flux_value = True

    # Reading data
    waves,fluxes,ivars,masks = coadd1d.load_1dspec_to_array(fnames,gdobj=gdobj,order=None,ex_value=ex_value,flux_value=flux_value)

    return waves, fluxes, ivars, masks

waves, fluxes, ivars, masks = read_gmos_stack()
#waves, fluxes, ivars, masks = read_lris_stack()

# Coadding
wave_stack, flux_stack, ivar_stack, mask_stack, scale_array = \
    coadd1d.long_comb(waves, fluxes, ivars, masks, wave_method='pixel', scale_method='poly', maxiter_reject = 5, \
                      qafile='J0252_gmos', outfile='J0252_gmos.fits', verbose=False, debug=True)

#