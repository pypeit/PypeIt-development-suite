import os
#from pypeit.core import coadd1d
from pypeit.core import coadd1d
from pypeit.core import load
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
            flux_arr = np.zeros((nspec, nfiles))
            wave_arr = np.zeros((nspec, nfiles))
            ivar_arr = np.zeros((nspec, nfiles))
        flux_arr[:, idx] = np.flip(flux)
        ivar_arr[:, idx] = np.flip(ivar)
        wave_arr[:, idx] = np.flip(wave)

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
    waves,fluxes,ivars,masks,header = load.load_1dspec_to_array(fnames,gdobj=gdobj,order=None,ex_value=ex_value,
                                                                flux_value=flux_value)

    return waves, fluxes, ivars, masks, header

def read_nires_stack(order=None):
    datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/NIRES/NIRES_May19/Science/')
    fnames = [datapath+'spec1d_flux_s190519_0037-J1007+2115_NIRES_2019May19T055221.895.fits',
              datapath+'spec1d_flux_s190519_0038-J1007+2115_NIRES_2019May19T055923.665.fits',
              datapath+'spec1d_flux_s190519_0041-J1007+2115_NIRES_2019May19T062048.865.fits',
              datapath+'spec1d_flux_s190519_0042-J1007+2115_NIRES_2019May19T062750.635.fits',
              datapath+'spec1d_flux_s190519_0045-J1007+2115_NIRES_2019May19T064943.885.fits',
              datapath+'spec1d_flux_s190519_0046-J1007+2115_NIRES_2019May19T065646.165.fits',
              datapath+'spec1d_flux_s190519_0049-J1007+2115_NIRES_2019May19T071920.215.fits',
              datapath+'spec1d_flux_s190519_0050-J1007+2115_NIRES_2019May19T072621.985.fits',
              datapath+'spec1d_flux_s190519_0053-J1007+2115_NIRES_2019May19T074819.315.fits',
              datapath+'spec1d_flux_s190519_0054-J1007+2115_NIRES_2019May19T075521.595.fits',
              datapath+'spec1d_flux_s190519_0057-J1007+2115_NIRES_2019May19T081918.265.fits',
              datapath+'spec1d_flux_s190519_0058-J1007+2115_NIRES_2019May19T082620.545.fits']

    gdobj = ['OBJ0001','OBJ0001','OBJ0001','OBJ0001',
             'OBJ0001','OBJ0001','OBJ0001','OBJ0001',
             'OBJ0001', 'OBJ0001','OBJ0001','OBJ0001']

    # parameters for load_1dspec_to_array
    ex_value = 'OPT'
    flux_value = True

    # Reading data
    waves,fluxes,ivars,masks,header = load.load_1dspec_to_array(fnames,gdobj=gdobj,order=order,ex_value=ex_value,
                                                                flux_value=flux_value)

    return waves, fluxes, ivars, masks,header

def read_xshooter_nir_stack(order=None):
    datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439_old/vlt_xshooter_nir/Science/')
    fnames = [datapath+'J0439_XSHOOTER_NIR_01.fits',datapath+'J0439_XSHOOTER_NIR_02.fits',datapath+'J0439_XSHOOTER_NIR_03.fits',
              datapath+'J0439_XSHOOTER_NIR_04.fits',datapath+'J0439_XSHOOTER_NIR_05.fits',datapath+'J0439_XSHOOTER_NIR_06.fits',
              datapath+'J0439_XSHOOTER_NIR_07.fits',datapath+'J0439_XSHOOTER_NIR_08.fits',datapath+'J0439_XSHOOTER_NIR_09.fits']

    gdobj = ['OBJ0001','OBJ0001','OBJ0001',
             'OBJ0001','OBJ0001','OBJ0001',
             'OBJ0001','OBJ0001','OBJ0001']

    # parameters for load_1dspec_to_array
    ex_value = 'OPT'
    flux_value = True

    # Reading data
    waves,fluxes,ivars,masks,header = load.load_1dspec_to_array(fnames,gdobj=gdobj,order=order,ex_value=ex_value,
                                                                flux_value=flux_value)

    return waves, fluxes, ivars, masks, header


def read_deimos_stack():
    datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/DEIMOS/Science/')
    fnames = [datapath + 'spec1d_F_DE.20170527.37601-P261_OFF_DEIMOS_2017May27T102635.318.fits',
              datapath + 'spec1d_F_DE.20170527.38872-P261_OFF_DEIMOS_2017May27T104746.349.fits',
              datapath + 'spec1d_F_DE.20170527.41775-P261_OFF_DEIMOS_2017May27T113608.093.fits',
              datapath + 'spec1d_F_DE.20170527.43045-P261_OFF_DEIMOS_2017May27T115718.864.fits',
              datapath + 'spec1d_F_DE.20170527.44316-P261_OFF_DEIMOS_2017May27T121830.586.fits']
    gdobj = ['SPAT0764-SLIT0000-DET07',
             'SPAT0764-SLIT0000-DET07',
             'SPAT0758-SLIT0000-DET07',
             'SPAT0758-SLIT0000-DET07',
             'SPAT0758-SLIT0000-DET07']

    # parameters for load_1dspec_to_array
    ex_value = 'OPT'
    flux_value = True

    # Reading data
    waves,fluxes,ivars,masks,header = load.load_1dspec_to_array(fnames,gdobj=gdobj,order=None,ex_value=ex_value,
                                                                flux_value=flux_value)

    return waves, fluxes, ivars, masks, header


#### Longslit coadd test
#waves, fluxes, ivars, masks, header = read_gmos_stack()
#waves, fluxes, ivars, masks = read_lris_stack()
#waves, fluxes, ivars, masks, header = read_nires_stack(order=4)
waves, fluxes, ivars, masks, header = read_xshooter_nir_stack(order=13)
#waves, fluxes, ivars, masks = read_deimos_stack()

# Generate a wave_grid
wave_grid = coadd1d.new_wave_grid(waves, wave_method='pixel')

# Coadding
wave_stack, flux_stack, ivar_stack, mask_stack, outmask, weights, scales, rms_sn = coadd1d.long_combspec(
    wave_grid, waves, fluxes, ivars, masks, qafile='P261_deimos', outfile='P261_deimos.fits', debug=True)
