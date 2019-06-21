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

def nires_fnames():
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

    return fnames, gdobj

def read_nires_stack(order=None):

    fnames, gdobj = nires_fnames()
    # parameters for load_1dspec_to_array
    ex_value = 'OPT'
    flux_value = True

    # Reading data
    waves,fluxes,ivars,masks,header = load.load_1dspec_to_array(fnames,gdobj=gdobj,order=order,ex_value=ex_value,
                                                                flux_value=flux_value)

    return waves, fluxes, ivars, masks,header

def xshooter_fnames():
    datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439_old/vlt_xshooter_nir/Science/')
    fnames = [datapath+'J0439_XSHOOTER_NIR_01.fits',datapath+'J0439_XSHOOTER_NIR_02.fits',datapath+'J0439_XSHOOTER_NIR_03.fits',
              datapath+'J0439_XSHOOTER_NIR_04.fits',datapath+'J0439_XSHOOTER_NIR_05.fits',datapath+'J0439_XSHOOTER_NIR_06.fits',
              datapath+'J0439_XSHOOTER_NIR_07.fits',datapath+'J0439_XSHOOTER_NIR_08.fits',datapath+'J0439_XSHOOTER_NIR_09.fits']

    gdobj = ['OBJ0001','OBJ0001','OBJ0001',
             'OBJ0001','OBJ0001','OBJ0001',
             'OBJ0001','OBJ0001','OBJ0001']

    return fnames, gdobj

def xshooter_fnames_newflux():
    datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439_old/vlt_xshooter_nir/Science/')
    fnames = [datapath+'J0439_XSHOOTER_NIR_01_flux.fits',datapath+'J0439_XSHOOTER_NIR_02_flux.fits',datapath+'J0439_XSHOOTER_NIR_03_flux.fits',
              datapath+'J0439_XSHOOTER_NIR_04_flux.fits',datapath+'J0439_XSHOOTER_NIR_05_flux.fits',datapath+'J0439_XSHOOTER_NIR_06_flux.fits',
              datapath+'J0439_XSHOOTER_NIR_07_flux.fits',datapath+'J0439_XSHOOTER_NIR_08_flux.fits',datapath+'J0439_XSHOOTER_NIR_09_flux.fits']

    gdobj = ['OBJ0001','OBJ0001','OBJ0001',
             'OBJ0001','OBJ0001','OBJ0001',
             'OBJ0001','OBJ0001','OBJ0001']

    return fnames, gdobj

def xshooter_fnames_newflux3():
    datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439_old/vlt_xshooter_nir/Science/')
    fnames = [datapath+'J0439_XSHOOTER_NIR_01_flux.fits',datapath+'J0439_XSHOOTER_NIR_02_flux.fits',
              datapath + 'J0439_XSHOOTER_NIR_03_flux.fits']
    gdobj = ['OBJ0001','OBJ0001','OBJ0001']

    return fnames, gdobj

def xshooter_fnames_newflux1():
    datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439_old/vlt_xshooter_nir/Science/')
    fnames = [datapath+'J0439_XSHOOTER_NIR_01_flux.fits']
    gdobj = ['OBJ0001']

    return fnames, gdobj


def xshooter_fnames_newflux1a():
    datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439_old/vlt_xshooter_nir/Science/')
    fnames = [datapath+'J0439_XSHOOTER_NIR_06_flux.fits']
    gdobj = ['OBJ0001']

    return fnames, gdobj

def feige110_xshooter_fnames():
    datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Science/')
    fnames = [datapath+'spec1d_XSHOO.2018-11-08T00:11:57.074-Feige110_XShooter_NIR_2018Nov08T001157.074._flux.fits',
              datapath+'spec1d_XSHOO.2018-11-08T00:16:56.583-Feige110_XShooter_NIR_2018Nov08T001656.583._flux.fits']

    gdobj = ['OBJ0001','OBJ0001']

    return fnames, gdobj

def J0020_xshooter_fnames():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0020-3653/NIR/Science/')
    fnames = [datapath + 'spec1d_VHSJ0020-3653OffsetstarB_XShooter_NIR_2017Dec17T024443.537_flux.fits',
              datapath + 'spec1d_VHSJ0020-3653OffsetstarB_XShooter_NIR_2017Dec17T030550.032_flux.fits',
              datapath + 'spec1d_VHSJ0020-3653OffsetstarB_XShooter_NIR_2017Oct26T001535.660_flux.fits',
              datapath + 'spec1d_VHSJ0020-3653OffsetstarB_XShooter_NIR_2017Oct26T002641.612_flux.fits']

    gdobj = ['OBJ0001','OBJ0001','OBJ0001','OBJ0001']
    return fnames, gdobj

def J0226_xshooter_fnames():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0226+0302/pypeit_nir/Science/')
    fnames = [datapath + 'spec1d_XSHOO.2017-12-17T03:50:47.125-PSOJ036.5078blindoffset_XShooter_NIR_2017Dec17T035047.125_flux.fits',
              datapath + 'spec1d_XSHOO.2017-12-17T04:11:13.716-PSOJ036.5078blindoffset_XShooter_NIR_2017Dec17T041113.716_flux.fits',
              datapath + 'spec1d_XSHOO.2018-01-14T02:12:34.014-PSOJ036.5078blindoffset_XShooter_NIR_2018Jan14T021234.014_flux.fits',
              datapath + 'spec1d_XSHOO.2018-01-14T02:33:00.603-PSOJ036.5078blindoffset_XShooter_NIR_2018Jan14T023300.603_flux.fits']

    gdobj = ['OBJ0001','OBJ0001','OBJ0001','OBJ0001']
    return fnames, gdobj


def J0224_xshooter_fnames():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0224-4711/pypeit_nir/Science/')
    fnames = [datapath + 'spec1d_XSHOO.2017-11-23T06:52:51.782-VDESJ0224-4711blindoffset_XShooter_NIR_2017Nov23T065251.782_flux.fits',
              datapath + 'spec1d_XSHOO.2017-11-23T07:13:18.374-VDESJ0224-4711blindoffset_XShooter_NIR_2017Nov23T071318.374_flux.fits',
              datapath + 'spec1d_XSHOO.2018-01-19T01:57:51.708-VDESJ0224-4711blindoffset_XShooter_NIR_2018Jan19T015751.708_flux.fits',
              datapath + 'spec1d_XSHOO.2018-01-19T02:18:18.297-VDESJ0224-4711blindoffset_XShooter_NIR_2018Jan19T021818.297_flux.fits']

    gdobj = ['OBJ0001','OBJ0001','OBJ0001','OBJ0001']
    return fnames, gdobj

def J0224_xshooter_fnamesa():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0224-4711/pypeit_nir/Science/')
    fnames = [datapath + 'spec1d_XSHOO.2017-11-23T06:52:51.782-VDESJ0224-4711blindoffset_XShooter_NIR_2017Nov23T065251.782_flux.fits',
              datapath + 'spec1d_XSHOO.2017-11-23T07:13:18.374-VDESJ0224-4711blindoffset_XShooter_NIR_2017Nov23T071318.374_flux.fits']
    gdobj = ['OBJ0001','OBJ0001']
    return fnames, gdobj

def J0224_xshooter_fnamesb():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0224-4711/pypeit_nir/Science/')
    fnames = [datapath + 'spec1d_XSHOO.2018-01-19T01:57:51.708-VDESJ0224-4711blindoffset_XShooter_NIR_2018Jan19T015751.708_flux.fits',
              datapath + 'spec1d_XSHOO.2018-01-19T02:18:18.297-VDESJ0224-4711blindoffset_XShooter_NIR_2018Jan19T021818.297_flux.fits']

    gdobj = ['OBJ0001','OBJ0001']
    return fnames, gdobj

def LTT_xshooter_frames():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0020-3653/NIR/Science/')
    fnames = [datapath + 'spec1d_STD,FLUX_XShooter_NIR_2017Dec17T081653.582_flux.fits',
              datapath + 'spec1d_STD,FLUX_XShooter_NIR_2017Dec17T082243.751_flux.fits']
    gdobj = ['OBJ0001','OBJ0001']

    return fnames, gdobj

def J1048_xshooter_fnames():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/OBS_DATA/XSHOOTER/NIR/ut20170202/Science/')
    fnames = [datapath + 'spec1d_XSHOO.2017-02-02T04:19:28.545-VIKJ1048m0109_XShooter_NIR_2017Feb02T041928.545_flux.fits',
              datapath + 'spec1d_XSHOO.2017-02-02T04:40:14.422-VIKJ1048m0109_XShooter_NIR_2017Feb02T044014.422_flux.fits',
              datapath + 'spec1d_XSHOO.2017-02-02T05:17:52.162-VIKJ1048m0109_XShooter_NIR_2017Feb02T051752.162_flux.fits']
    gdobj = ['OBJ0001', 'OBJ0001', 'OBJ0001']
    return fnames, gdobj

def read_xshooter_nir_stack(order=None):

    fnames, gdobj = xshooter_fnames()

    # parameters for load_1dspec_to_array
    ex_value = 'OPT'
    flux_value = True

    # Reading data
    waves,fluxes,ivars,masks,header = load.load_1dspec_to_array(fnames,gdobj=gdobj,order=order,ex_value=ex_value,
                                                                flux_value=flux_value)

    return waves, fluxes, ivars, masks, header


def deimos_fnames():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/DEIMOS/Science/')
    fnames = [datapath + 'spec1d_F_DE.20170527.37601-P261_OFF_DEIMOS_2017May27T102635.318.fits',
              datapath + 'spec1d_F_DE.20170527.38872-P261_OFF_DEIMOS_2017May27T104746.349.fits',
              datapath + 'spec1d_F_DE.20170527.41775-P261_OFF_DEIMOS_2017May27T113608.093.fits',
              datapath + 'spec1d_F_DE.20170527.43045-P261_OFF_DEIMOS_2017May27T115718.864.fits',
              datapath + 'spec1d_F_DE.20170527.44316-P261_OFF_DEIMOS_2017May27T121830.586.fits']
    objids = ['SPAT0764-SLIT0000-DET07',
              'SPAT0764-SLIT0000-DET07',
              'SPAT0758-SLIT0000-DET07',
              'SPAT0758-SLIT0000-DET07',
              'SPAT0758-SLIT0000-DET07']

    return fnames, objids

def read_deimos_stack():
    fnames, objids = deimos_fnames()
    # parameters for load_1dspec_to_array
    ex_value = 'OPT'
    flux_value = True

    # Reading data
    waves, fluxes, ivars, masks, header = load.load_1dspec_to_array(fnames,gdobj=objids,order=None,ex_value=ex_value,
                                                                    flux_value=flux_value)
    return waves, fluxes, ivars, masks, header



#### Longslit coadd test
#waves, fluxes, ivars, masks, header = read_gmos_stack()
#waves, fluxes, ivars, masks = read_lris_stack()
#waves, fluxes, ivars, masks, header = read_nires_stack(order=4)
#waves, fluxes, ivars, masks, header = read_xshooter_nir_stack(order=13)
#waves, fluxes, ivars, masks, header = read_deimos_stack()

# Generate a wave_grid
#wave_grid = coadd1d.new_wave_grid(waves, wave_method='pixel')

# Testing DEIMOS
#fnames, objids = deimos_fnames()
#wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.multi_combspec(fnames, objids, show=True, debug=True,
#                                                                        outfile='P261_coadd.fits')

# Testing NIRES
#fnames, objids = nires_fnames()
#wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames, objids, show=True)


# Test XSHOOTER
#sensfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Feige110_sens_tell.fits')
sensfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/NIR_Stack/Feige110_sens_tell_20190620.fits')

#fnames, objids = J0226_xshooter_fnames()
#outfile = 'J0226.fits'
#fnames, objids = xshooter_fnames_newflux3()
#outfile = 'J0439_3exp'
#fnames, objids = xshooter_fnames_newflux()
#fnames, objids = xshooter_fnames_newflux1a()
#fnames, objids = J0020_xshooter_fnames()
#outfile = 'J0020'
#fnames, objids = feige110_xshooter_fnames()
#fnames, objids = J0224_xshooter_fnames()
#outfile = 'J0224'
#fnames, objids = LTT_xshooter_frames()
#outfile = 'LTT3218'
fnames, objids = J1048_xshooter_fnames()
outfile = 'J1048'

wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames, objids, show=True, sensfile=sensfile,
                                                                      outfile=outfile)

# Coadding
#wave_stack, flux_stack, ivar_stack, mask_stack, outmask, weights, scales, rms_sn = coadd1d.combspec(
#    wave_grid, waves, fluxes, ivars, masks, debug=True)
