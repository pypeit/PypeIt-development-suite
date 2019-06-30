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

def TELL_xshooter_frames():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0224-4711/Test_tell/')
    fnames = [#datapath + 'spec1d_XSHOO.2017-11-23T07:44:02.747-STD,TELLURIC_XShooter_NIR_2017Nov23T074402.747_flux.fits',
              #datapath + 'spec1d_XSHOO.2017-11-23T07:44:57.633-STD,TELLURIC_XShooter_NIR_2017Nov23T074457.633_flux.fits',
              datapath + 'spec1d_XSHOO.2017-11-23T07:46:53.532-STD,TELLURIC_XShooter_NIR_2017Nov23T074653.532_flux.fits',
              datapath + 'spec1d_XSHOO.2017-11-23T07:47:33.917-STD,TELLURIC_XShooter_NIR_2017Nov23T074733.917_flux.fits']
    gdobj = ['OBJ0001','OBJ0001','OBJ0001','OBJ0001']

    return fnames, gdobj

def J1048_xshooter_fnames():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/OBS_DATA/XSHOOTER/NIR/ut20170202/Science/')
    fnames = [datapath + 'spec1d_XSHOO.2017-02-02T04:19:28.545-VIKJ1048m0109_XShooter_NIR_2017Feb02T041928.545_flux.fits',
              datapath + 'spec1d_XSHOO.2017-02-02T04:40:14.422-VIKJ1048m0109_XShooter_NIR_2017Feb02T044014.422_flux.fits',
              datapath + 'spec1d_XSHOO.2017-02-02T05:17:52.162-VIKJ1048m0109_XShooter_NIR_2017Feb02T051752.162_flux.fits']
    gdobj = ['OBJ0001', 'OBJ0001', 'OBJ0001']
    return fnames, gdobj

def pisco_xshooter_fnames():
    datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/Pypeit_files/PISCO_NIR/Science/')
    fnames = [datapath + 'spec1d_XSHOO.2017-06-28T23:51:39.115-PSOJ205p09_1_XShooter_NIR_2017Jun28T235139.115_flux.fits',
        datapath + 'spec1d_XSHOO.2017-06-29T00:12:24.992-PSOJ205p09_1_XShooter_NIR_2017Jun29T001224.992_flux.fits',
        datapath + 'spec1d_XSHOO.2018-02-15T08:07:00.164-PSOJ205p09_2_XShooter_NIR_2018Feb15T080700.164_flux.fits',
        datapath + 'spec1d_XSHOO.2018-02-15T08:27:44.713-PSOJ205p09_2_XShooter_NIR_2018Feb15T082744.713_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-10T05:08:40.055-PSOJ205p09_3_XShooter_NIR_2018Mar10T050840.055_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-10T05:29:25.933-PSOJ205p09_3_XShooter_NIR_2018Mar10T052925.933_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-10T05:58:33.093-PSOJ205p09_4_XShooter_NIR_2018Mar10T055833.093_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-10T06:19:22.960-PSOJ205p09_4_XShooter_NIR_2018Mar10T061922.960_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-11T08:41:01.429-PSOJ205p09_5_XShooter_NIR_2018Mar11T084101.429_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-11T09:01:47.309-PSOJ205p09_5_XShooter_NIR_2018Mar11T090147.309_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-18T06:26:29.402-PSOJ205p09_6_XShooter_NIR_2018Mar18T062629.402_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-18T06:47:15.945-PSOJ205p09_6_XShooter_NIR_2018Mar18T064715.945_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-18T08:32:32.302-PSOJ205p09_7_XShooter_NIR_2018Mar18T083232.302_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-18T08:53:18.845-PSOJ205p09_7_XShooter_NIR_2018Mar18T085318.845_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-20T05:22:01.367-PSOJ205p09_8_XShooter_NIR_2018Mar20T052201.367_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-20T05:42:48.575-PSOJ205p09_8_XShooter_NIR_2018Mar20T054248.575_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-21T06:06:09.751-PSOJ205p09_9_XShooter_NIR_2018Mar21T060609.751_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-21T06:26:55.631-PSOJ205p09_9_XShooter_NIR_2018Mar21T062655.631_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-21T07:13:08.797-PSOJ205p09_10_XShooter_NIR_2018Mar21T071308.797_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-21T07:33:56.008-PSOJ205p09_10_XShooter_NIR_2018Mar21T073356.008_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-22T05:00:30.229-PSOJ205p09_11_XShooter_NIR_2018Mar22T050030.229_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-22T05:21:17.439-PSOJ205p09_11_XShooter_NIR_2018Mar22T052117.439_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-22T06:01:30.465-PSOJ205p09_12_XShooter_NIR_2018Mar22T060130.465_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-22T06:22:17.675-PSOJ205p09_12_XShooter_NIR_2018Mar22T062217.675_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-23T04:45:05.305-PSOJ205p09_13_XShooter_NIR_2018Mar23T044505.305_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-23T05:05:50.520-PSOJ205p09_13_XShooter_NIR_2018Mar23T050550.520_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-23T05:43:50.910-PSOJ205p09_14_XShooter_NIR_2018Mar23T054350.910_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-23T06:04:37.452-PSOJ205p09_14_XShooter_NIR_2018Mar23T060437.452_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-23T06:37:00.347-PSOJ205p09_15_XShooter_NIR_2018Mar23T063700.347_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-23T06:57:45.561-PSOJ205p09_16_XShooter_NIR_2018Mar23T065745.561_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-23T07:36:36.220-PSOJ205p09_16_XShooter_NIR_2018Mar23T073636.220_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-23T07:57:23.428-PSOJ205p09_16_XShooter_NIR_2018Mar23T075723.428_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-24T05:20:10.198-PSOJ205p09_17_XShooter_NIR_2018Mar24T052010.198_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-24T05:40:56.740-PSOJ205p09_17_XShooter_NIR_2018Mar24T054056.740_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-24T06:17:18.332-PSOJ205p09_18_XShooter_NIR_2018Mar24T061718.332_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-24T06:38:05.542-PSOJ205p09_18_XShooter_NIR_2018Mar24T063805.542_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-24T07:13:35.465-PSOJ205p09_19_XShooter_NIR_2018Mar24T071335.465_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-24T07:34:22.671-PSOJ205p09_19_XShooter_NIR_2018Mar24T073422.671_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-25T05:18:51.790-PSOJ205p09_20_XShooter_NIR_2018Mar25T051851.790_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-25T05:39:38.333-PSOJ205p09_20_XShooter_NIR_2018Mar25T053938.333_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-26T06:36:32.911-PSOJ205p09_21_XShooter_NIR_2018Mar26T063632.911_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-26T06:57:20.118-PSOJ205p09_21_XShooter_NIR_2018Mar26T065720.118_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-27T06:57:29.269-PSOJ205p09_22_XShooter_NIR_2018Mar27T065729.269_flux.fits',
        datapath + 'spec1d_XSHOO.2018-03-27T07:18:15.146-PSOJ205p09_22_XShooter_NIR_2018Mar27T071815.146_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-08T04:59:57.666-PSOJ205p09_23_XShooter_NIR_2018Apr08T045957.666_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-08T05:20:44.206-PSOJ205p09_23_XShooter_NIR_2018Apr08T052044.206_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-08T05:57:54.622-PSOJ205p09_24_XShooter_NIR_2018Apr08T055754.622_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-08T06:18:42.492-PSOJ205p09_24_XShooter_NIR_2018Apr08T061842.492_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-09T04:59:20.273-PSOJ205p09_25_XShooter_NIR_2018Apr09T045920.273_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-09T05:20:06.817-PSOJ205p09_25_XShooter_NIR_2018Apr09T052006.817_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-09T05:58:52.879-PSOJ205p09_26_XShooter_NIR_2018Apr09T055852.879_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-09T06:19:39.420-PSOJ205p09_26_XShooter_NIR_2018Apr09T061939.420_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-10T05:07:58.959-PSOJ205p09_27_XShooter_NIR_2018Apr10T050758.959_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-10T05:28:46.165-PSOJ205p09_27_XShooter_NIR_2018Apr10T052846.165_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-10T06:05:42.360-PSOJ205p09_28_XShooter_NIR_2018Apr10T060542.360_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-10T06:26:30.231-PSOJ205p09_28_XShooter_NIR_2018Apr10T062630.231_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-15T04:58:27.291-PSOJ205p09_29_XShooter_NIR_2018Apr15T045827.291_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-15T05:19:15.167-PSOJ205p09_29_XShooter_NIR_2018Apr15T051915.167_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-15T05:59:39.364-PSOJ205p09_30_XShooter_NIR_2018Apr15T055939.364_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-15T06:20:25.908-PSOJ205p09_30_XShooter_NIR_2018Apr15T062025.908_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-16T04:23:25.898-PSOJ205p09_31_XShooter_NIR_2018Apr16T042325.898_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-16T04:44:12.441-PSOJ205p09_31_XShooter_NIR_2018Apr16T044412.441_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-16T05:26:58.617-PSOJ205p09_32_XShooter_NIR_2018Apr16T052658.617_flux.fits',
        datapath + 'spec1d_XSHOO.2018-04-16T05:47:45.160-PSOJ205p09_32_XShooter_NIR_2018Apr16T054745.160_flux.fits']
    gdobj = ['OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
             'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
             'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
             'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
             'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
             'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
             'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
             'OBJ0001']
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
fnames, objids = deimos_fnames()
wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.combspec(fnames, objids, show=True, debug=True,outfile='P261_coadd.fits')

# Testing NIRES
#fnames, objids = nires_fnames()
#wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames, objids, show=True)


# Test XSHOOTER
#sensfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/J0439/NIR/Feige110_sens_tell.fits')
#sensfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/NIR_Stack/Feige110_sens_tell_wang.fits')

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
#fnames, objids = J1048_xshooter_fnames()
#outfile = 'J1048'
#fnames, objids = pisco_xshooter_fnames()
#outfile = 'Pisco_all_box'
#fnames, objids = TELL_xshooter_frames()
#outfile = 'TELL_B8IV_V5p8'

#wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames, objids, show=True, sensfile=sensfile,
#                                                                      ex_value='OPT', outfile=outfile, debug=False)

# Coadding
#wave_stack, flux_stack, ivar_stack, mask_stack, outmask, weights, scales, rms_sn = coadd1d.combspec(
#    wave_grid, waves, fluxes, ivars, masks, debug=True)
