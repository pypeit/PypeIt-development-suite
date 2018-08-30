



from pypeit import msgs
import numpy as np
from IPython import embed
from pypeit import flatfield
from astropy.io import fits
import os
from pypeit.spectrographs.util import load_spectrograph


type = 'LRIS_red'
devpath = os.getenv('PYPEIT_DEV')

if type == 'LRIS_red':
    rawpath = devpath + '/RAW_DATA/Keck_LRIS_red/multi_400_8500_d560/'
    masterpath = devpath + '/REDUX_OUT/Keck_LRIS_red/multi_400_8500_d560/MF_keck_lris_red/'
    biasfile = masterpath + 'MasterBias_A_02_aa.fits'
    msbias = fits.getdata(biasfile)
    pixflat_image_files = np.core.defchararray.add(rawpath, ['r170320_2057.fits','r170320_2058.fits','r170320_2059.fits']).tolist()

    det = 2
    spectro_name = 'keck_lris_red'
    spectrograph = load_spectrograph(spectrograph=spectro_name)
    par = spectrograph.default_pypeit_par()
    flatField = flatfield.FlatField(spectrograph, file_list=pixflat_image_files,det=det, par=par['calibrations']['pixelflatframe']
                                    , msbias = msbias)
