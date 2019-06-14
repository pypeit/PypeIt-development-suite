import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#from pypeit.core import coadd1d
#from coadd1d_old import *
#from pypeit.core.coadd1d import *
from pypeit.core import load, save
from pypeit.core import coadd1d

## GMOS test
datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/GMOS/R400_Flux/')
fnames = [datapath + 'spec1d_flux_S20180903S0136-J0252-0503_GMOS-S_1864May27T160716.387.fits', \
          datapath + 'spec1d_flux_S20180903S0137-J0252-0503_GMOS-S_1864May27T160719.968.fits', \
          datapath + 'spec1d_flux_S20180903S0138-J0252-0503_GMOS-S_1864May27T160723.353.fits', \
          datapath + 'spec1d_flux_S20180903S0141-J0252-0503_GMOS-S_1864May27T160727.033.fits', \
          datapath + 'spec1d_flux_S20180903S0142-J0252-0503_GMOS-S_1864May27T160730.419.fits', \
          datapath + 'spec1d_flux_S20181015S0140-J0252-0503_GMOS-S_1864May27T185252.770.fits']
objids = ['SPAT1073-SLIT0001-DET03', 'SPAT1167-SLIT0001-DET03', 'SPAT1071-SLIT0001-DET03', 'SPAT1072-SLIT0001-DET03',
         'SPAT1166-SLIT0001-DET03', 'SPAT1073-SLIT0001-DET03']

# parameters for load_1dspec_to_array
ex_value = 'OPT'
flux_value = True

#wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.multi_combspec(fnames, objids, debug=False)


### NIRES test
datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/NIRES/NIRES_May19/Science/')
fnames = [datapath + 'spec1d_flux_s190519_0037-J1007+2115_NIRES_2019May19T055221.895.fits',
          datapath + 'spec1d_flux_s190519_0038-J1007+2115_NIRES_2019May19T055923.665.fits',
          datapath + 'spec1d_flux_s190519_0041-J1007+2115_NIRES_2019May19T062048.865.fits',
          datapath + 'spec1d_flux_s190519_0042-J1007+2115_NIRES_2019May19T062750.635.fits',
          datapath + 'spec1d_flux_s190519_0045-J1007+2115_NIRES_2019May19T064943.885.fits',
          datapath + 'spec1d_flux_s190519_0046-J1007+2115_NIRES_2019May19T065646.165.fits',
          datapath + 'spec1d_flux_s190519_0049-J1007+2115_NIRES_2019May19T071920.215.fits',
          datapath + 'spec1d_flux_s190519_0050-J1007+2115_NIRES_2019May19T072621.985.fits',
          datapath + 'spec1d_flux_s190519_0053-J1007+2115_NIRES_2019May19T074819.315.fits',
          datapath + 'spec1d_flux_s190519_0054-J1007+2115_NIRES_2019May19T075521.595.fits',
          datapath + 'spec1d_flux_s190519_0057-J1007+2115_NIRES_2019May19T081918.265.fits',
          datapath + 'spec1d_flux_s190519_0058-J1007+2115_NIRES_2019May19T082620.545.fits']

objids = ['OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
         'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001',
         'OBJ0001', 'OBJ0001', 'OBJ0001', 'OBJ0001']

wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames, objids,debug=False)
