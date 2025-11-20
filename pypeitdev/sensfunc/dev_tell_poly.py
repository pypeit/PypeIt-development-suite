import os
import numpy as np

import telluric
from pypeit.core.flux_calib import apply_sensfunc
from pypeit.core import coadd1d
from pypeit import msgs
show=True


#### Optical
telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_VIS_4900_11100_R25000.fits')

z_obj = 6.32
objname = 'J0100+2802'

#z_obj = 6.51
#objname = 'J0224-4711'

# run telluric.poly_telluric to get the final results
spec1dfluxfile = 'spec1d_stack_{:}_newcoadd.fits'.format(objname)
telloutfile = 'spec1d_stack_{:}_newcoadd_tellmodel.fits'.format(objname)
outfile = 'spec1d_stack_{:}_newcoadd_tellcorr.fits'.format(objname)

TelPoly = telluric.poly_telluric(spec1dfluxfile, telgridfile, telloutfile, outfile, z_obj=z_obj, polyorder=3,
                                 fit_region_min=[9200.0], fit_region_max=[9700.0], func='legendre', model='exp',
                                 mask_lyman_a=True, debug_init=True, debug=True, show=show)

##### NIR
telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
# run telluric.poly_telluric to get the final results
spec1dfluxfile = 'spec1d_stack_J0439_2exp.fits'
telloutfile = 'spec1d_stack_J0439_2exp_tellmodel.fits'
outfile = 'spec1d_stack_J0439_2exp_tellcorr.fits'

#TelPoly = telluric.poly_telluric(spec1dfluxfile, telgridfile, telloutfile, outfile, z_obj=z_obj, polyorder=3,
#                                 fit_region_min=[11150.0], fit_region_max=[11600.0], func='legendre', model='exp',
#                                 mask_lyman_a=True, debug_init=True, debug=True, show=show)

#TelPoly = telluric.poly_telluric(spec1dfluxfile, telgridfile, telloutfile, outfile, z_obj=z_obj, polyorder=3,
#                                 fit_region_min=[11150.0,19940.0], fit_region_max=[11600.0,20800.0], func='legendre', model='exp',
#                                 mask_lyman_a=True, debug_init=True, debug=True, show=show)


