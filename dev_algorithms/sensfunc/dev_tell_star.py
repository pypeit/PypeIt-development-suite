import os
import numpy as np

import telluric
from pypeit.core.flux_calib import apply_sensfunc
from pypeit.core import coadd1d
from pypeit import msgs
show=True

objname = 'Hip022840'
star_mag = 5.965 # V-band
star_type = 'B5' # B5/7II D spectroscopic binary
dev_path = os.getenv('PYPEIT_DEV')
datapath = os.path.join(dev_path, 'REDUX_OUT/VLT_XSHOOTER/NIR/Science')

#spec1dfiles = #['spec1d_XSHOO.2016-08-02T09:57:17.147-STD,TELLURIC_XShooter_NIR_2016Aug02T095717.147.fits',#
spec1dfiles =  ['spec1d_XSHOO.2016-08-02T09:57:56.903-STD,TELLURIC_XShooter_NIR_2016Aug02T095756.903.fits']

fnames = [os.path.join(datapath, file) for file in spec1dfiles]
nfiles = len(fnames)
objids = ['OBJ0001']*nfiles

sensfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/LTT3218_sens_tell.fits')

## Apply the sensfunc to the telluric star spectra
#apply_sensfunc(fnames, sensfile, extinct_correct=False, tell_correct=False, debug=False, show=True)

fnames_flux = [f.replace('.fits', '_flux.fits') for f in fnames]

# TODO: change the outfile to work with datapath. It's a hard coding on these names in coadd1d
#wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames_flux, objids, show=True, sensfile=sensfile,
#                                                                      ex_value='OPT', outfile=objname,
#                                                                      show_order_scale=True,
#                                                                      debug=True)

#
# run telluric.qso_telluric to get the final results
telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
spec1dfluxfile = 'spec1d_stack_{:}.fits'.format(objname)
telloutfile = 'spec1d_stack_{:}_tellmodel.fits'.format(objname)
outfile = 'spec1d_stack_{:}_tellcorr.fits'.format(objname)

TelStar = telluric.star_telluric(spec1dfluxfile, telgridfile, telloutfile, outfile, polyorder=7, debug_init=True, debug=True, show=show,
                                 star_mag=star_mag, star_type=star_type)
