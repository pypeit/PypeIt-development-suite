import os
import numpy as np

import telluric
from flux1d import apply_sensfunc
from pypeit.core import coadd1d

debug = False
show = True
z_qso = 7.085
npca = 8
ex_value = 'OPT'
pca_file = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux//qso_pca_1200_3100.pckl')
telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/TelFit_Paranal_NIR_9800_25000_R25000.fits')
sensfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/LTT3218_sens_tell.fits')

qsoname = 'J1120+0641'


datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/XSHOOTER/{:}/NIR/Science/'.format(qsoname))
#TODO: change the spec1dlist to the pypeit format and change the reader accordingly
spec1dfiles = np.genfromtxt(os.path.join(datapath,'spec1dlist'),dtype='str')

nfiles = len(spec1dfiles)
fnames = []
for ifile in range(nfiles):
    fnames.append(os.path.join(datapath,spec1dfiles[ifile]))
#TODO: the objids shoul be read in from the pypeit format file as noted above.
objids = ['OBJ0001']*nfiles

# apply the sensfunc to all spectra
apply_sensfunc(fnames, sensfile, extinct_correct=False, tell_correct=False, debug=debug, show=False)

# let's coadd all the fluxed spectra
wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames, objids, show=True, sensfile=sensfile,
                                                                      ex_value='OPT', outfile=qsoname, debug=debug)

# run telluric.qso_telluric to get the final results
spec1dfluxfile = 'spec1d_stack_{:}.fits'.format(qsoname)
telloutfile = 'spec1d_stack_{:}_tellmodel.fits'.format(qsoname)
outfile = 'spec1d_stack_{:}_tellcorr.fits'.format(qsoname)

TelQSO = telluric.qso_telluric(spec1dfluxfile, telgridfile, pca_file, z_qso, telloutfile, outfile,
                               create_bal_mask=None, debug=debug, show=show)