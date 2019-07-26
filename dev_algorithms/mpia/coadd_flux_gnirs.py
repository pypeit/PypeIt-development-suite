
import os
import numpy as np

from pypeit.core import telluric
from pypeit.core.flux_calib import apply_sensfunc
from pypeit.core import coadd1d
from pypeit import msgs
import glob

clobber=False
debug = False
show = True

z_qso = 7.54
npca = 8
ex_value = 'OPT'
qsoname = 'pisco'

mpia_path = os.path.join(os.getenv('HOME'), 'Dropbox/MPIA_PypeIt')
dev_suite_path =  os.path.join(mpia_path, 'work/gemini_gnirs_A/Science')
#dev_suite_path = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT/Gemini_GNIRS/GNIRS/Science')

# List of science files for pisco
file_list = glob.glob(dev_suite_path + '/spec1d*pisco*.fits')
spec1dfiles = []
for ifile in file_list:
    if '_flux' not in ifile:
        spec1dfiles.append(ifile)


nfiles = len(spec1dfiles)
objids = ['OBJ0001']*nfiles
# One telluric file
std1dfile = os.path.join(dev_suite_path, 'spec1d_cN20170331S0206-HIP62745_GNIRS_2017Mar31T083351.681.fits')
#According to Simbad HIP62745 is an A0 star with V=8.86
star_type = 'A0'
Vmag = 8.86
sensfile = os.path.join(mpia_path, 'HIP62745_sens_tell_gnirs.fits')
# telgridfile is a large grid of atmosphere models
telgridfile = os.path.join(mpia_path,'TelFit_MaunaKea_3100_26100_R20000.fits')

# TODO: set sensfile=None if you want to derive sensfunc from std1dfile
if not os.path.exists(sensfile) or clobber:
    # Run sensfun_teluric to get the sensitivity function
    TelSens = telluric.sensfunc_telluric(std1dfile, telgridfile, sensfile, star_type=star_type, star_mag=Vmag,
                                         mask_abs_lines=True, debug=True)

## Apply the sensfunc to all spectra. This is flux calibration only, but not telluric correction.
spec1dfiles_flux = [f.replace('.fits', '_flux.fits') for f in spec1dfiles]
if not os.path.exists(spec1dfiles_flux[0]):
    apply_sensfunc(spec1dfiles, sensfile, extinct_correct=False, tell_correct=False, debug=False, show=False)

# Now let's coadd the spectrum
spec1dfluxfile = os.path.join(mpia_path, 'spec1d_coadd_{:}.fits'.format(qsoname))
wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(spec1dfiles_flux, objids, show=show,
                                                                      show_exp=True,
                                                                      sensfile=sensfile, ex_value='OPT',
                                                                      outfile=spec1dfluxfile, debug=debug)

sys.exit(-1)
# This is a pickle file containing the PCA model for the QSOs
pca_file = os.path.join(mpia_path, 'qso_pca_1200_3100.pckl')

# run telluric.qso_telluric to get the final results
telloutfile = os.path.join(mpia_path, '{:}_tellmodel.fits'.format(qsoname))
outfile = os.path.join(mpia_path, 'spec1d_coadd_{:}_tellcorr.fits'.format(qsoname))

# TODO: add other modes here
TelQSO = telluric.qso_telluric(spec1dfluxfile, telgridfile, pca_file, z_qso, telloutfile, outfile,
                               create_bal_mask=None, debug=True, show=show)





