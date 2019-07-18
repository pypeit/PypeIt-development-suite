
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
mpia_outpath = os.path.join(mpia_path,'output')
dev_suite_path = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT/Gemini_GNIRS/GNIRS/Science')

# List of science files for pisco
spec1dfiles = glob.glob(dev_suite_path + '/spec1d*pisco*.fits')
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
spec1dfluxfile = 'spec1d_coadd_{:}.fits'.format(qsoname)
wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(spec1dfiles_flux, objids, show=show,
                                                                      sensfile=sensfile, ex_value='OPT',
                                                                      outfile=spec1dfluxfile, debug=debug)

# This is a pickle file containing the PCA model for the QSOs
pca_file = os.path.join(mpia_path, 'qso_pca_1200_3100.pckl')

# run telluric.qso_telluric to get the final results
telloutfile = '{:}_tellmodel.fits'.format(qsoname)
outfile = 'spec1d_coadd_{:}_tellcorr.fits'.format(qsoname)

# TODO: add other modes here
TelQSO = telluric.qso_telluric(spec1dfluxfile, telgridfile, pca_file, z_qso, telloutfile, outfile,
                               create_bal_mask=None, debug=True, show=show)





import os
import numpy as np

import telluric
from pypeit.core.flux_calib import apply_sensfunc
from pypeit.core import coadd1d
from pypeit import msgs

debug = False
show = True
do_sens = False

z_qso = 7.50
npca = 8
ex_value = 'OPT'
qsoname = 'BlueHawaii_GNIRS'

datapath = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/GNIRS/BlueHawaii/')

# TODO: change the spec1dlist to the pypeit format and change the reader accordingly
spec1dlist = 'spec1dlist'
spec1dfiles = np.genfromtxt(os.path.join(datapath, spec1dlist),dtype='str')
nfiles = len(spec1dfiles)
fnames = []
for ifile in range(nfiles):
    fnames.append(os.path.join(datapath,spec1dfiles[ifile]))

#TODO: the objids shoul be read in from the pypeit format file as noted above.
objids = ['OBJ0001']*nfiles

std1dfile = os.path.join(os.getenv('HOME'),
                          'Dropbox/PypeIt_Redux/GNIRS/BlueHawaii/spec1d_N20190531S0067-HIP56736_GNIRS_2019May31T073323.707.fits')

# get the pca pickle file and atmosphere model grid
pca_file = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/qso_pca_1200_3100.pckl')
telgridfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/TelFit_MaunaKea_3100_26100_R20000.fits')

# TODO: set sensfile=None if you want to derive sensfunc from std1dfile
sensfile = os.path.join(os.getenv('HOME'), 'Dropbox/PypeIt_Redux/GNIRS/HIP56736_sens_tell_gnirs.fits')
if do_sens:
    if std1dfile is None:
        msgs.error('You need either give a std1dfile to derive sensfunc')
    else:
        # run telluric.sensfunc_telluric to get the sensfile
        TelSens = telluric.sensfunc_telluric(std1dfile, telgridfile, sensfile, star_type='A0', star_mag=8.78,
                                             star_ra=None, star_dec=None, mask_abs_lines=True, debug=True)

## Apply the sensfunc to all spectra (only sensfunc but not tellluric)
# TODO: change show=False to show=show
apply_sensfunc(fnames, sensfile, extinct_correct=False, tell_correct=False, debug=debug, show=False)

fnames_flux = [f.replace('.fits', '_flux.fits') for f in fnames]

## Let's coadd all the fluxed spectra
# you should get a coadded spectrum named as 'spec1d_stack_{:}.fits'.format(qsoname)
#                a straight merge of individual order stacked spectra named as 'spec1d_merge_{:}.fits'.format(qsoname)
#                a individual order stacked spectra (multi-extension) named as 'spec1d_order_{:}.fits'.format(qsoname)
# TODO: change the outfile to work with datapath. It's a hard coding on these names in coadd1d
wave_stack, flux_stack, ivar_stack, mask_stack = coadd1d.ech_combspec(fnames_flux, objids, show=show, sensfile=sensfile,
                                                                      ex_value='OPT', outfile=qsoname, debug=debug)
# run telluric.qso_telluric to get the final results
spec1dfluxfile = 'spec1d_stack_{:}.fits'.format(qsoname)
telloutfile = 'spec1d_stack_{:}_tellmodel.fits'.format(qsoname)
outfile = 'spec1d_stack_{:}_tellcorr.fits'.format(qsoname)

# TODO: add other modes here
TelQSO = telluric.qso_telluric(spec1dfluxfile, telgridfile, pca_file, z_qso, telloutfile, outfile,
                               create_bal_mask=None, debug=True, show=show)


