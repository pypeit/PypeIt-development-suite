import os
import numpy as np
from pypeit import specobjs

boguspath = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J0020/calwebb/pypeit/Science/'


b_filenames= ['spec1d_bogus_jw01222012001_03102_00001_slit_S200A1.fits',
              'spec1d_bogus_jw01222012001_03102_00002_slit_S200A1.fits',
              'spec1d_bogus_jw01222012001_03102_00003_slit_S200A1.fits',
              'spec1d_bogus_jw01222012001_03104_00001_slit_S200A2.fits',
              'spec1d_bogus_jw01222012001_03104_00002_slit_S200A2.fits',
              'spec1d_bogus_jw01222012001_03104_00003_slit_S200A2.fits']

bogusfiles = [os.path.join(boguspath, file) for file in b_filenames]

# Mask the bogus spectra below this wavelength, since this is covered in the original non-doctored frames with F070
wave_bogus = 12699.12 # This is the red cutoff of the original F070/F140H frames, cutoff by calwebb wavelength solution

# Second order masking
# --------------------
#Technically we need to worry about second order contamination up to the red cutoff of the A2 slit, which is 18630,
#or wavelengths up to 9315A in first order. So wavelengths in the quasar from 7000A-9315.0 can contaminate from 14000-18630A.

#X-shooter: The X-shooter spectrum is completely absorbed in the 7000-9315A region. The only possible feaature is a
# tentative 3sigma transmission spike around 8189A (Ly-a at z = 5.7) in the X-shooter spectrum.

#JWST 140H: The JWST spectrum is completely absorbed in the 8100-9315A region. There appears to be a similar feature at
# 8188A in the JWST spectrum in all three exposures, but loooking at the 2d it does not look very convincing
# (cosmics in all three?). There is no obvious evidence for second order contamination at 2*8188=16376A in the JWST spectrum though.

#Masking strategy: Given that we see a transmission spike like feature in the X-shooter and all three JWST exposures we
# conservatively mask the transmissinon spike from 8186-8190 which as second order would be 16372-16380A.

wave_1st_order_bad_min, wave_1st_order_bad_max = 8186.0, 8190.0
wave_2nd_order_bad_min, wave_2nd_order_bad_max = 2.0*wave_1st_order_bad_min, 2.0*wave_1st_order_bad_max

for file in bogusfiles:
    sobjs = specobjs.SpecObjs.from_fitsfile(file)
    # Mask the transmission spike contaminated region.
    bad_mask_opt = (sobjs.OPT_WAVE >= wave_2nd_order_bad_min) & (sobjs.OPT_WAVE <= wave_2nd_order_bad_max)
    # Mask all below wave_bogus since these wavelengths are covered in the original F070/140H frames
    # Good wavelengths are > wave_bogus and not in the transmission spike contaminated region
    bogus_gpm_opt = (sobjs.OPT_WAVE > wave_bogus) & np.logical_not(bad_mask_opt)
    sobjs.OPT_MASK         *= bogus_gpm_opt
    sobjs.OPT_COUNTS       *= bogus_gpm_opt
    sobjs.OPT_COUNTS_IVAR  *= bogus_gpm_opt
    sobjs.OPT_COUNTS_NIVAR *= bogus_gpm_opt
    sobjs.OPT_COUNTS_SIG   *= bogus_gpm_opt
    #sobjs.OPT_WAVE         *= bogus_gpm_opt

    # Mask the transmission spike contaminated region.
    bad_mask_box = (sobjs.BOX_WAVE >= wave_2nd_order_bad_min) & (sobjs.BOX_WAVE <= wave_2nd_order_bad_max)
    # Mask all below wave_bogus since these wavelengths are covered in the original F070/140H frames
    # Good wavelengths are > wave_bogus and not in the transmission spike contaminated region
    bogus_gpm_box = (sobjs.BOX_WAVE > wave_bogus) & np.logical_not(bad_mask_box)
    sobjs.BOX_MASK         *= bogus_gpm_box
    sobjs.BOX_COUNTS       *= bogus_gpm_box
    sobjs.BOX_COUNTS_IVAR  *= bogus_gpm_box
    sobjs.BOX_COUNTS_NIVAR *= bogus_gpm_box
    sobjs.BOX_COUNTS_SIG   *= bogus_gpm_box
    #sobjs.BOX_WAVE         *= bogus_gpm_box

    outfile = file.replace('bogus', 'bogus_masked')
    sobjs.write_to_fits(sobjs.header, outfile, overwrite=True)


# For the normal reductions with F070 we need to mask above wave_bogus to not double count.

g_filenames= ['spec1d_jw01222012001_03102_00001_slit_S200A1.fits',
              'spec1d_jw01222012001_03102_00002_slit_S200A1.fits',
              'spec1d_jw01222012001_03102_00003_slit_S200A1.fits',
              'spec1d_jw01222012001_03104_00001_slit_S200A2.fits',
              'spec1d_jw01222012001_03104_00002_slit_S200A2.fits',
              'spec1d_jw01222012001_03104_00003_slit_S200A2.fits']

goodfiles = [os.path.join(boguspath, file) for file in g_filenames]

for file in goodfiles:
    sobjs = specobjs.SpecObjs.from_fitsfile(file)
    # Mask all below wave_bogus, above wave_bogus is good
    bogus_gpm_opt = (sobjs.OPT_WAVE <= wave_bogus)
    sobjs.OPT_MASK         *= bogus_gpm_opt
    sobjs.OPT_COUNTS       *= bogus_gpm_opt
    sobjs.OPT_COUNTS_IVAR  *= bogus_gpm_opt
    sobjs.OPT_COUNTS_NIVAR *= bogus_gpm_opt
    sobjs.OPT_COUNTS_SIG   *= bogus_gpm_opt
#    sobjs.OPT_WAVE         *= bogus_gpm_opt

    bogus_gpm_box = (sobjs.BOX_WAVE <= wave_bogus)
    sobjs.BOX_MASK         *= bogus_gpm_box
    sobjs.BOX_COUNTS       *= bogus_gpm_box
    sobjs.BOX_COUNTS_IVAR  *= bogus_gpm_box
    sobjs.BOX_COUNTS_NIVAR *= bogus_gpm_box
    sobjs.BOX_COUNTS_SIG   *= bogus_gpm_box
#    sobjs.BOX_WAVE         *= bogus_gpm_box
    outfile = file.replace('spec1d', 'spec1d_masked')
    sobjs.write_to_fits(sobjs.header, outfile, overwrite=True)


