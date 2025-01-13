import os
from pypeit import specobjs

boguspath = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J0411/calwebb/pypeit/Science/'


b_filenames= ['spec1d_bogus_jw01222002001_03102_00001_slit_S200A1.fits',
             'spec1d_bogus_jw01222002001_03102_00002_slit_S200A1.fits',
             'spec1d_bogus_jw01222002001_03102_00003_slit_S200A1.fits',
             'spec1d_bogus_jw01222002001_03104_00001_slit_S200A2.fits',
             'spec1d_bogus_jw01222002001_03104_00002_slit_S200A2.fits',
             'spec1d_bogus_jw01222002001_03104_00003_slit_S200A2.fits']

bogusfiles = [os.path.join(boguspath, file) for file in b_filenames]

# Mask the bogus spectra below this wavelength, since this is covered in the original non-doctored frames with F070
#wave_bogus = 12680.0
wave_bogus = 12480.0
wave_2ndorder = 2.0*7930.0 # MODS spectrum shows Ly-beta proximity zone emission down to lam =7930.0A

for file in bogusfiles:
    sobjs = specobjs.SpecObjs.from_fitsfile(file)
    # Mask all below wave_bogus, above wave_bogus is good
    bogus_gpm_opt = (sobjs.OPT_WAVE > wave_bogus) & (sobjs.OPT_WAVE < wave_2ndorder)
    sobjs.OPT_MASK         *= bogus_gpm_opt
    sobjs.OPT_COUNTS       *= bogus_gpm_opt
    sobjs.OPT_COUNTS_IVAR  *= bogus_gpm_opt
    sobjs.OPT_COUNTS_NIVAR *= bogus_gpm_opt
    sobjs.OPT_COUNTS_SIG   *= bogus_gpm_opt
    sobjs.OPT_WAVE         *= bogus_gpm_opt

    bogus_gpm_box = (sobjs.BOX_WAVE > wave_bogus) & (sobjs.OPT_WAVE < wave_2ndorder)
    sobjs.BOX_MASK         *= bogus_gpm_box
    sobjs.BOX_COUNTS       *= bogus_gpm_box
    sobjs.BOX_COUNTS_IVAR  *= bogus_gpm_box
    sobjs.BOX_COUNTS_NIVAR *= bogus_gpm_box
    sobjs.BOX_COUNTS_SIG   *= bogus_gpm_box
    sobjs.BOX_WAVE         *= bogus_gpm_box

    outfile = file.replace('bogus', 'bogus_masked')
    sobjs.write_to_fits(sobjs.header, outfile, overwrite=True)


# For the normal reductions with F070 we need to mask above wave_bogus to not double count.


g_filenames= ['spec1d_jw01222002001_03102_00001_slit_S200A1.fits',
              'spec1d_jw01222002001_03102_00002_slit_S200A1.fits',
              'spec1d_jw01222002001_03102_00003_slit_S200A1.fits',
              'spec1d_jw01222002001_03104_00001_slit_S200A2.fits',
              'spec1d_jw01222002001_03104_00002_slit_S200A2.fits',
              'spec1d_jw01222002001_03104_00003_slit_S200A2.fits']

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
    sobjs.OPT_WAVE         *= bogus_gpm_opt

    bogus_gpm_box = (sobjs.BOX_WAVE <= wave_bogus)
    sobjs.BOX_MASK         *= bogus_gpm_box
    sobjs.BOX_COUNTS       *= bogus_gpm_box
    sobjs.BOX_COUNTS_IVAR  *= bogus_gpm_box
    sobjs.BOX_COUNTS_NIVAR *= bogus_gpm_box
    sobjs.BOX_COUNTS_SIG   *= bogus_gpm_box
    sobjs.BOX_WAVE         *= bogus_gpm_box
    outfile = file.replace('spec1d', 'spec1d_masked')
    sobjs.write_to_fits(sobjs.header, outfile, overwrite=True)


