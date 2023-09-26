import os
from astropy.io import fits


#rawpath = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1222/level_12/01222/'
boguspath = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1222/level_12/01222_bogus_F100LP/'

exp_list = []
for detname in ['nrs1', 'nrs2']:
    # NIRSPEC 3-point dither dither center
    scifile1 = os.path.join(boguspath, 'bogus_jw01222002001_03102_00001_' + detname + '_rate.fits')
    scifile2 = os.path.join(boguspath, 'bogus_jw01222002001_03102_00002_' + detname + '_rate.fits')
    scifile3 = os.path.join(boguspath, 'bogus_jw01222002001_03102_00003_' + detname + '_rate.fits')
    scifile4 = os.path.join(boguspath, 'bogus_jw01222002001_03104_00001_' + detname + '_rate.fits')
    scifile5 = os.path.join(boguspath, 'bogus_jw01222002001_03104_00002_' + detname + '_rate.fits')
    scifile6 = os.path.join(boguspath, 'bogus_jw01222002001_03104_00003_' + detname + '_rate.fits')
    exp_list += [scifile1, scifile2, scifile3, scifile4, scifile5, scifile6]

for file in exp_list:
    fits.setval(file, 'FILTER', value='F100LP')


