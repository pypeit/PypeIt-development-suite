import os
from matplotlib import pyplot as plt
from astropy.io import fits

from pypeit.display import display


path = '/Users/joe/jwst_redux/NIRSPEC_ER/Raw/02736_ERO_SMACS0723_G395/level_2'

# NIRSPEC 3-point dither
file1 = os.path.join(path, 'jw02736007001_03103_00001_nrs1_rate.fits')
file2 = os.path.join(path, 'jw02736007001_03103_00002_nrs1_rate.fits')
file3 = os.path.join(path, 'jw02736007001_03103_00003_nrs1_rate.fits')
hdu1 = fits.open(file1)

image1, err1 = hdu1[1].data, hdu1[2].data
hdu2 = fits.open(file2)
image2, err2 = hdu2[1].data, hdu2[2].data
hdu3 = fits.open(file3)
image3, err3 = hdu3[1].data, hdu3[2].data

diff = image1 - (image2 + image3)/2.0

display.show_image(diff, chname='diff2d')