

from pypeit import ginga
from astropy.io import fits

lowrdx_path = '/Users/joe/pypeit_vs_lowredux/c17_60L/Red400/'
pypeit_path = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Keck_LRIS_red/multi_400_8500_d560/MF_keck_lris_red/'

flatfile_lowrdx = lowrdx_path + 'pixflat-r170320_2057.fits'
flatfile_pypeit = pypeit_path + 'MasterFlatField_A_01_aa.fits'

hdu = fits.open(flatfile_lowrdx)
flat_lowrdx = hdu[0].data
flat_lowrdx = flat_lowrdx[:,0:1024]

hdu = fits.open(flatfile_pypeit)
flat_pypeit = hdu[0].data

cuts = (0.8,1.5)
viewer, ch = ginga.show_image(flat_lowrdx, chname='LOWREDUX',cuts = cuts)
viewer, ch = ginga.show_image(flat_pypeit, chname='PYPEIT',cuts = cuts)