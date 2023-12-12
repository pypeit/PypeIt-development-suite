import numpy as np
import os
from pypeit.utils import inverse, fast_running_median
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy import constants as const
from matplotlib import pyplot as plt


coaddfile = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J0020/calwebb/pypeit/J0020_1dcoadd.fits'
#coaddfile = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J0411/calwebb/pypeit/J0411_1dcoadd.fits'
#coaddfile = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J1007/calwebb/pypeit/J1007_1dcoadd.fits'
#coaddfile = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J0313/calwebb/pypeit/J0313_1dcoadd.fits'
outfile = os.path.join(os.path.dirname(coaddfile), os.path.basename(coaddfile).replace('.fits', '_Flam.fits'))
#hdu = fits.open(coaddfile)
#spec_table = Table(hdu[1].data)
#coaddfile = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/02073_CLEAR_PRISM/J0252/pypeit/' \
#            'Science_coadd/spec1d_jw02073007001_03101_00001_nrs1_rate-jw02073007001_03101_00003_nrs1_rate-2073_8713.fits'

spec_table = Table.read(coaddfile, format='fits')


angle_unit = 1.0
#angle_unit = u.steradian
#angle_unit = u.steradian if 'J1007' in coaddfile else 1.0
flux_MJy = spec_table['flux']*u.MJy/angle_unit
wave = spec_table['wave_grid_mid']*u.angstrom
sigma_MJy = np.sqrt(inverse(spec_table['ivar']))*u.MJy/angle_unit
gpm = spec_table['mask'].astype(bool)

# Convert to F_lambda
# F_lambda = F_nu * c / lambda^2
#if 'J1007' in coaddfile:
#nu_to_flam_factor =const.c/np.square(wave)
#pixel_area = 4.795157053786361e-13*u.steradian
#factor = pixel_area*nu_to_flam_factor
#else:
nu_to_flam_factor =const.c/np.square(wave)
factor = nu_to_flam_factor

F_lam = (flux_MJy*factor).decompose().to(u.erg/u.s/u.cm**2/u.angstrom).value/1e-17
sigma_lam = (sigma_MJy*factor).decompose().to(u.erg/u.s/u.cm**2/u.angstrom).value/1e-17

F_lam_sm = fast_running_median(F_lam*gpm, 10)
sigma_lam_sm = fast_running_median(sigma_lam*gpm, 10)
ymin = -1.2*np.max(sigma_lam_sm)
ymax = 1.2*np.max(F_lam_sm)

plt.plot(wave, F_lam, drawstyle='steps-mid', label='JWST')
plt.ylim([ymin, ymax])
plt.ylabel(r'$F_\lambda$ (10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
plt.xlabel(r'$\lambda$ ($\AA$)')
plt.legend()
plt.show()

# Create an output table
spec_table['F_lam'] = F_lam
spec_table['sigma_lam'] = sigma_lam
# Write it out to disk
spec_table.write(outfile, overwrite=True)

