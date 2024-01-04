import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import ZScaleInterval, ImageNormalize
from scipy import signal
from pypeit.coadd3d import DataCube

filename = '//Users/joe/kcrm_redux/J0252_testbed/j0252_SNR_reshift.fits'
cube = DataCube.from_file(filename)
flux_cube = cube.flux  # Flux datacube
error_cube = cube.sig  # Errors associated with each voxel of the flux datacube
ivar_cube = cube.ivar  # Inverse variance cube
wave = cube.wave
wcs = cube.wcs
#nspec = flux_cube.shape[0]
#wave = cube.wcs.spectral.all_pix2world(np.arange(nspec), 0)[0] * 1.0E10


z_qso = 7.00
wave_cen =  (1.0 + z_qso)*1215.67
wave_min = wave_cen*(1.0 - 500.0/3.0e5)
wave_max = wave_cen*(1.0 + 500.0/3.0e5)
wave_psf_max = wave_cen*(1.0 + 10000.0/3.0e5)
wave_mask = (wave >= wave_min) & (wave <= wave_max)
psf_mask = (wave >= wave_max) & (wave <= wave_psf_max)
wave_slice = np.argmin(np.abs(wave-wave_cen))
sub_image = np.median(flux_cube[wave_mask,:,:], axis=0)
psf_image = np.median(flux_cube[psf_mask,:,:], axis=0)
# Normalize the psf_image
# Create 2d images of the x and y coordinates
x = np.arange(psf_image.shape[0])
y = np.arange(psf_image.shape[1])
ximg, yimg = np.meshgrid(x, y)
spat_mask = (ximg >= 20.0) & (ximg <= 70.0) & (yimg >= 20.0) & (yimg <= 70.0)
psf_max = np.max(psf_image*spat_mask)
lya_max = np.max(sub_image*spat_mask)

psf_model = (psf_image/psf_max)*lya_max

diff_image = sub_image - psf_model



norm = ImageNormalize(sub_image, interval=ZScaleInterval())
fig = plt.figure(figsize=(20,20))
fig.add_subplot(111, projection=wcs, slices=('x', 'y', wave_slice))
plt.imshow(sub_image, origin='lower', cmap=plt.cm.viridis, norm=norm)
plt.xlabel('RA')
plt.ylabel('Dec')
plt.show()