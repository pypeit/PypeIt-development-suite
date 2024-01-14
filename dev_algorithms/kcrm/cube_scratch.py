import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import ZScaleInterval, ImageNormalize
from pypeit.coadd3d import DataCube

filename = '/Users/joe/kcrm_redux/2023sep23/redux/J0252_test_datacube.fits'
cube = DataCube.from_file(filename)
flux_cube = cube.flux  # Flux datacube
error_cube = cube.sig  # Errors associated with each voxel of the flux datacube
ivar_cube = cube.ivar  # Inverse variance cube
wcs = cube.wcs
nspec = flux_cube.shape[0]
wave = cube.wcs.spectral.all_pix2world(np.arange(nspec), 0)[0] * 1.0E10

wave_cen = 8.008*1215.67
wave_min = wave_cen*(1.0 - 500.0/3.0e5)
wave_max = wave_cen*(1.0 + 500.0/3.0e5)
wave_mask = (wave >= wave_min) & (wave <= wave_max)
wave_slice = np.argmin(np.abs(wave-wave_cen))
sub_image = np.median(flux_cube[wave_mask,:,:], axis=0)
norm = ImageNormalize(sub_image, interval=ZScaleInterval())
fig = plt.figure(figsize=(20,20))
fig.add_subplot(111, projection=wcs, slices=('x', 'y', wave_slice))
plt.imshow(sub_image, origin='lower', cmap=plt.cm.viridis, norm=norm)
plt.xlabel('RA')
plt.ylabel('Dec')
plt.show()


