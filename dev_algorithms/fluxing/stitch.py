
from matplotlib import pyplot as plt
from pypeit import sensfunc



# Load in a sensfunc
sensfile = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/keck_mosfire/Y_long/sens_m191118_0064-GD71_MOSFIRE_2019Nov18T104704.507.fits'
wave, zeropoint, meta_table, out_table, header_sens = sensfunc.SensFunc.load(sensfile)
plt.plot(wave[:,0], zeropoint[:, 0])
plt.show()

# Read in multiple sensfuncs for the same grating (can be different stars) but different dispangle (central wavelength).
# Overpolot them, then we will code up the stitching together.

