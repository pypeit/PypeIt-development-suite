
import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table


# Get this from the esi structure
cbin = 1
echfile = '/Users/joe/Dropbox/Cowie_2002-02-17/FSpec/SDSS1030+0524a_ech.fits.gz'

spec2d_in = Table.read(echfile)
ordrs = np.array([0, 9],dtype=int)
nordrs = ordrs[1] - ordrs[0] + 1

bad = np.zeros((nordrs,nordrs, 2))
if cbin == 1.0:
    bad[0, 0, :] = [4330, 4400]  # UPPER EDGE
    bad[1, 0, :] = [4190, 4230]  # LOWER EDGE
    bad[1, 1, :] = [4645, 4685]  # HOTSPOT
    bad[2, 0, :] = [4440, 4500]  # LOWER EDGE
    bad[2, 1, :] = [5028, 5065]  # UPPER EDGE= + BAD COLUMN
    bad[3, 0, :] = [4675, 4860]  # LOWER EDGE
    bad[4, 0, :] = [5050, 5250]  # LOWER EDGE
    bad[5, 0, :] = [5600, 5710]  # LOWER EDGE
    bad[6, 0, :] = [6235, 6280]  # LOWER EDGE
elif cbin == 2:
    bad[0, 0, :] = [4325, 4347]  # UPPER EDGE
    bad[1, 0, :] = [4190, 4230]  # LOWER EDGE
    bad[1, 1, :] = [4651, 4676]  # HOTSPOT
    bad[2, 0, :] = [4440, 4480]  # LOWER EDGE
    bad[2, 1, :] = [5044, 5060]  # UPPER EDGE= + BAD COLUMN
    bad[3, 0, :] = [4675, 4831]  # LOWER EDGE
    bad[4, 0, :] = [5000, 5220]  # LOWER EDGE
    bad[5, 0, :] = [5600, 5700]  # LOWER EDGE
    bad[6, 0, :] = [6225, 6270.0]  # LOWER EDGE

for qq in np.arange(ordrs[0],ordrs[1]+1):
    ithis = spec2d['var'][:,qq] > 0.0
    wave_min = spec2d['wave'][ithis,qq].min()
    wave_max = spec2d['wave'][ithis,qq].max()
    fitmask = (spec2d['wave'][:,qq] >= wave_min) & (spec2d['wave'][:,qq] <= wave_max)
    nfit = np.sum(fitmask)
    loglam = np.log10(spec2d['wave'][fitmask,qq])



