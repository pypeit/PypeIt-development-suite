

import numpy as np
from matplotlib import pyplot as plt
from pypeit.core import telluric

tellfile = '/Users/joe/python/PypeIt/pypeit/data/telluric/atm_grids/TelFit_MaunaKea_3100_26100_R20000.fits'

tell_dict = telluric.read_telluric_grid(tellfile)


airmass_guess = 1.5
resln_guess = 2000.0

tell_guess = (np.median(tell_dict['pressure_grid']), np.median(tell_dict['temp_grid']), np.median(tell_dict['h2o_grid']),
              airmass_guess, resln_guess, 0.0, 1.0)

ind_lower =0
ind_upper = tell_dict['wave_grid'].size -3
tell_model = telluric.eval_telluric(tell_guess, tell_dict, ind_lower=ind_lower, ind_upper=ind_upper)
wave_grid = tell_dict['wave_grid'][ind_lower:ind_upper +1]

plt.plot(wave_grid,tell_model)