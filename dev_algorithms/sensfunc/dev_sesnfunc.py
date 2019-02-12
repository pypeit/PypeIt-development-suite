


import numpy as np
import matplotlib.pyplot as plt
import os
from pypeit.core import flux
from pypeit.core import load

iord = 8
spec1dfile = '/Users/joe/Dropbox/PypeIt_Redux/XSHOOTER/Pypeit_files/PISCO_nir_REDUCED/Science_coadd/spec1d_STD,FLUX.fits'
sobjs, head = load.load_specobjs(spec1dfile)
exptime = head['EXPTIME']

wave = sobjs[iord].optimal['WAVE_GRID']
wave_mask = sobjs[iord].optimal['WAVE_GRID'] > 0.0
counts = sobjs[iord].optimal['COUNTS']
counts_ivar = sobjs[iord].optimal['COUNTS_IVAR']

# Create copy of the arrays to avoid modification and convert to electrons / s
wave_star = wave.copy()
flux_star = counts.copy() / exptime
ivar_star = counts_ivar.copy() * exptime ** 2
std_dict = flux.get_standard_spectrum(star_type=star_type, star_mag=star_mag, ra=ra, dec=dec)


