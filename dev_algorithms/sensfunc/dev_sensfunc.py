


import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import astropy.units as u
from astropy.io import fits
from pypeit.core import flux
from pypeit.core import load
from pypeit import utils
PYPEIT_FLUX_SCALE = 1e-17

iord = 7
spec1dfile = '/Users/fdavies/Dropbox/PypeIt_Redux/XSHOOTER/Pypeit_files/PISCO_nir_REDUCED/Science_coadd/spec1d_STD,FLUX.fits'
sobjs, head = load.load_specobjs(spec1dfile)
exptime = head['EXPTIME']
airmass = head['AIRMASS']

wave_mask = sobjs[iord].optimal['WAVE'] > 10000.0*u.Angstrom
wave = sobjs[iord].optimal['WAVE'][wave_mask]
counts = sobjs[iord].optimal['COUNTS'][wave_mask]
counts_ivar = sobjs[iord].optimal['COUNTS_IVAR'][wave_mask]

# Create copy of the arrays to avoid modification and convert to electrons / s
wave_star = wave.copy()
counts_ps = counts.copy()/exptime
counts_ps_ivar = counts_ivar.copy() * exptime ** 2

xshooter_file = '/Users/fdavies/PypeIt-development-suite/dev_algorithms/sensfunc/xshooter_standards/fEG274.dat'
output = np.loadtxt(xshooter_file)
wave_true = output[:,0]
flux_true = output[:,1]/PYPEIT_FLUX_SCALE
# Create a fake std_dict for EG274
std_dict = {}
std_dict['std_ra'] = head['RA']
std_dict['std_dec'] = head['DEC']
std_dict['exptime'] = exptime
std_dict['airmass'] = head['AIRMASS']
std_dict['std_name'] = ['EG274']
std_dict['cal_file'] = ['EG274']
std_dict['wave'] = wave_true
std_dict['flux'] = flux_true

def sensfunc(theta, wave_star, counts_ps, counts_ps_ivar, wave_min, wave_max, tellmodel, flux_true):

    sensmodel = utils.func_val(theta, wave_star, 'legendre', minx=wave_min, maxx=wave_max)
    if np.sum(sensmodel) < 1e-6:
        return np.inf
    else:
        chi_vec = (sensmodel != 0.0)*(tellmodel*flux_true/(sensmodel + (sensmodel == 0.0)) - counts_ps)*counts_ps_ivar
        chi2 = np.sum(np.square(chi_vec))
    return chi2

#std_dict = flux.get_standard_spectrum(star_type=star_type, star_mag=star_mag, ra=ra, dec=dec)

# Interpolate standard star spectrum onto the data wavelength grid
flux_star = scipy.interpolate.interp1d(std_dict['wave'], std_dict['flux'], bounds_error=False,fill_value='extrapolate')(wave_star)

# Load in the telluric grid
tel_wave, tel_model, pg, tg, hg, ag = read_telluric_grid('TelFit_Paranal_NIR_AM1.03_R5500.fits')

seed = np.fmin(int(np.abs(np.sum(counts_ps[np.isfinite(counts_ps)]))), 2 ** 32 - 1)
random_state = np.random.RandomState(seed=seed)
bounds = [(shift_cc + nspec * shift_mnmx[0], shift_cc + nspec * shift_mnmx[1]), stretch_mnmx]
# TODO Can we make the differential evolution run faster?
result = scipy.optimize.differential_evolution(sensfunc, args=(y1, y2), tol=1e-4,bounds=bounds, disp=False, polish=True, seed=random_state)
