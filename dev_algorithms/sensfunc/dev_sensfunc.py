


import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import astropy.units as u
from astropy.io import fits
from pypeit.core import flux
from pypeit.core import load
from pypeit.core import coadd2d
from pypeit import utils
PYPEIT_FLUX_SCALE = 1e-17
from astropy.io import fits



def read_telluric_grid(filename):
    hdul = fits.open(filename)
    wave_grid = hdul[1].data
    model_grid = hdul[0].data

    pg = hdul[0].header['PRES0']+hdul[0].header['DPRES']*np.arange(0,hdul[0].header['NPRES'])
    tg = hdul[0].header['TEMP0']+hdul[0].header['DTEMP']*np.arange(0,hdul[0].header['NTEMP'])
    hg = hdul[0].header['HUM0']+hdul[0].header['DHUM']*np.arange(0,hdul[0].header['NHUM'])
    if hdul[0].header['NAM'] > 1:
        ag = hdul[0].header['AM0']+hdul[0].header['DAM']*np.arange(0,hdul[0].header['NAM'])
    else:
        ag = hdul[0].header['AM0']+1*np.arange(0,1)

    return 10.0*wave_grid, model_grid, pg, tg, hg, ag


def interp_telluric_grid(theta,pg,tg,hg,ag,model_grid):

    press,temp,hum,airmass = theta
    if len(pg) > 1:
        p_ind = int(np.round((press-pg[0])/(pg[1]-pg[0])))
    else:
        p_ind = 0
    if len(tg) > 1:
        t_ind = int(np.round((temp-tg[0])/(tg[1]-tg[0])))
    else:
        t_ind = 0
    if len(hg) > 1:
        h_ind = int(np.round((hum-hg[0])/(hg[1]-hg[0])))
    else:
        h_ind = 0
    if len(ag) > 1:
        a_ind = int(np.round((airmass-ag[0])/(ag[1]-ag[0])))
    else:
        a_ind = 0

    return model_grid[p_ind,t_ind,h_ind,a_ind]


def sensfunc(theta, wave_star, counts_ps, counts_ps_ivar, wave_min, wave_max, flux_true, order, tell_dict):

    theta_sens = theta[:order+1]
    theta_tell = theta[order+1:]
    sensmodel = utils.func_val(theta_sens, wave_star, 'legendre', minx=wave_min, maxx=wave_max)
    tellmodel = interp_telluric_grid(theta_tell,tell_dict['pg'],tell_dict['tg'],tell_dict['hg'],
                                     tell_dict['ag'], tell_dict['tell_model'])
    if np.sum(sensmodel) < 1e-6:
        return np.inf
    else:
        chi_vec = (sensmodel != 0.0)*(tellmodel*flux_true/(sensmodel + (sensmodel == 0.0)) - counts_ps)*np.sqrt(counts_ps_ivar)
        chi2 = np.sum(np.square(chi_vec))
    return chi2


iord = 7
spec1dfile = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/XSHOOTER/Pypeit_files/PISCO_nir_REDUCED/Science_coadd/spec1d_STD,FLUX.fits')
sobjs, head = load.load_specobjs(spec1dfile)
exptime = head['EXPTIME']
airmass = head['AIRMASS']

wave_mask = sobjs[iord].optimal['WAVE_GRID'] > 0.0
wave = sobjs[iord].optimal['WAVE_GRID'][wave_mask]
counts = sobjs[iord].optimal['COUNTS'][wave_mask]
counts_ivar = sobjs[iord].optimal['COUNTS_IVAR'][wave_mask]

# Create copy of the arrays to avoid modification and convert to electrons / s
wave_star = wave.copy()
counts_ps = counts.copy()/exptime
counts_ps_ivar = counts_ivar.copy() * exptime ** 2
counts_ps_mask = sobjs[iord].optimal['MASK'][wave_mask]

dev_path = os.getenv('PYPEIT_DEV')
xshooter_file = os.path.join(dev_path,'dev_algorithms/sensfunc/xshooter_standards/fEG274.dat')
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

#std_dict = flux.get_standard_spectrum(star_type=star_type, star_mag=star_mag, ra=ra, dec=dec)

# Interpolate standard star spectrum onto the data wavelength grid
flux_star = scipy.interpolate.interp1d(std_dict['wave'], std_dict['flux'], bounds_error=False,fill_value='extrapolate')(wave_star)

# Load in the telluric grid
telgridfile = os.path.join(dev_path,'dev_algorithms/sensfunc/TelFit_Paranal_NIR_AM1.03_R7000.fits')
tell_wave_grid, tell_model_grid, pg, tg, hg, ag = read_telluric_grid(telgridfile)
ind_lower, ind_upper = coadd2d.get_wave_ind(tell_wave_grid, np.min(wave_star), np.max(wave_star))
tell_wave_grid = tell_wave_grid[ind_lower:ind_upper]
tell_model_grid = tell_model_grid[:,:,:,:,ind_lower:ind_upper]

tell_guess = (750.0,0.0,50.0,airmass)
tell_model1 = interp_telluric_grid(tell_guess,pg,tg,hg,ag,tell_model_grid)
tell_dict = dict(pg=pg,tg=tg,hg=hg,ag=ag,tell_model=tell_model_grid)

sensguess = tell_model1*flux_star/(counts_ps + (counts_ps < 0.0))
inmask = counts_ps_mask & (counts_ps > 0.0) &  np.isfinite(sensguess) & (counts_ps_ivar > 0.0)
order = 5
func = 'legendre'
wave_min = wave_star.min()
wave_max = wave_star.max()

mask, coeff = utils.robust_polyfit_djs(wave_star, sensguess, order, function = func, minx = wave_min, maxx = wave_max,
                                 inmask=inmask, lower=3.0, upper=3.0,use_mad=True)
sensfit_guess = utils.func_val(coeff, wave_star, func, minx=wave_min, maxx=wave_max)
delta_coeff = (0.7, 1.3)

seed = np.fmin(int(np.abs(np.sum(counts_ps[np.isfinite(counts_ps)]))), 2 ** 32 - 1)
random_state = np.random.RandomState(seed=seed)
bounds_poly = [(this_coeff*delta_coeff[0], this_coeff*delta_coeff[1]) for this_coeff in coeff]
bounds_tell = [(tell_dict['pg'].min(), tell_dict['pg'].max()),
               (tell_dict['tg'].min(), tell_dict['tg'].max()),
               (tell_dict['hg'].min(), tell_dict['hg'].max()),
               (tell_dict['ag'].min(), tell_dict['ag'].max())]
bounds_poly.extend(bounds_tell)
# TODO Can we make the differential evolution run faster?
result = scipy.optimize.differential_evolution(sensfunc, args=(wave_star, counts_ps, counts_ps_ivar, wave_min, wave_max,
                                                               flux_star, order, tell_dict), tol=1e-4,
                                                               bounds=bounds_poly,
                                                               popsize=40,recombination=0.6,
                                                               disp=True, polish=True, seed=random_state)

coeff_out = result.x[:order+1]
tell_out = result.x[order+1:]

tellfit = interp_telluric_grid(tell_out,pg,tg,hg,ag,tell_model_grid)
sensfit = utils.func_val(coeff_out, wave_star, func, minx=wave_min, maxx=wave_max)

plt.plot(wave_star,counts_ps*sensfit)
plt.plot(wave_star,counts_ps*sensfit/(tellfit + (tellfit == 0.0)))
plt.plot(wave_star,flux_star)
plt.ylim(-0.1*flux_star.max(),1.5*flux_star.max())
plt.show()
