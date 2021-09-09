import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from pypeit.core.coadd import multi_combspec
from IPython import embed


def ivarsmooth(flux, ivar, window):
    '''
    Boxcar smoothign of width window with ivar weights
    Args:
        flux:
        ivar:
        window:
    Returns:
    '''
    nflux = (flux.shape)[0]
    halfwindow = int(np.floor((np.round(window) - 1) / 2))
    shiftarr = np.zeros((nflux, 2 * halfwindow + 1))
    shiftivar = np.zeros((nflux, 2 * halfwindow + 1))
    shiftindex = np.zeros((nflux, 2 * halfwindow + 1))
    indexarr = np.arange(nflux)
    indnorm = np.outer(indexarr, (np.zeros(2 * halfwindow + 1) + 1))
    for i in np.arange(-halfwindow, halfwindow + 1, dtype=int):
        shiftarr[:, i + halfwindow] = np.roll(flux, i)
        shiftivar[:, i + halfwindow] = np.roll(ivar, i)
        shiftindex[:, i + halfwindow] = np.roll(indexarr, i)
    wh = (np.abs(shiftindex - indnorm) > (halfwindow + 1))
    shiftivar[wh] = 0.0
    outivar = np.sum(shiftivar, axis=1)
    nzero, = np.where(outivar > 0.0)
    zeroct = len(nzero)
    smoothflux = np.sum(shiftarr * shiftivar, axis=1)
    if (zeroct > 0):
        smoothflux[nzero] = smoothflux[nzero] / outivar[nzero]
    else:
        smoothflux = np.roll(flux, 2 * halfwindow + 1)  # kill off NAN’s
    return (smoothflux, outivar)


# spec_file = "/home/sbechtel/Documents/DEIMOS_Light_Echo/Targets/J1438A/det_all/setup_Both/Science_coadd/spec1d_DE.20190605.30172-DE.20190605.35227-J1438A.fits"
# Extens = 9, 30


# spec_file = "/home/sbechtel/Documents/DEIMOS_Light_Echo/Targets/J1438A/det_all/setup_Star/Science_coadd/spec1d_DE.20190605.30172-DE.20190605.35227-J1438A.fits"
# Extens = 10, 25

spec_file = "/home/sbechtel/Documents/DEIMOS_Light_Echo/Targets/J1630A/det_all/setup_Both/Science_coadd/spec1d_DE.20190705.25097-DE.20190705.36380-J1630A.fits"
# Extens = 6, 51

spec_hdu = fits.open(spec_file)
blue_data = spec_hdu[6].data
red_data = spec_hdu[51].data

blue_waves = blue_data['OPT_WAVE']
blue_flux = blue_data['OPT_COUNTS']
blue_ivars = blue_data['OPT_COUNTS_IVAR']

red_waves = red_data['OPT_WAVE']
red_flux = red_data['OPT_COUNTS']
red_ivars = red_data['OPT_COUNTS_IVAR']

if len(blue_waves) > len(red_waves):

    red_waves_hold = np.zeros_like(blue_waves)
    red_flux_hold = np.zeros_like(blue_waves)
    red_ivars_hold = np.zeros_like(blue_waves)

    diff = len(blue_waves) - len(red_waves)

    red_waves_hold[diff:] = red_waves
    red_flux_hold[diff:] = red_flux
    red_ivars_hold[diff:] = red_ivars

    red_waves = red_waves_hold
    red_flux = red_flux_hold
    red_ivars = red_ivars_hold

else:

    blue_waves_hold = np.zeros_like(red_waves)
    blue_flux_hold = np.zeros_like(red_waves)
    blue_ivars_hold = np.zeros_like(red_waves)

    diff = len(red_waves) - len(blue_waves)

    blue_waves_hold[diff:] = blue_waves
    blue_flux_hold[diff:] = blue_flux
    blue_ivars_hold[diff:] = blue_ivars

    blue_waves = blue_waves_hold
    blue_flux = blue_flux_hold
    blue_ivars = blue_ivars_hold

red_mask = red_waves > 10
blue_mask = blue_waves > 10

waves = np.zeros((len(blue_waves), 2))
fluxes = np.zeros((len(blue_waves), 2))
ivars = np.zeros((len(blue_waves), 2))
masks = np.zeros((len(blue_waves), 2), dtype=bool)

waves[:, 0] = blue_waves
waves[:, 1] = red_waves
fluxes[:, 0] = blue_flux
fluxes[:, 1] = red_flux
ivars[:, 0] = blue_ivars
ivars[:, 1] = red_ivars
masks[:, 0] = blue_mask
masks[:, 1] = red_mask

wgmax = np.max(red_waves)
wgmin = np.min(blue_waves[blue_waves > 10])
new_waves, new_flux, new_ivars, new_masks = multi_combspec(waves, fluxes, ivars, masks, wave_grid_max=wgmax,
                                                           wave_grid_min=wgmin)

'''plt.plot(new_waves,new_flux,'k')
plt.plot(blue_waves,blue_flux,'--')
plt.plot(red_waves,red_flux,'r--')
plt.show()'''

# Need to create some routine to match the blue/red continuum.
# Try smoothing both sides, fitting blue with linear regression (weighted by ivar), then compare red side fit to
# extrapolation from this fit?

# print(blue_waves[-300])
num = 50

blue_hold_waves = blue_waves[blue_waves > 10][-num:]
blue_smooth_flux, blue_smooth_ivar = ivarsmooth(blue_flux[blue_waves > 10][-num:],
                                                blue_ivars[blue_waves > 10][-num:], 5)
blue_smooth_flux = blue_smooth_flux.reshape((num, 1))
blue_smooth_waves = blue_hold_waves.reshape((num, 1))

blue_regr = LinearRegression()
blue_regr.fit(blue_smooth_waves, blue_smooth_flux, blue_smooth_ivar)

plt.plot(blue_waves, blue_flux, '--')
plt.plot(blue_smooth_waves[:, 0], blue_smooth_flux, 'g')
plt.plot(blue_smooth_waves[:, 0], blue_regr.predict(blue_smooth_waves), 'r')
plt.show()

red_hold_waves = red_waves[red_waves > 10][:num]
red_smooth_flux, red_smooth_ivar = ivarsmooth(red_flux[red_waves > 10][:num], red_ivars[red_waves > 10][:num], 5)
red_smooth_flux = red_smooth_flux.reshape((num, 1))
red_smooth_waves = red_hold_waves.reshape((num, 1))

red_regr = LinearRegression()
red_regr.fit(red_smooth_waves, red_smooth_flux, red_smooth_ivar)

plt.plot(blue_waves, blue_flux, 'b--', alpha=0.2)
plt.plot(red_waves, red_flux, 'r--', alpha=0.2)
plt.plot(blue_smooth_waves[:, 0], blue_regr.predict(blue_smooth_waves), 'b')
plt.plot(red_smooth_waves[:, 0], blue_regr.predict(red_smooth_waves), 'b')
plt.plot(blue_smooth_waves[:, 0], red_regr.predict(blue_smooth_waves), 'r')
plt.plot(red_smooth_waves[:, 0], red_regr.predict(red_smooth_waves), 'r')
plt.show()

# Just looked at overlap average for smoothed flux vs regression flux for star and ratios were identical to 0.0015.
# Don't need regression?

# For non-stellar spectrum, difference was large at 0.904 smooth vs 1.034 regr. Need regression?

# Would think regression is better since variance is used as weight, but by eye smooth is better. Only 7 overlap values,
# perhaps this is too small an overlap. Try another?

# Also need to try varying number of samples for regression, width of smoothing, and weights for regression.

lower_mask = blue_smooth_waves >= red_smooth_waves[0]
upper_mask = red_smooth_waves <= blue_smooth_waves[-1]

smooth_ratio = blue_smooth_flux[lower_mask].mean() / red_smooth_flux[upper_mask].mean()

blue_regr_mean = blue_regr.predict(blue_smooth_waves[lower_mask].reshape(lower_mask.sum(), 1)).mean()
red_regr_mean = red_regr.predict(red_smooth_waves[upper_mask].reshape(upper_mask.sum(), 1)).mean()

regr_ratio = blue_regr_mean / red_regr_mean

embed()