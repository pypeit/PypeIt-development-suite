

from astropy.io import fits
from pypeit.display import display
import numpy as np
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import gaussian_filter



qso_obj_file = '/Users/joe/Downloads/MF.20200528.40601.fits'
qso_sky_file = '/Users/joe/Downloads/MF.20200528.40505.fits'
star_obj_file = '/Users/joe/Downloads/MF.20200528.40413.fits'
star_sky_file = '/Users/joe/Downloads/MF.20200528.40360.fits'

# Read in data
star_obj = fits.getdata(star_obj_file)
star_sky = fits.getdata(star_sky_file)
star_diff = star_obj - star_sky
qso_obj = fits.getdata(qso_obj_file)
qso_sky = fits.getdata(qso_sky_file)
qso_diff = qso_obj - qso_sky

fwhm = 3.0
star_smooth = gaussian_filter(star_diff, sigma=fwhm)
qso_smooth = gaussian_filter(qso_diff, sigma=fwhm)

display.connect_to_ginga(raise_err=True, allow_new=True)
# Display images
star_mean, star_med, star_sigma = sigma_clipped_stats(star_diff[1062:1069, 1034:1052], sigma_lower=5.0, sigma_upper=5.0)
qso_mean, qso_med, qso_sigma = sigma_clipped_stats(qso_diff[1062:1069, 1034:1052], sigma_lower=5.0, sigma_upper=5.0)
viewer, ch_star = display.show_image(star_diff,'Star', clear=True)#, cuts = (-5.0*star_sigma, 5.0*star_sigma), clear=True)
viewer, ch_qso = display.show_image(qso_diff,'QSO')#, cuts = (-5.0*qso_sigma, 5.0*qso_sigma))
out_star = ch_star.cut_levels((star_mean -5.0*star_sigma, star_mean + 10.0*star_sigma))
out_qso = ch_qso.cut_levels((qso_mean -5.0*qso_sigma, qso_mean + 7.0*qso_sigma))

# After displaying all the images sync up the images with WCS_MATCH
shell = viewer.shell()
shell.start_global_plugin('WCSMatch')
shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [ch_star], {})



