from astropy.io import fits
from pypeit.display import display
import numpy as np
from photutils import RectangularAperture, aperture_photometry
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import gaussian_filter
from photutils import DAOStarFinder

ZP = {
    'Y':28.13,
    'J':27.89,
    'H':27.54,
    'Ks':26.93,
    'J2':26.93,
    'J3':26.88,
    'H1':26.88,
    'H2':26.85
}


qso_obj_file = '/Users/joe/Downloads/MF.20200528.40601.fits'
qso_sky_file = '/Users/joe/Downloads/MF.20200528.40505.fits'
star_obj_file = '/Users/joe/Downloads/MF.20200528.40413.fits'
star_sky_file = '/Users/joe/Downloads/MF.20200528.40360.fits'

#qso_obj_file = '/home/riccardo/Downloads/MF.20200528.40601.fits'
#qso_sky_file = '/home/riccardo/Downloads/MF.20200528.40505.fits'
#star_obj_file = '/home/riccardo/Downloads/MF.20200528.40413.fits'
#star_sky_file = '/home/riccardo/Downloads/MF.20200528.40360.fits'

# Read in data
star_obj = fits.getdata(star_obj_file)
star_hdr=fits.getheader(star_obj_file)
star_sky = fits.getdata(star_sky_file)
star_diff = star_obj - star_sky
qso_obj = fits.getdata(qso_obj_file)
qso_hdr=fits.getheader(qso_obj_file)
qso_sky = fits.getdata(qso_sky_file)
qso_diff = qso_obj - qso_sky


plate_scale = 0.18 # Put in correct value
centroid=[star_hdr['CRPIX1']+15,star_hdr['CRPIX2']+20]
lengthx=round(2/plate_scale)
lengthy=round(3.5/plate_scale)

# Determine acquisition window location # Find the centroid of the image
sky_mean, sky_med, sky_sigma = sigma_clipped_stats(qso_sky[centroid[1]-lengthy:centroid[1]+lengthy, centroid[0]-lengthx:centroid[0]+lengthx],sigma_lower=5.0, sigma_upper=5.0)
acq_mask = (qso_sky > (sky_med -2.0*sky_sigma)) & (qso_sky < (sky_med + 2.0*sky_sigma))
# qso_mean, qso_med, qso_sigma = sigma_clipped_stats(qso_diff[centroid[0]-lengthx:centroid[0

#acq_mask = np.zeros_like(star_sky, dtype=bool)


fwhm = 3.0
star_smooth = gaussian_filter(star_diff, sigma=fwhm)
qso_smooth = gaussian_filter(qso_diff, sigma=fwhm)

# Find the centroid of the image
star_mean, star_med, star_sigma = sigma_clipped_stats(star_diff[centroid[0]-lengthx:centroid[0]+lengthx, centroid[1]-lengthy:centroid[1]+lengthy],
                                                      sigma_lower=5.0, sigma_upper=5.0)
qso_mean, qso_med, qso_sigma = sigma_clipped_stats(qso_diff[centroid[0]-lengthx:centroid[0]+lengthx, centroid[1]-lengthy:centroid[1]+lengthy],
                                                   sigma_lower=5.0, sigma_upper=5.0)
daofind = DAOStarFinder(fwhm=5.0, threshold=5.*star_sigma)
sources = daofind(star_diff)

# Define coordinates and apertures
aperture = RectangularAperture([sources['xcentroid'], sources['ycentroid']], lengthx*2, lengthy)
box_mask = aperture.to_mask(method='center')
back_aperture = RectangularAperture([sources['xcentroid'], sources['ycentroid']], lengthx*2, lengthy*2)

# Estimate background and aperture photometry for the star
phot_star = aperture_photometry(star_diff, aperture)
phot_back_star = aperture_photometry(star_diff, back_aperture)
counts_star=phot_star['aperture_sum']
counts_back_star=(phot_back_star['aperture_sum']-phot_star['aperture_sum'])/(back_aperture.area-aperture.area)*aperture.area
f_star=(counts_star-counts_back_star)/star_hdr['ELAPTIME']
m_star=-2.5 * np.log10(f_star) + ZP[star_hdr['FILTER']]
print(f_star,m_star)

# Estimate background and aperture photometry for the QSO
phot_qso = aperture_photometry(qso_diff, aperture)
phot_back_qso = aperture_photometry(qso_diff, back_aperture)
counts_qso=phot_qso['aperture_sum']
counts_back_qso=(phot_back_qso['aperture_sum']-phot_qso['aperture_sum'])/(back_aperture.area-aperture.area)*aperture.area
f_qso=(counts_qso-counts_back_qso)/qso_hdr['ELAPTIME']
m_qso=-2.5 * np.log10(f_qso) + ZP[qso_hdr['FILTER']]
print(f_qso,m_qso)
'''
# Display images
display.connect_to_ginga(raise_err=True, allow_new=True)
viewer, ch_star = display.show_image(star_diff,'Star', clear=True, cuts = (-5.0*star_sigma, 5.0*star_sigma))
viewer, ch_qso = display.show_image(qso_diff,'QSO', cuts = (-5.0*qso_sigma, 5.0*qso_sigma))
out_star = ch_star.cut_levels((star_mean -5.0*star_sigma, star_mean + 10.0*star_sigma))
out_qso = ch_qso.cut_levels((qso_mean -5.0*qso_sigma, qso_mean + 7.0*qso_sigma))

# After displaying all the images sync up the images with WCS_MATCH
shell = viewer.shell()
shell.start_global_plugin('WCSMatch')
shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [ch_star], {})
'''