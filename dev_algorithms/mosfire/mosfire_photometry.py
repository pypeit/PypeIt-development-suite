from astropy.io import fits
from pypeit.display import display
import numpy as np
from photutils import RectangularAperture, aperture_photometry
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import gaussian_filter
from photutils import DAOStarFinder
from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph


# The following code is from MOSFIRE.CSU
def mm_to_pix(millimeters):
    demagnification = 7.24254
    return millimeters/demagnification/0.018

# following code comes from old trunk of MOSFIRE code tree:
# http://code.google.com/p/mosfire/source/browse/MOSFIRE/CSU.py?r=823640967a0497bd7bafbb0f8147228c3d3993bd

def csu_mm_to_pix(x_mm, slitno):
    """
    Convert a slit's position into a pixel value. This is a linear approximation to a sixth order polynomial fit by ccs from cooldown 8 data.
    """

    center_pix = (1042.986, 1035.879)
    numslits = 46
    rotation = np.radians(0.2282)
    tempscale = 0.99646
    mm = 1

    # _kfp is keck focal plane, the x_kfp has a "fudge factor"
    x_kfp = x_mm * tempscale - 5.0 * 1.2 * mm
    y_kfp = 5.8*mm*tempscale*(numslits - slitno + 0.70)
    # _mfp is Mosfire focal plane, convert the mms into pixels
    x_mfp = 2047 - mm_to_pix(x_kfp)
    y_mfp = mm_to_pix(y_kfp)


    # Rotate around the center
    x_mfp -= center_pix[0]
    y_mfp -= center_pix[1]

    x_mfp = np.cos(rotation)*x_mfp - np.sin(rotation)*y_mfp
    y_mfp = np.sin(rotation)*x_mfp + np.cos(rotation)*y_mfp
    x_mfp += center_pix[0]
    y_mfp += center_pix[1]

    return np.array([x_mfp, y_mfp])


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

# Object is J21.3, offset is J18.99
#qso_obj_file = '/Users/joe/Downloads/MF.20200528.40601.fits'
#qso_sky_file = '/Users/joe/Downloads/MF.20200528.40505.fits'
#star_obj_file = '/Users/joe/Downloads/MF.20200528.40413.fits'
#star_sky_file = '/Users/joe/Downloads/MF.20200528.40360.fits'

qso_obj_file = '/home/riccardo/Downloads/m201023_0093.fits'
qso_sky_file = '/home/riccardo/Downloads/m201023_0092.fits'
star_obj_file = '/home/riccardo/Downloads/m201023_0091.fits'
star_sky_file = '/home/riccardo/Downloads/m201023_0090.fits'


# Build Science image
spectrograph = load_spectrograph('keck_mosfire')
det = 1
parset = spectrograph.default_pypeit_par()
parset['scienceframe']['process']['use_illumflat'] = False
parset['scienceframe']['process']['use_pixelflat'] = False

skyImg = buildimage.buildimage_fromlist(
    spectrograph, det, parset['scienceframe'], [qso_sky_file], ignore_saturation=False)
sys.exit(-1)
# Read in data
#star_obj = fits.getdata(star_obj_file)
#star_hdr = fits.getheader(star_obj_file)
#star_sky = fits.getdata(star_sky_file)
#star_sky_hdr = fits.getheader(star_sky_file)
#star_diff = star_obj - star_sky
#qso_obj = fits.getdata(qso_obj_file)
#qso_hdr =fits.getheader(qso_obj_file)
#qso_sky = fits.getdata(qso_sky_file)
qso_sky_hdr =fits.getheader(qso_sky_file)
qso_diff = qso_obj - qso_sky

# Get bar information from header
xpix = np.zeros(92)
ypix = np.zeros(92)
for islit in range(1, 93):
    slit = int(islit + 1) / 2
    pos = qso_sky_hdr["B%0.2iPOS" % islit]
    xpix[islit-1], ypix[islit-1] = csu_mm_to_pix(pos, slit)

display.connect_to_ginga(raise_err=True, allow_new=True)
viewer, ch_sky = display.show_image(qso_diff,'QSO_SKY') #, cuts = (-5.0*qso_sigma, 5.0*qso_sigma))
display.show_points(viewer, ch_sky, ypix, xpix)

sys.exit(-1)

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