

import os
import numpy as np
# data models
from jwst import datamodels
from IPython import embed

os.environ['CRDS_PATH'] = '/Users/joe/crds_cache/jwst_pub'
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds-pub.stsci.edu'
from matplotlib import pyplot as plt
from astropy.io import fits

from pypeit.display import display
from jwst.pipeline import Spec2Pipeline

rawpath_level2 = '/Users/joe/jwst_redux/NIRSPEC/Raw/02736_ERO_SMACS0723_G395MG235M/level_2'

output_dir = '/Users/joe/jwst_redux/NIRSPEC/redux/calwebb/output'
asn_file = '/Users/joe/jwst_redux/NIRSPEC/redux/calwebb/Raw/jw02736-o007_20220712t161855_spec2_006_asn.json'

runflag = False
if runflag:
    spec2 = Spec2Pipeline()
    spec2.save_results = True
    spec2.output_dir = output_dir
    result = spec2(asn_file)

# take a look at the results - open the level 2b files
calfile = os.path.join(output_dir, 'jw02736007001_03103_00001_nrs1_cal.fits')
s2dfile = calfile.replace('_cal.fits', '_s2d.fits')
x1dfile = calfile.replace('_cal.fits', '_x1d.fits')
cal = datamodels.open(calfile)  # this contains the calibrated unrectified 2D spectra
s2d = datamodels.open(s2dfile)  # this contains the calibrated *rectified* 2D spectra
x1d = datamodels.open(x1dfile)  # # this contains the aperture-extracted 1D spectra
for i, slit in enumerate(cal.slits):
    print(slit.name)
    calsci = slit.data  # contains the pixel data from the cal file (SCI extension)
    s2dsci = s2d.slits[i].data  # contains the pixel data from the s2d file

    # determine the wavelength scale of the s2d data for plotting purposes
    # get the data model WCS object
    wcsobj = s2d.slits[i].meta.wcs
    y, x = np.mgrid[:s2dsci.shape[0], : s2dsci.shape[1]]  # grid of pixel x,y indices
    det2sky = wcsobj.get_transform('detector','world')  # the coordinate transform from detector space (pixels) to sky (RA, DEC in degrees)
    ra, dec, s2dwave = det2sky(x, y)  # RA, Dec, wavelength (microns) for each pixel
    s2dwaves = s2dwave[0, :]  # only need a single row of values since this is the rectified spectrum
    xtint = np.arange(100, s2dsci.shape[1], 100)
    xtlab = np.round(s2dwaves[xtint], 2)  # wavelength labels for the x-axis

    # get wavelength & flux from the x1d data model
    x1dwave = x1d.spec[i].spec_table.WAVELENGTH
    x1dflux = x1d.spec[i].spec_table.FLUX

    # plot the unrectified calibrated 2D spectrum
    display.show_image(calsci, chname='unrectified')

    #show_image(calsci, -0.01, 0.01, aspect=5., scale='linear', units='MJy')

    # plot the rectified 2D spectrum
    display.show_image(s2dsci, chname='rectified')
    #show_image(s2dsci, -0.01, 0.01, aspect=5., scale='linear', units='MJy')
    #plt.xticks(xtint, xtlab)
    #plt.xlabel('wavelength (microns)')

    embed()
    # plot the 1D extracted spectrum
    fig = plt.figure(figsize=(19, 8))
    plt.plot(x1dwave, x1dflux)
    plt.xlabel('wavelength (microns)')
    plt.ylabel('flux (Jy)')
    plt.show()

