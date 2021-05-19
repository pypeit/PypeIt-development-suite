from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
from pypeit import specobj
from pypeit import specobjs
from pypeit import utils
from pypeit.spectrographs.util import load_spectrograph
# Task for you is to read in the standard star reduction from Extract directory and make some plots of the
# of the counts vs wavelength for each order. The relevant tags are wave, fx, var (sqrt var is sigma)

# Extra credit (you probably won't get this)

# Try to create a SpecObjs in echelle format that contains the data from that standard file. Look at echelle in Dev Suite forPypeit for examples. 

# Read in a Cowie xidl output for a single exposure
xidl_file = '/Users/joe/HIRES_redux/J0100+2802/Cowie/Extract/Obj_0147R.fits.gz'
head2d = fits.getheader(xidl_file, ext=1)

hdu = fits.open(xidl_file)
data = hdu[1].data

# TODO Implement all detectors
norders =len(data)
sobjs = specobjs.SpecObjs()
for iord in range(norders-1,-1,-1):
    thisobj = specobj.SpecObj('Echelle', DET=1, OBJTYPE='science', ECH_ORDERINDX=iord, ECH_ORDER=data[iord]['ORDER'])
    # Set boxcar extracted flux
    thisobj.BOX_WAVE = data[iord]['BOX_WV']
    thisobj.BOX_COUNTS = data[iord]['BOX_FX']
    thisobj.BOX_COUNTS_IVAR = utils.inverse(data[iord]['BOX_VAR'])
    thisobj.BOX_COUNTS_SIG = np.sqrt(data[iord]['BOX_VAR'])
    thisobj.BOX_MASK = data[iord]['BOX_VAR'] > 0.0
    # Do the same things for the opt flags
    thisobj.OPT_WAVE = data[iord]['WAVE']
    thisobj.OPT_COUNTS = data[iord]['FX']
    thisobj.OPT_COUNTS_IVAR = utils.inverse(data[iord]['VAR'])
    thisobj.OPT_COUNTS_SIG = np.sqrt(data[iord]['VAR'])
    thisobj.TRACE_SPAT = data[iord]['TRACE']
    thisobj.FHWM = data[iord]['SPATIAL_FWHM'] # want this?
    thisobj.NAME = data[iord]['FIELD'] # want this?
    thisobj.OPT_COUNTS_SKY = data[iord]['SKY'] # want this?

    thisobj.set_name()
    sobjs.add_sobj(thisobj)


spectrograph = load_spectrograph('keck_hires')
# Somehow hack the fitstbl or whatever
# keywords from pypeit.core.meta.define_core_meta
row_fitstbl = hdu[0].header
row_fitstbl['ra'] = None
row_fitstbl['dec'] = None
row_fitstbl['target'] = None
row_fitstbl['dispname'] = None
row_fitstbl['decker'] = None
row_fitstbl['binning'] = None
row_fitstbl['mjd'] = None
row_fitstbl['airmass'] = None
row_fitstbl['exptime'] = None

subheader = spectrograph.subheader_for_spec(row_fitstbl, head2d)
outfile = 'test.fits'
sobjs.write_to_fits(subheader, outfile)

