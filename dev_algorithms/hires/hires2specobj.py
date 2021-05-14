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
    thisobj = specobj.SpecObj('Echelle', 1, OBJTYPE='science', ECH_ORDERINDX=iord, ECH_ORDER=data[iord]['ORDER'])
    # Set boxcar extracted flux
    thisobj.BOX_WAVE = data[iord]['BOX_WV']
    thisobj.BOX_COUNTS = data[iord]['BOX_FX']
    thisobj.BOX_COUNTS_IVAR = utils.inverse(data[iord]['BOX_VAR'])
    thisobj.BOX_COUNTS_SIG = np.sqrt(data[iord]['BOX_VAR'])
    thisobj.BOX_MASK = data[iord]['BOX_VAR'] > 0.0
    # Do the same things for the opt flags
    thisobj.set_name()
    sobjs.add_sobj(thisobj)


spectrograph = load_spectrograph('keck_hires')
# Somehow hack the fitstbl or whatever
subheader = spectrograph.subheader_for_spec(row_fitstbl, head2d)
outfile = 'test.fits'
sobjs.write_to_fits(subheader, outfile)

