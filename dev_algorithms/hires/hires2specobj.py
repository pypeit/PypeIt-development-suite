import os
from matplotlib import pyplot as plt
from astropy.io import fits
from pypeit.metadata import PypeItMetaData
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
#xidl_file = '/Users/suksientie/Research/J0100+2802/Cowie/Extract/Obj_0147R.fits.gz'
head2d = fits.getheader(xidl_file, ext=1)

hdu = fits.open(xidl_file)
data = hdu[1].data

spectrograph = load_spectrograph('keck_hires_red') # apparently spectrographs.keck_hires.py exists, but name is set to keck_hires_red at the moment
detector = spectrograph.get_detector_par(None, 1)

# TODO Implement all detectors
norders =len(data)
sobjs = specobjs.SpecObjs()
for iord in range(norders-1,-1,-1):
    thisobj = specobj.SpecObj('Echelle', DET=1, OBJTYPE='science', ECH_ORDERINDX=iord, ECH_ORDER=data[iord]['ORDER'])
    # Set boxcar extracted flux
    thisobj.BOX_WAVE = data[iord]['BOX_WV']
    thisobj.BOX_COUNTS = data[iord]['BOX_FX'].astype(float) # otherwise TypeError
    thisobj.BOX_COUNTS_IVAR = utils.inverse(data[iord]['BOX_VAR'].astype(float))
    thisobj.BOX_COUNTS_SIG = np.sqrt(data[iord]['BOX_VAR'].astype(float))
    thisobj.BOX_MASK = data[iord]['BOX_VAR'] > 0.0
    # Do the same things for the opt flags
    thisobj.OPT_WAVE = data[iord]['WAVE']
    thisobj.OPT_COUNTS = data[iord]['FX'].astype(float)
    thisobj.OPT_COUNTS_IVAR = utils.inverse(data[iord]['VAR'].astype(float))
    thisobj.OPT_COUNTS_SIG = np.sqrt(data[iord]['VAR'].astype(float))
    thisobj.OPT_MASK = data[iord]['VAR'] > 0.0

    thisobj.TRACE_SPAT = data[iord]['TRACE'].astype(float)
    #thisobj.FHWM = data[iord]['SPATIAL_FWHM'] # want this?
    thisobj.NAME = data[iord]['FIELD'] # want this?
    thisobj.OPT_COUNTS_SKY = data[iord]['SKY'].astype(float) # want this?
    thisobj.ECH_FRACPOS = 0.5 # This a bogus value to ensure that we get a decent name
    thisobj.DETECTOR = detector
    thisobj.SPAT_PIXPOS = np.median(data[iord]['TRACE'][ data[iord]['VAR']> 0.0])
    thisobj.WAVE_RMS = 0.0
    thisobj.set_name()
    sobjs.add_sobj(thisobj)





# Hack in the detector container nonsense


#par = spectrograph.default_pypeit_par()
# Trying to hack fitstbl
#fitstbl = PypeItMetaData(spectrograph, par)


# Somehow hack the fitstbl or whatever
# keywords from pypeit.core.meta.define_core_meta
row_fitstbl = hdu[0].header
row_fitstbl['filename'] = None

#row_fitstbl['ra'] = None
#row_fitstbl['dec'] = None
#row_fitstbl['target'] = None
#row_fitstbl['dispname'] = None
#row_fitstbl['decker'] = None
#row_fitstbl['binning'] = None
#row_fitstbl['mjd'] = None
#row_fitstbl['airmass'] = None
#row_fitstbl['exptime'] = None


subheader = spectrograph.subheader_for_spec(row_fitstbl, head2d, allow_missing=True)
sci_path = 'Science'
outfile = os.path.join(sci_path, 'spec1d_0147R.fits')
outfile_txt = outfile.replace('.fits', '.txt')
sobjs.write_to_fits(subheader, outfile)
sobjs.write_info(outfile_txt, spectrograph.pypeline)
