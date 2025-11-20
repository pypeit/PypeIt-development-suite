

from astropy.table import Table
from astropy.io import fits
from pypeit.core import arc
import numpy as np
import scipy

#infile  = 'XSHOOTER_modelsky.fits'
#outfile = 'OH_XSHOOTER_lines.dat'

infile  = 'GNIRS_modelsky.fits'
outfile = 'OH_GNIRS_lines.dat'


hdu = fits.open(infile)
spec = hdu[0].data
wave = 1e4*hdu[1].data
pixvec = np.arange(spec.size)

tampl, tcent, twid, centerr, w, yprep, nsig = arc.detect_lines(spec, nfitpix=7, sigdetect = 2.0, FWHM = 12, cont_samp = 10, debug = True)
peaks_good = tcent[w]
ampl_good = tampl[w]
wave_peak = scipy.interpolate.interp1d(pixvec, wave, bounds_error=False, fill_value='extrapolate')(peaks_good)
npeak = len(wave_peak)
ion = npeak*['OH']
NIST = npeak*[1]
Instr = npeak*[32]
Source = npeak*['IDL code']

dat = Table([wave_peak, ion, NIST, Instr, ampl_good, Source], names=('wave', 'ion','NIST','Instr','amplitude','Source'))
dat.write(outfile,format='ascii.fixed_width')