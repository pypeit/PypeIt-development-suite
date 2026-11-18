import os, sys
import numpy as np
from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.edgetrace import EdgeTraceSet
mode=sys.argv[1]
instrument = load_spectrograph('keck_deimos')
par = instrument.default_pypeit_par()
files=[os.path.join(os.getenv('PYPEIT_DEV'),'RAW_DATA','keck_deimos','1200G_Cooper',f) for f in ['d0115_0023.fits.gz','d0115_0024.fits.gz','d0115_0025.fits.gz']]
tI = buildimage.buildimage_fromlist(instrument,1,par['calibrations']['traceframe'],files)
par = instrument.config_specific_par(tI.files[0])
e = EdgeTraceSet(tI, instrument, par['calibrations']['slitedges'], auto=True)
s=e.get_slits()
np.savez(f'mdm4_{mode}.npz', left=s.left_init[s.nspec//2], right=s.right_init[s.nspec//2], mid=s.maskdef_id)
print(mode, s.nslits)
