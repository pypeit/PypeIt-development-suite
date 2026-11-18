import os, sys
import numpy as np
from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.edgetrace import EdgeTraceSet
mode=sys.argv[1]
instrument = load_spectrograph('keck_mosfire')
par = instrument.default_pypeit_par()
files=[os.path.join(os.getenv('PYPEIT_DEV'),'RAW_DATA','keck_mosfire','J_multi',f) for f in ['m191015_0002.fits','m191015_0003.fits','m191015_0004.fits']]
tI = buildimage.buildimage_fromlist(instrument,1,par['calibrations']['traceframe'],files)
par = instrument.config_specific_par(tI.files[0])
e = EdgeTraceSet(tI, instrument, par['calibrations']['slitedges'], auto=False)
e.initial_trace(); e.centroid_refine(); e.fit_refine(); e.merge_traces()
if e.par['auto_pca'] and e.can_pca():
    e.pca_refine(); e.peak_refine(rebuild_pca=True)
r=e.edge_fit[e.pca.reference_row,:]
print(mode, e.ntrace, round(r[0]), round(r[-1]), round(r[9]), round(r[10]))
