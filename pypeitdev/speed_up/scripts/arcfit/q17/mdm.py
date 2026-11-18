import os, sys, contextlib
import numpy as np
from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.edgetrace import EdgeTraceSet
from pypeit.core import arc
mode = sys.argv[1]
ctx = arc.arcfit_options(match_curvefit=True) if mode=='A' and hasattr(arc,'arcfit_options') else contextlib.nullcontext()
instrument = load_spectrograph('keck_deimos')
par = instrument.default_pypeit_par()
files=[os.path.join(os.getenv('PYPEIT_DEV'),'RAW_DATA','keck_deimos','1200G_Cooper',f) for f in ['d0115_0023.fits.gz','d0115_0024.fits.gz','d0115_0025.fits.gz']]
tI = buildimage.buildimage_fromlist(instrument,1,par['calibrations']['traceframe'],files)
par = instrument.config_specific_par(tI.files[0])
with ctx:
    e = EdgeTraceSet(tI, instrument, par['calibrations']['slitedges'], auto=False)
    e.initial_trace(); e.centroid_refine(); e.fit_refine(); e.merge_traces()
    if e.par['auto_pca'] and e.can_pca():
        e.pca_refine(); e.peak_refine(rebuild_pca=True)
r = e.edge_fit[e.pca.reference_row,:]
np.save(f'mdm_{mode}.npy', np.vstack([r, e.traceid]))
print(mode, e.ntrace, np.round(r[:14],2), np.round(r[-12:],2))
