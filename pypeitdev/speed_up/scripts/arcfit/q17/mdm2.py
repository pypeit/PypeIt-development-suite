import os, sys, pickle
import numpy as np
from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.edgetrace import EdgeTraceSet
from pypeit.core import arc
calls=[]
orig=arc.detect_lines
def rec(censpec,*a,**k):
    out=orig(censpec,*a,**k)
    calls.append((np.array(censpec,copy=True),a,k))
    return out
arc.detect_lines=rec
instrument = load_spectrograph('keck_deimos')
par = instrument.default_pypeit_par()
files=[os.path.join(os.getenv('PYPEIT_DEV'),'RAW_DATA','keck_deimos','1200G_Cooper',f) for f in ['d0115_0023.fits.gz','d0115_0024.fits.gz','d0115_0025.fits.gz']]
tI = buildimage.buildimage_fromlist(instrument,1,par['calibrations']['traceframe'],files)
par = instrument.config_specific_par(tI.files[0])
e = EdgeTraceSet(tI, instrument, par['calibrations']['slitedges'], auto=False)
e.initial_trace(); e.centroid_refine(); e.fit_refine(); e.merge_traces()
e.pca_refine(); e.peak_refine(rebuild_pca=True)
pickle.dump(calls, open('mdm_calls.pkl','wb'))
print(len(calls), [ (c[0].shape, {k:v for k,v in c[2].items() if k!='debug'}) for c in calls])
