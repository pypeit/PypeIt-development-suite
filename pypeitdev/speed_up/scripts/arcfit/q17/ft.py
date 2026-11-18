"""Re-run full_template on one stored arc spectrum.
usage: ft.py SETUP KEY SLITIDX VARIANT [FWHM_FROM_RUN] [ARC_FROM_RUN]
VARIANT: native (whatever code is on PYTHONPATH), B (arcfit code, B also in get_xcorr_arc),
         A (arcfit code, A everywhere)."""
import sys, glob, contextlib, types, json
import numpy as np
from pypeit.par import PypeItPar
from pypeit.wavecalib import WaveCalib
from pypeit.core.wavecal import autoid, wvutils
from pypeit.core import arc
setup, key, si, variant = sys.argv[1:5]; si=int(si)
fwrun = sys.argv[5] if len(sys.argv)>5 else 'develop'
arcrun = sys.argv[6] if len(sys.argv)>6 else fwrun
R='/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_PR/'
pfp=sorted(glob.glob(R+'arcfitBxa/'+setup+'/*_UTC_*.par'))[0]
par=PypeItPar.from_cfg_lines(cfg_lines=[l for l in open(pfp).read().splitlines() if 'fast_line_fit' not in l])
wp=par['calibrations']['wavelengths']
wa=WaveCalib.from_file(R+arcrun+'/'+setup+f'/Calibrations/WaveCalib_{key}.fits', chk_version=False)
wf=WaveCalib.from_file(R+fwrun+'/'+setup+f'/Calibrations/WaveCalib_{key}.fits', chk_version=False)
spec=wa.arc_spectra[:,[si]]
fwhm=wf.wv_fits[si].fwhm
ctx=contextlib.nullcontext()
if variant=='B':
    class P:
        def __getattr__(self,n): return getattr(arc,n)
        def arcfit_options(self,**k): return contextlib.nullcontext()
    wvutils.arc=P()
elif variant=='A':
    ctx=arc.arcfit_options(match_curvefit=True)
from pypeit.images import buildimage
from pypeit.core import parse
akey=key.rsplit('_',1)
am=buildimage.ArcImage.from_file(R+arcrun+'/'+setup+f'/Calibrations/Arc_{key}.fits', chk_version=False)
binspec,_=parse.parse_binning(am.detector.binning)
nonlin=am.detector.saturation*am.detector.nonlinear
LOG=[]
import os
from pypeit.core.wavecal import waveio
if os.environ.get('FIXTPL', '0') == '1':
    _lt=waveio.load_template
    def lt(*a,**k):
        r=list(_lt(*a,**k)); w=r[0]
        seams=np.where(np.diff(w)<=0)[0]
        if seams.size: r[0]=w[:seams[0]+1]; r[1]=r[1][:seams[0]+1]
        LOG.append(('template', len(w), '->', len(r[0])))
        return tuple(r)
    waveio.load_template=lt
_xs=wvutils.xcorr_shift
def xs(*a,**k):
    r=_xs(*a,**k); LOG.append(('xcorr_shift',round(float(r[0]),3),round(float(r[1]),4))); return r
wvutils.xcorr_shift=xs
_xss=wvutils.xcorr_shift_stretch
def xss(*a,**k):
    r=_xss(*a,**k); LOG.append(('shift_stretch',)+tuple(np.round(np.atleast_1d(np.asarray(r[:4],dtype=float)),4))); return r
wvutils.xcorr_shift_stretch=xss
_rid=autoid.reidentify
def rid(*a,**k):
    r=_rid(*a,**k); ids=r[2].get('IDs',np.array([])) if isinstance(r[2],dict) else np.array([])
    LOG.append(('reidentify ndet',len(r[0]),'nID',int(np.sum(np.asarray(ids)>0)))); return r
autoid.reidentify=rid
with ctx:
    wvc=autoid.full_template(spec, wp['lamps'], wp, np.array([0]), 1, binspec, nsnippet=wp['nsnippet'],
                              slit_ids=np.array([wa.spat_ids[si]]), measured_fwhms=np.array([fwhm]), nonlinear_counts=nonlin)
wvc=wvc[0]
d=wvc['0'] if '0' in wvc else list(wvc.values())[0]
x=np.arange(spec.shape[0])
w=d.pypeitfit.eval(x/(spec.shape[0]-1)) if d is not None and d.pypeitfit is not None else np.zeros_like(x)
wref=wf.wv_fits[si]
out=dict(variant=variant, fwrun=fwrun, arcrun=arcrun, fwhm=float(fwhm), nlin=int(len(d.pixel_fit)) if d is not None else 0, rms=float(d.rms) if d is not None else -1)
np.save(os.environ.get("OUTD",".")+f"/ft_{os.environ.get('FIXTPL','0')}_{setup.replace('/','_')}_{si}_{variant}_{fwrun}_{arcrun}.npy", dict(out, wave=w, pix=None if d is None else d.pixel_fit, wfit=None if d is None else d.wave_fit))
[print('  ',l) for l in LOG]
print('FT', json.dumps(out), 'wave[0,mid,-1]', np.round([w[0], w[len(w)//2], w[-1]],1))
