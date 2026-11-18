"""full_template on every slit of a setup (one process).  usage: ftall.py SETUP KEY TAG [FIXTPL=0/1]
Uses the arc spectra / fwhm of the run matching the code (develop->develop, arcfit->arcfitBxa)."""
import sys, glob, os, json, numpy as np
from pypeit.par import PypeItPar
from pypeit.wavecalib import WaveCalib
from pypeit.core.wavecal import autoid, waveio
from pypeit.images import buildimage
from pypeit.core import parse
setup,key,tag,fix=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]=='1'
R='/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_PR/'
run='develop' if tag=='dev' else 'arcfitBxa'
pfp=sorted(glob.glob(R+'arcfitBxa/'+setup+'/*_UTC_*.par'))[0]
par=PypeItPar.from_cfg_lines(cfg_lines=[l for l in open(pfp).read().splitlines() if 'fast_line_fit' not in l])
wp=par['calibrations']['wavelengths']
w=WaveCalib.from_file(R+run+'/'+setup+f'/Calibrations/WaveCalib_{key}.fits', chk_version=False)
wd=WaveCalib.from_file(R+'develop/'+setup+f'/Calibrations/WaveCalib_{key}.fits', chk_version=False)
am=buildimage.ArcImage.from_file(R+run+'/'+setup+f'/Calibrations/Arc_{key}.fits', chk_version=False)
binspec,_=parse.parse_binning(am.detector.binning); nonlin=am.detector.saturation*am.detector.nonlinear
if fix:
    _lt=waveio.load_template
    def lt(*a,**k):
        r=list(_lt(*a,**k)); ww=r[0]; s=np.where(np.diff(ww)<=0)[0]
        if s.size: r[0]=ww[:s[0]+1]; r[1]=r[1][:s[0]+1]
        return tuple(r)
    waveio.load_template=lt
ok=np.array([i for i,f in enumerate(w.wv_fits) if f is not None and f.fwhm is not None])
fw=np.array([f.fwhm if f is not None else None for f in w.wv_fits],dtype=object)
wvc,_=autoid.full_template(w.arc_spectra, wp['lamps'], wp, ok, 1, binspec, nsnippet=wp['nsnippet'],
                         slit_ids=w.spat_ids, measured_fwhms=fw, nonlinear_counts=nonlin)
ns=w.arc_spectra.shape[0]; x=np.arange(ns)/(ns-1)
rows=[]
for i in range(len(w.wv_fits)):
    d=wvc.get(str(i)); fd=wd.wv_fits[i]
    if d is None or fd is None or fd.pypeitfit is None: rows.append(dict(i=i,ok=False)); continue
    wv=d.pypeitfit.eval(x); wref=fd.pypeitfit.eval(x)
    disp=np.median(np.abs(np.diff(wref)))
    rows.append(dict(i=i,ok=True,nlin=len(d.pixel_fit),rms=float(d.rms),nlin_dev=len(fd.pixel_fit),rms_dev=float(fd.rms),
                     dpix=float(np.max(np.abs(wv-wref))/disp), mono=bool(np.all(np.diff(wv)>0)), fwhm=float(fw[i])))
json.dump(rows,open(f"ftall_{setup.replace('/','_')}_{tag}_fix{int(fix)}.json",'w'))
for r in rows:
    if r['ok']: print('ROW %2d nlin %3d rms %7.3f | dev nlin %3d rms %6.3f | dpix %8.2f mono %s'%(r['i'],r['nlin'],r['rms'],r['nlin_dev'],r['rms_dev'],r['dpix'],r['mono']))
