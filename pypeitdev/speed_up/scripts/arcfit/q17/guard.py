import glob, os, re, json, numpy as np
from pypeit.wavecalib import WaveCalib
R='/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_PR/'
out=[]
for run in ['develop','arcfitBxa']:
    for pf in sorted(glob.glob(R+run+'/*/*/*_UTC_*.par')):
        txt=open(pf).read()
        m=re.search(r'\[\[wavelengths\]\](.*?)\n    \[\[', txt, re.S)
        if not m: continue
        blk=m.group(1)
        meth=re.search(r'\n +method = (\S+)',blk).group(1)
        if meth!='full_template': continue
        if re.search(r'\n +echelle = True',blk): continue
        thr=float(re.search(r'\n +rms_thresh_frac_fwhm = (\S+)',blk).group(1))
        sdir=os.path.dirname(pf)
        for wf in sorted(glob.glob(sdir+'/Calibrations/WaveCalib_*.fits')):
            try: w=WaveCalib.from_file(wf, chk_version=False)
            except Exception as e: continue
            ns=w.arc_spectra.shape[0]; x=np.arange(ns)/(ns-1)
            for i,f in enumerate(w.wv_fits):
                if f is None or f.pypeitfit is None: continue
                wv=f.pypeitfit.eval(x)
                out.append(dict(run=run,setup=sdir.split(run+'/')[1],key=os.path.basename(wf)[9:-5],i=i,sid=int(w.spat_ids[i]),
                                rms=float(f.rms),fwhm=float(f.fwhm) if f.fwhm else None,thr=thr,
                                mono=bool(np.all(np.diff(wv)>0) or np.all(np.diff(wv)<0)),nlin=len(f.pixel_fit)))
json.dump(out,open('guard.json','w'))
for run in ['develop','arcfitBxa']:
    r=[x for x in out if x['run']==run]
    bad_rms=[x for x in r if x['fwhm'] and x['rms']>x['thr']*x['fwhm']]
    bad_mono=[x for x in r if not x['mono']]
    print(run,'full_template slits',len(r),'setups',len({x['setup'] for x in r}),'| rms>thr',len(bad_rms),'| non-monotonic',len(bad_mono),'| either',len({(x['setup'],x['key'],x['i']) for x in bad_rms+bad_mono}))
    for x in bad_mono: print('   nonmono', x['setup'],x['key'],x['i'],'rms %.2f thr %.2f nlin %d'%(x['rms'],x['thr']*x['fwhm'],x['nlin']))
