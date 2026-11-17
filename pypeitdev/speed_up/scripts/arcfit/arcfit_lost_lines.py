import sys, glob, warnings; warnings.simplefilter('ignore')
import numpy as np
from pypeit.core import arc
from pypeit.wavecalib import WaveCalib
specs=[]
for f in sorted(glob.glob(sys.argv[1])):
    wc = WaveCalib.from_file(f, chk_version=False); specs += [wc.arc_spectra[:, i] for i in range(wc.arc_spectra.shape[1])]
calls=[]; orig=arc.fit_arcspec
arc.fit_arcspec = lambda *a: (calls.append(a), orig(*a))[1]
for s in specs: arc.detect_lines(s, sigdetect=5., fwhm=4.)
arc.fit_arcspec = orig
pk_l=[]; pk_all=[]; flat=[]; flat_all=[]; edge=[]
for (x,y,pixt,fitp),s in zip(calls,specs):
    n=orig(x,y,pixt,fitp); arc.USE_VECTORIZED_ARCFIT=False; o=orig(x,y,pixt,fitp); arc.USE_VECTORIZED_ARCFIT=True
    gl=lambda r:(r[2]>0)&(r[2]<12/2.35)&(np.abs(r[1]-pixt)<3)&(r[1]>0)
    lost=gl(o)&~gl(n)
    win=y[np.clip(pixt[:,None]+np.arange(-3,4),0,y.size-1)]
    f=np.sort(win,1)[:,-3]/win.max(1)   # 3rd-highest / max: ~1 for a plateau
    pk_l+=list(s[pixt[lost]]); pk_all+=list(s[pixt]); flat+=list(f[lost]); flat_all+=list(f)
    edge+=list(np.minimum(pixt[lost], y.size-1-pixt[lost]))
print('lost n=%d; raw peak counts median lost %.0f vs all %.0f; frac lost >30k: %.2f'%(len(pk_l),np.median(pk_l),np.median(pk_all),np.mean(np.array(pk_l)>3e4)))
print('plateau metric (3rd/max) median lost %.3f vs all %.3f; lost with metric>0.95: %.2f'%(np.median(flat),np.median(flat_all),np.mean(np.array(flat)>0.95)))
print('distance to array edge, lost: min %d median %d'%(min(edge),np.median(edge)))
