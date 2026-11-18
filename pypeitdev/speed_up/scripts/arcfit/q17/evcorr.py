import sys, pickle, numpy as np, contextlib
from pypeit.core import arc
caps=pickle.load(open('corrcaps.pkl','rb'))
s,k=caps[0]; k=dict(k); k.pop('min_frac_med_width',None)
modes=['legacy'] if not hasattr(arc,'arcfit_options') else ['A','B']
for m in modes:
    ctx=contextlib.nullcontext() if m=='legacy' else arc.arcfit_options(match_curvefit=(m=='A'))
    with ctx:
        k2=dict(k); k2['nfind']=None
        r=arc.detect_lines(s,**k2)
        r1=arc.detect_lines(s,**k)
    tampl_true,tampl,cen,wid,cerr,ww,_,nsig=r
    o=np.argsort(-tampl)[:6]
    print(m,'nfind=1 ->', np.round(r1[2],2), '| top peaks (cen,ampl,wid,nsig,good):', [(round(cen[i],1),round(tampl[i],3),round(wid[i],2),round(nsig[i],1), i in ww) for i in o])
