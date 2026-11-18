import pickle, numpy as np
from pypeit.core import arc
calls=pickle.load(open('mdm_calls.pkl','rb'))
res={}
for i,(s,a,k) in enumerate(calls):
    k=dict(k); k.pop('debug',None)
    for m in ['A','B']:
        with arc.arcfit_options(match_curvefit=(m=='A')):
            out=arc.detect_lines(s,*a,**k)
        res[(i,m)]=out
    tA=res[(i,'A')]; tB=res[(i,'B')]
    # out: tampl, tampl_cont, tcent, twid, centerr, w, arc, nsig
    cA=tA[2][tA[5]]; cB=tB[2][tB[5]]
    print('call',i,'nA',len(cA),'nB',len(cB))
    onlyB=[c for c in cB if np.min(np.abs(cA-c))>0.5]; onlyA=[c for c in cA if np.min(np.abs(cB-c))>0.5]
    print(' onlyB',np.round(onlyB,2),' onlyA',np.round(onlyA,2))
    # details near 303/306 in rectified coords unknown -> print all close pairs differences
    for c in onlyB+onlyA:
        for m,t in [('A',tA),('B',tB)]:
            j=np.argmin(np.abs(t[2]-c))
            print('  ',m,'cent %.2f wid %.3f ampl %.1f nsig %.1f good %s'%(t[2][j],t[3][j],t[0][j],t[7][j],t[5][j] if t[5].dtype==bool else j in t[5]))
pickle.dump(res,open('mdm_res.pkl','wb'))
