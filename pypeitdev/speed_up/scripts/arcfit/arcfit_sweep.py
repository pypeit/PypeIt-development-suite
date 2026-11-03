import time, warnings; warnings.simplefilter('ignore')
import numpy as np
from pypeit.core import arc, fitting
from arcfit_load import load
kw = dict(sigdetect=5., fwhm=4.)
def acc(r, pixt): return (r[2]>0)&(r[2]<12/2.35)&(r[1]>0)&(np.abs(r[1]-pixt)<3)
for which in ['sky','deimos']:
    specs = load(which); calls=[]; orig=arc.fit_arcspec
    arc.fit_arcspec = lambda *a: (calls.append(a), orig(*a))[1]
    for s in specs: arc.detect_lines(s, **kw)
    arc.fit_arcspec = orig
    arc.USE_VECTORIZED_ARCFIT=False
    t=time.perf_counter(); old=[orig(*c) for c in calls]; to=time.perf_counter()-t
    arc.USE_VECTORIZED_ARCFIT=True
    P=np.concatenate([c[2] for c in calls]); O=[np.concatenate([o[i] for o in old]) for i in range(4)]; ao=acc(O,P)
    print(f'== {which}: {P.size} lines; legacy {to:.2f}s; legacy-accepted {ao.sum()}')
    for gn, rs in [(0,0),(3,0),(3,20),(10,20)]:
        arc.ARCFIT_GN_ITER=gn; arc.ARCFIT_RESCUE_ITER=rs
        ncf=[0]; ofg=fitting.fit_gauss
        def cnt(*a,**k): ncf[0]+=1; return ofg(*a,**k)
        fitting.fit_gauss=cnt
        best=1e9
        for _ in range(3):
            t=time.perf_counter(); new=[orig(*c) for c in calls]; best=min(best,time.perf_counter()-t)
        fitting.fit_gauss=ofg
        N=[np.concatenate([n[i] for n in new]) for i in range(4)]; an=acc(N,P)
        both=an&ao; d=np.abs(N[1][both]-O[1][both]); r=np.abs(N[2][both]/O[2][both]-1)
        print(f' GN={gn:2d} R={rs:3d}: {best:.3f}s ({to/best:5.1f}x) curve_fit {ncf[0]//3:4d} | accepted {an.sum()} (lost {np.sum(ao&~an)}, gained {np.sum(an&~ao)}) '
              f'| both: |dc| med {np.median(d):.1e} p99 {np.percentile(d,99):.1e} max {d.max():.2f}, n>0.02 {np.sum(d>0.02)} | |dw/w| p99 {np.percentile(r,99):.1e} | NaN {sum(np.isnan(x).sum() for x in N)}')
arc.ARCFIT_GN_ITER=3; arc.ARCFIT_RESCUE_ITER=100
