import sys, pickle, runpy, numpy as np
from pypeit.core import arc
from pypeit.core.wavecal import wvutils
caps=[]
orig=arc.detect_lines
def dl(spec,*a,**k):
    if k.get('nfind')==1: caps.append((np.array(spec),dict(k)))
    return orig(spec,*a,**k)
arc.detect_lines=dl
sys.argv=['ft.py']+sys.argv[1:]
try:
    runpy.run_path('ft.py', run_name='__main__')
finally:
    pickle.dump(caps, open('corrcaps.pkl','wb'))
