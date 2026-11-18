import sys, glob, numpy as np
from pypeit.par import PypeItPar
from pypeit.images import buildimage
from pypeit.slittrace import SlitTraceSet
from pypeit.core.wavecal import autoid
from pypeit.core import arc
from pypeit.spectrographs.util import load_spectrograph
R='/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_PR/'
C='gemini_gnirs_echelle/32_SB_SXD/'
pfp=sorted(glob.glob(R+'arcfitBxa/'+C+'gemini_gnirs_echelle_32_sb_sxd_UTC_*.par'))[0]
par=PypeItPar.from_cfg_lines(cfg_lines=[l for l in open(pfp).read().splitlines() if 'fast_line_fit' not in l])
wp=par['calibrations']['wavelengths']
print('wp fwhm',wp['fwhm'],'sigdetect',wp['sigdetect'],'orders',wp['fwhm_spec_order'],wp['fwhm_spat_order'],'method',wp['method'])
msarc=buildimage.ArcImage.from_file(R+'arcfitBxa/'+C+'Calibrations/Arc_A_0_DET01.fits', chk_version=False)
slits=SlitTraceSet.from_file(R+'arcfitBxa/'+C+'Calibrations/Slits_A_all_DET01.fits.gz', chk_version=False)
left,right,_=slits.select_edges(initial=True)
sm=slits.slit_img(initial=True)
spec=load_spectrograph('gemini_gnirs_echelle')
gpm=np.logical_not(msarc.select_flag()) if hasattr(msarc,'select_flag') else np.ones(msarc.image.shape,bool)
if wp['method']!='full_template':
    gpm &= msarc.image < spec.get_detector_par(1)['saturation']*spec.get_detector_par(1)['nonlinear']
out={}
for mode in ['A','B']:
    rec=[]
    orig=arc.detect_lines
    def dl(*a,**k):
        r=orig(*a,**k); rec.append((r[2][r[5]],r[3][r[5]],r[7][r[5]])); return r
    arc.detect_lines=dl
    with arc.arcfit_options(match_curvefit=(mode=='A')):
        fm=autoid.map_fwhm(msarc.image,gpm,left,right,sm,nsample=10,slit_bpm=np.zeros(slits.nslits,bool),
                           specord=wp['fwhm_spec_order'],spatord=wp['fwhm_spat_order'],sigdetect=wp['sigdetect'],fwhm=wp['fwhm'])
    arc.detect_lines=orig
    vals=[f.eval(msarc.image.shape[0]//2,0.5) for f in fm]
    print(mode, 'fwhm at centre per order', np.round(vals,3))
    # order 8 = slit 0 -> first 10 recorded calls
    cen=np.concatenate([r[0] for r in rec[:10]]); wid=np.concatenate([r[1] for r in rec[:10]]); ns=np.concatenate([r[2] for r in rec[:10]])
    samp=np.concatenate([np.full(r[0].size,i) for i,r in enumerate(rec[:10])])
    out[mode]=dict(cen=cen,fwhm=wid*2.3548,nsig=ns,samp=samp,gpm=fm[0].bool_gpm,coef=fm[0].fitc, val=vals[0], curve=np.array([fm[0].eval(r,0.5) for r in range(msarc.image.shape[0])]))
    print('  order8 nlines', cen.size, 'used', fm[0].bool_gpm.sum(), 'fit eval at rows 511/600/800/1000', [round(float(fm[0].eval(r,0.5)),2) for r in [511,600,800,1000]])
import os; np.save(os.environ.get('OUTD','.')+'/gn7.npy', out, allow_pickle=True)
a,b=out['A'],out['B']
print('A cen range', a['cen'].min(), a['cen'].max()); print('B cen range', b['cen'].min(), b['cen'].max())
# lines only in B
for m,x,y in [('onlyB',b,a),('onlyA',a,b)]:
    sel=[i for i in range(x['cen'].size) if not np.any((y['samp']==x['samp'][i]) & (np.abs(y['cen']-x['cen'][i])<0.5))]
    print(m, len(sel), np.round(x['cen'][sel],1), np.round(x['fwhm'][sel],2), 'used', x['gpm'][sel].astype(int))
