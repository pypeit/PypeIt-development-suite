"""Rerun BuildWaveTilts on GNIRS 32_SB_SXD A_0 with chosen fit and fwhm per order."""
import sys, json, numpy as np
from pypeit.inputfiles import PypeItFile
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import PypeItPar
from pypeit.images import buildimage
from pypeit.slittrace import SlitTraceSet
from pypeit.wavecalib import WaveCalib
from pypeit import wavetilts
from pypeit.core import arc
calrun, fitmode, fwhmsrc = sys.argv[1:4]   # calib images from run; fit A|B; fwhm from develop|arcfitBxa
R='/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_PR/'
C='gemini_gnirs_echelle/32_SB_SXD/'
pf=PypeItFile.from_file(R+'arcfitBxa/'+C+'gemini_gnirs_echelle_32_sb_sxd.pypeit')
spec=load_spectrograph('gemini_gnirs_echelle')
mst=buildimage.TiltImage.from_file(R+calrun+'/'+C+'Calibrations/Tiltimg_A_0_DET01.fits', chk_version=False)
spar=spec.config_specific_par(mst.files[0] if hasattr(mst,'files') and mst.files else None) if False else spec.default_pypeit_par()
import glob
pf_par=sorted(glob.glob(R+'arcfitBxa/'+C+'gemini_gnirs_echelle_32_sb_sxd_UTC_*.par'))[0]
par=PypeItPar.from_cfg_lines(cfg_lines=[l for l in open(pf_par).read().splitlines() if 'fast_line_fit' not in l])
slits=SlitTraceSet.from_file(R+calrun+'/'+C+'Calibrations/Slits_A_all_DET01.fits.gz', chk_version=False)
wv=WaveCalib.from_file(R+fwhmsrc+'/'+C+'Calibrations/WaveCalib_A_0_DET01.fits', chk_version=False)
fwhms=[f.fwhm for f in wv.wv_fits]
if len(sys.argv)>4: fwhms[0]=float(sys.argv[4])
wp=par['calibrations']['wavelengths']
if 'fast_line_fit' in wp.keys(): wp['fast_line_fit']=(fitmode=='B')
sm=slits.slit_img(initial=False)
print('mask flags', slits.mask)
import os
if os.environ.get('SPECORD8'): tp=par['calibrations']['tilts']; tp['spec_order']=[int(os.environ['SPECORD8'])]+[4]*5
b=wavetilts.BuildWaveTilts(mst, slits, spec, par['calibrations']['tilts'], wp, det=1, measured_fwhms=fwhms)
t=b.run(doqa=False)
print("coeffs0", t.coeffs[:,:,0].ravel()[:4], "slitmask after", b.slits.mask, "t.spat_id", t.spat_id, "uniq sm", np.unique(sm)[:8])
tilts=t.fit2tiltimg(sm)
print("tilts nonzero", np.sum(tilts!=0))
ns=tilts.shape[0]; i=0; sid=slits.spat_id[i]
cen=np.round((slits.left_init[:,i]+slits.right_init[:,i])/2).astype(int)
v=tilts[np.arange(ns),cen]*(ns-1)
msk=sm[np.arange(ns),cen]==sid
rows=np.where(msk)[0]
d=v[rows]-rows
print(f'RESULT cal={calrun} fit={fitmode} fwhm={fwhmsrc} fwhm8={fwhms[0]:.3f} rows {rows.min()}-{rows.max()} max|tilt-row| {np.max(np.abs(d)):.2f} at row {rows[np.argmax(np.abs(d))]} tilt(520)-520={v[520]-520:.2f}')
ln=b.all_trace_dict[0] if hasattr(b,'all_trace_dict') else None
if ln is not None:
    print('  lines used', np.round(ln['tilts_spec'][ln['tilts_spec'].shape[0]//2][ln['use_tilt']],1) if 'tilts_spec' in ln else '')
np.save(os.environ.get('OUTD','.')+f'/gn6_fwhm{fwhms[0]:.2f}_so{os.environ.get("SPECORD8","4")}.npy', np.vstack([np.arange(ns), v, msk]))
td=b.all_trace_dict[0] if hasattr(b,'all_trace_dict') and b.all_trace_dict[0] is not None else b.trace_dict
ts=td['tilts_spec']; use=td['use_tilt']
print('LINES fwhm=%s fit=%s: lines_spec'%(fwhmsrc,fitmode), np.round(np.median(ts,axis=0),1), 'use', use.astype(int))
