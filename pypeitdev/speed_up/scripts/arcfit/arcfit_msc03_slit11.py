"""
Task 6 item 1: reproduce the keck_deimos/600ZD_M_6500 MSC03 slit 11 (spat_id 1063)
difference between legacy and variant A by re-running autoid.full_template on the
stored arc spectrum, crossing (fwhm from legacy run, fwhm from A run) x
(USE_VECTORIZED_ARCFIT False/True).
Usage: python arcfit_msc03_slit11.py [REDUX_OUT_C]
"""
import sys, warnings
from pathlib import Path
warnings.simplefilter('ignore')
import numpy as np
from pypeit import log
log.setLevel('ERROR')
from pypeit.wavecalib import WaveCalib
from pypeit.images.buildimage import ArcImage
from pypeit.par.pypeitpar import WavelengthSolutionPar
from pypeit.core import arc
from pypeit.core.wavecal import autoid

root = Path(sys.argv[1] if len(sys.argv) > 1 else
            '/home/xavier/Projects/PypeIt/PypeIt-development-suite/REDUX_OUT_C')
rel = 'keck_deimos/600ZD_M_6500/Calibrations/WaveCalib_A_0_MSC03.fits'
islit = 11
W = {v: WaveCalib.from_file(str(root / v / rel), chk_version=False) for v in ['legacy', 'A']}
msarc = ArcImage.from_file(str(root / 'A' / rel.replace('WaveCalib', 'Arc')), chk_version=False)
nlc = np.max([d.nonlinear_counts() for d in msarc.detector.detectors])
par = WavelengthSolutionPar.from_dict(dict(W['A'].par))
spec = W['A'].arc_spectra[:, islit]
assert np.array_equal(spec, W['legacy'].arc_spectra[:, islit])
fw = {v: W[v].wv_fits[islit].fwhm for v in W}
print(f'nonlinear_counts {nlc}; fwhm legacy {fw["legacy"]!r} A {fw["A"]!r}')
for fv in ['legacy', 'A']:
    for vec in [False, True]:
        arc.USE_VECTORIZED_ARCFIT = vec
        out, _ = autoid.full_template(spec, par['lamps'], par, [0], 1, 1, slit_ids=[1063],
                                      measured_fwhms=[fw[fv]], nonlinear_counts=nlc,
                                      nsnippet=par['nsnippet'], x_percentile=par['cc_percent_ceil'])
        f = out['0']
        g = f.pypeitfit.gpm.astype(bool)
        print(f'fwhm={fv:6s} vectorized={vec!s:5s}: rms {f.rms:.4f} nfit {len(f.pixel_fit)} ngood {g.sum()} '
              f'bluest good line {np.min(np.asarray(f.pixel_fit)[g]):.3f} wave(0) {f.wave_soln[0]:.2f}')
arc.USE_VECTORIZED_ARCFIT = True
