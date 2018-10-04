''' Script to develop wavelengths for GMOS
'''
from __future__ import absolute_import, division, print_function

import pdb

import numpy as np

from astropy.table import vstack

from linetools import utils as ltu
from linetools.spectra.xspectrum1d import XSpectrum1D

from pypeit.core.wavecal import autoid
from pypeit.core.wavecal import waveio
from pypeit.core import arc
from pypeit.spectrographs import gemini_gmos

# Load the spectra
jdict = ltu.loadjson('GMOS_R400_blue.json.gz')
spectrograph = gemini_gmos.GeminiGMOSNE2VSpectrograph()
arcparam = {}
spectrograph.setup_arcparam(arcparam,disperser='R400')
arcparam['n_first'] = 2
arcparam['n_final'] = 3
arcparam['func'] = 'legendre'
arcparam['nsig_rej'] = 2.
arcparam['nsig_rej_final'] = 3.
arcparam['disp'] = 0.67*2.
arcparam['match_toler'] = 3.
arcparam['disp_toler'] = 0.1
arcparam['Nstrong'] = 13

arccen = np.array(jdict['arccen'])

spec = arccen[:,1]
# Show me
if False:
    xspec = XSpectrum1D.from_tuple((np.arange(len(spec)), spec))
    xspec.plot(xspec=True)
dummy = np.zeros((1024,10))

IDpixels = [801.7, 636.01, 487.28, 322.49, 31.8]
IDwaves = [5189.191, 4966.465, 4766.197, 4546.3258, 4159.762]

# Line list
CuI = waveio.load_line_list('CuI', use_ion=True, NIST=True)
ArI = waveio.load_line_list('ArI', use_ion=True, NIST=True)
ArII = waveio.load_line_list('ArII', use_ion=True, NIST=True)
llist = vstack([CuI, ArI, ArII])
arcparam['llist'] = llist


# Simple calibration
final_fit = arc.simple_calib(dummy, arcparam, spec, IDpixels=IDpixels, IDwaves=IDwaves, nfitpix=9, sigdetect=4.)
arc.arc_fit_qa(None, final_fit, None, outfile='GMOS_R400_wave.png')
jdict = ltu.jsonify(final_fit)
ltu.savejson('GMOS_R400_blue_1_fit.json', jdict)
pdb.set_trace()

# Arc fitter
lines = ['CuI','ArI','ArII']
min_ampl = 1000.
arcfitter = autoid.General(spec.reshape((spec.size, 1)), lines, min_ampl=min_ampl,
                           rms_threshold=0.2, nonlinear_counts=80000.)

pdb.set_trace()

