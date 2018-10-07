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
chip = 2
if chip == 1:
    jdict = ltu.loadjson('GMOS_R400_blue.json.gz')
elif chip == 2:
    jdict = ltu.loadjson('MasterWaveCalib_A_02_aa.json')

outroot = 'GMOS_R400_long_'
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


spec = np.array(jdict['0']['spec'])

# Show me
if True:
    xspec = XSpectrum1D.from_tuple((np.arange(len(spec)), spec))
    xspec.plot(xspec=True)
dummy = np.zeros((1024,10))

pdb.set_trace()
#
if chip == 2:
    IDpixels = [1001.4, 942.6, 565.5, 355.3, 130.0]
    IDwaves = [7726.33, 7637.208, 7069.167, 6754.698, 6418.081]
    outfile = outroot+'2_fit.json'

# Line list
CuI = waveio.load_line_list('CuI', use_ion=True, NIST=True)
ArI = waveio.load_line_list('ArI', use_ion=True, NIST=True)
ArII = waveio.load_line_list('ArII', use_ion=True, NIST=True)
llist = vstack([CuI, ArI, ArII])
arcparam['llist'] = llist


# Simple calibration
final_fit = arc.simple_calib(dummy, arcparam, spec, IDpixels=IDpixels, IDwaves=IDwaves, nfitpix=9)#, sigdetect=5.) #sigdetect=7.)
arc.arc_fit_qa(None, final_fit, None, outfile='GMOS_R400_wave.png')
jdict = ltu.jsonify(final_fit)
ltu.savejson(outfile, jdict, overwrite=True)
pdb.set_trace()

# Arc fitter
lines = ['CuI','ArI','ArII']
min_ampl = 1000.
arcfitter = autoid.General(spec.reshape((spec.size, 1)), lines, min_ampl=min_ampl,
                           rms_threshold=0.2, nonlinear_counts=80000.)

pdb.set_trace()

