''' Script to develop wavelengths for GMOS
'''
from __future__ import absolute_import, division, print_function

import pdb

import numpy as np

from linetools import utils as ltu

from pypeit.core.wavecal import autoid

# Load the spectra
jdict = ltu.loadjson('GMOS_R400_blue.json.gz')

arccen = np.array(jdict['arccen'])

# Arc fitter
spec = arccen[:,1]
lines = ['CuI','ArI','ArII']
min_ampl = 1000.
arcfitter = autoid.General(spec.reshape((spec.size, 1)), lines, min_ampl=min_ampl,
                           rms_threshold=0.2, nonlinear_counts=80000.)

pdb.set_trace()

