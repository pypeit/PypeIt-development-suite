

import numpy as np
from astropy import table

from astropy.io import ascii

data = ascii.read('gnirs_ar_argon.dat')

wave_peak = data['col1'].data
nlines = wave_peak.shape[0]
ion = nlines*['Ar']
NIST = nlines*[1]
Instr = nlines*[32]
Source = nlines*['gemini_webpage']
ampl_good = data['col2'].data
b_indx = ampl_good == 'B'
ampl_good[b_indx] = 1.0

file_root_name = 'Ar_IR_GNIRS'
dat = table.Table([wave_peak, ion, NIST, Instr, ampl_good, Source],
            names=('wave', 'ion', 'NIST', 'Instr', 'amplitude', 'Source'))
dat.write(file_root_name + '_lines.dat', format='ascii.fixed_width')
