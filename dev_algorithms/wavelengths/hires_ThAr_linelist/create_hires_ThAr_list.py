import sys

import numpy as np
from astropy import table
from astropy import units
from pypeit.core.wave import airtovac
import pypeit.data.arc_lines.convert_NIST_to_lists as nist_to_list
from IPython import embed


# load merged line list from https://www2.keck.hawaii.edu/inst/hires/observing.html
# (only for 3007-7954A air)
lines = table.Table.read('ThAr_list_from_keck_3007-7954A_air.txt', format='ascii')

# convert to vacuum
wave_vac = airtovac(lines['wave'] * units.AA)
lines['wave_vac'] = wave_vac

# work on NIST lines
# convert NIST lines to list
thar_nist = nist_to_list.load_line_list('ThAr')
thar_nist.remove_columns(['Unc.', 'Unc._1', 'Aki'])

# match keck to NIST
thar_matched_nist = nist_to_list.init_line_list()
idx_matched = np.zeros(len(lines), dtype=bool)
for i,w in enumerate(lines['wave_vac']):
    imatch = np.where(np.abs(w - thar_nist['wave'].data) < 0.005)[0]
    if len(imatch) > 0:
        n = imatch[0]
        thar_matched_nist.add_row([thar_nist['Spectrum'][n], thar_nist['wave'][n], 1, 0, thar_nist['Rel.'][n], 'Keck+NIST'])
        idx_matched[i] = True

# print unmatched lines
print('Unmatched lines:')
for i,w in enumerate(lines['wave_vac']):
    if not idx_matched[i]:
        print(lines['wave_vac'][i])

thar_matched_nist.remove_row(0)
thar_matched_nist.sort('wave')

# write to file
outfile = 'ThAr_HIRES_lines.dat'
nist_to_list.write_line_list(thar_matched_nist, outfile, overwrite=False)

