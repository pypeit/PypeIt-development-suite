import sys
from astropy import table
from astropy import units
from pypeit.core.wave import airtovac
import pypeit.data.arc_lines.convert_NIST_to_lists as nist_to_list


# load line list from https://www2.keck.hawaii.edu/inst/lris/arc_calibrations.html
lines = table.Table.read('fene.dat', format='ascii')

# convert to vacuum
wave_vac = airtovac(lines['wave'] * units.AA)
lines['wave_vac'] = wave_vac

# work on FeI lines
fe = lines['ion'] == 'FeI'
# convert FeI NIST to list
sys.argv = ['../PypeIt/PypeIt/pypeit/data/arc_lines/convert_NIST_to_lists.py', 'FeI', '--skip_stop', '-r -1']
fe_nist = nist_to_list.main()

# match to FeI in lines
fe_matched_nist = nist_to_list.init_line_list()
for w in lines[fe]['wave_vac']:
    for n,nist_w in enumerate(fe_nist['wave']):
        if abs(nist_w - w)<0.0005:
            fe_matched_nist.add_row([fe_nist['ion'][n], fe_nist['wave'][n], fe_nist['NIST'][n], fe_nist['Instr'][n], fe_nist['amplitude'][n], fe_nist['Source'][n]])
fe_matched_nist.remove_row(0)
fe_matched_nist.sort('wave')

# work on NeI lines
neI = lines['ion'] == 'NeI'
# convert NeI NIST to list
sys.argv = ['../PypeIt/PypeIt/pypeit/data/arc_lines/convert_NIST_to_lists.py', 'NeI', '--skip_stop', '-r -1']
neI_nist = nist_to_list.main()

# match to NeI in lines
neI_matched_nist = nist_to_list.init_line_list()
for w in lines[neI]['wave_vac']:
    for n,nist_w in enumerate(neI_nist['wave']):
        if abs(nist_w - w)<0.0005:
            neI_matched_nist.add_row([neI_nist['ion'][n], neI_nist['wave'][n], neI_nist['NIST'][n], neI_nist['Instr'][n], neI_nist['amplitude'][n], neI_nist['Source'][n]])
neI_matched_nist.remove_row(0)
neI_matched_nist.sort('wave')


# work on NeII lines
neII = lines['ion'] == 'NeII'
# convert NeII NIST to list
sys.argv = ['../PypeIt/PypeIt/pypeit/data/arc_lines/convert_NIST_to_lists.py', 'NeII', '--skip_stop', '-r -1']
neII_nist = nist_to_list.main()

# match to NeII in lines
neII_matched_nist = nist_to_list.init_line_list()
for w in lines[neII]['wave_vac']:
    for n,nist_w in enumerate(neII_nist['wave']):
        if abs(nist_w - w)<0.05:
            neII_matched_nist.add_row([neII_nist['ion'][n], neII_nist['wave'][n], neII_nist['NIST'][n], neII_nist['Instr'][n], neII_nist['amplitude'][n], neII_nist['Source'][n]])
neII_matched_nist.remove_row(0)
neII_matched_nist.sort('wave')


# stack all matched
all_matched_nist = table.vstack([fe_matched_nist, neI_matched_nist, neII_matched_nist])
all_matched_nist.sort('wave')

# write to file
outfile = 'FeAr_lines.dat'
nist_to_list.write_line_list(all_matched_nist, outfile, overwrite=False)

