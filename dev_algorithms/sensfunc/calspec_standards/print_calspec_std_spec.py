# This creates a table with the spectral info (wave range and delta wave) for the calspec standard stars in pypeit/data/standards/calspec


from astropy.io import fits
from astropy.table import Table
import os
from pathlib import Path
import numpy as np
from IPython import embed




calspec_folder = '/Users/dpelliccia/PypeIt/PypeIt/pypeit/data/standards/calspec'
cpath = Path(calspec_folder).resolve()



files = np.sort(sorted(cpath.glob('*.fits.gz')))
exclude = ['README', 'calspec_info.txt', 'hz44_stis_plus_50000k_ext.fits.gz', 'bd_33d2642_004.fits.gz']

fnames = []
names = []
opt_wave_range = []
opt_dwave_range = []

ir_wave_range = []
ir_dwave_range = []


targname = []


for file in files:
	if file.name in exclude:
		continue
	hdu = fits.open(file)
	tab = Table(hdu[1].data)
	wave = tab['WAVELENGTH']
	fnames.append(file.name)
	names.append(file.name.split('_stis')[0].split('_nic')[0].split('_mod')[0])
	opt_wave = wave <= 11000.
	opt_wave_range.append(f'{round(wave[opt_wave].min())} - {round(wave[opt_wave].max())}' if np.any(opt_wave) else '')
	opt_dwave_range.append(f'{round(np.diff(wave[opt_wave]).min(),2)} - {round(np.diff(wave[opt_wave]).max(),2)}' if np.any(opt_wave) else '')

	ir_wave = (wave > 11000.) & (wave <=25000)
	ir_wave_range.append(f'{round(wave[ir_wave].min())} - {round(wave[ir_wave].max())}' if np.any(ir_wave) else '')
	ir_dwave_range.append(f'{round(np.diff(wave[ir_wave]).min(),2)} - {round(np.diff(wave[ir_wave]).max(),2)}' if np.any(ir_wave) else '')
		

# make print table
print_tab = Table()
print_tab['filename'] = fnames
# print_tab['names'] = names
print_tab['opt_wave_range(Angstrom)'] = opt_wave_range
print_tab['opt_dwave_range(Angstrom)'] = opt_dwave_range

print_tab['ir_wave_range(Angstrom)'] = ir_wave_range
print_tab['ir_dwave_range(Angstrom)'] = ir_dwave_range

print('')
print_tab.pprint_all()
print('')
# print_tab.write('new_calspec_files_info_nomod_NEW.csv', format='ascii.csv')

embed()