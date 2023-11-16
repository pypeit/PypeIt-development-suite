from pypeit import wavecalib, msgs
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from astropy.table import vstack, Table, join
import argparse
from IPython import embed


def load_all_wave_info(redux, input_spec=None):
	"""
	Load all wavecalib files and return a table with all the diagnostic info
	Args:
		redux (str): path to the redux folder
		input_spec (str): name of the spectrograph to be analyzed

	Returns:
		`astropy.table.Table`_: table with all the diagnostic info for all
		the datasets in the redux folder

	"""
	redux_path = Path(redux).resolve()
	spectrographs = sorted(redux_path.glob('*')) if input_spec is None else sorted(redux_path.glob(f'{input_spec}'))
	wtable = Table()
	# loop over spectrographs
	for s in spectrographs:
		if s.name in ['QL_CALIB']:
			continue
		datasets = sorted(s.glob('*'))
		# loop over datasets
		for d in datasets:
			wavefiles = sorted(d.glob('Calibrations/WaveCalib*'))

			# loop over wavefiles
			for w in wavefiles:
				wfile = str(w.resolve())
				try:
					waveCalib = wavecalib.WaveCalib.from_file(wfile, chk_version=False)
				except:
					msgs.warn(f'file {wfile} could not be loaded')

				tab = waveCalib.wave_diagnostics()
				tab['file'] = f'{s.name}/{d.name}/Calibrations/{w.name}'
				wtable = vstack([wtable, tab])

	return wtable


def get_values_for_histo(wtable):
	"""
	Get values, xlabels and bins for the histograms
	Args:
		wtable (`astropy.table.Table`_): table with all the diagnostic info for all

	Returns:
		tuple: values, xlabels, bins
	"""

	# RMS
	rms_diff = wtable['RMS_2'].data.data - wtable['RMS_1'].data.data
	rms_bins = 50 if np.all(rms_diff == 0) else np.arange(rms_diff.min(), rms_diff.max(), 0.05)
	rms_xlabels = r'$\Delta$ RMS pixel (new - old)'
	# Nline
	nlin_diff = wtable['Nlin_2'].data.data - wtable['Nlin_1'].data.data
	nlin_bins = 50 if np.all(nlin_diff == 0) else np.arange(nlin_diff.min(), nlin_diff.max(), 5)
	nlin_xlabels = r'$\Delta$ N lines (new - old)'
	# dWave
	dwave_diff = wtable['dWave_2'].data.data - wtable['dWave_1'].data.data
	dwave_bins = 50 if np.all(dwave_diff == 0) else np.arange(dwave_diff.min(), dwave_diff.max(), 0.05)
	dwave_xlabels = r'$\Delta$ dWave Ang (new - old)'
	# Wave_cen
	wavecen_diff = wtable['Wave_cen_2'].data.data - wtable['Wave_cen_1'].data.data
	wavecen_bins = 50 if np.all(wavecen_diff == 0) else np.arange(wavecen_diff.min(), wavecen_diff.max(), 1.)
	wavecen_xlabels = r'$\Delta$ Wave_cen Ang (new - old)'
	# measured_fwhm
	fwhm_diff = wtable['mesured_fwhm_2'].data.data - wtable['mesured_fwhm_1'].data.data
	fwhm_bins = 50 if np.all(fwhm_diff == 0) else np.arange(fwhm_diff.min(), fwhm_diff.max(), 0.5)
	fwhm_xlabels = r'$\Delta$ measured_fwhm pix (new - old)'

	# put everything in lists
	values = [rms_diff, nlin_diff, dwave_diff, wavecen_diff, fwhm_diff]
	xlabels = [rms_xlabels, nlin_xlabels, dwave_xlabels, wavecen_xlabels, fwhm_xlabels]
	bins = [rms_bins, nlin_bins, dwave_bins, wavecen_bins, fwhm_bins]

	return values, xlabels, bins


def get_values_for_plt(wtable):
	"""
	Get values, xlabels and ylabels for the plots

	Args:
		wtable (`astropy.table.Table`_): table with all the diagnostic info for all

	Returns:
		tuple: old_values, new_values, xlabels, ylabels
	"""

	# RMS
	rms_xlabels = r'OLD RMS pixel'
	rms_ylabels = r'NEW RMS pixel'
	# Nline
	nlin_xlabels = r'OLD N lines'
	nlin_ylabels = r'NEW N lines'
	# dWave
	dwave_xlabels = r'OLD dWave Ang'
	dwave_ylabels = r'NEW dWave Ang'
	# Wave_cen
	wavecen_xlabels = r'OLD Wave_cen Ang'
	wavecen_ylabels = r'NEW Wave_cen Ang'
	# measured_fwhm
	fwhm_xlabels = r'OLD measured_fwhm pix'
	fwhm_ylabels = r'NEW measured_fwhm pix'

	# put everything in lists
	old_values = [wtable['RMS_1'].data.data, wtable['Nlin_1'].data.data, wtable['dWave_1'].data.data,
				  wtable['Wave_cen_1'].data.data, wtable['mesured_fwhm_1'].data.data]
	new_values = [wtable['RMS_2'].data.data, wtable['Nlin_2'].data.data, wtable['dWave_2'].data.data,
				  wtable['Wave_cen_2'].data.data, wtable['mesured_fwhm_2'].data.data]
	xlabels = [rms_xlabels, nlin_xlabels, dwave_xlabels, wavecen_xlabels, fwhm_xlabels]
	ylabels = [rms_ylabels, nlin_ylabels, dwave_ylabels, wavecen_ylabels, fwhm_ylabels]

	return old_values, new_values, xlabels, ylabels


def get_parser():
	parser = argparse.ArgumentParser(description='Compare wavelength calibration diagnostics between two PypeIt runs')
	parser.add_argument('old_redux', type=str, help='Path to the old redux folder')
	parser.add_argument('new_redux', type=str, help='Path to the new redux folder')
	parser.add_argument('--spec', default=None, type=str, help='Name of the spectrograph to be compared')
	parser.add_argument('--print_tab', default=False, help='Print combined table in the terminal', action='store_true')
	parser.add_argument('--embed', default=False, action='store_true')
	return parser.parse_args()


def main(args):
	"""

	Args:
		args: arguments from the command line

	Returns:
		`astropy.table.Table`_: table of the combined wavelength diagnostics from the two PypeIt runs

	"""

	# check that the two redux folders exist
	if not Path(args.old_redux).exists():
		msgs.error(f'Folder "{args.old_redux}" does not exist')
	if not Path(args.new_redux).exists():
		msgs.error(f'Folder "{args.new_redux}" does not exist')

	wtable_old = load_all_wave_info(args.old_redux, input_spec=args.spec)
	wtable_new = load_all_wave_info(args.new_redux, input_spec=args.spec)

	combined_wtable = join(wtable_old, wtable_new, join_type='left', keys=['file', 'SpatID'])
	if args.print_tab:
		combined_wtable.pprint_all()

	masked_rows = combined_wtable.mask['RMS_2'] == True
	msgs.info(f'\n      **{np.sum(masked_rows)} slits are not present in the new run**')
	notmasked_combined_wtable = combined_wtable[np.logical_not(masked_rows)]

	# values for the histograms
	values_histo, xlabels_histo, bins_histo = get_values_for_histo(notmasked_combined_wtable)
	# values for the plots
	old_values, new_values, xlabels_plt, ylabels_plt = get_values_for_plt(notmasked_combined_wtable)

	# 1) plot
	nrows, ncols = 1, len(values_histo)
	fig = plt.figure(figsize=(5 * ncols, 6 * nrows))
	gs = gridspec.GridSpec(nrows, ncols)
	for ii in range(nrows * ncols):
		ax = plt.subplot(gs[ii // ncols, ii % ncols])
		ax.minorticks_on()
		ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
		try:
			y, x, _ = ax.hist(values_histo[ii], bins=bins_histo[ii], histtype='stepfilled', align='left')
			ax.axvline(0, color='k', linestyle='--', lw=1.5)
			ax.set_xlabel(xlabels_histo[ii])
			ax.set_ylabel('N')
			ax.set_ylim(0, y.max() * 0.9)
		except:
			continue

	fig.tight_layout(pad=1, h_pad=0.0, w_pad=1)

	# 2) plot
	nrows2, ncols2 = 1, len(old_values)
	fig2 = plt.figure(figsize=(5 * ncols2, 6 * nrows2))
	gs2 = gridspec.GridSpec(nrows2, ncols2)
	for jj in range(nrows2 * ncols2):
		ax2 = plt.subplot(gs2[jj // ncols2, jj % ncols2])
		ax2.minorticks_on()
		ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True)
		line = np.arange(0., np.max([old_values[jj].max(), new_values[jj].max()]), 0.1)
		ax2.plot(line, line, color='k', linestyle='--', lw=1.5)
		ax2.scatter(old_values[jj], new_values[jj], s=10)
		ax2.set_xlabel(xlabels_plt[jj])
		ax2.set_ylabel(ylabels_plt[jj])

	fig2.tight_layout(pad=1, h_pad=0.0, w_pad=1)

	plt.show()

	if args.embed:
		embed()

	return combined_wtable


if __name__ == '__main__':
	args = get_parser()
	main(args)


