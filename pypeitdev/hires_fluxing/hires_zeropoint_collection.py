# compute the sensitivity functions for a large number of reduced HIRES stars
# located in the Google Drive, and save them in the working directory (`working_dir`)
# the script `hires_plot_zeropoints.py` will collect all the existing sensitivity functions
# and plot the zeropoints in a single figure.

import shutil
import numpy as np
from pathlib import Path
from pypeit import msgs
from pypeit.scripts import sensfunc
from IPython import embed

skip_existing = True  # skip computing the sensitivity function if the output file already exists
boxcar = False  # use a boxcar extraction for the sensitivity function computation


# list of reduced stars to collect zeropoints from
star_list = ['Feige34', 'Feige110', 'G191B2B']

# path to the HIRES reduction of those stars on the Google Drive
hredux_path = Path('/Volumes/GoogleDrive/Shared drives/PypeIt ADAP 2020/backups/HIRES')

# working directory (change this to your working directory)
working_dir = Path('/Users/dpelliccia/Desktop/adap2020')

# check if the path to the HIRES reduction directory exists
if not hredux_path.is_dir():
    msgs.error('Path to the HIRES reduction directory does not exist.')

# loop over the stars
skipped = []
reasons = []
for star in star_list:
    msgs.info(f'\n\nWorking on star {star}')
    # get the path to the star reduction on the Google Drive
    star_path = hredux_path / star
    if not star_path.is_dir():
        msgs.error(f'Path to the star directory {star} does not exist.')
    # set the path to the star reduction on the working directory
    workstar_dir = working_dir / star
    # check if the path to the star reduction on the working directory exists otherwise create it
    if not workstar_dir.is_dir():
        workstar_dir.mkdir()
    # get the paths to the reduction directories for this star
    redux_paths = list(star_path.glob('*'))
    # loop over the reduction directories
    for redux_path in redux_paths:
        if not redux_path.is_dir():
            continue
        # get each run of the reduction for this star
        runs = list(redux_path.glob('*/complete/reduce/keck_hires_*'))
        # loop over the runs
        for run in runs:
            if not run.is_dir():
                continue
            # get the spec1d files for this run and check if they exist
            spec1ds = list(run.glob('Science/spec1d*.fits'))
            if len(spec1ds) == 0:
                msgs.warn(f'No spec1d files found in {run}.')
                skipped.append(str(run.parent.parent.parent))
                reasons.append('No spec1d files')
                continue
            # use the first spec1d file to compute the sensitivity function
            spec1d = spec1ds[0]

            # set the path to the pypeit_sensfunc run on the working directory
            fname = '_'.join(str(runs[0]).split(f'{star}/')[1].split('/')[0:2])
            workrun_dir = workstar_dir / fname
            # check if the path to the pypeit_sensfunc run on the working directory exists otherwise create it
            if not workrun_dir.is_dir():
                workrun_dir.mkdir()
            elif len(list(workrun_dir.glob('sens*.fits'))) > 0 and skip_existing:
                msgs.info(f'Sensitivity function for {workrun_dir.name} already exists. Skipping.')
                continue

            # set the path to the sensitivity function file
            out_sensfile = workrun_dir / spec1d.name.replace('spec1d','sens')
            # set the path to the sensitivity function parameters file
            par_outfile = workrun_dir / 'sensfunc.par'
            msgs.info(f'\n\nComputing sensitivity function for {workrun_dir.name}')
            try:
                # compute the sensitivity function
                parser = sensfunc.SensFunc.get_parser()
                if boxcar:
                    msgs.info('Using BOXCAR extraction for the sensitivity function computation.')
                    args = parser.parse_args([str(spec1d), '-f', '-o', str(out_sensfile),
                                              '--par_outfile', str(par_outfile), '--extr', 'BOX'])
                else:
                    args = parser.parse_args([str(spec1d), '-f', '-o', str(out_sensfile),
                                              '--par_outfile', str(par_outfile)])
                sensfunc.SensFunc.main(args)
                msgs.info(f'Sensitivity function for {workrun_dir.name} computed.')
            except Exception as e:
                msgs.warn(f'Error computing sensitivity function for {workrun_dir.name}.\n{e}')
                skipped.append(str(run.parent.parent.parent))
                reasons.append('Error computing sensitivity function')
            if not out_sensfile.is_file():
                msgs.warn(f'Sensitivity function file {out_sensfile} not created. '
                          f'Deleting the folder {workrun_dir.name}.')
                shutil.rmtree(workrun_dir)


# print the list of skipped sensitivity functions
msgs.info('Skipped sensitivity functions:')
for s, r in zip(skipped, reasons):
    msgs.info(f'{s}:   {r}')


