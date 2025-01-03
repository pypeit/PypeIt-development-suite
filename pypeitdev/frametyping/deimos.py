"""
Script trolls the Keck DEIMOS raw data directories in the DevSuite and
runs PypeItSetup on them.

The frametypes are then compared to those in the pre-cooked pypeit
files.  Any differences cause the script th fault with a ValueError.

Script can be executed from anywhere.  However, note it:

    - requires the PYPEIT_DEV environmental variable that points to the
      top-level directory of the DevSuite
    - will remove and overwrite a directory called `output` within the
      current working directory.

"""

import os
import glob
import warnings
import numpy
import shutil

from IPython import embed

from pypeit.pypeitsetup import PypeItSetup
from pypeit.par.util import parse_pypeit_file

def main():

    # Raw DEIMOS directory
    raw_dir = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos')

    # Get the list of setup directories
    setups = glob.glob(os.path.join(raw_dir, '*'))

    # Set the output path and *remove if* if it already exists
    output_path = os.path.join(os.getcwd(), 'output')
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)

    # Iterate through the setups
    for setup in setups:
 
        # Find the relevant pypeit file constructed by hand.
        by_hand_pypeit = os.path.join(os.getenv('PYPEIT_DEV'), 'pypeit_files',
                                      'keck_deimos_{0}.pypeit'.format(
                                        os.path.split(setup)[1].lower()))

        if not os.path.isfile(by_hand_pypeit):
            # It doesn't exist, so assume there is no by-hand pypeit
            # file to compare to
            warnings.warn('No by-hand pypeit file for DEIMOS setup: {0}'.format(setup))
            continue

        # Run pypeit_setup
        ps = PypeItSetup.from_file_root(setup, 'keck_deimos', output_path=output_path)
        ps.run(setup_only=True)
        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        # Read the frame types from both the by-hand and automated
        # pypeit files
        _, _, by_hand_frametypes, _, _ = parse_pypeit_file(by_hand_pypeit, file_check=False)
        _, _, auto_frametypes, _, _ = parse_pypeit_file(pypeit_files[0], file_check=False)

        # For each file in the by-hand list, check that the frame types
        # in the automatically generated pypeit file are identical
        for f in by_hand_frametypes.keys():
            type_list = numpy.sort(by_hand_frametypes[f].split(','))
            if 'science' in type_list or 'standard' in type_list:
                # Only ensuring that calibrations are correctly typed
                continue
            if f not in auto_frametypes.keys():
                raise KeyError('Frame {0} not automatically parsed for setup {1}.'.format(f, setup))
            if not numpy.array_equal(type_list, numpy.sort(auto_frametypes[f].split(','))):
                raise ValueError('Frame types differ for file {0} in setup {1}\n'.format(f, setup)
                                 + '    By-hand types: {0}'.format(by_hand_frametypes[f])
                                 + '    Automated types: {0}'.format(auto_frametypes[f]))

        # Clean up after every setup
        shutil.rmtree(output_path)


if __name__ == '__main__':
    main()

