
import os
import glob
import warnings
import numpy
import shutil

from IPython import embed

from pypeit.pypeitsetup import PypeItSetup
from pypeit.par.util import parse_pypeit_file

def main():

    raw_dir = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos')

    setups = glob.glob(os.path.join(raw_dir, '*'))
    output_path = os.path.join(os.getcwd(), 'output')
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)

    for setup in setups:
        
        by_hand_pypeit = os.path.join(os.getenv('PYPEIT_DEV'), 'pypeit_files',
                                      'keck_deimos_{0}.pypeit'.format(
                                        os.path.split(setup)[1].lower()))
        if not os.path.isfile(by_hand_pypeit):
            warnings.warn('No by-hand pypeit file for DEIMOS setup: {0}'.format(setup))
            continue

        ps = PypeItSetup.from_file_root(setup, 'keck_deimos', output_path=output_path)
        ps.run(setup_only=True)
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        _, _, by_hand_frametypes, _, _ = parse_pypeit_file(by_hand_pypeit, file_check=False)
        _, _, auto_frametypes, _, _ = parse_pypeit_file(pypeit_files[0], file_check=False)

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

