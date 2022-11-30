#
# See top-level LICENSE.rst file for Copyright information
#
# -*- coding: utf-8 -*-
"""
Classes and utility functions for running individual pypeit tests.
"""


import os.path
import shutil
import subprocess
import datetime
import traceback
import glob
from abc import ABC, abstractmethod

from astropy.table import Table
import numpy as np

from pypeit import inputfiles

from IPython import embed

_COVERAGE_ARGS = ["--source", "pypeit", "--omit", "*PypeIt/pypeit/tests/*,*PypeIt/pypeit/deprecated/*", "--parallel-mode"] 

class PypeItTest(ABC):
    """Abstract base class for classes that run pypeit tests and hold the results from those tests."""


    def __init__(self, setup, pargs, description, log_suffix):
        """
        Constructor

        Args:
            setup (:obj:`TestSetup`) Test setup containing the test.
            description (str): A description of the type of test performed.
            log_suffix (str): The suffix to use for the log file name for the type of test.
        """

        self.setup = setup
        self.description = description
        self.log_suffix = log_suffix
        self.coverage = pargs.coverage is not None
        self.env = os.environ
        """ :obj:`Mapping`: OS Environment to run the test under."""

        self.passed = None
        """ bool: True if the test passed, False if the test failed, None if the test is in progress"""

        self.command_line = None
        """ :obj:`list` of :obj:`str`: The command and arguments that were run for the test."""

        self.error_msgs = []
        """ :obj:`list` of :obj:`str`: List of any error messages generated by the test."""

        self.logfile = None
        """ str: The log file for the test """

        self.pid = None
        """ int: The process id of the child process than ran the test"""

        self.start_time = None
        """ :obj:`datetime.datetime`: The date and time the test started."""

        self.end_time = None
        """ :obj:`datetime.datetime`: The date and time the test finished."""


    def __str__(self):
        """Return a summary of the test and the status.

           Example: 'shane_kast_red/600_7500_d57 pypeit (without masters)'
        """
        return f"{self.setup} {self.description}"

    def get_logfile(self):
        """Return a unique logifle name for the test"""
        # Get a unique log file to prevent a test from overwriting the log from a previous test
        name = '{0}_{1}.{2}.log'.format(self.setup.instr.lower(), self.setup.name.lower(), self.log_suffix)
        return get_unique_file(os.path.join(self.setup.rdxdir, name))



    @abstractmethod
    def build_command_line(self):
        pass

    def run(self):
        """Run a test in a child process."""

        try:
            # Open a log for the test
            child = None
            self.logfile = self.get_logfile()            
            self.command_line = self.build_command_line()

            with open(self.logfile, "a") as f:
                try:
                    if self.coverage:
                        # Coverage will need the full path to the script
                        full_path_to_command = shutil.which(self.command_line[0])
                        if full_path_to_command is not None:
                            self.command_line[0] = full_path_to_command
                        else:
                            raise RuntimeError(f"Could not find full path for {self.command_line[0]}")

                        self.command_line = ["coverage", "run"] + _COVERAGE_ARGS + self.command_line
                    if self.start_time is None:
                        # If a subclass sets the start time or calls run multiple times,
                        # (see deimos QL) use the first value as the start rather than overwriting it.
                        self.start_time = datetime.datetime.now()
                        
                    child = subprocess.Popen(self.command_line, stdout=f, stderr=f, env=self.env, cwd=self.setup.rdxdir)
                    self.pid = child.pid
                    child.wait()
                    self.end_time = datetime.datetime.now()
                    self.passed = (child.returncode == 0)
                finally:
                    # Kill the child if the parent script exits due to a SIGTERM or SIGINT (Ctrl+C)
                    if child is not None:
                        child.terminate()

        except Exception:
            # An exception occurred while running the test
            # Make sure it's marked as failed and document that an exception ocurred
            self.error_msgs.append(f"Exception in Test {self}:")
            self.error_msgs.append(traceback.format_exc())
            self.passed = False

        return self.passed

    def check_for_missing_files(self):
        """Return a list of any missing files the test requires. This is called before testing begins, so
        files generated during testing should be included"""
        return []


class PypeItSetupTest(PypeItTest):
    """Test subclass that runs pypeit_setup"""
    def __init__(self, setup, pargs):
        super().__init__(setup, pargs, "pypeit_setup", "setup")
        setup.generate_pyp_file = True

    def run(self):

        if super().run():
            # Check for the pypeit file after running the test
            rdxdir = os.path.join(self.setup.rdxdir, self.setup.instr.lower() + '_A')
            pyp_file = os.path.join(
                rdxdir, f'{self.setup.instr}_A.pypeit')
            if not os.path.isfile(pyp_file):
                self.error_msgs.append(f"Could not find expected pypeit file {pyp_file}")
                self.passed = False
            else:
                # If the pypeit file was created, put it's location and the new output
                # directory into the setup object for subsequent tests to use
                pyp_file = os.path.split(pyp_file)[1]

                self.setup.pyp_file = pyp_file
                self.setup.rdxdir = rdxdir

        return self.passed

    def build_command_line(self):
        return ['pypeit_setup', '-r', self.setup.rawdir, '-s',
                self.setup.instr, '-c all', '-o', '--output_path', self.setup.rdxdir]


class PypeItReduceTest(PypeItTest):
    """Test subclass that runs run_pypeit"""

    def __init__(self, setup, pargs, ignore_masters=None, std=False):

        self.ignore_masters = ignore_masters if ignore_masters is not None else pargs.do_not_reuse_masters

        description = f"pypeit {'standards ' if std else ''}{'(ignore masters)' if self.ignore_masters else ''}"
        super().__init__(setup, pargs, description, "test")

        self.std = std
        # If the pypeit file isn't being created by pypeit_setup, copy it and update it's path
        if not self.setup.generate_pyp_file:
            self.pyp_file = template_pypeit_file(self.setup.dev_path,
                                                 self.setup.instr,
                                                 self.setup.name,
                                                 self.std)

            self.pyp_file = fix_pypeit_file_directory(self.pyp_file,
                                                      self.setup.dev_path,
                                                      self.setup.rawdir,
                                                      self.setup.instr,
                                                      self.setup.name,
                                                      self.setup.rdxdir,
                                                      self.std)

    def build_command_line(self):
        if self.setup.generate_pyp_file:
            self.pyp_file = self.setup.pyp_file

        command_line = ['run_pypeit', self.pyp_file, '-o']
        if self.ignore_masters:
            command_line += ['-m']

        return command_line

    def check_for_missing_files(self):
        if not self.setup.generate_pyp_file and not os.path.isfile(self.pyp_file):
            return [self.pyp_file]
        else:
            return []

class PypeItSensFuncTest(PypeItTest):
    """Test subclass that runs pypeit_sensfunc"""
    def __init__(self, setup, pargs, std_file, sens_file=None):
        super().__init__(setup, pargs, "pypeit_sensfunc", "test_sens")
        self.std_file = std_file
        self.sens_file = sens_file

        if self.sens_file is not None:
            self.sens_file = os.path.join(setup.dev_path, 'sensfunc_files', self.sens_file)

    def run(self):

        search_pattern = os.path.join(self.setup.rdxdir, "Science", self.std_file)
        files = glob.glob(search_pattern)

        if len(files) == 0:
            self.error_msgs.append(f"Could not find std file matching!: {search_pattern}")
            self.passed = False
        elif len(files) > 1:
            self.error_msgs.append(f"Found more than one std file matching!: {search_pattern}")
            self.passed = False
        else:
            self.std_file = files[0]
            return super().run()

    def build_command_line(self):
        command_line = ['pypeit_sensfunc', self.std_file]

        if self.sens_file is not None:
            command_line += ['-s', self.sens_file]

        return command_line

    def check_for_missing_files(self):
        if self.sens_file is not None and not os.path.exists(self.sens_file):
            return [self.sens_file]
        else:
            return []


class PypeItFluxSetupTest(PypeItTest):
    """Test subclass that runs pypeit_flux_setup"""
    def __init__(self, setup, pargs):
        super().__init__(setup, pargs, "pypeit_flux_setup", "test_flux_setup")

    def build_command_line(self):
        return ['pypeit_flux_setup', os.path.join(self.setup.rdxdir, 'Science')]


class PypeItFluxTest(PypeItTest):
    """Test subclass that runs pypeit_flux_calib"""
    def __init__(self, setup, pargs):
        super().__init__(setup, pargs, "pypeit_flux", "test_flux")

        self.flux_file = os.path.join(self.setup.dev_path, 'fluxing_files',
                                      '{0}_{1}.flux'.format(self.setup.instr.lower(), self.setup.name.lower()))

    def build_command_line(self):
        return ['pypeit_flux_calib', self.flux_file]

    def check_for_missing_files(self):
        if not os.path.exists(self.flux_file):
            return [self.flux_file]
        else:
            return []

class PypeItFlexureTest(PypeItTest):
    """Test subclass that runs pypeit_deimos_flexure"""
    def __init__(self, setup, pargs):
        super().__init__(setup, pargs, "pypeit_multislit_flexure", "test_flexure")

        self.flexure_file = os.path.join(self.setup.dev_path, 'flexure_files',
                                      '{0}_{1}.flex'.format(self.setup.instr.lower(), self.setup.name.lower()))

    def build_command_line(self):
        return ['pypeit_multislit_flexure', self.flexure_file, 'testing_']

    def check_for_missing_files(self):
        if not os.path.exists(self.flexure_file):
            return [self.flexure_file]
        else:
            return []
class PypeItCoadd1DTest(PypeItTest):
    """Test subclass that runs pypeit_coadd_1dspec"""

    def __init__(self, setup, pargs):
        super().__init__(setup, pargs, "pypeit_coadd_1dspec", "test_1dcoadd")

        self.coadd_file = os.path.join(self.setup.dev_path, 'coadd1d_files',
                                       '{0}_{1}.coadd1d'.format(self.setup.instr.lower(), self.setup.name.lower()))

    def build_command_line(self):

        # Double check the object ids in the coadd file to see if they are slightly off.
        # Correct them if they are

        file_to_obj = dict()          # Mapping of files to actual objects within them
        orig_coadd1d_lines = []       # original coadd1d lines
        corrected_filenames = []  # corrected coadd1d files 
        corrected_objid = []  # corrected objids


        # Read the dev-suite coadd file, and also gather the list of files
        orig_coadd1dFile = inputfiles.Coadd1DFile.from_file(self.coadd_file)
        cfg_lines = orig_coadd1dFile.cfg_lines

        unique_files = set()
        for row in orig_coadd1dFile.data:
            (filename, object) = row['filename'], row['obj_id'] 
            orig_coadd1d_lines.append((filename, object))
            unique_files.add(filename)
        
        # Build a mapping of files to objects using the spec1d txt info file
        for unique_file in unique_files:
            txt_filename = os.path.join(self.setup.rdxdir, os.path.splitext(unique_file)[0] + ".txt")
            txt_info = Table.read(txt_filename, format='ascii.fixed_width')
            file_to_obj[unique_file] = txt_info['name']

        # Go through each original coadd1d line and correct the object names
        # if needed
        for (filename, object) in orig_coadd1d_lines:
            if object in file_to_obj[filename]:
                # No correction needed
                corrected_filenames.append(filename)
                corrected_objid.append(object)
            else:
                # Find the spatial position within object name, which
                # can start with either SPAT or OBJ
                object_parts = object.split('-')
                if object_parts[0].startswith("SPAT"):
                    object_prefix_len = 4
                    object_spat_pos = float(object_parts[0][4:])
                elif object_parts[0].startswith("OBJ"):
                    object_prefix_len = 3
                    object_spat_pos = float(object_parts[0][3:])
                else:
                    # Unrecognized format, don't try to correct
                    with open(self.logfile, "a") as f:
                        print(f"WARNING: Could not correct coadd1d file, unrecognized object id format: {object}", file=f)
                    corrected_filenames.append(filename)
                    corrected_objid.append(object)
                    continue
                
                # Build an array of distances between the actual objects
                # and the original object, and find the closest one
                distances = np.fabs(np.array([float(x.split('-')[0][object_prefix_len:]) for x in file_to_obj[filename]]) - object_spat_pos)
                sorted_dist = np.argsort(distances)
                if distances[sorted_dist][0] <=2:
                    # If the closest one is within two pixels, choose it
                    closest_obj = file_to_obj[filename][sorted_dist][0]
                    corrected_filenames.append(filename)
                    corrected_objid.append(closest_obj)
                else:
                    with open(self.logfile, "a") as f:
                        print(f"WARNING: Could not correct coadd1d file, closest object '{file_to_obj[filename][sorted_dist][0]}' is more than 2 pixels away from: {object}", file=f)
                    corrected_filenames.append(filename)
                    corrected_objid.append(object)

        
        # Generate the new coadd1d file
        final_coadd_file = get_unique_file(os.path.join(self.setup.rdxdir, f"{self.setup.instr.lower()}_{self.setup.name.lower()}_merged.coadd1d"))
        data_block = Table()
        data_block['filename'] = corrected_filenames
        data_block['obj_id'] = corrected_objid
        new_coadd1dFile = inputfiles.Coadd1DFile(config=cfg_lines,
                                             data_table=data_block)
        new_coadd1dFile.write(final_coadd_file) 

        return ['pypeit_coadd_1dspec', final_coadd_file]

    def check_for_missing_files(self):
        if not os.path.exists(self.coadd_file):
            return [self.coadd_file]
        else:
            return []

class PypeItCoadd2DTest(PypeItTest):
    """Test subclass that runs pypeit_coadd_2dspec"""
    def __init__(self, setup, pargs, coadd_file=None, obj=None):
        super().__init__(setup, pargs, "pypeit_coadd_2dspec", "test_2dcoadd")
        self.obj = obj

        if coadd_file:
            self.coadd_file = template_coadd2d_file(setup.dev_path, setup.instr, setup.name)
        else:
            self.coadd_file = None

        if self.coadd_file is None and self.obj is None:
            raise ValueError('Must provide coadd2d file or object name.')

    def build_command_line(self):
        command_line = ['pypeit_coadd_2dspec']
        command_line += ['--obj', self.obj] if self.coadd_file is None else ['--file', self.coadd_file]
        return command_line

    def check_for_missing_files(self):
        if self.coadd_file and not os.path.exists(self.coadd_file):
            return [self.coadd_file]
        else:
            return []


class PypeItTelluricTest(PypeItTest):
    """Test subclass that runs pypeit_tellfit"""

    def __init__(self, setup, pargs, coadd_file, tell_file):
        super().__init__(setup, pargs, "pypeit_tellfit", 'test_tellfit')
        self.coadd_file = coadd_file

        self.tell_file = os.path.join(self.setup.dev_path, 'tellfit_files',
                                f'{self.setup.instr.lower()}_{self.setup.name.lower()}.tell') \
                            if tell_file else None

    def build_command_line(self):
        command_line = ['pypeit_tellfit', os.path.join(self.setup.rdxdir, self.coadd_file)]
        command_line += ['-t', f'{self.tell_file}']
        return command_line

class PypeItCollate1DTest(PypeItTest):
    """Test subclass that runs pypeit_collate_1d"""
    def __init__(self, setup, pargs, files, **options):
        super().__init__(setup, pargs, "pypeit_collate", "test_collate")
        self.files = files
        self.options = options

        # Cleanup some files so tests are repeatable (collate normally appends to these files)
        report_file = os.path.join(setup.rdxdir, "collate_report.dat")
        warnings_file = os.path.join(setup.rdxdir, "collate_warnings.txt")
        if os.path.exists(report_file):
            os.remove(report_file)
        if os.path.exists(warnings_file):
            os.remove(warnings_file)


    def build_command_line(self):

        expanded_files = []
        for file_pattern in self.files:
            expanded_files += glob.glob(os.path.join(self.setup.rdxdir, file_pattern))
            
        command_line = ['pypeit_collate_1d', '--spec1d_outdir', self.setup.rdxdir, '--spec1d_files'] + expanded_files

        for option in self.options:
            if self.options[option] is None:
                command_line += [option]
            else:
                command_line += [option, str(self.options[option])]

        return command_line


class PypeItQuickLookTest(PypeItTest):
    """Test subclass that runs quick look tests.
       The specific test script run depends on the instrument type.
    """

    def __init__(self, setup, pargs, files:list,  test_name:str=None, **options):
        super().__init__(setup, pargs, "pypeit_ql", "test_ql")
        self.files = files
        self.options = options
        self.redux_dir = os.path.abspath(pargs.outputdir)
        self.pargs = pargs
        self.test_name = test_name
        # Place the masters into REDUX_DIR/QL_MASTERS directory.
        self.output_dir = os.path.join(self.redux_dir, 'QL_MASTERS')

    def build_command_line(self):

        if self.setup.instr == 'keck_mosfire':
            command_line = ['pypeit_ql_jfh_multislit', 'keck_mosfire']
            if self.pargs.quiet:
                command_line += ['--no_gui', '--writefits']
        elif self.setup.instr == 'keck_lris_red_mark4':
            command_line = ['pypeit_ql_jfh_multislit', 'keck_lris_red_mark4']
            if self.pargs.quiet:
                command_line += ['--no_gui', '--writefits']
        else:
            # Redux folder
            redux_path = os.path.join(self.redux_dir,
                    self.setup.instr, self.setup.name)
            last_folder = 'QL'
            if self.test_name is not None:
                last_folder += '_' + self.test_name
            redux_path = os.path.join(redux_path, last_folder)
                    
            command_line = [
                'pypeit_ql', self.setup.instr,
                '--redux_path', redux_path]

        if self.setup.instr in ['keck_mosfire', 'keck_lris_red_mark4']:  # TESTING USING JFH QL
            command_line += [self.setup.rawdir] + self.files
        else:
            command_line += ['--full_rawpath', self.setup.rawdir, 
                             '--rawfiles'] + self.files
        # Aditional options
        for option in self.options:
            if self.options[option] is None:
                command_line += [option]
            elif self.options[option] in ['USE_MASTERS_DIR', 
                                          'USE_CALIB_DIR']:
                idir = os.path.join(self.redux_dir,
                    self.setup.instr, self.setup.name)
                # Masters?
                if self.options[option] == 'USE_MASTERS_DIR':
                    # Crazy hack for pypeit_setup 
                    if self.setup.instr in ['shane_kast_blue'] and \
                        self.setup.name in ['600_4310_d55']:
                            idir = os.path.join(idir,
                                f'{self.setup.instr}_A')
                    idir = os.path.join(idir, 'Masters')
                # Finally
                command_line += [option, idir]
            else:
                command_line += [option, str(self.options[option])]

        return command_line

    def run(self):
        """Generate any required quick look masters before running the quick look test"""

        if self.setup.instr == 'keck_nires' or (self.setup.instr == 'keck_mosfire' and self.setup.name == 'Y_long') or \
                (self.setup.instr == 'keck_lris_red_mark4' and self.setup.name == 'long_600_10000_d680'):
            try:

                # Build the masters with the output going to a log file
                logfile = get_unique_file(os.path.join(self.setup.rdxdir, "build_ql_masters_output.log"))
                with open(logfile, "w") as log:
                    result = subprocess.run([os.path.join(self.setup.dev_path, 'build_ql_masters'),
                                             self.setup.instr, "-s", self.setup.name, '--output_dir', self.output_dir, '--redux_dir', self.redux_dir, '--force_copy'],
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

                    print(result.stdout if isinstance(result.stdout, str) else result.stdout.decode(errors='replace'), file=log)

                if result.returncode != 0:
                    self.error_msgs.append("Failed to generate QL masters.")
                    self.passed = False
                    return False
            except Exception:
                # Prevent any exceptions from escaping the "run" method
                self.error_msgs.append("Exception building QL masters:")
                self.error_msgs.append(traceback.format_exc())
                self.passed = False
                return False

        # Run the quick look test via the parent's run method, setting the environment
        # to use the newly generated masters
        self.env = os.environ.copy()
        self.env['QL_MASTERS'] = self.output_dir
        return super().run()


def pypeit_file_name(instr, setup, std=False):
    base = '{0}_{1}'.format(instr.lower(), setup.lower())
    return '{0}_std.pypeit'.format(base) if std else '{0}.pypeit'.format(base)


def template_pypeit_file(dev_path, instr, setup, std=False):
    return os.path.join(dev_path, 'pypeit_files', pypeit_file_name(instr, setup, std=std))


def coadd2d_file_name(instr, setup):
    return '{0}_{1}.coadd2d'.format(instr.lower(), setup.lower())


def template_coadd2d_file(dev_path, instr, setup):
    return os.path.join(dev_path, 'coadd2d_files', coadd2d_file_name(instr, setup))


def fix_pypeit_file_directory(pyp_file, dev_path, raw_data, instr, setup, rdxdir, std=False):
    """
    Use template pypeit file to write the pypeit file relevant to the
    exising directory structure.

    Returns:
        str: The path to the corrected pypeit file.
    """
    # Read the pypeit file
    if not os.path.isfile(pyp_file):
        raise FileNotFoundError(pyp_file)
    with open(pyp_file, 'r') as infile:
        lines = infile.readlines()

    # Replace the default path with the local one
    for kk, iline in enumerate(lines):
        if 'data read' in iline:
            old_path = lines[kk+1].strip().split(' ')[1] if 'path' in lines[kk+1] \
                            else lines[kk+1].strip()
            subdir = ''
            newdpth = ' path ' if 'path' in lines[kk+1] else ' '
            newdpth += os.path.join(raw_data, subdir)
            newdpth += '\n'
            lines[kk+1] = newdpth
        elif 'flatfield' in iline and 'pixelflat_file' in lines[kk+1]:
            newcpth = os.path.join(dev_path, 'CALIBS', os.path.split(lines[kk+1])[1])
            lines[kk+1] = '        pixelflat_file = {0}'.format(newcpth)

    # Write the pypeit file
    pyp_file = os.path.join(rdxdir, pypeit_file_name(instr, setup, std=std))
    with open(pyp_file, 'w') as ofile:
        ofile.writelines(lines)
    return pyp_file

def get_unique_file(file):
    """Ensures a file name is unique on the file system, modifying it if neccessary.

    Args:
        file (str): The full pathname of the file.

    Return:
        A unique version of the passed in file name. If the file doesn't already exist it is returned
        unchanged, otherwise a number is added to name (before the file extension) to make it unique.
    """
    file_num = 2
    (file_base, file_ext) = os.path.splitext(file)
    while os.path.exists(file):
        file = f'{file_base}.{file_num}{file_ext}'
        file_num += 1

    return file


