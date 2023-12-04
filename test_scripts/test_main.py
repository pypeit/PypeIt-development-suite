#!/usr/bin/env python3
#
# See top-level LICENSE.rst file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script runs the PypeIt development suite of tests
"""

import sys
import os
import os.path
import subprocess
from queue import PriorityQueue, Empty
from threading import Thread, Lock
import traceback
import datetime
from pathlib import Path
import textwrap

import numpy as np
import pypeit 
# Stop logging from pypeit.par.utils when reading/writing coadd1d files
pypeit.msgs.reset(verbosity=0) 


from .test_setups import TestPhase, all_tests, all_setups
from .pypeit_tests import get_unique_file, _COVERAGE_ARGS

test_run_queue = PriorityQueue()
""":obj:`queue.Queue`: FIFO queue for test setups to be run."""

class TestPriorityList(object):
    """A class for reading and updating the order test setups are are tested.
    
    The order that test setups should be tested in is stored as a list of instr/setup keys in a file named 
    'test_priority_list' in $PYPEIT_DEV. This class is responsible for reading, writing and updating the file, and
    for translating the order into a priority that can be used in a PriorityQueue.
    
    To determine the order to run tests, the slowest tests are run first in order to optimize CPU efficiency when
    running with multiple threads.  Once the test suite is completed, the order is recalculated based on how long
    the tests actually took.  If the order is changed, the new ordering is written to the test_priority_list file.
    
    Attributes:
        _prority_map (:obj:`dict` of str to int): Maps instr/setup name identifiers to an integer priority.

        _updated (boolean): Whether the _priority_map has been modified and should be written to disk.
        _file (str):        The file name to write the priority list to.
    """

    def __init__(self, file):
        """Reads the priority list from a file."""
        self._priority_map = dict()
        self._updated = False
        self._file = file

        count = 0
        if os.path.exists(file):
            with open(file, "r") as f:
                while True:
                    line = f.readline()
                    if len(line) == 0:
                        break
                    setup_key = line.strip()
                    if len(setup_key) > 0:
                        self._priority_map[setup_key] = count
                        count += 1

    def __len__(self):
        """Return how many test setups are in the priority list"""
        return len(self._priority_map)

    def set_test_setup_priority(self, setup):
        """Set the test priority of a TestSetup object."""

        if setup.key in self._priority_map:
            setup.priority = self._priority_map[setup.key]
        else:
            # If the test setup is currently known, assume it is short and put it at the end of the priority list
            # by assigning it a large number.
            setup.priority = sys.maxsize

    def update_priorities(self, setups):
        """Update the test priorities based on the actual runtimes of the test setups"""
        setup_durations = []
        for setup in setups:
            duration = sum([test.end_time - test.start_time for test in setup.tests], datetime.timedelta())
            setup_durations.append((setup.key, duration))

        count = 0
        new_priority_map = dict()
        for setup_key in [item[0] for item in sorted(setup_durations, key=lambda x: x[1])]:
            new_priority_map[setup_key] = count
            count += 1

        if new_priority_map != self._priority_map:
            self._priority_map = new_priority_map
            self._updated = True

    def write(self):
        """Write the test priority list to a file if it has changed."""
        if self._updated:
            with open(self._file, "w") as f:
                for item in sorted(self._priority_map.items(), key=lambda x: x[1], reverse=True):
                    print(item[0], file=f)

            self.updated = False


class TestSetup(object):
    """Representation of a test setup within the pypeit development suite.

    Attributes:
        instr (str):        The instrument for the test setup
        name (str):         The name of the test setup
        key (str):          The "instrument/setup name" key that identifies this setup in the all_test data structure
                            in test_setups.py.
        rawdir (str):       The directory with the raw data for the test setup
        rdxdir (str):       The output directory for the test setup. This can be changed as tests are run.
        dev_path (str):     The path of the Pypeit-development-suite repository
        pyp_file (str):     The .pypeit file used for the test. This may be created by a PypeItSetupTest.
        std_pyp_file (str): The standards .pypeit file used for some tests.
        priority (int):     The priority of the TestSetup. Used in a PriorityQueue to determine the order used to
                            run test setups.

        generate_pyp_file (boolean): Set to true if this setup will generate it's own .pypeit file wint pypeit_setup.

        tests (:obj:`list` of :obj:`PypeItTest`): The list of tests to run in this test setup. The tests will be run
                                                  in sequential order

        missing_files (:obj:`list` of str): List of missing files preventing the test setup from running.

    """
    def __init__(self, instr, name, rawdir, rdxdir, dev_path):
        self.instr = instr
        self.name = name
        self.key = f'{self.instr}/{self.name}'
        self.rawdir = rawdir
        self.rdxdir = rdxdir
        self.dev_path = dev_path
        self.pyp_file = None
        self.std_pyp_file = None
        self.priority = 0
        self.generate_pyp_file = False
        self.tests = []
        self.missing_files = []

    def __str__(self):
        """Return a string representation of this setup of the format "instr/name"""""
        return self.key

    def __lt__(self, other):
        """
        Compare this test setup with another using the "priority" attribute. This is included so that a
        TestSetup object can be placed into a PriorityQueue
        """
        return self.priority < other.priority


def red_text(text):
    """Utiltiy method to wrap text in the escape sequences to make it appear red in a terminal"""
    return f'\x1B[1;31m{text}\x1B[0m'

def green_text(text):
    """Utiltiy method to wrap text in the escape sequences to make it appear green in a terminal"""
    return f'\x1B[1;32m{text}\x1B[0m'


class TestReport(object):
    """Class for reporting on the status and results of testing.

    The test_started, test_skipped, and test_completed methods of this class are called during testing to keep track of
    status and provide incremental status updates to the console as tests are running.

    detailed_report and summary report are called after testing to report on the results.

    Attributes:

    pargs (:obj:`Namespace`:) The parsed command line arguments to pypeit_test
    setups (:obj:`list` of :obj:`TestSetup`):  The test setups that are being tested.

    start_time (:obj:`datetime.datetime`): The date and time testing started.
    end_time (:obj:`datetime.datetime`):   The date and time testing finished.

    num_tests (int):   The total number of tests in all of the test setups.
    num_passed (int):  The number of tests that have passed.
    num_failed (int):  The number of tests that have failed.
    num_skipped (int): The number of tests that were skipped because they depended on the results of a failed tests.
    num_active (int):  The number of tests that are currently in progress.

    failed_tests (:obj:`list` of str):  List of names of tests that have failed
    skipped_tests (:obj:`list` of str): List of names of tests that have been skipped

    testing_complete (bool): Whether testing has completed.
    lock (:obj:`threading.Lock`): Lock used to synchronize access when multiple threads are reporting status. This
                                  prevents scrambled output being sent to stdout.
    """
    def __init__(self, pargs):
        self.pargs = pargs
        self.test_setups = []
        self.end_time = None

        self.num_tests = 0
        self.num_passed = 0
        self.num_failed = 0
        self.num_skipped = 0
        self.num_active = 0
        self.failed_tests = []
        self.skipped_tests = []
        self.lock = Lock()
        self.testing_complete = False
        self.start_time = datetime.datetime.now()

        self.pytest_results=dict()

        if pargs.report is not None and os.path.exists(pargs.report):
            # Remove any old report files if we've been asked to overwrite it
            if not pargs.quiet:
                print(f"Overwriting existing report {pargs.report}",flush=True)
            os.unlink(pargs.report)


    def _get_test_counts(self):
        """Helper method to create a string with the current test counts"""
        verbose_info = f'{self.num_active:2} active/' if self.pargs.verbose else ''
        return f'{verbose_info}{self.num_passed:2} passed/{self.num_failed:2} failed/{self.num_skipped:2} skipped'

    def setup_testing_started(self,setups):
        """Called once test setup testing has started.
        
        Args:
            setups (list of :obj:`TestSetup`):
                The list of test setups being run. These are used to generate the detailed report once
                testing is complete.

        """
        with self.lock:
            self.test_setups = setups
            # Create the report file (if needed) and write the header to it
            if self.pargs.report:
                try:
                    with open(self.pargs.report, "a") as report_file:
                        self.detailed_report_header(output=report_file)

                except Exception as e:
                    print(f"Could not open report file {self.pargs.report}", file=sys.stderr)
                    traceback.print_exc()
                    sys.exit(1)



    def test_started(self, test):
        """Called when a test has started executing"""
        with self.lock:
            self.num_tests += 1
            self.num_active += 1


            if not self.pargs.quiet:
                verbose_info = ''
                if self.pargs.verbose:
                    verbose_info = f' at {datetime.datetime.now().ctime()}'

                print(f'{self._get_test_counts()} STARTED {test}{verbose_info}', flush=True)

    def test_skipped(self, test):
        """Called when a test has been skipped because a test before it has failed"""
        with self.lock:
            self.num_skipped += 1
            self.skipped_tests.append(test)

            if not self.pargs.quiet:
                print(f'{self._get_test_counts()} {red_text("SKIPPED")} {test}', flush=True)

    def test_completed(self, test):
        """Called when a test has finished executing."""
        with self.lock:
            self.num_active -= 1
            if test.passed:
                self.num_passed += 1
            else:
                self.num_failed += 1
                self.failed_tests.append(test)

            if not self.pargs.quiet:
                verbose_info = ''
                if self.pargs.verbose:
                    if test.end_time is not None and test.start_time is not None:
                        duration = test.end_time-test.start_time
                    else:
                        duration = 'n/a'

                    verbose_info = f' with pid {test.pid} at {datetime.datetime.now().ctime()} Duration {duration}'

                if test.passed:
                    print(f'{self._get_test_counts()} {green_text("PASSED")}  {test}{verbose_info}', flush=True)
                else:
                    print(f'{self._get_test_counts()} {red_text("FAILED")}  {test}{verbose_info}', flush=True)
                    self.report_on_test(test, flush=True)

    def test_setup_completed(self, test_setup):
        """Called once all of the tests in a test setup have completed"""
        if self.pargs.report is not None:
            with self.lock:
                with open(self.pargs.report, "a") as report_file:
                    self.report_on_setup(test_setup, report_file)            

    def testing_completed(self):
        """Called once all test setups have complete"""
        self.end_time = datetime.datetime.now()
        if self.pargs.report is not None:
            with open(self.pargs.report, "a") as report_file:
                self.summary_report(report_file)

    def pytest_started(self, test_descr):
        """Called when a set of pytest tests have started.
        
        Args:
            test_descr (str): A short description of the pytest test suite that
                              will be displayed in the test report and uniquely
                              identify the test suite. For example:
                              "Unit Tests".
        """

        if not self.pargs.quiet:
            print(f"Running {test_descr}", flush=True)

        if self.pargs.report is not None:
            with open(self.pargs.report, "a") as report_file:
                print(f"{test_descr} Results:", file=report_file)
                print ("-------------------------", file=report_file)                

    def pytest_line(self, test_descr, line):
        """Called for each line ouptut from a pytest run. Each line is echoed to
        stdout and, if requested, a report file.
        
        Args:
            test_descr (str): A short description of the pytest test suite that
                              will be displayed in the test report and uniquely
                              identify the test suite. For example:
                              "Unit Tests".

            line (str):       A line from the stdout of pytest.
        
        """
        if not self.pargs.quiet:
            print(line, flush=True)

        if self.pargs.report is not None:
            with open(self.pargs.report, "a") as report_file:
                print(line, file=report_file)
        
        # Save any summary lines found for reporting later.
        if "warnings" in line or "passed" in line or "failed" in line:
            self.pytest_results[test_descr] = line.replace("=", "")


    def detailed_report(self, output=sys.stdout):
        """Display a detailed report on testing to the given output stream"""

        self.detailed_report_header(output)

        for setup in self.test_setups:
            self.report_on_setup(setup, output)
        print ("-------------------------", file=output)
        self.summary_report(output)

    def detailed_report_header(self, output):
        """Display the header information of a detaile report"""

        # Report header information
        print('Reduced data for the following setups:', file=output)
        for setup in self.test_setups:
            print(f'    {setup}', file=output)
        print('', file=output)

        if self.pargs.threads > 1:
            print(f'Ran tests in {self.pargs.threads} parallel processes\n', file=output)

    def summarize_setup_tests(self, output=sys.stdout):
        """Display a summary of the PypeIt setup tests"""

        calib_text = '(Existing calibrations ignored)' if self.pargs.do_not_reuse_calibs else ''
        if self.num_tests == self.num_passed:
            print("\x1B[" + "1;32m" +
                  "--- PYPEIT DEVELOPMENT SUITE PASSED {0}/{1} TESTS {2} ---".format(
                      self.num_passed, self.num_tests, calib_text)
                  + "\x1B[" + "0m" + "\r", file=output)
        else:
            print("\x1B[" + "1;31m" +
                  "--- PYPEIT DEVELOPMENT SUITE FAILED {0}/{1} TESTS {2} ---".format(
                      self.num_failed, self.num_tests, calib_text)
                  + "\x1B[" + "0m" + "\r", file=output)
            print('Failed tests:', file=output)
            for t in self.failed_tests:
                print('    {0}'.format(t), file=output)
            print('Skipped tests:', file=output)
            for t in self.skipped_tests:
                print('    {0}'.format(t), file=output)

    def summarize_pytest_results(self, test_descr, output=sys.stdout):
        """Display a summary of a pytest run."""
        if test_descr not in self.pytest_results:
            # The tests weren't run, so no results
            return
        results = self.pytest_results[test_descr]
        if "failed" in results:
            print("\x1B[" + "1;31m" + f"--- PYTEST {test_descr.upper()} FAILED " + "\x1B[" + "0m"
                  + results +  "\x1B[" + "1;32m" + "---" + "\x1B[" + "0m" + "\r", file=output)
        else:
            print("\x1B[" + "1;32m" + f"--- PYTEST {test_descr.upper()} PASSED " + "\x1B[" + "0m"
                  + results +  "\x1B[" + "1;32m" + "---" + "\x1B[" + "0m" + "\r", file=output)

    def performance_results(self, output):
        """Display performance statistics on PypeIt tests."""
        print("Setup,Test Type,Start Time,End Time,Duration(s),Memory Usage (bytes),Duration (D:H:M:S), Memory Usage (MiB)", file=output)
        for setup in self.test_setups:
            for test in setup.tests:
                if test.start_time is not None and test.end_time is not None:
                    duration = test.end_time - test.start_time
                    duration_secs = duration.total_seconds()
                else:
                    duration = ""
                    duration_secs = ""

                if test.max_mem is None:
                    mem_usage = ""
                    mem_usage_megs = ""
                else:
                    mem_usage = test.max_mem
                    mem_usage_megs = test.max_mem / (2**20)

                print(f'{test.setup},{test.description},{test.start_time},{test.end_time},{duration_secs},{mem_usage},{duration},{mem_usage_megs}', file=output)


    def print_tail(self, file, num_lines, output=sys.stdout, flush=False):
        """Print the last num_lines of a file."""
        result = subprocess.run(['tail', f'-{num_lines}', file], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(result.stdout.decode(), file=output, flush=flush)


    def summary_report(self, output=sys.stdout):
        """Display a summary report on the results of testing to the given output stream"""

        print ("\nTest Summary\n--------------------------------------------------------", file=output)
        self.summarize_pytest_results("PypeIt Unit Tests", output)
        self.summarize_pytest_results("Unit Tests", output)
        self.summarize_pytest_results("Vet Tests", output)
        self.summarize_setup_tests(output)

        if self.pargs.coverage is not None:
            print(f"Coverage results:", file=output)
            self.print_tail(self.pargs.coverage, 1, output)

        print(f"Testing Started at {self.start_time.isoformat()}", file=output)
        print(f"Testing Completed at {self.end_time.isoformat()}", file=output)
        print(f"Total Time: {self.end_time - self.start_time}", file=output)


    def report_on_test(self, test, output=sys.stdout, flush=False):
        """Print a detailed report on the status of a test to the given output stream."""

        if test.passed:
            result = green_text('--- PASSED')
        elif test.passed is None:
            result = red_text('--- SKIPPED')
        else:
            result = red_text('--- FAILED')

        if test.start_time is not None and test.end_time is not None:
            duration = test.end_time - test.start_time
        else:
            duration = None

        print("----", file=output, flush=flush)
        print(f"{test.setup} {test.description} Result: {result}\n", file=output, flush=flush)
        print(f'Logfile:    {test.logfile}', file=output, flush=flush)
        print(f'Process Id: {test.pid}', file=output, flush=flush)
        print(f'Start time: {test.start_time.ctime() if test.start_time is not None else "n/a"}', file=output, flush=flush)
        print(f'End time:   {test.end_time.ctime() if test.end_time is not None else "n/a"}', file=output, flush=flush)
        print(f'Duration:   {duration}', file=output, flush=flush)
        print(f'Mem Usage:  {test.max_mem}', file=output, flush=flush)
        print(f"Command:    {' '.join(test.command_line) if test.command_line is not None else ''}", file=output, flush=flush)
        print('', file=output, flush=flush)
        print('Error Messages:', file=output, flush=flush)

        for msg in test.error_msgs:
            print(msg, file=output, flush=flush)

        print('', file=output, flush=flush)
        print("End of Log:", file=output, flush=flush)
        if test.logfile is not None and os.path.exists(test.logfile):
            self.print_tail(test.logfile, 3, output, flush)

        print('\n', file=output, flush=flush)

    def report_on_setup(self, setup, output=sys.stdout):
        """Print a detailed report on the status of a test setup and the tests within it to the given output stream."""

        print ("-------------------------", file=output)
        print (f"Test Setup: {setup}\n", file=output)
        print ("-------------------------", file=output)
        print("Directories:", file=output)
        print(f"         Raw data: {setup.rawdir}", file=output)
        print(f"    PypeIt output: {setup.rdxdir}", file=output)
        print("Files:", file=output)
        print(f"     .pypeit file: {setup.pyp_file}", file=output)
        print(f" Std .pypeit file: {setup.std_pyp_file}", file=output)

        print("Tests:", file=output)

        for t in setup.tests:
            self.report_on_test(t, output)

def clear_coverage_data(redux_out):
    """Clear any leftover coverage data that may be left over form an interrupted prior run."""
    path = Path(redux_out)
    for file in path.rglob(".coverage*"):
        file.unlink(missing_ok = True)

def run_pytest(pargs, test_descr, test_dir, test_report, 
               redux_out=None):
    """Run pytest on a directory of test files.
    
    Args:
        pargs (:obj:`argparse.Namespace`): The arguments to pypeit_test, as returned by argparse.

        test_descr (str): 
            A short description of the pytest test suite that will be displayed in 
            the test report and uniquely identify the test suite. For example:
            "Unit Tests".

        test_dir (str): 
            The directory where the pytest tests to run are stored.

        test_report (:obj:`TestReport`): 
            The test report object for this test run. Output from the pytest run
            will be passed to this object.

        redux_out (str):
            The location of the output of this dev-suite run. Optional, only required
            if the pytest suite requires the output from the dev-suite (i.e vet_tests).
    """
    abs_test_dir = os.path.abspath(test_dir)

    test_report.pytest_started(test_descr)

    # Tests not written to run in parallel yet
    #args = ["pytest", "-v", "--color=yes", "-n", str(pargs.threads)]

    # Run pytest using coverage if requested
    if pargs.coverage is not None:
        args = ["coverage", "run"] + _COVERAGE_ARGS + ["-m", "pytest", "-v", "--color=yes"]
    else:
        args = ["pytest", "-v", "--color=yes"]

    if not pargs.show_warnings:
        args.append("--disable-warnings")
    
    if redux_out is not None:
        args += ["--redux_out", redux_out]

    args.append(abs_test_dir)

    # Run pytest, sending the outpu to the test report.
    # We change the current directory so that the coverage output goes to the outputdir
    with subprocess.Popen(args,stderr=subprocess.STDOUT, stdout=subprocess.PIPE,cwd=pargs.outputdir) as p:
        while(p.poll() is None):
            test_report.pytest_line(test_descr, p.stdout.readline().decode().strip())

def generate_coverage_report(pargs):

    # Find the coverage files
    coverage_files = [str(path) for path in Path(pargs.outputdir).rglob(".coverage.*")]
    
    if len(coverage_files) > 0:
        if not pargs.quiet:
            print("Combining coverage files...", flush=True)
    else:
        with open(pargs.coverage, "w") as f:
            print("Couldn't find coverage files to combine.", file=f)
        return

    # Combine them, changing to the output dir to keep the coverage 
    # data files there
    process = subprocess.run(["coverage", "combine"] + coverage_files, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=pargs.outputdir)
    if process.returncode != 0:
        if not pargs.quiet:
            print("Failed to combine coverage data.", flush=True)
        with open(pargs.coverage, "w") as f:
            print("Failed to combine coverage files. Output:", file=f)
            print(process.stdout, file=f)
        return

    # Generate the report.
    with open(pargs.coverage, "w") as f:
        process = subprocess.run(["coverage", "report", "-m"], stdout=f, stderr=subprocess.STDOUT, cwd=pargs.outputdir)

def raw_data_dir():
    return os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA')


def available_data():
    return [x.name for x in Path(raw_data_dir()).glob('*') if x.is_dir()]

def available_setups(raw_data, instr):
    """Return the setups available in RAW_DATA directory for a given instrument."""
    return [x.name for x in Path(raw_data).joinpath(instr).glob('*') if x.is_dir()]    


def parser(options=None):
    import argparse

    dirs = available_data()
    all_tests = np.unique([ [d.split('_')[0], d.split('_')[1]] for d in dirs ])

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Run pypeit tests on a set of instruments.  '
                                                 'Typical call for testing pypeit when developing '
                                                 'new code is `./pypeit_test all`.  Execution '
                                                 'requires you to have a PYPEIT_DEV environmental '
                                                 'variable, pointing to the top-level directory '
                                                 'of the dev-suite repository (typically the '
                                                 'location of this script).  Raw data for testing '
                                                 'is expected to be at ${PYPEIT_DEV}/RAW_DATA.  '
                                                 'To run all tests for the supported instruments, '
                                                 'use \'all\'.  To only run the basic '
                                                 'reductions, use \'reduce\'.  To only run the '
                                                 'tests that use the results of the reductions, '
                                                 'use \'afterburn\'\'. Use \'list\' to view all '
                                                 'supported setups.')

    parser.add_argument('tests', type=str, nargs='+', default=None,
                        help='Which test types to run. Options are:  '
                             'pypeit_tests, unit, reduce, afterburn, ql, vet, or all. Use list to show all supported instruments and setups.')
    parser.add_argument('-o', '--outputdir', type=str, default='REDUX_OUT',
                        help='Output folder.')
    parser.add_argument('-i', '--instruments', type=str, nargs='+', 
                        help='One or more instruments to run tests for. Use "pypeit_test list" to see all supported instruments.')
    parser.add_argument('-s', '--setups', type=str, nargs='+', 
                        help='One or more setups to run tests for. Use "pypeit_test list" to see all supported setups.')
    parser.add_argument('--debug', default=False, action='store_true',
                        help='Debug using only blue setups')
    parser.add_argument('-p', '--prep_only', default=False, action='store_true',
                        help='Only prepare to execute run_pypeit, but do not actually run it.')
    parser.add_argument('-m', '--do_not_reuse_calibs', default=False, action='store_true',
                        help='run pypeit without using any existing processed calibration frames')
    parser.add_argument('-t', '--threads', default=1, type=int,
                        help='Run THREADS number of parallel tests.')
    parser.add_argument('-q', '--quiet', default=False, action='store_true',
                        help='Supress all output to stdout. If -r is not a given, a report file will be '
                             'written to <outputdir>/pypeit_test_results.txt')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Output additional detailed information while running the tests and output a '
                             'detailed report at the end of testing. This has no effect if -q is given')
    parser.add_argument('--coverage', default=None, type=str, 
                        help='Collect code coverage information. and write it to the given file.')
    parser.add_argument('-r', '--report', default=None, type=str,
                        help='Write a detailed test report to REPORT.')
    parser.add_argument('-c', '--csv', default=None, type=str,
                        help='Write performance numbers to a CSV file.')
    parser.add_argument('-w', '--show_warnings', default=False, action='store_true',
                        help='Show warnings when running unit tests and vet tests.')
    return parser.parse_args() if options is None else parser.parse_args(options)

def show_setup_list():
    """Show a list of all instruments and the setups they support."""
    print("All instruments and test setups supported by the dev-suite:\n")
    for instrument in all_setups.keys():
        print(instrument)
        # Print an indented line wrapped list of setups beneath the instrument.
        # "break_long_words=False" prevents it from breaking up setup names at an underscore,
        # which looks bad.
        setups = " ".join(all_setups[instrument])
        for line in textwrap.wrap(setups, width=80, initial_indent = "    ", 
                                  subsequent_indent="    ", break_long_words=False):
            print(line)

def thread_target(test_report):
    """Thread target method for running tests."""
    while not test_report.testing_complete:
        try:
            test_setup = test_run_queue.get(timeout=2)
        except Empty:
            # queue is currently empty but another thread may put another test on it, so try again
            continue

        passed = True
        for test in test_setup.tests:

            if not passed:
                test_report.test_skipped(test)
            else:
                test_report.test_started(test)
                passed = test.run()
                test_report.test_completed(test)

        test_report.test_setup_completed(test_setup)

        # Count the test setup as done. This needs to be done to allow the join() call in main to return when
        # all of the tests have been completed
        test_run_queue.task_done()


def main():

    # ---------------------------------------------------------------------------
    # Parse command line arguments

    pargs = parser()

    if 'list' in pargs.tests:
        show_setup_list()
        return 0

    if pargs.threads <=0:
        raise ValueError("Number of threads must be >= 1")
    elif pargs.threads > 1:
        # Set the OMP_NUM_THREADS to 1 to prevent numpy multithreading from competing for resources
        # with the multiple processes started by this script
        os.environ['OMP_NUM_THREADS'] = '1'

    raw_data = raw_data_dir()
    if not os.path.isdir(raw_data):
        raise NotADirectoryError('No directory: {0}'.format(raw_data))

    if not os.path.exists(pargs.outputdir):
        os.mkdir(pargs.outputdir)
    pargs.outputdir = os.path.abspath(pargs.outputdir)

    # Make sure we can create a writable report file (if needed) before starting the tests
    if pargs.report is None and pargs.quiet:
        # If there's no report file specified in the command line, but we're in quiet mode,
        # make up a report file name
        pargs.report = get_unique_file(os.path.join(pargs.outputdir, "pypeit_test_results.txt"))

    # ---------------------------------------------------------------------------
    # Determine which tests to run

    flg_pypeit_tests = False
    flg_unit = False
    flg_reduce = False
    flg_after = False
    flg_ql = False
    flg_vet = False

    write_priorities = False

    for test in pargs.tests:
        if test == "all":
            flg_pypeit_tests = True
            flg_unit = True
            flg_reduce = True
            flg_after = True
            flg_ql = True
            flg_vet = True

            # Write the test priority file if all tests are being run
            if pargs.instruments is None and pargs.setups is None and pargs.debug == False:
                write_priorities = True

        elif test == "pypeit_tests":
            flg_pypeit_tests = True
        elif test == "unit":
            flg_unit = True
        elif test == "reduce":
            flg_reduce = True
        elif test == "after" or test=="afterburn":
            flg_after = True
        elif test == "ql":
            flg_ql = True
        elif test == "vet":
            flg_vet = True
        else:
            print("\x1B[" + "1;31m" + "\nERROR - " + "\x1B[" + "0m" +
                  "Invalid test selected: {}\n\n".format(test) +
                  "Consult the help (pypeit_test -h)")
            return 1
            

    # ---------------------------------------------------------------------------
    # Determine which instruments will be tested


    unsupported = []
    instruments = []
    all_instruments = all_setups.keys()
    if pargs.instruments is not None and len(pargs.instruments) > 0:
        for instr in pargs.instruments:
            if instr in all_instruments:
                instruments.append(instr) 
            else:
                unsupported.append(instr)

    # Setups may be specified with a "instr/setup" syntax, parse those out
    # and make sure the instruments are included
    argument_setup_names = []
    if pargs.setups is not None and len(pargs.setups) > 0:
        for setup in pargs.setups:
            if "/" in setup:
                (instr, setup_name) = setup.split("/")
                if instr not in instruments and instr in all_instruments:
                    instruments.append(instr)
                argument_setup_names.append(setup_name)
            else:
                argument_setup_names.append(setup)

    if len(unsupported) > 0:
        print("\x1B[" + "1;33m" + "\nWARNING - " + "\x1B[" + "0m" +
                "The following instruments are not supported: {0}\n\n".format(
                unsupported))
        return 1

    # If no instruments were supplied by either the instruments or
    # setups arguments, test all instruments
    if len(instruments) == 0:
        instruments = all_instruments

    # Report
    if not pargs.quiet:
        if "all" in pargs.tests:
            print('Running all tests.')
        else:
            if flg_pypeit_tests:
                print('Running unit tests in pypeit/tests')
            if flg_unit:
                print('Running dev suite unit tests')
            if flg_reduce:
                print('Running reduce tests.')
            if flg_after:
                print('Running afterburner tests')
            if flg_ql is True:
                print('Running quick look tests')
            if flg_vet is True:
                print('Running vet tests')

    # Clean up prior coverage results that could be left over from an
    # interrupted dev suite run
    if pargs.coverage is not None:
        clear_coverage_data(pargs.outputdir)
 
    # Start Unit Tests
    test_report = TestReport(pargs)

    # For coverage testing, run the PypeIt unit tests too
    if flg_pypeit_tests and not pargs.prep_only:
        pypeit_tests_dir = Path(pypeit.__file__).parent.joinpath("tests")
        run_pytest(pargs, "PypeIt Unit Tests", str(pypeit_tests_dir), test_report)

    dev_path = os.getenv('PYPEIT_DEV')
    if flg_unit is True and not pargs.prep_only:
        run_pytest(pargs, "Unit Tests", os.path.join(dev_path, "unit_tests"), test_report)


    if flg_reduce or flg_after or flg_ql:
        # ---------------------------------------------------------------------------
        # Build the TestSetup and PypeItTest objects for testing

        # Load test setup priority from file
        priority_list = TestPriorityList('test_priority_list')
        if not pargs.quiet and pargs.verbose:
            print(f'Loaded {len(priority_list)} setup priorities')

        # Report on instruments
        if not pargs.quiet:
            print('Running tests on the following instruments:')
            for instr in instruments:
                print('    {0}'.format(instr))
            print('')


        setups = []
        missing_files = []
        for instr in instruments:
            # Only do blue instruments
            if pargs.debug and instr != 'shane_kast_blue':
                continue

            # Setups        
            if len(argument_setup_names) > 0:
                setup_names = [name for name in argument_setup_names if name in all_setups[instr]]

                # No setups for this instrument specified, so run all setups
                if len(setup_names)==0:
                    setup_names = all_setups[instr]                    
            elif pargs.debug:
                setup_names = ['600_4310_d55']
            else:
                setup_names = all_setups[instr]

            # Build test setups, check for missing files, and run any prep work
            for setup_name in setup_names:

                setup = build_test_setup(pargs, instr, setup_name, flg_reduce, flg_after,
                                        flg_ql)
                missing_files += setup.missing_files

                # set setup priority from file
                priority_list.set_test_setup_priority(setup)

                setups.append(setup)

            print('Reducing data from {0} for the following setups:'.format(instr))
            for name in setup_names:
                print('    {0}'.format(name))
            print('')

        # ---------------------------------------------------------------------------
        # Check all the data and relevant files exist before starting!
        if len(missing_files) > 0:
            raise ValueError('Missing the following files:\n    {0}'.format(
                            '\n    '.join(missing_files)))


        # ---------------------------------------------------------------------------
        # Run the tests
        test_report.setup_testing_started(setups)
        # Add tests to the test_run_queue
        for setup in setups:
            if len(setup.tests) == 0:
                continue
            test_run_queue.put(setup)

        # Start threads to run the tests
        if not pargs.quiet and pargs.threads > 1:
            print(f'Running tests in {pargs.threads} parallel processes')

        thread_pool = []
        for i in range(pargs.threads):
            new_thread = Thread(target=thread_target, args=[test_report])
            thread_pool.append(new_thread)
            new_thread.start()

        # Wait for the tests to finish
        test_run_queue.join()

        # Set the test status to complete and then wait for the threads to finish.
        # We don't run the threads as daemon threads so that main() can be called multiple times
        # in unit tests
        test_report.testing_complete = True
        for thread in thread_pool:
            thread.join()

        if not pargs.quiet:
            test_report.summarize_setup_tests()

    # Run the vet tests
    if flg_vet is True:
        run_pytest(pargs, "Vet Tests", os.path.join(dev_path, "vet_tests"), test_report, redux_out=pargs.outputdir)


    # ---------------------------------------------------------------------------
    # Build test priority list for next time, but only if all tests succeeded
    # and all tests were being run
    if write_priorities and test_report.num_passed == test_report.num_tests:
        priority_list.update_priorities(setups)
        priority_list.write()
        if not pargs.quiet and pargs.verbose:
            print(f'Wrote {len(priority_list)} setup priorities')

    if pargs.coverage is not None:
        generate_coverage_report(pargs)

    # ---------------------------------------------------------------------------
    # Finish up the report on the test results
    test_report.testing_completed()

    if pargs.csv is not None:
        with open(pargs.csv, "w") as f:
            test_report.performance_results(f)

    if not pargs.quiet:
        if pargs.verbose:
            test_report.detailed_report()
        else:
            test_report.summary_report()

    return test_report.num_failed


def build_test_setup(pargs, instr, setup_name, flg_reduce, flg_after, flg_ql):
    """
    Builds a TestSetup object including the tests that it will run

    Args:
        pargs (:obj:`argparse.Namespace`): 
            The arguments to pypeit_test, as returned by argparse.
        
        instr (str): 
            The instrument the setup is for.

        setup_name (str): 
            The name of the test setup.

        flg_reduce (bool): 
            Whether or not reduce tests are being run.

        flg_after (bool): 
            Whether or not afterburner tests are being run.

        flg_ql (bool): 
            Whether or not quick look tests are being run.

    Returns:
        :obj:`TestSetup`:
            A TestSetup object representing the test setup.
    """

    dev_path = os.getenv('PYPEIT_DEV')
    raw_data = raw_data_dir()

    # Directory with raw data
    rawdir = os.path.join(raw_data, instr, setup_name)

    # Directory for reduced data
    rdxdir = os.path.join(pargs.outputdir, instr, setup_name)
    if not os.path.exists(rdxdir):
        # Make the directory
        os.makedirs(rdxdir)

    # Create the test setup and set it's priority
    setup = TestSetup(instr, setup_name, rawdir, rdxdir, dev_path)

    # Go through each test type and add it to this setup if it's applicable and
    # selected by the command line arguments
    for test_descr in all_tests:

        # Check instruments
        if setup.instr not in test_descr['setups']:
            continue
        # Check setup
        if setup_name not in test_descr['setups'][setup.instr]:
            continue

        for kwargs in test_descr['setups'][setup.instr][setup_name]:
            # Create the test, this will also run any prep_only steps in the
            # __init__ method
            try:
                test = test_descr['factory'](setup, pargs, **kwargs)
            except FileNotFoundError as e:
                # If the test prep work found a missing file, just record it in the
                # list of missing files. This allows all missing files for the
                # test suite to be reported at once
                setup.missing_files.append(str(e))
                continue

            # Check for any missing files
            missing_files = test.check_for_missing_files()
            if len(missing_files) > 0:
                setup.missing_files += missing_files
                continue

            # Skip the test if it wasn't selected by the command line
            if pargs.prep_only and test_descr['type'] != TestPhase.PREP:
                continue

            if not flg_reduce and test_descr['type'] == TestPhase.REDUCE:
                continue

            if not flg_after and test_descr['type'] == TestPhase.AFTERBURN:
                continue

            if not flg_ql and test_descr['type'] == TestPhase.QL:
                continue

            setup.tests.append(test)

    return setup



