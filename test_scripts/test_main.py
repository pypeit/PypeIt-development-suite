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
from IPython import embed

import numpy as np

from .test_setups import TestPhase, all_tests, develop_setups, supported_instruments
from .test_setups import cooked_setups, ql_setups
from .pypeit_tests import get_unique_file

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
    def __init__(self, pargs, setups):
        self.pargs = pargs
        self.test_setups = setups
        self.start_time = None
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

        self.testing_started()


    def _get_test_counts(self):
        """Helper method to create a string with the current test counts"""
        verbose_info = f'{self.num_active:2} active/' if self.pargs.verbose else ''
        return f'{verbose_info}{self.num_passed:2} passed/{self.num_failed:2} failed/{self.num_skipped:2} skipped'

    def testing_started(self):
        """Called once testing has started"""
        with self.lock:
            self.start_time = datetime.datetime.now()

            # Create the report file and write the header to it
            if self.pargs.report:
                try:
                    with open(self.pargs.report, "w") as report_file:
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
        """Called once all testing is complete"""
        self.end_time = datetime.datetime.now()
        if self.pargs.report is not None:
            with open(self.pargs.report, "a") as report_file:
                print ("-------------------------", file=report_file)
                self.summary_report(report_file)

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

    def summary_report(self, output=sys.stdout):
        """Display a summary report on the results of testing to the given output stream"""

        masters_text = 'reused' if self.pargs.masters else 'ignored'
        if self.num_tests == self.num_passed:
            print("\n" + "\x1B[" + "1;32m" +
                  "--- PYPEIT DEVELOPMENT SUITE PASSED {0}/{1} TESTS (Masters {2}) ---".format(
                      self.num_passed, self.num_tests, masters_text)
                  + "\x1B[" + "0m" + "\r", file=output)
        else:
            print("\n" + "\x1B[" + "1;31m" +
                  "--- PYPEIT DEVELOPMENT SUITE FAILED {0}/{1} TESTS (Masters {2}) ---".format(
                      self.num_failed, self.num_tests, masters_text)
                  + "\x1B[" + "0m" + "\r", file=output)
            print('Failed tests:', file=output)
            for t in self.failed_tests:
                print('    {0}'.format(t), file=output)
            print('Skipped tests:', file=output)
            for t in self.skipped_tests:
                print('    {0}'.format(t), file=output)

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
        print(f"Command:    {' '.join(test.command_line) if test.command_line is not None else ''}", file=output, flush=flush)
        print('', file=output, flush=flush)
        print('Error Messages:', file=output, flush=flush)

        for msg in test.error_msgs:
            print(msg, file=output, flush=flush)

        print('', file=output, flush=flush)
        print("End of Log:", file=output, flush=flush)
        if test.logfile is not None and os.path.exists(test.logfile):
            result = subprocess.run(['tail', '-3', test.logfile], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            print(result.stdout.decode(), file=output, flush=flush)

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


def raw_data_dir():
    return os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA')


def available_data():
    walk = os.walk(raw_data_dir())
    return next(walk)[1]


def parser(options=None):
    import argparse

    dirs = available_data()
    all_tests = np.unique([ [d.split('_')[0], d.split('_')[1]] for d in dirs ])

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Run pypeit tests on a set of instruments.  '
                                                 'Typical call for testing pypeit when developing '
                                                 'new code is `./pypeit_test develop`.  Execution '
                                                 'requires you to have a PYPEIT_DEV environmental '
                                                 'variable, pointing to the top-level directory '
                                                 'of the dev-suite repository (typically the '
                                                 'location of this script).  Raw data for testing '
                                                 'is expected to be at ${PYPEIT_DEV}/RAW_DATA.  '
                                                 'To run all tests for the supported instruments, '
                                                 'use \'develop\'.  To only run the basic '
                                                 'reductions, use \'reduce\'.  To only run the '
                                                 'tests that use the results of the reductions, '
                                                 'use \'afterburn\'.  To run all possible tests '
                                                 '(beware!), use \'all\'.')

    parser.add_argument('tests', type=str, default=None,
                        help='Instrument or test to run.  For instrument-specific tests, you '
                             'can provide the telescope or the spectrograph, but beware of '
                             'non-unique matches.  E.g. \'mage\' selects all the magellan '
                             'instruments, not just \'magellan_mage\'.  Options include: '
                             'develop, reduce, afterburn, all, ql, {0}'.format(', '.join(all_tests)))
    parser.add_argument('-o', '--outputdir', type=str, default='REDUX_OUT',
                        help='Output folder.')
    # TODO: Why is this an option?
    parser.add_argument('-i', '--instrument', type=str, help="Restrict to input instrument")
    parser.add_argument('-s', '--setup', type=str, help="Single out a setup to run")
    parser.add_argument('--debug', default=False, action='store_true',
                        help='Debug using only blue setups')
    parser.add_argument('-p', '--prep_only', default=False, action='store_true',
                        help='Only prepare to execute run_pypeit, but do not actually run it.')
    parser.add_argument('-m', '--masters', default=False, action='store_true',
                        help='run pypeit using any existing masters')
    parser.add_argument('-t', '--threads', default=1, type=int,
                        help='Run THREADS number of parallel tests.')
    parser.add_argument('-q', '--quiet', default=False, action='store_true',
                        help='Supress all output to stdout. If -r is not a given, a report file will be '
                             'written to <outputdir>/pypeit_test_results.txt')
    parser.add_argument('--no_gui', default=False, action='store_true',
                        help='Supress any GUIs displayed by any tests.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Output additional detailed information while running the tests and output a '
                             'detailed report at the end of testing. This has no effect if -q is given')
    parser.add_argument('-r', '--report', default=None, type=str,
                        help='Write a detailed test report to REPORT.')
    return parser.parse_args() if options is None else parser.parse_args(options)

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
    outputdir = os.path.abspath(pargs.outputdir)

    # Make sure we can create a writable report file (if needed) before starting the tests
    if pargs.report is None and pargs.quiet:
        # If there's no report file specified in the command line, but we're in quiet mode,
        # make up a report file name
        pargs.report = get_unique_file(os.path.join(outputdir, "pypeit_test_results.txt"))

    if pargs.report:
        try:
            report_file = open(pargs.report, "w")
        except Exception as e:
            print(f"Could not open report file {pargs.report}", file=sys.stderr)
            traceback.print_exc()
            sys.exit(1)

    # Only write the test priority file if all develop tests are being run
    if pargs.tests == "develop" and pargs.instrument is None and pargs.setup is None and pargs.debug is False:
        write_priorities = True
    else:
        write_priorities = False

    # ---------------------------------------------------------------------------
    # Determine which instruments and setups will be tested

    all_instruments = available_data()
    flg_after = False
    flg_ql = False

    # Development instruments (in returned dictonary keys) and setups
    devsetups = develop_setups
    tests_that_only_use_dev_setups = ['develop', 'reduce', 'afterburn']

    # Cooked instruments (in returned dictonary keys) and setups
    cooksetups = cooked_setups

    # Setup
    unsupported = []
    if pargs.tests == 'all':
        instruments = np.array([item for item in all_instruments
                                        for inst in supported_instruments
                                            if inst.lower() in item.lower()])
    elif pargs.tests in tests_that_only_use_dev_setups:
        instruments = np.array(list(devsetups.keys())) if pargs.instrument is None \
                        else np.array([pargs.instrument])
        if pargs.tests == 'afterburn':
            # Only do the flux-calibration and coadding tests
            flg_after = True
    elif pargs.tests.lower() == 'cooked':
        instruments = np.array(list(cooksetups.keys())) if pargs.instrument is None \
                        else np.array([pargs.instrument])
    elif pargs.tests.lower() == 'ql':      
        flg_ql = True    
        instruments = np.array(list(ql_setups.keys())) if pargs.instrument is None \
                        else np.array([pargs.instrument])
    else:
        instruments = np.array([item for item in all_instruments 
                                    if pargs.tests.lower() in item.lower()])
        unsupported = [item for item in instruments 
                            if not np.any([inst.lower() in item.lower()
                                for inst in supported_instruments]) ]

    # Check that instruments is not blank
    if len(instruments) == 0:
        print("\x1B[" + "1;31m" + "\nERROR - " + "\x1B[" + "0m" +
              "Invalid test selected: {0:s}\n\n".format(pargs.tests) +
              "Consult the help (pypeit_test -h) or select one of the " +
              "available RAW_DATA directories: {0}".format(', '.join(all_instruments)))
        return 1

    if len(unsupported) > 0:
        if not pargs.quiet:
            print("\x1B[" + "1;33m" + "\nWARNING - " + "\x1B[" + "0m" +
                  "The following tests have not been validated and may not pass: {0}\n\n".format(
                  unsupported))

    # Report
    if not pargs.quiet:
        print('Running tests on the following instruments:')
        for instr in instruments:
            print('    {0}'.format(instr))
        print('')

    # ---------------------------------------------------------------------------
    # Build the TestSetup and PypeItTest objects for testing

    # Load test setup priority from file
    priority_list = TestPriorityList('test_priority_list')
    if not pargs.quiet and pargs.verbose:
        print(f'Loaded {len(priority_list)} setup priorities')

    setups = []
    missing_files = []
    for instr in instruments:
        # Only do blue instruments
        if pargs.debug and 'blue' not in instr:
            continue

        # Setups
        setup_names = next(os.walk(os.path.join(raw_data, instr)))[1]
        if pargs.setup is not None and pargs.setup not in setup_names:
            # No setups selected
            continue
        # Limit to a single setup
        if pargs.setup is not None:
            setup_names = [ pargs.setup ]
        # Limit to development setups
        elif pargs.tests in tests_that_only_use_dev_setups:
            setup_names = devsetups[instr]
        elif pargs.tests.lower() == 'cooked':
            setup_names = cooksetups[instr]
        elif pargs.tests.lower() == 'ql':
            setup_names = ql_setups[instr]

        # Build test setups, check for missing files, and run any prep work
        for setup_name in setup_names:

            setup = build_test_setup(pargs, instr, setup_name, flg_after, flg_ql)
            missing_files += setup.missing_files

            # set setup priority from file
            priority_list.set_test_setup_priority(setup)

            setups.append(setup)

        # Report
        if not pargs.quiet:
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
    test_report = TestReport(pargs, setups)


    # Add tests to the test_run_queue
    for setup in setups:
        test_run_queue.put(setup)

    # Start threads to run the tests
    if not pargs.quiet and pargs.threads > 1:
        print(f'Running tests in {pargs.threads} parallel processes')

    test_report = TestReport(pargs, setups)
    test_report.start_time = datetime.datetime.now()

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

    # ---------------------------------------------------------------------------
    # Build test priority list for next time, but only if all tests succeeded
    # and all develop tests were being run
    if write_priorities and test_report.num_passed == test_report.num_tests:
        priority_list.update_priorities(setups)
        priority_list.write()
        if not pargs.quiet and pargs.verbose:
            print(f'Wrote {len(priority_list)} setup priorities')


    # ---------------------------------------------------------------------------
    # Finish up the report on the test results
    test_report.testing_completed()

    if not pargs.quiet:
        if pargs.verbose:
            test_report.detailed_report()
        else:
            test_report.summary_report()

    return test_report.num_failed


def build_test_setup(pargs, instr, setup_name, flg_after, flg_ql):
    """Builds a TestSetup object including the tests that it will run"""

    dev_path = os.getenv('PYPEIT_DEV')
    raw_data = raw_data_dir()
    outputdir = os.path.abspath(pargs.outputdir)

    # Directory with raw data
    rawdir = os.path.join(raw_data, instr, setup_name)

    # Directory for reduced data
    rdxdir = os.path.join(outputdir, instr, setup_name)
    if not os.path.exists(rdxdir):
        # Make the directory
        os.makedirs(rdxdir)

    # Create the test setup and set it's priority
    setup = TestSetup(instr, setup_name, rawdir, rdxdir, dev_path)

    # Go through each test for this setup and add it to the setup if it's
    # selected by the command line arguments
    for test_descr in all_tests:
        if '*' in test_descr['setups']:
            key = '*'
        elif setup.instr in test_descr['setups']:
            key = setup.instr
        elif setup.key in test_descr['setups']:
            key = setup.key
        else:
            continue

        if isinstance(test_descr['setups'], dict):
            kwargs = test_descr['setups'][key]
        else:
            kwargs = dict()

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

        if pargs.tests == "reduce" and test_descr['type'] not in (TestPhase.PREP, TestPhase.REDUCE):
            continue

        if flg_after and test_descr['type'] not in (TestPhase.PREP, TestPhase.AFTERBURN):
            continue

        if flg_ql and test_descr['type'] not in (TestPhase.PREP, TestPhase.QL):
            continue

        setup.tests.append(test)

    return setup



