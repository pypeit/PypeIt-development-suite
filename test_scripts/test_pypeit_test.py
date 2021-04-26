#
# See top-level LICENSE.rst file for Copyright information
#
# -*- coding: utf-8 -*-

"""
Pytest based test functions and Mock classes for testing the PypeIt dev suite scripts.

Rather than a series of bottom up style unit tests, these tests are intended to
cover as much code as possible to avoid finding out there's a syntax error in a script
at the end of a 10 hour test run. To accomplish this the tests call the pypeit_test.main() function
with mock objects and functions to replace calling the actual dev suite tests.

Although the tests attempt to verify the results of pypeit_test.main(), it is usually helpful to check the stdout
from the tests. This can be done with:

    python -m  pytest --capture=no -v test_scripts/test_pypeit_test.py

To get the coverage of the tests:

    pip install coverage
    coverage run -m pytest test_scripts/test_pypeit_test.py
    coverage report -m

"""

import pytest
import subprocess
import sys
import os
import random
from test_scripts import test_main
from test_scripts.pypeit_tests import PypeItReduceTest
import time


class MockPopen(object):
    """
    Mock of the Popen objects returned by subprocess.Popen.

    Attributes:

        pid (int):           Randomly generate integer to simulate a process id
        returncode (int):    Simulated return code from the child process. 0 unless failure_case is True
        failure_case (bool): If set to True by __init__, causes the simulated test to appear to fail. Defaults to False.
    """

    def __init__(self, failure_case=False):
        self.pid = random.randint(1,99999)
        self.returncode = 0
        self.failure_case = failure_case

    def wait(self):
        t = random.randint(1, 3)
        time.sleep(t)
        if self.failure_case:
            self.returncode = 1
        return 0

    def terminate(self):
        pass

class MockCompletedProcess(object):
    """
    Mock of the CompletedProcess objects returned by subprocess.run.  Used to simulate the results of tailing the logs
    from the pypeit test scripts
    """

    def __init__(self, returncode=0):
        self.stdout = b"This is \nSample log output\nFor unit testing\n"
        self.returncode = returncode

def mock_popen(*args, **kwargs):
    """
    Mock function for subprocess.Popen()
    """
    return MockPopen()

def mock_run(*args, **kwargs):
    """
    Mock function for subprocess.run()
    """
    return MockCompletedProcess()

def mock_failed_run(*args, **kwargs):
    """
    Mock function for subprocess.run() that has a returncode not equal to zero
    """
    return MockCompletedProcess(1)

def mock_raises_run(*args, **kwargs):
    """
    Mock function for subprocess.run() that raises an exception
    """
    raise RuntimeError("Unit testing Exception")


def create_dummy_files(base_path, files):
    """
    Utility function to create dummy files for testing. The files allow tests that check for the existence of files
    created by previous tests to pass.

    Args:
        base_path (:obj:`pathlib.Path`): The base path to create the files in
        files (sequence):                The relative paths of the files to create.

    """
    for file in files:
        full_path = base_path / file
        full_path.parent.mkdir(parents=True, exist_ok=True)
        with open(full_path, 'w') as f:
            print("dummy content", file=f)


def change_dir(new_directory):
    """
    Utility method to change the local directory with a context manager that will restore
       the original directory.  It is used to prevent a test from affecting other tests
       when changing the directory
    """
    class cd_context_manager:
        def __init__(self, dir):
            self.new_dir = dir
            self.old_dir = None

        def __enter__(self):
            self.old_dir = os.getcwd()
            os.chdir(self.new_dir)

        def __exit__(self, exc_type, exc_val, exc_tb):
            if self.old_dir is not None:
                os.chdir(self.old_dir)

            return False # Indicate an exception should propagate beyond the context

    return cd_context_manager(new_directory)

def test_main_with_missing_file_failures(monkeypatch, tmp_path):
    """
    Test test_main.main() when there are failures due to missing files
    """
    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        # Test failure to generate pypeit file
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '1',
                                          '-i', 'shane_kast_blue', '-s', '600_4310_d55', 'shane_kast_blue'])
        assert test_main.main() == 1

        # Test missing input for sensfunc
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '1',
                                          '-i', 'gemini_gnirs', '-s', '32_SB_SXD', 'gemini_gnirs'])
        assert test_main.main() == 1

        # Test not being able to find unique input for sensfunc
        missing_files = ['gemini_gnirs/32_SB_SXD/Science/spec1d_cN20170331S0206-HIP62745_GNIRS_2017Mar31T083351.681.fits',
                         'gemini_gnirs/32_SB_SXD/Science/spec1d_cN20170331S0206-HIP62745_GNIRS_2017Mar31T083351.682.fits']
        create_dummy_files(tmp_path, missing_files)

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '1',
                                          '-i', 'gemini_gnirs', '-s', '32_SB_SXD', 'gemini_gnirs'])
        os.chdir(tmp_path)
        assert test_main.main() == 1

        # Make sure the test_priority_list doesn't exist if there's a failure
        assert not os.path.exists(tmp_path / "test_priority_list")


def test_main_with_build_command_failure(monkeypatch, tmp_path):
    """
    Test test_main.main() when there are exceptions building commands
    """

    with monkeypatch.context() as m2:

        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        # Test error building command
        def mock_build_command_line(self):
            raise RuntimeError("Unit testing Exception")

        monkeypatch.setattr(PypeItReduceTest, "build_command_line", mock_build_command_line)
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', '-v',
                                          '-i', 'gemini_gnirs', 'develop'])

        assert test_main.main() == 2

def test_main_with_test_failure(monkeypatch, tmp_path):
    """
    Test test_main.main() when a pypeit test script returns a failure.
    """

    with monkeypatch.context() as m3:

        # Test error running command, with verbose and multiple threads
        def mock_failure_popen(*args, **kwargs):
            return MockPopen(True)

        monkeypatch.setattr(subprocess, "Popen", mock_failure_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', '-v', 'shane_kast_blue'])
        assert test_main.main() == 3


def test_main_develop_without_failures(monkeypatch, tmp_path):
    """
    Test test_main.main() on a simulated dev suite run with no errors.
    """

    test_order = ['p200_dbsp_blue/600_4000_d55',
                  'shane_kast_blue/600_4310_d55',
                  'keck_lris_blue/multi_300_5000_d680',
                  'keck_lris_blue_orig/long_600_4000_d500',
                  'keck_lris_blue/multi_600_4000_d560',
                  'shane_kast_blue/452_3306_d57',
                  'keck_lris_blue/long_600_4000_d560',
                  'keck_lris_blue/long_400_3400_d560'
                  ]

    priority_list = tmp_path / 'test_priority_list'

    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', 'develop'])

        missing_files = ['gemini_gmos/GS_HAM_R400_860/Science/spec1d_S20181219S0316-GD71_GMOS-S_1864May27T230832.356.fits',
                         'gemini_gnirs/32_SB_SXD/Science/spec1d_cN20170331S0206-HIP62745_GNIRS_2017Mar31T083351.681.fits',
                         'shane_kast_blue/600_4310_d55/shane_kast_blue_A/shane_kast_blue_A.pypeit',
                         'shane_kast_blue/600_4310_d55/shane_kast_blue_A/Science/spec1d_b24-Feige66_KASTb_2015May20T041246.960.fits',
                         'keck_deimos/900ZD_LVM_5500/Science/spec1d_DE.20110729.54545-Feige110_DEIMOS_2011Jul29T150856.803.fits',
                         'keck_mosfire/Y_long/Science/spec1d_m191118_0064-GD71_MOSFIRE_2019Nov18T104704.507.fits']

        create_dummy_files(tmp_path, missing_files)

        # Change to the temp path so that the test_priority_list is written there
        with change_dir(tmp_path):
            # Write out a short test priority list to exercise the code to detect changes from a full run
            with open(priority_list, "w") as f:
                for setup_key in test_order:
                    print(setup_key, file=f)

            # Get the size of it to make sure it changes
            first_stat_info = priority_list.stat()

            assert test_main.main() == 0

        assert priority_list.exists()
        second_stat_info = priority_list.stat()
        assert first_stat_info.st_size < second_stat_info.st_size

def test_main_debug_with_verbose_and_report(monkeypatch, tmp_path):
    """
    Test test_main.main() with the --debug option, verbose output, and an external report file
    """
    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        # This will include a failure and skipped tests as well as passed tests
        report_path = tmp_path / 'test_output.txt'
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', '--debug',
                                          '-v', '-r', str(report_path), 'develop'])
        assert test_main.main() == 1
        # Verify the report exists and contains data
        stat_result = report_path.stat()
        assert stat_result.st_size > 0


def test_main_debug_priority_list(monkeypatch, tmp_path, capsys):
    """
    Test test_main.main() with the --debug option, and make sure tests run in the order given by the
    test_priority list
    """

    test_order = ['p200_dbsp_blue/600_4000_d55',
                  'shane_kast_blue/600_4310_d55',
                  'keck_lris_blue/multi_300_5000_d680',
                  'keck_lris_blue_orig/long_600_4000_d500',
                  'keck_lris_blue/multi_600_4000_d560',
                  'shane_kast_blue/452_3306_d57',
                  'keck_lris_blue/long_600_4000_d560',
                  'keck_lris_blue/long_400_3400_d560'
                  ]

    # The last test won't get written to the test priority list file so the code to assign the default
    # priority executes
    last_test = 'shane_kast_blue/830_3460_d46'


    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        # Run from the tmp_path so the test priority list is separated from other tests
        with change_dir(tmp_path):

            # Write out a test priority list
            with open('test_priority_list', "w") as f:
                for setup_key in test_order:
                    print(setup_key, file=f)

            # Add the last test back to the test_order now that the priority list is written
            test_order.append(last_test)

            # Prevent debug suite from failing
            missing_files = ['shane_kast_blue/600_4310_d55/shane_kast_blue_A/shane_kast_blue_A.pypeit',
                             'shane_kast_blue/600_4310_d55/shane_kast_blue_A/Science/spec1d_b24-Feige66_KASTb_2015May20T041246.960.fits']
            create_dummy_files(tmp_path, missing_files)

            monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '--debug',
                                              '-v', 'develop'])
            assert test_main.main() == 0

            # Use captured stdout to make sure test setups were run in the order in test_order
            # We verify this by splitting the stdout into words and looking for "STARTED" and the test setup
            # listed after it. This setup is compared against the next expected setup in test_order, but care is taken
            # to deal with test setups that start multiple tests
            output = capsys.readouterr().out
            output_words = output.split()
            current_pos = 0
            found_started = False
            last_setup_found = ""
            for word in output_words:
                if word == 'STARTED':
                    found_started = True
                    continue

                if found_started == True:
                    found_started = False
                    if word == last_setup_found:
                        # Don't check against the next test setup if this is a test in the previous test setup
                        continue
                    # If more tests have been added to the debug tests, don't get an index out of
                    # bound exception
                    if current_pos >= len(test_order):
                        break
                        
                    assert word == test_order[current_pos]
                    last_setup_found = word
                    current_pos += 1

def test_main_with_quiet(monkeypatch, tmp_path, capsys):
    """
    Test test_main.main() with the --quiet option
    """

    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)


        # This will include a failure and skipped tests as well as passed tests
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', '--debug',
                                          '-q', 'develop'])
        assert test_main.main() == 1

        # Verify quiet mode rally is quiet
        captured = capsys.readouterr()
        assert len(captured.out) == 0
        assert len(captured.err) == 0

        # Verify quiet mode created a test report in default location
        default_report_path = tmp_path / 'pypeit_test_results.txt'
        stat_result = default_report_path.stat()
        assert stat_result.st_size > 0

def test_main_with_all(monkeypatch, tmp_path):
    """
    Test test_main.main() with the "all" argument.
    """
    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)


        # This will fail due to missing files
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', 'all'])

        with pytest.raises(ValueError):
            test_main.main()

def test_main_in_steps(monkeypatch, tmp_path):
    """
    Test test_main.main() prep-only, reduce, afterburn, and ql options.
    """

    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        missing_files = ['shane_kast_blue/600_4310_d55/shane_kast_blue_A/shane_kast_blue_A.pypeit',
                         'shane_kast_blue/600_4310_d55/shane_kast_blue_A/Science/spec1d_b24-Feige66_KASTb_2015May20T041246.960.fits']
        create_dummy_files(tmp_path, missing_files)

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4',
                                          '-i', 'shane_kast_blue', '--prep_only', 'develop'])
        assert test_main.main() == 0
        assert (tmp_path / "shane_kast_blue" / "452_3306_d57" / "shane_kast_blue_452_3306_d57.pypeit").exists()

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4',
                                          '-i', 'shane_kast_blue', '-s', '600_4310_d55', 'reduce'])
        assert test_main.main() == 0

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4',
                                          '-i', 'shane_kast_blue', '-s', '600_4310_d55', 'afterburn'])
        assert test_main.main() == 0

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4',
                                          '-i', 'shane_kast_blue', '-s', '600_4310_d55', 'ql'])
        assert test_main.main() == 0

def test_quick_look_build_masters_failure(monkeypatch, tmp_path):
    """
    Test test_main.main() ql tests that fail building masters.
    """

    # Failure from the build script returning non-zero
    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_failed_run)
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-i', 'keck_nires', 'ql'])

        assert test_main.main() == 1

    # Failure from an exception raised when running the build script
    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_raises_run)
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-i', 'keck_nires', 'ql'])

        assert test_main.main() == 1

def test_quick_look_masters_env(monkeypatch, tmp_path):
    """
    Test test_main.main() ql tests with and without a NIRES_MASTERS environment variable.
    """

    # Generate masters without the environment variable set
    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)
        monkeypatch.delenv('NIRES_MASTERS', raising=False)
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-i', 'keck_nires', 'ql'])

        assert test_main.main() == 0

    # Generate masters with the environment variable set
    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)
        monkeypatch.setenv('NIRES_MASTERS', str(tmp_path))
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-i', 'keck_nires', 'ql'])

        assert test_main.main() == 0
