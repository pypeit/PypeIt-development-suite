#
# See top-level LICENSE.rst file for Copyright information
#
# -*- coding: utf-8 -*-

"""
Test driver for the pypeit_test script.

Rather than a series of bottom up style unit tests, these tests are intended to
cover as much code as possible to avoid finding out there's a syntax error in pypeit_test
at the end of a 10 hour test run. To accomplish this the tests call the main() function in pypeit_test
with the subprocess calls to run tests replaced by mock objects.

Although the tests attempt to verify the results of pypeit_test.main(), it is usually helpful to check the stdout
from the tests. This can be done with:

    pytest --capture=no test_pypeit_test.py

To get the coverage of the tests:

    pip install coverage
    coverage run -m pytest test_pypeit_test.py
    coverage report -m

"""

import pytest

import subprocess
import sys
import random
import pypeit_test

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

    def __init__(self):
        self.stdout = b"This is \nSample log output\nFor unit testing\n"

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
        full_path.parent.mkdir(parents=True)
        with open(full_path, 'w') as f:
            print("dummy content", file=f)

def test_main_with_missing_file_failures(monkeypatch, tmp_path):
    """
    Test pypeit_test.main() when there are failures due to missing files
    """
    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        # Test failure to generate pypeit file
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '1',
                                          '-i', 'shane_kast_blue', '-s', '600_4310_d55', 'shane_kast_blue'])
        assert pypeit_test.main() == 1

        # Test missing input for sensfunc
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '1',
                                          '-i', 'gemini_gnirs', '-s', '32_SB_SXD', 'gemini_gnirs'])
        assert pypeit_test.main() == 1

def test_main_with_build_command_failure(monkeypatch, tmp_path):
    """
    Test pypeit_test.main() when there are exceptions building commands
    """

    with monkeypatch.context() as m2:

        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        # Test error building command
        def mock_build_command_line(self):
            raise RuntimeError("Unit testing Exception")

        monkeypatch.setattr(pypeit_test.PypeItReduceTest, "build_command_line", mock_build_command_line)
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', '-v',
                                          '-i', 'gemini_gnirs', 'develop'])

        assert pypeit_test.main() == 2

def test_main_with_test_failure(monkeypatch, tmp_path):
    """
    Test pypeit_test.main() when a pypeit test script returns a failure.
    """

    with monkeypatch.context() as m3:

        # Test error running command, with verbose and multiple threads
        def mock_failure_popen(*args, **kwargs):
            return MockPopen(True)

        monkeypatch.setattr(subprocess, "Popen", mock_failure_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', '-v', 'shane_kast_blue'])
        assert pypeit_test.main() == 3


def test_main_develop_without_failures(monkeypatch, tmp_path):
    """
    Test pypeit_test.main() on a simulated dev suite run with no errors.
    """

    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', 'develop'])

        missing_files = ['gemini_gmos/GS_HAM_R400_860/Science/spec1d_S20181219S0316-GD71_GMOS-S_1864May27T230832.356.fits',
                         'gemini_gnirs/32_SB_SXD/Science/spec1d_cN20170331S0206-HIP62745_GNIRS_2017Mar31T083351.681.fits',
                         'shane_kast_blue/600_4310_d55/shane_kast_blue_A/shane_kast_blue_A.pypeit',
                         'shane_kast_blue/600_4310_d55/shane_kast_blue_A/Science/spec1d_b24-Feige66_KASTb_2015May20T041246.960.fits']
        create_dummy_files(tmp_path, missing_files)
        assert pypeit_test.main() == 0

def test_main_debug_with_verbose_and_report(monkeypatch, tmp_path):
    """
    Test pypeit_test.main() with the --debug option, verbose output, and an external report file
    """
    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        # This will include a failure and skipped tests as well as passed tests
        report_path = tmp_path / 'test_output.txt'
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', '--debug',
                                          '-v', '-r', str(report_path), 'develop'])
        assert pypeit_test.main() == 1
        # Verify the report exists and contains data
        stat_result = report_path.stat()
        assert stat_result.st_size > 0

def test_main_with_quiet(monkeypatch, tmp_path, capsys):
    """
    Test pypeit_test.main() with the --quiet option
    """

    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)


        # This will include a failure and skipped tests as well as passed tests
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', '--debug',
                                          '-q', 'develop'])
        assert pypeit_test.main() == 1

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
    Test pypeit_test.main() with the "all" argument.
    """
    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)


        # This will fail due to missing files
        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4', 'all'])

        with pytest.raises(FileNotFoundError):
            pypeit_test.main()

def test_main_in_steps(monkeypatch, tmp_path):
    """
    Test pypeit_test.main() prep-only, reduce, afterburn, and ql options.
    """

    with monkeypatch.context() as m:
        monkeypatch.setattr(subprocess, "Popen", mock_popen)
        monkeypatch.setattr(subprocess, "run", mock_run)

        missing_files = ['shane_kast_blue/600_4310_d55/shane_kast_blue_A/shane_kast_blue_A.pypeit',
                         'shane_kast_blue/600_4310_d55/shane_kast_blue_A/Science/spec1d_b24-Feige66_KASTb_2015May20T041246.960.fits']
        create_dummy_files(tmp_path, missing_files)

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4',
                                          '-i', 'shane_kast_blue', '--prep_only', '-s', '600_4310_d55', 'develop'])
        assert pypeit_test.main() == 0

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4',
                                          '-i', 'shane_kast_blue', '-s', '600_4310_d55', 'reduce'])
        assert pypeit_test.main() == 0

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4',
                                          '-i', 'shane_kast_blue', '-s', '600_4310_d55', 'afterburn'])
        assert pypeit_test.main() == 0

        monkeypatch.setattr(sys, "argv", ['pypeit_test', '-o', str(tmp_path), '-t', '4',
                                          '-i', 'shane_kast_blue', '-s', '600_4310_d55', 'ql'])
        assert pypeit_test.main() == 0
