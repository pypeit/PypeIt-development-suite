# The Python Spectroscopic Data Reduction Pipeline (PypeIt) development suite

This repository provides data and input files used for testing
[PypeIt](https://github.com/pypeit/PypeIt).

## Installation

To install, run:

`git clone https://github.com/pypeit/PypeIt-development-suite.git`

This pulls the versioned files into a directory called
`PypeIt-development-suite`.  Both the testing script in this repo and
the tests that are included with the PypeIt code expect an environmental
variable that points to this directory.  For example,

 - If you use c-shell, add the following to your `~/.cshrc`:

   ```
   setenv PYPEIT_DEV ${HOME}/PypeIt-development-suite
   ```

 - If you use bash, add the following to your `~/.bashrc` or
   `~/.bash_profile`:

   ```
   export PYPEIT_DEV=${HOME}/PypeIt-development-suite
   ```

## Data Access

Given its volume, this repo does not contain the raw data.  Instead the
raw test data are hosted by in a 
[Google Drive](https://drive.google.com/drive/folders/1oh19siB1-F0jjmY-F_jr73eA-TQYEiFW?usp=sharing).  

The easiest way to access the
development suite data is download and use [Google File
Stream](https://support.google.com/drive/answer/7329379?hl=en).  Google
File Stream is a Dropbox-like application that syncs your Google Drive
with a local directory on your machine.  

If you cannot or would rather
not use Google File Stream, you can simply download the appropriate
files directly using the Google Drive web interface (although this can
be a bit onerous and does not keep in sync with the remote Team Drive).

Using Google File Stream, the PypeIt team drive you will be able to
access the development suite data at the path: 

```
/Volumes/GoogleDrive/Team\ Drives/PHYS-GP-Hennawi/PypeIt/PypeIt-development-suite/
```

The Team Drive contains three directories that should be accessible for
testing PypeIt (see below): `RAW_DATA`, and `CALIBS`.

  - If you download the data directly from the Team Drive, place them in
    your `PypeIt-development-suite` directory.  **Make sure that you *do
    not* add these directories to the repo!**

  - If you're using Google File Stream, add symlinks to your
    `PypeIt-development-suite` directory as follows (be sure to include
    the \ in the Team\ Drives otherwise the space in "Team Drives" will
    cause problems):

    ```
    cd $PYPEIT_DEV
    ln -s /Volumes/GoogleDrive/Team\ Drives/PHYS-GP-Hennawi/PypeIt/PypeIt-development-suite/RAW_DATA  RAW_DATA
    ln -s /Volumes/GoogleDrive/Team\ Drives/PHYS-GP-Hennawi/PypeIt/PypeIt-development-suite/CALIBS  CALIBS
    ```

## Testing PypeIt

This suite includes several examples of the PypeIt data reduction process 
for some of the supported instruments.

If you have [PypeIt](https://github.com/pypeit/PypeIt) installed, you
should be able to run the tests in this directory using the
`pypeit_test` script:

```
$ $PYPEIT_DEV/pypeit_test all
```
The above runs all pypeit tests, including the unit tests in the PypeIt 
installation. To run a subset of tests you can specify one or more test type. For example:

```
# Run the non-pytest tests, similar to develop in prior versions
$ $PYPEIT_DEV/pypeit_test reduce afterburn ql

# Run reduce tests followed by the vet pytest tests
$ $PYPEIT_DEV/pypeit_test reduce vet
```

All of the supported test test types are shown in the table below.

| Test Name   | Description                                                    |
|-------------|----------------------------------------------------------------|
|pypeit_tests | Runs the pytest tests installed with PypeIt in pypeit/tests. These tests are self contained and can run in CI.                              |
|unit         | Runs the pytest tests installed in $PYPEIT_DEV/unit. These tests require the RAW_DATA directory.                                          |
|reduce       | Runs the reduction tests that call run_pypeit directly.        |
|afterburn    | Runs the PypeIt tests that directly call PypeIt post-reduction  tools (e.g. flux calibration, coadding, etc).                                  |
|ql           | Runs the Quick Look tests.                                     |
|vet          | Runs the pytest tests that verify the results from earlier PypeIt tests.                                                                  |
|all          | Runs all of the above, in the order listed above.              |
|list         | This does not run any tests, instead it lists all of the supported instruments and setups. (See below).                                 |


## Running Unit and Vet tests separately

The unit and vet tests can also be run directly using pytest. For example:

```
$ cd $PYPEIT_DEV

# Run all dev-suite unit tests
$ pytest unit_tests  

# Run all dev-suite vet tests
$ pytest vet_tests   

# Run the script tests in both unit_tests and vet_tests
$ pytest unit_tests/test_scripts.py vet_tests/test_scripts.py
```

See the [pytest docs](https://docs.pytest.org/) for more information on running pytest.

## Selecting test setups and instruments to test

The ``-i`` and ``-s`` options to ``pypeit_test`` can be used to select multiple
instruments and setups to test. Setups can be specified in conjunction with
an instrument or by using a ``/`` between instrument and setup name. For example:

```
# Run all of pytest tests and all of the shane_kast_blue and shane_kast_red tests
cd $PYPEIT_DEV
$ ./pypeit_test all -i shane_kast_blue shane_kast_red

# Run one reduce test from shane_kast_blue and shane_kast_red respectively
$ ./pypeit_test reduce -i shane_kast_blue shane_kast_red -s 600_4310_d55 600_7500_d57

# Run the same tests as above, using the / syntax
$ ./pypeit_test reduce -s shane_kast_blue/600_4310_d55 shane_kast_red/600_7500_d57

```

Run ``pypeit_test list`` to see a list of all supported instruments and setups.

## Test Reports

A test report can be generated using the ``-r`` option, for example:

```
$PYPEIT_DEV/pypeit_test all -r test_report.txt
```

The contents of the report contain additional information about the test
setups run. Below is an example of for one test setup:

```
-------------------------
Test Setup: keck_mosfire/Y_multi

-------------------------
Directories:
         Raw data: /PypeIt-development-suite/RAW_DATA/keck_mosfire/Y_multi
    PypeIt output: /tmp/REDUX_OUT/keck_mosfire/Y_multi
Files:
     .pypeit file: None
 Std .pypeit file: None
Tests:
----
keck_mosfire/Y_multi pypeit (without masters) Result: --- PASSED

Logfile:    /tmp/REDUX_OUT/keck_mosfire/Y_multi/keck_mosfire_y_multi.test.log
Process Id: 300
Start time: Tue May 17 17:12:43 2022
End time:   Tue May 17 17:47:25 2022
Duration:   0:34:41.349185
Command:    run_pypeit /tmp/REDUX_OUT/keck_mosfire/Y_multi/keck_mosfire_y_multi.pypeit -o

Error Messages:

End of Log:
[INFO]    :: run_pypeit.py 146 main() - Generating QA HTML
Wrote: /tmp/REDUX_OUT/keck_mosfire/Y_multi/QA/MF_A.html
Wrote: /tmp/REDUX_OUT/keck_mosfire/Y_multi/QA/MF_A.html
```

## Code coverage

The dev suite can also collect coverage data for ``PypeIt`` using [Coverage](https://coverage.readthedocs.io/) . To do this add ``--coverage <coverage report file>`` 
to the ``pypeit_test`` command:

```
$ cd $PYPEIT_DEV
$ ./pypeit_test all --coverage coverage_report_file.txt
```
The coverage report contains a file by file list of the coverage information, including missed lines. It ends with a summary of the total code coverage.
The unit tests, and deprecated sections of the ``PypeIt`` code base are omitted.
For example:

```
Name                                                               Stmts   Miss  Cover   Missing
------------------------------------------------------------------------------------------------
/home/dusty/work/PypeIt/pypeit/__init__.py                            20      5    75%   45-49
...
/home/dusty/work/PypeIt/pypeit/wavemodel.py                          330    286    13%   50-79, 108-137, 149-152, 195-211, 227-233, 247-257, 315-418, 459-535, 568-602, 627-648, 686-715, 765-788, 841-863
/home/dusty/work/PypeIt/pypeit/wavetilts.py                          240    104    57%   90, 204-205, 297, 311-319, 435-525, 562-581, 593, 598-600, 611-614, 627-630, 635, 661-681, 688, 714-717, 722-729
------------------------------------------------------------------------------------------------
TOTAL                                                              41785  22139    47%
```

## Parallel Testing
The development suite currently takes over 12 hours to run. This can be sped up by running parallel tests:
```
./pypeit_test -t 2 all
```

The number of threads that can be run depends on the amount of memory available. As a rule of thumb 1 thread per 16G of
available memory should be safe.  For systems with more virtual CPUs than physical CPU cores (i.e. Hyperthreading) the 
number of threads should not exceed the number of physical cores, or else there could be a performance hit as threads 
compete for resources.  

To keep all cpus active as long as possible ''pypeit_test'' runs the slowest tests first. To do this it needs the ``test_priority_list`` file which contains a list of all the test setups ordered from slowest to fastest. This file is re-written everytime a run of the full test suite passes, and should be 
kept up to date by periodically pushing it to git.

The pytest portion of the dev-suite currently cannot be run in parallel.

## Headless Testing
Some of the tests in the dev-suite will start GUI applications. To run in a
headless environment where this isn't possible, QT must still be installed.
To do so, first install ``PypeIt`` with QT5 support, as documented at https://pypeit.readthedocs.io/en/latest/installing.html. Next install the correct
QT package for the OS.


| OS Distribution   | Package name     |
|-------------------|------------------|
| Ubuntu 21.04      | qt5-default      |
| Ubuntu 22.04      | qtbase5-dev      |
| Centos 7          | qt5-qtbase-devel |

Finally, set the ``QT_QPA_PLATFORM`` environment variable to ``offscreen``. There is a convenience file
named ``source_headless_test.sh`` that will do this. For example:

```
source $PYPEIT_DEV/source_headless_test.sh
```

## Running in Nautilus
The dev-suite can be run in the [Nautilus cluster](https://ucsd-prp.gitlab.io/).
To generate the YAML for a dev-suite job, use ``gen_kube_devsuite``.  If needed,
a specific branch of both the PypeIt repository and the Pypeit-development-suite
repository can be chosen using ``-p`` and ``-d`` respectively.  These default to
``develop`` if not specified. The YAML file can then be ran using ``kubectl``.

For example to run using the pypeit_branch on the PypeIt repo and the devsuite_branch on the PypeIt-development-suite repo:

```
$ cd $PYPEIT_DEV/nauitilus
$ ./gen_kube_devsuite devsuite-job-name devsuite_job_file.yml -p pypeit_branch -d devsuite_branch
$ kubectl create -f devsuite_job_file.yml
```

The results of the dev-suite are copied to Nautilus S3 under ``s3://pypeit/Reports/``. And can be retrieved using the AWS CLI as follows:

```
aws --endpoint  https://s3-west.nrp-nautilus.io s3 cp s3://pypeit/Reports/devsuite-job-name.report .
```

``gen_kube_devsuite`` has additional code for generating coverage information and the test priority list. If ``--coverage`` and ``--priority_list`` are used, these files are also copied to S3:

```
$ ./gen_kube_devsuite coverage-job-name coverage_job_file.yml --coverage 
$ ./gen_kube_devsuite priority-job-name priority_job_file.yml --priority_list
$ kubectl create -f coverage_job_file.yml
$ kubectl create -f priority_job_file.yml
$ # Wait several hours 
$ export ENDPOINT=https://s3-west.nrp-nautilus.io 
$ aws --endpoint $ENDPOINT s3 cp s3://pypeit/Reports coverage-job-name.report .
$ aws --endpoint $ENDPOINT s3 cp s3://pypeit/Reports/priority-job-name.report .
$ aws --endpoint $ENDPOINT s3 cp s3://pypeit/Reports/coverage-job-name.coverage.report .
$ aws --endpoint $ENDPOINT s3 cp s3://pypeit/Reports/priority-job-name.test_priority_list .
```

Notice that ``--coverage`` can affect the performance of tests, so it's best
not to run it and ``--priority_list`` together.

To monitor a test in Nautilus as it is running, the logs can be tailed:

```
$ kubectl get pods

NAME                         READY    STATUS    RESTARTS
devsuite-job-name--1-fpjxw   1/1      RUNNING   0

$ kubectl logs -f devsuite-job-name--1-fpjxw
```

Additionally they can be monitored with the [Nautilus Grafana page](https://grafana.nrp-nautilus.io/?orgId=1).


## Additional Options

```
$ $PYPEIT_DEV/pypeit_test -h
usage: pypeit_test [-h] [-o OUTPUTDIR] [-i INSTRUMENTS [INSTRUMENTS ...]]
                   [-s SETUPS [SETUPS ...]] [--debug] [-p] [-m] [-t THREADS]
                   [-q] [-v] [--coverage COVERAGE] [-r REPORT] [-w]
                   tests [tests ...]

Run pypeit tests on a set of instruments. Typical call for testing pypeit when
developing new code is `./pypeit_test all`. Execution requires you to have a
PYPEIT_DEV environmental variable, pointing to the top-level directory of the
dev-suite repository (typically the location of this script). Raw data for
testing is expected to be at ${PYPEIT_DEV}/RAW_DATA. To run all tests for the
supported instruments, use 'all'. To only run the basic reductions, use
'reduce'. To only run the tests that use the results of the reductions, use
'afterburn''. Use 'list' to view all supported setups.

positional arguments:
  tests                 Which test types to run. Options are: pypeit_tests,
                        unit, reduce, afterburn, ql, vet, or all. Use list to
                        show all supported instruments and setups.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Output folder. (default: REDUX_OUT)
  -i INSTRUMENTS [INSTRUMENTS ...], --instruments INSTRUMENTS [INSTRUMENTS ...]
                        One or more instruments to run tests for. Use
                        "pypeit_test list" to see all supported instruments.
                        (default: None)
  -s SETUPS [SETUPS ...], --setups SETUPS [SETUPS ...]
                        One or more setups to run tests for. Use "pypeit_test
                        list" to see all supported setups. (default: None)
  --debug               Debug using only blue setups (default: False)
  -p, --prep_only       Only prepare to execute run_pypeit, but do not
                        actually run it. (default: False)
  -m, --do_not_reuse_masters
                        run pypeit without using any existing masters
                        (default: False)
  -t THREADS, --threads THREADS
                        Run THREADS number of parallel tests. (default: 1)
  -q, --quiet           Supress all output to stdout. If -r is not a given, a
                        report file will be written to
                        <outputdir>/pypeit_test_results.txt (default: False)
  -v, --verbose         Output additional detailed information while running
                        the tests and output a detailed report at the end of
                        testing. This has no effect if -q is given (default:
                        False)
  --coverage COVERAGE   Collect code coverage information. and write it to the
                        given file. (default: None)
  -r REPORT, --report REPORT
                        Write a detailed test report to REPORT. (default:
                        None)
  -w, --show_warnings   Show warnings when running unit tests and vet tests.
                        (default: False)
```

