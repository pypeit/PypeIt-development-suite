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
testing PypeIt (see below): `RAW_DATA`, `Cooked`, and `CALIBS`.

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
    ln -s /Volumes/GoogleDrive/Team\ Drives/PHYS-GP-Hennawi/PypeIt/PypeIt-development-suite/Cooked  Cooked
    ln -s /Volumes/GoogleDrive/Team\ Drives/PHYS-GP-Hennawi/PypeIt/PypeIt-development-suite/CALIBS  CALIBS
    ```

## Testing PypeIt

This suite includes several examples of the PypeIt data reduction process 
for some of the supported instruments.

If you have [PypeIt](https://github.com/pypeit/PypeIt) installed, you
should be able to run the tests in this directory using the
`pypeit_test` script:

```
$ $PYPEIT_DEV/pypeit_test -h
usage: pypeit_test [-h] [-o OUTPUTDIR] [-i INSTRUMENT] [-s SETUP] [--debug] [-p] [-m] [-t THREADS] [-q] [-v] [-r REPORT] tests

Run pypeit tests on a set of instruments. Typical call for testing pypeit when developing new code is `./pypeit_test develop`. Execution requires you to have a PYPEIT_DEV environmental variable, pointing to the top-level directory of the dev-suite repository
(typically the location of this script). Raw data for testing is expected to be at ${PYPEIT_DEV}/RAW_DATA. To run all tests for the supported instruments, use 'develop'. To only run the basic reductions, use 'reduce'. To only run the tests that use the results of
the reductions, use 'afterburn'. To run all possible tests (beware!), use 'all'.

positional arguments:
  tests                 Instrument or test to run. For instrument-specific tests, you can provide the telescope or the spectrograph, but beware of non-unique matches. E.g. 'mage' selects all the magellan instruments, not just 'magellan_mage'. Options include:
                        develop, reduce, afterburn, all, ql, alfosc, apf, binospec, dbsp, deimos, esi, fire, flamingos2, fors2, gemini, gmos, gnirs, hires, isis, kast, kcwi, keck, lbt, levy, lris, luci, mage, magellan, mdm, mmirs, mmt, mods, mosfire, nires,
                        nirspec, not, osmos, p200, shane, tspec, vlt, wht, xshooter

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Output folder. (default: REDUX_OUT)
  -i INSTRUMENT, --instrument INSTRUMENT
                        Restrict to input instrument (default: None)
  -s SETUP, --setup SETUP
                        Single out a setup to run (default: None)
  --debug               Debug using only blue setups (default: False)
  -p, --prep_only       Only prepare to execute run_pypeit, but do not actually run it. (default: False)
  -m, --masters         run pypeit using any existing masters (default: False)
  -t THREADS, --threads THREADS
                        Run THREADS number of parallel tests. (default: 1)
  -q, --quiet           Supress all output to stdout. If -r is not a given, a report file will be written to <outputdir>/pypeit_test_results.txt (default: False)
  -v, --verbose         Output additional detailed information while running the tests and output a detailed report at the end of testing. This has no effect if -q is given (default: False)
  -r REPORT, --report REPORT
                        Write a detailed test report to REPORT. (default: None)
```

For example, to run all the DEIMOS setups:

```
./pypeit_test deimos
```

or just the `1200G_M_7750` setup specifically:

```
./pypeit_test deimos -s 1200G_M_7750
```

or the development suite for full tests

```
./pypeit_test develop
```

## Parallel Testing
The development suite currently takes over 12 hours to run. This can be sped up by running parallel tests:
```
./pypeit_test -t 2 develop
```

The number of threads that can be run depends on the amount of memory available. As a rule of thumb 1 thread per 16G of
available memory should be safe.  For systems with more virtual CPUs than physical CPU cores (i.e. Hyperthreading) the 
number of threads should not exceed the number of physical cores, or else there could be a performance hit as threads 
compete for resources.  

