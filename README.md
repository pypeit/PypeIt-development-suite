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

    * If you use c-shell, add the following to your `~/.cshrc`:

```
setenv PYPEIT_DEV ${HOME}/PypeIt-development-suite
```

    * If you use bash, add the following to your `~/.bashrc` or
      `~/.bash_profile`:

```
export PYPEIT_DEV=${HOME}/PypeIt-development-suite
```

## Data Access

Given its volume, this repo does not contain the raw data.  Instead the
raw test data are hosted by in a Google Team Drive.  Please contact Joe
Hennawi <joe@physics.ucsb.edu> for access.

Once you have the correct permissions, the easiest way to access the
development suite data is download and use [Google File
Stream](https://support.google.com/drive/answer/7329379?hl=en).  Google
File Stream is a Dropbox-like application that syncs your Google Drive
with a local directory on your machine.  If you cannot or would rather
not use Google File Stream, you can simply download the appropriate
files directly using the Google Drive web interface (although this can
be a bit onerous and does not keep in sync with the remote Team Drive).

Using Google File Stream, the PypeIt team drive you will be able to
access the development suite data at the path: 

/Volumes/GoogleDrive/Team\ Drives/PHYS-GP-Hennawi/PypeIt/PypeIt-development-suite/

The Team Drive contains three directories that should be accessible for
testing PypeIt (see below): `RAW_DATA`, `Cooked`, and `CALIBS`.

    - If you download the data directly from the Team Drive, place them
      in your `PypeIt-development-suite` directory.  If you ever submit
      a PR to this dev suite, make sure that you **do not** add these
      directories to the repo!

    - If you're using Google File Stream, add symlinks to you
      `PypeIt-development-suite` directory as follows (be sure to
      include the \ in the Team\ Drives otherwise the space in "Team
      Drives" will cause problems):

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
usage: pypeit_test [-h] [-outputdir OUTPUTDIR] [-setup SETUP] tests

positional arguments:
  tests                 Tests to run. Options include: all, kast, lris,
                        deimos, nires, keck

optional arguments:
  -h, --help            show this help message and exit
  -outputdir OUTPUTDIR  Output folder. Default is ./REDUX_OUT (default: None)
  -setup SETUP          Single out a setup to run (default: None)
```

For example, to run all the DEIMOS setups:

```
./pypeit_test deimos
```

or just the `1200G_M_7750` setup specifically:

```
./pypeit_test deimos -s 1200G_M_7750
```


