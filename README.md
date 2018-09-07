# PYPIT-development-suite v0.4.1

The Python Spectroscopic Data Reduction Pipeline development suite

This repository contains a suite of test data and PYPIT input files. 
This suite includes several examples of the PYPIT data reduction process 
for some of the supported instruments.

If you have PYPIT installed, you should be able to run
the tests from this directory with:

./pypit_test all

Current allowed arguments include all, kast, isis, deimos, wht, and lris.

Use -outputdir to choose a different output directory.

The datasets are stored on a Google Team Drive. Contact Joe Hennawi at joe@physics.ucsb.edu if
you would like to have access to the drive.

The easiest way to use the RAW_DATA in the development suite is download Google File Stream
for your architecture at (unfortunately not yet available for Linux): 

https://support.google.com/drive/answer/7329379?hl=en

File Stream is a Dropbox like application that syncs your Google Drive with a local directory
on your machine. Then once you are added to the PypeIt team drive you will be able to access
the development suite data at the path: 

/Volumes/GoogleDrive/Team\ Drives/PHYS-GP-Hennawi/PypeIt/PypeIt-development-suite/RAW_DATA

If you don't want to deal with File Stream, then just download the RAW_DATA directory
above. However note that you will need to re-download it or parts of it if new datasets are
added (whereas File Stream would just sync with it automatically).

Once you have the RAW_DATA directory on your machine you need to set the PYPEIT_DEV environment
variable to this location, i.e. if you use c-shell in your .cshrc add

setenv PYPEIT_DEV /Volumes/GoogleDrive/Team\ Drives/PHYS-GP-Hennawi/PypeIt/PypeIt-development-suite/RAW_DATA/

or for bash

env PYPEIT_DEV=/Volumes/GoogleDrive/Team\ Drives/PHYS-GP-Hennawi/PypeIt/PypeIt-development-suite/RAW_DATA/

The pypeit_test script should now run and will reduce the raw data in this RAW_DATA directory. 