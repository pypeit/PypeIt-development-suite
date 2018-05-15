# PYPIT-development-suite v0.4.0

The Python Spectroscopic Data Reduction Pipeline development suite

This repository contains a suite of test data and PYPIT input files. 
This suite includes several examples of the PYPIT data reduction process 
for some of the supported instruments.

If you have PYPIT installed, you should be able to run
the tests from this directory with:

./pypit_test

Current allowed arguments include all, kast, isis, deimos, wht, and lris.

Many of the datasets are included as part of the Repository,
but several exceed GitHub capacity.  To include testing on
those, you will need to download them.

* DEIMOS
    1. Download the tarball from here: https://drive.google.com/open?id=1frZk1VsF9PMvwkMQEGaKtwkvD5gRXJBx
    1. cd RAW_DATA
    1. Unpack the tarball here
    1. Now the tests for DEIMOS will be enabled
* Cooked files [i.e. at least partly processed by PYPIT]
    1. Download the tarball from here: https://drive.google.com/open?id=1ZyZrk58N2peFvA7n_EmWIb42mmf33OIz
    1. cd Cooked
    1. Unpack the tarball here


Use -outputdir to choose a different output directory.
