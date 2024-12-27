""" For the development and testing of 2D ARCS
"""

from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

# General imports
import numpy as np
import scipy
import json
import matplotlib.pyplot as plt
from scipy.io import readsav
from astropy.io import ascii
from astropy.stats import sigma_clip

# PYPEIT imports
from pypeit.core import pydl
from pypeit.core import arc

###############################################################
# Porting XIDL code x_fit2darc to python

# This version: Ema Farina 2018.11.23

# PURPOSE of the XIDL code:
#  To fit the arc lines identified in x_fitarc as a function of
#  their y-centroid and order number. The main routine is in
#  x_fit2darc. The fit is a simple least-squares with one round 
#  of rejection.
#
# Feige Wang runned the code on his GNRIS data. These will be used
# test that the PYPEIT code will arrive to the same outputs of
# XIDL
#
# Test runs on GNRIS and NIRES. There is some hacking to make 
# to fix the XIDL format PypeIt compatible.

debug = False

# spec = 'GNRIS'
spec = 'NIRES'

###########################################################################
###########################################################################
###########################################################################

if spec is 'GNRIS':

    print("Test for GNRIS spectrograph")

    # Reading in the output from XIDL for GNRIS.
    # Order vector
    order = [3, 4, 5, 6, 7, 8]

    # Number of identified lines per order
    pixid = np.zeros_like(order)

    # Read pixels and wavelengths from sv_lines_clean.txt
    # this is just reading the file that Feige gave me.
    f = open('./sav_files/sv_lines_clean.txt', 'r')
    PIXWL_str = f.readlines()
    full = {}
    index = 0
    for line in PIXWL_str:
        full[index] = np.fromstring(line, dtype=float, sep=' ')
        index = index+1
    all_pix = {}
    all_wv = {}
    index = 0
    for ii in np.arange(0, 12, 2):
        all_pix[index] = full[ii][np.nonzero(full[ii])]
        all_wv[index] = full[ii+1][np.nonzero(full[ii+1])]
        pixid[index] = len(all_pix[index])
        index = index+1

    # Now I have a dict with pixels [all_pix] , one with
    # corresponding wavelengths [all_wl], and one vector with 
    # the orders [order].
    # I'm creating now the vectors resambling those in XIDL.

    all_pix_pypeit = []
    t_pypeit = []
    all_wv_pypeit = []
    npix_pypeit = []
    index = 0

    for ii in all_pix.keys():
        all_pix_pypeit = np.concatenate((all_pix_pypeit,
                                         np.array(all_pix[ii])))
        t_tmp = np.full_like(np.array(all_pix[ii]), np.float(order[index]))
        t_pypeit = np.concatenate((t_pypeit, t_tmp))
        all_wv_pypeit = np.concatenate((all_wv_pypeit,
                                         np.array(all_wv[ii])))
        npix_tmp = np.full_like(np.array(all_pix[ii]), np.size(all_pix[ii]))
        npix_pypeit = np.concatenate((npix_pypeit, npix_tmp))
        index = index + 1

    # Not the real number but a good approximation
    nspec =int(np.max(all_pix_pypeit))

###########################################################################
###########################################################################
###########################################################################

if spec is 'NIRES':

    print("Test for NIRES spectrograph")

    # Reading in the output from XIDL for NIRES.
    # Order vector
    order = [3, 4, 5, 6, 7]

    # Number of identified lines per order
    pixid = np.zeros_like(order)

    # Read pixels and wavelengths from sv_lines_clean.txt
    # this is just reading the jason file that Feige gave me.
    

    with open('./json_files/nires_wavecalib.json') as json_nires:
        nires_data = json.load(json_nires)

    all_pix_pypeit = []
    t_pypeit = []
    all_wv_pypeit = []
    npix_pypeit = []
    index = 0

    nspec =int(nires_data[str(0)]['xnorm'])

    for ii in nires_data.keys():
        if ii != 'arcparam' and ii != 'steps':
            all_pix_pypeit = np.concatenate((all_pix_pypeit,
                                             np.array(nires_data[ii]['xfit'])*(np.array(nires_data[ii]['xnorm'])-1.)))
            t_tmp = np.full_like(np.array(nires_data[ii]['xfit']), np.float(nires_data[ii]['norder']))
            # t_tmp = np.full_like(np.array(nires_data[ii]['xfit']), np.float(ii)+1.)
            t_pypeit = np.concatenate((t_pypeit, t_tmp))
            all_wv_pypeit = np.concatenate((all_wv_pypeit,
                                             np.array(nires_data[ii]['yfit'])))

###########################################################################
###########################################################################
###########################################################################


# Run the actual fit:
fit_dict = arc.fit2darc(all_wv_pypeit,
                        all_pix_pypeit,
                        t_pypeit,
                        nspec,
                        debug=debug)

# Plots
arc.fit2darc_global_qa(fit_dict)
arc.fit2darc_orders_qa(fit_dict)


