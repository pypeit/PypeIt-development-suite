""" For the development and testing of 2D ARCS
"""

from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

# General imports
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.io import readsav
from astropy.io import ascii
from astropy.stats import sigma_clip

# PYPEIT imports
from pypeit.core import pydl

###############################################################
# Porting XIDL code x_fit2darc to python

# PURPOSE of the XIDL code:
#  To fit the arc lines identified in x_fitarc as a function of
#  their y-centroid and order number. The main routine is in
#  x_fit2darc. The fit is a simple least-squares with one round 
#  of rejection.

# Feige runned the code on his GNRIS data. I will use this to
# test that the PYPEIT code will arrive to the same outputs of
# XIDL

debug = False
debug = True

# Reading in the output from XIDL for GNRIS.

# The relevant vectors are:
# - xidl['pix_nrm_xidl']
# - xidl['t_nrm_xidl']
xidl = ascii.read('./sav_files/data.xidl')
pix_nrm_xidl = np.array(xidl['pix_nrm_xidl'].data)
t_nrm_xidl = np.array(xidl['t_nrm_xidl'].data)

# These are the .sav files from Feige contiaining
# alpha, beta, p, res, work2d, work2dim, and wv_mod
# for the 1st and the 2nd round of fitting.

alpha_1_dict = readsav('./sav_files/alpha1.sav', python_dict=True)
alpha_1_xidl = alpha_1_dict['alpha']
alpha_2_dict = readsav('./sav_files/alpha2.sav', python_dict=True)
alpha_2_xidl = alpha_2_dict['alpha']

beta_1_dict = readsav('./sav_files/beta1.sav', python_dict=True)
beta_1_xidl = beta_1_dict['beta']
beta_2_dict = readsav('./sav_files/beta2.sav', python_dict=True)
beta_2_xidl = beta_2_dict['beta']

p_1_dict = readsav('./sav_files/p1.sav', python_dict=True)
p_1_xidl = p_1_dict['p']
p_2_dict = readsav('./sav_files/p2.sav', python_dict=True)
p_2_xidl = p_2_dict['p']

res_1_dict = readsav('./sav_files/res1.sav', python_dict=True)
res_1_xidl = res_1_dict['res']
res_2_dict = readsav('./sav_files/res2.sav', python_dict=True)
res_2_xidl = res_2_dict['res']

work2d_1_dict = readsav('./sav_files/work2d1.sav', python_dict=True)
work2d_1_xidl = work2d_1_dict['work2d']
work2d_2_dict = readsav('./sav_files/work2d2.sav', python_dict=True)
work2d_2_xidl = work2d_2_dict['work2d']

work2di_1_dict = readsav('./sav_files/work2di1.sav', python_dict=True)
work2di_1_xidl = work2di_1_dict['work2di']
work2di_2_dict = readsav('./sav_files/work2di2.sav', python_dict=True)
work2di_2_xidl = work2di_2_dict['work2di']

wv_mod_1_dict = readsav('./sav_files/wv_mod_1.sav', python_dict=True)
wv_mod_1_xidl = wv_mod_1_dict['wv_mod']
wv_mod_2_dict = readsav('./sav_files/wv_mod_2.sav', python_dict=True)
wv_mod_2_xidl = wv_mod_2_dict['wv_mod']

# I also load the vectors t and all_wv
t_dict = readsav('./sav_files/t.sav', python_dict=True)
t_xidl = t_dict['t']

all_wv_dict = readsav('./sav_files/all_wv.sav', python_dict=True)
all_wv_xidl = all_wv_dict['all_wv']

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

# Setting the same format of XIDL
all_wv_pypeit = all_wv_pypeit * t_pypeit

# NORMALIZE PIX
mnx = np.min(all_pix_pypeit)
mxx = np.max(all_pix_pypeit)
nrm = np.array([0.5 * (mnx + mxx), mxx - mnx])
pix_nrm_pypeit = 2. * (all_pix_pypeit - nrm[0])/nrm[1]

# NORMALIZE ORDERS
mnx = np.min(t_pypeit)
mxx = np.max(t_pypeit)
nrmt = np.array([0.5 * (mnx + mxx), mxx - mnx])
t_nrm_pypeit = 2. * (t_pypeit - nrmt[0])/nrmt[1]


if debug: 
    # Check that the vectors match the XIDL ones.
    print(" ")
    print(" t PYPEIT vs. XIDL:")
    print("  - min={}".format(np.min(t_pypeit-t_xidl)))
    print("  - max={}".format(np.max(t_pypeit-t_xidl)))
    print(" ")
    print(" all_wv PYPEIT vs. XIDL:")
    print("  - min={}".format(np.min(all_wv_pypeit-all_wv_xidl)))
    print("  - max={}".format(np.max(all_wv_pypeit-all_wv_xidl)))
    print(" ")
    print(" pix_nrm PYPEIT vs. XIDL:")
    print("  - min={}".format(np.min(pix_nrm_pypeit-pix_nrm_xidl)))
    print("  - max={}".format(np.max(pix_nrm_pypeit-pix_nrm_xidl)))
    print(" ")
    print(" t_nrm PYPEIT vs. XIDL:")
    print("  - min={}".format(np.min(t_nrm_pypeit-t_nrm_xidl)))
    print("  - max={}".format(np.max(t_nrm_pypeit-t_nrm_xidl)))
    print(" ")
    # Plot pixels vs. order with wavelengths as colorbar
    print(" ")
    print(" Plot of wavelengths as a function of normalized orders and")
    print(" pixels.")
    print(" ")
    plt.figure()
    cm = plt.cm.get_cmap('RdYlBu_r')
    sc = plt.scatter(t_nrm_pypeit, pix_nrm_pypeit,
                     c=all_wv_pypeit/t_pypeit/10000., cmap=cm)
    cbar = plt.colorbar(sc)
    cbar.set_label(r'WAVELENGTHS [$\mu$m]', rotation=270,
                   labelpad=20)
    plt.xlabel(r'NORMALIZED ORDERS')
    plt.ylabel(r'NORMALIZED PIX')
    plt.show()

# The following stricktly resamble what is happening from
# line 148 of x_fit2darc.pro

# pixid is a vector that cointains the number of lines id in
# each order. This means that sum(pixid) is the total number
# of line identified.
invvar = np.ones(np.sum(pixid), dtype=np.float64)
if debug: 
    print(" ")
    print(" The shape of invvar is: {}".format(invvar.shape))
    print(" ")

# SETUP THE FUNCTIONS
# these are input values of the function. nycoeff is along pixel
# direction for each order. This should be 3 since the maximum 
# 1-d fits for any of the orders is 3.
nycoeff = 3
# nocoeff is the order direction. 5 seems to give better rms.
nocoeff = 5

work2d = np.zeros((nycoeff*nocoeff, np.sum(pixid)), dtype=np.float64)
worky = pydl.flegendre(pix_nrm_pypeit, nycoeff)
workt = pydl.flegendre(t_nrm_pypeit, nocoeff)

for i in range(nocoeff):
    for j in range(nycoeff):
        work2d[j*nocoeff+i,:] = worky[j,:] * workt[i,:]


if debug: 
    print(" ")
    print(" The shape of work2d is: {}".format(work2d.shape))
    print(" while from XIDL the shape is: {}".format(work2d_1_xidl.shape))
    print(" ")
    print(" The shape of worky is: {}".format(worky.shape))
    print(" The shape of workt is: {}".format(workt.shape))
    print(" ")
    print(" work2d PYPEIT vs. XIDL:")
    print("  - min={}".format(np.min(work2d-work2d_1_xidl)))
    print("  - max={}".format(np.max(work2d-work2d_1_xidl)))
    print(" ")

# Do the matrix algebra
work2di = np.transpose(work2d * np.outer(np.ones(nocoeff*nycoeff,
                                         dtype=float),
                                         invvar))
if debug: 
    print(" ")
    print(" The shape of work2di is: {}".format(work2di.shape))
    print(" while from XIDL the shape is: {}".format(work2di_1_xidl.shape))
    print(" ")
    print(" The shape of all_wv_pypeit.T is: {}".format(all_wv_pypeit.T.shape))
    print(" while from XIDL the shape of all_wv_xidl is: {}".format(all_wv_xidl.shape))
    print(" ")
    print(" work2di PYPEIT vs. XIDL:")
    print("  - min={}".format(np.min(work2di-work2di_1_xidl)))
    print("  - max={}".format(np.max(work2di-work2di_1_xidl)))
    print(" ")

    # Plot to check that everything works fine
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    im1 = ax1.imshow((work2d_1_xidl-work2d)*1e8)
    cbar1 = fig.colorbar(im1)
    cbar1.set_label(r'(WORK2D XIDL - WORK2D PYTHON)$\times$10$^{-8}$',
                     rotation=270, labelpad=20)
    ax1.set_aspect('auto')
    ax2 = fig.add_subplot(122)
    ax2.hist(np.ndarray.flatten(work2d_1_xidl-work2d)*1e8)
    ax2.set_xlabel(r'(WORK2D XIDL - WORK2D PYTHON)$\times$10$^{-8}$')
    ax2.set_aspect('auto')
    plt.show()
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    im1 = ax1.imshow((work2di_1_xidl-work2di)*1e8)
    cbar1 = fig.colorbar(im1)
    cbar1.set_label(r'(WORK2DI XIDL - WORK2DI PYTHON)$\times$10$^{-8}$',
                     rotation=270, labelpad=20)
    ax1.set_aspect('auto')
    ax2 = fig.add_subplot(122)
    ax2.hist(np.ndarray.flatten(work2di_1_xidl-work2di)*1e8)
    ax2.set_xlabel(r'(WORK2DI XIDL - WORK2DI PYTHON)$\times$10$^{-8}$')
    ax2.set_aspect('auto')
    plt.show()

alpha = work2d.dot(work2di)
beta = all_wv_pypeit.dot(work2di)

solve = np.linalg.solve(alpha,beta)
wv_mod = solve.dot(work2d)

if debug: 
    plt.figure()
    plt.scatter(all_wv_pypeit / t_pypeit,
                wv_mod / t_pypeit - all_wv_pypeit / t_pypeit,
                marker="v",
                label='PypeIt')
    plt.scatter(all_wv_xidl / t_xidl,
                wv_mod_1_xidl/t_xidl - all_wv_xidl / t_xidl,
                marker="^",
                label='XIDL')
    plt.legend()
    plt.title(r'1st ITERATION')
    plt.xlabel(r'WAVELENGTHS [$\AA$]')
    plt.ylabel(r'RESIDUALS [$\AA$]')
    plt.show()

# Mask Values
# msk = True means a bad value
msk = sigma_clip(wv_mod-all_wv_pypeit).mask
if np.any(msk):
    print("Rejecting: {} of {} lines.".format(len(msk[np.where(msk == True)]),len(msk)))
    invvar[msk] = 0.

# Do the matrix algebra
work2di = np.transpose(work2d * np.outer(np.ones(nocoeff*nycoeff,
                                         dtype=float),
                                         invvar))
alpha = work2d.dot(work2di)
beta = all_wv_pypeit.dot(work2di)

solve = np.linalg.solve(alpha,beta)
wv_mod = solve.dot(work2d)

plt.figure()
plt.scatter(all_wv_pypeit / t_pypeit,
            wv_mod / t_pypeit - all_wv_pypeit / t_pypeit,
            marker="v",
            label='PypeIt')
plt.scatter(all_wv_xidl / t_xidl,
            wv_mod_2_xidl/t_xidl - all_wv_xidl / t_xidl,
            marker="^",
            label='XIDL')
plt.legend()
plt.title(r'2nd ITERATION')
plt.xlabel(r'WAVELENGTHS [$\AA$]')
plt.ylabel(r'RESIDUALS [$\AA$]')
plt.show()

# Finsih
gd_wv = invvar > 0.
resid = (wv_mod[gd_wv] / t_pypeit[gd_wv] - all_wv_pypeit[gd_wv] / t_pypeit[gd_wv])
fin_rms = np.sqrt(np.mean(resid**2))

print("RMS: {} Ang*Order#.".format(fin_rms))

plt.figure()
cm = plt.cm.get_cmap('RdYlBu_r')
sc = plt.scatter(t_pypeit, all_pix_pypeit,
                 c=wv_mod/t_pypeit/10000., cmap=cm)
cbar = plt.colorbar(sc)
cbar.set_label(r'WAVELENGTHS [$\mu$m]', rotation=270,
               labelpad=20)
plt.title(r'RESULT OF THE FIT')
plt.xlabel(r'ORDERS')
plt.ylabel(r'PIXELS')
plt.show()

