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

# Reading in the output from XIDL for GNRIS.

# The relevant vectors are:
# - xidl['pix_nrm_xidl']
# - xidl['t_nrm_xidl']
xidl = ascii.read('data.xidl')
pix_nrm_xidl = np.array(xidl['pix_nrm_xidl'].data)
t_nrm_xidl = np.array(xidl['t_nrm_xidl'].data)

# These are the .sav files from Feige contiaining
# alpha, beta, p, res, work2d, and work2di for the
# 1st and the 2nd round of fitting.

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

# Order vector
order = [3, 4, 5, 6, 7, 8]

# Number of identified lines per order
pixid = np.zeros_like(order)

# Read pixels and wavelengths from sv_lines_clean.txt
# this is just reading the file that Feige gave me.
f = open('sv_lines_clean.txt', 'r')
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

# Plot to check that everything works fine
plt.figure()
for jj in np.arange(0, 5):
    plt.scatter(all_pix[jj], all_wv[jj],
                label='Order {}'.format(str(order[jj])))
plt.legend()
plt.xlabel(r'pixels')
plt.ylabel(r'wavelength [$\AA$]')
plt.show()

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


# Check that the vectors match the XIDL ones.
print(" ")
print(" pix_nrm PYPEIT vs. XIDL:")
print("  - min={}".format(np.min(pix_nrm_pypeit-pix_nrm_xidl)))
print("  - max={}".format(np.max(pix_nrm_pypeit-pix_nrm_xidl)))
print(" ")
print(" t_nrm PYPEIT vs. XIDL:")
print("  - min={}".format(np.min(t_nrm_pypeit-t_nrm_xidl)))
print("  - max={}".format(np.max(t_nrm_pypeit-t_nrm_xidl)))
print(" ")

# The following stricktly resamble what is happening from
# line 148 of x_fit2darc.pro

# pixid is a vector that cointains the number of lines id in
# each order. This means that sum(pixid) is the total number
# of line identified.
invvar = np.ones(np.sum(pixid), dtype=np.float64)
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

print(" ")
print(" The shape of work2d is: {}".format(work2d.shape))
print(" while from XIDL the shape is: {}".format(work2d_1_xidl.shape))
print(" ")
print(" The shape of worky is: {}".format(worky.shape))
print(" The shape of workt is: {}".format(workt.shape))
print(" ")

for i in range(nocoeff):
    for j in range (nycoeff):
        work2d[j*nocoeff+i,:] = worky[j,:] * workt[i,:]

print(" ")
print(" work2d PYPEIT vs. XIDL:")
print("  - min={}".format(np.min(work2d-work2d_1_xidl)))
print("  - max={}".format(np.max(work2d-work2d_1_xidl)))
print(" ")

# Do the matrix algebra
work2di = np.transpose(work2d * np.outer(np.ones(nocoeff*nycoeff,
                                         dtype=float),
                                         invvar))
print(" ")
print(" The shape of work2di is: {}".format(work2di.shape))
print(" while from XIDL the shape is: {}".format(work2di_1_xidl.shape))
print(" ")

print(" ")
print(" work2di PYPEIT vs. XIDL:")
print("  - min={}".format(np.min(work2di-work2di_1_xidl)))
print("  - max={}".format(np.max(work2di-work2di_1_xidl)))
print(" ")

# Plot to check that everything works fine
plt.figure()
plt.scatter(np.ndarray.flatten(work2d_1_xidl),
            np.ndarray.flatten((work2d_1_xidl-work2d) * 1e8),
            label='work2d')
plt.scatter(np.ndarray.flatten(work2di_1_xidl),
            np.ndarray.flatten((work2di_1_xidl-work2di) * 1e8),
            label='work2di')
plt.legend()
plt.xlabel(r'RESULTS from XIDL')
plt.ylabel(r'(XIDL - PYTHON)$\times$10$^{-8}$')
plt.show()

plt.figure()
plt.imshow((work2d_1_xidl-work2d)*1e8)
cbar = plt.colorbar()
cbar.set_label('(WORK2D XIDL - WORK2D PYTHON)x10${^-8}$', rotation=270)
plt.show()

plt.figure()
plt.imshow((work2di_1_xidl-work2di)*1e8)
cbar = plt.colorbar()
cbar.set_label('(WORK2DI XIDL - WORK2DI PYTHON)x10${^-8}$', rotation=270)
plt.show()

alpha = np.matmul(work2di.T, work2d.T)
beta = np.matmul(work2di.T, all_wv_pypeit)

# Plot to check that everything works fine
plt.figure()
plt.scatter(np.ndarray.flatten(alpha_1_xidl),
            np.ndarray.flatten((alpha_1_xidl-alpha)),
            label='alpha')
plt.legend()
plt.xlabel(r'RESULTS from XIDL')
plt.ylabel(r'(XIDL - PYTHON)')
plt.show()

plt.figure()
plt.scatter(np.ndarray.flatten(beta_1_xidl),
            np.ndarray.flatten((beta_1_xidl-beta)),
            label='beta')
plt.legend()
plt.xlabel(r'RESULTS from XIDL')
plt.ylabel(r'(XIDL - PYTHON)')
plt.show()

plt.figure()
plt.imshow(alpha_1_xidl-alpha)
cbar = plt.colorbar()
cbar.set_label('ALPHA XIDL - ALPHA PYTHON', rotation=270)
plt.show()

# Brute force test
alpha_test = np.zeros((15,15))
for row in range(15):
    for column in range(15):
        alpha_test[row,column] = np.sum(work2di[:,column]*work2d[row,:])
filt_bf = np.abs(alpha_test-alpha) > 5e-8

filt_xidl = np.abs(alpha_1_xidl-alpha) > 5e-8

print(" ")
print(" Compare np.matmul with a brute force columns times row")
print(" multiplication. This is the opertation that is expected")
print(" to be performed by IDL")
print(" Are the two operations different? {}".format(np.any(filt_bf)))
print(" Now compare the result from PYTHON with those from IDL:")
print(" Are the two operations different? {}".format(np.any(filt_xidl)))
print(" ")

alpha_test2 = np.zeros((15,15))
for row in range(15):
    for column in range(15):
        alpha_test2[column,row] = np.sum(work2di[:,row]*work2d[column,:])
filt2 = np.abs(alpha_test2-alpha_1_xidl) > 1e-8

print(" ")
print(" The shape of alpha is: {}".format(alpha.shape))
print(" while from XIDL the shape is: {}".format(alpha_1_xidl.shape))
print(" The shape of beta is: {}".format(beta.shape))
print(" while from XIDL the shape is: {}".format(beta_1_xidl.shape))
print(" ")

print(" ")
print(" alpha PYPEIT vs. XIDL:")
print("  - min={}".format(np.min(alpha-alpha_1_xidl)))
print("  - max={}".format(np.max(alpha-alpha_1_xidl)))
print(" beta PYPEIT vs. XIDL:")
print("  - min={}".format(np.min(beta-beta_1_xidl)))
print("  - max={}".format(np.max(beta-beta_1_xidl)))
print(" ")


#####     p = scipy.linalg.cholesky(alpha, lower=False)
#####     res = scipy.linalg.cho_solve((p, False), beta)
#####     
#####     print(res)
#####     print(res_1_xidl)
#####     
#####     
#####     # Plot to check that everything works fine
#####     plt.figure()
#####     plt.scatter(res_1_xidl, res)
#####     plt.xlabel(r'RESULTS from XIDL')
#####     plt.ylabel(r'RESULTS from PYTHON')
#####     plt.show()
#####     
#####################################################
#####################################################
#####################################################
#####################################################

### 
### 
### 
### res = scipy.linalg.cho_solve((p, False), beta)
### wv_mod = np.matmul(work2d, res)
### 

'''
   alpha = work2di # work2d
   beta = work2di # all_wv[*]
;   beta = work2di # (alog10(all_wv[*]))
   choldc, alpha, p, /double
   res = cholsol(alpha,p,beta, /double)
   wv_mod = dblarr(npix)
   wv_mod[*] = work2d # res
'''

# Get RMS


### inoder = [slitpix == 3]


'''
   ;; Get RMS
   gd_wv = where(invvar GT 0.0, ngd)
   msk = bytarr(npix)
   msk[gd_wv] = 1B
   ;; RESID
;   resid = (wv_mod[gd_wv] - all_wv[gd_wv]);/t[gd_wv]  ; Ang
;   fin_rms = sqrt( total( resid^2 ) / float(ngd))
;   print, 'x_fit2darc: RMS = ', fin_rms, ' Ang*Order#'
;   if keyword_set(CHKRES) then x_splot, all_wv[gd_wv], resid, psym1=1, /block
;   if keyword_set(CHKRES) then x_splot, resid, /block

   ;; REJECT
;   djs_iterstat, (wv_mod-alog10(all_wv)), sigrej=2.5, mask=msk
   djs_iterstat, (wv_mod-all_wv), sigrej=sigrej, mask=msk
   gd = where(msk EQ 1B, complement=bad, ncomplement=nbad)

   ;; RESET invvar
   if nbad NE 0 then begin
       print, 'x_fit2darc_work: Rejecting ', $
         nbad, ' of ', n_elements(all_wv), ' lines'
       invvar[bad] = 0.
   endif

   ;; Do the matrix algebra
   work2di = transpose(work2d * (invvar[*] # replicate(1,nocoeff*nycoeff)))
   alpha = work2di # work2d
   beta = work2di # (all_wv[*])
;   beta = work2di # (alog10(all_wv[*]))
   choldc, alpha, p
   res = cholsol(alpha,p,beta, /double)
   wv_mod = all_wv * 0.0
   wv_mod[*] = work2d # res


   ;; Finish
   gd_wv = where(invvar GT 0.0, ngd)
   resid = (wv_mod[gd_wv] - all_wv[gd_wv]);/t[gd_wv]  ; Ang
;   resid = 10^wv_mod[gd_wv] - all_wv[gd_wv]
   fin_rms = sqrt( total( resid^2 ) / float(ngd))
   print, 'x_fit2darc: RMS = ', fin_rms, ' Ang*Order#'
;   if keyword_set(CHKRES) then x_splot, all_wv[gd_wv], resid, psym1=1, /block
   if keyword_set(CHKRES) then x_splot, resid, /block

   ;; OUTPUT Structure
   out_str = { $
               nrm: nrm, $
               nrmt: nrmt, $
               ny: nycoeff, $
               no: nocoeff, $
               res: res }
'''
