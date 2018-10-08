
""" For the development and testing of 2D ARCS
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

# General imports
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import ascii

# import sys
# from astropy.table import Table
# from astropy.io import fits
# from scipy.io import readsav
# from scipy.special import ndtr
# from IPython import embed

# PYPEIT imports
from pypeit.core import pydl

###############################################################
# Porting XIDL code x_fit2darc to python

# PURPOSE:
#  To fit the arc lines identified in x_fitarc as a fucntion of
#  their y-centroid and order number.  The main routine is in
#  x_fit2darc.  The fit is a simple least-squares with one
#  round of rejection.

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

# Order vector
order = [3, 4, 5, 6, 7, 8]

# Number of identified lines per order
npix = np.zeros_like(order)

# Read pixels and wavelengths from sv_lines_clean.txt
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
    npix[index] = len(all_pix[index])
    index = index+1

plt.figure()
for jj in np.arange(0,5):
    plt.scatter(all_pix[jj],all_wv[jj],
                label='Order {}'.format(str(order[jj])))
plt.legend()
plt.xlabel(r'pixels')
plt.ylabel(r'wavelength [$\AA$]')
plt.show()

# Now I have a dict with pixels [all_pix] , one with
# corresponding wavelengths [all_wl], and one vector with 
# the orders [order].
# I'm creating now the vectors resampling those in XIDL.

all_pix_pypeit = []
t_nrm_pypeit = []



# The following stricktly resample what is happening from
# line 135 of x_fit2darc.pro

# NORMALIZE PIX
mnx = 999999999.
mxx = -mnx
for ii in all_pix.keys():
    if mnx > np.min(all_pix[ii]):
        mnx = np.min(all_pix[ii])
    if mxx < np.max(all_pix[ii]):
        mxx = np.max(all_pix[ii])
nrm = np.array([0.5 * (mnx + mxx), mxx - mnx])
pix_nrm = {}
for ii in all_pix.keys():
    pix_nrm[ii] = 2. * (all_pix[ii] - nrm[0])/nrm[1]

# NORMALIZE ORDER
mnx = np.min(order)
mxx = np.max(order)
nrmt = np.array([0.5 * (mnx + mxx), mxx - mnx])
t_nrm = 2. * (order - nrmt[0])/nrmt[1]

# Do we actually need this?

## form XIDL   invvar = replicate(1., npix)

# SETUP THE FUNCTIONS
# these are input values of the function
# nycoeff is along pixel direction for each order. This should
# be 3 since the maximum 1-d fits for any of the orders is 3.
nycoeff = 3
# nocoeff is the order direction. 5 seems to give better rms.
nocoeff = 5

# python works with y,x instead of x,y

work2d = np.zeros((nycoeff*nocoeff,np.sum(npix)),dtype=np.float64)

worky = pydl.flegendre(pix_nrm, nycoeff)
workt = pydl.flegendre(t_nrm, nocoeff)



'''
   
   for i=0,nocoeff-1 do begin
       for j=0,nycoeff-1 do begin
           work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
       endfor
   endfor
'''

'''
from pydl
def cholesky_band(l, mininf=0.0):
    """Compute Cholesky decomposition of banded matrix.
def cholesky_solve(a, bb):
    """Solve the equation Ax=b where A is a Cholesky-banded matrix.
'''


'''

;  Setup the Functions

   ;; Do the matrix algebra
   work2di = transpose(work2d * (invvar[*] # replicate(1,nocoeff*nycoeff)))
   alpha = work2di # work2d
   beta = work2di # all_wv[*]
;   beta = work2di # (alog10(all_wv[*]))
   choldc, alpha, p, /double
   res = cholsol(alpha,p,beta, /double)
   wv_mod = dblarr(npix)
   wv_mod[*] = work2d # res

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
