import numpy as np
import os

from astropy.table import Table
from astropy.io import ascii,fits

import deimos_tools
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt

import scipy.ndimage as scipynd
from scipy.optimize import curve_fit


# NEED TO REMOVE THIS SOMEHOW
DEIMOS_DROPBOX = '/Users/mgeha/Dropbox/DEIMOS/'


########################################
def gaussian(x,*p) :
    # A gaussian peak with:
    #   Constant Background          : p[0]
    #   Peak height above background : p[1]
    #   Central value                : p[2]
    #   Standard deviation           : p[3]
    return p[0]+p[1]*np.exp(-1.*(x-p[2])**2/(2.*p[3]**2))


#######################################################
# SPLIT PYPEIT SPEC!D EXTENSION NAME
#######################################################
def parse_spec1d_ext(pext_str):

        # GET DETECTOR
        tmp = pext_str.split('DET')
        det = int(tmp[1])
        
        # GET SPAT NUMBER
        tmp2 = pext_str.split('-')
        tmp3 = tmp2[0].split('T')
        spat = int(tmp3[1])
        
        # GET SLIT NUMBER
        tmp3 = tmp2[1].split('T')
        slit = int(tmp3[1])

        
        # HACKED XPOSITION ON MASK
        xslit_pos = int(spat) + 2048.*(int(det)-1.)

        return det, spat, slit, xslit_pos


#######################################################
# CALC SN
#######################################################
def calc_rb_SN(r1,b1, hdu):

    # TRY SPECTRUM AND CALCULATE SN
    try:
        rSN = np.median(hdu[r1].data['OPT_COUNTS'] * np.sqrt(hdu[r1].data['OPT_COUNTS_IVAR']))
        bSN = np.median(hdu[b1].data['OPT_COUNTS'] * np.sqrt(hdu[b1].data['OPT_COUNTS_IVAR']))
        aa=hdu[b1].data['OPT_WAVE'] * hdu[r1].data['OPT_WAVE']        
    except:
        rSN=0
        bSN=0

    return rSN, bSN                




#######################################################
#  GENERAL SCRIPT TO CALL 
def load_spectrum(single_slit,hdu,vacuum=0):

    r = single_slit['rname']
    b = single_slit['bname']

    try:

        tmp_wave = np.concatenate((hdu[b].data['OPT_WAVE'],hdu[r].data['OPT_WAVE']),axis=None)
        all_flux = np.concatenate((hdu[b].data['OPT_COUNTS'],hdu[r].data['OPT_COUNTS']),axis=None)
        all_sky = np.concatenate((hdu[b].data['OPT_COUNTS_SKY'],hdu[r].data['OPT_COUNTS_SKY']),axis=None)
        all_ivar = np.concatenate((hdu[b].data['OPT_COUNTS_IVAR'],hdu[r].data['OPT_COUNTS_IVAR']),axis=None)

        fitwave  = single_slit['fit_slope']*tmp_wave + single_slit['fit_b']
        vwave = tmp_wave - fitwave

        # CONVERT PYPEIT OUTPUT WAVELENGTHS FROM VACUUM TO AIR
        all_wave = vwave
        if (vacuum == 0):
            all_wave = vwave / (1.0 + 2.735182E-4 + 131.4182 / vwave**2 + 2.76249E8 / vwave**4)


        # TRIM ENDS
        all_wave=all_wave[5:-15]
        all_flux=all_flux[5:-15]
        all_ivar=all_ivar[5:-15]
        all_sky=all_sky[5:-15]

        # REMOVE CRAZY 500-SIGMA VALUES
        cmask = (all_flux > np.percentile(all_flux,0.1)) & (all_flux < np.percentile(all_flux,99.9))

        m=np.median(all_flux[cmask])
        s=np.std(all_flux[cmask])
        mm = (all_flux > 500.*s + m) | (all_flux < m-50.*s)
        all_flux[mm] = m
        all_ivar[mm] = 1e6
        if (np.sum(mm) > 10):
            print('Removing more than 10 pixels of data')
        
    except:
        print('no data!')
        n=8172
        all_wave=np.zeros(n)
        all_flux=np.zeros(n)
        all_ivar=np.zeros(n)
        all_sky=np.zeros(n)

    return all_wave,all_flux,all_ivar,all_sky


def load_old_templates():
    c = 299792.5
    old_tmpl = DEIMOS_DROPBOX+'/templates/idl_templates/deimos-052715.fits'
    thdu = fits.open(old_tmpl)
    h   = thdu[0].header

    log_wave = 10**(h['COEFF0']+np.arange(h['NAXIS1'])*h['COEFF1'])

    tmpl = thdu[0].data
    m = (log_wave >= 8400) & (log_wave <= 8720)
    t1  = tmpl[1,m]
    t2  = tmpl[3,m]
    t3  = tmpl[6,m]

    lwave= log_wave[m]
    vscale = np.mean([c*(lwave[i + 1] - lwave[i])/lwave[i] for i in range(len(lwave)-1)]) 

    return lwave,vscale,t1,t2,t3


def filename_telluric(h2o,o2):
    dir = '/Users/mgeha/Dropbox/DEIMOS/templates/tellurics/telluric_0.01A_h2o_'
    filename = dir+'{:0.0f}_o2_{:0.2f}_.fits'.format(h2o,o2)
    return filename

def parse_tfile(tfile):

    spl = tfile.split('_')
    h2o = np.float(spl[3])
    o2  = np.float(spl[5])
    return o2,h2o

def find_chi2_roots(values, chi2):

    z = np.polyfit(values,chi2,2)
    p = np.poly1d(z)

    p2 = np.roots(p)
    min_w    = p2[0].real
    min_chi2 = p(min_w)
    
    return min_w,min_chi2


