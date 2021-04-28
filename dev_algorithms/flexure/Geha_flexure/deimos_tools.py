import numpy as np
import os

from astropy import table
from astropy.table import Table, Column
from astropy.io import ascii
from astropy.time import Time


from astropy import units as u
from astropy.coordinates import SkyCoord,EarthLocation

DEIMOS_DROPBOX = '/Users/mgeha/Dropbox/DEIMOS/'


######################################################
# HELIOCENTRIC CORRECTION PER MASK
def zspec_helio(zspec):
    i=0
    t = Time(zspec['MJD'],format='mjd')
    m=zspec['RA'] != 0
    r = np.median(zspec['RA'][m])
    d = np.median(zspec['DEC'][m])
    sc = SkyCoord(r,d, unit=(u.deg,u.deg))

    keck = EarthLocation.from_geodetic(lat=19.8283*u.deg, lon=-155.4783*u.deg, height=4160*u.m)
    heliocorr = sc.radial_velocity_correction('heliocentric',obstime=t, location=keck)

    vhelio = heliocorr.to(u.km/u.s) * (u.s/u.km)

    return vhelio

######################################################
# READ DEIMOS OBJECT GOOGLE DOCUMENT
def deimos_google():
    key = '1V2aVg1QghpQ70Lms40zNUjcCrycBF2bjgs-mrp6ojI8'
    gid=1906496323
    url = 'https://docs.google.com/spreadsheets/d/{0}/export?format=csv&gid={1}'.format(key, gid)
    masklist = Table.read(url, format='csv')

    gid =0
    url = 'https://docs.google.com/spreadsheets/d/{0}/export?format=csv&gid={1}'.format(key, gid)
    objlist = ascii.read(url, format='csv')
    
    return objlist,masklist

######################################################
# READ DEIMOS OBJECT GOOGLE DOCUMENT
def deimos_google_M31():
    key = '1V2aVg1QghpQ70Lms40zNUjcCrycBF2bjgs-mrp6ojI8'
    gid=237917397
    url = 'https://docs.google.com/spreadsheets/d/{0}/export?format=csv&gid={1}'.format(key, gid)
    masklist = Table.read(url, format='csv')

    gid =180140955
    url = 'https://docs.google.com/spreadsheets/d/{0}/export?format=csv&gid={1}'.format(key, gid)
    objlist = ascii.read(url, format='csv')
    
    return objlist,masklist


######################################################
# FIX BAD RA/DEC VALUES IN ZSPEC FILES
def fix_zspec(zspec_file):
  

    z=Table.read(zspec_file)
    #print(z['DEC'])
    dnew = []
    for dec in z['DEC']:
        d=dec.split(':')
        print(dec)
        if (d[0].strip() != '0.0'):
            if (int(d[1]) > 40):
                d[0]='-31'
            if (int(d[1]) < 30):
                d[0]='-32'
            dnew.append(d[0]+':'+d[1]+':'+d[2])

        else:
            dnew.append(dec)
            
    z['DEC'] = dnew

    z.write('new.fits',format='fits')
    

######################################################
# FILL A COLUMN
def filled_column(name, fill_value, size):
    """
    Tool to allow for large strings
    """
    return Column([fill_value]*int(size), name)


######################################################
def calc_na_EW(spec,wvl,ivar,z,plot=False):

    wave = wvl / (1.0+z)  

    # 21AA window centered on 8190AA
    wline = [8178., 8200.5]
    wred  = [8203., 8230.]
    wblue = [8155., 8175.]
    waver = (wred[0] + wred[1])/2.
    waveb = (wblue[0] + wblue[1])/2.

    mnaI  = (wave > wline[0]) & (wave < wline[1])
    mred = (wave > wred[0]) & (wave < wred[1])
    mblue = (wave > wblue[0]) & (wave < wblue[1])

    # INITIALIZE QUANTITIES
    na_EW = 0
    na_EW_err = -99

    # DETERMINE WEIGHTED MEAN OF BLUE/RED PSEUDO-CONTINUUM BAND
    # DON"T CALCULATE IF DATA DOESN"T EXIST
    if (np.sum(mblue) != 0) & (np.sum(mred) != 0): 
        sum1 = np.sum(spec[mblue] * ivar[mblue]**2 )
        sum2 = np.sum( ivar[mblue]**2 )
        bcont = sum1 / sum2

        sum1 = np.sum(spec[mred] * ivar[mred]**2 )
        sum2 = np.sum( ivar[mred]**2 )
        rcont = sum1 / sum2


        # DEFINE CONTINUUM LINE BETWEEN RED/BLUE PASSBAND (y=mx+b)
        mline = (rcont - bcont) / (waver - waveb)
        bline = rcont - (mline * waver)
        continuum = (mline * wave) + bline

        
        # CALCULATE DELTA LAMBAs, ASSUME NON-LINEAR BINNING
        dlambda = np.zeros(np.size(wave))
        for i,item in enumerate(wave):
            if (i != np.size(wave)-1):
                dlambda[i] = wave[i+1]-wave[i]            
                
        # CALCULATE NaI EW
        na_EW = np.sum((1 - (spec[mnaI] / continuum[mnaI])) * dlambda[mnaI])
        na_EW_var = np.sum(1./np.sum(ivar[mnaI] * (continuum[mnaI]/dlambda[mnaI])))
        na_EW_err = np.sqrt(na_EW_var)
   

        # PLOT IF REQUESTED
        if (plot == 1):
            plt.figure()
            na = [8183,8195]
            plt.axvline(na[0],color='r')
            plt.axvline(na[1],color='r')
            plt.title(na_EW)
            plt.plot(wave,spec)
            plt.plot(wave,continuum)
            plt.xlim(8170,8230)

    return na_EW, na_EW_err


######################################################
# ROUGH MEMEBERSHIP
def membership_CMD(zspec,obj):

    # GET ISOCHRONE PROPERTIES    
    EBV = obj['EBV_SF11']
    dist= obj['Dist_kpc']
    iso = obj['iso_guess']

    r,gr = plot_isochrone_padova(dist*1e3,EBV,iso)

    r_hb,gr_hb = plot_isochrone_HB(dist*1e3,EBV)    

        
    dmin = []
    emin = []
    for star in zspec:
        err = np.sqrt(star['GMAG_ERR']**2 + star['RMAG_ERR']**2)

        d = (r - star['RMAG'])**2 + (gr - (star['GMAG'] - star['RMAG']))**2
        d2 = (r_hb - star['RMAG'])**2 + (gr_hb - (star['GMAG'] - star['RMAG']))**2

        tmp = np.min(d)
        tmp2=np.min(d2)        
        if tmp2 < tmp:
            tmp=tmp2
        dmin.append(tmp)  
        emin.append(err)
    
    emin = np.array(emin)
    m    =  (emin > 1) 
    emin[m] = 0
    m    =  (emin > 0.3) & (star['RMAG_ERR'] < 22)   #. account for bad errors in UMa2    
    emin[m] = 0
    
    #*********************************
    # CMD THRESHOLD == 0.1 PLUS ERRORS
    mem = np.array(dmin) < 0.1 + emin   # SET THRESHOLD PLUS PHOTOMETRIC ERROR
    
    return mem



######################################################
def plot_isochrone_padova(dist,EBV,iso):

    iso = ascii.read(DEIMOS_DROPBOX+'/Photometry/isochrones/iso_t12_z'+str(iso)+'.dat')

    #A(g)/(E(B-V)) = 3.793    
    Ag = 3.793 * EBV
    Ar = 2.751 * EBV

    r_iso = iso['rmag'] + 5.*np.log10(dist) - 5. + Ar
    gr_iso = iso['gmag'] - iso['rmag'] + (Ag - Ar)
    
    return r_iso,gr_iso


######################################################
def plot_isochrone_HB(dist,EBV):

    iso = Table.read(DEIMOS_DROPBOX+'/Photometry/isochrones/M92_fiducial.dat',format='ascii',guess=False)
    hb = iso['typ'] == 1
    
    #A(g)/(E(B-V)) = 3.793    
    Ag = 3.793 * EBV
    Ar = 2.751 * EBV

    r_iso = iso['rmag'][hb] + 5.*np.log10(dist) - 5. + Ar
    gr_iso = iso['gmr'][hb] + (Ag - Ar) + 0.2
    
    return r_iso,gr_iso 
