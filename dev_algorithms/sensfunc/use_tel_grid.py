import numpy as np
from astropy.io import fits

def read_telluric_grid(filename):
    hdul = fits.open(filename)
    wave_grid = hdul[1].data
    model_grid = hdul[0].data

    pg = hdul[0].header['PRES0']+hdul[0].header['DPRES']*np.arange(0,hdul[0].header['NPRES'])
    tg = hdul[0].header['TEMP0']+hdul[0].header['DTEMP']*np.arange(0,hdul[0].header['NTEMP'])
    hg = hdul[0].header['HUM0']+hdul[0].header['DHUM']*np.arange(0,hdul[0].header['NHUM'])
    if hdul[0].header['NAM'] > 1:
        ag = hdul[0].header['AM0']+hdul[0].header['DAM']*np.arange(0,hdul[0].header['NAM'])
    else:
        ag = hdul[0].header['AM0']+1*np.arange(0,1)

    return wave_grid, model_grid, pg, tg, hg, ag

def interp_telluric_grid(theta,pg,tg,hg,ag,model_grid):

    press,temp,hum,airmass = theta
    if len(pg) > 1:
        p_ind = np.mod(int(np.round((press-pg[0])/(pg[1]-pg[0]))))
    else:
        p_ind = 0
    if len(tg) > 1:
        t_ind = np.mod(int(np.round((temp-tg[0])/(tg[1]-tg[0]))))
    else:
        t_ind = 0
    if len(hg) > 1:
        h_ind = np.mod(int(np.round((hum-hg[0])/(hg[1]-hg[0]))))
    else:
        h_ind = 0
    if len(ag) > 1:
        a_ind = np.mod(int(np.round((airmass-ag[0])/(ag[1]-ag[0]))))
    else:
        a_ind = 0

    return model_grid[p_ind,t_ind,h_ind,a_ind]



