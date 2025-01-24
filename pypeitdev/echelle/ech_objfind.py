


import numpy as np
from pypeit.core import extract
from pypeit import msgs
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from matplotlib import pyplot as plt

from pypeit import ginga
from pypeit.core import pydl
from pypeit import utils
from pypeit.core import pixels
from sklearn.decomposition import PCA
from pypeit import specobjs
from pypeit.core import extract
from astropy.stats import SigmaClip
#from pydl.pydlutils.spheregroup import spheregroup
from pypeit.core.pydl import spheregroup

# HIRES
spectro = 'MAGE'
user = '/Users/joe/'
if spectro == 'HIRES':
    hdu = fits.open(user + '/Dropbox/hires_fndobj/f_hires0184B.fits.gz')
    objminsky =hdu[2].data
    ivar  = hdu[1].data
    mask = (ivar > 0.0)
    order_str = Table.read(user + '/Dropbox/hires_fndobj/OStr_B_04.fits')
    slit_left = (order_str['LHEDG']).T
    slit_righ = (order_str['RHEDG']).T
    plate_scale = 0.36
    std_trace = None
elif spectro == 'MAGE':
    from scipy.io import readsav
    hdu = fits.open(user + 'Dropbox/hires_fndobj/MAGE/f_mage2013.fits.gz')
    objminsky = hdu[0].data - hdu[2].data
    ivar  = hdu[4].data
    #ivar = utils.calc_ivar(var)
    mask = (ivar > 0.0)
    slit_left = readsav(user + 'Dropbox/hires_fndobj/MAGE/left_edge.sav',python_dict=False)['left_edge'].T
    slit_righ = readsav(user + 'Dropbox/hires_fndobj/MAGE/right_edge.sav',python_dict=False)['right_edge'].T
    plate_scale = 0.3
    hdu_std = fits.open(user + 'Dropbox/hires_fndobj/MAGE/ObjStr1046.fits')
    std_trace = hdu_std[1].data['trace'].T
    # Standard is the ObjStr1046.fits exten = 1
elif spectro == 'ESI':
    # ESI
    hdu = fits.open(user + '/Dropbox/Cowie_2002-02-17/Final/fringe_ES.20020217.35453.fits.gz')
    sciimg = hdu[0].data
    var  = hdu[1].data
    ivar = utils.calc_ivar(var)
    mask = (var > 0.0)
    skyimg = hdu[2].data
    objminsky = sciimg - skyimg
    hdu_sedg = fits.open(user + '/Dropbox/Cowie_2002-02-17/Flats/SEdgECH75_1x1.fits')
    data = hdu_sedg[0].data
    slit_left = data[0,:,:].T
    slit_righ = data[1,:,:].T
    plate_scale = 0.149
    hdu_std =fits.open(user + '/Dropbox/Cowie_2002-02-17/Extract/Obj_ES.20020217.20037.fits.gz')
    std_trace = hdu_std[1].data['trace'][:,0:slit_left.shape[0]].T
    std_trace[:,-1] = std_trace[:,-2]+ (np.median(std_trace[:200,-1])-np.median(std_trace[:200,-2]))
    #The structure is exten = 1, and the xpos,ypos are the standard trace
elif spectro == 'NIRES':
    from linetools import utils as ltu
    jdict = ltu.loadjson(user + '/Dropbox/hires_fndobj/tilt_nires.json')
    slit_left = np.array(jdict['lcen'])
    slit_righ = np.array(jdict['rcen'])
    hdu = fits.open(user + 'Dropbox/hires_fndobj/spec2d_J1724+1901_NIRES_2018Jun04T130207.856.fits')
    objminsky = hdu[1].data - hdu[3].data
    ivar = hdu[2].data
    mask = (ivar>0)
    plate_scale = 0.123
    std_trace = None
elif spectro == 'GNIRS':
    from scipy.io import readsav
    #hdu = fits.open(user + '/Dropbox/hires_fndobj/sci-N20170331S0216-219.fits')
    #hdu = fits.open(user + '/Dropbox/hires_fndobj/GNIRS/J021514.76+004223.8/Science/J021514.76+004223.8_1/sci-N20170927S0294-297.fits')
    #hdu = fits.open(user + 'Dropbox/hires_fndobj/GNIRS/J005424.45+004750.2/Science/J005424.45+004750.2_7/sci-N20171021S0264-267.fits')
    hdu = fits.open(user + 'Dropbox/hires_fndobj/GNIRS/J002407.02-001237.2/Science/J002407.02-001237.2_5/sci-N20171006S0236-239.fits')
    obj = hdu[0].data
    #objminsky = obj - hdu[1].data
    objminsky = hdu[1].data - obj # test negative trace
    ivar  = hdu[2].data
    #ivar = utils.calc_ivar(var)
    mask = (ivar > 0.0)
    #slit_left = readsav(user + '/Dropbox/hires_fndobj/left_edge.sav', python_dict=False)['left_edge'].T
    #slit_righ = readsav(user + '/Dropbox/hires_fndobj/right_edge.sav', python_dict=False)['right_edge'].T
    #slit_left = readsav(user + '/Dropbox/hires_fndobj/GNIRS/J021514.76+004223.8/left_edge_J0215.sav',python_dict=False)['left_edge'].T
    #slit_righ = readsav(user + '/Dropbox/hires_fndobj/GNIRS/J021514.76+004223.8/right_edge_J0215.sav',python_dict=False)['right_edge'].T
    #slit_left = readsav(user + '/Dropbox/hires_fndobj/GNIRS/J005424.45+004750.2/left_edge_J0054.sav',python_dict=False)['left_edge'].T
    #slit_righ = readsav(user + '/Dropbox/hires_fndobj/GNIRS/J005424.45+004750.2/right_edge_J0054.sav',python_dict=False)['right_edge'].T
    slit_left = readsav(user + 'Dropbox/hires_fndobj/GNIRS/J002407.02-001237.2/left_edge_J0024.sav',python_dict=False)['left_edge'].T
    slit_righ = readsav(user + 'Dropbox/hires_fndobj/GNIRS/J002407.02-001237.2/right_edge_J0024.sav',python_dict=False)['right_edge'].T
    plate_scale = 0.15
    #Use std. Telluric directory. Extension 4 is the object structure, and you want xpos and ypos. Xpos is the spatial trace that you want
    # to pass in. ypos is just np.arange(nspec)
    hdu_std =fits.open(user + 'Dropbox/hires_fndobj/GNIRS/J002407.02-001237.2/Science/HIP117774_5/tel-N20171006S0228-231.fits')
    std_xpos = hdu_std[4].data['XPOS'].T
    std_trace = np.zeros(np.shape(slit_left))
    for i in range(6):
        std_trace[:,i] = std_xpos[:,2*i]


ordermask = pixels.slit_pixels(slit_left, slit_righ, objminsky.shape[1])

#viewer, ch = ginga.show_image(objminsky)
#ginga.show_slits(viewer,ch, slit_left, slit_righ)

#yvec = np.outer(np.ones(norders), spec_vec)

#tset_left = pydl.xy2traceset(yvec, slit_left.T, ncoeff=ncoeff)
#slit_left_fit = tset_left.yfit.T

#tset_righ = pydl.xy2traceset(yvec, slit_righ.T, ncoeff=ncoeff)
#slit_righ_fit = tset_righ.yfit.T

#ginga.show_slits(viewer,ch, slit_left_fit, slit_righ_fit)

slit_mid = (slit_left + slit_righ)/2.0
#usepca = np.zeros(slit_mid.shape[1],dtype=bool)
#usepca[0:8] = True
#pca_out = pca_trace(slit_mid, npca = 4, npoly_cen = 3)
#sys.exit(-1)
#viewer, ch = ginga.show_image(ordermask > 0)
#ginga.show_slits(viewer,ch, slit_mid, pca_out)

# create the ouptut images skymask and objmask
#skymask = np.zeros_like(objminsky, dtype=bool)
#objmask = np.zeros_like(objminsky, dtype=bool)
image = objminsky.copy()
inmask = mask.copy()
#Routine starts here.
#------
ncoeff = 5
box_radius = 2.0  # arcseconds
min_snr = 0.2
nabove_min_snr = 2
npca = None
pca_explained_var = 99.8
sig_thresh = 3.0
box_radius = 2.0
ncoeff = 5
trim_edg = (5,5)


sobjs_final = extract.ech_objfind(image, ivar, ordermask, slit_left, slit_righ,inmask=inmask,plate_scale=plate_scale,
                                  std_trace =  std_trace, pca_explained_var=pca_explained_var, min_snr=min_snr,
                                  nabove_min_snr=nabove_min_snr,box_radius=box_radius,sig_thresh=sig_thresh,
                                  show_peaks=False,show_fits=False,show_trace=True,debug=True)
