



from pypeit import msgs
import numpy as np
from IPython import embed
from pypeit import flatfield
from astropy.io import fits
import os
from pypeit.spectrographs.util import load_spectrograph
from pypeit import ginga
from pypeit import traceslits
from pypeit.core import parse
from pypeit.core import pixels
from pypeit import scienceimage

type = 'LRIS_red'
devpath = os.getenv('PYPEIT_DEV')

if type == 'LRIS_red':
    det = 1
    sdet = parse.get_dnum(det, prefix=False)
    rawpath = devpath + '/RAW_DATA/Keck_LRIS_red/multi_400_8500_d560/'
    masterpath = devpath + '/REDUX_OUT/Keck_LRIS_red/multi_400_8500_d560/MF_keck_lris_red/'

    # Read in the msbias for bias subtraction
    biasfile = masterpath + 'MasterBias_A_' + sdet + '_aa.fits'
    msbias = fits.getdata(biasfile)
    # Read in and process flat field images
    pixflat_image_files = np.core.defchararray.add(rawpath, ['r170320_2057.fits','r170320_2058.fits','r170320_2059.fits']).tolist()
    det = 2
    spectro_name = 'keck_lris_red'
    spectrograph = load_spectrograph(spectrograph=spectro_name)
    par = spectrograph.default_pypeit_par()
    flatField = flatfield.FlatField(spectrograph, file_list=pixflat_image_files,det=det, par=par['calibrations']['pixelflatframe']
                                    , msbias = msbias)
    flat = flatField.build_pixflat()
    # Read in the tilts
    tiltsfile = masterpath + 'MasterTilts_A_' + sdet + '_aa.fits'
    mstilts = fits.getdata(tiltsfile)
    # Read in the tslits_dict
    traceslitsroot = masterpath + 'MasterTrace_A_01_aa'
    Tslits = traceslits.TraceSlits.from_master_files(traceslitsroot)
    tslits_dict = {}
    tslits_dict['lcen']=Tslits.lcen
    tslits_dict['rcen']=Tslits.rcen
    tslits_dict['slitpix'] = pixels.slit_pixels(tslits_dict['lcen'],tslits_dict['rcen'], flat.shape, Tslits.par['pad'])
elif type == 'ESI':
    flatfile = '/Users/joe/REDUX/esi_redux/Mar_2008/Flats/FlatECH10_1x1_D.fits.gz'
    piximgfile = '/Users/joe/REDUX/esi_redux/Mar_2008/Final/f_esi1044.fits.gz'
    waveimgfile = '/Users/joe/REDUX/esi_redux/Mar_2008/Arcs/ArcECH10_1x1IMG.fits.gz'
    sedg_file = '/Users/joe/REDUX/esi_redux/Mar_2008/Flats/SEdgECH10_1x1.fits.gz'
    flat = fits.getdata(flatfile)
    (nspec, nspat) = flat.shape
    piximg = fits.getdata(piximgfile, 3)
    mstilts = piximg/nspec
    slit_edges = fits.getdata(sedg_file,0)
    tslits_dict = {}
    tslits_dict['lcen']=slit_edges[0,:,:].T
    tslits_dict['rcen']=slit_edges[0,:,:].T
    tslits_dict['slitpix'] = pixels.slit_pixels(tslits_dict['lcen'],tslits_dict['rcen'], flat.shape, 0)


maskslits = np.zeros(tslits_dict['lcen'].shape[1], dtype=bool)
gdslits = np.where(~maskslits)[0]

# Loop on slits
for slit in gdslits:
    msgs.info("Computing flat field image for slit: {:d}".format(slit + 1))
    slit_left = tslits_dict['lcen'][:, slit]
    slit_righ = tslits_dict['rcen'][:, slit]
    thismask = (tslits_dict['slitpix'] == slit + 1)

    # Function compute_flats(flat, mstilts, slit_left, sligt_righ, thismask,
    # spec_samp_fine = 1.2, spec_samp_coarse = 50, spat_samp =  5.0, spat_illum_thresh = 0.03
    ximg, edgmask = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg=trim_edg)

