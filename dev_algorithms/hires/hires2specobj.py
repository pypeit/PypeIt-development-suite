import os
from matplotlib import pyplot as plt
from astropy.io import fits
from pypeit.metadata import PypeItMetaData
import numpy as np
from pypeit import specobj
from pypeit import specobjs
from pypeit import utils
from pypeit.spectrographs.util import load_spectrograph
from IPython import embed
import IPython


def read_in_hires_xidl(path, name, outfile, colors=['R', 'G', 'B'], file_type='.fits.gz'):
    """

    :param path: the path to the xidl outputs
    :param name: the name (number) of the exposure
    :param colors: the colors of the three detectors (right now this only works for these three in this order...)
    :param file_type: the ending to the xidl output file to read in
    :return:
    """

    extract_file_0 = os.path.join(path, 'Extract', f'Obj_{name}' + colors[0] + file_type)
    extract_file_1 = os.path.join(path, 'Extract', f'Obj_{name}' + colors[1] + file_type)
    extract_file_2 = os.path.join(path, 'Extract', f'Obj_{name}' + colors[2] + file_type)

    final_file_0 = os.path.join(path, 'Final', f'f_hires{name}' + colors[0] + file_type)
    head_final = fits.getheader(final_file_0)

    # get all the headers
    head2d_0 = fits.getheader(extract_file_0, ext=1)
    head2d_1 = fits.getheader(extract_file_1, ext=1)
    head2d_2 = fits.getheader(extract_file_2, ext=1)

    #TODO I don't think you need this below with the final header
    # find the difference in the headers
    diff_01 = fits.diff.HeaderDiff(head2d_0, head2d_1)
    diff_02 = fits.diff.HeaderDiff(head2d_0, head2d_2)
    diff_12 = fits.diff.HeaderDiff(head2d_2, head2d_1)

    # check to see if any header varies more than the 'NAXIS2' parameter
    if np.any([len(diff_01.diff_keyword_values),
               len(diff_02.diff_keyword_values),
               len(diff_12.diff_keyword_values)] > [1, 1, 1]):
        print('more than just differences in number of orders in the header')


    # define the final header and add the number of rows together for the final header
    head2d = head2d_0.copy()
    head2d['NAXIS2'] = head2d_0['NAXIS2'] + head2d_1['NAXIS2'] + head2d_2['NAXIS2']

    spectrograph = load_spectrograph('keck_hires_red') # apparently spectrographs.keck_hires.py exists, but name is set to keck_hires_red at the moment
    detector = spectrograph.get_detector_par(None, 1)

    n_orders = head2d['NAXIS2']

    cut_0 = head2d_0['NAXIS2']
    cut_1 = head2d_0['NAXIS2'] + head2d_1['NAXIS2']

    sobjs = specobjs.SpecObjs()
    for iord in range(n_orders-1, -1, -1):
        if iord < cut_0:
            hdu = fits.open(extract_file_0)
            order_idx_use = iord
        elif iord < cut_1:
            hdu = fits.open(extract_file_1)
            order_idx_use = iord - cut_0
        else:
            hdu = fits.open(extract_file_2)
            order_idx_use = iord - cut_1

        data = hdu[1].data

        thisobj = specobj.SpecObj('Echelle', DET=1, OBJTYPE='science', ECH_ORDERINDX=iord, ECH_ORDER=data[order_idx_use]['ORDER'])
        # Set boxcar extracted flux
        # TODO same kludge as below for boxcar
        thisobj.BOX_WAVE = data[order_idx_use]['BOX_WV']
        thisobj.BOX_COUNTS = data[order_idx_use]['BOX_FX'].astype(float) # otherwise TypeError
        thisobj.BOX_COUNTS_IVAR = utils.inverse(data[order_idx_use]['BOX_VAR'].astype(float))
        thisobj.BOX_COUNTS_SIG = np.sqrt(data[order_idx_use]['BOX_VAR'].astype(float))
        thisobj.BOX_MASK = data[order_idx_use]['BOX_VAR'] > 0.0
        # Do the same things for the opt flags
        # Find non-unique wavelengths

        # This is a kludge to deal with a bug in xidl that produces multiple values at the same wavelength
        uniq_vals, uniq_ind, uniq_counts = np.unique(data[order_idx_use]['WAVE'], return_index=True, return_counts=True)
        dupl_wave = (uniq_counts > 1) & (uniq_vals != 0.0)
        if np.any(dupl_wave):
            data[order_idx_use]['WAVE'][uniq_ind[dupl_wave] + 1] = 0.0
            data[order_idx_use]['VAR'][uniq_ind[dupl_wave] + 1]  = 0.0

        thisobj.OPT_WAVE = data[order_idx_use]['WAVE']
        thisobj.OPT_COUNTS = data[order_idx_use]['FX'].astype(float)
        thisobj.OPT_COUNTS_IVAR = utils.inverse(data[order_idx_use]['VAR'].astype(float))
        thisobj.OPT_COUNTS_SIG = np.sqrt(data[order_idx_use]['VAR'].astype(float))
        thisobj.OPT_MASK = data[order_idx_use]['VAR'] > 0.0

        thisobj.TRACE_SPAT = data[order_idx_use]['TRACE'].astype(float)
        #thisobj.FHWM = data[order_idx_use]['SPATIAL_FWHM'] # want this?
        thisobj.NAME = data[order_idx_use]['FIELD'] # want this?
        thisobj.OPT_COUNTS_SKY = data[order_idx_use]['SKY'].astype(float) # want this?
        thisobj.ECH_FRACPOS = 0.5 # This a bogus value to ensure that we get a decent name
        thisobj.DETECTOR = detector
        thisobj.SPAT_PIXPOS = np.median(data[order_idx_use]['TRACE'][data[order_idx_use]['VAR']> 0.0])
        thisobj.WAVE_RMS = 0.0
        thisobj.set_name()
        sobjs.add_sobj(thisobj)


    # Hack in the detector container nonsense

    #par = spectrograph.default_pypeit_par()
    # Trying to hack fitstbl
    #fitstbl = PypeItMetaData(spectrograph, par)

    # Somehow hack the fitstbl or whatever
    # keywords from pypeit.core.meta.define_core_meta
    row_fitstbl = hdu[0].header  # this is all empty for each so don't need to combine all 3
    row_fitstbl['filename'] = None
    row_fitstbl['ra'] = head_final['RA']
    row_fitstbl['dec'] = head_final['DEC']
    row_fitstbl['airmass'] = float(head_final['AIRMASS'])
    row_fitstbl['exptime'] = float(head_final['EXPTIME'])
    row_fitstbl['target'] = head_final['TARGNAME']
    row_fitstbl['mjd'] = head_final['MJD']
    row_fitstbl['decker'] = head_final['DECKNAME']
    row_fitstbl['dispname'] = 'HIRES-R'


    # row_fitstbl['ra'] = None
    # row_fitstbl['dec'] = None
    # row_fitstbl['target'] = None
    # row_fitstbl['dispname'] = None
    # row_fitstbl['decker'] = None
    # row_fitstbl['binning'] = None
    # row_fitstbl['mjd'] = None
    # row_fitstbl['airmass'] = None
    # row_fitstbl['exptime'] = None

    subheader = spectrograph.subheader_for_spec(row_fitstbl, head2d, allow_missing=True)

    outfile_txt = outfile.replace('.fits', '.txt')
    sobjs.write_to_fits(subheader, outfile)
    sobjs.write_info(outfile_txt, spectrograph.pypeline)


if __name__ == '__main__':
    # Task for you is to read in the standard star reduction from Extract directory and make some plots of the
    # of the counts vs wavelength for each order. The relevant tags are wave, fx, var (sqrt var is sigma)

    # Extra credit (you probably won't get this)

    # Try to create a SpecObjs in echelle format that contains the data from that standard file. Look at echelle in Dev Suite forPypeit for examples.

    # Read in a Cowie xidl output for a single exposure
    # xidl_file = '/Users/joe/HIRES_redux/J0100+2802/Cowie/Extract/Obj_0147R.fits.gz'
    # xidl_file = '/Users/suksientie/Research/J0100+2802/Cowie/Extract/Obj_0147R.fits.gz'
    # xidl_file = '/home/molly/projects/OBS/HIRES/J0100+2802/Cowie/Extract/Obj_0147R.fits.gz'

    outpath = './'
    xidl_path = '/Users/joe/HIRES_redux/J0100+2802/Cowie/'
    obj_name = '0143'
    outfile = os.path.join(outpath, f'spec1d_{obj_name}.fits')

    read_in_hires_xidl(xidl_path, obj_name, outfile)
