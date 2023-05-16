import os
from astropy.io import fits
import numpy as np
import jwst_1overf_brammer as b
import jwst_1overf_unfold as unfold
from pypeit.display import display
import jwst_mosaic_slits as jms

def run_brammer(ratefile, bpm=None, fix_rows=False, in_place=False, writeout=False, make_plot=False):

    im = fits.open(ratefile)
    im_data = im['SCI'].data.T

    fig, mod, data = b.exposure_oneoverf_correction(ratefile, rot_pypeit_fmt=True, bpm=bpm, in_place=in_place, deg_pix=256, make_plot=make_plot) # default
    mod[~np.isfinite(mod)] = 0

    if fix_rows:
        axis = 1
        fig, fix_row_mod, data = b.exposure_oneoverf_correction(ratefile, rot_pypeit_fmt=True, bpm=bpm, axis=axis, in_place=in_place, deg_pix=2048, make_plot=make_plot)

        fix_row_mod[~np.isfinite(fix_row_mod)] = 0
        total_model = mod + fix_row_mod

    else:
        total_model = mod

    corr_data = im_data - total_model

    if writeout:
        rate1overffile = ratefile.replace('_rate.fits', '_rate_1overf.fits')
        im_data -= total_model
        im.writeto(rate1overffile, overwrite=True)

    if make_plot:
        display.connect_to_ginga(raise_err=True, allow_new=True)
        display.show_image(im_data, chname='raw data')
        display.show_image(corr_data, chname='outimg')
        display.show_image(total_model, chname='total model')

    return corr_data, total_model

def old_run_brammer(ratefile, fix_rows=False, in_place=False, writeout=False, savefig=False):

    im = fits.open(ratefile)
    axis=0 # default
    fig, mod, data = b.exposure_oneoverf_correction(ratefile, erode_mask=False, in_place=in_place, axis=axis, deg_pix=256) # default
    if savefig:
        fig.savefig(ratefile.split('.fits')[0] + f'_onef_axis{axis}.png')

    mod[~np.isfinite(mod)] = 0
    im['SCI'].data -= mod

    if fix_rows:
        axis=1
        fig, mod, data = b.exposure_oneoverf_correction(ratefile, erode_mask=False, in_place=in_place, axis=axis, deg_pix=2048)

        fig.savefig(ratefile.split('.fits')[0] + f'_onef_axis{axis}.png')
        mod[~np.isfinite(mod)] = 0
        im['SCI'].data -= mod

    if writeout:
        rate1overffile = ratefile.replace('_rate.fits', '_rate_1overf.fits')
        im.writeto(rate1overffile, overwrite=True)

    return im['SCI'].data, mod

def get_objmask(ratefile, calfile, flatfile, slitname=None, show=False):
    nrs_calib, gdslits = jms.get_one_calib_obj(calfile, flatfile)
    data_t = fits.open(ratefile)['SCI'].data.T
    objmask = np.zeros(np.shape(data_t)).astype(bool)

    if slitname is not None:
        i = np.argwhere(np.array(gdslits) == slitname).squeeze()
        objmask[nrs_calib[i].slit_slice] = True

    else:
        for i in range(len(nrs_calib)):
            objmask[nrs_calib[i].slit_slice] = True

    if show:
        display.connect_to_ginga(raise_err=True, allow_new=True)
        display.show_image(data_t, chname='data')
        for i in range(len(nrs_calib)):
            display.show_image(data_t[nrs_calib[i].slit_slice], chname='%s' % gdslits[i])

    return objmask

def run_unfold(ratefile, bpm=None, namp=4, evenOdd=True, skip_col=False, show=False):
    data = fits.open(ratefile)['SCI'].data.T
    outimg, modelimg = unfold.fnoise_sub(data, bpm=bpm, namp=namp, sub_bkg=False, mask_brightstar=False, evenOdd=evenOdd, skip_col=skip_col, show=show)

    if show:
        display.connect_to_ginga(raise_err=True, allow_new=True)
        display.show_image(data, chname='raw data')
        display.show_image(outimg, chname='outimg')
        display.show_image(modelimg, chname='modelimg')

    return outimg, modelimg

def run_display_b_unfold(ratefile, calfile, flatfile, show=False):

    objmask_all = get_objmask(ratefile, calfile, flatfile, show=False)

    # run and compare the results of brammer and unfold_jwst's corrections
    data = fits.open(ratefile)['SCI'].data.T
    err = fits.open(ratefile)['ERR'].data

    b_corrdata, b_mod = run_brammer(ratefile, bpm=objmask_all, fix_rows=True, make_plot=show)
    unfold_corrdata, unfold_mod = run_unfold(ratefile, bpm=objmask_all, namp=4, evenOdd=True, skip_col=False, show=show)
    out = b_corrdata, b_mod, unfold_corrdata, unfold_mod, data

    return out

    allim_cuts = [-0.1, 0.1]
    display.connect_to_ginga(raise_err=True, allow_new=True)
    display.show_image(data, chname='data', cuts=allim_cuts)
    display.show_image(b_corrdata, chname='brammer corr', cuts=allim_cuts)
    display.show_image(unfold_corrdata, chname='unfold corr', cuts=allim_cuts)
    display.show_image(b_mod, chname='brammer model', cuts=allim_cuts)
    display.show_image(unfold_mod, chname='unfold model', cuts=allim_cuts)

    """
    # looking at regions far from the slit
    allim_cuts = [-0.1, 0.1]
    display.show_image(data[1250:1980, 570:1450], chname='data (sub)', cuts=allim_cuts)
    display.show_image(b_corrdata[1250:1980, 570:1450], chname='brammer corr (sub)', cuts=allim_cuts)
    display.show_image(unfold_corrdata[1250:1980, 570:1450], chname='unfold corr (sub)', cuts=allim_cuts)
    display.show_image(b_mod[1250:1980, 570:1450], chname='brammer model (sub)', cuts=allim_cuts)
    display.show_image(unfold_mod[1250:1980, 570:1450], chname='unfold model (sub)', cuts=allim_cuts)
    """
    return out

if __name__ == "__main__":


    #disperser = 'PRISM_02073'
    disperser = '235H'

    detectors = ['nrs1', 'nrs2']
    exp_list = []

    reduce_slits = ['S200A1']

    target = 'J0313'
    mode = 'FS'
    #mode = 'MSA'
    detectors = ['nrs1', 'nrs2']
    exp_list = []
    for detname in detectors:
        if disperser == 'PRISM_02073':
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/NIRSPEC_2073/level_12/02073/'
            redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/02073_CLEAR_PRISM/calwebb'
            output_dir = os.path.join(redux_dir, 'output')
            pypeit_output_dir = os.path.join(redux_dir, 'pypeit')

            # NIRSPEC 3-point dither
            # dither center
            scifile1 = os.path.join(rawpath_level2, 'jw02073008001_03101_00001_' + detname + '_rate.fits')
            scifile2 = os.path.join(rawpath_level2, 'jw02073008001_03101_00002_' + detname + '_rate.fits')
            scifile3 = os.path.join(rawpath_level2, 'jw02073008001_03101_00003_' + detname + '_rate.fits')
        elif 'J0313' == target:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1764/level_12/01764/'
            output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J0313/calwebb/output'
            pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J0313/calwebb/pypeit'
            if disperser == '235H':
                # NIRSPEC 3-point dither dither center
                if reduce_slits[0] == 'S200A1':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764014001_03102_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764014001_03102_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764014001_03102_00003_' + detname + '_rate.fits')
                elif reduce_slits[0] == 'S200A2':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764014001_03104_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764014001_03104_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764014001_03104_00003_' + detname + '_rate.fits')
            elif disperser == '395H':
                # NIRSPEC 3-point dither dither center
                if reduce_slits[0] == 'S200A2':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764014001_03106_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764014001_03106_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764014001_03106_00003_' + detname + '_rate.fits')
                elif reduce_slits[0] == 'S200A1':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764014001_03108_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764014001_03108_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764014001_03108_00003_' + detname + '_rate.fits')
            elif disperser == '140H':
                # NIRSPEC 3-point dither dither center
                if reduce_slits[0] == 'S200A1':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764014001_0310a_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764014001_0310a_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764014001_0310a_00003_' + detname + '_rate.fits')
                elif reduce_slits[0] == 'S200A2':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764014001_0310c_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764014001_0310c_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764014001_0310c_00003_' + detname + '_rate.fits')
        elif 'J1007' == target:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1764/level_12/01764/'
            output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J1007/calwebb/output'
            pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/J1007/calwebb/pypeit'
            if disperser == '235H':
                # NIRSPEC 3-point dither dither center
                if reduce_slits[0] == 'S200A1':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764006001_04102_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764006001_04102_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764006001_04102_00003_' + detname + '_rate.fits')
                elif reduce_slits[0] == 'S200A2':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764006001_04104_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764006001_04104_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764006001_04104_00003_' + detname + '_rate.fits')
            elif disperser == '395H':
                # NIRSPEC 3-point dither dither center
                if reduce_slits[0] == 'S200A2':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764006001_04106_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764006001_04106_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764006001_04106_00003_' + detname + '_rate.fits')
                elif reduce_slits[0] == 'S200A1':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764006001_04108_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764006001_04108_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764006001_04108_00003_' + detname + '_rate.fits')
            elif disperser == '140H':
                # NIRSPEC 3-point dither dither center
                if reduce_slits[0] == 'S200A1':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764006001_0410a_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764006001_0410a_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764006001_0410a_00003_' + detname + '_rate.fits')
                elif reduce_slits[0] == 'S200A2':
                    scifile1 = os.path.join(rawpath_level2, 'jw01764006001_0410c_00001_' + detname + '_rate.fits')
                    scifile2 = os.path.join(rawpath_level2, 'jw01764006001_0410c_00002_' + detname + '_rate.fits')
                    scifile3 = os.path.join(rawpath_level2, 'jw01764006001_0410c_00003_' + detname + '_rate.fits')

        exp_list.append([scifile1, scifile2, scifile3])

    scifiles_1 = exp_list[0]
    scifiles_2 = exp_list[1] if len(exp_list) > 1 else []
    scifiles = [scifiles_1, scifiles_2]
    scifiles_all = scifiles_1 + scifiles_2
    nexp = len(scifiles_1)

    basenames = []
    basenames_1 = []
    basenames_2 = []
    for sci1, sci2 in zip(scifiles_1, scifiles_2):
        b1 = os.path.basename(sci1).replace('_rate.fits', '')
        b2 = os.path.basename(sci2).replace('_rate.fits', '')
        basenames_1.append(b1)
        basenames_2.append(b2)
        basenames.append(b1.replace('_nrs1', ''))


    # Output file names
    intflat_output_files_1 = []
    msa_output_files_1 = []
    cal_output_files_1 = []

    intflat_output_files_2 = []
    msa_output_files_2 = []
    cal_output_files_2 = []

    for base1, base2 in zip(basenames_1, basenames_2):
        if mode == 'MSA':
            msa_output_files_1.append(os.path.join(output_dir, base1 + '_msa_flagging.fits'))
            msa_output_files_2.append(os.path.join(output_dir, base2 + '_msa_flagging.fits'))
        else:
            msa_output_files_1.append(os.path.join(output_dir, base1 + '_assign_wcs.fits'))
            msa_output_files_2.append(os.path.join(output_dir, base2 + '_assign_wcs.fits'))

        intflat_output_files_1.append(os.path.join(output_dir, base1 + '_interpolatedflat.fits'))
        cal_output_files_1.append(os.path.join(output_dir, base1 + '_cal.fits'))
        intflat_output_files_2.append(os.path.join(output_dir, base2 + '_interpolatedflat.fits'))
        cal_output_files_2.append(os.path.join(output_dir, base2 + '_cal.fits'))



    # Testing 1/f
    bpm_obj = get_objmask(scifiles_2[0], cal_output_files_2[0], intflat_output_files_2[0], slitname=None, show=True)
    run_brammer(scifiles_2[0], bpm=bpm_obj, fix_rows=False, in_place=False, writeout=False, make_plot=True)
    run_unfold(scifiles_2[0], bpm=bpm_obj, namp=4, evenOdd=True, skip_col=False, show=True)

    # Read in multi exposure calwebb outputs
    msa_multi_list_1 = []
    intflat_multi_list_1 = []
    final_multi_list_1 = []
    msa_multi_list_2 = []
    intflat_multi_list_2 = []
    final_multi_list_2 = []
    #nslits_1 = np.zeros(nexp, dtype=int)
    #nslits_2 = np.zeros(nexp, dtype=int)
    #t_eff = np.zeros(nexp, dtype=float)

    #ndetectors = 2
    # Create arrays to hold JWST spec2, but only load the files when they're needed
    #msa_data = np.empty((ndetectors, nexp), dtype=object)
    #flat_data = np.empty((ndetectors, nexp), dtype=object)
    #cal_data = np.empty((ndetectors, nexp), dtype=object)

    #dither_offsets = np.zeros((ndetectors,nexp), dtype=float)
    # TODO: This probably isn't correct.  I.e., need to know offsets and slit
    # position angle.
    #for iexp in range(nexp):
    #    with fits.open(scifiles_1[iexp]) as hdu:
    #        dither_offsets[0,iexp] = hdu[0].header['YOFFSET']
    #for idet in range(1,ndetectors):
    #    dither_offsets[idet] = dither_offsets[0]
    #dither_offsets_pixels = dither_offsets.copy()
    #for idet in range(ndetectors):
    #    dither_offsets_pixels[idet] /= det_container_list[idet].platescale
    # NOTE: Sign convention requires this calculation of the offset
    #dither_offsets_pixels = dither_offsets_pixels[:,0,None] - dither_offsets_pixels
    #print(dither_offsets_pixels)



