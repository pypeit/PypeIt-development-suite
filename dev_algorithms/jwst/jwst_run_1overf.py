from astropy.io import fits
import numpy as np
import jwst_1overf_brammer as b
import jwst_1overf_unfold as unfold
from pypeit.display import display

def run_brammer(ratefile, fix_rows=False, in_place=False, writeout=False, savefig=False):

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

def run_unfold(ratefile):
    data = fits.open(ratefile)['SCI'].data
    outimg, modelimg = unfold.fnoise_sub(data)

    return outimg, modelimg

def run_display_b_unfold(ratefile):
    # run and compare the results of brammer and unfold_jwst's corrections
    data = fits.open(ratefile)['SCI'].data
    err = fits.open(ratefile)['ERR'].data

    b_corrdata, b_mod = run_brammer(ratefile)
    unfold_corrdata, unfold_mod = run_unfold(ratefile)
    out = b_corrdata, b_mod, unfold_corrdata, unfold_mod, data

    allim_cuts = [-0.1, 0.1]
    display.connect_to_ginga(raise_err=True, allow_new=True)
    display.show_image(data, chname='data', cuts=allim_cuts)
    display.show_image(b_corrdata, chname='brammer corr', cuts=allim_cuts)
    display.show_image(unfold_corrdata, chname='unfold corr', cuts=allim_cuts)
    display.show_image(b_mod, chname='brammer model', cuts=allim_cuts)
    display.show_image(unfold_mod, chname='unfold model', cuts=allim_cuts)

    # looking at regions far from the slit
    allim_cuts = [-0.1, 0.1]
    display.show_image(data[1250:1980, 570:1450], chname='data (sub)', cuts=allim_cuts)
    display.show_image(b_corrdata[1250:1980, 570:1450], chname='brammer corr (sub)', cuts=allim_cuts)
    display.show_image(unfold_corrdata[1250:1980, 570:1450], chname='unfold corr (sub)', cuts=allim_cuts)
    display.show_image(b_mod[1250:1980, 570:1450], chname='brammer model (sub)', cuts=allim_cuts)
    display.show_image(unfold_mod[1250:1980, 570:1450], chname='unfold model (sub)', cuts=allim_cuts)

    return out