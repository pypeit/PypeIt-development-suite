
import numpy as np
from ginga.util import grc
import subprocess
from astropy.io import fits
from astropy.stats import sigma_clipped_stats


def clear_all():
    """
    Clear all of the ginga canvasses

    """
    viewer = connect_to_ginga()
    shell = viewer.shell()
    chnames = shell.get_channel_names()
    for ch in chnames:
        shell.delete_channel(ch)

def connect_to_ginga(host='localhost', port=9000, raise_err=False, allow_new=False):
    """
    Connect to a RC Ginga.

    Args:
        host (:obj:`str`, optional):
            Host name.
        port (:obj:`int`, optional):
            Probably should remain at 9000
        raise_err (:obj:`bool`, optional):
            Raise an error if no connection is made, otherwise just
            raise a warning and continue
        allow_new (:obj:`bool`, optional):
            Allow a subprocess to be called to execute a new ginga
            viewer if one is not already running.

    Returns:
        RemoteClient: connection to ginga viewer.
    """
    # Start
    viewer = grc.RemoteClient(host, port)
    # Test
    sh = viewer.shell()
    try:
        tmp = sh.get_current_workspace()
    except:
        if allow_new:
            subprocess.Popen(['ginga', '--modules=RC'])

            # NOTE: time.sleep(3) is now insufficient. The loop below
            # continues to try to connect with the ginga viewer that
            # was just instantiated for a maximum number of iterations.
            # If the connection is remains unsuccessful, an error is
            # thrown stating that the connection timed out.
            maxiter = int(1e6)
            for i in range(maxiter):
                try:
                    viewer = grc.RemoteClient(host, port)
                    sh = viewer.shell()
                    tmp = sh.get_current_workspace()
                except:
                    continue
                else:
                    break
            if i == maxiter-1:
                print('Timeout waiting for ginga to start.  If window does not appear, type '
                '`ginga --modules=RC` on the command line.  In either case, wait for '
                'the ginga viewer to open and try the pypeit command again.')
            return viewer

        if raise_err:
            raise ValueError
        else:
            print('Problem connecting to Ginga.  Launch an RC Ginga viewer and '
                      'then continue: \n    ginga --modules=RC')

    # Return
    return viewer


def show_image(inp, chname='Image', waveimg=None, mask=None, exten=0, cuts=None,
               clear=False, wcs_match=False):
    """
    Display an image using Ginga.

    .. todo::
        - implement instrument specific reading
        - use the `mask` as a boolean mask if `bitmask` is not provided.

    Args:
        inp (:obj:`str`, numpy.ndarray):
            The image to view.  If a string is provided, it must be the
            name of a fits image that can be read by `astropy.io.fits`.
        chname (:obj:`str`, optional):
            The name of the ginga channel to use.
        waveimg (:obj:`str`, optional):
            The name of a FITS image with the relevant WCS coordinates
            in its header, mainly for wavelength array.  If None, no WCS
            is used.
        exten (:obj:`int`, optional):
            The extension of the fits file with the image to show.  This
            is only used if the input is a file name.
        cuts (array-like, optional):
            Initial cut levels to apply when displaying the image.  This
            object must have a length of 2 with the lower and upper
            levels, respectively.
        clear (:obj:`bool`, optional):
            Clear any existing ginga viewer and its channels.
        wcs_match(:obj:`bool`, optional):
            Use this as a reference image for the WCS and match all
            image in other channels to it.

    Returns:
        ginga.util.grc.RemoteClient, ginga.util.grc._channel_proxy: The
        ginga remote client and the channel with the displayed image.

    Raises:
        ValueError:
            Raised if `cuts` is provided and does not have two elements
            or if bitmask is provided but mask is not.
    """
    # Input checks
    if cuts is not None and len(cuts) != 2:
        raise ValueError('Input cuts must only have two elements, the lower and upper cut.')

    # Instantiate viewer
    viewer = connect_to_ginga()
    # Read or set the image data.  This will fail if the input is a
    # string and astropy.io.fits cannot read the image.
    img = io.fits_open(inp)[exten].data if isinstance(inp, str) else inp

    if clear:
        clear_all()

    ch = viewer.channel(chname)
    # Header
    header = {}
    header['NAXIS1'] = img.shape[1]
    header['NAXIS2'] = img.shape[0]

    # Giddy up
#    waveimg = None
    if waveimg is not None:
        sh = viewer.shell()
        args = [chname, chname, grc.Blob(img.tobytes()), img.shape, img.dtype.name, header,
                grc.Blob(waveimg.tobytes()), waveimg.dtype.name, {}]
        sh.call_global_plugin_method('SlitWavelength', 'load_buffer', args, {})
    else:
        ch.load_np(chname, img, 'fits', header)

    # These commands set up the viewer. They can be found at
    # ginga/ginga/ImageView.py
    canvas = viewer.canvas(ch._chname)
    out = canvas.clear()
    out = ch.set_color_map('ramp')
    out = ch.set_intensity_map('ramp')
    out = ch.set_color_algorithm('linear')
    out = ch.restore_contrast()
    out = ch.restore_cmap()
    if cuts is not None:
        out = ch.cut_levels(cuts[0], cuts[1])

    # WCS Match this to other images with this as the reference image?
    if wcs_match:
        # After displaying all the images since up the images with WCS_MATCH
        shell = viewer.shell()
        out = shell.start_global_plugin('WCSMatch')
        out = shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [chname], {})


    # TODO: I would prefer to change the color map to indicate these
    # pixels rather than overplot points. Because for large numbers of
    # masked pixels, this is super slow. Need to ask ginga folks how to
    # do that.

    return viewer, ch

connect_to_ginga(raise_err=True, allow_new=True)
# Generate a fake image
img1 = fits.getdata('/Users/joe/ginga_test.fits')
img_mask = (fits.getdata('/Users/joe/ginga_mask.fits')).astype(bool)
mean, med, sigma = sigma_clipped_stats(img1[img_mask], sigma_lower=5.0, sigma_upper=5.0)
cuts = (med - 4.0 * sigma, med + 4.0 * sigma)

waveimg = np.repeat(np.arange(img1.shape[0])[:,np.newaxis], img1.shape[1], axis=1)
#viewer, ch1 = show_image(img1, waveimg=waveimg, cuts = (-0.5, 0.5), chname='IMG1', wcs_match = True, clear=True)
viewer, ch1 = show_image(img1, chname='IMG_WITHOUT_WAVE_OR_CUTS', clear=True)
viewer, ch1 = show_image(img1, cuts= (-0.3313194903428431, 0.3306333975688511), chname='IMG_WITHOUT_WAVE_CUTS')
viewer, ch1 = show_image(img1, cuts=cuts, chname='IMG_WITHOUT_WAVE_CUTS2')
viewer, ch1 = show_image(img1, waveimg=waveimg, chname='IMG_WITH_WAVE_NOCUTS')
viewer, ch1 = show_image(img1, waveimg=waveimg, cuts= (-0.3313194903428431, 0.3306333975688511), chname='IMG_WITH_WAVE_CUTS1')
viewer, ch1 = show_image(img1, waveimg=waveimg, cuts= cuts, chname='IMG_WITH_WAVE_CUTS2')
