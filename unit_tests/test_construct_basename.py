"""
Module to test pypeit.outputfiles.construct_basename using real raw science
frames from a range of spectrographs, chosen to cover the different raw file
extensions supported by PypeIt (see Spectrograph.allowed_extensions).
"""
import os

from IPython import embed

import pytest

from pypeit.spectrographs.util import load_spectrograph
from pypeit import outputfiles


def grab_basename(specstr, rawfile):
    spec = load_spectrograph(specstr)
    headarr = spec.get_headarr(rawfile)
    target = spec.get_meta_value(headarr, 'target')
    mjd = spec.get_meta_value(headarr, 'mjd')
    return outputfiles.construct_basename(rawfile, target, spec.camera, mjd,
                                          spec.allowed_extensions)


def test_construct_basename_aat_uhrf():
    # Raw file extension: '.FTS'
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'aat_uhrf', '3875',
                         'RUN0034.FTS')
    try:
        basename = grab_basename('aat_uhrf', ifile)
    except:
        pytest.fail('AAT UHRF construct_basename failed.')
    assert basename == 'RUN0034-ZETAOPHCN_UHRF_19930511T152129.000'
    assert not basename.endswith('.FTS'), 'Extension not stripped'


def test_construct_basename_wht_isis_blue():
    # Raw file extension: '.fit.gz'
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'wht_isis_blue',
                         'long_R300B_d5300', 'r2324653.fit.gz')
    try:
        basename = grab_basename('wht_isis_blue', ifile)
    except:
        pytest.fail('WHT ISIS Blue construct_basename failed.')
    assert basename == 'r2324653-MMT0155-0821_ISISb_20160209T201715.544'
    assert '.fit' not in basename, 'Extension not stripped'


def test_construct_basename_gemini_gnirs():
    # Raw file extension: '.fits'
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'gemini_gnirs_echelle',
                         '10_LB_SXD', 'N20160717S0114.fits')
    try:
        basename = grab_basename('gemini_gnirs_echelle', ifile)
    except:
        pytest.fail('Gemini GNIRS construct_basename failed.')
    assert basename == 'N20160717S0114-SDSS-J212329.47-005052.9_GNIRS_20160717T092320.135'


def test_construct_basename_gemini_gmos_fitsgz():
    # Raw file extension: '.fits.gz'
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'gemini_gmos',
                         'GS_HAM_R400_700', 'S20181005S0079.fits.gz')
    try:
        basename = grab_basename('gemini_gmos_south_ham', ifile)
    except:
        pytest.fail('Gemini GMOS (.fits.gz) construct_basename failed.')
    assert basename == 'S20181005S0079-FRB180924_GMOS-S_20181005T023228.667'


def test_construct_basename_gemini_gmos_fitsbz2():
    # Raw file extension: '.fits.bz2'
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'gemini_gmos',
                         'GS_HAM_R400_795_SENS', 'S20221207S0183.fits.bz2')
    try:
        basename = grab_basename('gemini_gmos_south_ham', ifile)
    except:
        pytest.fail('Gemini GMOS (.fits.bz2) construct_basename failed.')
    assert basename == 'S20221207S0183-LTT3218_GMOS-S_20221207T073027.350'


def test_construct_basename_soar_goodman():
    # Raw file extension: '.fits.fz'
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'soar_goodman_red', 'M2',
                         '0320_FRB210320_host_05-04-2021.fits.fz')
    try:
        basename = grab_basename('soar_goodman_red', ifile)
    except:
        pytest.fail('SOAR Goodman construct_basename failed.')
    assert basename == '0320_FRB210320_host_05-04-2021-FRB210320_host_red_20210405T050429.443'
    assert '.fits' not in basename, 'Extension not fully stripped'


def test_construct_basename_mdm_modspec():
    # Raw file extension: '.fit'
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'mdm_modspec', 'Echelle',
                         'MDM_Science_Test.fit')
    try:
        basename = grab_basename('mdm_modspec', ifile)
    except:
        pytest.fail('MDM Modspec construct_basename failed.')
    assert basename == 'MDM_Science_Test-KIC9653110_Echelle_20170611T051040.800'
