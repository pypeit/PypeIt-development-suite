"""
Module to test spectrograph read functions
"""
from pathlib import Path
import os

from IPython import embed

from pypeit import spectrographs
from pypeit.spectrographs.util import load_spectrograph


# test that we can find the correct files
def generic_file_search_test(spec, root, nfiles):
    devroot = os.getenv('PYPEIT_DEV')
    assert devroot is not None, 'PYPEIT_DEV not defined'
    _root = Path(devroot).absolute() / 'RAW_DATA' / spec / root
    assert _root.is_dir(), f'{_root} is not a directory'
    files = load_spectrograph(spec).find_raw_files(_root)
    assert len(files) == nfiles, f'Found {len(files)} {spec} files in {_root}, not {nfiles}'

def test_filesearch():
    #         Spectrograph name, setup directory, expected number of files
    args = [['aat_uhrf', '3875', 5],
            ['bok_bc', '300', 8],
            ['gemini_flamingos2', 'HK_HK', 9],
            ['gemini_gnirs_echelle', '32_SB_SXD', 28],
            ['gemini_gnirs_ifu', 'LR_IFU_32mm', 23],
            ['gtc_osiris', 'R1000B', 9],
            ['keck_deimos', '1200B_LVM_5200', 8],
            ['keck_esi', '1arc_1x1', 10],
            ['keck_hires', 'HS1700+6416_H45aH_RED_B2_ECH_0.00_XD_-0.00_1x2', 13],
            ['keck_kcrm', 'low_rl', 23],
            ['keck_kcwi', 'medium_bl', 9],
            ['keck_lris_blue', 'long_400_3400_d560', 6],
            ['keck_lris_blue_orig', 'long_600_4000_d500', 7],
            ['keck_lris_red', 'long_150_7500_d560', 4],
            ['keck_lris_red_mark4', 'long_400_8500_d560', 12],
            ['keck_lris_red_orig', 'long_1200_7500_d560', 3],
            ['keck_mosfire', 'J2_long', 16],
            ['keck_nires', 'ABBA_wstandard', 10],
            ['keck_nirspec_high', 'Hband', 10],
            ['ldt_deveny', 'DV1', 10],
            ['magellan_fire', 'FIRE', 13],
            ['magellan_mage', '1x1', 51],
            ['mdm_modspec', 'Echelle', 13],
            ['mmt_binospec', 'Longslit_G600', 5],
            ['mmt_bluechannel', '300l', 15],
            ['mmt_mmirs', 'HK_zJ', 5],
            ['not_alfosc', 'grism10', 10],
            ['ntt_efosc2', 'gr4', 9],
            ['p200_dbsp_blue', '1200_5000_d55', 8],
            ['p200_dbsp_red', '1200_7100_d68', 6],
            ['p200_tspec', 'TSPEC', 16],
            ['shane_kast_blue', '452_3306_d57', 11],
            ['shane_kast_red', '1200_5000_d57', 4],
            ['soar_goodman_blue', 'M1', 5],
            ['tng_dolores', 'LRB', 8],
            ['vlt_fors2', '300I', 30],
            ['vlt_sinfoni', 'K_0.8', 14],
            ['wht_isis_blue', 'long_R300B_d5300', 16],
    ]
    for a in args:
        generic_file_search_test(*a)


def test_mdm_osmos():
    s = spectrographs.mdm_osmos.MDMOSMOSMDM4KSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'mdm_osmos',
                                'MDM4K', 'MDM_science.fits')
    assert os.path.isfile(example_file), 'Could not find example file for MDM OSMOS read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (1016, 4128)
    assert bpm.shape == (4060, 1008)

def test_gemini_flamingos():
    s = spectrographs.gemini_flamingos.GeminiFLAMINGOS2Spectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'gemini_flamingos2',
                                'HK_HK', 'S20161108S0072.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Gemini Flamingos2 read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (2048, 2048)
    assert bpm.shape == (2048, 2048)

def test_gemini_gnirs():
    s = spectrographs.gemini_gnirs.GeminiGNIRSSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'gemini_gnirs_echelle',
                                '32_SB_SXD', 'cN20170331S0246.fits')
    assert os.path.isfile(example_file), 'Could not find example file for Gemini GNIRS read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (1022, 1024)
    assert bpm.shape == (1022, 1024)

def test_lbt_luci_ii():
    s = spectrographs.lbt_luci.LBTLUCI2Spectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'lbt_luci',
                                'LUCI-II', 'luci2.20181122.0110.fits')
    assert os.path.isfile(example_file), 'Could not find example file for LBT Luci-II read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (2048,2048)
    assert bpm.shape == (2040, 2040)

def test_lbt_luci_i():
    s = spectrographs.lbt_luci.LBTLUCI1Spectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'lbt_luci',
                                'LUCI-I', 'luci1.20181124.0034.fits')
    assert os.path.isfile(example_file), 'Could not find example file for LBT Luci-I read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (2048,2048)
    assert bpm.shape == (2040, 2040)

def test_gemini_gmos_gmossham():
    s = spectrographs.gemini_gmos.GeminiGMOSSHamSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'gemini_gmos',
                                'GS_HAM_R400_700', 'S20181005S0085.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Gemini GMOS-S Ham read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (512, 1152)
    assert bpm.shape == (1024, 512)


def test_gemini_gmos_gmossham_mosaic():
    s = spectrographs.gemini_gmos.GeminiGMOSSHamSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'gemini_gmos',
                                'GS_HAM_R400_700', 'S20181005S0085.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Gemini GMOS-S Ham read.'
    det = (1,2,3)
    mosaic, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (3, 512, 1152)
    assert bpm.shape == (3, 1024, 512)


def test_magellanfire_echelle():
    s = spectrographs.magellan_fire.MagellanFIREEchelleSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'magellan_fire',
                                'FIRE_Echelle', 'fire_0048.fits')
    assert os.path.isfile(example_file), 'Could not find example file for MagellanFire Echelle read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (2048,2048)
    assert bpm.shape == (2040, 2040)

def test_magellanfire_long():
    s = spectrographs.magellan_fire.MagellanFIRELONGSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'magellan_fire',
                                'FIRE_Long', 'fire_0015.fits')
    assert os.path.isfile(example_file), 'Could not find example file for Magellan Fire Long read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (2048,2048)
    assert bpm.shape == (2040, 351)

def test_keckdeimos():
    s = spectrographs.keck_deimos.KeckDEIMOSSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos',
                                '830G_L_8400', 'd0914_0002.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck DEIMOS read.'
    det = 2
    #data, _ = s.load_raw_frame(example_file, det=det)
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    #
    bpm = s.bpm(example_file, det) #shape=shape) # filename=example_file)
    assert data.shape == (4096,2128)
    assert bpm.shape == (4096,2048)


def test_kecklrisblue():
    s = spectrographs.keck_lris.KeckLRISBSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_blue',
                                'long_400_3400_d560', 'LB.20160109.14149.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck LRIS blue read.'
    det = 2
    #data, _ = s.load_raw_frame(example_file, det=det)
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    #
    bpm = s.bpm(example_file, det) #shape=shape) # filename=example_file)
    assert data.shape == (2048,1154)
    assert bpm.shape == (2048,1024)


def test_kecklrisred():
    s = spectrographs.keck_lris.KeckLRISRSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_red',
                                'long_600_7500_d560', 'LR.20160216.05529.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck LRIS red read.'
    det = 1
    #data, _ = s.load_raw_frame(example_file, det=det)
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    #
    bpm = s.bpm(example_file, det)#, debug=True) #shape=shape) # filename=example_file)
    assert data.shape == (2068,1110)
    assert bpm.shape == (2048,1024)


def test_kecknires():
    s = spectrographs.keck_nires.KeckNIRESSpectrograph()
    # TODO: Any Keck NIRES files to read?


def test_kecknirspec():
    s = spectrographs.keck_nirspec.KeckNIRSPECLowSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_nirspec_low',
                                'LOW_NIRSPEC-1', 'NS.20160414.02637.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck NIRSPEC read.'
    #data, _ = s.load_raw_frame(example_file)
    det=1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'



def test_shanekastredret():
    s = spectrographs.shane_kast.ShaneKastRedRetSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'shane_kast_red',
                                '600_7500_d55_ret', 'r112.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Shane Kast red read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (250, 1232)
    assert bpm.shape == (1200, 250)


def test_shanekastred():
    spectrographs.shane_kast.ShaneKastRedSpectrograph()
    # TODO: Any Shane Kast Red files to read?


def test_tngdolores():
    s = spectrographs.tng_dolores.TNGDoloresSpectrograph()
    # TODO: Any TNG Dolores files to read?


def test_vltxshooteruvb():
    s = spectrographs.vlt_xshooter.VLTXShooterUVBSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'vlt_xshooter',
                                'UVB_1x1', 'XSHOO.2010-04-28T05:34:32.723.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for VLT Xshooter UVB read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (3000, 2144)
    assert bpm.shape == (3000, 2048)


def test_vltxshootervis():
    s = spectrographs.vlt_xshooter.VLTXShooterVISSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'vlt_xshooter',
                                'VIS_1x1', 'XSHOO.2010-04-28T05:34:37.853.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for VLT Xshooter VIS read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (4000, 2106)
    assert bpm.shape == (4000, 2048)


def test_vltxshooternir():
    s = spectrographs.vlt_xshooter.VLTXShooterNIRSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'vlt_xshooter',
                                'NIR', 'XSHOO.2016-08-02T08:45:49.494.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for VLT Xshooter NIR read.'
    det = 1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (1100,2048)
    assert bpm.shape == (2045, 1097)

