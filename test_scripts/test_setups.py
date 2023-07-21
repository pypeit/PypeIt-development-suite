# See top-level LICENSE.rst file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This modules defines all contents of the PypeIt dev suite.  It is intended to be a single place to edit when adding new
instruments, setups, tests or test types.

To add a new instrument and/or setup to the dev suite:

1) Make sure the necessary raw data, pypeit files, and other test files are in the PypeIt-development-suite repo and
   the development suite google drive. See https://pypeit.readthedocs.io/en/latest/development.html#development-suite
2) Add the instrument and setup(s) to the all_setups dict
3) If this is a new instrument, add the instrument to the supported_instruments list.
3) If additional tests are desired add 'instrument/setup' to the desired test attribute.

To add a new test for an existing setup:

1) Find the appropriate XXX_tests attribute for the test type. For example, to add a quick look test, edit
   quick_look_tests.
2) If the test type does not require any arguments, the attribute will be a set. Simply add the setup.
3) If the test type requires arguments, they will be passed to the PypeItTest subclasses constructor as keyword
   arguments. Put the arguments into a dict using the keyword argument as a key. For example:

   telluric_tests = {'gemini_gnirs_echelle/32_SB_SXD':
                         {'coadd_file': 'pisco_coadd.fits', 'redshift': 7.52, 'objmodel': 'qso'},

To add a new type of test:

1) Edit pypeit_tests.py to add a new subclass to run the new test type as a child process.
2) Add a new test list with at least one instrument/setup that runs the new test.
3) Add the test to the all_tests list. This list defines the order the test types are run for a test setup, the
   PypeItTest subclass that runs the test, and the test phase (prep, reduce, afterburn, quicklook).

Attributes:
    _reduce_setups:          The test setups that support reduction. A dict of instruments to the supported test
                             setups for the instrument.
                             Each setup should have data in $PYPEIT_DEV/RAW_DATA/instrument/setup
    all_setups:              All of test setups that comprise the "reduce", "afterburn", and "ql" tests in the dev
                             suite. Effectively all of the test that are not run by pytest.
    all_tests:               A list of the test types supported by the dev suite and which test setups they are run
                             on.  The test types are listed in the order they run in so that tests can depend on the
                             result of previous tests.

                             Each test type is represented by a dict with the following keys:

                             'factory': A callable that builds a PypeItTest subclasss to run the test.
                             This callable will be passed a TestSetup object, the command line arguments to
                             pypeit_test (as returned by argparse.ArgumentParser), and any keyword arguments included
                             in the 'setups' key.

                             'type': A TestPhase enum that is either PREP, REDUCE, AFTERBURN, or QL.

                             'setups': Which setups should run the test along with any arguments needed to run the test.

                             The setup can also be specified as an instrument name to indicate every setup for the
                             instrument should run the test type, or as 'instrument/setup' to indicate only a specific
                             setup should run the test type.

                             The value for 'setups' can be a list of setups, or a dict. If a dict is used, it maps
                             a setup specifier to the keyword arguments that will be passed to the PypeItTest __init__
                             method.  The key word arguments are represented as a dict and passed to the __init__
                             method with the ** operator.

                            The 'setups' values are taken from the following private attributes of this module:

                            _pypeit_setup:      Test setups that run pypeit_setup to generate a .pypeit file.
                            _additional_reduce: Test setups that run additional reduce tests beyond the default test.
                            _sensfunc:          Test setups that run pypeit_sensfunc.
                            _flux_setup:        Test setups that run pypeit_flux_setup.
                            _flux:              Test setups that run pypeit_flux_calib.
                            _flexure:           Test setups that run pypeit_deimos_flexure.
                            _coadd1d:           Test setups that run pypeit_coadd_1dspec.
                            _coadd2d:           Test setups that run pypeit_coadd_2dspec.
                            _telluric:          Test setups that run pypeit_tellfit.
                            _quick_look:        Test setups that run quick look script. The actual script run is chosen
                                                based on the instrument.

"""

from . import pypeit_tests
from enum import Enum, IntEnum, auto
import copy

class TestPhase(Enum):
    """Enumeration for specifying the test phase that a test runs in.

    Values:

    PREP
    REDUCE
    AFTERBURN
    QL
    UNIT
    """
    PREP      = auto()
    REDUCE    = auto()
    AFTERBURN = auto()
    QL        = auto()



# This dict specifies all the instruments and setups that are supported by the dev suite.
#  The keys are the instruments and the values are a list of the supported setups.
all_setups  = {
    'bok_bc': ['300','600'],
    'gemini_gnirs_echelle': ['32_SB_SXD', '10_LB_SXD'],
    'gemini_gnirs_ifu': ['LR_IFU_32mm'],
    'gemini_gmos': ['GS_HAM_R400_700', 'GS_HAM_R400_860',
                    'GN_HAM_R400_885', 'GN_HAM_NS_B600_620',
                    'GS_HAM_MULTI_R400_700', 'GN_E2V_MULTI_R400_600',
                    'GS_HAM_B600_MOS'],
    'gemini_flamingos2': ['HK_HK', 'JH_JH'],
    'gtc_osiris': ['R1000B', 'R1000BMOS', 'R1000RMOS', 'R2500R', 'R2500V'],
    'gtc_osiris_plus': ['R1000R', 'R300B'],
    'keck_esi': ['Ech_1x1', 'Ech_2x1'],
    'keck_deimos': ['600ZD_M_6500', '600ZD_tilted', '1200G_M_7750', '830G_LVM_8400', '830G_M_8100_26',
                    '830G_M_8500', '830G_L_8100', '1200B_M_5200', '1200G_M_5500',
                    '900ZD_M_6000', '1200B_LVM_5200', '900ZD_LVM_5500',
                    '830G_M_9000_dither'],
#    'keck_hires': ['J0100+2802_RED_C1_ECH_-0.82_XD_1.62_1x2',
#                   'J0100+2802_RED_C1_ECH_-0.91_XD_1.46_1x2',
#                   #'J0100+2802_RED_C1_ECH_0.88_XD_1.46_1x2',
#                   'J0100+2802_RED_C2_ECH_0.74_XD_1.39_1x3',
#                   'HS1700+6416_H45aH_RED_B2_ECH_0.00_XD_-0.00_1x2',
#                   'J0306+1853_U074_RED_C2_ECH_0.72_XD_1.42_1x3',
#                   'J0306+1853_U074_RED_C2_ECH_-0.86_XD_1.31_1x3',
#                   'J1723+2243_W241_RED_C5_ECH_0.08_XD_0.90_2x2',
#                   'J1723+2243_W241_RED_C5_ECH_-0.15_XD_0.90_2x2',
#                   ],
    'keck_kcwi': ['bh2_4200', 'bl'],
    'keck_nires': ['ABBA_wstandard', 'ABBA_nostandard', 'ABC_nostandard', 'ABpat_wstandard', 'ABBA_nostandard_faint'],
    'keck_nirspec': ['LOW_NIRSPEC-1'],
    'keck_mosfire': ['Y_long', 'J_multi', 'K_long', 'Y_multi', 'long2pos1_H', 'longslit_3x0.7_H', 'mask1_K_with_continuum', 'mask1_J_with_continuum', 'J2_long'],
    'keck_lris_blue': ['multi_600_4000_d560', 'long_400_3400_d560', 'long_600_4000_d560',
                      'multi_300_5000_d680', 'multi_600_4000_slitmask'],
    'keck_lris_blue_orig': ['long_600_4000_d500'],
    'keck_lris_red': ['long_600_7500_d560', 'multi_1200_9000_d680_1x2',
                      'multi_600_5000_d560', 'multi_1200_9000_d680',
                      'multi_400_8500_d560', 'long_600_10000_d680',
                      'long_400_8500_longread'],  # Longslit read-out mode
    'keck_lris_red_orig': ['long_300_5000'],
    'keck_lris_red_mark4': ['long_400_8500_d560', 'long_600_10000_d680', 'multi_600_10000_slitmask'],
    'lbt_luci': ['LUCI-I', 'LUCI-II'],
    'lbt_mods': ['MODS1R_Longslit', 'MODS2R_Longslit'],
    'ldt_deveny': ['DV1', 'DV2', 'DV3', 'DV4', 'DV5', 'DV6', 'DV7', 'DV8', 'DV9'],
    'magellan_fire': ['FIRE_Echelle', 'FIRE_Long'],
    'magellan_mage': ['1x1'],
    'mdm_osmos': ['MDM4K'],
    'mdm_modspec': ['Echelle'],
    'mmt_binospec': ['Longslit_G600', 'Multislit_G270', 'Longslit_G1000'],
    'mmt_mmirs': ['HK_zJ', 'J_zJ', 'K_K'],
    'mmt_bluechannel': ['300l', '500GPM', '800GPM', '832GPM_1st', '832GPM_2nd', '1200GPM'],
    'ntt_efosc2': ['gr5', 'gr6'],
    'not_alfosc': ['grism3', 'grism4', 'grism5', 'grism7', 'grism10', 'grism11', 'grism17', 'grism18', 'grism19', 'grism20', 'grism4_nobin'],
    'p200_dbsp_blue': ['600_4000_d55', '600_4000_d68', '1200_5000_d68'],
    'p200_dbsp_red': ['316_7500_d55', '600_10000_d55', '1200_7100_d68'],
    'p200_tspec': ['TSPEC'],
    'shane_kast_blue': ['452_3306_d57', '600_4310_d55', '830_3460_d46'],
    'shane_kast_red': ['300_7500_Ne', '600_7500_d55_ret', '600_7500_d57', '600_5000_d46', '1200_5000_d57'],
    'soar_goodman_red': ['M1', 'M2', '600red'],
    'soar_goodman_blue': ['M1'],
    'tng_dolores': ['LRB'],
    'vlt_fors2': ['300I', '600Z'],
    'vlt_sinfoni': ['K_0.8'],
    'vlt_xshooter': ['VIS_1x1', 'VIS_2x1', 'VIS_2x2', 'VIS_manual', 'NIR', 'UVB_1x1', 'UVB_1x1_Feige110', 'VIS_1x1_Feige110', 'NIR_Feige110'],
}

# Tests for full reductions
_reduce_setups = {}
for instr in all_setups:
    _reduce_setups[instr] = {}
    for setup in all_setups[instr]:
        _reduce_setups[instr][setup] = [{}]

_pypeit_setup = {
    'shane_kast_blue': {
        '600_4310_d55': [{}]}}

_additional_reduce = {
    'keck_lris_red': {
        'long_600_7500_d560': [dict(ignore_calibs=True)]},
    'gemini_gmos': {
        'GS_HAM_R400_860': [dict(std=True)]},
    }

_sensfunc = {
    'shane_kast_blue': {
        '600_4310_d55': [dict(std_file='spec1d_*Feige66*.fits')]},
    'shane_kast_red': {
        '600_7500_d55_ret': [dict(std_file='spec1d_*G191b2b*.fits',
                                  sens_file="shane_kast_red_600_7500_d55_ret.sens")]},
    'gemini_gnirs_echelle': {
        '32_SB_SXD': [dict(std_file='spec1d_*S0206-HIP62745*.fits',
                           sens_file='gemini_gnirs_echelle_32_sb_sxd.sens')]},
    'gemini_gmos': {
        'GS_HAM_R400_860': [dict(std_file='spec1d_**GD71*.fits',
                                 sens_file='gemini_gmos_gs_ham_r400_860.sens')],
        'GS_HAM_R400_700': [dict(std_file='spec1d_**LTT7379*.fits',
                                 sens_file='gemini_gmos_gs_ham_r400_700.sens')]},
    'keck_deimos': {
        '900ZD_LVM_5500': [dict(std_file='spec1d_*Feige110*.fits',
                                sens_file='keck_deimos_900zd_lvm_5500.sens')]},
    'keck_mosfire': {
        'Y_long': [dict(std_file='spec1d_*0064-GD71*.fits',
                        sens_file='keck_mosfire_Y_long.sens')]},
    'keck_lris_red_mark4': {
        'long_600_10000_d680': [dict(std_file='spec1d_*00127-GD153*.fits',
                                     sens_file='keck_lris_red_mark4_long_600_10000_d680.sens')]
        },
    'ldt_deveny': {
        'DV2': [dict(std_file='spec1d_**BD+33d2642**.fits',
                     sens_file='ldt_deveny_dv2.sens')],
        'DV6': [dict(std_file='spec1d**G191-B2B**.fits',
                     sens_file='ldt_deveny_dv6.sens')]
        },
    }


_flux_setup = {
    'shane_kast_blue': {
        '600_4310_d55': [{}]},
    'gemini_gnirs_echelle': {
        '32_SB_SXD': [{}]},
    'gemini_gmos': {
        'GS_HAM_R400_860': [{}]},
    }

_flux = {
    'shane_kast_blue': {
        '600_4310_d55': [{}]},
    'shane_kast_red': {
        '600_7500_d55_ret': [{}]},
    'gemini_gnirs_echelle': {
        '32_SB_SXD': [{}]},
    'gemini_gmos': {
        'GS_HAM_R400_860': [{}],
        'GS_HAM_R400_700': [{}]},
    'keck_deimos': {
        '900ZD_LVM_5500': [{}],
        '600ZD_M_6500': [{}]},
    'ldt_deveny': {
        'DV2': [{}],
        'DV6': [{}]}
    }

_flexure ={
    'keck_deimos': {
        '830G_M_8500': [{}]}
    }

_coadd1d = {
    'shane_kast_blue': {
        '600_4310_d55': [{}]},
    'gemini_gnirs_echelle': {
        '32_SB_SXD': [{}]},
    'gemini_gmos': {
        'GS_HAM_R400_860': [{}]},
    'gemini_gmos': {
        'GS_HAM_R400_700': [{}]},
    }

_coadd2d = {
    'gemini_gnirs_echelle': {
        '32_SB_SXD': [dict(coadd_file=True)]},
    'keck_lris_blue': {
        'multi_600_4000_d560': [dict(coadd_file=True)]},
    'vlt_xshooter': {
        'VIS_manual': [dict(coadd_file=True)]},
    'keck_deimos': {
        '830G_M_9000_dither': [dict(coadd_file=True)]},
    'keck_mosfire': {
        'long2pos1_H': [dict(coadd_file=True)]},
    'keck_nires': {
        'ABBA_wstandard': [dict(coadd_file=True)]},
    'keck_nires': {
        'ABBA_nostandard_faint': [dict(coadd_file=True)]}
    }

_telluric = {
    'gemini_gnirs_echelle': {
        '32_SB_SXD': [dict(coadd_file='pisco_coadd.fits', tell_file=True)]},
    'gemini_gmos': {
        'GS_HAM_R400_700': [dict(coadd_file='FRB180924_opt.fits',
                                 tell_file=True)]},
    }

_collate1d = {
    'keck_deimos': {
        '830G_M_8500':
                   [{'files': ['Science/spec1d_*DE.20100913.22358*.fits'],
                    '--refframe': 'heliocentric',
                    '--wv_rms_thresh': 0.1,
                    '--flux': None}]}
        }

_quick_look = {
    'shane_kast_blue': {
        '600_4310_d55':  [
            # (1) Basic execution with no pre-existing directories.
            #   - All output is to QL_std/
            #   - Calibrations are generated and placed in
            #     QL_std/shane_kast_blue_A/Calibrations
            #   - Science results are in QL_std/b27/Science
            #   - The QL_std/b27/Calibrations directory is a symlink to
            #     QL_std/shane_kast_blue_A/Calibrations
            {'test_name': 'std',
             'files' : ['b1.fits.gz', 'b10.fits.gz', 'b27.fits.gz'],
            },
            # (2) Use the pre-existing calibrations generated during the
            # "reduce" run for this setup.  The QL script is pointed directly to
            # the directory with the calibrations to use.
            #   - All output is to QL_cooked/
            #   - Calibrations are *not* generated
            #   - Science results are in QL_cooked/b27/Science
            #   - The QL_cooked/b27/Calibrations directory is a symlink to
            #     shane_kast_blue_A/Calibrations
            {'test_name': 'cooked',
             'files': ['b1.fits.gz', 'b10.fits.gz', 'b27.fits.gz'],
             '--setup_calib_dir': 'USE_CALIB_DIR',
            },
            # (3) Same as test 2, except that the QL script is pointed to the
            # top-level directory that could potentially house calibrations for
            # multiple setups.  The script has to match the setup used for the
            # science frames to the correct calibrations (although only one
            # setup exists).
            {'test_name': 'match',
             'files': ['b1.fits.gz', 'b10.fits.gz', 'b27.fits.gz'],
             '--parent_calib_dir': 'USE_ARCHIVE_CALIB_DIR',
            },
            # (4) Same as test 2, but process two stacked science frames
            {'test_name': 'multi',
             'files': ['b1.fits.gz', 'b10.fits.gz', 'b27.fits.gz', 'b28.fits.gz'],
             '--setup_calib_dir': 'USE_CALIB_DIR',
            },
            # (5) Same as test 2, but change the boxcar extaction width
            {'test_name': 'boxcar',
             'files': ['b1.fits.gz', 'b10.fits.gz', 'b27.fits.gz'],
             '--setup_calib_dir': 'USE_CALIB_DIR',
             '--boxcar_radius': 2.,
            },
        ],
    },
    'shane_kast_red': {
        '600_7500_d57': [
            {'files': ['r122.fits'],
              '--setup_calib_dir': 'USE_CALIB_DIR',
            },
        ]
    },
    'keck_lris_red': {
        'long_600_7500_d560': [
            {'test_name': 'det',
                'files': ['LR.20160216.40478.fits.gz'],
              '--det': 2,
              '--setup_calib_dir': 'USE_CALIB_DIR',
            },
        ],
    },
    'keck_deimos': {
        '600ZD_M_6500': [
            {'test_name': 'maskID', # Run with maskID
                'files': ['d1010_0056.fits.gz'],
              '--maskID': 958454,
              '--setup_calib_dir': 'USE_CALIB_DIR',
            },
            {'test_name': 'slitspatnum', # Run with slitspatnum
                'files': ['d1010_0056.fits.gz'],
              '--slitspatnum': 'MSC02:368,MSC02:452',
              '--setup_calib_dir': 'USE_CALIB_DIR',
            },
        ]
    },
    'keck_mosfire': {
        'J_multi': [
            {'files': ['m191014_0170.fits'],
              '--setup_calib_dir': 'USE_CALIB_DIR',
            }
        ],
        'Y_long': [
            # (1) Run without using archived calibrations
            # NOTE: This takes ~8min, so not really quick...
#            {'test_name': 'raw',
#             'files': ['m191119_0027.fits', 'm191119_0037.fits', 'm191120_0043.fits',
#                       'm191120_0044.fits', 'm191120_0045.fits', 'm191120_0046.fits'],
#             '--coadd': None, '--spec_samp_fact': 2.0, '--spat_samp_fact': 2.0,
#            },
            # (2) Run with archived calibrations
            # NOTE: Takes ~2min
            {'test_name': 'arc',
             'files': ['m191120_0043.fits', 'm191120_0044.fits',
                       'm191120_0045.fits', 'm191120_0046.fits'],
             '--parent_calib_dir': 'USE_ARCHIVE_CALIB_DIR',
             '--coadd': None, '--spec_samp_fact': 2.0, '--spat_samp_fact': 2.0,
            },
        ]
    },
    'keck_lris_red_mark4': {
        'long_600_10000_d680': [
            {'files': ['r220127_00123.fits', 'r220127_00124.fits'],
              '--coadd': None, '--spec_samp_fact': 2.0, '--spat_samp_fact': 2.0,
              '--setup_calib_dir': 'USE_CALIB_DIR',
            }
        ],
    },
    'keck_nires': {
        'ABpat_wstandard': [
            # (1) Generate the calibrations from scratch and just reduce a
            # single frame
            {'test_name': 'one',
             'files': ['NR.20191211.07572.fits',    # flat
                       'NR.20191211.07688.fits',    # lamp-off flat
                       'NR.20191211.26257.fits'],   # science (arc,tilt)
            },
            # (2) Use existing calibrations and just reduce a single frame
            {'test_name': 'cooked',
             'files': ['NR.20191211.26257.fits'],
             '--setup_calib_dir': 'USE_CALIB_DIR',
            },
            # (3) Use exising calibrations and reduce both a standard and a
            # single science frame
            {'test_name': 'std',
             'files': ['NR.20191211.26257.fits',    # science (arc,tilt)
                       'NR.20191211.27199.fits'],   # standard
             '--setup_calib_dir': 'USE_CALIB_DIR',
            },
            # (4) Use exising calibrations and reduce both a standard and an
            # AB dither sequence
            {'test_name': 'ABstd',
             'files': ['NR.20191211.26257.fits',    # science: A
                       'NR.20191211.26611.fits',    # science: B
                       'NR.20191211.27199.fits'],   # standard
             '--snr_thresh': 5,
             '--setup_calib_dir': 'USE_CALIB_DIR',
            },
            # (5) Use existing global calibrations
            {'test_name': 'ABarc',
             'files': ['NR.20191211.26257.fits',    # science: A
                       'NR.20191211.26611.fits'],   # science: B
             '--snr_thresh': 5,
             '--parent_calib_dir': 'USE_ARCHIVE_CALIB_DIR',
            },
        ]
    },
    }


# The order of these tests in all_tests determine the order they run
# in for the setup. So that tests that depend on previous tests must
# be in the right order. e.g. PypeItSetupTest must come before
# PypeItReduceTest and PypeItSensFuncTest must come before
# PypeItFluxTest.
#
all_tests = [{'factory': pypeit_tests.PypeItSetupTest,
              'type':    TestPhase.PREP,
              'setups':  _pypeit_setup},
             {'factory': pypeit_tests.PypeItReduceTest,
              'type':    TestPhase.REDUCE,
              'setups':  _reduce_setups},
             {'factory': pypeit_tests.PypeItReduceTest,
              'type':    TestPhase.REDUCE,
              'setups':  _additional_reduce},
             {'factory': pypeit_tests.PypeItSensFuncTest,
              'type':    TestPhase.AFTERBURN,
              'setups':  _sensfunc},
             {'factory': pypeit_tests.PypeItFluxSetupTest,
              'type':    TestPhase.AFTERBURN,
              'setups':  _flux_setup},
             {'factory': pypeit_tests.PypeItFluxTest,
              'type':    TestPhase.AFTERBURN,
              'setups':  _flux},
             {'factory': pypeit_tests.PypeItFlexureTest,
              'type':    TestPhase.AFTERBURN,
              'setups':  _flexure},
             {'factory': pypeit_tests.PypeItCollate1DTest,
              'type':    TestPhase.AFTERBURN,
              'setups':  _collate1d},
             {'factory': pypeit_tests.PypeItCoadd1DTest,
              'type':    TestPhase.AFTERBURN,
              'setups':  _coadd1d},
             {'factory': pypeit_tests.PypeItCoadd2DTest,
              'type':    TestPhase.AFTERBURN,
              'setups':  _coadd2d},
             {'factory': pypeit_tests.PypeItTelluricTest,
              'type':    TestPhase.AFTERBURN,
              'setups':  _telluric},
             {'factory': pypeit_tests.PypeItQuickLookTest,
              'type':    TestPhase.QL,
              'setups':  _quick_look},
             ]
