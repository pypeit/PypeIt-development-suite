
import sys
import argparse
import copy
import os

import numpy as np

from pypeit.sensfunc import SensFunc
from pypeit.core import fitting
from pypeit.core.wavecal import wvutils
import pypeit.io

from pypeit import utils
from matplotlib import pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages


def stitch_with_polyfit(wave_1, zp_fit_1, gpm_1,
                        wave_2, zp_fit_2, gpm_2,
                        fitting_wave, fitting_zp_fit, order,
                        use_polyfit_min, use_polyfit_max):

    fitc = np.polynomial.polynomial.polyfit(fitting_wave, fitting_zp_fit, order)
    fit_info = [np.min(fitting_wave), np.max(fitting_wave), fitc]

    if use_polyfit_max is None or use_polyfit_min is None:
        # Stitch without polyfit, but still return the fit_info
        # This is used to generate a plot from which the min/max can be chosen
        combined_wave = np.concatenate((wave_1, wave_2))

        combined_zp_fit = np.concatenate((zp_fit_1, zp_fit_2))
    
        combined_zp_fit_gpm = np.concatenate((gpm_1, gpm_2))

    else:
        use_polyfit_1_mask = wave_1 > use_polyfit_min
        use_polyfit_2_mask = wave_2 < use_polyfit_max

        wave_grid_min = wave_1[use_polyfit_1_mask][0]
        wave_grid_max = wave_2[use_polyfit_2_mask][-1]

        joining_wave = wvutils.get_wave_grid(fitting_wave,
                                            wave_method='linear',
                                            wave_grid_min=wave_grid_min,
                                            wave_grid_max=wave_grid_max,
                                            spec_samp_fact=1.0)[0]

        joining_zp_fit = np.polynomial.polynomial.polyval(joining_wave, fitc)

        joining_gpm = np.ones_like(joining_wave, dtype=bool)

        combined_wave = np.concatenate((wave_1[np.logical_not(use_polyfit_1_mask)],
                                    joining_wave,
                                    wave_2[np.logical_not(use_polyfit_2_mask)]))

        combined_zp_fit = np.concatenate((zp_fit_1[np.logical_not(use_polyfit_1_mask)],
                                        joining_zp_fit,
                                        zp_fit_2[np.logical_not(use_polyfit_2_mask)]))
    
        combined_zp_fit_gpm = np.concatenate((gpm_1[np.logical_not(use_polyfit_1_mask)],
                                            joining_gpm,
                                            gpm_2[np.logical_not(use_polyfit_2_mask)]))

    return combined_wave, combined_zp_fit, combined_zp_fit_gpm, fit_info


def stitch_sensfunc(args, sflist):
    if args.grating == '1200G':
        return stitch_1200G_sensfunc(sflist)
    elif args.grating == '1200B':
        return stitch_1200B_sensfunc(sflist)
    elif args.grating == '600ZD':
        return stitch_600ZD_sensfunc(sflist)
    elif args.grating == '900ZD':
        return stitch_900ZD_sensfunc(sflist)
    elif args.grating == '830G':
        return stitch_830G_sensfunc(sflist)
    else:
        raise ValueError(f"{args.grating} not Implemented (yet)")


def stitch_1200G_sensfunc(sflist):
    """
    Stitch together sensfunc files generated from throughput data from Gret Wriths scripts. Specifically
    The following files were used:
        sens_from_gdw_data/extract/1200G/sens_2005aug27_d0827_0046.fits
        sens_from_gdw_data/extract/1200G/sens_2005aug27_d0827_0047.fits
        sens_from_gdw_data/extract/1200G/sens_2005aug27_d0827_0048.fits

    """


    # Sanity check the first sensfunc file by masking any 0 values and ensuring there are no NaNs or +/- inf values
    file1_det3_nzpm = (sflist[0].sens['SENS_ZEROPOINT_FIT'][0] != 0) & (sflist[0].sens['SENS_WAVE'][0] != 0)
    file1_det7_nzpm = (sflist[0].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[0].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm]).all()

    # Take everything in the first sensfunc, nudging the det 7 sensfunc up a bit to
    # match det 3

    combined_wave = np.concatenate((sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm], sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm]))
    combined_zp_fit = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm], sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm]+0.0354))
    combined_zp_fit_gpm = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][0][file1_det3_nzpm], sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][1][file1_det7_nzpm]))

    # Sanity check the second sensfunc file by masking any 0 values and ensuring there are no NaNs or +/- inf values
    file2_det3_nzpm = (sflist[1].sens['SENS_ZEROPOINT_FIT'][0] != 0) & (sflist[1].sens['SENS_WAVE'][0] != 0)
    file2_det7_nzpm = (sflist[1].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[1].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_ZEROPOINT_FIT'][0][file2_det3_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm]).all()

    # For detector 3 in the second file, take everything that doesn't overlap detector 7 in the first
    non_overlap = sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm] > np.max(sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm])
    combined_wave = np.concatenate((combined_wave, sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm][non_overlap]))
    combined_zp_fit = np.concatenate((combined_zp_fit, sflist[1].sens['SENS_ZEROPOINT_FIT'][0][file2_det3_nzpm][non_overlap]+0.0983))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm, sflist[1].sens['SENS_ZEROPOINT_FIT_GPM'][0][file2_det3_nzpm][non_overlap]))

    # For detecter 7 in the second file, nudge it down to match the  previous sensfunc
    file2_det7_offset = -0.0232
    combined_wave = np.concatenate((combined_wave, sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm]))
    combined_zp_fit = np.concatenate((combined_zp_fit, sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm]+file2_det7_offset))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm, sflist[1].sens['SENS_ZEROPOINT_FIT_GPM'][1][file2_det7_nzpm]))

    # Sanity check the third sensfunc file detector 7 by masking any 0 values and ensuring there are no NaNs or +/- inf values
    file3_det7_nzpm = (sflist[2].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[2].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm]).all()
    assert np.isfinite(sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm]).all()

    # Combine the non overlapping part of file 3 det 7's sensfunc. 
    # This uses polynomial fitting to smooth the discontinuity between the two sensfuncs.
    # File 3 det 7's sensfunc is also nudged up a bit to match the previous sensfunc
    # better
    non_overlap = sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm] > np.max(combined_wave)
    file3_det7_offset = 0.0899
    fitting_wave = np.concatenate((sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm], sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm][non_overlap]))
    fitting_zp_fit =  np.concatenate((sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm]+file2_det7_offset, sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm][non_overlap]+file3_det7_offset))

    # The region to overwite with the polyfit values
    polyfit1_write_min = 8074.5
    polyfit1_write_max = 8709.5

    (combined_wave, combined_zp_fit, 
     combined_zp_fit_gpm, fit_info) = stitch_with_polyfit(
                                        combined_wave,
                                        combined_zp_fit,
                                        combined_zp_fit_gpm,
                                        sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm][non_overlap],
                                        sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm][non_overlap] + file3_det7_offset,
                                        sflist[2].sens['SENS_ZEROPOINT_FIT_GPM'][1][file3_det7_nzpm][non_overlap],
                                        fitting_wave, fitting_zp_fit, 30, polyfit1_write_min, polyfit1_write_max)

    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    polyfit_areas = (combined_wave > polyfit1_write_min) & (combined_wave < polyfit1_write_max)

    return (combined_wave, combined_zp_fit, combined_zp_fit_gpm, polyfit_areas)



def stitch_1200B_sensfunc(sflist):
    """
    Stitch together sensfunc files generated from throughput data from Gret Wriths scripts. Specifically
    The following files were used:
        sens_from_gdw_data/extract/1200B/sens_2017oct03_d1003_0055.fits
        sens_from_gdw_data/extract/1200B/sens_2017oct03_d1003_0056.fits
        sens_from_gdw_data/extract/1200B/sens_2017oct03_d1003_0057.fits

    """

    # Sanity check both sensors from the first sensfunc file

    file1_det3_nzpm = (sflist[0].sens['SENS_ZEROPOINT_FIT'][0] != 0) & (sflist[0].sens['SENS_WAVE'][0] != 0)
    file1_det7_nzpm = (sflist[0].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[0].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm]).all()

    ##### Numbers derived from plotting the sensfuncs:

    # The discontinuity around the det 3/7 edge in sensfunc file 1
    polyfit1_det_min = 5050.
    polyfit1_det_max = 5450.

    # The region to overwrite with the polyfit values
    polyfit1_write_min = 4905.35
    polyfit1_write_max = 5554.7

    # Use polynomial fitting to smooth the join between file 1 det 3 and 7.  
    # The fitting uses a fitting_wave/zp_fit that masks off the region around
    # the detector boundry to get a polynomial that fits the join area better
    file1_det3_fitting_mask = sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm] < polyfit1_det_min
    file1_det7_fitting_mask = sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm] > polyfit1_det_max
    fitting_wave = np.concatenate((sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm][file1_det3_fitting_mask],
                                   sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm][file1_det7_fitting_mask]))
    fitting_zp_fit = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm][file1_det3_fitting_mask], 
                                     sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm][file1_det7_fitting_mask]))


    (combined_wave, combined_zp_fit, 
     combined_zp_fit_gpm, fit_info) = stitch_with_polyfit(
                                        sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm], 
                                        sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][0][file1_det3_nzpm],
                                        sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][1][file1_det7_nzpm],
                                        fitting_wave, fitting_zp_fit, 15, polyfit1_write_min, polyfit1_write_max)


    # Sanity check both sensors from the second sensfunc file

    file2_det3_nzpm = (sflist[1].sens['SENS_ZEROPOINT_FIT'][0] != 0) & (sflist[1].sens['SENS_WAVE'][0] != 0)
    file2_det7_nzpm = (sflist[1].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[1].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_ZEROPOINT_FIT'][0][file2_det3_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm]).all()

    # For detector 3 in the second, take everything that doesn't overlap detector 7 in the first
    # Move it up a bit and it's a smooth enought fit that no polynomial fitting is needed
    
    non_overlap = sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm] > np.max(sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm])
    combined_wave = np.concatenate((combined_wave, sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm][non_overlap]))
    combined_zp_fit = np.concatenate((combined_zp_fit, sflist[1].sens['SENS_ZEROPOINT_FIT'][0][file2_det3_nzpm][non_overlap]+0.0474))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm, sflist[1].sens['SENS_ZEROPOINT_FIT_GPM'][0][file2_det3_nzpm][non_overlap]))

    # Take all of detecter 7 in the second file. We nudge it down a bit match det 3 better and use a polynomial fit
    # to smooth the disonctinuity 

    # The discontinuity around the det 3/7 edge in sensfunc file 2
    polyfit2_det_min = 6850.
    polyfit2_det_max = 7200.

    # The region to overwrite with the polyfit values
    polyfit2_write_min = 6747.3
    polyfit2_write_max = 7289.15

    file2_det7_offset = -0.0845

    file2_det3_fitting_mask = sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm][non_overlap] < polyfit2_det_min
    file2_det7_fitting_mask = sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm] > polyfit2_det_max
    fitting_wave = np.concatenate((sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm][non_overlap][file2_det3_fitting_mask],
                                   sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm][file2_det7_fitting_mask]))
    fitting_zp_fit = np.concatenate((sflist[1].sens['SENS_ZEROPOINT_FIT'][0][file2_det3_nzpm][non_overlap][file2_det3_fitting_mask]+0.0474, 
                                     sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm][file2_det7_fitting_mask] + file2_det7_offset))


    (combined_wave, combined_zp_fit, 
     combined_zp_fit_gpm, fit_info) = stitch_with_polyfit(
                                        combined_wave,
                                        combined_zp_fit,
                                        combined_zp_fit_gpm,
                                        sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm],
                                        sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm] + file2_det7_offset,
                                        sflist[1].sens['SENS_ZEROPOINT_FIT_GPM'][1][file2_det7_nzpm],
                                        fitting_wave, fitting_zp_fit, 15, polyfit2_write_min, polyfit2_write_max)

    

    # For the third sensfunc file, ignore detector 3 because it overlaps everything.

    # Sanity check detector 7 from the third file
    file3_det7_nzpm = (sflist[2].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[2].sens['SENS_WAVE'][1] != 0)
    assert np.isfinite(sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm]).all()
    assert np.isfinite(sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm]).all()


    # Take the non overlapping parts of detector 7 and do a polynomail fit to smooth over the join

    non_overlap = sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm] > np.max(sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm])


    # The discontinuity around the joining of file2 det 7 and file 3 det 7
    polyfit3_det_min = 8250
    polyfit3_det_max = 8450

    # The region to overwrite with the polyfit values
    polyfit3_write_min = 8211.1
    polyfit3_write_max = 8544.1

    file2_det7_fitting_mask = sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm] < polyfit3_det_min
    file3_det7_fitting_mask = sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm][non_overlap] > polyfit3_det_max
    fitting_wave = np.concatenate((sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm][file2_det7_fitting_mask],
                                   sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm][non_overlap][file3_det7_fitting_mask]))
    fitting_zp_fit = np.concatenate((sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm][file2_det7_fitting_mask]-0.0845,
                                     sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm][non_overlap][file3_det7_fitting_mask]))


    (combined_wave, combined_zp_fit, 
     combined_zp_fit_gpm, fit_info) = stitch_with_polyfit(
                                        combined_wave,
                                        combined_zp_fit,
                                        combined_zp_fit_gpm,
                                        sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm][non_overlap],
                                        sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm][non_overlap],
                                        sflist[2].sens['SENS_ZEROPOINT_FIT_GPM'][1][file3_det7_nzpm][non_overlap],
                                        fitting_wave, fitting_zp_fit, 15, polyfit3_write_min, polyfit3_write_max)


    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    combined_edges = ((combined_wave >= polyfit1_write_min) & (combined_wave <= polyfit1_write_max)) | ((combined_wave >= polyfit2_write_min) & (combined_wave <= polyfit2_write_max)) | ((combined_wave >= polyfit3_write_min) & (combined_wave <= polyfit3_write_max))
    return (combined_wave, combined_zp_fit, combined_zp_fit_gpm, combined_edges)

def stitch_900ZD_sensfunc(sflist):

    """
    Stitch together sensfunc files generated from throughput data from Gret Wriths scripts. Specifically
    The following files were used:
        sens_from_gdw_data/extract/900ZD/sens_2010oct06_d1006_0138.fits
        sens_from_gdw_data/extract/900ZD/sens_2010oct06_d1006_0139.fits
        sens_from_gdw_data/extract/900ZD/sens_2010oct06_d1006_0140.fits
    """

    # Sanity check the first sensfunc file by masking any 0 values and ensuring there are no NaNs or +/- inf values
    file1_det3_nzpm = (sflist[0].sens['SENS_ZEROPOINT_FIT'][0] != 0) & (sflist[0].sens['SENS_WAVE'][0] != 0)
    file1_det7_nzpm = (sflist[0].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[0].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm]).all()

    # Combine file 1 detectors 3 and 7 with a polynomial fit to smooth out the discontinuity.

    # Mask the fitting around the detector boundary and towards the outer ends of the sensfuncs to get
    # the best fit
    fitting_wave = np.concatenate((sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm], sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm]))
    fitting_zp_fit = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm], sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm]))

    fitting_mask = ((fitting_wave > 4629.4) & (fitting_wave < 5324.5)) | ((fitting_wave > 5530.4) & (fitting_wave < 9000))  


    # The region to overwrite with the polyfit values
    polyfit1_write_min = 4910.2
    polyfit1_write_max = 5805.2

    (combined_wave, combined_zp_fit, 
     combined_zp_fit_gpm, fit_info) = stitch_with_polyfit(
                                        sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm], 
                                        sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][0][file1_det3_nzpm],
                                        sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][1][file1_det7_nzpm],
                                        fitting_wave[fitting_mask], fitting_zp_fit[fitting_mask], 15, 
                                        polyfit1_write_min, polyfit1_write_max)

    
    
    # Now add the parts of file 3, det3 that don't overlap file 1 det 7.
    # We're skipping the second file entirely because of overlaps
    # We also translate up file 3 det 3 up to match file 1 det 7
    file3_det3_nzpm = (sflist[2].sens['SENS_ZEROPOINT_FIT'][0] != 0) & (sflist[2].sens['SENS_WAVE'][0] != 0)
    assert np.isfinite(sflist[2].sens['SENS_WAVE'][0][file3_det3_nzpm]).all()
    assert np.isfinite(sflist[2].sens['SENS_ZEROPOINT_FIT'][0][file3_det3_nzpm]).all()


    non_overlap = sflist[2].sens['SENS_WAVE'][0][file3_det3_nzpm] > np.max(sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm])

    file3_det3_offset = 0.10533590352139655
    combined_wave = np.concatenate((combined_wave, sflist[2].sens['SENS_WAVE'][0][file3_det3_nzpm][non_overlap]))
    combined_zp_fit = np.concatenate((combined_zp_fit, sflist[2].sens['SENS_ZEROPOINT_FIT'][0][file3_det3_nzpm][non_overlap]+file3_det3_offset))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm, sflist[2].sens['SENS_ZEROPOINT_FIT_GPM'][0][file3_det3_nzpm][non_overlap]))

    # Now add file 3 det 7, moving it down to match det 3 and using a polynomial fit to smooth the 
    # discontinuities between det 3 and 7

    file3_det7_nzpm = (sflist[2].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[2].sens['SENS_WAVE'][1] != 0)
    assert np.isfinite(sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm]).all()
    assert np.isfinite(sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm]).all()

    file3_det7_offset = -0.02687150393769855

    fitting_wave=np.concatenate((sflist[2].sens['SENS_WAVE'][0][file3_det3_nzpm][non_overlap], sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm]))
    fitting_mask = (fitting_wave > 7600) & (fitting_wave < 8400) 
    fitting_zp_fit = np.concatenate((sflist[2].sens['SENS_ZEROPOINT_FIT'][0][file3_det3_nzpm][non_overlap]+file3_det3_offset, sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm]+file3_det7_offset))

    polyfit2_write_min = 7876.
    polyfit2_write_max = 8055.7

    (combined_wave, combined_zp_fit, 
     combined_zp_fit_gpm, fit_info) = stitch_with_polyfit(
                                        combined_wave,
                                        combined_zp_fit, 
                                        combined_zp_fit_gpm,
                                        sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm],
                                        sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm]+file3_det7_offset,
                                        sflist[2].sens['SENS_ZEROPOINT_FIT_GPM'][1][file3_det7_nzpm],
                                        fitting_wave[fitting_mask], fitting_zp_fit[fitting_mask], 30, polyfit2_write_min, polyfit2_write_max)

    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    polyfit_areas = ((combined_wave > polyfit1_write_min) & (combined_wave < polyfit1_write_max)) | ((combined_wave > polyfit2_write_min) & (combined_wave < polyfit2_write_max))


    return (combined_wave, combined_zp_fit, combined_zp_fit_gpm, polyfit_areas)

def stitch_600ZD_sensfunc(sflist):
    """
    Stitch together sensfunc files generated from throughput data from Gret Wriths scripts. Specifically
    The following files were used:
    sens_from_gdw_data/extract/600ZD/sens_2010sep24_d0924_0008.fits
    sens_from_gdw_data/extract/600ZD/sens_2010sep24_d0924_0009.fits
    sens_from_gdw_data/extract/600ZD/sens_2010sep24_d0924_0010.fits

    """

    # Stitch file1 det3 and det7 with a polynomial fit to smooth the discontinuity between them

    # Eliminate any 0 values in SENS_WAVE or SENS_ZEROPOINT_FIT
    file1_det3_nzpm = (sflist[0].sens['SENS_ZEROPOINT_FIT'][0] != 0) & (sflist[0].sens['SENS_WAVE'][0] != 0)
    file1_det7_nzpm = (sflist[0].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[0].sens['SENS_WAVE'][1] != 0)

    # Ensure there are no NaNs or +/- inf values
    assert np.isfinite(sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm]).all()

    # The polynomial fit will be done on a subset of the sensfuncs to improve the fit and eliminate
    # odd behavior near the detector boundary. These numbers were chosen by examining plots of the
    # individual sensfuncs
    file1_det3_fitting_mask = (sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm] > 5000.) & (sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm] < 5712.)
    
    fitting_wave = np.concatenate((sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm][file1_det3_fitting_mask], sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm]))
    fitting_zp_fit = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm][file1_det3_fitting_mask], sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm]))

    # The region to overwrite with the polyfit values
    polyfit1_write_min = 5379.55
    polyfit1_write_max = 6336.25

    (combined_wave, combined_zp_fit, 
     combined_zp_fit_gpm, fit_info) = stitch_with_polyfit(
                                        sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm], 
                                        sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][0][file1_det3_nzpm],
                                        sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][1][file1_det7_nzpm],
                                        fitting_wave, fitting_zp_fit, 15, polyfit1_write_min, polyfit1_write_max)


    # Append the sensfunc from file 2 det 7.
    # The two sensfuncs are stitched together at a point that avoids odd behavior at the detector boundary
    # in file 1 det 7's sensfunc.  File 2 det 7's sensfunc is then translated up a bit to be a better match

    file2_det7_nzpm = (sflist[1].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[1].sens['SENS_WAVE'][1] != 0)
    assert np.isfinite(sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm]).all()

    combined_wave_idx = combined_wave < 7558.1
    file2_det_7_wave_idx = sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm] > 7558.1
    combined_wave = np.concatenate((combined_wave[combined_wave_idx], sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm][file2_det_7_wave_idx]))
    combined_zp_fit = np.concatenate((combined_zp_fit[combined_wave_idx], sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm][file2_det_7_wave_idx] + 0.0156))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm[combined_wave_idx], sflist[1].sens['SENS_ZEROPOINT_FIT_GPM'][1][file2_det7_nzpm][file2_det_7_wave_idx]))


    # Append the portion of file3 detector 7 that does not overlap file 2 detector 7. It is also 
    # translated up to match
    file3_det7_nzpm = (sflist[2].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[2].sens['SENS_WAVE'][1] != 0)
    assert np.isfinite(sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm]).all()
    assert np.isfinite(sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm]).all()

    file3_det_7_wave_idx = sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm] > np.max(combined_wave)

    combined_wave = np.concatenate((combined_wave, sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm][file3_det_7_wave_idx]))
    combined_zp_fit = np.concatenate((combined_zp_fit, sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm][file3_det_7_wave_idx]+0.0663))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm, sflist[2].sens['SENS_ZEROPOINT_FIT_GPM'][1][file3_det7_nzpm][file3_det_7_wave_idx]))

    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    polyfit_areas = (combined_wave > polyfit1_write_min) & (combined_wave < polyfit1_write_max)
    
    return (combined_wave, combined_zp_fit, combined_zp_fit_gpm, polyfit_areas)

def stitch_830G_sensfunc(sflist):

    """
    Stitch together sensfunc files generated from throughput data from Gret Wriths scripts. Specifically
    The following files were used:
        sens_from_gdw_data/extract/830G/sens_2010oct06_d1006_0135.fits
        sens_from_gdw_data/extract/830G/sens_2010oct06_d1006_0136.fits
        sens_from_gdw_data/extract/830G/sens_2010oct06_d1006_0137.fits
    """

    # Sanity check the first sensfunc file by masking any 0 values and ensuring there are no NaNs or +/- inf values
    file1_det3_nzpm = (sflist[0].sens['SENS_ZEROPOINT_FIT'][0] != 0) & (sflist[0].sens['SENS_WAVE'][0] != 0)
    file1_det7_nzpm = (sflist[0].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[0].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm]).all()



    # Take everything in the first sensfunc but mask data around the transition between 3 to 7 to 
    # let the polynomial fit smooth it out.  
    fitting_wave = np.concatenate((sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm], sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm]))
    fitting_zp_fit = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm], sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm]))
    fitting_mask = np.logical_not((fitting_wave > 5305.8) & (fitting_wave < 5528.8))

    # The region to overwrite with the polyfit values
    polyfit1_write_min = 5196.4
    polyfit1_write_max = 5782.5

    (combined_wave, combined_zp_fit, 
     combined_zp_fit_gpm, fit_info) = stitch_with_polyfit(
                                        sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm], 
                                        sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][0][file1_det3_nzpm],
                                        sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm],
                                        sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][1][file1_det7_nzpm],
                                        fitting_wave[fitting_mask], fitting_zp_fit[fitting_mask], 15, 
                                        polyfit1_write_min, polyfit1_write_max)

    # Sanity check file 2, det 7
    file2_det7_nzpm = (sflist[1].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[1].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm]).all()

    # Add file 2, det 7, using a polynomial fit to join the sensfuncs

    fitting_wave = np.concatenate((sflist[0].sens['SENS_WAVE'][1][file1_det7_nzpm], sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm]))
    fitting_zp_fit = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT'][1][file1_det7_nzpm], sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm]))
    fitting_mask = (fitting_wave > 6000) & (fitting_wave < 7500)

    # Truncate file 2 det 7's sensfunc at around the intersection with the next sensfunc
    file2_det7_trunc = 8742.8
    file2_det7_trunc_mask = sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm] < file2_det7_trunc

    # The region that will be overwritten by the polynomial fit
    polyfit2_write_min = 6660.2
    polyfit2_write_max = 7160.2

    (combined_wave, combined_zp_fit, 
     combined_zp_fit_gpm, fit_info) = stitch_with_polyfit(
                                        combined_wave,
                                        combined_zp_fit, 
                                        combined_zp_fit_gpm,
                                        sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm][file2_det7_trunc_mask],
                                        sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm][file2_det7_trunc_mask],
                                        sflist[1].sens['SENS_ZEROPOINT_FIT_GPM'][1][file2_det7_nzpm][file2_det7_trunc_mask],
                                        fitting_wave[fitting_mask], fitting_zp_fit[fitting_mask], 15, 
                                        polyfit2_write_min, polyfit2_write_max)

    # For the third sens func, ignore detector 3 because it overlaps everything.
    # Take the non overlapping parts of detector 7 and translate it up 
    file3_det7_nzpm = (sflist[2].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[2].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm]).all()
    assert np.isfinite(sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm]).all()

    non_overlap = sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm] > file2_det7_trunc
    combined_wave = np.concatenate((combined_wave, sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm][non_overlap]))
    combined_zp_fit = np.concatenate((combined_zp_fit, sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm][non_overlap]+0.0483))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm, sflist[2].sens['SENS_ZEROPOINT_FIT_GPM'][1][file3_det7_nzpm][non_overlap]))

    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    polyfit_areas = ((combined_wave > polyfit1_write_min) & (combined_wave < polyfit1_write_max)) | ((combined_wave > polyfit2_write_min) & (combined_wave < polyfit2_write_max))

    return (combined_wave, combined_zp_fit, combined_zp_fit_gpm, polyfit_areas)

def write_stitched_sensfunc(sflist, grating, combined_wave, combined_zp_fit,combined_zp_fit_gpm):
    
    newsf = copy.deepcopy(sflist[0])
    newsf.sens =  None
    newsens = SensFunc.empty_sensfunc_table(1, len(combined_wave))
    newsens['SENS_WAVE'] = combined_wave
    newsens['SENS_ZEROPOINT_FIT'] = combined_zp_fit
    newsens['SENS_ZEROPOINT_FIT_GPM'] = combined_zp_fit_gpm
    newsf.sens = newsens

    newsf.spec1df = None
    newsf.std_name = None
    newsf.std_cal = None
    newsf.std_ra = None
    newsf.std_dec = None
    newsf.airmass = None
    newsf.telluric = None
    newsf.wave = np.empty((combined_wave.size,1))
    newsf.wave[:,0] = combined_wave
    newsf.zeropoint = np.empty((combined_zp_fit.size,1))
    newsf.zeropoint[:,0] = combined_zp_fit
    newsf.throughput = None
    file_name = f"keck_deimos_{grating}_sensfunc.fits"

    hdr = pypeit.io.initialize_header(primary=True)
    hdr['HISTORY'] = "This DEIMOS sensfunc was stitched together from the following source files."
    # Trim off spec1d_ from base of each spec1d file
    hdr['HISTORY'] = os.path.basename(sflist[0].spec1df)[7:]
    hdr['HISTORY'] = os.path.basename(sflist[1].spec1df)[7:]
    hdr['HISTORY'] = os.path.basename(sflist[2].spec1df)[7:]
    
    newsf.to_file(file_name, primary_hdr = hdr, overwrite=True)
    return (newsf, file_name)

def build_figure():

    utils.pyplot_rcparams()
    fig = plot.figure(figsize=(12,8))
    axis = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    return fig, axis

def create_plot(axis, color, label, x, y, linewidth=2.5, marker = '', linestyle='solid'):

    axis.plot(x, y, color=color, linewidth=linewidth, marker=marker, linestyle=linestyle, label=label)
    xmin = (0.98*x.min())
    #xmin = 3000
    xmax = (1.02*x.max())
    #ymin = 14.0
    ymin = 16.0
    ymax = (1.05*y.max())
    #ymax = 22

    return xmin, xmax, ymin, ymax


def setup_axes(axis,xmin, ymin, xmax, ymax, title):
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    axis.legend()
    axis.set_xlabel('Wavelength (Angstroms)')
    axis.set_ylabel('Zeropoint Fit (AB mag)')
    axis.set_title('PypeIt stitched SensFunc ' + title)

def plot_stitch_results(sf, polyfit_areas, sflist, filename):

    fig, axis = build_figure()

    sens_gpm = sf.sens['SENS_ZEROPOINT_FIT_GPM'][0]
    bpm = np.logical_not(sens_gpm)

    num_plots = 2

    if np.sum(bpm) != 0: # Only plot bad pixels if there are some
        num_plots += 1        

    if sflist is not None:
        num_plots += (2*len(sflist))
    xmin = np.zeros(num_plots)
    ymin = np.zeros(num_plots)
    xmax = np.zeros(num_plots)
    ymax = np.zeros(num_plots)

    i = 0

    # Plot the translated/stitched zero points without the "bad pixels" marked when generating sensfuncs
    x = sf.wave[sens_gpm]
    y = sf.zeropoint[sens_gpm]
    (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, (1, .7, .7, .5), f'Translated/stitched SENS_ZEROPOINT_FIT', x,y, marker='.', linestyle='none')
    i+=1

    # Plot the bad pixels
    if np.sum(bpm) != 0: # Only plot bad pixels if there are some
        x = sf.wave[bpm]
        y = sf.zeropoint[bpm]
        (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, "red", f'Inverse SENS_ZEROPOINT_FIT_GPM', x, y, marker='.', linestyle='none')
        i+=1

    # Plot the areas of the zeropoint that came from polynomial fit in blue
    if polyfit_areas is not None:
        x = sf.wave[polyfit_areas]
        y = sf.zeropoint[polyfit_areas]
        (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, "blue", f'Polynomial fit', x, y, marker='.', linestyle='none')
        i+=1

    # Plot the original sensfuncs in light gray in the background
    for bksf in sflist:
        for det in [0, 1]:
            gpm = bksf.sens['SENS_ZEROPOINT_FIT_GPM'][det]
            x = bksf.sens['SENS_WAVE'][det][gpm]
            y = bksf.sens['SENS_ZEROPOINT_FIT'][det][gpm]
            if i == 3:
                label = "Original sensfuncs"
            else:
                label = None
            (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, (.5, .5, .5), label, x,y, linewidth = 1, linestyle='solid')#linestyle='dashed')
            i+=1


    setup_axes(axis, np.min(xmin), np.min(ymin), np.max(xmax),np.max(ymax), filename)


    plot.show()

    pdf_filename = os.path.splitext(filename)[0] + ".pdf"
    with PdfPages(pdf_filename) as pdf:
        pdf.savefig(fig)



def parse_args(options=None, return_parser=False):
    parser = argparse.ArgumentParser(description='Combine DEIMOS sensitivity functions to create a general purpose one.')
    parser.add_argument("grating", type=str, choices=['1200G', '1200B', '600ZD', '830G', '900ZD'])
    parser.add_argument("filelist", type=str, nargs="+",
                        help="List of sensfunc files to combine.")
            

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """ Executes sensitivity function computation.
    """
    sflist = []
    for file in args.filelist:
        sf = SensFunc.from_file(file)
        sf.sensfile = file
        sflist.append(sf)


    sflist.sort(key=lambda x: x.sens['WAVE_MIN'][0])

    if len(sflist)!=3:
        print("Failed. Currently require exaclty 3 sensfunc files.")
        sys.exit(1)


    (combined_wave,  combined_zp_fit, combined_zp_fit_gpm, polyfit_areas) = stitch_sensfunc(args, sflist)

    (newsf, newfile) = write_stitched_sensfunc(sflist, args.grating, combined_wave, combined_zp_fit, combined_zp_fit_gpm)

    plot_stitch_results(newsf, polyfit_areas, sflist, newfile)



def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
