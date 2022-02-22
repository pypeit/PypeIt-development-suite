import numpy as np
from pypeit.core.wavecal import wvutils

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


def stitch_sensfunc(grating, sflist):
    if grating == '1200G':
        return stitch_1200G_sensfunc(sflist)
    elif grating == '1200B':
        return stitch_1200B_sensfunc(sflist)
    elif grating == '600ZD':
        return stitch_600ZD_sensfunc(sflist)
    elif grating == '900ZD':
        return stitch_900ZD_sensfunc(sflist)
    elif grating == '830G':
        return stitch_830G_sensfunc(sflist)
    else:
        raise ValueError(f"{grating} not Implemented (yet)")


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
    file2_det3_nzpm = (sflist[1].sens['SENS_ZEROPOINT_FIT'][0] != 0) & (sflist[1].sens['SENS_WAVE'][0] != 0)

    assert np.isfinite(sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm]).all()
    assert np.isfinite(sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_ZEROPOINT_FIT'][0][file2_det3_nzpm]).all()



    # Take the sensfuncs from first detector in file 1 and in file 2, merging them at the point where they intersect
    
    file1_det3_mask = sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm] < 5255.7
    file2_det3_mask = sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm] > 5255.7

    combined_wave = np.concatenate((sflist[0].sens['SENS_WAVE'][0][file1_det3_nzpm][file1_det3_mask], sflist[1].sens['SENS_WAVE'][0][file2_det3_nzpm][file2_det3_mask])) 
    combined_zp_fit = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT'][0][file1_det3_nzpm][file1_det3_mask], sflist[1].sens['SENS_ZEROPOINT_FIT'][0][file2_det3_nzpm][file2_det3_mask])) 
    combined_zp_fit_gpm = np.concatenate((sflist[0].sens['SENS_ZEROPOINT_FIT_GPM'][0][file1_det3_nzpm][file1_det3_mask], sflist[1].sens['SENS_ZEROPOINT_FIT_GPM'][0][file2_det3_nzpm][file2_det3_mask])) 

    # Now translate up the sensfunc from det 3 in file 3 to match the joined sensfunc and stitch it on at around 6500

    file3_det3_nzpm = (sflist[2].sens['SENS_ZEROPOINT_FIT'][0] != 0) & (sflist[2].sens['SENS_WAVE'][0] != 0)

    assert np.isfinite(sflist[2].sens['SENS_WAVE'][0][file3_det3_nzpm]).all()
    assert np.isfinite(sflist[2].sens['SENS_ZEROPOINT_FIT'][0][file3_det3_nzpm]).all()

    combined_wave_mask = combined_wave < 6500
    file3_det3_mask = sflist[2].sens['SENS_WAVE'][0][file3_det3_nzpm] > np.max(combined_wave[combined_wave_mask])

    combined_wave = np.concatenate((combined_wave[combined_wave_mask], sflist[2].sens['SENS_WAVE'][0][file3_det3_nzpm][file3_det3_mask]))
    combined_zp_fit = np.concatenate((combined_zp_fit[combined_wave_mask], sflist[2].sens['SENS_ZEROPOINT_FIT'][0][file3_det3_nzpm][file3_det3_mask]+0.08742193265786469))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm[combined_wave_mask], sflist[2].sens['SENS_ZEROPOINT_FIT_GPM'][0][file3_det3_nzpm][file3_det3_mask]))

    # Now translate sensfunc 7 from file 2 down to match the end of the combined sensfunc and stitch them to gether at around 7750

    file2_det7_nzpm = (sflist[1].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[1].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm]).all()
    assert np.isfinite(sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm]).all()

    combined_wave_mask = combined_wave < 7750
    file2_det7_mask = sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm] > np.max(combined_wave[combined_wave_mask])

    combined_wave = np.concatenate((combined_wave[combined_wave_mask], sflist[1].sens['SENS_WAVE'][1][file2_det7_nzpm][file2_det7_mask]))
    combined_zp_fit = np.concatenate((combined_zp_fit[combined_wave_mask], sflist[1].sens['SENS_ZEROPOINT_FIT'][1][file2_det7_nzpm][file2_det7_mask]-0.03777170821255993))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm[combined_wave_mask], sflist[1].sens['SENS_ZEROPOINT_FIT_GPM'][1][file2_det7_nzpm][file2_det7_mask]))

    # Finish the stitched sensfunc by adding the non-overlapping portions of det 7 from file 3.
    file3_det7_nzpm = (sflist[2].sens['SENS_ZEROPOINT_FIT'][1] != 0) & (sflist[2].sens['SENS_WAVE'][1] != 0)

    assert np.isfinite(sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm]).all()
    assert np.isfinite(sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm]).all()

    file3_det7_mask = sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm] > np.max(combined_wave)

    combined_wave = np.concatenate((combined_wave, sflist[2].sens['SENS_WAVE'][1][file3_det7_nzpm][file3_det7_mask]))
    combined_zp_fit = np.concatenate((combined_zp_fit, sflist[2].sens['SENS_ZEROPOINT_FIT'][1][file3_det7_nzpm][file3_det7_mask]-0.003240994416337628))
    combined_zp_fit_gpm = np.concatenate((combined_zp_fit_gpm, sflist[2].sens['SENS_ZEROPOINT_FIT_GPM'][1][file3_det7_nzpm][file3_det7_mask]))

    polyfit_areas = None

    return (combined_wave, combined_zp_fit, combined_zp_fit_gpm, polyfit_areas)



def stitch_830G_sensfunc_old(sflist):

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

