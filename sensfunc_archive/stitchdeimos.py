import numpy as np
from stitchutils import sanity_check_sf, stitch_sf_polyfit, stitch_sf_polyfit_old, stitch_sf, truncate_sf, translate_sf, gradient_stitch_sf

from pypeit import msgs

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
        sens_from_gdw_data/extract/1200G/sens_2023jan17_d0117_0102.fits
        sens_from_gdw_data/extract/1200G/sens_2023jan17_d0117_0103.fits
        sens_from_gdw_data/extract/1200G/sens_2023jan17_d0117_0105.fits

    """

    # Sanity check sensfuncs
    # Note we use file 4 because the original set of data had 4 files, although we don't use the third.
    file1_det3_sf = sanity_check_sf(sflist[0].sens[0])
    file1_det7_sf = sanity_check_sf(sflist[0].sens[1])
    file2_det3_sf = sanity_check_sf(sflist[1].sens[0])
    file2_det7_sf = sanity_check_sf(sflist[1].sens[1])
    file4_det3_sf = sanity_check_sf(sflist[2].sens[0])
    file4_det7_sf = sanity_check_sf(sflist[2].sens[1])

    # Combine file 1 detectors 3 and 7
    # Translate detector 7 up to meet detector 3 and stitch at the end of det 3.
    # There's a gap between the sensors. We don't try to fill that gap, we just
    # use the last slope in detector 3 to predict where it's next point would be
    # using the detector 7 wavelengths, and translate the detector 7 sensfunc up to match that.
    det3_end = file1_det3_sf['SENS_WAVE'][-1]
    det3_end_zp = file1_det3_sf['SENS_ZEROPOINT_FIT'][-1]
    det3_slope = (det3_end_zp - file1_det3_sf['SENS_ZEROPOINT_FIT'][-2]) / (det3_end - file1_det3_sf['SENS_WAVE'][-2])
    det7_start = file1_det7_sf['SENS_WAVE'][0]
    det7_start_zp = file1_det7_sf['SENS_ZEROPOINT_FIT'][0]
    det3_to_7_distance = det7_start - det3_end
    new_det7_start_zp = det3_end_zp + (det3_slope*det3_to_7_distance)
    offset = new_det7_start_zp - det7_start_zp

    combined_sf1 = stitch_sf(file1_det3_sf, translate_sf(file1_det7_sf, offset))

    # For the next stitch, the translate file 2, det3 up to meet the previous 
    # file and truncated the overlap. 
    trunc_file2_det3_sf = truncate_sf(file2_det3_sf, min_wave=combined_sf1['SENS_WAVE'][-1])
    offset = combined_sf1['SENS_ZEROPOINT_FIT'][-1] - trunc_file2_det3_sf['SENS_ZEROPOINT_FIT'][0]
    combined_sf2 = stitch_sf(combined_sf1, translate_sf(file2_det3_sf, offset))

    # Translate file 2 detector 7 up to match, again dealing with a gap between sensfuncs
    cf2_end = combined_sf2['SENS_WAVE'][-1]
    cf2_end_zp = combined_sf2['SENS_ZEROPOINT_FIT'][-1]
    cf2_slope = (cf2_end_zp - combined_sf2['SENS_ZEROPOINT_FIT'][-2]) / (cf2_end - combined_sf2['SENS_WAVE'][-2])
    det7_start = file2_det7_sf['SENS_WAVE'][0]
    det7_start_zp = file2_det7_sf['SENS_ZEROPOINT_FIT'][0]
    det3_to_7_distance = det7_start - cf2_end
    new_det7_start_zp = cf2_end_zp + (cf2_slope*det3_to_7_distance)
    offset = new_det7_start_zp - det7_start_zp
    combined_sf3 = stitch_sf(combined_sf2, translate_sf(file2_det7_sf, offset))

    # Truncate and stitch on file 4 detector 3
    trunc_file4_det3_sf = truncate_sf(file4_det3_sf, min_wave=combined_sf3['SENS_WAVE'][-1])
    offset = combined_sf3['SENS_ZEROPOINT_FIT'][-1] - trunc_file4_det3_sf['SENS_ZEROPOINT_FIT'][0]

    combined_sf4 = stitch_sf(combined_sf3, translate_sf(trunc_file4_det3_sf, offset))

    # For file 4 detector 7, we need to do a polynomial fit, but first we translate it up to better
    # match the previous sensfunc. 

    # There's a gap in the wavelengths, use the slope at the end of the previous combined sensfunc
    # to get a better offset for translation
    cf4_end_slope = ( (combined_sf4['SENS_ZEROPOINT_FIT'][-1] - combined_sf4['SENS_ZEROPOINT_FIT'][-2])
                    / (combined_sf4['SENS_WAVE'][-1] - combined_sf4['SENS_WAVE'][-2]))
    cf4_to_f4d7_distance = file4_det7_sf['SENS_WAVE'][0] - combined_sf4['SENS_WAVE'][-1]
    new_f4d7_zp = combined_sf4['SENS_ZEROPOINT_FIT'][-1] + cf4_end_slope * cf4_to_f4d7_distance
    offset = new_f4d7_zp - file4_det7_sf['SENS_ZEROPOINT_FIT'][0]

    # Only use specic regions for fitting
    areas_for_fitting = [(8400, 8900), (9200, 9800)]
    cf4_fitting_mask = (combined_sf4['SENS_WAVE'] >= areas_for_fitting[0][0]) & (combined_sf4['SENS_WAVE'] <= areas_for_fitting[0][1])
    f4d7_fitting_mask = (file4_det7_sf['SENS_WAVE'] >= areas_for_fitting[1][0]) & (file4_det7_sf['SENS_WAVE'] <= areas_for_fitting[1][1])

    # The polyfit region is an approximation, stitch_polyfit_sf finds the closest match
    # between the polynomial fit and the two sensfuncs within this range
    polyfit_region = (8800, 9300)


    # We'll let the polynomial fit fill in the gap in the wavelengths using the telgrid file
    combined_sf5, fit_info = stitch_sf_polyfit(combined_sf4,     translate_sf(file4_det7_sf, offset),
                                               cf4_fitting_mask, f4d7_fitting_mask,
                                               5, # order for polynomial fit
                                               polyfit_region,
                                               gen_wave_grid=True, 
                                               telgrid_file = "TelFit_MaunaKea_3100_26100_R20000.fits")

    msgs.info(f"Region overwritten with polyfit: {fit_info[0]} to {fit_info[1]}")

    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    polyfit_areas = (combined_sf5['SENS_WAVE'] > fit_info[0]) & (combined_sf5['SENS_WAVE'] < fit_info[1])

    return (combined_sf5, polyfit_areas, fit_info)



def stitch_1200B_sensfunc(sflist):
    """
    Stitch together sensfunc files generated from throughput data from Gret Wriths scripts. Specifically
    The following files were used:
        sens_from_gdw_data/extract/1200B/sens_2017oct03_d1003_0055.fits
        sens_from_gdw_data/extract/1200B/sens_2017oct03_d1003_0056.fits
        sens_from_gdw_data/extract/1200B/sens_2017oct03_d1003_0057.fits

    """

    # Sanity check both sensors from the first sensfunc file

    file1_det3_sf = sanity_check_sf(sflist[0].sens[0])
    file1_det7_sf = sanity_check_sf(sflist[0].sens[1])

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
    file1_det3_fitting_mask = file1_det3_sf['SENS_WAVE'] < polyfit1_det_min
    file1_det7_fitting_mask = file1_det7_sf['SENS_WAVE'] > polyfit1_det_max

    fitting_wave   = np.concatenate((file1_det3_sf['SENS_WAVE'][file1_det3_fitting_mask],
                                     file1_det7_sf['SENS_WAVE'][file1_det7_fitting_mask]))
    fitting_zp_fit = np.concatenate((file1_det3_sf['SENS_ZEROPOINT_FIT'][file1_det3_fitting_mask],
                                     file1_det7_sf['SENS_ZEROPOINT_FIT'][file1_det7_fitting_mask]))

    (combined_sf, fit_info) = stitch_sf_polyfit_old(file1_det3_sf,
                                                    file1_det7_sf,
                                                    fitting_wave,fitting_zp_fit,
                                                    15, 
                                                    use_polyfit_min=polyfit1_write_min, 
                                                    use_polyfit_max=polyfit1_write_max)


    # Sanity check both sensors from the second sensfunc file
    file2_det3_sf = sanity_check_sf(sflist[1].sens[0])
    file2_det7_sf = sanity_check_sf(sflist[1].sens[1])

    # For detector 3 in the second, take everything that doesn't overlap detector 7 in the first
    # Move it up a bit and it's a smooth enought fit that no polynomial fitting is needed
    max_file1_det7_wave = np.max(file1_det7_sf['SENS_WAVE'])
    trans_file2_det3 = translate_sf(truncate_sf(file2_det3_sf, min_wave=max_file1_det7_wave), 0.0474)
    combined_sf = stitch_sf(combined_sf, trans_file2_det3)
    
    # Take all of detecter 7 in the second file. 
    # We nudge it down a bit match det 3 better and use a polynomial fit
    # to smooth the disonctinuity 

    # The discontinuity around the det 3/7 edge in sensfunc file 2
    polyfit2_det_min = 6850.
    polyfit2_det_max = 7200.

    # The region to overwrite with the polyfit values
    polyfit2_write_min = 6747.3
    polyfit2_write_max = 7289.15

    file2_det7_offset = -0.0845
    trans_file2_det7 = translate_sf(file2_det7_sf, file2_det7_offset)

    file2_det3_fitting_mask = trans_file2_det3['SENS_WAVE'] < polyfit2_det_min
    file2_det7_fitting_mask = trans_file2_det7['SENS_WAVE'] > polyfit2_det_max
    fitting_wave   = np.concatenate((trans_file2_det3['SENS_WAVE'][file2_det3_fitting_mask],
                                     trans_file2_det7['SENS_WAVE'][file2_det7_fitting_mask]))
    fitting_zp_fit = np.concatenate((trans_file2_det3['SENS_ZEROPOINT_FIT'][file2_det3_fitting_mask],
                                     trans_file2_det7['SENS_ZEROPOINT_FIT'][file2_det7_fitting_mask]))

    (combined_sf, fit_info) = stitch_sf_polyfit_old(
                                        combined_sf,
                                        trans_file2_det7,
                                        fitting_wave, fitting_zp_fit,
                                        15, 
                                        use_polyfit_min=polyfit2_write_min,
                                        use_polyfit_max=polyfit2_write_max)

    

    # For the third sensfunc file, ignore detector 3 because it overlaps everything.

    # Sanity check detector 7 from the third file
    file3_det7_sf = sanity_check_sf(sflist[2].sens[1])

    # Take the non overlapping parts of detector 7 and do a polynomail fit to smooth over the join
    trunc_file3_det7 = truncate_sf(file3_det7_sf, min_wave = np.max(combined_sf['SENS_WAVE']))


    # The discontinuity around the joining of file2 det 7 and file 3 det 7
    polyfit3_det_min = 8250
    polyfit3_det_max = 8450

    # The region to overwrite with the polyfit values
    polyfit3_write_min = 8211.1
    polyfit3_write_max = 8544.1

    file2_det7_fitting_mask = trans_file2_det7['SENS_WAVE'] < polyfit3_det_min
    file3_det7_fitting_mask = trunc_file3_det7['SENS_WAVE'] > polyfit3_det_max
    fitting_wave   = np.concatenate((trans_file2_det7['SENS_WAVE'][file2_det7_fitting_mask],
                                     trunc_file3_det7['SENS_WAVE'][file3_det7_fitting_mask]))
    fitting_zp_fit = np.concatenate((trans_file2_det7['SENS_ZEROPOINT_FIT'][file2_det7_fitting_mask],
                                     trunc_file3_det7['SENS_ZEROPOINT_FIT'][file3_det7_fitting_mask]))


    (combined_sf, fit_info) = stitch_sf_polyfit_old(
                                        combined_sf,
                                        trunc_file3_det7,
                                        fitting_wave, fitting_zp_fit,
                                        15, 
                                        use_polyfit_min=polyfit3_write_min,
                                        use_polyfit_max=polyfit3_write_max)


    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    combined_edges = (   ((combined_sf['SENS_WAVE'] >= polyfit1_write_min) & (combined_sf['SENS_WAVE'] <= polyfit1_write_max))
                       | ((combined_sf['SENS_WAVE'] >= polyfit2_write_min) & (combined_sf['SENS_WAVE'] <= polyfit2_write_max))
                       | ((combined_sf['SENS_WAVE'] >= polyfit3_write_min) & (combined_sf['SENS_WAVE'] <= polyfit3_write_max)) )
    return (combined_sf, combined_edges, fit_info)

def stitch_900ZD_sensfunc(sflist):

    """
    Stitch together sensfunc files generated from throughput data from Gret Wriths scripts. Specifically
    The following files were used:
        sens_from_gdw_data/extract/900ZD/sens_2010oct06_d1006_0138.fits
        sens_from_gdw_data/extract/900ZD/sens_2010oct06_d1006_0139.fits
        sens_from_gdw_data/extract/900ZD/sens_2010oct06_d1006_0140.fits
    """

    # Sanity check the first sensfunc file by masking any 0 values and ensuring there are no NaNs or +/- inf values
    file1_det3_sf = sanity_check_sf(sflist[0].sens[0])
    file1_det7_sf = sanity_check_sf(sflist[0].sens[1])

    # Combine file 1 detectors 3 and 7 with a polynomial fit to smooth out the discontinuity.

    # Mask the fitting around the detector boundary and towards the outer ends of the sensfuncs to get
    # the best fit
    file1_det3_fitting_mask = (file1_det3_sf['SENS_WAVE'] > 4629.4) & (file1_det3_sf['SENS_WAVE'] < 5324.5)
    file1_det7_fitting_mask = (file1_det7_sf['SENS_WAVE'] > 5530.4) & (file1_det7_sf['SENS_WAVE'] < 9000)

    fitting_wave   = np.concatenate((file1_det3_sf['SENS_WAVE'][file1_det3_fitting_mask],
                                     file1_det7_sf['SENS_WAVE'][file1_det7_fitting_mask]))
    fitting_zp_fit = np.concatenate((file1_det3_sf['SENS_ZEROPOINT_FIT'][file1_det3_fitting_mask],
                                     file1_det7_sf['SENS_ZEROPOINT_FIT'][file1_det7_fitting_mask]))

    # The region to overwrite with the polyfit values
    polyfit1_write_min = 4910.2
    polyfit1_write_max = 5805.2

    (combined_sf, fit_info) = stitch_sf_polyfit_old(
                                        file1_det3_sf, file1_det7_sf,
                                        fitting_wave, fitting_zp_fit,
                                        15, 
                                        use_polyfit_min=polyfit1_write_min, use_polyfit_max=polyfit1_write_max)

    
    
    # Now add the parts of file 3, det3 that don't overlap file 1 det 7.
    # We're skipping the second file entirely because of overlaps
    file3_det3_sf = sanity_check_sf(sflist[2].sens[0])

    trunc_file3_det3 = truncate_sf(file3_det3_sf, min_wave=np.max(file1_det7_sf['SENS_WAVE']))

    # We also translate up file 3 det 3 up to match file 1 det 7
    file3_det3_offset = 0.10533590352139655
    trans_file3_det3 = translate_sf(trunc_file3_det3, file3_det3_offset)

    combined_sf = stitch_sf(combined_sf, trans_file3_det3)

    # Now add file 3 det 7, moving it down to match det 3 and using a polynomial fit to smooth the 
    # discontinuities between det 3 and 7

    file3_det7_sf = sanity_check_sf(sflist[2].sens[1])
    file3_det7_offset = -0.02687150393769855
    trans_file3_det7 = translate_sf(file3_det7_sf, file3_det7_offset)

    file3_det3_fitting_mask = trans_file3_det3['SENS_WAVE'] > 7600
    file3_det7_fititng_mask = trans_file3_det7['SENS_WAVE'] < 8400

    fitting_wave   = np.concatenate((trans_file3_det3['SENS_WAVE'][file3_det3_fitting_mask],
                                     trans_file3_det7['SENS_WAVE'][file3_det7_fititng_mask]))
    fitting_zp_fit = np.concatenate((trans_file3_det3['SENS_ZEROPOINT_FIT'][file3_det3_fitting_mask],
                                     trans_file3_det7['SENS_ZEROPOINT_FIT'][file3_det7_fititng_mask]))
    
    polyfit2_write_min = 7876.
    polyfit2_write_max = 8055.7

    (combined_sf, fit_info) = stitch_sf_polyfit_old(
                                        combined_sf, trans_file3_det7,
                                        fitting_wave, fitting_zp_fit,
                                        30,
                                        use_polyfit_min=polyfit2_write_min, use_polyfit_max=polyfit2_write_max)

    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    polyfit_areas = (  ((combined_sf['SENS_WAVE'] > polyfit1_write_min) & (combined_sf['SENS_WAVE'] < polyfit1_write_max)) 
                     | ((combined_sf['SENS_WAVE'] > polyfit2_write_min) & (combined_sf['SENS_WAVE'] < polyfit2_write_max)))


    return (combined_sf, polyfit_areas, fit_info)

def stitch_600ZD_sensfunc(sflist):
    """
    Stitch together sensfunc files generated from throughput data from Gret Wriths scripts. Specifically
    The following files were used:
        sens_from_gdw_data/extract/600ZD/sens_2023jan17_d0117_0098.fits
        sens_from_gdw_data/extract/600ZD/sens_2023jan17_d0117_0099.fits
        sens_from_gdw_data/extract/600ZD/sens_2023jan17_d0117_0101.fits

    """


    # Sanity check the sensfuncs
    file1_det3_sf = sanity_check_sf(sflist[0].sens[0])
    file1_det7_sf = sanity_check_sf(sflist[0].sens[1])
    file2_det7_sf = sanity_check_sf(sflist[1].sens[1])

    # We call this "file4" because there were originally four exposures taken, but we aren't using the
    # third
    file4_det7_sf = sanity_check_sf(sflist[2].sens[1])

    # Stitch file1 det3 and det7 with a polynomial fit to smooth the discontinuity between them

    # The polynomial fit will be done on a subset of the sensfuncs to improve the fit and eliminate
    # odd behavior near the detector boundary. These numbers were chosen by examining plots of the
    # individual sensfuncs
    areas_for_fitting = [(5500, 5850), (6100, 6400)]
    f1d3_fitting_mask = (file1_det3_sf['SENS_WAVE'] >= areas_for_fitting[0][0]) & (file1_det3_sf['SENS_WAVE'] <= areas_for_fitting[0][1])
    f1d7_fitting_mask = (file1_det7_sf['SENS_WAVE'] >= areas_for_fitting[1][0]) & (file1_det7_sf['SENS_WAVE'] <= areas_for_fitting[1][1])

    # An approximate region to overwrite with the polyfit values.
    # stitch_sf_polyfit will figure out the exact best place to use within this region
    approx_polyfit_region = [5800, 6200]

    # Experiment showed that order 5 produces a good fit
    order = 5

    combined_sf1, fit_info = stitch_sf_polyfit(file1_det3_sf, file1_det7_sf,
                                               f1d3_fitting_mask, f1d7_fitting_mask,
                                               order, approx_polyfit_region,
                                               gen_wave_grid=True, telgrid_file = "TelFit_MaunaKea_3100_26100_R20000.fits")


    # Now add the sensfunc from file 2 det 7.

    # The stitching point was chosen by looking at the plot of file1 det7's sensfunc and picking a
    # wavelength before it curved upwards at the detector boundary
    stitch_point = 7600

    # Truncate the combined sensfunc at that point
    combined_sf1_trunc = truncate_sf(combined_sf1, max_wave=stitch_point)

    # File 2 det 7's sensfunc needs to be translated up a bit to be a better match
    # Find the offset by subtracting the two functions at the stitching point
    combined_sf1_last_zp = combined_sf1_trunc['SENS_ZEROPOINT_FIT'][-1]
    combined_sf1_last_wave = combined_sf1_trunc['SENS_WAVE'][-1]
    file2_det7_mask = file2_det7_sf['SENS_WAVE'] >= combined_sf1_last_wave
    file2_det7_first_zp = file2_det7_sf['SENS_ZEROPOINT_FIT'][file2_det7_mask][0]
    offset = combined_sf1_last_zp - file2_det7_first_zp

    # Move file2_det7 up a bit
    trans_file2_det7 = translate_sf(file2_det7_sf, offset)

    # Stitch
    combined_sf2 = stitch_sf(combined_sf1_trunc, trans_file2_det7)

    # Choose file 4 detector 7 for last sensfunc to join. We translate it up to meet the combined sf. 
    # This time we subtract the values of the two sensitivity functions at the stitch point 
    # to find the offset instead of finidng one by hand.

    combined_sf2_last_zp = combined_sf2['SENS_ZEROPOINT_FIT'][-1]
    combined_sf2_last_wave = combined_sf2['SENS_WAVE'][-1]
    file4_det7_mask = file4_det7_sf['SENS_WAVE'] >= combined_sf2_last_wave
    file4_det7_first_zp = file4_det7_sf['SENS_ZEROPOINT_FIT'][file4_det7_mask][0]
    offset = combined_sf2_last_zp - file4_det7_first_zp

    trans_file4_det7 = translate_sf(file4_det7_sf, offset)
    combined_sf3 = stitch_sf(combined_sf2, trans_file4_det7)

    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    polyfit_areas = (combined_sf3['SENS_WAVE'] > fit_info[0]) & (combined_sf3['SENS_WAVE'] < fit_info[1])
    
    return combined_sf3, polyfit_areas, fit_info


def stitch_830G_sensfunc(sflist):

    file1_det3_sf = sanity_check_sf(sflist[0].sens[0])
    file2_det3_sf = sanity_check_sf(sflist[1].sens[0])
    file3_det3_sf = sanity_check_sf(sflist[2].sens[0])
    file2_det7_sf = sanity_check_sf(sflist[1].sens[1])
    file3_det7_sf = sanity_check_sf(sflist[2].sens[1])
    file4_det7_sf = sanity_check_sf(sflist[3].sens[1])

    ### Take the sensfuncs from first detector in file 1 and in file 2, merging them at the point with the closest slope

    area_to_stitch = (5100, 5400)
    combined_sf1, stitch_point, offset = gradient_stitch_sf(file1_det3_sf, file2_det3_sf, area_to_stitch)

    ### Now translate up the sensfunc from det 3 in file 3 to match the joined sensfunc and 
    ### stitch it on at around 6500
    area_to_stitch = (6300, 6700)

    # Compare gradients to find the point where the two functions are closest in slope
    combined_sf2, stitch_point, offset = gradient_stitch_sf(combined_sf1, file3_det3_sf, area_to_stitch)

    ### Now translate sensfunc 7 from file 2 down to match the end of the combined sensfunc and stitch them 
    ### together at around 7750
    area_to_stitch = (7700, 7900)
    combined_sf3, stitch_point, offset = gradient_stitch_sf(combined_sf2, file2_det7_sf, area_to_stitch)

    area_to_stitch = (8300, 8600)
    combined_sf4, stitch_point, offset = gradient_stitch_sf(combined_sf3, file3_det7_sf, area_to_stitch)

    area_to_stitch=(9800, 9925)
    combined_sf5, stitch_point, offset = gradient_stitch_sf(combined_sf4, file4_det7_sf, area_to_stitch)


    return (combined_sf5, None, None)