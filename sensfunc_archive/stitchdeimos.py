import numpy as np
from pypeit.core.wavecal import wvutils
from stitchutils import sanity_check_sf, stitch_sf_polyfit, stitch_sf, truncate_sf, translate_sf

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

    # The below area is overwritten by the polyfit results
    area_to_use_polyfit = (8817, 9275)


    # We'll let the polynomial fit fill in the gap in the wavelengths using the telgrid file
    combined_sf5, fit_info = stitch_sf_polyfit(combined_sf4,     translate_sf(file4_det7_sf, offset),
                                               cf4_fitting_mask, f4d7_fitting_mask,
                                               5, # order for polynomial fit
                                               use_polyfit_min = area_to_use_polyfit[0],
                                               use_polyfit_max = area_to_use_polyfit[1],
                                               gen_wave_grid=True, 
                                               telgrid_file = "TelFit_MaunaKea_3100_26100_R20000.fits")

    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    polyfit_areas = (combined_sf5['SENS_WAVE'] > area_to_use_polyfit[0]) & (combined_sf5['SENS_WAVE'] < area_to_use_polyfit[1])

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


    (combined_sf, fit_info) = stitch_sf_polyfit(file1_det3_sf,
                                                file1_det7_sf,
                                                file1_det3_fitting_mask,
                                                file1_det7_fitting_mask,
                                                15, 
                                                use_polyfit_min=polyfit1_write_min, 
                                                use_polyfit_max=polyfit1_write_max)


    # Sanity check both sensors from the second sensfunc file
    file2_det3_sf = sanity_check_sf(sflist[1].sens[0])
    file2_det7_sf = sanity_check_sf(sflist[1].sens[1])

    # For detector 3 in the second, take everything that doesn't overlap detector 7 in the first
    # Move it up a bit and it's a smooth enought fit that no polynomial fitting is needed
    combined_sf = stitch_sf(combined_sf, file2_det3_sf, translate2=0.0474)
    
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

    prev_combined_fitting_mask = combined_sf['SENS_WAVE'] < polyfit2_det_min
    file2_det7_fitting_mask = file2_det7_sf['SENS_WAVE'] > polyfit2_det_max

    (combined_sf, fit_info) = stitch_sf_polyfit(
                                        combined_sf,
                                        file2_det7_sf,
                                        prev_combined_fitting_mask,
                                        file2_det7_fitting_mask,
                                        15, 
                                        use_polyfit_min=polyfit2_write_min,
                                        use_polyfit_max=polyfit2_write_max,
                                        translate2=file2_det7_offset)

    

    # For the third sensfunc file, ignore detector 3 because it overlaps everything.

    # Sanity check detector 7 from the third file
    file3_det7_sf = sanity_check_sf(sflist[2].sens[1])

    # Take the non overlapping parts of detector 7 and do a polynomail fit to smooth over the join
    non_overlap = file3_det7_sf['SENS_WAVE'] > np.max(combined_sf['SENS_WAVE'])


    # The discontinuity around the joining of file2 det 7 and file 3 det 7
    polyfit3_det_min = 8250
    polyfit3_det_max = 8450

    # The region to overwrite with the polyfit values
    polyfit3_write_min = 8211.1
    polyfit3_write_max = 8544.1

    prev_combined_fitting_mask = combined_sf['SENS_WAVE'] < polyfit3_det_min
    file3_det7_fitting_mask = non_overlap & (file3_det7_sf['SENS_WAVE'] > polyfit3_det_max)

    (combined_sf, fit_info) = stitch_sf_polyfit(
                                        combined_sf,
                                        file3_det7_sf,
                                        prev_combined_fitting_mask,
                                        file3_det7_fitting_mask,
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


    # The region to overwrite with the polyfit values
    polyfit1_write_min = 4910.2
    polyfit1_write_max = 5805.2

    (combined_sf, fit_info) = stitch_sf_polyfit(
                                        file1_det3_sf, file1_det7_sf,
                                        file1_det3_fitting_mask, file1_det7_fitting_mask,
                                        15, 
                                        use_polyfit_min=polyfit1_write_min, use_polyfit_max=polyfit1_write_max)

    
    
    # Now add the parts of file 3, det3 that don't overlap file 1 det 7.
    # We're skipping the second file entirely because of overlaps
    # We also translate up file 3 det 3 up to match file 1 det 7
    file3_det3_sf = sanity_check_sf(sflist[2].sens[0])
    file3_det3_offset = 0.10533590352139655

    combined_sf = stitch_sf(combined_sf, file3_det3_sf, translate2=0.10533590352139655)

    # Now add file 3 det 7, moving it down to match det 3 and using a polynomial fit to smooth the 
    # discontinuities between det 3 and 7

    file3_det7_sf = sanity_check_sf(sflist[2].sens[1])
    file3_det7_offset = -0.02687150393769855

    prev_combined_fitting_mask = combined_sf['SENS_WAVE'] > 7600
    file3_det7_fitting_mask = file3_det7_sf['SENS_WAVE'] < 8400
    
    polyfit2_write_min = 7876.
    polyfit2_write_max = 8055.7

    (combined_sf, fit_info) = stitch_sf_polyfit(
                                        combined_sf, file3_det7_sf,
                                        prev_combined_fitting_mask, file3_det7_fitting_mask,
                                        30, translate2=file3_det7_offset,
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

    # Stitch file1 det3 and det7 with a polynomial fit to smooth the discontinuity between them

    # Eliminate any 0 values in SENS_WAVE or SENS_ZEROPOINT_FIT
    file1_det3_sf = sanity_check_sf(sflist[0].sens[0])
    file1_det7_sf = sanity_check_sf(sflist[0].sens[1])

    # The polynomial fit will be done on a subset of the sensfuncs to improve the fit and eliminate
    # odd behavior near the detector boundary. These numbers were chosen by examining plots of the
    # individual sensfuncs
    file1_det3_fitting_mask = (file1_det3_sf['SENS_WAVE'] > 5500.) & (file1_det3_sf['SENS_WAVE'] < 5850.)
    file1_det7_fitting_mask = (file1_det7_sf['SENS_WAVE'] > 6100.) & (file1_det7_sf['SENS_WAVE'] < 6400.)

    # The region to overwrite with the polyfit values
    polyfit1_write_min = 5818.8
    polyfit1_write_max = 6137.2

    (combined_sf, fit_info) = stitch_sf_polyfit(file1_det3_sf, file1_det7_sf,
                                                file1_det3_fitting_mask, file1_det7_fitting_mask,
                                                5,
                                                use_polyfit_min=polyfit1_write_min, use_polyfit_max=polyfit1_write_max)


    # Append the sensfunc from file 2 det 7.
    # The two sensfuncs are stitched together at a point that avoids odd behavior at the detector boundary
    # in file 1 det 7's sensfunc.  File 2 det 7's sensfunc is then translated up a bit to be a better match
    file2_det7_sf = sanity_check_sf(sflist[1].sens[1])
    combined_sf = stitch_sf(combined_sf, file2_det7_sf, stitch_point = 7599.5, translate2 = 0.01047)


    # Append the portion of file3 detector 7 that does not overlap file 2 detector 7. It is also 
    # translated up to match
    file3_det7_sf = sanity_check_sf(sflist[2].sens[1])

    combined_sf = stitch_sf(combined_sf, file3_det7_sf, translate2 = 0.0568)

    # Create a mask for the areas that were filled in with a polynomial fit. This is used for 
    # plotting
    polyfit_areas = (combined_sf['SENS_WAVE'] > polyfit1_write_min) & (combined_sf['SENS_WAVE'] < polyfit1_write_max)
    
    return combined_sf, polyfit_areas, fit_info


def stitch_830G_sensfunc_v4(sflist): # v4

    # Sanity check our input
    file1_sf = sanity_check_sf(sflist[0].sens[0])
    file2_sf = sanity_check_sf(sflist[1].sens[0])
    file3_sf = sanity_check_sf(sflist[2].sens[0])
    file4_sf = sanity_check_sf(sflist[3].sens[0])

    # Start by stitching file 1 with file 2's sensfunc.
    # We do this by picking the area where we want to stitch visually.
    # In this case it's the region between the start of file 2's sensfunc and 5600 angstroms.
    area_to_stitch = (file2_sf['SENS_WAVE'][0], 5600)

    # Then we find the point within that range that the two sensfuncs are closest in gradient. 
    # This should have the smoothest connection between the two sensfuncs.
    # We compute this point rather than hardcoding it can change slightly each time they're generated, 
    # but should remain consistent enough that this will work.
    diff_mask1 = (file1_sf['SENS_WAVE'] > area_to_stitch[0]) & (file1_sf['SENS_WAVE'] < area_to_stitch[1])
    diff_mask2 = (file2_sf['SENS_WAVE'] > area_to_stitch[0]) & (file2_sf['SENS_WAVE'] < area_to_stitch[1])
    g1 = np.gradient(file1_sf['SENS_ZEROPOINT_FIT'], file1_sf['SENS_WAVE'])
    g2 = np.gradient(file2_sf['SENS_ZEROPOINT_FIT'], file2_sf['SENS_WAVE'])

    min_index = np.argmin(np.fabs(g1[diff_mask1] - g2[diff_mask2]))
    stitch_point = file1_sf['SENS_WAVE'][diff_mask1][min_index]
    
    # We then translate the second sensfunc down slightly to match the first at this point, 
    # and then combine them.
    offset = file1_sf['SENS_ZEROPOINT_FIT'][diff_mask1][min_index] - file2_sf['SENS_ZEROPOINT_FIT'][diff_mask2][min_index]
    msgs.info(f"Stitching file 1 and file2 at {stitch_point} with offset {offset}")
    combined_sf1 = stitch_sf(file1_sf, translate_sf(file2_sf, offset), stitch_point=stitch_point)

    # We do the same to add file 3 to the combined sensfunc
    area_to_stitch = (6150, 7000)

    diff_mask1 = (combined_sf1['SENS_WAVE'] > area_to_stitch[0]) & (combined_sf1['SENS_WAVE'] < area_to_stitch[1])
    diff_mask2 = (file3_sf['SENS_WAVE'] > area_to_stitch[0]) & (file3_sf['SENS_WAVE'] < area_to_stitch[1])
    g1 = np.gradient(combined_sf1['SENS_ZEROPOINT_FIT'], combined_sf1['SENS_WAVE'])
    g2 = np.gradient(file3_sf['SENS_ZEROPOINT_FIT'], file3_sf['SENS_WAVE'])

    min_index = np.argmin(np.fabs(g1[diff_mask1] - g2[diff_mask2]))
    stitch_point = combined_sf1['SENS_WAVE'][diff_mask1][min_index]
    
    offset = combined_sf1['SENS_ZEROPOINT_FIT'][diff_mask1][min_index] - file3_sf['SENS_ZEROPOINT_FIT'][diff_mask2][min_index]
    msgs.info(f"Stitching file3 at {stitch_point} with offset {offset}")
    combined_sf2 = stitch_sf(combined_sf1, translate_sf(file3_sf, offset), stitch_point=stitch_point)

    # And again for file 4
    area_to_stitch=(9400,9900)

    diff_mask1 = (combined_sf2['SENS_WAVE'] > area_to_stitch[0]) & (combined_sf2['SENS_WAVE'] < area_to_stitch[1])
    diff_mask2 = (file4_sf['SENS_WAVE'] > area_to_stitch[0]) & (file4_sf['SENS_WAVE'] < area_to_stitch[1])
    g1 = np.gradient(combined_sf2['SENS_ZEROPOINT_FIT'], combined_sf2['SENS_WAVE'])
    g2 = np.gradient(file4_sf['SENS_ZEROPOINT_FIT'], file4_sf['SENS_WAVE'])

    min_index = np.argmin(np.fabs(g1[diff_mask1] - g2[diff_mask2]))
    stitch_point = combined_sf2['SENS_WAVE'][diff_mask1][min_index]
    
    offset = combined_sf2['SENS_ZEROPOINT_FIT'][diff_mask1][min_index] - file4_sf['SENS_ZEROPOINT_FIT'][diff_mask2][min_index]
    msgs.info(f"Stitching file4 at {stitch_point} with offset {offset}")
    combined_sf3= stitch_sf(combined_sf2, translate_sf(file4_sf, offset), stitch_point=stitch_point)

    return combined_sf3, None, None

def stitch_830G_sensfunc_v3(sflist): # v3

    # Sanity check our input
    file1_sf = sanity_check_sf(sflist[0].sens[0])
    file2_sf = sanity_check_sf(sflist[1].sens[0])
    file3_sf = sanity_check_sf(sflist[2].sens[0])
    file4_sf = sanity_check_sf(sflist[3].sens[0])

    #oops deleted, need to re-create from jupyter

def stitch_830G_sensfunc_v2(sflist):

    # Sanity check our input
    file1_sf = sanity_check_sf(sflist[0].sens[0])
    file3_sf = sanity_check_sf(sflist[2].sens[0])
    file4_sf = sanity_check_sf(sflist[3].sens[0])

    # Start by stitching file 1 with file 2's sensfunc.

    # We stitch starting at the largest value in file 1's sensfunc
    file1_max_idx = np.argmax(file1_sf['SENS_ZEROPOINT_FIT'])
    file1_max_zp = file1_sf['SENS_ZEROPOINT_FIT'][file1_max_idx]
    wave_at_max_file1_zp = file1_sf['SENS_WAVE'][file1_max_idx]

    # We do this with a polynomial fit to smooth the transition between the two sensfuncs.
    # We know from visual inspection/trial and error that we only want to polyfit a specific range, defined
    # by fitting masks.
    file1_fitting_mask = (file1_sf['SENS_WAVE'] >= 5500) & (file1_sf['SENS_WAVE'] <=  wave_at_max_file1_zp)
    file3_fitting_mask = (file3_sf['SENS_ZEROPOINT_FIT'] >= file1_max_zp) & (file3_sf['SENS_WAVE'] <= 8000)

    # First we get the entire polynomial fit, and then find where to join it to the existing sensfuncs
    polyfit_sf, fi = stitch_sf_polyfit(file1_sf, file3_sf, file1_fitting_mask, file3_fitting_mask, 8)         

    # We join the polynomial fit to the two sensitivity functions by finding the closest point within a known
    # wavelength range found by looking at plots of all functions. We do it this way because sensitivity functions can change slightly 
    # each time they're generated.

    # Find the closest point between file 1 and the polynomial fit
    diff_mask1 = (file1_sf['SENS_WAVE']>6220) & (file1_sf['SENS_WAVE']<6240)
    diff_mask2 = (polyfit_sf['SENS_WAVE']>6220) & (polyfit_sf['SENS_WAVE']<6240)

    # Find the closest point
    min_index = np.argmin(np.fabs(file1_sf['SENS_ZEROPOINT_FIT'][diff_mask1] - polyfit_sf['SENS_ZEROPOINT_FIT'][diff_mask2]))

    use_polyfit_min1 = file1_sf['SENS_WAVE'][diff_mask1][min_index]

    # Do the same for file 3
    diff_mask1 = (file3_sf['SENS_WAVE']>7110) & (file3_sf['SENS_WAVE']<7130)
    diff_mask2 = (polyfit_sf['SENS_WAVE']>7110) & (polyfit_sf['SENS_WAVE']<7130)

    # Find the closest point
    min_index = np.argmin(np.fabs(file3_sf['SENS_ZEROPOINT_FIT'][diff_mask1] - polyfit_sf['SENS_ZEROPOINT_FIT'][diff_mask2]))

    use_polyfit_max1 = file3_sf['SENS_WAVE'][diff_mask1][min_index]

    # Now stitch using this min/max
    combined_sf1, fi = stitch_sf_polyfit(file1_sf, file3_sf, file1_fitting_mask, file3_fitting_mask, 8, gen_wave_grid=False, 
                                         use_polyfit_min=use_polyfit_min1, use_polyfit_max=use_polyfit_max1)         

    # We now do the same process to stitch the combined sensfunc to file4
    pre_comb_fitting_mask = (combined_sf1['SENS_WAVE'] > 8800) & (combined_sf1['SENS_WAVE'] < 9200)
    file4_fitting_mask = (file4_sf['SENS_WAVE']>9450) & (file4_sf['SENS_WAVE'] < 9700 )

    # Get a temporary polynomial fit
    polyfit_sf, fi = stitch_sf_polyfit(combined_sf1, file4_sf, pre_comb_fitting_mask, file4_fitting_mask, 5)


    # Find the closest point between the fit and the previous combined sensfunc
    diff_mask1 = (combined_sf1['SENS_WAVE']>9100) & (combined_sf1['SENS_WAVE']<9200)
    diff_mask2 = (polyfit_sf['SENS_WAVE']>9100) & (polyfit_sf['SENS_WAVE']<9200)

    min_index = np.argmin(np.fabs(combined_sf1['SENS_ZEROPOINT_FIT'][diff_mask1] - polyfit_sf['SENS_ZEROPOINT_FIT'][diff_mask2]))

    use_polyfit_min2 = combined_sf1['SENS_WAVE'][diff_mask1][min_index]

    # Find the closest point between the fit and file 4's sensfunc
    diff_mask1 = (file4_sf['SENS_WAVE']>9560) & (file4_sf['SENS_WAVE']<9600)
    diff_mask2 = (polyfit_sf['SENS_WAVE']>9560) & (polyfit_sf['SENS_WAVE']<9600)
    min_index = np.argmin(np.fabs(file4_sf['SENS_ZEROPOINT_FIT'][diff_mask1] - polyfit_sf['SENS_ZEROPOINT_FIT'][diff_mask2]))

    use_polyfit_max2 = file4_sf['SENS_WAVE'][diff_mask1][min_index]

    # Use this range to stitch the sensfuncs. The stitch point of 9400 is chosen to be between use_polyfit_min and use_polyfit_max
    combined_sf2, fi =stitch_sf_polyfit(combined_sf1, file4_sf,
                                        pre_comb_fitting_mask, file4_fitting_mask,
                                        5, stitch_point=9400, gen_wave_grid=False,
                                        use_polyfit_min=use_polyfit_min2, use_polyfit_max=use_polyfit_max2)


    polyfit_areas = (((combined_sf2['SENS_WAVE'] >= use_polyfit_min1) & (combined_sf2['SENS_WAVE'] <= use_polyfit_max1)) | 
                     ((combined_sf2['SENS_WAVE'] >= use_polyfit_min2) & (combined_sf2['SENS_WAVE'] <= use_polyfit_max2)))
    return combined_sf2, polyfit_areas, None

def stitch_830G_sensfunc_v1(sflist):

    # Sanity check our input
    file1_sf = sanity_check_sf(sflist[0].sens[0])
    file2_sf = sanity_check_sf(sflist[1].sens[0])
    file4_sf = sanity_check_sf(sflist[3].sens[0])

    # Start by stitching file 1 with file 2's sensfunc.
    # We do this by picking the area where we want to stitch visually.
    # In this case it's the region between 5200 and 5600 angstroms.
    # Then we find the cloest point between the two sensfuncs along this range.
    # This will be the stitching point. We do it this way because sensfuncs
    # can change slightly each time they're generated, but this should remain consistent. 
    # We then translate the second sensfunc down slightly to match the first at this point, 
    # and then combine them.

    # Find a range within which to subtract the two sensfuncs
    diff_mask1 = (file1_sf['SENS_WAVE'] > 5200) & (file1_sf['SENS_WAVE'] < 5600)
    diff_mask2 = (file2_sf['SENS_WAVE'] > 5200) & (file2_sf['SENS_WAVE'] < 5600)

    # Make sure the start of each range is equal
    assert file1_sf['SENS_WAVE'][diff_mask1][0] == file2_sf['SENS_WAVE'][diff_mask2][0]

    # Make sure the ranges are equal in length
    assert len(file1_sf['SENS_WAVE'][diff_mask1]) == len(file2_sf['SENS_WAVE'][diff_mask2]), "Need to pick equal length subsets of sensfuncs to subtract"

    # Find the closest point
    min_index = np.argmin(file1_sf['SENS_ZEROPOINT_FIT'][diff_mask1] - file2_sf['SENS_ZEROPOINT_FIT'][diff_mask2])

    # Find the offset needed to get file2's sensfunc to
    # match file1's at this point
    offset = file1_sf['SENS_ZEROPOINT_FIT'][diff_mask1][min_index] - file2_sf['SENS_ZEROPOINT_FIT'][diff_mask2][min_index]

    stitch_point = file1_sf['SENS_WAVE'][diff_mask1][min_index]
    msgs.info(f"830G stitch point 1 is at {stitch_point}")
    msgs.info(f"830G offset for sensfunc file 2 is {offset}")

    translated_file2_sf = translate_sf(file2_sf, offset)
    combined_sf1 = stitch_sf(file1_sf, translated_file2_sf, stitch_point=stitch_point)

    # Next we add file4's sensfunc to this, using the same technique as above
    diff_mask1 = (combined_sf1['SENS_WAVE'] > 7500) & (combined_sf1['SENS_WAVE'] < 8000)
    diff_mask2 = (file4_sf['SENS_WAVE'] > 7500) & (file4_sf['SENS_WAVE'] < 8000)

    # Make sure the start of each range is the same
    assert combined_sf1['SENS_WAVE'][diff_mask1][0] == file4_sf['SENS_WAVE'][diff_mask2][0]

    # Make sure the ranges are equal in length
    assert len(combined_sf1['SENS_WAVE'][diff_mask1]) == len(file4_sf['SENS_WAVE'][diff_mask2]), "Need to pick equal length subsets of sensfuncs to subtract"

    # Find the closest point
    min_index = np.argmin(combined_sf1['SENS_ZEROPOINT_FIT'][diff_mask1] - file4_sf['SENS_ZEROPOINT_FIT'][diff_mask2])

    # Use this as the point to stitch the sesnfuncs together, and find the offset needed to get file2's sensfunc to
    # match file1's
    offset = combined_sf1['SENS_ZEROPOINT_FIT'][diff_mask1][min_index] - file4_sf['SENS_ZEROPOINT_FIT'][diff_mask2][min_index]
    stitch_point = combined_sf1['SENS_WAVE'][diff_mask1][min_index]
    msgs.info(f"830G stitch point 2 is {stitch_point}")
    msgs.info(f"830G offset for sensfunc file 4 is {offset}")

    translated_file4_sf = translate_sf(file4_sf, offset)
    combined_sf2 = stitch_sf(combined_sf1, translated_file4_sf, stitch_point=stitch_point)

    return combined_sf2, None, None


def stitch_830G_sensfunc_old(sflist):

    """
    Stitch together sensfunc files generated from throughput data from Gret Wriths scripts. Specifically
    The following files were used:
        sens_from_gdw_data/extract/830G/sens_2010oct06_d1006_0135.fits
        sens_from_gdw_data/extract/830G/sens_2010oct06_d1006_0136.fits
        sens_from_gdw_data/extract/830G/sens_2010oct06_d1006_0137.fits
    """

    # Sanity check the sensfuncs by masking any 0 values and ensuring there are no NaNs or +/- inf values
    file1_det3_sf = sanity_check_sf(sflist[0].sens[0])
    file2_det3_sf = sanity_check_sf(sflist[1].sens[0])

    # Take the sensfuncs from first detector in file 1 and in file 2, merging them at the point where they intersect
    stitch_point = 5225.7

    combined_sf = stitch_sf(file1_det3_sf, file2_det3_sf, stitch_point)

    # Now translate up the sensfunc from det 3 in file 3 to match the joined sensfunc and stitch it on at around 6500

    file3_det3_sf = sanity_check_sf(sflist[2].sens[0])
    stitch_point = 6500
    file3_det3_offset = +0.08742193265786469

    combined_sf = stitch_sf(combined_sf, translate_sf(file3_det3_sf, file3_det3_offset), stitch_point=stitch_point)
    
    # Now translate sensfunc 7 from file 2 down to match the end of the combined sensfunc and stitch them together at around 7750

    file2_det7_sf = sanity_check_sf(sflist[1].sens[1])

    stitch_point = 7750
    file2_det7_offset = -0.03777170821255993

    combined_sf = stitch_sf(combined_sf, translate_sf(file2_det7_sf, file2_det7_offset), stitch_point=stitch_point)

    # Finish the stitched sensfunc by adding the non-overlapping portions of det 7 from file 3.
    file3_det7_sf = sanity_check_sf(sflist[2].sens[1])
    file3_det7_offset = -0.003240994416337628

    combined_sf = stitch_sf(combined_sf, translate_sf(file3_det7_sf, file3_det7_offset))

    return (combined_sf, None, None)

def gradient_stitch(sf1, sf2, area_to_stitch):
    sf1_g = np.gradient(sf1['SENS_ZEROPOINT_FIT'], sf1['SENS_WAVE'])
    sf2_g = np.gradient(sf2['SENS_ZEROPOINT_FIT'], sf2['SENS_WAVE'])

    diff_mask1 = (sf1['SENS_WAVE'] > area_to_stitch[0]) & (sf1['SENS_WAVE'] < area_to_stitch[1])
    diff_mask2 = (sf2['SENS_WAVE'] > area_to_stitch[0]) & (sf2['SENS_WAVE'] < area_to_stitch[1])
    min_index = np.argmin(np.fabs(sf1_g[diff_mask1] - sf2_g[diff_mask2]))
    stitch_point = sf1['SENS_WAVE'][diff_mask1][min_index]
    offset = sf1['SENS_ZEROPOINT_FIT'][diff_mask1][min_index] - sf2['SENS_ZEROPOINT_FIT'][diff_mask2][min_index]
    combined_sf = stitch_sf(sf1, translate_sf(sf2, offset), stitch_point)
    return combined_sf, stitch_point, offset

def stitch_830G_sensfunc_v0(sflist): #_v0

    file1_det3_sf = sanity_check_sf(sflist[0].sens[0])
    file2_det3_sf = sanity_check_sf(sflist[1].sens[0])
    file3_det3_sf = sanity_check_sf(sflist[2].sens[0])
    file2_det7_sf = sanity_check_sf(sflist[1].sens[1])
    file3_det7_sf = sanity_check_sf(sflist[2].sens[1])

    ### Take the sensfuncs from first detector in file 1 and in file 2, merging them at the point where they intersect

    area_to_stitch = (5100, 5400)
    # Compute closest point, within the stitching area
    diff_mask1 = (file1_det3_sf['SENS_WAVE'] > area_to_stitch[0]) & (file1_det3_sf['SENS_WAVE'] < area_to_stitch[1])
    diff_mask2 = (file2_det3_sf['SENS_WAVE'] > area_to_stitch[0]) & (file2_det3_sf['SENS_WAVE'] < area_to_stitch[1])
    min_index = np.argmin(np.fabs(file1_det3_sf['SENS_ZEROPOINT_FIT'][diff_mask1] - file2_det3_sf['SENS_ZEROPOINT_FIT'][diff_mask2]))
    stitch_point = file1_det3_sf['SENS_WAVE'][diff_mask1][min_index]
    offset = file1_det3_sf['SENS_ZEROPOINT_FIT'][diff_mask1][min_index] - file2_det3_sf['SENS_ZEROPOINT_FIT'][diff_mask2][min_index]

    combined_sf1 = stitch_sf(file1_det3_sf, file2_det3_sf, stitch_point)

    ### Now translate up the sensfunc from det 3 in file 3 to match the joined sensfunc and 
    ### stitch it on at around 6500
    area_to_stitch = (6300, 6700)

    # Compare gradients to find the point where the two functions are closest in slope
    combined_sf2, stitch_point, offset = gradient_stitch(combined_sf1, file3_det3_sf, area_to_stitch)

    ### Now translate sensfunc 7 from file 2 down to match the end of the combined sensfunc and stitch them 
    ### together at around 7750
    area_to_stitch = (7700, 7900)
    combined_sf3, stitch_point, offset = gradient_stitch(combined_sf2, file2_det7_sf, area_to_stitch)

    # Finish the stitched sensfunc by adding the non-overlapping portions of det 7 from file 3.
    stitch_point = combined_sf3['SENS_WAVE'][-1]
    overlap_mask = file3_det7_sf['SENS_WAVE'] > stitch_point
    offset = combined_sf3['SENS_ZEROPOINT_FIT'][-1] - file3_det7_sf['SENS_ZEROPOINT_FIT'][overlap_mask][0]

    combined_sf4 = stitch_sf(combined_sf3, translate_sf(file3_det7_sf, offset))
    return (combined_sf4, None, None)

def stitch_830G_sensfunc(sflist): #_v5

    file1_det3_sf = sanity_check_sf(sflist[0].sens[0])
    file2_det3_sf = sanity_check_sf(sflist[1].sens[0])
    file3_det3_sf = sanity_check_sf(sflist[2].sens[0])
    file2_det7_sf = sanity_check_sf(sflist[1].sens[1])
    file3_det7_sf = sanity_check_sf(sflist[2].sens[1])
    file4_det7_sf = sanity_check_sf(sflist[3].sens[1])

    ### Take the sensfuncs from first detector in file 1 and in file 2, merging them at the point with the closest slope

    area_to_stitch = (5100, 5400)
    combined_sf1, stitch_point, offset = gradient_stitch(file1_det3_sf, file2_det3_sf, area_to_stitch)

    ### Now translate up the sensfunc from det 3 in file 3 to match the joined sensfunc and 
    ### stitch it on at around 6500
    area_to_stitch = (6300, 6700)

    # Compare gradients to find the point where the two functions are closest in slope
    combined_sf2, stitch_point, offset = gradient_stitch(combined_sf1, file3_det3_sf, area_to_stitch)

    ### Now translate sensfunc 7 from file 2 down to match the end of the combined sensfunc and stitch them 
    ### together at around 7750
    area_to_stitch = (7700, 7900)
    combined_sf3, stitch_point, offset = gradient_stitch(combined_sf2, file2_det7_sf, area_to_stitch)

    area_to_stitch = (8300, 8600)
    combined_sf4, stitch_point, offset = gradient_stitch(combined_sf3, file3_det7_sf, area_to_stitch)

    area_to_stitch=(9800, 9925)
    combined_sf5, stitch_point, offset = gradient_stitch(combined_sf4, file4_det7_sf, area_to_stitch)


    return (combined_sf5, None, None)