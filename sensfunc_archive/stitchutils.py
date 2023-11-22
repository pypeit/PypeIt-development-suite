"""Utility functions to aid stitching sensfuncs together."""
import numpy as np
from pypeit.core.telluric import read_telluric_grid
from pypeit import utils
from pypeit import msgs
from pypeit.core.wavecal import wvutils

from matplotlib import pyplot as plot

def sanity_check_sf(sf):
    """Sanity check a sensitivity function and return a 
    sensfunc-like dictionary useable by the other functions in stitchutils.py
    
    Args:
        sf (astropy.table.Table): A PypeIt sensitivity function table, as found in
                                  the "sens" value in the SensFunc data container.

    Return (dict): A dictionary that contains the "SENS_WAVE", "SENS_ZEROPOINT_FIT",
                   and "SENS_ZEROPOINT_FIT_GPM" data from the passed in columns.
                   The data arrays are masked to avoid zero values.

    Raises: AssertionError - Raised if there are infinite or NaN values in the sensfunc.
    """
    # Return a sensfunc-like dict with the zero pixels masked
    mask = (sf['SENS_ZEROPOINT_FIT'] != 0) & (sf['SENS_WAVE'] != 0)

    masked_sf = {'SENS_WAVE':              sf['SENS_WAVE'][mask],
                 'SENS_ZEROPOINT_FIT':     sf['SENS_ZEROPOINT_FIT'][mask],
                 'SENS_ZEROPOINT_FIT_GPM': sf['SENS_ZEROPOINT_FIT_GPM'][mask],
                }
    
    # Sanity check our inputs
    assert np.isfinite(masked_sf['SENS_WAVE']).all()
    assert np.isfinite(masked_sf['SENS_ZEROPOINT_FIT']).all()

    return masked_sf

def truncate_sf(sf, min_wave=None, max_wave=None):
    """Truncates a sensfunc between a min and max wavelength.
    
    Args:
        min_wave (float): The minimum wavelength to truncate at. Values <= to this
                          value are not included in the resulting sensfunc. If
                          None no truncation is done at the start of the sensfunc.
                          Defaults to None.
        min_wave (float): The maximum wavelength to truncate at. Values >= to this
                          value are not included in the resulting sensfunc. If
                          None no truncation is done at the end of the sensfunc.
                          Defaults to None.

    Returns (dict):  A dictionary with the truncated sensfunc data.
    """
    if min_wave is None and max_wave is None:
        msgs.warn("Truncate called with none wavelengths, returning original sensfunc")
        return sf
    
    if min_wave is not None:
        truncate_mask = sf['SENS_WAVE'] > min_wave
    else:
        truncate_mask = np.ones_like(sf['SENS_WAVE'], dtype=bool)

    if max_wave is not None:
        truncate_mask = truncate_mask & (sf['SENS_WAVE'] < max_wave)

    truncated_sf = {'SENS_WAVE':              sf['SENS_WAVE'][truncate_mask],
                    'SENS_ZEROPOINT_FIT':     sf['SENS_ZEROPOINT_FIT'][truncate_mask],
                    'SENS_ZEROPOINT_FIT_GPM': sf['SENS_ZEROPOINT_FIT_GPM'][truncate_mask],}
    return truncated_sf

def translate_sf(sf, offset):
    """Translate a sensfunc by adding (or subtracting) an offset from it's zeropoint
    values.
    
    Args:
        offset (float): The offset to add to the sensfunc.

    Returns (dict):  A dictionary with the translate sensfunc data.
    
    """
    translated_sf =  {'SENS_WAVE':              sf['SENS_WAVE'],
                      'SENS_ZEROPOINT_FIT':     sf['SENS_ZEROPOINT_FIT']+offset,
                      'SENS_ZEROPOINT_FIT_GPM': sf['SENS_ZEROPOINT_FIT_GPM']}
    return translated_sf

def stitch_sf(sf1, sf2, stitch_point=None):
    """Stitches (appends) two sensitivity function together.
    If no stitch_point is given, the non-overlapping portions of sf2 are
    appended to sf1. Otherwise values < the stitch point from sf1 are
    combined with values >= the stitch point from sf2.
    Args:
        sf1 (dict): The first sensitivity function.
        sf2 (dict): The second sensitivity function.
        stitch_point (float): The wavelength value to do the stitching at. Defaults to None.

    Returns (dict): A new sensfunc with the combined data from sf1 and sf2.
    """
    if stitch_point is None:
        # Just append the portion of sf2 that doesn't overlap sf1
        sf1_mask = np.ones_like(sf1['SENS_WAVE'], dtype=bool)
        sf2_mask = sf2['SENS_WAVE'] > np.max(sf1['SENS_WAVE'])
    else:
        sf1_mask = sf1['SENS_WAVE'] < stitch_point
        sf2_mask = sf2['SENS_WAVE'] >= stitch_point

    combined_wave       = np.concatenate((sf1['SENS_WAVE'][sf1_mask],              sf2['SENS_WAVE'][sf2_mask]))
    combined_zp_fit     = np.concatenate((sf1['SENS_ZEROPOINT_FIT'][sf1_mask],     sf2['SENS_ZEROPOINT_FIT'][sf2_mask]))
    combined_zp_fit_gpm = np.concatenate((sf1['SENS_ZEROPOINT_FIT_GPM'][sf1_mask], sf2['SENS_ZEROPOINT_FIT_GPM'][sf2_mask]))

    return {'SENS_WAVE':              combined_wave,
            'SENS_ZEROPOINT_FIT':     combined_zp_fit,
            'SENS_ZEROPOINT_FIT_GPM': combined_zp_fit_gpm }

def gradient_stitch_sf(sf1, sf2, area_to_stitch):
    """Stitches two sensitivity function together.
    It searches a given wavelength wave for the point where the gradient of
    the two sensfuncs are closest, in an attempt to get a smooth joining.

    Args:
        sf1 (dict): The first sensitivity function.
        sf2 (dict): The second sensitivity function.
        area_to_stitch (tuple of float): The min/max wavelenghts of the area to search for a stitching point.

    Returns (dict): A new sensfunc with the combined data from sf1 and sf2.
    """

    sf1_g = np.gradient(sf1['SENS_ZEROPOINT_FIT'], sf1['SENS_WAVE'])
    sf2_g = np.gradient(sf2['SENS_ZEROPOINT_FIT'], sf2['SENS_WAVE'])

    diff_mask1 = (sf1['SENS_WAVE'] > area_to_stitch[0]) & (sf1['SENS_WAVE'] < area_to_stitch[1])
    diff_mask2 = (sf2['SENS_WAVE'] > area_to_stitch[0]) & (sf2['SENS_WAVE'] < area_to_stitch[1])
    min_index = np.argmin(np.fabs(sf1_g[diff_mask1] - sf2_g[diff_mask2]))
    stitch_point = sf1['SENS_WAVE'][diff_mask1][min_index]
    offset = sf1['SENS_ZEROPOINT_FIT'][diff_mask1][min_index] - sf2['SENS_ZEROPOINT_FIT'][diff_mask2][min_index]
    combined_sf = stitch_sf(sf1, translate_sf(sf2, offset), stitch_point)
    return combined_sf, stitch_point, offset


def stitch_sf_polyfit_old(sf1, sf2,
                          fitting_wave, fitting_zp_fit, order,
                          use_polyfit_min, use_polyfit_max):
    """Stitches two sensitivity function together using a polynomial fit
    to smooth the join. 
    
    This version requires exact values for where the polynomial fit values
    overwrite the original sensfunc data. It also doesn't use telgrid files
    to generate a wavelength grid over the stitched region.

    Args:
        sf1 (dict): The first sensitivity function.
        sf2 (dict): The second sensitivity function.

        fitting_wave (np.array):   The wavelength values to use for the polynomial fit.
                                   Typically this is truncated and has any detector 
                                   boundaries removed to get the closesit fit.

        fitting_zp_fit (np.array): The zeropoint values to use for the polynomial fit.
                                   Typically this is truncated and has any detector 
                                   boundaries removed to get the closesit fit.

        order (int):               The order to use for the fitting polynomial.

        use_polyfit_min (float):   The minimum wavelength in the first sensfunc to overwrite
                                   with the polynomial fit.

        use_polyfit_max (float):   The maximum wavelength in the second sensfunc to
                                   overwrite with the polynomial fit.
                                 

    Returns (dict):  A new sensfunc with the combined data from sf1 and sf2 and the polynomial fit.
            (tuple): A tuple containing the min and max wavelenghts of the polynomial fit,
                     and an array of constants for the polynomial used.
    """

    wave_1, zp_fit_1, gpm_1 = sf1['SENS_WAVE'], sf1['SENS_ZEROPOINT_FIT'], sf1['SENS_ZEROPOINT_FIT_GPM']
    wave_2, zp_fit_2, gpm_2 = sf2['SENS_WAVE'], sf2['SENS_ZEROPOINT_FIT'], sf2['SENS_ZEROPOINT_FIT_GPM']

    fitc = np.polynomial.polynomial.polyfit(fitting_wave, fitting_zp_fit, order)

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

    combined_sf = {
        'SENS_WAVE':              combined_wave,
        'SENS_ZEROPOINT_FIT':     combined_zp_fit,
        'SENS_ZEROPOINT_FIT_GPM': combined_zp_fit_gpm
    }


    return combined_sf, (use_polyfit_min, use_polyfit_max, fitc)

def stitch_sf_polyfit(sf1,sf2, 
                      sf1_fitting_mask, sf2_fitting_mask,
                      order,
                      polyfit_region,
                      gen_wave_grid=False, telgrid_file=None):
    """Stitches two sensitivity function together using a polynomial fit
    to smooth the join. 
    
    This version only requires an approximate region to replace with the polynomial fit.
    It chooses values within that region where the fit is closest to the original data.
    It also can generate a wavelength grid for the polynomial fit region using a telgrid file.

    Args:
        sf1 (dict): The first sensitivity function.
        sf2 (dict): The second sensitivity function.

        sf1_fitting_mask (np.array):   
            A mask for sf1 to use for the polynomial fit.
            Typically this is truncated and has any detector 
            boundaries removed to get the closesit fit.

        sf2_fitting_mask (np.array):   
            A mask for sf2 to use for the polynomial fit.
            Typically this is truncated and has any detector 
            boundaries removed to get the closesit fit.

        order (int):
            The order to use for the fitting polynomial.

        polyfit_region (tuple of float):
            The approximate wavelength range within the resulting sensfunc to replace
            with the polynomial fit.

        gen_wave_grid (bool, Optional):   
            Whether to generate a new wavelength grid for the polynomial fit region.
            Defaults to False
                                 
        telgrid_file (str, Optional):
            The telluric file to use when generating a wavelength grid. Defaults to None.
            Must be given if gen_wave_grid is True.

    Returns (dict):  A new sensfunc with the combined data from sf1 and sf2 and the polynomial fit.
            (tuple): A tuple containing the min and max wavelengths within polyfit_region that
                     were replaced by the polynomial fit, and an array of
                     constants for the polynomial used.
    """

    if gen_wave_grid is True:
        if telgrid_file is None:
            raise ValueError("Must specify telgrid_file if gen_wave_grid is True")


    fitting_wave = np.concatenate((sf1['SENS_WAVE'][sf1_fitting_mask], sf2['SENS_WAVE'][sf2_fitting_mask]))
    fitting_zp_fit = np.concatenate((sf1['SENS_ZEROPOINT_FIT'][sf1_fitting_mask], sf2['SENS_ZEROPOINT_FIT'][sf2_fitting_mask]))

    fitc = np.polynomial.polynomial.polyfit(fitting_wave, fitting_zp_fit, order)

    if gen_wave_grid:
        wave_grid_min = sf1['SENS_WAVE'][-1]
        wave_grid_max = sf2['SENS_WAVE'][0]
        msgs.info(f"Generating wave grid between: {wave_grid_min} and {wave_grid_max}")
        tell_dict = read_telluric_grid(telgrid_file, wave_min=wave_grid_min, wave_max=wave_grid_max)
        # Sometimes read_telluric grid returns a larger grid
        grid_mask = (tell_dict['wave_grid'] > wave_grid_min) & (tell_dict['wave_grid'] < wave_grid_max)
        joining_wave = tell_dict['wave_grid'][grid_mask]
        combined_wave = np.concatenate((sf1['SENS_WAVE'], joining_wave, sf2['SENS_WAVE']))
    else:
        combined_wave = np.concatenate((sf1['SENS_WAVE'], sf2['SENS_WAVE']))

    fit_region_mask = (combined_wave >= polyfit_region[0]) & (combined_wave <= polyfit_region[1])    
    polyfit_zp = np.polynomial.polynomial.polyval(combined_wave[fit_region_mask], fitc)

    sf1_fit_diff_mask = combined_wave[fit_region_mask] <= sf1['SENS_WAVE'][-1]
    sf1_diff_mask = sf1['SENS_WAVE'] >= polyfit_region[0]
    closest_idx = np.argmin(np.fabs(sf1['SENS_ZEROPOINT_FIT'][sf1_diff_mask] - polyfit_zp[sf1_fit_diff_mask]))
    use_polyfit_min = sf1['SENS_WAVE'][sf1_diff_mask][closest_idx]

    sf2_fit_diff_mask = combined_wave[fit_region_mask] >= sf2['SENS_WAVE'][0]
    sf2_diff_mask = sf2['SENS_WAVE'] <= polyfit_region[1]
    closest_idx = np.argmin(np.fabs(sf2['SENS_ZEROPOINT_FIT'][sf2_diff_mask] - polyfit_zp[sf2_fit_diff_mask]))
    use_polyfit_max = sf2['SENS_WAVE'][sf2_diff_mask][closest_idx]

    use_sf1_mask = sf1['SENS_WAVE'] < use_polyfit_min
    use_sf2_mask = sf2['SENS_WAVE'] > use_polyfit_max
    use_polyfit_mask = (combined_wave[fit_region_mask] >= use_polyfit_min) & (combined_wave[fit_region_mask] <= use_polyfit_max)

    combined_zp_fit = np.concatenate((sf1['SENS_ZEROPOINT_FIT'][use_sf1_mask],
                                      polyfit_zp[use_polyfit_mask],
                                      sf2['SENS_ZEROPOINT_FIT'][use_sf2_mask]))

    combined_zp_fit_gpm = np.concatenate((sf1['SENS_ZEROPOINT_FIT_GPM'][use_sf1_mask],
                                          np.ones_like(polyfit_zp[use_polyfit_mask], dtype=bool),
                                          sf2['SENS_ZEROPOINT_FIT_GPM'][use_sf2_mask]))

    combined_sf = {
        'SENS_WAVE':              combined_wave,
        'SENS_ZEROPOINT_FIT':     combined_zp_fit,
        'SENS_ZEROPOINT_FIT_GPM': combined_zp_fit_gpm
    }


    return combined_sf, (use_polyfit_min, use_polyfit_max, fitc)

