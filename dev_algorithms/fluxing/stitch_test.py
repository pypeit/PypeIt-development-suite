# Meant to be used in ipython when developing stitched sensfuncs.  
# To use:
# 1 Copy stitch sensfunc funtion from sensfunc_archive/stitchdeimos.py
# 2 Edit the files_XXXX arrays to contain the correct files.
# 3 Change the stitch function for the new stitched sensfunc. Change it to return None
#   as the fith return value.
# 4 Paste this file into ipython
# 5 To stitch and graph:
#    sflist = prep(files_600ZD, 'sens/extract')
#    targetsf, polyfit_areas, fit_info = stitch(sflist, "600ZD")
#    graph(targetsf, polyfit_areas=polyfit_areas, fit_info=fit_info, sflist=sflist
# 6 When determining how to do a polynomial fit, pass in None as use_polyfit_min/use_polyfit_max to
#   stitch_with_polyfit. Then change the stitch function to return None
#   for "polyfit_areas", and fit_info from stitch_with_polyfit as the fith return value. Use this graph
#   to tweak the polyfit and determine what range to overwrite with the fit. 
# 7 Once the polynomial fit is working, pass the desired range to fit using use_polyfit_min/use_polyfit_max, 
#   return polyfit_areas as the fourth return value, and None as the fit_info (fifth return value)

from math import comb
import sys
import os

import numpy as np
from numpy.lib.polynomial import poly

from pypeit.sensfunc import SensFunc
from pypeit import utils
from matplotlib import pyplot as plot
import copy
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


files_1200G = ['1200G/sens_2023jan17_d0117_0102.fits',
               '1200G/sens_2023jan17_d0117_0103.fits',
               '1200G/sens_2023jan17_d0117_0104.fits',
               '1200G/sens_2023jan17_d0117_0105.fits']

files_1200B = ['1200B/sens_2017oct03_d1003_0055.fits',
               '1200B/sens_2017oct03_d1003_0056.fits',
               '1200B/sens_2017oct03_d1003_0057.fits']

files_900ZD = ['900ZD/sens_2010oct06_d1006_0138.fits',
               '900ZD/sens_2010oct06_d1006_0139.fits',
               '900ZD/sens_2010oct06_d1006_0140.fits']

files_600ZD = ['600ZD/sens_2023jan17_d0117_0098.fits',
               '600ZD/sens_2023jan17_d0117_0099.fits',
               '600ZD/sens_2023jan17_d0117_0101.fits']


files_830G =  ['830G/sens_2010oct06_d1006_0135.fits',
               '830G/sens_2010oct06_d1006_0136.fits',
               '830G/sens_2010oct06_d1006_0137.fits']

def prep(files, fileroot=''):
    """ Builds a list of sensfunc objects from files..
    """
    filelist = [os.path.join(fileroot, file) for file in files]
    sflist = []
    for file in filelist:
        sf = SensFunc.from_file(file)
        sf.sensfile = file
        sflist.append(sf)

    return sflist

def stitch(sflist, grating):

    (combined_wave, combined_zp_fit, combined_zp_fit_gpm, polyfit_areas, fit_info) = stitch_sensfunc(grating, sflist)


    newsens = SensFunc.empty_sensfunc_table(1, len(combined_wave))
    newsens['SENS_WAVE'] = combined_wave
    newsens['SENS_ZEROPOINT_FIT'] = combined_zp_fit
    newsens['SENS_ZEROPOINT_FIT_GPM'] = combined_zp_fit_gpm

    targetsf = copy.deepcopy(sflist[0])
    targetsf.sens = newsens
    targetsf.wave = np.empty((combined_wave.size,1))
    targetsf.wave[:,0] = combined_wave
    targetsf.zeropoint = np.empty((combined_zp_fit.size,1))
    targetsf.zeropoint[:,0] = combined_zp_fit

    return targetsf, polyfit_areas, fit_info

def fit(sf, invar, order, bound, weights=None):
    combined_wave = sf.sens['SENS_WAVE'][0]
    combined_zp_fit = sf.sens['SENS_ZEROPOINT_FIT'][0]

    invvar = np.full_like(combined_wave, invar)
    # 5050 - 6500
    #if weights is None:
    #    weights = np.ones(combined_wave.size, dtype=float)
    #from IPython import embed
    #embed()
    for i in range(len(combined_wave)):
        if i != 0:
            if combined_wave[i-1] > combined_wave[i]:
                raise ValueError(f"Blah {i}")
    
    #zpfpf = fitting.robust_fit(combined_wave, combined_zp_fit, order, 
    #                           use_mad=False, maxiter=30, upper=bound, lower=bound, invvar=invvar, 
    #                           weights=weights, sticky=False)
    zpfpf = np.polynomial.polynomial.polyfit(combined_wave, combined_zp_fit, order)
    #zpfpf.to_file(args.output + f"_zpfit_fit_{order}.fits", overwrite=True)

    return zpfpf, invvar

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

def create_err_plot(axis, color, label, x, y, yerr, linewidth=2.5, linestyle='solid'):

    axis.errorbar(x, y, yerr=yerr, color=color, linewidth=linewidth, linestyle=linestyle, label=label)
    xmin = (0.98*x.min())
    xmax = (1.02*x.max())
    #ymin = 14.0
    ymin = 16.0
    ymax = (1.05*y.max())

    return xmin, xmax, ymin, ymax


def setup_axes(axis,xmin, ymin, xmax, ymax, title):
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    axis.legend()
    axis.set_xlabel('Wavelength (Angstroms)')
    axis.set_ylabel('Zeropoint Fit (AB mag)')
    axis.set_title('PypeIt Zeropoint Fit for ' + title)

def graph(sf, polyfit_areas, sflist = None, fit_info=None, label = "stitch test"):

    fig, axis = build_figure()

    sens_gpm = sf.sens['SENS_ZEROPOINT_FIT_GPM'][0]
    bpm = np.logical_not(sens_gpm)

    num_plots = 1

    if np.sum(bpm) != 0: # Only plot bad pixels if there are some
        num_plots += 1
        
    if polyfit_areas is not None:
        num_plots +=1

    if fit_info is not None:
        num_plots +=1

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

    # Plot the whole polynomial fit in cyan
    if fit_info is not None:
        xmask = np.logical_and(sf.sens['SENS_WAVE'][0] >= fit_info[0], sf.sens['SENS_WAVE'][0] <= fit_info[1])
        x = sf.sens['SENS_WAVE'][0][xmask]
        fitted_y = np.polynomial.polynomial.polyval(x, fit_info[2])# fit.eval(x)
        (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, "cyan", f'entire polyfit', x, fitted_y, marker='.', linestyle='none')
        i+=1



    # Plot the original sensfuncs in light gray in the background
    if sflist is not None:
        for bksf in sflist:
            for det in [0, 1]:
                gpm = bksf.sens['SENS_ZEROPOINT_FIT_GPM'][det]
                x = bksf.sens['SENS_WAVE'][det][gpm]
                y = bksf.sens['SENS_ZEROPOINT_FIT'][det][gpm]
                if i == 3:
                    label = "Original sensfuncs"
                else:
                    label = ""
                (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, (.5, .5, .5), label, x,y, linewidth = 1, linestyle='solid')#linestyle='dashed')
                i+=1

    setup_axes(axis, np.min(xmin), np.min(ymin), np.max(xmax),np.max(ymax), label)
    plot.show()


def old_graph(sf, fit_info, polyfit_, sflist=None, grating="1200G"):#, sens_gpm):

    fig, axis = build_figure()

    sens_gpm = sf.sens['SENS_ZEROPOINT_FIT_GPM'][0]
    bpm = np.logical_not(sens_gpm)
    #bpm_x = x[bpm
    num_plots = 2 if fit_info is not None else 1
    if combined_edge_mask is not None:
        num_plots += 1

    if np.sum(bpm) != 0: # Only plot bad pixels if there are some
        num_plots += 1        

    #invvar = np.full_like(combined_wave, invar)
    if sflist is not None:
        num_plots += (2*len(sflist))
    xmin = np.zeros(num_plots)
    ymin = np.zeros(num_plots)
    xmax = np.zeros(num_plots)
    ymax = np.zeros(num_plots)

    i = 0
    #sigma = 1.0/np.sqrt(invvar[0])
    #(xmin[i], xmax[i], ymin[i], ymax[i]) = create_err_plot(axis, (1, .7, .7, 0.5), f'error bars, sigma {sigma:.4f}', x,sf.sens['SENS_ZEROPOINT_FIT'][0], yerr=sigma)#, linestyle='dashed') )
    #i+=1

    x = sf.sens['SENS_WAVE'][0][sens_gpm]
    (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, (1, .7, .7, .5), f'Translated/stitched SENS_ZEROPOINT_FIT', x,sf.sens['SENS_ZEROPOINT_FIT'][0][sens_gpm], marker='.', linestyle='none')
    i+=1


    if np.sum(bpm) != 0: # Only plot bad pixels if there are some
        x = sf.sens['SENS_WAVE'][0][bpm]
        y = sf.sens['SENS_ZEROPOINT_FIT'][0][bpm]
        (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, "red", f'Inverse SENS_ZEROPOINT_FIT_GPM', x, y, marker='.', linestyle='none')
        i+=1

    if combined_edge_mask is not None:
        x = sf.sens['SENS_WAVE'][0][combined_edge_mask]
        y = sf.sens['SENS_ZEROPOINT_FIT'][0][combined_edge_mask]
        (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, "blue", f'Polynomial fit', x, y, marker='.', linestyle='none')
        i+=1


    if fit_info is not None:
        xmask = np.logical_and(sf.sens['SENS_WAVE'][0] >= fit_info[0], sf.sens['SENS_WAVE'][0] <= fit_info[1])
        x = sf.sens['SENS_WAVE'][0][xmask]
        fitted_y = np.polynomial.polynomial.polyval(x, fit_info[2])# fit.eval(x)
        # Old color was .5, .5, .5
        (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, (0, 0, 0), f'fitted zeropoint', x, fitted_y, linestyle='dashed')
        i+=1

    
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

    #x = sf.sens['SENS_WAVE'][0][np.logical_not(sens_gpm)]
    #y = sf.sens['SENS_ZEROPOINT_FIT'][0][np.logical_not(sens_gpm)]
    #(xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, "yellow", f'SensFunc zeropoint fit gpm', x, y, marker='^', linestyle='none')
    #i+=1

    #if num_plots == 4:
    #    y = sf.sens['SENS_ZEROPOINT_FIT'][0][bpm]
    #    (xmin[i], xmax[i], ymin[i], ymax[i]) = create_plot(axis, "orange", f'PypeItFit gpm', bpm_x, y, marker='o', linestyle='none')
    #    i+=1

    setup_axes(axis, np.min(xmin), np.min(ymin), np.max(xmax),np.max(ymax), grating)


    plot.show()

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

