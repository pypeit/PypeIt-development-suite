import numpy as np
import os

from astropy.table import Table
from astropy.io import ascii,fits

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

from pypeit.core import fitting
from pypeit.core.wavecal import wvutils
#from astropy.modeling import models, fitting

import dmost_utils, dmost_slit_matching

from IPython import embed

#DEIMOS_DROPBOX = '/Users/mgeha/Dropbox/DEIMOS/'


########################################
def gaussian(x,*p) :
    # A gaussian peak with:
    #   Constant Background          : p[0]
    #   Peak height above background : p[1]
    #   Central value                : p[2]
    #   Standard deviation           : p[3]
    return p[0]+p[1]*np.exp(-1.*(x-p[2])**2/(2.*p[3]**2))


########################################
def gauss_guess(x,y):
    norm = np.median(np.percentile(y,50))
    w=np.mean(x)
    N_guess   = np.max(y) - np.min(y)
    sig_guess = 0.5
    p0 = [norm,N_guess,w,sig_guess]

    return p0


########################################
def spec_normalize(wave,spec,ivar,nlow=75,npoly=1):
    
    m = (spec > np.percentile(spec,nlow)) & (spec < np.percentile(spec,99))

    z = np.polyfit(wave[m],spec[m],npoly)
    p = np.poly1d(z)
    fit = p(wave)
    
    # NORMALIZE SPECTRUM
    nwave = wave
    nspec = spec/fit
    nivar = ivar*fit**2

    return nwave,nspec,nivar


#####################################################
# CALCULATE SKY EMISSION LINE 
#
def sky_em_residuals(wave,flux,ivar,plot=0, orig=True, toler=0.3):

    # READ SKY LINES -- THESE ARE VACUUM WAVELENGTHS
    #sky_file = DEIMOS_DROPBOX+'Other_data/sky_single_mg.dat'
    # TODO -- Path in PypeIt
    sky_file = 'sky_single_mg.dat'
    sky=ascii.read(sky_file)
    
    dwave = []
    diff = []
    diff_err  = []
    los = []
    los_err= []

    if orig:
        for line in sky['Wave']:
            wline = [line-5.,line+5.] 
            mw    = (wave > wline[0]) & (wave < wline[1])
            
            p=[0,0,0,0]
            if np.sum(mw) > 20:
                p0 = gauss_guess(wave[mw],flux[mw])
                try:
                    p, pcov = curve_fit(gaussian,wave[mw],flux[mw], 
                                        sigma = 1./np.sqrt(ivar[mw]), p0=p0)
                    perr = np.sqrt(np.diag(pcov))
                except:
                    p=p0
                    p[2] = -99
                    perr=p0
                gfit = gaussian(wave[mw],*p)
                d = p[2] - line

                if (plot==1):
                    plt.figure(figsize=(8,3)) 
                    plt.plot(wave[mw],gfit,'g')
                    plt.plot(wave[mw],flux[mw])
                    plt.title('{} {:0.2f} diff= {:0.3f}'.format(line,p[3],d))

                if ~np.isfinite(perr[2]):
                    perr[2] = 1000.
                dwave = np.append(dwave,line)
                diff = np.append(diff,d)
                diff_err = np.append(diff_err,perr[2])
                los = np.append(los,p[3])
                los_err = np.append(los_err,perr[3])
    # New approach
    all_tcent, all_ecent, cut_tcent, icut, arc_cont_sub, \
        all_twid, all_tampl = wvutils.arc_lines_from_spec(
            flux, sigdetect=5., return_more=True)

    # Linearly interpolate wavelengths 
    # Match up to known lines
    dwave2 = []
    diff2 = []
    diff_err2  = []
    los2 = []
    f_wv = interp1d(np.arange(wave.size), wave)
    gd_wave = f_wv(all_tcent[icut])
    for line in sky['Wave']:
        if np.absmin(line-gd_wave) < toler:
            imin = np.argmin(np.abs(np.absline-gd_wave))
            dwave2.append(line)
            diff2.append(gd_wave[imin]-line)
            diff_err2.append(all_ecent[icut][imin])
            los2.append(all_twid[icut][imin])


    embed(header='compare these on 113 of dmost_flexure')
            
    m=(diff_err < 0.1) & (diff_err > 0.0)
    return dwave[m],diff[m],diff_err[m],los[m],los_err[m]


#######################################################
#
def fit_sky_linear(wlines,wdiff,wdiff_err):
    
    z = np.polyfit(wlines,wdiff,w = 1./wdiff_err,deg= 1)
    p = np.poly1d(z)

    return p

#######################################################
#  
#
def measure_sky_lines(slits, nslits, hdu):

    for i in np.arange(0,nslits,1):

        if (slits['rSN'][i] > 1.) & (slits['bSN'][i] > 1.):
            r = slits['rname'][i]
            b = slits['bname'][i]


            # SKY LINES FIRST
            r_sky_line, r_sky_diff,r_sky_ediff,r_los,r_elos = sky_em_residuals(hdu[r].data['OPT_WAVE'], \
                                                    hdu[r].data['OPT_COUNTS_SKY'],\
                                                    hdu[r].data['OPT_COUNTS_IVAR'])

            b_sky_line, b_sky_diff,b_sky_ediff,b_los,b_elos = sky_em_residuals(hdu[b].data['OPT_WAVE'], \
                                                    hdu[b].data['OPT_COUNTS_SKY'],\
                                                    hdu[b].data['OPT_COUNTS_IVAR'])


            sky_diff  = np.concatenate((r_sky_diff,b_sky_diff),axis=None)
            sky_lines = np.concatenate((r_sky_line,b_sky_line),axis=None)
            sky_ediff = np.concatenate((r_sky_ediff,b_sky_ediff),axis=None)
            sky_los   = np.concatenate((r_los,b_los),axis=None)

            # FIT SINGLE SLIT SKY LINES WITH A LINE           
            fitted_line = fit_sky_linear(sky_lines,sky_diff,sky_ediff)

            
            slits['fit_slope'][i] = fitted_line[1]
            slits['fit_b'][i]     = fitted_line[0]
            slits['fit_los'][i]   = np.median(sky_los)


    return slits

#######################################################
#  UPDATE SLITS FITS WITH ALL_MASK FIT
#
def update_flexure_fit(slits, nslits, hdu, pmodel_m,pmodel_b,pmodel_los):

    fslits = slits

    # UPDATE FITS
#    fslits['fit_slope'] = pmodel_m(slits['xpos'],slits['ypos'])
#    fslits['fit_b']     = pmodel_b(slits['xpos'],slits['ypos'])
#    fslits['fit_los']   = pmodel_los(slits['xpos'],slits['ypos'])

    fslits['fit_slope'] = pmodel_m(slits['objra'],slits['objdec'])
    fslits['fit_b']     = pmodel_b(slits['objra'],slits['objdec'])
    fslits['fit_los']   = pmodel_los(slits['objra'],slits['objdec'])

    # CALCULATE RESIDUALS FROM FIT
    resid_sky = []
    for f in fslits:

        if f['rSN'] > 0:
            all_wave,all_flux,all_ivar,all_sky = dmost_utils.load_spectrum(f,hdu,vacuum = 1)

            dwave,diff,diff_err,los,elos = sky_em_residuals(all_wave,all_sky,all_ivar,plot=0)
            m=np.isfinite(diff)
            sky_mean = np.average(np.abs(diff[m]), weights = 1./diff_err[m]**2)
            resid_sky = np.append(resid_sky,sky_mean)
        if f['rSN'] <= 0:
            resid_sky = np.append(resid_sky,-1)

    fslits['resid_sky'] = resid_sky

    return fslits


#######################################################
def fit_mask_surfaces(slits):

    m = (slits['rSN'] > 1.) & (slits['bSN'] > 1.)

    mu =  np.median(slits['fit_slope'][m])
    sd =  np.std(slits['fit_slope'][m])
    mu2 =  np.median(slits['fit_b'][m])
    sd2 =  np.std(slits['fit_b'][m])


    mgood=(np.abs(slits['fit_slope']-mu) < 2.*sd)  & (np.abs(slits['fit_b']-mu2) < 2.*sd2) & \
                            (slits['rSN'] > 1.) & (slits['bSN'] > 1.)

    
    # FIT ALL SURFACES WITH 3D POLYNOMIAL
    p_init = models.Polynomial2D(degree=3)
    fit_p = fitting.LevMarLSQFitter()

    # FIT FOR SLOPES, INTERCEPTS, LOS
    
    pmodel_m = fit_p(p_init, slits['objra'][mgood], slits['objdec'][mgood], slits['fit_slope'][mgood])
    pmodel_b = fit_p(p_init, slits['objra'][mgood], slits['objdec'][mgood], slits['fit_b'][mgood])
    pmodel_los = fit_p(p_init, slits['objra'][mgood], slits['objdec'][mgood], slits['fit_los'][mgood])

    
    return pmodel_m,pmodel_b,pmodel_los



#######################################################

def qa_flexure_plots(plot_dir, nslits, slits, fslits, hdu):

    header = hdu[0].header
    mask = header['TARGET'].strip()
    fnames = header['FILENAME'].split('.')
    pdf2 = matplotlib.backends.backend_pdf.PdfPages(plot_dir+'QA/flex_slits_'+mask+'_'+fnames[2]+'.pdf')
    plt.rcParams.update({'figure.max_open_warning': 0})
    for i in np.arange(0,nslits,1):

        if (slits['rSN'][i] > 0.) & (slits['bSN'][i] > 0.):
            r = slits['rname'][i]
            b = slits['bname'][i]


            # SKY LINES FIRST
            r_sky_line, r_sky_diff,r_sky_ediff,r_los,r_elos = sky_em_residuals(hdu[r].data['OPT_WAVE'], \
                                                    hdu[r].data['OPT_COUNTS_SKY'],\
                                                    hdu[r].data['OPT_COUNTS_IVAR'])

            b_sky_line, b_sky_diff,b_sky_ediff,b_los,b_elos = sky_em_residuals(hdu[b].data['OPT_WAVE'], \
                                                    hdu[b].data['OPT_COUNTS_SKY'],\
                                                    hdu[b].data['OPT_COUNTS_IVAR'])

            fig, (ax1,ax2) = plt.subplots(1, 2,figsize=(20,4))
            ax1.plot(r_sky_line,r_sky_diff,'ro',alpha=0.8,label='Red chip: Sky Emission')
            ax1.plot(b_sky_line,b_sky_diff,'bo',alpha=0.8,label='Blue chip: Sky Emission')
            ax1.errorbar(b_sky_line,b_sky_diff,yerr=b_sky_ediff,fmt='none',ecolor='b',alpha=0.5)
            ax1.errorbar(r_sky_line,r_sky_diff,yerr=r_sky_ediff,fmt='none',ecolor='r',alpha=0.5)
            ax1.text(6320,0,'{}'.format(b),fontsize=11)
            ax1.text(8500,0,'{}'.format(r),fontsize=11)
            ax1.set_ylim(-0.45,0.45)

            x=np.arange(6000,9000,1)
            l1 = slits['fit_slope'][i]*x + slits['fit_b'][i]
            l2 = fslits['fit_slope'][i]*x + fslits['fit_b'][i]
            ax1.plot(x,l1,'-')
            ax1.plot(x,l2,'--')
            ax1.axhline(linewidth=1, color='grey',alpha=0.5)
            ax1.set_ylabel('Wavelength offset (AA)')
            ax1.set_xlabel('Wavelength (AA)')
            ax1.set_xlim(6300,9100)
            t = 'Sky Line Fits , resid = {:0.4f} AA, arc = {:0.2f}'.format(slits['resid_sky'][i],0.32*slits['rms_arc_r'][i])
            ax1.set_title(t)

            sky_diff  = np.concatenate((r_sky_diff,b_sky_diff),axis=None)
            sky_lines = np.concatenate((r_sky_line,b_sky_line),axis=None)
            sky_ediff = np.concatenate((r_sky_ediff,b_sky_ediff),axis=None)
            sky_los   = np.concatenate((r_los,b_los),axis=None)


            ax2.plot(r_sky_line,r_los,'ro',alpha=0.8,label='Red chip: Sky Emission')
            ax2.plot(b_sky_line,b_los,'bo',alpha=0.8,label='Blue chip: Sky Emission')
            ax2.errorbar(r_sky_line,r_los,yerr=r_elos,fmt='none',ecolor='r',alpha=0.5)
            ax2.errorbar(b_sky_line,b_los,yerr=b_elos,fmt='none',ecolor='b',alpha=0.5)
            ax2.axhline(fslits['fit_los'][i],linewidth=1, color='grey',alpha=0.5)

            ax2.set_title('Line widths')
            ax2.set_xlabel('Wavelength (AA)')
            ax2.set_ylim(0.3,0.8)
            ax2.set_xlim(6300,9100)

            pdf2.savefig()
    pdf2.close()
    plt.close('all')

    #########################################################################
    # CREATE FULL MASK FITS
    pdf = matplotlib.backends.backend_pdf.PdfPages(plot_dir+'QA/flex_mask_'+mask+'_'+fnames[2]+'.pdf')
    xslit = slits['objra']
    yslit = slits['objdec']
    t=2.

    mu =  np.median(slits['fit_slope'])
    sd =  np.std(slits['fit_slope'])
    mu2 =  np.median(slits['fit_b'])
    sd2 =  np.std(slits['fit_b'])
    mu3 =  np.median(slits['fit_los'])
    sd3 =  np.std(slits['fit_los'])

    # PLOT FITTED VALUES
    fig, (ax1,ax2,ax3) = plt.subplots(1, 3,figsize=(22,5))
 
    mm1=-0.00005
    mm2=0.00005
    print(mu-t*sd,mu+t*sd)
    ax1.scatter(xslit,yslit,c=slits['fit_slope'],cmap="cool",vmin = mm1,vmax=mm2 )# mu-t*sd,vmax=mu+t*sd)
    ax1.set_ylabel('Dec [deg]')
    ax1.set_xlabel('RA [deg]')
    ax1.set_title('Wave MEASURE: line slope')
    #cax, _ = matplotlib.colorbar.make_axes(ax1)
    #normalize = matplotlib.colors.Normalize(vmin = mu-t*sd,vmax=mu+t*sd)
    #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)


    ax2.scatter(xslit,yslit,c=slits['fit_b'],cmap="summer",vmin = mu2-t*sd2,vmax=mu2+t*sd2)
    ax2.set_ylabel('Dec [deg]')
    ax2.set_xlabel('RA [deg]')
    ax2.set_title('Wave MEASURE: line intercept')
    cax, _ = matplotlib.colorbar.make_axes(ax2)
    normalize = matplotlib.colors.Normalize(vmin = mu2-t*sd2,vmax=mu2+t*sd2)
    #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='summer',norm=normalize)


    ax3.scatter(xslit,yslit,c=slits['fit_los'],cmap="cool",vmin = mu3-t*sd3,vmax=mu3+t*sd3)
    ax3.set_ylabel('Dec [deg]')
    ax3.set_xlabel('RA [deg]')
    ax3.set_title('Wave MEASURE: line width')
    cax, _ = matplotlib.colorbar.make_axes(ax3)
    normalize = matplotlib.colors.Normalize(vmin = mu3-t*sd3,vmax=mu3+t*sd3)
    #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)

    pdf.savefig()
    
    #######################
    # PLOT MEASURED VALUES
    fig, (ax1,ax2,ax3) = plt.subplots(1, 3,figsize=(22,5))
 
    ax1.scatter(xslit,yslit,c=fslits['fit_slope'],cmap="cool",vmin = mu-t*sd,vmax=mu+t*sd)

    ax1.set_ylabel('Dec [deg]')
    ax1.set_xlabel('RA [deg]')
    ax1.set_title('Wave fit: line slope')
    cax, _ = matplotlib.colorbar.make_axes(ax1)
    normalize = matplotlib.colors.Normalize(vmin = mu-t*sd,vmax=mu+t*sd)
    #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)


    ax2.scatter(xslit,yslit,c=fslits['fit_b'],cmap="summer",vmin = mu2-t*sd2,vmax=mu2+t*sd2)
    ax2.set_ylabel('Dec [deg]')
    ax2.set_xlabel('RA [deg]')
    ax2.set_title('Wave fit: line intercept')
    cax, _ = matplotlib.colorbar.make_axes(ax2)
    normalize = matplotlib.colors.Normalize(vmin = mu2-t*sd2,vmax=mu2+t*sd2)
    #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='summer',norm=normalize)


    ax3.scatter(xslit,yslit,c=fslits['fit_los'],cmap="cool",vmin = mu3-t*sd3,vmax=mu3+t*sd3)
    ax3.set_ylabel('Dec [deg]')
    ax3.set_xlabel('RA [deg]')
    ax3.set_title('Wave fit: line width')
    cax, _ = matplotlib.colorbar.make_axes(ax3)
    normalize = matplotlib.colors.Normalize(vmin = mu3-t*sd3,vmax=mu3+t*sd3)
    #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)

    
    pdf.close()




#######################################################
def flexure_correct(hdu, data_dir,clobber=0):


    filename = hdu.filename()
    tmp = filename.split('spec1d')
    
    txt = filename.replace('.fits','.txt')
    
    slit_table_file = data_dir + '/dmost/dmost'+tmp[1]

    # IF FILE DOESN"T EXIST GENERATE
    if (not os.path.isfile(slit_table_file)) | (clobber == 1):

        # CREATE SLIT TABLE
        slits, nslits = dmost_slit_matching.create_slit_table(hdu,data_dir,txt)

        # INITIAL SKY LINE STUFF
        slits = measure_sky_lines(slits, nslits,hdu)

        # FIT SURFACES
        pmodel_m, pmodel_b,pmodel_los = fit_mask_surfaces(slits)

     
        # ADD TO TABLE
        fslits = update_flexure_fit(slits,nslits, hdu, pmodel_m, pmodel_b,pmodel_los)

        # REFIT FOR QA PLOTS
        qa_flexure_plots(data_dir,nslits,slits,fslits,hdu)

        fslits.write(slit_table_file,overwrite=True)

    # ELSE READ IT IN
    if os.path.isfile(slit_table_file):
        fslits = Table.read(slit_table_file)

    return fslits

    
