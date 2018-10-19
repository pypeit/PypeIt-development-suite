""" 2D ARCS
"""

from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

# General imports
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.io import readsav
from astropy.io import ascii
from astropy.stats import sigma_clip

# PYPEIT imports
from pypeit.core import pydl
from pypeit import msgs

def prettyplot():
    # set some plotting parameters
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams["ytick.direction"] = 'in'
    plt.rcParams["xtick.direction"] = 'in'
    plt.rcParams["xtick.major.size"] = 6
    plt.rcParams["ytick.major.size"] = 6
    plt.rcParams["xtick.minor.size"] = 3
    plt.rcParams["ytick.minor.size"] = 3
    plt.rcParams["xtick.major.width"] = 1
    plt.rcParams["ytick.major.width"] = 1
    plt.rcParams["xtick.minor.width"] = 1
    plt.rcParams["ytick.minor.width"] = 1
    plt.rcParams["axes.linewidth"] = 1
    plt.rcParams["lines.linewidth"] = 5.5
    plt.rcParams["lines.markeredgewidth"] = 2.0
    plt.rcParams["patch.linewidth"] = 3
    plt.rcParams["hatch.linewidth"] = 3.0
    plt.rcParams["font.size"] = 13
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["legend.handletextpad"] = 1.0


def fit2darc(all_wv,
             all_pix,
             t,
             nycoeff=3,
             nocoeff=5,
             sigmarjct=3.0,
             debug=True):
    """Routine to obtain the 2D wavelength solution for an
    echelle spectrograph. This is calculated from the y-centroid
    and the order number of identified arc lines. The fit is a 
    simple least-squares with one round of rejections.
    This is a direct porting of the XIDL code: x_fit2darc.pro

    Parameters
    ----------
    all_wv: np.array
     wavelength of the identified lines
    all_pix: np.array
      y-centroid position of the identified lines
    t: np.array
      order number of the identified lines
    nycoeff : np.int
      order of the fitting along the pixel direction for each order
    nocoeff : np.int
      order of the fitting in the order direction
    sigmarjct: np.float
      sigma level for the rejection
    debug: boolean
      Extra plots to check the status of the procedure

    Returns:
    -------
    """

    # To use the legendre polynomial pixels and orders
    # need to be normalized in the -1,+1 range
    # Normalize pixels
    mnx = np.min(all_pix)
    mxx = np.max(all_pix)
    nrmp = np.array([0.5 * (mnx + mxx), mxx - mnx])
    pix_nrm = 2. * (all_pix - nrmp[0])/nrmp[1]
    # Normalize orders
    mnx = np.min(t)
    mxx = np.max(t)
    nrmt = np.array([0.5 * (mnx + mxx), mxx - mnx])
    t_nrm = 2. * (t - nrmt[0])/nrmt[1]

    if debug:
        # set some plotting parameters
        prettyplot()
        plt.figure(figsize=(7,5))
        msgs.info("Plot identified lines")
        cm = plt.cm.get_cmap('RdYlBu_r')
        sc = plt.scatter(t_nrm, pix_nrm,
                         c=all_wv/10000., cmap=cm)
        cbar = plt.colorbar(sc)
        cbar.set_label(r'Wavelength [$\mu$m]', rotation=270,
                       labelpad=20)
        plt.xlabel(r'Normalized Orders')
        plt.ylabel(r'Normalized Pixels')
        plt.title(r'Location of the identified lines')
        plt.show()

    msgs.info("First iteration")
    # all lines have the same weight
    invvar = np.ones(len(all_wv), dtype=np.float64)
    all_wv_order = all_wv * t
    work2d = np.zeros((nycoeff*nocoeff, len(all_wv)), dtype=np.float64)
    worky = pydl.flegendre(pix_nrm, nycoeff)
    workt = pydl.flegendre(t_nrm, nocoeff)
    for i in range(nocoeff):
        for j in range(nycoeff):
            work2d[j*nocoeff+i,:] = worky[j,:] * workt[i,:]
    work2di = np.transpose(work2d * np.outer(np.ones(nocoeff*nycoeff,
                                             dtype=np.float64),
                                             invvar))
    alpha = work2d.dot(work2di)
    beta = all_wv_order.dot(work2di)
    res = np.linalg.solve(alpha,beta)
    wv_mod = res.dot(work2d)
    if debug:
        # set some plotting parameters
        prettyplot()
        plt.figure(figsize=(7,5))
        plt.axhline(y=np.average(wv_mod / t - all_wv),
                    color='r', linestyle='--')
        plt.axhline(y=+np.std(wv_mod / t - all_wv),
                    color='r', linestyle=':')
        plt.axhline(y=-np.std(wv_mod / t - all_wv),
                    color='r', linestyle=':')
        plt.scatter(all_wv/10000.,
                    wv_mod / t - all_wv,
                    marker="v")
        plt.text(np.min(all_wv/10000), np.average(wv_mod/t-all_wv),
                 r'Average={0:.1f}$\AA$'.format(np.average(wv_mod/t-all_wv)),
                 ha="left", va="bottom",
                 bbox=dict(boxstyle="square",
                           ec=(1., 0.5, 0.5),
                           fc=(1., 0.8, 0.8),
                           alpha=0.7,
                           ))
        plt.text(np.max(all_wv/10000), np.std(wv_mod/t-all_wv),
                 r'Sigma={0:.1f}$\AA$'.format(np.std(wv_mod/t-all_wv)),
                 ha="right", va="bottom",
                 bbox=dict(boxstyle="square",
                           ec=(1., 0.5, 0.5),
                           fc=(1., 0.8, 0.8),
                           alpha=0.7,
                           ))
        plt.title(r'Residuals after 1st iteration')
        plt.xlabel(r'Wavelength [$\mu$m]')
        plt.ylabel(r'Residuals [$\AA$]')
        plt.show()

    msgs.info("Second iteration")
    # Mask Values
    # msk = True means a bad value
    msk = sigma_clip(wv_mod-all_wv_order, sigma=sigmarjct, cenfunc=np.ma.mean).mask
    if np.any(msk):
        msgs.info("Rejecting: {} of {} lines.".format(len(msk[np.where(msk == True)]),len(msk)))
        invvar[msk] = 0.
        work2di = np.transpose(work2d * np.outer(np.ones(nocoeff*nycoeff,
                                                 dtype=np.float64),
                                                 invvar))
        alpha = work2d.dot(work2di)
        beta = all_wv_order.dot(work2di)
        res = np.linalg.solve(alpha,beta)
        wv_mod = res.dot(work2d)
        if debug:
            prettyplot()
            plt.figure(figsize=(7,5))
            plt.axhline(y=np.average(wv_mod[~msk] / t[~msk] - all_wv[~msk]),
                        color='r', linestyle='--')
            plt.axhline(y=+np.std(wv_mod[~msk] / t[~msk] - all_wv[~msk]),
                        color='r', linestyle=':')
            plt.axhline(y=-np.std(wv_mod[~msk] / t[~msk] - all_wv[~msk]),
                        color='r', linestyle=':')
            plt.scatter(all_wv[msk]/10000.,
                        wv_mod[msk] / t[msk] - all_wv[msk],
                        marker="v",
                        label=r'Rejected values')
            plt.scatter(all_wv[~msk]/10000.,
                        wv_mod[~msk] / t[~msk] - all_wv[~msk],
                        marker="v",
                        label=r'Good values')
            plt.text(np.min(all_wv/10000), np.average(wv_mod[~msk]/t[~msk]-all_wv[~msk]),
                     r'Average={0:.1f}$\AA$'.format(np.average(wv_mod[~msk]/t[~msk]-all_wv[~msk])),
                     ha="left", va="bottom",
                     bbox=dict(boxstyle="square",
                               ec=(1., 0.5, 0.5),
                               fc=(1., 0.8, 0.8),
                               alpha=0.7,
                               ))
            plt.text(np.max(all_wv/10000), np.std(wv_mod[~msk]/t[~msk]-all_wv[~msk]),
                     r'Sigma={0:.1f}$\AA$'.format(np.std(wv_mod[~msk]/t[~msk]-all_wv[~msk])),
                     ha="right", va="bottom",
                     bbox=dict(boxstyle="square",
                               ec=(1., 0.5, 0.5),
                               fc=(1., 0.8, 0.8),
                               alpha=0.7,
                               ))
            plt.legend()
            plt.title(r'Residuals after 2nd iteration')
            plt.xlabel(r'Wavelength [$\mu$m]')
            plt.ylabel(r'Residuals [$\AA$]')
            plt.show()
    else:
        msgs.info("No line rejected")

    # Check quality
    gd_wv = invvar > 0.
    resid = (wv_mod[gd_wv]-all_wv_order[gd_wv])
    fin_rms = np.sqrt(np.mean(resid**2))
    msgs.info("RMS: {0:.5f} Ang*Order#".format(fin_rms))

    # Plot QA

    all_pix_qa = np.arange(np.min(all_pix),np.max(all_pix),1)
    pix_nrm_qa = 2. * (all_pix_qa - nrmp[0])/nrmp[1]
    worky_qa = pydl.flegendre(pix_nrm_qa, nycoeff)
    mn, mx = np.min(wv_mod/t), np.max(wv_mod/t)
    order = np.arange(np.min(t),np.max(t)+1,1)

    prettyplot()

    plt.figure(figsize=(7,5))
    plt.title(r'Arc 2D FIT, nx={0:.0f}, ny={1:.0f}, RMS={2:.5f} Ang*Order#'.format(nocoeff, nycoeff,fin_rms))
    plt.xlabel(r'Wavelength [$\AA$]')
    plt.ylabel(r'Row [pixel]')

    for ii in order:
        # define the color
        rr = (ii-np.max(order))/(np.min(order)-np.max(order))
        gg = 0.0
        bb = (ii-np.min(order))/(np.max(order)-np.min(order))
        tsub = np.ones_like(len(all_pix_qa),dtype=np.float64) * ii
        t_nrm_qa = 2. * (tsub - nrmt[0])/nrmt[1]
        work2d_qa = np.zeros((nycoeff*nocoeff, len(all_pix_qa)), dtype=np.float64)
        workt_qa = pydl.flegendre(t_nrm_qa, nocoeff)
        for i in range(nocoeff):
            for j in range(nycoeff):
                work2d_qa[j*nocoeff+i,:] = worky_qa[j,:] * workt_qa[i,:]
        wv_mod_qa = res.dot(work2d_qa)
        plt.plot(wv_mod_qa/ii, all_pix_qa,
                 color=(rr,gg,bb), linestyle='-')
        # Residuals
        print(t[t == ii])
        resid_qa = (wv_mod[t == ii]-all_wv_order[t == ii])/t[t == ii]
        plt.scatter(wv_mod[t == ii]/t[t == ii]+100*resid_qa, all_pix[t == ii],
                    color=(rr,gg,bb))
    plt.text(mx,np.max(all_pix),
             r'residuals $\times$100',
             ha="right", va="top",)
    plt.show()


'''

       

           wv = dblarr(npix)
           wv[*] = work2d # out_str.res / ordr_str[jj].order
           
           mn = min(wv, max=mx)
           
           oplot, wv, all_pix, color=clr.blue
           
           ;; Resid
           pts = where(t EQ ii, npts)
           if npts NE 0 then begin
               nlin = sv_lines[jj].nlin
               sres = (wv_mod[pts] - all_wv[pts])/ordr_str[jj].order
               oplot, wv_mod[pts]/t[pts] + sres*500., $
                 sv_lines[jj].pix[0:nlin-1],  $
                 psym=1, color=clr.black
               ;; RMS
           
           ;; Label
           xyouts, 0.5, 0.96, 'Arc 2D FIT (Res x500) nx='+strtrim(nocoeff,2)+ $
       endfor
           
   ;;;;;;;;;;
       ;; Individual plots
       nordr = n_elements(ordr_str)
       !p.multi = [0, 3, 2]
       pixrms = fltarr(nordr)
       for jj=0L,nordr-1 do begin
           ;; NORMALIZE ORDER
           ii = ordr_str[jj].order
           tsub = replicate(float(ii), npix)
           t_nrm = 2. * (tsub - nrmt[0])/nrmt[1]
           
           ;; work2d and wv
           work2d = dblarr(npix,nycoeff*nocoeff)
           workt = flegendre(t_nrm[*], nocoeff)
           
           for i=0,nocoeff-1 do begin
               for j=0,nycoeff-1 do begin
                   work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
               endfor
           endfor
           
           wv = dblarr(npix)
;           wv[*] = work2d # out_str.res
           wv[*] = work2d # out_str.res / tsub
           
;           mn = min(10^wv, max=mx)
;           plot, all_pix, 10^wv, color=clr.black, $
           mn = min(wv, max=mx)
           plot, all_pix, wv, color=clr.black, $
             background=clr.white, charsize=1.5, yrange=[mn, mx], $
             xrange=[0.,sz[1]], xstyle=1, ystyle=1, xtitle='Row', $
             ytitle='Wavelength', xmargin=[11,2], ymargin=[5,1], /nodata
           
           ;; Fit
           oplot, all_pix, wv, color=clr.blue
;           oplot, all_pix, 10^wv, color=clr.blue
           
           ;; Resid
           pts = where(t EQ ii, npts)
           rms = 9.99
           if npts NE 0 then begin
               nlin = sv_lines[jj].nlin
;               sres = 10^wv_mod[pts] - all_wv[pts] 
;               oplot, sv_lines[jj].pix[0:nlin-1],  10^wv_mod[pts] + sres*100., $
               sres = (wv_mod[pts] - all_wv[pts])/t[pts]
               oplot, sv_lines[jj].pix[0:nlin-1],  wv_mod[pts]/t[pts] $
                 + sres*100., $
                 psym=1, color=clr.black
               ;; Rej
               rej = where(invvar[pts] LE 0., nrej, ncomplement=nnorej, $
                           complement=norej)
               if nrej NE 0 then $
                 oplot, [(sv_lines[jj].pix[0:nlin-1])[rej]],  $
                        [(all_wv[pts]/t[pts] + sres*100)[rej]], $
                        psym=2, color=clr.red
;                       [(10^wv_mod[pts] + sres*100.)[rej]], $
               ;; RMS
               if nnorej NE 0 then $
                 rms = sqrt( total( sres[norej]^2 ) / float(nnorej-1)) $
               else rms = 9.99
           endif
           
           ;; Label
;           if side EQ 2 then ylbl = mn + (mx-mn)*0.08*(findgen(5)+1) $
;           else ylbl = mx - (mx-mn)*0.08*(findgen(5)+1) 
           ylbl = mx - (mx-mn)*0.08*(findgen(5)+1) 
           
           xyouts, sz[1]*0.05, ylbl[0], 'Order = '+strtrim(ii,2), $
             color=clr.black, charsize=1.5
           dwv = abs(wv[0]-wv[npix-1])/float(all_pix[npix-1])
;           dwv = abs(10^wv[0]-10^wv[npix-1])/float(all_pix[npix-1])
           xyouts, sz[1]*0.05, ylbl[2], '!9Dl!X = '+string(dwv,format='(f6.4)'), $
             color=clr.black, charsize=1.5
           xyouts, sz[1]*0.05, ylbl[1], 'RMS(pix) = '+$
             string(rms/dwv,format='(f4.2)'), color=clr.black, charsize=1.5
           pixrms[jj] = rms/dwv
       endfor
           
       x_psclose
       !p.multi=[0,1,1]

       replace_title = '"' + '%%Title: '+qafil + ' ' +systime() + '"'
       ps_replacetitle, replace_title, qafil

       spawn, 'gzip -f '+qafil
   endif

   return, 'Success'

end

'''
