;-----------------------------------------------------------------------
pro deimos_throughput_grating_summary_plots, params, summary, $
  OUTDIR=outdir, $
  CLOBBER=clobber, $
  VERBOSE=verbose
;-----------------------------------------------------------------------

;; declarations...
csize = 1.75
template = '/tmp/temp.%.ps'
angstrom = STRING(197B)

;; generate filenames...
root = OUTDIR+'/eff_vs_time'
summary.eff_vs_time_bintab = root + '.fits'
summary.eff_vs_time_pdf    = root + '.pdf'
summary.eff_vs_time_plot   = root + '.png'
summary.eff_vs_time_tab    = root + '.txt'

root = OUTDIR+'/eff_vs_wave'
summary.eff_vs_wavelength_bintab = root + '.fits'
summary.eff_vs_wavelength_pdf    = root + '.pdf'
summary.eff_vs_wavelength_plot   = root + '.png'
summary.eff_vs_wavelength_tab    = root + '.txt'

root = OUTDIR+'/eff_latest'
summary.efficiency_current_bintab = root + '.fits'
summary.efficiency_current_pdf    = root + '.pdf'
summary.efficiency_current_plot   = root + '.png'
summary.efficiency_current_tab    = root + '.txt'

;;----------------------------------------
;; current efficiency...
;;----------------------------------------

;; generate data for current efficiency by combining data from the
;; most recent dataset...
last_date = params[0].date

;; locate other data from this date...
good = where( params.date eq last_date, count)
if count eq 0 then message, 'no data matching date '+last_date

;; get array of cenlam...
cenlam = params[good].cenlam
cenlam_uniq = cenlam[UNIQ( cenlam, SORT(cenlam))]
n_cenlam_uniq = n_elements(cenlam_uniq)

keepers = lonarr(n_cenlam_uniq)
dw_median = fltarr(n_cenlam_uniq)

;; select the last exposure for each cenlam and determine data range...
for i=0,n_cenlam_uniq-1 do begin
    good = where( params.date eq last_date and $
                  params.cenlam eq cenlam_uniq[i], $
                  count)
    if count lt 1 then message, 'no data matching date and cenlam'
    j = good[0]
    keepers[i] = j

    ;; Read in the THRUPUT data
    thru_str = xmrdfits( params[j].infile, 2, /silent)

    ;; determine min/max values...
    xminval = min(thru_str.wav, max=xmaxval)
    yminval = min(thru_str.eff, max=ymaxval)

    if i eq 0 then begin
        wav_min = xminval
        eff_min = yminval
        wav_max = xmaxval
        eff_max = ymaxval
    endif else begin
        wav_min = min( [wav_min, xminval])
        wav_max = max( [wav_max, xmaxval])
        eff_min = min( [eff_min, yminval])
        eff_max = max( [eff_max, ymaxval])
    endelse 

    ;; determine pixel spacing...
    npix = n_elements(thru_str.wav)
    dw_temp = thru_str.wav[1:npix-1] - thru_str.wav[0:npix-2]
    dw_median[i] = median(dw_temp)

endfor

;; create a new grid for the overall throughput data...
dw = 10.*avg(dw_median)
dw2 = dw/2.
npix = fix(float(wav_max - wav_min)/dw)
wav = wav_min + findgen(npix)*dw + dw2

;; allocate array to hold latest re-interpolated eff curves...
interp = fltarr(npix,n_cenlam_uniq)

;; re-interpolate data onto grid...
for i=0,n_cenlam_uniq-1 do begin
    j = keepers[i]
    thru_str = xmrdfits( params[j].infile, 2, /silent)
    eff = INTERPOL( thru_str.eff, thru_str.wav, wav)

    ;; screen out-of-bounds pixels...
    w1 = min(thru_str.wav, max=w2)
    bad = where( wav lt w1 or wav gt w2, count)
    if count gt 0 then eff[bad] = !values.f_nan

    ;; insert data into array...
    interp[0,i] = eff
endfor

;; now compute the maximum at each data point...
if n_cenlam_uniq gt 1 then eff = max( interp, dim=2, /nan)

;; start plot...
title = 'Latest DEIMOS Efficiency Data for Grating ' $
        +params[0].grating $
        +' ('+last_date+')'
psfile = mktemp(template)
psland, psfile
DEVICE, SET_FONT='Helvetica', /TT_FONT  
DEVICE, /ISOLATIN1
plot, wav, eff, $
      xrange=[wav_min,wav_max], $
      yrange=[eff_min,eff_max], $
      charsize=csize, $
      xstyle=1, ystyle=1, $
      thick=10, xthick=1, ythick=1, $
      xtitle='Wavelength ['+Angstrom+']', $
      ytitle='End-to-end Efficiency', $
      title=title, font=0, $
      xmargin=[2,1], ymargin=[1,1]

;; add lines to plot...
for i=0,n_cenlam_uniq-1 do begin
    j = keepers[i]
    thru_str = xmrdfits( params[j].infile, 2, /silent)
    oplot, thru_str.wav, thru_str.eff, linestyle=1
endfor 

;; complete plot...
id, font=1
device, /close
ps2other, psfile, $
          png=summary.efficiency_current_plot, $
          pdf=summary.efficiency_current_pdf, $
          verbose=verbose, /delete

;;----------------------------------------
;; time evolution...
;;----------------------------------------

;; compute extrema...
time_min = min(params.jd, max=time_max, /nan)
eff_min  = min(params.efficiency, max=eff_max, /nan)

;; define date labels...
dummy = LABEL_DATE(DATE_FORMAT=['%Y']) 

;; initialize plot...
title = 'DEIMOS Efficiency in Various Passbands for Grating ' +params[0].grating 
psfile = mktemp(template)
psland, psfile
plot, [0], [0], $
      xrange=[time_min, time_max], $
      yrange=[eff_min, eff_max], $
      charsize=csize, $
      xstyle=1, ystyle=1, $
      thick=2, xthick=1, ythick=1, $
      xtitle='Time', $
      ytitle='End-to-end Efficiency', $
      title=title, font=0, $
      xmargin=[2,1], ymargin=[1,1], $
      XTICKUNITS = ['Time'], $  
      XTICKFORMAT='LABEL_DATE', XTICKINTERVAL=1

;; loop over passbands...
n_bands = n_elements( params[0].lambda_eff)
for i=0,n_bands-1 do begin

    ;; locate good values...
    good = where( finite(params.efficiency[i]), count)

    ;; skip if fewer than 2 data points...
    if count lt 2 then continue

    ;; select good data...
    x = params[good].jd
    y = params[good].efficiency[i]

    order = sort(x)
    x = x[order]
    y = y[order]

    oplot, x, y, psym=-4, linestyle=i+1

    ;;for j=0,count-1 do begin
    ;;k = [good[j]]
    ;;print, 'jd=', x[k], ' jd2=', params[k].jd, $
    ;;' date=', params[k].date, ' set=', params[k].dataset
    ;;endfor 
    
endfor

;; complete plot...
id, font=1
device, /close
ps2other, psfile, $
          png=summary.eff_vs_time_plot, $
          pdf=summary.eff_vs_time_pdf, $
          verbose=verbose, /delete

end
