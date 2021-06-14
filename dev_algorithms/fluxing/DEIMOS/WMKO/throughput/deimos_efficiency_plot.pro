;-----------------------------------------------------------------------
pro showEvent, date1, date2, label
;-----------------------------------------------------------------------

x1 = xdatestring2jd( date1)
x2 = xdatestring2jd( date2)

y1 = !y.crange[0]
y2 = !y.crange[1]

xbox = [x1,x2,x2,x1]
ybox = [y1,y1,y2,y2]
ymin = 0.5*(y1+y2)
polyfill, xbox, ybox, color=192

space = '  '
xyouts, x1, y2, label+space, align=1.0, orient=90, charsize=1.5

end

;-----------------------------------------------------------------------
pro read_deimos_eff_sav_file, infile, extract, fextract, sensstd, efficiency
;-----------------------------------------------------------------------
extract = 0B
restore, infile
if size( extract, /tname) ne 'STRUCT' then  $
    message, 'WARNING: there is no struct named EXTRACT in file '+file
if size( fextract, /tname) ne 'STRUCT' then  $
    message, 'WARNING: there is no struct named FEXTRACT in file '+file
if size( sensstd, /tname) ne 'STRUCT' then  $
    message, 'WARNING: there is no struct named SENSSTD in file '+file
if size( efficiency, /tname) ne 'STRUCT' then  $
    message, 'WARNING: there is no struct named EFFICIENCY in file '+file

end

;-----------------------------------------------------------------------
pro deimos_efficiency_plot, infiles, psfile, $
  GRANAME=graname, WAVELEN=wavelen, LAMBDA_EFF=lambda_eff, VERBOSE=verbose
;-----------------------------------------------------------------------
; Purpose:
;       Plot efficiency results for combination of grating and
;       wavelength.
;-----------------------------------------------------------------------

;; help, graname, wavelen, lambda_eff
device, true=24
device, decomp=0

;; count files...
n = n_elements(infiles)

;; allocate arrays...
dates = fltarr(n)
eff = fltarr(n)

;; loop over files...
ngood = 0
nbad = 0
for i=0,n-1 do begin

    breakname, infiles[i], dirname, rootname, extn

    ;; read data....
    read_deimos_eff_sav_file, infiles[i], extract, fextract, $
      sensstd, efficiency

    ;; check params...
    if extract.meta.grating ne graname then begin
        if keyword_set(VERBOSE) then $
          print, rootname, ': wrong grating (', extract.meta.grating, ')'
        continue
    endif 

    if keyword_set(WAVELEN) then begin
        if nint(extract.meta.central_wave) ne nint(wavelen) then begin
            if keyword_set(VERBOSE) then $
              print, rootname, ': wrong wavelength (', extract.meta.central_wave, ')'
            continue
        endif
    endif 

    ;; loop over tabulated efficiency values...
    n_passbands = n_elements(efficiency)
    for j=0,n_passbands-1 do begin
        if nint(efficiency[j].lambda_eff) eq nint(lambda_eff) then begin
            if finite(efficiency[j].efficiency) then begin
                if keyword_set(VERBOSE) then $
                  print, rootname, ': good data for date ', extract.meta.date
                dates[ngood] = xdatestring2jd( extract.meta.date)
                eff[ngood] = efficiency[j].efficiency
                ngood += 1
            endif else begin
                if keyword_set(VERBOSE) then $
                  print, rootname, ': NaN data'
                nbad += 1
            endelse 
        endif 
    endfor 
endfor 

help, ngood, nbad

;; check for not enough good data
if ngood lt 2 then begin
    message, 'not enough good data', /info
    return
endif 

;; truncate arrays...
dates = dates[0:ngood-1]
eff   = eff[0:ngood-1]

date1 = xdatestring2jd( '01JAN2002')
date2 = xdatestring2jd( '01JAN2014')

psopen, psfile, /land, /helvetica, /color, $
        xsize=11.0, ysize=8.5, /inches
loadct, 0
;; buf = label_date( date_format='%M!C%Y')
Angstrom = STRING(197B)
if keyword_set(WAVELEN) then begin
    format = '("DEIMOS Grating=",a," Cenlam=",i4,a," Passband=",i4,a)'
    title = string( format=format, graname, nint(wavelen), Angstrom, $
                    nint(lambda_eff), Angstrom)
endif else begin
    format = '("DEIMOS Grating=",a," Passband=",i4,a)'
    title = string( format=format, graname, $
                    nint(lambda_eff), Angstrom)
endelse 
;; plot, dates, eff, XTICKFORMAT = 'LABEL_DATE', psym=6, title=title, $
!p.font = 1
plot, dates, eff, psym=6, title=title, $
      xtitle='Date', ytitle='Efficiency', $
      yrange=[0.,0.3], $
      xrange=[date1,date2], xstyle=1, $
      xticks=6, xtickunits='Years', xminor=4, $
      charsize=2, /nodata
id, psfile

;; add underlay...
showEvent, '22MAY2002', '13NOV2002', 'M1'
showEvent, '04MAY2004', '20OCT2004', 'M1'
showEvent, '14AUG2007', '14SEP2007', 'M2'
showEvent, '01SEP2009', '23SEP2010', 'M1+M2+M3'
showEvent, '05NOV2010', '05DEC2010', 'M3'
showEvent, '08MAY2012', '31OCT2012', 'M1+M2'
showEvent, '25MAR2005', '25APR2005', 'Tent'
showEvent, '10JUL2005', '10AUG2005', 'Coll'

;; overplot...
plot, dates, eff, psym=6, title=title, $
      xtitle='Date', ytitle='Efficiency', $
      yrange=[0.,0.3], $
      xrange=[date1,date2], xstyle=1, $
      xticks=6, xtickunits='Years', xminor=4, $
      charsize=2, /noeras

psclose
print, 'wrote psfile ', psfile

end 
