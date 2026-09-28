pro thruput

  if not keyword_set( PSFILE ) then psfile = 'thruput.ps'

  if not keyword_set( LSZ ) then lsz = 1.
  c = x_constants()

  ;; Input
  strcfil = ['/raid/Keck/ESI/data/2000-Apr-07/esi_2000apr07.fits', $
             '/raid/Keck/ESI/data/2005-Aug-04/esi_2005aug04.fits' $
            ]
  frame = [52L,84L]

  nplt = n_elements(frame)

  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
  clrs = x_setclrs(nplt)
  !p.multi=[0,1,1]

  ;; Overall
  yrng = [16.5, 20]
  plot, [0.], [0.], color=clr.black, $
    background=clr.white, charsize=2.2,$
    xmargin=[7,1], ymargin=[4,1], xtitle='Wavelength', $
    ytitle='Zero point (AB Mag)', /nodata, xthick=5, ythick=5, xstyle=1, ystyle=1, $
    xr=[4000., 1e4], yr=yrng

  ;; Loop
  for qq=0L,nplt-1 do begin
      ;; ESI structure
      esi = esi_ar(strcfil[qq])
      idx = where(esi.frame EQ frame[qq])

      ;; Obj structure
      pos = strpos(strcfil[qq],'esi_')
      cf = ''
      if frame[qq] LT 10 then cf = '00'
      if frame[qq] GE 10 and frame[qq] LT 100 then cf = '0'
      objfil = strmid(strcfil[qq],0,pos)+'Extract/Obj_esi0'+cf+ $
               strtrim(frame[qq],2)+'.fits'
      obj = xmrdfits(objfil,1,/sile)

      ;; Sensitivity function
      sensfil = strmid(strcfil[qq],0,pos)+'Extract/sens_esi0'+cf+ $
                strtrim(frame[qq],2)+'.idl'
      restore, sensfil

      for jj=0L,9 do begin
          apx = where(obj[jj].box_var GT 0.,na)
          gpx = apx[0] + lindgen(apx[na-1]-apx[0]+1)
          
          sens = x_calcfit(obj[jj].wave[gpx], fitstr=tot_fit[jj])
          ;; 
          dwv = obj[jj].wave[gpx] - shift(obj[jj].wave[gpx],1) ;; Ang
          dwv[0] = dwv[1]
          
          ;; Zero point (flamb)
;          flam = 1. / sens * dwv / esi[idx].exp * 1e-16
          flam = 1. / sens * dwv * 1e-16
          fnu = flam / c.c * obj[jj].wave[gpx]^2 * 1e-8 ; cgs
          mag = -2.5*alog10(fnu) - 48.6
          print, median(mag)

          ;; Plot
          oplot, obj[jj].wave[gpx], mag, color=clrs[qq]

      endfor
      ;; Label
      xyouts, 8000., yrng[1]-0.2*(qq+1), strmid(strcfil[qq],pos)+' '+$
              string(esi[idx].slit,format='(f4.2)')+' '+esi[idx].Obj, $
              color=clrs[qq], charsiz=lsz
  endfor

  ;; Label
;  xyouts

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  return
end
