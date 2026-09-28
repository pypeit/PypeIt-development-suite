; NAME:
; mkhtml_specthru
;    Version 1.0
;
; PURPOSE:
;    Build a web page (and figures) summarizing a set of spectroscopic
;    thruput data
;
; CALLING SEQUENCE:
;  mkhtml_specthru, thruput_files, TITLE=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /MAUNAKEA -- Use values for MK
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Jan-2008 Written by JXP
;-
;------------------------------------------------------------------------------
pro mkhtml_specthru_mktimefig, thru_files, wvtime, fig_timeeff, OUTPTH=outpth

  if not keyword_set( CSIZE ) then csize = 2.
  if not keyword_set( LSZ ) then lsz = 2.5
  if not keyword_set( OUTPTH ) then outpth = ''

  nwav = n_elements(wvtime)
  nfil = n_elements(thru_files)

  all_eff = fltarr(nwav, nfil)
  xdate = fltarr(nfil)

  i2000 = x_setjdate('01JAN2000')
  ;; Loop for values
  for jj=0L,nfil-1 do begin
      ;; Read in 
      meta = xmrdfits(thru_files[jj],1)
      thru_str = xmrdfits(thru_files[jj],2)

      ;; Evaluate
      linterp, thru_str.wav, thru_str.eff, wvtime, ans
      all_eff[*,jj] = ans

      ;; Date
      jd = x_setjdate(meta.date)
      xdate[jj] = 2000. + (jd - i2000) / 365 + ((jd-i2000) MOD 365) / 365.
  endfor


  xrng = [min(xdate, max=xmax)-0.1, xmax+0.1]
  yrng = [0., max(all_eff)*1.05]

  ;; Smooth
  npt = 80L
  xstp = 0.25
  xsmth = 2000. + xstp * findgen(npt)
  ysmth = fltarr(nwav,npt)
  for ss=0L,npt-2 do begin
      gd = where(xdate GT xsmth[ss] and xdate LT xsmth[ss+1], ngd)
      case ngd of 
          0: 
          1: ysmth[*,ss] = all_eff[*,gd]
          else: ysmth[*,ss] = djs_median( all_eff[*,gd], 2)
      endcase
  endfor

  ;; Plot
  for xx=0,1 do begin
      if xx EQ 0 then begin
          x_psopen, OUTPTH+fig_timeeff+'.ps', /maxs
          !p.multi=[0,1,1]
      endif

      clr = getcolor(/load)
      xclr = x_setclrs(nwav)
      plot, [0], [0], xrange=xrng, $
            yrange=yrng, xtickn=xspaces, xmargin=[8,1], ymargin=[4,1], $
            charsize=csize, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1, xthick=5, ythick=5, /nodata, $
            xtitle='Epoch (yr)', $
            ytitle='Efficiency'
      
      ;; Individual
      for ss=0L,nwav-1 do oplot, xdate, all_eff[ss,*], color=xclr[ss], psym=1
      
      ;; Smoothed
      gsm = where(ysmth[0,*] GT 0.)
      for ss=0L,nwav-1 do oplot, xsmth[gsm]+xstp/2, ysmth[ss,gsm], color=xclr[ss], $
        thick=4
      
      ;; Label
      for ss=0L,nwav-1 do $
             xyouts, xrng[0] + (xrng[1]-xrng[0])*(0.1*(ss+1)), $
                     0.01, string(round(wvtime[ss]),format='(i5)'), $
                     color=xclr[ss], charsiz=1.6
      
      
      if xx EQ 0 then begin
          x_psclose
          spawn, 'ps2pdf '+OUTPTH+fig_timeeff+'.ps '+OUTPTH+fig_timeeff+'.pdf'
          spawn, '\rm '+OUTPTH+fig_timeeff+'.ps'
      endif else begin
          img = tvrd(/true)
          write_jpeg, OUTPTH+fig_timeeff+'.jpg', img, true=1
      endelse
  endfor
  return
end

pro mkhtml_specthru_mkindivfigs, meta, thru_str, $
  fig_zppix, fig_zpang, fig_eff, LAST=last, OUTPTH=outpth


  if not keyword_set( CSIZE ) then csize = 2.
  if not keyword_set( LSZ ) then lsz = 2.5
  if not keyword_set( OUTPTH ) then outpth = ''

  ;; ZP Pixel
  x_psopen, OUTPTH+fig_zppix+'.ps', /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  xrng = [min(thru_str.wav, max=xmax), xmax]
  yrng = [min(thru_str.zp_pix, max=ymax), ymax*1.05]
;  ymax = max(thru_str.zp_pix)
;  medy = median(thru_str.zp_pix)
;  diff = ymax - medy
;  yrng = [ymax-2*diff, ymax]

  plot, [0], [0], xrange=xrng, $
        yrange=yrng, xtickn=xspaces, xmargin=[8,1], ymargin=[4,1], $
        charsize=csize, background=clr.white, color=clr.black, $
        xstyle=1, ystyle=1, xthick=5, ythick=5, /nodata, $
        xtitle='Wavelength (Ang)', $
        ytitle='Zero Point (AB Mag) for 1e!u-!N/s/pix'
  oplot, thru_str.wav, thru_str.zp_pix, color=clr.black
  xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.18, $
          meta.std_name, color=clr.black, charsiz=lsz
  xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.1, $
          meta.date, color=clr.black, charsiz=lsz
  xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.02, $
          'AM='+string(meta.airmass,format='(f4.2)')+'   Bin='+$
          string(meta.spec_bin,format='(i2)')+' Grating: '+$
          string(meta.grating,format='(a12)'), $
          color=clr.black, charsiz=lsz
  x_psclose
  spawn, 'ps2pdf '+OUTPTH+fig_zppix+'.ps '+OUTPTH+fig_zppix+'.pdf'
  spawn, '\rm '+OUTPTH+fig_zppix+'.ps'

  if keyword_set(LAST) then begin
      clr = getcolor(/load)
      plot, [0], [0], xrange=xrng, $
            yrange=yrng, xtickn=xspaces, xmargin=[8,1], ymargin=[4,1], $
            charsize=1.5, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1, thick=3, xthick=3, ythick=3, /nodata, $
            xtitle='Wavelength (Ang)', $
            ytitle='Zero Point (AB Mag) for 1e!u-!N/s/pix'
      oplot, thru_str.wav, thru_str.zp_pix, color=clr.black, thick=3
      xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.18, $
              meta.std_name, color=clr.black, charsiz=lsz
      xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.1, $
              meta.date, color=clr.black, charsiz=lsz
      xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.02, $
              'AM='+string(meta.airmass,format='(f4.2)')+'   Bin='+$
              string(meta.spec_bin,format='(i2)')+' Grating: '+$
              string(meta.grating,format='(a12)'), $
              color=clr.black, charsiz=lsz
      img = tvrd(/true)
      write_jpeg, OUTPTH+fig_zppix+'.jpg', img, true=1
  endif

  ;; ZP Ang
  x_psopen, OUTPTH+fig_zpang+'.ps', /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  xrng = [min(thru_str.wav, max=xmax), xmax]
  yrng = [min(thru_str.zp_ang, max=ymax), ymax*1.05]
;  ymax = max(thru_str.zp_ang)
;  medy = median(thru_str.zp_ang)
;  diff = ymax - medy
;  yrng = [ymax-2*diff, ymax]

  plot, [0], [0], xrange=xrng, $
        yrange=yrng, xtickn=xspaces, xmargin=[8,1], ymargin=[4,1], $
        charsize=csize, background=clr.white, color=clr.black, $
        xstyle=1, ystyle=1, xthick=5, ythick=5, /nodata, $
        xtitle='Wavelength (Ang)', $
        ytitle='Zero Point (AB Mag) for 1e!u-!N/s/Ang'
  oplot, thru_str.wav, thru_str.zp_ang, color=clr.black
  xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.18, $
          meta.std_name, color=clr.black, charsiz=lsz
  xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.1, $
          meta.date, color=clr.black, charsiz=lsz
  xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.02, $
          'AM='+string(meta.airmass,format='(f4.2)')+'  Grating: '+$
          string(meta.grating,format='(a12)'), $
          color=clr.black, charsiz=lsz
  x_psclose
  spawn, 'ps2pdf '+OUTPTH+fig_zpang+'.ps '+OUTPTH+fig_zpang+'.pdf'
  spawn, '\rm '+OUTPTH+fig_zpang+'.ps'

  ;; Efficiency
  x_psopen, OUTPTH+fig_eff+'.ps', /maxs
  !p.multi=[0,1,1]
  clr = getcolor(/load)

  xrng = [min(thru_str.wav, max=xmax), xmax]
  yrng = [min(thru_str.eff, max=ymax), ymax*1.05]
;  ymax = max(thru_str.eff)
;  medy = median(thru_str.eff)
;  diff = ymax - medy
;  yrng = [ymax-2*diff, ymax]

  plot, [0], [0], xrange=xrng, $
        yrange=yrng, xtickn=xspaces, xmargin=[8,1], ymargin=[4,1], $
        charsize=csize, background=clr.white, color=clr.black, $
        xstyle=1, ystyle=1, xthick=5, ythick=5, /nodata, $
        xtitle='Wavelength (Ang)', $
        ytitle='End-to-end Efficiency'
  oplot, thru_str.wav, thru_str.eff, color=clr.black
  x_psclose
  spawn, 'ps2pdf '+OUTPTH+fig_eff+'.ps '+OUTPTH+fig_eff+'.pdf'
  spawn, '\rm '+OUTPTH+fig_eff+'.ps'
  if keyword_set(LAST) then begin
      clr = getcolor(/load)
      plot, [0], [0], xrange=xrng, $
            yrange=yrng, xtickn=xspaces, xmargin=[8,1], ymargin=[4,1], $
            charsize=csize, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1, xthick=3, ythick=3, /nodata, $
            xtitle='Wavelength (Ang)', thick=3, $
            ytitle='End-to-end Efficiency'
      oplot, thru_str.wav, thru_str.eff, color=clr.black, thick=3
      xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.18, $
              meta.std_name, color=clr.black, charsiz=lsz
      xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.1, $
              meta.date, color=clr.black, charsiz=lsz
      xyouts, xrng[0]+(xrng[1]-xrng[0])*0.2, yrng[0]+(yrng[1]-yrng[0])*0.02, $
              'AM='+string(meta.airmass,format='(f4.2)')+'  Grating: '+$
              string(meta.grating,format='(a12)'), $
              color=clr.black, charsiz=lsz
      img = tvrd(/true)
      write_jpeg, OUTPTH+fig_eff+'.jpg', img, true=1
  endif

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mkhtml_specthru, thru_files, OUTFIL=outfil, TITLE=title, GRATING=grating, $
                     WVTIME=wvtime, OUTPTH=outpth, OUTDIR=outdir, LICK=lick, $
                     ALL_GRATING=all_grating
                     
                     
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mkhtml_specthru, thru_files, OUTFIL=, OUTPTH=, /ALL_GRATING [v1.0]'
      return
  endif 


  if not keyword_set( OUTPTH ) then OUTPTH = ''
  if not keyword_set( OUTDIR ) then OUTDIR = 'THRU_DIR/'
  a = findfile(OUTPTH+OUTDIR+'/..', count=count)
  if count EQ 0 then file_mkdir, OUTPTH
  if not keyword_set( OUTFIL ) then outfil = OUTPTH+'thru_summ.html'
  if not keyword_set( TITLE ) then title = 'Thruput Summary'

  nfil = n_elements(thru_files)

  ;; Sort by date
  if nfil GT 0 then begin
      all_dates = dblarr(nfil)
      all_grating = strarr(nfil)
      for jj=0L,nfil-1 do begin
          meta = xmrdfits(thru_files[jj],1)
          ;; Convert to Julian
          all_dates[jj] = x_setjdate(meta.date)
          all_grating[jj] = strtrim(meta.grating,2)
          print, meta.date, all_dates[jj]
      endfor
      ;; Check on unique gratings
      uni = uniq(all_grating)
      if n_elements(uni) GT 1 then begin
          if not keyword_set(ALL_GRATING) then begin
              print, 'mkhtml_specthru: More than 1 grating not allowed'
              print, 'mkhtml_specthru: Rerun and specify the grating you want'
              return
          endif
          ;; DEAL
          stop
      endif
      ;; Sort
      srt = sort(all_dates)
      thru_files = thru_files[srt]
  endif
  
  ;; Tables
  close, /all
  openw, 1, outfil
  
  ;; Header
  printf, 1,  ' <html><head><meta name="Generator" content="Manual">'
  printf, 1, '<title>'+title+'</title><style>'
  printf, 1, 'body'
  printf, 1, '{'
  printf, 1, 'background-color:#FFFFEE;'
  printf, 1, 'font-family:sans-serif;'
  printf, 1, 'margin-top:0;'
  printf, 1, '}'
  printf, 1, 'a:hover { color:#FF0000; font-weight:bolder }'
  printf, 1, '</style></head>'
  printf, 1, '<body bgcolor="#ffffee">'
  printf, 1, '<basefont size="3" face="ARIAL"> <br><br><br>'

  printf, 1, '<table width="800" border="0" cellspacing="0" cellpadding="0">'
  printf, 1, '<tr><td align="center">'
  printf, 1, '<h1>'+title+'</h1></tr></td>'
  printf, 1, '</table>'

  printf, 1, '<hr>'


  ;; Table begin
  printf, 1, '<h2>'
  printf, 1, '<table style="font-size: 18px; line-height: 100%;" width="800"'
  printf, 1, 'border="1" cellpadding="1" cellspacing="1"><tbody>' + $
    '<tr bgcolor="#ddddff">'
  printf, 1, '<td align="center">Summary</td>'
  printf, 1, '</tr></tbody></table>'
  printf, 1, '<table style="font-size: 12px; line-height: 100%;" width="800"'
  printf, 1, 'border="1" cellpadding="1" cellspacing="1"><tbody>' + $
    '<tr bgcolor="#ddddff">'
  printf, 1, '<th width="12%">Date</th>'
  printf, 1, '<th width="12%">Standard</th>'
  printf, 1, '<th width="6%">RA (J2000)</th>'
  printf, 1, '<th width="6%">DEC (J2000)</th>'
  printf, 1, '<th width="12%">Grating</th>'
  printf, 1, '<th width="12%">Blocking</th>'
  printf, 1, '<th width="12%">Cent. Wave</th>'
  printf, 1, '<th width="12%">Slit</th>'
  printf, 1, '<th width="12%">Airmass</th>'
  printf, 1, '<th width="32%">Conditions</th>'
  printf, 1, '<th width="10%">Figures</th>'
  printf, 1, '<th width="10%">Tables</th>'

  for q=0L,nfil-1 do begin

      ;; Root
      istrt = strpos(thru_files[q], '/', /reverse_search) > 0L
      iend = strpos(thru_files[q], '.fits')
      root = strmid(thru_files[q], istrt+1, iend-istrt-1)

      ;; Read in the META and THRUPUT data
      meta = xmrdfits(thru_files[q],1)
      meta.grating = strtrim(meta.grating,2)
      thru_str = xmrdfits(thru_files[q],2)

      ;; Color
      if (q MOD 4) EQ 0 then $
        printf, 1, '<tr bgcolor="Yellow">' $
      else printf, 1, '<tr bgcolor="White">'

      ;; Date
      printf, 1, '<td align="center">'+meta.date+'</td>'

      ;; Name
      printf, 1, '<td align="center">'+meta.std_name+'</td>'

      ;; RA DEC
      printf, 1, '<td align="center">'+meta.ra+'</td>'
      printf, 1, '<td align="center">'+meta.dec+'</td>'

      ;; Grating
      printf, 1, '<td align="center">'+strtrim(meta.grating,2)+'</td>'

      ;; Blocking filter
      printf, 1, '<td align="center">'+meta.blocking+'</td>'

      ;; Central Wavelength
      printf, 1, '<td align="center">'+string(round(meta.central_wave),format='(i5)')+'</td>'

      ;; Slit
      printf, 1, '<td align="center">'+string(meta.slit_width,format='(f8.2)')+ $
              '</td>'

      ;; Airmass
      printf, 1, '<td align="center">'+string(meta.airmass,format='(f4.2)')+'</td>'

      ;; Conditions
      printf, 1, '<td align="center">'+meta.conditions+'</td>'


      ;; Figures
      printf, 1, '<td align="left">'
      fig_zppix = OUTDIR+'/ZPPIX_'+root
      fig_zpang = OUTDIR+'/ZPANG_'+root
      fig_eff = OUTDIR+'/EFF_'+root
      mkhtml_specthru_mkindivfigs, meta, thru_str, fig_zppix, $
        fig_zpang, fig_eff, $
        LAST=(q EQ (nfil-1)), OUTPTH=outpth
;      spawn, 'gzip -f '+fig_zppix
;      spawn, 'gzip -f '+fig_zpang
      printf, 1, '<a href="'+fig_zppix+'.pdf">ZP_pix</a> '
      printf, 1, '<a href="'+fig_zpang+'.pdf">ZP_Ang</a> '
      printf, 1, '<a href="'+fig_eff+'.pdf">Efficiency</a> '
      printf, 1, '</td>'
      

      ;; ZP Tables
      printf, 1, '<td align="left">'
      outtab = OUTDIR+'/ZPPIX_'+root+'.dat'
      writecol, OUTPTH+outtab, thru_str.wav, thru_str.zp_pix
      printf, 1, '<a href="'+outtab+'">ZP_pix</a> '
      outtab = OUTDIR+'/ZPANG_'+root+'.dat'
      printf, 1, '<a href="'+outtab+'">ZP_Ang</a> '
      writecol, OUTPTH+outtab, thru_str.wav, thru_str.zp_ang
      outtab = OUTDIR+'/EFF'+root+'.dat'
      printf, 1, '<a href="'+outtab+'">Efficiency</a> '
      writecol, OUTPTH+outtab, thru_str.wav, thru_str.eff
      printf, 1, '</td>'

      ;; Print
      printf, 1, '</tr>'
  endfor

  printf, 1, '</table>'
  printf, 1, '<br>'


  printf, 1, '<table width="800" border="0" cellspacing="0" cellpadding="0">'
  printf, 1, '<tr><td align="left">'
  printf, 1, 'The above table summarizes the set of standard stars that have'
  printf, 1, 'been processed.  They should have all been observed in photometric'
  printf, 1, 'conditions.  You may note, however, that the slitwidth and airmass'
  printf, 1, 'vary consdierably.  The efficiency calculations have made an estimated'
  printf, 1, 'correction for the airmass, but NO correction for slit losses.<br>'
  if keyword_set(LICK) then obsvty = 'Lick' else obsvty = 'Keck'
  printf, 1, 'If you have a standard star to contribute, please contact the '+obsvty
  printf, 1, ' staff or JXP.  <br><br> Here is a summary of the standards observed:'
  printf, 1, '</tr></td>'

  printf, 1, '<tr><td align="left">'
  printf, 1, '<UL>'
  printf, 1, '<LI> ZP_pix:  Zero point (AB magnitude) to give one electron per second'
  printf, 1, 'per pixel.  This does NOT include an airmass correction.'
  printf, 1, '<LI> ZP_Ang:  Zero point (AB magnitude) to give one electron per second'
  printf, 1, 'per Angstrom.  This does NOT include an airmass correction.'
  printf, 1, '<LI> Efficiency:  Total efficiency of the telescope+instrument. This'
  printf, 1, 'does include an airmass correction.'
  printf, 1, '</UL>'
  printf, 1, '</tr></td>'


  printf, 1, '</tr>'
  printf, 1, '</tbody></table><br>'
  printf, 1, '</h2>
  printf, 1, '<hr>'
  printf, 1, '<br>'

  ;; Current ZP and Efficiency
  printf, 1, '<table width="800" border="0" cellspacing="0" cellpadding="0">'
  printf, 1, '<tr><td align="center"><h2>Current ZP and Efficiency Curves</h2></td></tr>'

  ;; Text
  printf, 1, '<tr><td align="left">'
  printf, 1, 'The following figures show the most recent measurements of the '
  printf, 1, 'ZeroPoint and efficiency measurements observed through a '
  printf, 1, 'sufficiently large slit (usually >2 arcsec).'
  printf, 1, '</tr></td>'

  ;; Zeropoint
  printf, 1, '<tr><td align="center"><img src="'+fig_zppix+ $
          '.jpg" width="700" align="center"/> </td></tr>'

  ;; Efficiency
  printf, 1, '<tr><td align="center"><img src="'+fig_eff+ $
          '.jpg" width="700" align="center"/> </td></tr>'

  printf, 1, '</table>'
  printf, 1, '<hr>'

  ;; Time evolution
  if keyword_set(WVTIME) and nfil NE 1 then begin
      ;; Make the figures
      fig_timeeff = OUTDIR+'/time_efficiency'
      mkhtml_specthru_mktimefig, thru_files, wvtime, fig_timeeff, OUTPTH=outpth
      

      ;; Start Table
      printf, 1, '<table width="800" border="0" cellspacing="0" cellpadding="0">'
      printf, 1, '<tr><td align="center"><h2>Time Evolution</h2></td></tr>'

      ;; Text
      printf, 1, '<tr><td align="left">'
      printf, 1, 'The following figure shows the time history of the total efficiency'
      printf, 1, '(telescopes + instrument) evaluated at a few wavelenths as a function'
      printf, 1, 'of time.  Note that this analysis is restricted to wider slits.'
      printf, 1, '</tr></td>'

      printf, 1, '<tr><td align="center"><a href="'+fig_timeeff+ $
              '.pdf"><img src="'+fig_timeeff+ $
              '.jpg" width="700" align="center"/></a></td></tr>'

      ;; End Table
      printf, 1, '</table>'
      printf, 1, '<hr>'
  endif


  ;; Closing
  printf, 1, '<font size="1">'
  printf, 1, 'Last modified: '+systime(/utc)
  printf, 1, '</font>'
  printf, 1, '</body></html>'
  

  close, /all
  ;; All done
  print, 'mk_webtab: All done'
  return
end
  
