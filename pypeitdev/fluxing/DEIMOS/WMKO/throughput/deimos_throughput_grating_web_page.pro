;------------------------------------------------------------------------
pro deimos_throughput_grating_web_page, $
  OUTFILE=outfile, $
  TITLE=title, $
  INDEX=index, $
  GRATING=grating, $
  FIG_EFF=fig_eff, $
  FIG_ZP_ANG=fig_zp_ang, $
  FIG_TIME=fig_time
;------------------------------------------------------------------------
; Purpose:
;       Generate web page for a single grating
;------------------------------------------------------------------------

;; parse keywords...
if ~ keyword_set( OUTFILE ) then outfile = OUTPATH+'thru_summ.html'
if ~ keyword_set( TITLE )   then title = 'Thoughput Summary'
if ~ keyword_set( GRATING)  then grating = ''

extn1 = '.png'
extn2 = '.pdf'

;; Current ZP and Efficiency...
body = ''
if keyword_set( FIG_EFF) then begin
    pdf = strsub( fig_eff, extn1, extn2)
    body += '<h2>Current DEIMOS '+grating+' Efficiency</h2>'
    body += '<a href="'+pdf+'">'
    body += '<img src="'+fig_eff+'" alt="DEIMOS '+grating+' efficiency plot">'
    body += '</a>'
endif 

if keyword_set( FIG_ZP_ANG) then begin
    pdf = strsub( fig_zp_ang, extn1, extn2)
    body += '<h2>Current DEIMOS '+grating+' Zero Point</h2>'
    body += '<a href="'+pdf+'">'
    body += '<img src="'+fig_zp_ang+'" alt="DEIMOS '+grating+' zeropoint plot">'
    body += '</a>'
endif 

if keyword_set( FIG_TIME) then begin
    pdf = strsub( fig_time, extn1, extn2)
    body += '<h2>DEIMOS '+grating+' Efficiency vs. Time</h2>'
    body += '<a href="'+pdf+'">'
    body += '<img src="'+fig_time+'" alt="DEIMOS '+grating+' efficiency vs. time">'
    body += '</a>'
endif 

;; write the file...
deimos_write_web_page, outfile, body, title=title, index=index

end

