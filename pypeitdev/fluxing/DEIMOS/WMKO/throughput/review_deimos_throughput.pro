;---------------------------------------------------------------------
pro read_sensstd_sav_file, filename, $
  EXTRACT=extract, $
  FEXTRACT=fextract, $
  SENSSTD=sensstd, $
  STATUS=status
;-----------------------------------------------------------------------
; Given filename, return extract, fextract, sensstd
;-----------------------------------------------------------------------

status = 1B
extract = 0B
fextract = 0B
sensstd = 0B

;; read the input structures (meta, spec, diag0, diag1)...
restore, filename

;; confirm that we got the expected things...
if size( extract, /tname) ne 'STRUCT' then begin
    status = 0B
    message, 'WARNING: there is no struct named EXTRACT in file '+file, /info
endif 
if size( fextract, /tname) ne 'STRUCT' then begin
    status = 0B
    message, 'WARNING: there is no struct named FEXTRACT in file '+file, /info
endif 
if size( sensstd, /tname) ne 'STRUCT' then begin
    status = 0B
    message, 'WARNING: there is no struct named SENSSTD in file '+file, /info
endif 

end

;-----------------------------------------------------------------------
pro display_extract, EXTRACT=extract, FIDUCIAL=fiducial, XPOS=xpos, $
                     WINDOW=window
;-----------------------------------------------------------------------

wset, window

;; scale Y axis from sky to peak with 10% buffer
dy = extract.peak - extract.sky
y_min = extract.sky  - 0.1 * dy
y_max = extract.peak + 0.1 * dy

;; generate axes...
plot, extract.x, extract.flux, $
      xtitle='Column [px]', ytitle='Sky-subtracted Flux', $
      yrange=[y_min,y_max], ystyle=1, title='Extraction Profile'

;; overplot extraction region...
oplot, [extract.x1,extract.x2], [extract.peak,extract.peak], psym=-1

;; mark left sky region...
oplot, [extract.xleft1, extract.xleft2], [extract.sky,extract.sky], psym=-1

;; mark right sky region...
oplot, [extract.xright1, extract.xright2], [extract.sky,extract.sky], psym=-1

;; indicate sky level as dotted horizontal line...
oplot, !x.crange, [extract.sky,extract.sky], linestyle=1

;; plot dotted vertical line at the measured center...
oplot, [extract.center,extract.center], !y.crange, linestyle=1
xyouts, extract.center, extract.sky, $
        ' Center='+stringify(extract.center,'(f15.1)'), $
        orient=90.

;; annotate...
nxyouts, xpos, 0.55, $
         'Row='+stringify(extract.row)
nxyouts, xpos, 0.50, $
         'Sky='+stringify(nint(extract.sky),'(i10)')

;; add fiducial...
if keyword_set(FIDUCIAL) then begin
    oplot, [fiducial.center,fiducial.center], !y.crange, color=gdwcolor('red')
endif 

end

;------------------------------------------------------------------------
pro display_profile, PROFILE=profile, FIDUCIAL=fiducial, XPOS=xpos, $
                     WINDOW=window
;------------------------------------------------------------------------

wset, window

plot, profile.u, profile.flux, $
      xtitle='Column [px]', ytitle='Sky-subtracted Flux', $
      xstyle=1, $
      title='Profile at Rows '+stringify(profile.v1)+ $
      '-'+stringify(profile.v2)
oplot, [profile.u_center,profile.u_center], !y.crange, linestyle=1
xyouts, profile.u_center, !y.crange[0], $
        ' Center='+stringify(profile.u_center,'(f15.1)'), $
        orient=90.
nxyouts, xpos, 0.10, 'Sky='+stringify(nint(profile.sky),'(i10)')
nxyouts, xpos, 0.05, 'S/N='+stringify(profile.s2n,'(f15.1)')

if keyword_set(FIDUCIAL) then begin
    oplot, fiducial.u, fiducial.flux, color=gdwcolor('red')
endif 

end

;------------------------------------------------------------------------
pro display_trace, TRACE=trace, FIDUCIAL=fiducial, PROFILE=profile, $
                   XPOS=xpos, WINDOW=window
;------------------------------------------------------------------------
; display the trace...
;------------------------------------------------------------------------

wset, window
x = [trace.good.v, trace.bad.v]
y = [trace.good.u, trace.bad.u]
x1 = min(x[where(x ne 0)], max=x2)
y1 = min(y[where(y ne 0)], max=y2)

y1f = min(fiducial.fit.u, max=y2f)
ymin = y1 < y1f
ymax = y2 > y2f

plot, [0], [0], /nodata, xtitle='Y [px]', ytitle='X [px]', /yno, $
      xstyle=1, title='Trace', xrange=[x1,x2], yrange=[ymin,ymax]

;; plot good data...
if trace.good.count gt 0 then $
  oplot, trace.good.v, trace.good.u, psym=6, symsize=0.25

;; plot bad data...
if trace.bad.count gt 0 then $
  oplot, trace.bad.v, trace.bad.u, psym=7, symsize=0.5
nxyouts, xpos, 0.55, stringify(trace.bad.count)+' bad points'
nxyouts, xpos, 0.5, 'Sky='+stringify(nint(profile.sky),'(i10)')

;; add trace...
if ( trace.fit.status ) then begin
    oplot, trace.fit.v, trace.fit.u
endif else begin
    nxyouts, 0.5, 0.5, 'WARNING: BAD TRACE', align=0.5
endelse 

;; add fiducial...
if keyword_set(FIDUCIAL) then begin
    if ( fiducial.fit.status ) then begin
        oplot, fiducial.fit.v, fiducial.fit.u, color=gdwcolor('red')
    endif 
endif 

end

;------------------------------------------------------------------------
pro display_spectrum, SPECTRUM=spectrum, FIDUCIAL=fiducial, WINDOW=window
;------------------------------------------------------------------------
; display the spectrum...
;------------------------------------------------------------------------

wset, window

;; get min/max values...
ymin1 = min( spectrum.flux, max=ymax1)
ymin2 = min( fiducial.flux, max=ymax2)
ymin = ymin1 < ymin2
ymax = ymax1 > ymax2
print, 'Spectrum range is ', ymin1, ymax1
print, 'Fiducial range is ', ymin2, ymax2

plot, spectrum.pixel, spectrum.flux, xtitle='Row', ytitle='flux [e-/px]', $
      xstyle=1, title='Raw spectrum', yrange=[ymin, ymax]
if keyword_set(FIDUCIAL) then begin
    message, 'overplotting fiducial', /info
    oplot, fiducial.pixel, fiducial.flux, color=gdwcolor('red')
endif 
end 

;-----------------------------------------------------------------------
pro update_generic_list, filename, v1
;-----------------------------------------------------------------------

;; if file doesn't exist, then just create it...
if ~ file_test(filename) then begin
    openw, ounit, filename, /get_lun
    printf, ounit, v1
    free_lun, ounit

endif else begin

    ;; file already exists, so update it...

    ;; get current contents...
    readcol, filename, format='(a)', c1, /silent
    match = where( c1 eq v1, count)
    
    if count eq 0 then begin
        
        ;; add new value to file...
        openw, ounit, filename, /get_lun, /append
        printf, ounit, v1
        free_lun, ounit
        
    endif else begin
        message, 'WARNING: '+v1+' already appears in file '+filename, /info
    endelse 
endelse 

end

;-----------------------------------------------------------------------
pro update_commented_list, filename, v1, string
;-----------------------------------------------------------------------

;; if file doesn't exist, then just create it...
if ~ file_test(filename) then begin

    write_csv, filename, [v1], [string]

endif else begin

    ;; file already exists, so update it...

    ;; get current contents...
    result = read_csv(filename)
    
    ;; check for empty file...
    if size(result, /n_dim) eq 0 then begin
        write_csv, filename, [v1], [string]
    endif else begin
        
        c1 = result.field1
        c2 = result.field2
        match = where( c1 eq v1, count)
        
        if count eq 0 then begin
            
            ;; add new value to file...
            write_csv, filename, [c1,v1], [c2,string]
            
        endif else begin
            message, 'WARNING: '+v1+' already appears in file '+filename, /info
        endelse 
    endelse 
endelse 

end

;-----------------------------------------------------------------------
pro update_center_list, filename, infile, center
;-----------------------------------------------------------------------

;; if file doesn't exist, then just create it...
if ~ file_test(filename) then begin

    write_csv, filename, [infile], [center]

endif else begin

    ;; file already exists, so update it...

    ;; get current contents...
    result = read_csv(filename)
    c1 = result.field1
    c2 = result.field2
    match = where( c1 eq infile, count)
    
    if count eq 0 then begin
        
        ;; add new value to file...
        c1 = [c1, infile]
        c2 = [c2, center]
        
    endif else if count eq 1 then begin

        ;; replace existing value...
        c1[match] = infile
        c2[match] = center

    endif else begin
        message, 'ERROR: Center file '+filename+' already has mutiple entries for image '+infile, /info
        return
    endelse 

    ;; save results..
    write_csv, filename, c1, c2
endelse 

end

;;------------------------------------------------------------------------
pro view_trace, infile, chip, trace, FIDUCIAL=fiducial
;;------------------------------------------------------------------------
;; display the image...
;;------------------------------------------------------------------------

;; read image...
image = deimos_read_chip( infile, chip, header=header)

;; get rootname...
breakname, infile, dirname, rootname, extn
message, '  Displaying '+rootname+' in ATV...', /info

atv, image
message, '  Overplotting trace points in green', /info
if trace.good.count gt 0 then $
  atvplot, trace.good.u, trace.good.v, color='green', psym=6
if trace.bad.count gt 0 then $
  atvplot, trace.bad.u, trace.bad.v, color='red', psym=7
if ( trace.fit.status ) then begin
    color = 'cyan'
    message, '  Overplotting trace fit in '+color, /info
    atvplot, trace.fit.u, trace.fit.v, color=color
endif else begin
    message, '  Warning: bad trace fit!', /info
endelse

;; display fiducial...
if keyword_set(fiducial) then begin
    atvplot, fiducial.fit.u, fiducial.fit.v, color='yellow'
endif

;; allow time for reflection...
message, '  Please press Q in the ATV window to continue...', /info
atv_activate

end

;-----------------------------------------------------------------------
pro review_deimos_throughput, sav_file, $
  FORCE=force, $
  FIDUCIAL=fiducial, $
  EXIT_STATUS=exit_status
;-----------------------------------------------------------------------
;+
; NAME:
;	REVIEW_DEIMOS_THROUGHPUT
;
; PURPOSE:
;	This procedure is used to review DEIMOS throughput
;	measurements and trigger re-reduction of data as desired.
;
; CATEGORY:
;	Data analysis
;
; CALLING SEQUENCE:
;       REVIEW_DEIMOS_THROUGHPUT, sav_file
;
; INPUTS:
;	sav_file:	name of SAV file generated by deimos_sensstd_gdw
;
; KEYWORDS (INPUT):
;       FORCE: set this keyword to force review of object regardless
;       of current status
;
;       FIDUCIAL: set this keyword to overplot data from the fiducial
;       dataset for this configuration
;
; KEYWORDS (OUTPUT):
;       EXIT_STATUS: returns 0 for normal exit or 1 if the observer
;       selects "QUIT"
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2011-Jun-21	GDW	Original version
;-
; TODO:
;       - fix plotting problem
;       - allow fiducials
;       - add parsing of center list to reductions
;       - add parsing of reject list to reductions
;-----------------------------------------------------------------------

xwinsize = 500
ywinsize = 400

;; define constants...
xpos=0.025
default='A'
null = ''
exit_status = 0

;; define files...
thru_dir = getenv('DEIMOS_THRU_DIR')
if thru_dir eq null then $
  message, 'DEIMOS_THRU_DIR is not defined -- abort!'
data_dir = thru_dir + '/dat/'
center_list = data_dir + 'centers.lst'

;; read review database...
review_database = data_dir + 'review.fits'
review = mrdfits( review_database, 1)

;; check for file...
if ~ file_test( sav_file) then begin
    message, 'WARNING: no save file '+sav_file+' for infile '+infile+'; add to reject list and skip', /info
    reason = 'auto extraction failed'
    update_review_database, review_database, review, infile, status='bad', $
      comment=reason
    return
endif 

;; read the input structures...
read_sensstd_sav_file, sav_file, $
  EXTRACT=extract, $
  FEXTRACT=fextract, $
  SENSSTD=sensstd, $
  STATUS=status
if ~ status then return
infile = extract.meta.filename

;; get the name of the original FITS file...

;; get existing record for this target...
record = parse_review_database( review, infile)

if size( record, /tname) eq 'STRUCT' then begin

    ;; skip objects already considered "bad"...
    if record.status eq 'bad' and ~ keyword_set(FORCE) then begin
        message, 'WARNING: file '+infile+' previously rejected with this comment: '+ $
                 record.comment+'; skip', /info
        return
    endif
    
    ;; skip objects already considered good...
    if record.status eq 'good' and ~ keyword_set(FORCE) then begin
        message, 'WARNING: file '+infile+' previously approved; skip', /info  
        return
    endif

    oldcomment = record.comment
endif else begin
    oldcomment = ''
endelse 

;; define windows...
profile_window = 0
extract_window = 1
trace_window = 2
spec_window = 3

;; create windows...
window, profile_window, xsize=xwinsize, ysize=ywinsize
window, extract_window, xsize=xwinsize, ysize=ywinsize
window, trace_window, xsize=xwinsize, ysize=ywinsize
window, spec_window, xsize=xwinsize, ysize=ywinsize

;; set graphics settings...
!p.font = 1
!p.multi = 0
!p.charsize = 1.5

;; loop over extensions...
for ext=0,1 do begin

    ;; read the input structures...
    read_sensstd_sav_file, sav_file, $
      EXTRACT=extract, $
      FEXTRACT=fextract, $
      SENSSTD=sensstd, $
      STATUS=status
    if ~ status then return

    ;; extract components...
    meta = extract.meta
    if ext eq 0 then begin
        diag  = extract.diag_blue
        diag2 = fextract.diag_blue
    endif else begin
        diag  = extract.diag_red
        diag2 = fextract.diag_red
    endelse 
    
    ;; give feedback...
    print, '  file      = ', meta.filename
    print, '  fiducial  = ', fextract.meta.filename
    print, '  target    = ', meta.std_name
    print, '  grating   = ', meta.grating
    print, '  cenlam    = ', meta.central_wave
    print, '  dwfilnam  = ', meta.blocking
    print, '  extension = ', ext
    print, '  comment   = ', oldcomment

    ;; extract components from diag...
    trace = diag.trace
    profile = diag.profile
    spectrum = diag.spectrum
    extract = diag.extract
    chip = diag.chip
    
    ;; extract components from diag...
    trace2 = diag2.trace
    profile2 = diag2.profile
    spectrum2 = diag2.spectrum
    extract2 = diag2.extract
    chip2 = diag2.chip

    ;; display the trace...
    display_trace, TRACE=trace, FIDUCIAL=TRACE2, PROFILE=profile, $
                   XPOS=xpos, WINDOW=trace_window

    ;; display the extracted profile...
    display_extract, EXTRACT=extract, FIDUCIAL=extract2, XPOS=xpos, $
                     window=extract_window

    ;; display the spectrum...
    display_spectrum, SPECTRUM=spectrum, FIDUCIAL=spectrum2, $
                      WINDOW=spec_window

    ;; display the initial profile...
    display_profile, PROFILE=profile, FIDUCIAL=profile2, XPOS=xpos, $
                     WINDOW=profile_window

    ;; display image and trace...
    view_trace, infile, chip, trace, FIDUCIAL=trace2

    ;; determine next action...
    while 1 do begin
        prompt = '(A)ccept, (R)eject, (C)enter, (F)ollowup, (V)iew in ATV, or (Q)uit ['+default+']: '
        answer = null
        read, answer, prompt=prompt
        if answer eq null then answer = default
        answer = strupcase( answer)
        
        if answer eq 'Q' then begin
            print, 'Quitting...'
            exit_status=1
            return

        endif else if answer eq 'A' then begin
            print, 'Results accepted.'

        endif else if answer eq 'R' then begin

            prompt = 'Please enter comment (or none for no change): '
            reason = null
            read, reason, prompt=prompt
            if answer eq null then answer = oldcomment
            print, 'Adding ', infile, ' to the list of rejects.'
            update_review_database, review_database, review, infile, $
              status='bad', comment=reason
            return

        endif else if answer eq 'F' then begin

            prompt = 'Please enter comment: '
            reason = null
            read, reason, prompt=prompt
            if answer eq null then answer = oldcomment
            print, 'Adding ', infile, ' to the list of followups.'
            update_review_database, review_database, review, infile, $
              status='review', comment=reason
            return

        endif else if answer eq 'C' then begin

            ;; select and redisplay the profile to ensure that we get
            ;; the mapping correct...
            wset, profile_window
            display_profile, PROFILE=profile, XPOS=xpos
            print, 'Please place cursor at desired center and click the mouse...'
            cursor, xc, yc, /up, /data
            print, 'Got new center at X=', xc
            update_center_list, center_list, infile, xc

            ;; trigger re-analysis of spectrum...
            print, 'Re-extracting spectrum...'
            do_deimos_throughput, input=[infile], /clobber

            ;; force re-do of this extension...
            ext--

        endif else if answer eq 'V' then begin
            atv_activate
            continue

        endif else begin
            print, 'Sorry, "', answer, '" is not a valid option.'
            continue

        endelse
        
        ;; if we got here, then we got a valid response...
        break
    endwhile

    ;; if we got here then the target is approved...
    if ext eq 1 then begin
        message, 'Adding target to approved list.', /info
        update_review_database, review_database, review, infile, $
          status='good'
    endif

endfor

end
