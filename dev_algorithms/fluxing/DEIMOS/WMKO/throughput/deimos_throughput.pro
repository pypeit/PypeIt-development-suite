;-----------------------------------------------------------------------
pro deimos_throughput, infile, wavelength, counts, $
  OUTFILE=outfile,     $
  SAVFILE=savfile,     $
  META=meta,           $
  TRACEFILE=tracefile, $
  TRACE_CENTER=trace_center, $
  STATUS=status,       $
  ERRMSG=errmsg,       $
  VERBOSE=verbose,      $
  LOGUNIT=logunit
;-----------------------------------------------------------------------
;+
; NAME:
;	DEIMOS_THROUGHPUT
;
; PURPOSE:
;       Given a standard star spectrum acquired with DEIMOS, this
;       procedure will generate a throughput curve.
;
; CATEGORY:
;	Data analysis
;
; CALLING SEQUENCE:
;       DEIMOS_THROUGHPUT, infile, wavelengths, counts
;
; INPUTS:
;       infile = name of the FITS image containing the spectrum
;
; KEYWORD PARAMETERS:
;       META    = name of output keyword to contain metadata
;       OUTFILE = name of optional output file containing data.
;       SAVFILE = if set, save variables to named IDL SAVE file 
;       TRACEFILE = name of file to write trace data into
;       TRACE_CENTER = if set, use specified column as trace center
;       VERBOSE = if set, provide wordy output on terminal
;
; OUTPUT KEYWORDS:
;       ERRMSG  = output error message (blank if no error)
;       STATUS  = output status (1=good, 0=bad)
;
; OUTPUTS:
;       wavelengths = set of wavelengths at which throughput is
;               computed [Angstroms]
;       counts = array of detected electrons received at 
;               corresponding wavelengths [e-]
;
; SAVEFILE CONTENTS:
;       extract:
;               meta: metadata
;               spec { wavelength, counts }
;               diag_blue { 
;                       infile: name of input raw FITS file
;                       chip:   CCD number
;                       color: R or B
;                       trace {
;                               good: { count, u, v}
;                               bad:  { count, u, v}
;                               fit:  { status, u, v}
;                               }
;                       profile {
;                               u = column number
;                               flux = flux at column
;                               v1 = starting row number for sum
;                               v2 = ending row number for sum
;                               sky = sky value
;                               s2n = signal to noise at peak
;                               u_center = measured center
;                               }
;                       extract { # parameters from extract_spectrum
;                               row:i, $     ; row number
;                               center:v[i], $ ; trace center at this row
;                               radius:radius, $ ; extraction radius
;                               background_width:background_width, $ ; background width
;                               buffer:1, $   ; buffer
;                               start:start, $ ; start of plot window
;                               stop:stop, $ ; end of plot window
;                               peak:peak, $ ; peak value within window
;                               sky:sky, $   ; sky value
;                               x1:x1, $     ; start of extraction window
;                               x2:x2, $     ; end of extraction window
;                               xleft1:xleft1, $ ; start of left sky window
;                               xleft2:xleft2, $ ; end of left sky window
;                               xright1:xright1, $ ; start of right sky window
;                               xright2:xright2, $ ; end of right sky window
;                               x:x[start:stop], $ ; pixel numbers
;                               flux:line[start:stop] } ; data values
;                               }
;                       spectrum {
;                               pixel
;                               flux
;                               }
;
;               diag_red: same as diag_blue
;
; INVOKED SUBROUTINES:
;       deimos_read_chip
;       trace_spectrum
;       extract_spectrum
;       deimos_get_wavelengths
;
; PROCEDURE:
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2005-Aug-30	GDW	[v0.0] Original version
;       2008-May-08     GDW     [v0.1] Extensive revision
;-
;-----------------------------------------------------------------------

;; define constants...
min_s2n = 1.
u_window = 90
v_window = 51
max_shift = 1.0
dispaxis = 2
n_truncate = 4
version = '0.1'
status = 0B                     ; bad by default
errmsg = ''                     ; error message is empty by default

;; get rootname...
breakname, infile, dirname, rootname, extn

;; get image header parameters...
header   = headfits( infile)
slmsknam = strtrim(sxpar( header, 'slmsknam'),2)
mosmode  = strtrim(sxpar( header, 'mosmode'),2)
targname = strtrim(sxpar( header, 'targname'),2)
airmass  = sxpar( header, 'airmass')
gratepos = sxpar( header, 'gratepos')
keyword  = string( format='(a,i1,a)', 'g', gratepos, 'tltwav')
gXtltwav = sxpar( header, keyword)
cenlam   = nint(  gXtltwav)
exptime  = sxpar( header, 'ttime')
gratenam = strtrim( sxpar( header, 'gratenam'),2)
dwfilnam = strtrim( sxpar( header, 'dwfilnam'),2)

;; skip pathological cases...
if mosmode ne 'Spectral' then begin
    errmsg = 'not in Spectral readout mode -- skipping'
    plog, logunit, 'WARNING: image '+rootname + ' ' + errmsg
    return
endif
if slmsknam ne 'None' then begin
    errmsg = 'not in slitless mode -- skipping'
    plog, logunit, 'WARNING: image '+rootname + ' ' + errmsg
    return
endif

if keyword_set(VERBOSE) then begin
    print, 'file     = ', infile
    print, 'target   = ', targname
    print, 'grating  = ', gratenam
    print, 'cenlam   = ', cenlam
    print, 'dwfilnam = ', dwfilnam
endif 

;; set up metadata struct...
meta = { $
       observatory: ' ', $
       instrument: ' ',        $ ; ID number
       date: ' ',   $           ; DDMMYYYY (e.g. 20JAN2004)
       grating: ' ',   $        ; grating name
       gratepos: 0, $           ; grating position
       central_wave: 0., $      ; Central wavelength of grating
       blocking: '', $          ; Blocking filter
       airmass: 0., $
       spec_bin: 0, $
       slit_width: 0., $
       std_name: ' ', $
       ra: ' ',    $            ; RA (J2000)
       dec: ' ',    $           ; DEC (J2000)
       conditions: ' ',    $    ; sky conditions
       filename: ' ', $
       frame: 0L,   $
       program: '', $
       version: '' }

;; define metadata from keywords...
meta.program = 'deimos_throughput'
meta.version = version
meta.observatory = 'Keck'
meta.instrument = 'DEIMOS'
meta.frame = sxpar( header, 'frameno')
meta.filename = infile
meta.grating = strtrim(gratenam,2)
meta.central_wave = gXtltwav
meta.blocking = dwfilnam
meta.airmass = sxpar( header, 'airmass')
binning = sxpar( header, 'binning')
meta.spec_bin = long( (strsplit( binning, ',', /extract))[1])
meta.slit_width = 99.
meta.ra = sxpar( header, 'targra')
meta.dec = sxpar( header, 'targdec')
meta.gratepos = sxpar( header, 'gratepos')
dateobs = sxpar( header, 'DATE-OBS')
date = strmid(dateobs,8,2)+$
       x_getmonth(long(strmid(dateobs,5,2)))+$
       strmid(dateobs,0,4)
meta.date = date
meta.std_name = targname
title = targname + ' ' + gratenam + '@' + stringify(cenlam) $
        + ' ' + dwfilnam + ' ' + infile

;; define wavelength limits...
filters = ['BAL12','GG400','GG455','GG495','OG550']
cut_on  = [3500., 4000., 4400., 4950., 5500.]
cut_off = [8000., 8000., 8800., 10000., 11000.]
good = where( filters eq meta.blocking, count)
if count ne 1 then begin
    errmsg = 'incorrect number of filters matching DWFILNAM ('+meta.blocking+')'
    plog, logunit, 'WARNING: image '+rootname + ' ' + errmsg
    return
endif 
index = good[0]
good_wavelength_min = cut_on[index]
good_wavelength_max = cut_off[index]

if keyword_set(VERBOSE) then $
  plog, logunit, 'Processing file '+rootname

;; define starting point for trace.  Note that we do not use the
;; center row because the spectrum may not reach that far.  On the
;; blue side we start at 10% below the cutoff and in red 10% above.
nrows = 4096
v0 = nint([0.9*nrows,0.1*nrows])

;; loop over chips...
chips = [3,7]
colors = [ 'B', 'R']
n_chips = 2
n_iter = 0

for i=0,n_chips-1 do begin

    n_iter += 1
    chip = chips[i]
    color = colors[i]

    ;; read image...
    plog, logunit, '  Reading '+color
    image = deimos_read_chip( infile, chip, header=header)

    ;; get size...
    buf = size(image)
    nx = buf[1]
    ny = buf[2]

    ;; if trace center not set then set to Nan...
    if ~ keyword_set(trace_center) then trace_center=!values.f_nan

    ;; trace spectrum...
    plog, logunit, '  Tracing '+color
    trace_status = trace_spectrum( image, u, v, $
                                   U0=trace_center, $
                                   V0=v0[i], $
                                   U_WINDOW=u_window, $
                                   V_WINDOW=v_window, $
                                   U_GOOD=u_good, V_GOOD=v_good, $
                                   DISPAXIS=dispaxis, $
                                   MAX_SHIFT=max_shift, $
                                   MIN_S2N=min_s2n, PROFILE=profile, BPM=bpm)

    ;; fit the trace...
    if ( trace_status ) then begin
        plog, logunit, '  Good trace in '+color
        fit_struct = x_setfitstrct( FLGREJ=1)
        buf = x_fitrej( v_good, u_good, FITSTR=fit_struct)
        v_trace = findgen( ny)
        u_trace = x_calcfit( v_trace, FITSTR=fit_struct)
    endif else begin
        plog, logunit, '  WARNING: bad trace in '+color
        v_trace = findgen(ny)
        u_trace = replicate( !values.f_nan, ny)
    endelse 

    ;; compile diagnostics...
    select = where( bpm ne 0B, ngood)
    if ngood gt 0 then begin
        good = {count:ngood, u:u[select], v:v[select]}
    endif else begin
        good = {count:ngood, u:[0], v:[0]}
    endelse 

    select = where( bpm eq 0B, nbad)
    if nbad gt 0 then begin
        bad = {count:nbad, u:u[select], v:v[select]}
    endif else begin
        bad = {count:nbad, u:[0], v:[0]}
    endelse 
    plog, logunit, '  Good points='+stringify(ngood)
    plog, logunit, '  Bad points ='+stringify(nbad)

    fit = {status:trace_status, u:u_trace, v:v_trace}
    trace = {good:good, bad:bad, fit:fit}

;    ;; display trace...
;    if keyword_set(PSFILE) then begin
;        if keyword_set(VERBOSE) then $
;          message, '  Plotting trace in plot window', /info

;        ;; generate plot area...
;        plot, v, u, /nodata, xtitle='Y [px]', ytitle='X [px]', /yno, $
;              xstyle=1, charsize = 1.25

;        ;; plot good data...
;        if good.count gt 0 then oplot, good.v, good.u, psym=6, symsize=0.25

;        ;; plot bad data...
;        if bad.count gt 0 then oplot, bad.v, bad.u, psym=7, symsize=0.5
;        nxyouts, xpos, 0.55, stringify(bad.count)+' bad points'
;        nxyouts, xpos, 0.5, 'Sky='+stringify(nint(profile.sky),'(i10)')

;        ;; add trace...
;        if ( fit.status ) then begin
;            oplot, fit.v, fit.u
;        endif else begin
;            nxyouts, 0.5, 0.5, 'WARNING: BAD TRACE', align=0.5
;        endelse 

;        ;; annotate...
;        nxyouts, xpos, ypos, 'Trace'

;    endif 

    if fit.status then begin

        ;; generate optional trace file...
        if keyword_set(TRACEFILE) then begin
            tfile = tracefile + '_' + stringify(chip) + '_trace.region'
            openw, tunit, tfile, /get_lun
            printf, tunit, format='("# image = ", a)', infile
            printf, tunit, format='("# chip = ", i)', chip
            npts = n_elements( u_trace)
            for j=0,npts-1 do begin
                printf, tunit, format='(f,f)', $
                        u_trace[j], v_trace[j]
            endfor 
            free_lun, tunit            
        endif

    endif else begin
        errmsg = 'WARNING: bad trace for '+rootname
        plog, logunit, errmsg

;        if keyword_Set(PSFILE)  then begin

;            ;; plot initial profile...
;            plot, profile.u, profile.flux, $
;                  xtitle='Column [px]', ytitle='Sky-subtracted Flux', $
;                  xstyle=1, charsize = 1.25
;            oplot, [profile.u_center,profile.u_center], !y.crange, linestyle=1
;            nxyouts, xpos, ypos, 'Profile'
;            nxyouts, xpos, 0.5, $
;                'Profile at Rows '+stringify(profile.v1)+'-'+stringify(profile.v2)
;            xyouts, profile.u_center, !y.crange[0], $
;                    ' Center='+stringify(profile.u_center,'(f15.1)'), $
;                    orient=90.
;            nxyouts, xpos, 0.10, 'Sky='+stringify(nint(profile.sky),'(i10)')
;            nxyouts, xpos, 0.05, 'S/N='+stringify(profile.s2n,'(f15.1)')
;            if keyword_set(PSFILE) then psclose
;        endif 
        ;; restore original graphics state...
        return
    endelse

;    ;; optional feedback...
;    if keyword_set(INTERACTIVE) then begin
;        if keyword_set(VERBOSE) then $
;          message, '  Displaying '+rootname+' in ATV...', /info
;        atv, image
;        if keyword_set(VERBOSE) then $
;          message, '  Overplotting trace points in green', /info
;        atvplot, u, v, color='green', psym=2
;        if ( trace_status ) then begin
;            if keyword_set(VERBOSE) then $
;              message, '  Overplotting trace fit in red', /info
;            atvplot, u_trace, v_trace, color='red'
;        endif else begin
;            message, '  Warning: bad trace fit!', /info
;        endelse 
;        message, '  Please press Q in the ATV window to continue...', /info
;        atv_activate
;    endif 

    ;; get spectrum...
    if keyword_set(VERBOSE) then $
      plog, logunit, '  Extracting spectrum in ' + color
    extract_spectrum, image, v_trace, u_trace, spec, extract_params, $
                      DISPAXIS=dispaxis, $
                      DISPLAY=0B, VERBOSE=0, PARAM_ROW=v0[i]
    n_pix = n_elements(spec)
    spectrum = { pixel:v_trace, flux:spec }

;    ;; show profile...    
;    if keyword_set(PSFILE) then begin

;        ;; scale Y axis from sky to peak with 10% buffer
;        dy = extract_params.peak - extract_params.sky
;        y_min = extract_params.sky  - 0.1 * dy
;        y_max = extract_params.peak + 0.1 * dy

;        ;; generate axes...
;        plot, extract_params.x, extract_params.flux, $
;              xtitle='Column [px]', ytitle='Sky-subtracted Flux', $
;              charsize=1.25, yrange=[y_min,y_max], ystyle=1

;        ;; overplot extraction region...
;        oplot, [extract_params.x1,extract_params.x2], $
;               [extract_params.peak,extract_params.peak], psym=-1

;        ;; mark left sky region...
;        oplot, [extract_params.xleft1, extract_params.xleft2], $
;               [extract_params.sky,extract_params.sky], psym=-1

;        ;; mark right sky region...
;        oplot, [extract_params.xright1, extract_params.xright2], $
;               [extract_params.sky,extract_params.sky], psym=-1

;        ;; indicate sky level as dotted horizontal line...
;        oplot, !x.crange, [extract_params.sky,extract_params.sky], linestyle=1
            
;        ;; plot dotted vertical line at the measured center...
;        oplot, [extract_params.center,extract_params.center], !y.crange, $
;               linestyle=1
;        xyouts, extract_params.center, extract_params.sky, $
;                ' Center='+stringify(extract_params.center,'(f15.1)'), $
;                orient=90.

;        ;; annotate...
;        nxyouts, xpos, ypos, 'Profile'
;        nxyouts, xpos, 0.55, $
;                 'Row='+stringify(extract_params.row)
;        nxyouts, xpos, 0.50, $
;                 'Sky='+stringify(nint(extract_params.sky),'(i10)')
;    endif

;    ;; more optional feedback...
;    if keyword_set(PSFILE) then begin
;        plot, v_trace, spec, xtitle='row', ytitle='flux [e-/px]', $
;              xstyle=1, charsize = 1.25
;        nxyouts, xpos, ypos, 'Raw spectrum'
;    endif

    ;; generate wavelengths...
    if i eq 0 then begin
        wave = deimos_get_wavelengths( gratenam, cenlam, n_pix, /BLUE)
    endif else begin
        wave = deimos_get_wavelengths( gratenam, cenlam, n_pix, /RED)
    endelse

    ;; concatenate spectra...
    if i eq 0 then begin
        wavelength = wave[0:n_pix-n_truncate-1]
        counts     = spec[0:n_pix-n_truncate-1]
    endif else begin
        wavelength = [wavelength,wave[0:n_pix-1]]
        counts     = [counts,spec[0:n_pix-1]]
    endelse

    ;; compile diagnostics...
    diag = {infile:infile, $
            chip:chip, $
            color:color, $
            trace:trace, $
            profile:profile, $
            extract:extract_params, $
            spectrum:spectrum}
    if i eq 0 then begin
        diag_blue = diag
    endif else begin
        diag_red = diag
    endelse 

endfor

;if keyword_set(PSFILE) then begin
;    if keyword_set(VERBOSE) then $
;      message, '  Plotting spectrum', /info
;    if keyword_set(PSFILE) then $
;      !p.multi = [1,1,3]        ; plot on bottom 1/3 of chart
;    plot, wavelength, counts, linestyle=1, $
;          xtitle='Wavelength ['+Angstrom+']', ytitle='Flux [e-/px]', $
;          xstyle=1, charsize = 1.25
;    nxyouts, xpos, ypos, 'Spectrum'
;endif

;; remove bad wavelengths...
good = where( wavelength ge good_wavelength_min and $
              wavelength le good_wavelength_max, count)
if count lt 1 then begin
    errmsg = 'no good wavelengths in spectrum'
    plog, logunit, errmsg
    return
endif 
wavelength = wavelength[good]
counts = counts[good]

if keyword_set(PSFILE) then begin
    oplot, wavelength, counts
endif

;; divide the counts by the exposure time...
counts = counts/exptime

;; generate a structure with the spectrum...
n = n_elements(wavelength)
line = {wavelength:0., counts:0.}
spec = REPLICATE(line,n)
spec.wavelength = wavelength
spec.counts     = counts

;; optional output file...
if keyword_set( OUTFILE) then begin

    ;; dump meta data to file...
    sxaddpar, outheader, 'EXTNAME', 'METADATA', 'Name of this extension'
    mwrfits, meta, outfile, outheader, /create

    sxaddpar, outheader, 'EXTNAME', 'SPECTRUM', 'Name of this extension'
    mwrfits, spec, outfile, outheader
    if keyword_set(VERBOSE) then plog, logunit, 'Writing output file '+outfile

endif 

if keyword_set( SAVFILE) then begin
    plog, logunit, 'Writing SAV file '+savfile
    extract = {meta:meta, $
               spec:spec, $
               diag_blue:diag_blue, $
               diag_red:diag_red }
    save, extract, filename=savfile, /verbose
endif 

;; mark as "good"...
status = 1B

end
