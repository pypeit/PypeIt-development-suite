;-----------------------------------------------------------------------
pro run_review_deimos_throughput, TARGNAME=targname, GRATENAM=gratenam, $
  WAVELEN=wavelen, DWFILNAM=dwfilnam, TEST=test, QUERY=query, $ 
  REVIEW=review
;-----------------------------------------------------------------------
;+
; NAME:
;	RUN_REVIEW_DEIMOS_THROUGHPUT
;
; PURPOSE:
;       This procedure will select DEIMOS throughput data meeting
;       various user-specified criteria and will run the
;       review_deimos_throughput program to analyze them.
;
; CATEGORY:
;       IPM
;
; CALLING SEQUENCE:
;       run_review_deimos_throughput
;
; KEYWORD PARAMETERS:
;	TARGNAME: Set this keyword to select spectra of the specified
;       target name (header keyword TARGNAME)
;
;       GRATENAM: Set this keyword to select spectra with the
;       specified grating name (header keyword GRATENAM)
;
;       WAVELEN: Set this keyword to select spectra with the specified
;       wavelength (header keyword WAVELEN)
;
;       DWFILNAM: Set this keyword to select spectra with the specified
;       filter (header keyword DWFILNAM)
;
;       QUERY: set this keyword to query the user before reviewing
;       each file
;
;       TEST: set this keyword to print the files to be reviewed
;       without running the review program
;
; OUTPUTS:
;       See review_deimos_throughput
;
; RESTRICTIONS:
;	The input file list 'input.lst' should be updated before
;	running this program by executing the program
;	build_throughput_list.
;
; EXAMPLE:
;       1) Review all available spectra:
;              > run_review_deimos_throughput
;
;       2) Review all available spectra of target 'Feige 34':
;              > run_review_deimos_throughput, targname='Feige 34'
;
;       3) Review all available spectra with grating '600ZD':
;              > run_review_deimos_throughput, gratenam='600ZD'
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2011-Oct-31	GDW	Original version
;-
;-----------------------------------------------------------------------

debug = 0B

;; define files...
thru_dir = getenv('DEIMOS_THRU_DIR')
if thru_dir eq '' then $
  message, 'DEIMOS_THRU_DIR is not defined -- abort!'
data_dir = thru_dir + '/dat/'
input_list   = data_dir + 'input.lst'
review_list  = data_dir + 'review.lst'

;; get input lists...
result = read_csv( input_list)
    
;; extract fields...
files       = result.field1
targnames   = result.field2
gratenames  = result.field3
wavelengths = result.field4
dwfilnames  = result.field5
dateobss    = result.field6

n = n_elements(files)

;; read review database...
review_database = data_dir + 'review.fits'
db = mrdfits( review_database, 1)

;; construct status array...
status = strarr(n)
for i=0,n-1 do begin
    
    ;; retrieve record for this file...
    record = parse_review_database( db, files[i])

    ;; insert into status array...
    if size(record, /type) ne 8 then begin
        status[i] = ''
    endif else begin
        status[i] = record.status
    endelse 
    
endfor 

;; create array of flags and set all to "on"
select = bytarr(n) + 1B

;; allocate comment array...
comment = strarr(n)

;; cut on target name...
if keyword_set( TARGNAME) then begin
    if targname ne '' then begin
        if debug then print, 'select targname=', targname
        bad = where( targnames ne targname, count)
        if count gt 0 then select[bad] = 0B
    endif else begin
        message, 'ERROR: null TARGNAME'
    endelse 
endif 

;; cut on gratenam...
if keyword_set( GRATENAM) then begin
    if gratenam ne '' then begin
        if debug then print, 'select gratenam=', gratenam
        bad = where( gratenames ne gratenam, count)
        if count gt 0 then select[bad] = 0B
    endif else begin
        message, 'ERROR: null GRATENAM'
    endelse 
endif 

;; cut on dwfilnam...
if keyword_set( DWFILNAM) then begin
    if dwfilnam ne '' then begin
        if debug then print, 'select dwfilnam=', dwfilnam
        bad = where( dwfilnames ne dwfilnam, count)
        if count gt 0 then select[bad] = 0B
    endif else begin
        message, 'ERROR: null DWFILNAM'
    endelse 
endif 

;; cut on wavelen...
if keyword_set( WAVELEN) then begin
    if wavelen gt 0 then begin
        if debug then print, 'select wavelen=', wavelen
        bad = where( wavelengths ne wavelen, count)
        if count gt 0 then select[bad] = 0B
    endif else begin
        message, 'ERROR: illegal WAVELEN'
    endelse 
endif 

;; cut on review...
if keyword_set( REVIEW) then begin

    ;; read list...
    if debug then print, 'select status=review'

    bad = where( status ne 'review', count)
    if count gt 0 then select[bad] = 0B

endif

;; execute reviewing program...
good = where( select eq 1B, count)
print, 'Found ', count, ' files matching criteria'
space = ' '
for j=0,count-1 do begin
    i = good[j]

    BREAKNAME, files[i], Dirname, Rootname, Extn
    if j eq 0 then begin
        print, 'i', 'rootname', $
               'targname', $
               'gratenam', $
               'wavelen', $
               'dwfilnam', $
               'dateobs', $
               'status', $
               format='(a5,a12,a15,a10,a8,a9,a12,a7)'
        print, '-', '--------', $
               '--------', $
               '--------', $
               '-------', $
               '--------', $
               '-------', $
               '------', $
               format='(a5,a12,a15,a10,a8,a9,a12,a7)'
    endif 
    
    print, i, rootname, $
           targnames[i], $
           gratenames[i], $
           wavelengths[i], $
           dwfilnames[i], $
           dateobss[i], $
           status[i], $
           format='(i5,a12,a15,a10,i8,a9,a12,a7)'
    
    if keyword_set(QUERY) then begin
        answer = ''
        read, 'Review file '+files[i]+'? (y/n): ', answer
        if answer eq '' or answer eq 'y' or answer eq 'Y' then begin
            print, 'OK'
        endif else begin
            print, 'skipping!'
            continue
        endelse 
        
    endif
    
    ;; review...
    if ~ keyword_set(TEST) then begin
        review_deimos_throughput, files[i], EXIT=exit
        if exit ne 0 then return
    endif 
endfor 

print, ''
print, 'n    =', n
print, 'good =', count

end
