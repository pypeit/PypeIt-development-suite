;-----------------------------------------------------------------------
function fiducial_file, META=meta, $
  GRATING=grating, $
  WAVELENGTH=wavelength, $
  FILTER=filter, $
  TARGET=target, $
  EXTRACT=extract, $
  ERRMSG=errmsg
;-----------------------------------------------------------------------
;+
; NAME:
;	FIDUCIAL_FILE
;
; PURPOSE:
;	Given the grating, wavelength, filter, and target name of an
;	observation, this function will parse the fiducial file list
;	and will return the name of the corresponding fiducial file,
;	or an empty string if no fiducial exists.  By default, the
;	returned value is the name of the RAW file for the fiducial.
;
; CATEGORY:
;	IPM
;
; CALLING SEQUENCE:
;	result = fiducial_file()
;
; INPUTS:
;       meta: structure with metadata
;
;	    grating: name of the grating (e.g., '1200G')
;
;	    wavelength: integer central wavelength (e.g., 7000)
;
;       filter: name of the blocking filter (e.g., 'OG550')
;
;       target: name of the standard star (e.g., 'HZ44')
;
;       EXTRACT: if set, return the filename of the extracted SAV
;       file instead of the raw fits file.
;
; OUTPUTS:
;	Function returns the full disk name of the fiducial file.
;
; PROCEDURE:
;       - read fiducial file table to retrieve directory and root name
;       - assemble file name
;       - verify existence of file
;
; EXAMPLE:
;       ffile = fiducial_file( grating='1200G', wavelength=7000, \
;               filter='GG455', target='HZ44', dir='2013apr01')
;
; AUTHOR:
;	Lysha Matsunobu, WMKO
;
; MODIFICATION HISTORY:
; 	2012-Jan-31	LM	Original version
;       2014-Feb-18     GDW     Add META keyword
;       2017-Dec-15     CA      Use variable $DEIMOS_THRU_DIR
;                               instead of harcoded data directory
;                               name.
;-      2018-Nov-15     TL      Add dir for fiducial duplicity check
;-----------------------------------------------------------------------

null = ''
errmsg = null

thru_dir = getenv('DEIMOS_THRU_DIR')

if keyword_set(META) then begin

    ;; get info from metadata structure...
    target     = meta.std_name
    grating    = meta.grating
    wavelength = nint(meta.central_wave)
    filter     = meta.blocking
    dir        = strmid(meta.date, 5,4) + strmid(strlowcase(meta.date), 2,3) $
                 + strmid(meta.date, 0,2)

endif else begin

    ;; verify keywords are set...
    if ~ keyword_set( GRATING)    then message, 'ERROR: GRATING keyword undefined'
    if ~ keyword_set( WAVELENGTH) then message, 'ERROR: WAVELENGTH keyword undefined'
    if ~ keyword_set( FILTER) then message, 'ERROR: FILTER keyword undefined'
    if ~ keyword_set( TARGET) then message, 'ERROR: TARGET keyword undefined'
    if ~ keyword_set( DATE) then message, 'ERROR: DATE keyword undefined'

endelse 

help, target, grating, wavelength, filter

;; verify existence of file...
infile = getenv('DEIMOS_THRU_DIR') + '/dat/fiducials.csv'
if ~ file_test( infile) then begin
    errmsg = 'ERROR: required fiducial file "' + infile + '" not found -- abort!'
    message, errmsg
endif 

;; parse input file...
result = read_csv( infile)

;; decompose struct...
gratings    = result.field1
wavelengths = result.field3
filters     = result.field2
targets     = result.field4
dirs        = result.field5
filename    = result.field6

;; find lines which match input...
match = where( gratings eq grating and $
               wavelengths eq wavelength and $
               filters eq filter and $
               targets eq target and $
               dirs eq dir, count)

;; verify that we got the correct number of returns...
pathname = null                 ; default is null string
if count eq 1 then begin

    ffile = thru_dir + '/raw/' $
            + dirs[match] + '/' + $
            filename[match] + '.fits'
    gzfile = ffile+'.gz'

    ;; check existence of file...
    if file_test( ffile) then begin
        pathname = ffile
    endif else if file_test( gzfile) then begin
        pathname = gzfile
    endif else begin
        errmsg = 'WARNING: expected fiducial file "'+ffile+'" not found!'
        message, errmsg, /info
        pathname = null
    endelse 

endif else if count eq 0 then begin
    sep = '/'
    errmsg = 'WARNING: no fiducial file found for ' $
             +target+sep+grating+sep+stringify(wavelength)+sep+filter
    message, errmsg, /info
endif else begin
    sep = '/'
    errmsg = 'WARNING: multiple fiducial files for ' $
             +target+sep+grating+sep+stringify(wavelength)+sep+filter
    message, errmsg, /info
endelse 

;; if the extracted file was requested, then return it...
if keyword_set(EXTRACT) and pathname ne null then begin
    ffile = raw2extract(pathname) + '.sav'
    gzfile = ffile+'.gz'

    ;; check existence of file...
    if file_test( ffile) then begin
        pathname = ffile
    endif else if file_test( gzfile) then begin
        pathname = gzfile
    endif else begin
        errmsg = 'WARNING: expected fiducial file "'+ffile+'" not found!'
        message, errmsg, /info
        pathname = null
    endelse 

endif 

return, pathname

end
