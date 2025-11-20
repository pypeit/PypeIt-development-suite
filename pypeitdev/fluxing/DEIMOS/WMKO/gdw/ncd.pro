pro ncd, search_dir
;-----------------------------------------------------------------------
; ncd.cl / GDW / 28 Jul 93
;
; Procedure to nicely change directories.  Takes a string as input,
; and changes directories to the first subdirectory whose name
; matches the string.
;
; Parameters:
;
; dir [string]
;	search string indicating directory to seek
;
; ncddat [string]
;	file containing directory listings.  See
; 	~gwirth/bin/mkncd to create this file.
;
;-----------------------------------------------------------------------

val = strtrim( search_dir,2)
ncddat = '~/.ncddat'            ; file listing NCD info

; initialize variables...
null = ''                       ; null value
s1 = null                       ; matching directory
s2 = null                       ; current line of ncddat file

np = n_params()

; if a search string was supplied, then look for it...
if np gt 0 then begin
    if val ne null then begin

        openr, unit, ncddat, /get_lun
        
        ; scan through the input file searching for substring...
        while not eof(unit) and s1 eq null do begin
            readf, unit, s2
            if strpos( s2, val) gt -1 then $
              s1 = s2
        endwhile
        
        ; print error message if none found...
        if s1 eq null then begin
            message, 'no matching subdirectory found', /informational
            close, unit
            return
        endif
        
        close, unit
    endif
endif

; effect the change requested...
cd, s1
print, s1

end
