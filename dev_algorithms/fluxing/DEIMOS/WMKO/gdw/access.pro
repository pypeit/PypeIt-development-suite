;-----------------------------------------------------------------------
function access, filename, mode
;-----------------------------------------------------------------------
;+
; NAME:
;	ACCESS
;
; PURPOSE:
;	This function determines whether the specified file exists and
;	can be opened in the requested mode.  The mode is presumed to
;	be "readonly" by default. 
;
; CATEGORY:
;	File I/O
;
; CALLING SEQUENCE:
;	Result = ACCESS( Filename, Mode)
;
; INPUTS:
;	Filename [string]: Name of the file to check accessibility of
;
; OPTIONAL INPUTS:
;	Mode [char]:       Access mode to be verified.  Legal values:
;                            "R" to test for "read" access [default]
;                            "W" to test for "write" access
;	
; OUTPUTS:
;	This function returns the "true" value (1L) is the file exists 
;	and can be opened in the specified mode.  Otherise, it returns 
;	"false" (OL).
;
; SIDE EFFECTS:
;	The specified file is opened for read or write access, then
;	closed again after verification.
;
;
; EXAMPLES:
;	1) Test the existence of the file "foo.dat":
;		status = ACCESS( 'foo.dat')
;	2) Test whether we have write access to the file "foo.dat":
;		status = ACCESS( 'foo.dat','w')
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
;	1999-Jul-27  [GDW] Original version
;	1999-Jul-28  [GDW] Changed "close" statements to "free_lun" in
;	order to prevent exhaustion of logical unit numbers.
;-
;-----------------------------------------------------------------------

; default value is false...
value = 0L

; check for missing filename...
if( n_params() lt 1) then begin
    message, 'filename argument not specified'
    return, value
endif

; check for default mode...
if( n_params() lt 2) then mode = 'r'
lmode = strlowcase( mode)

if lmode eq 'r' or lmode eq 'w' then begin

    ; open read-only mode...
    openr, unit, filename, /get_lun, error=err
    if (err eq 0) then begin
        free_lun, unit
        value=1L
    endif

    ; if read succeeds, then we've established existence.  Now check
    ; write access...
    if value and lmode eq 'w' then begin
        value = 0L
        openw, unit, filename, /get_lun, error=err
        if (err eq 0) then begin
            free_lun, unit
            value=1L
        endif
    endif
endif else begin
    message, 'specified mode not recognized.  Use "R" or "W".', /inf
endelse

return, value

end
