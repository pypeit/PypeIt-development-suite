;-----------------------------------------------------------------------
function stringify, value, format
;-----------------------------------------------------------------------
;+
; NAME:
;	STRINGIFY
;
; PURPOSE:
;	This function converts a real value into a string based on the 
;	specified format.
;
; CATEGORY:
;	Formatting.
;
; CALLING SEQUENCE:
;	Result = STRINGIFY( Value, Format)
;
; INPUTS:
;	Value	Value which is to be converted.
;
; OPTIONAL INPUTS:
;	Format	Optional format statement
;	
; OUTPUTS:
;	This function returns a string version of the input array with 
;	leading and trailing spaces removed.
;
; RESTRICTIONS:
;	None
;
; EXAMPLE:
;	1) format a number into a string:
;		label = "X=" + stringify(x)
;
;	1) format a number into a string, with formatting:
;		label = "X=" + stringify(x,'(f15.3)')
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	1999-Dec-14	GDW	Original version
;-
;-----------------------------------------------------------------------

; verify input...
np = n_params()
if np lt 1 or np gt 2 then message, 'wrong number of parameters'

; format the value into an appropriate string...
if np gt 1 then begin
    text = string( format=format, value)
endif else begin
    text = string( value)
endelse

; remove leading and trailing spaces...
return, strtrim( text, 2)
end
