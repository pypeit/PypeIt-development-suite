;-----------------------------------------------------------------------
function cdf_frac, array, frac, epsilon, MAX_ITER=max_iter, VERBOSE=verbose
;-----------------------------------------------------------------------
;+
; NAME:
;	CDF_FRAC
;
; PURPOSE:
;	This function estimates the value which lies at the frac*100
;	percentile distribution of the data.  For example, if frac=0.5 
;	then this routine returns an estimate of the median.
;
; CATEGORY:
;	Statistics
;
; CALLING SEQUENCE:
;	Result = CDF_FRAC( Array, Frac)
;
; INPUTS:
;	Array:	Data to be analyzed
;	Frac:	Fractional point in the cumulative distribution
;		function to be determined.  Legal values: 0<Frac<1
;
; OPTIONAL INPUTS:
;	Epsilon: maximum error in the returned estimate.
;		Positive values are interpreted as absolute limits, and
;		negative values imply a fractional accuracy (i.e.,
;		Frac=0.5, EPSILON=-0.01 will return a value which is
;		guaranteed to be between the 0.49 and 0.51 percentile).
;		[-0.01]
;
; KEYWORD PARAMETERS:
;	MAX_ITER:	Maximum number of iterations to complete [100]
;	VERBOSE:	Print diagnostic updates? [n]
;
; OUTPUTS:
;	This function returns an estimate of the value which is higher 
;	than the specified fraction of the values in the array.
;	
; EXAMPLE:
;	1. To estimate the median value of an array of data to within 10
;	data units:
;		median = CDF_FRAC( Array, 0.5, 10.)
;
;	2. To estimate the 75th percentile value of an array of data
;	to within 1% of the data range:
;		median = CDF_FRAC( Array, 0.5)
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	1999-Nov-4	GDW	Original version
;-
;-----------------------------------------------------------------------

; verify parameters...
np = n_params()
if( np lt 2 or np gt 3) then begin
    message, 'Required parameters missing.'
endif

; defaults...
if( not keyword_set( MAX_ITER)) then max_iter=100
verbose = keyword_set(VERBOSE)
if( np lt 3) then epsilon=-0.01

; verify values...
if( frac lt 0 or frac gt 1) then begin
    message, 'Frac value is out of range: ' + string(frac)
endif

if( max_iter lt 0) then begin
    message, 'max_iter is out of range: ' + string(max_iter)
endif

; determine the type of data being used...
if( size( array, /type) le 3) then begin
    type_int = 1L
endif else begin
    type_int = 0L
endelse

; find endpoints...
min_val = min( array, max=max_val)

; special cases...
if( frac eq 0) then return, min_val
if( frac eq 1) then return, max_val

; set starting endpoints...
lo = min_val
hi = max_val

; define the number of points below each limit...
n = n_elements( array)
n_lo = 0
n_hi = n-1
n_target = nint(frac*n)
if( epsilon le 0) then n_stop = -epsilon*n > 1

; loop...
for i=1,max_iter do begin

    ; define a new midpoint...
    mid = lo + 0.5*(hi-lo)

    if verbose then begin
        print, format='(i4,$)', i
        print, format='(g12.5,$)', lo
        print, format='(g12.5,$)', mid
        print, format='(g12.5,$)', hi
    endif

    ; can we quit?
    if( type_int and hi-lo lt 1) then begin
        if verbose then print, "stopped"
        return, mid
    endif

    ; stop if:
    ; 1) we're in fractional epsilon mode
    ; 2) the high and low limits have converged
    ; 3) the number of points between the high and low limits is less than 
    ;    the number defining the precision of the estimate
    if( epsilon le 0) then begin
        if verbose then begin
            print, format='(i9,$)', n_hi-n_lo
            print, format='(i9,$)', n_stop
            print, ""
        endif
        if hi-lo le 0 or n_hi - n_lo le n_stop then return, mid
    ; stop if:
    ; 1) we're in absolte epsilon mode
    ; 2) the difference between the high and low limits 
    ;    is less than epsilon
    endif else begin
        if verbose then begin
            print, format='(g12.5,$)', hi-lo
            print, format='(g12.5,$)', epsilon
            print, ""
        endif
        if hi - lo le epsilon then return, mid
    endelse

    ; count points below value of mid...
    buf = where( array lt mid, n_mid)

    ; shift the appropriate limit...
    if( n_mid eq n_target) then begin
        return, mid
    endif else if( n_mid gt n_target) then begin
        hi = mid
        n_hi = n_mid
    endif else begin
        lo = mid
        n_lo = n_mid
    endelse

endfor

message, 'maximum iterations exceeded', /continue
return, mid

end
