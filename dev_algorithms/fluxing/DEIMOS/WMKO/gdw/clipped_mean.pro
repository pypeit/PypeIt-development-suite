;-----------------------------------------------------------------------
function clipped_mean, array, mean, sigma
;-----------------------------------------------------------------------
;+
; NAME:
;	CLIPPED_MEAN
;
; PURPOSE:
;	This function uses iterative sigma clipping to estimate the
;	mean and standard deviation of the "background" pixels in an image.
;
; CATEGORY:
;	Data analyis
;
; CALLING SEQUENCE:
;	Result = CLIPPED_MEAN( Array, Mean, Sigma)
;
; INPUTS:
;	Array:	Image data
;
; OPTIONAL INPUTS:
;	
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	This function returns 'true' (1L) on success and 'false' (0L)
;	on failure.
;
; OPTIONAL OUTPUTS:
;	Describe optional outputs here.  If the routine doesn't have any, 
;	just delete this section.
;
; PROCEDURE:
;
; EXAMPLE:
;	Status = CLIPPED_MEAN( Array, Mean, Sigma)
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	1999-Jul-24	GDW	Original version
;-
;
; LOCAL VARIABLES:
;    real(2) result     array containing results from moment()
;    real    mean       average pixel value of sample
;    real    sigma      standard deviation of sample
;    real    toplimit   maximum pixel value allowed in sample
;    real    topsigma   number of sigma above mean to clip
;    real    oldsigma   value of sigma on previous iteration
;    real    delta2     squared fraction change in sigma
;    real    mindelta   minimum fraction change in stddev for
;                       more iterations 
;    real    mindelta2  squared value of mindelta
;    int     maxiter    maximum number of slipping iterations
;    int     k          iteration counter
;
; INTERNAL PROCEDURES USED:
;    moment, where
;-----------------------------------------------------------------------

debug = keyword_set( DEBUG)

; define constants...
topsigma = 3                    ; sigma clipping factor above mean
maxiter  = 10                   ; maximum number of clipping iterations
mindelta = 0.1                  ; minimum fraction change in stddev for 
                                ; continuing iterations

; initialize...
oldsigma = 0.
mindelta2 = mindelta^2

; Now truncate the pixels at the requested number of sigmas
; above the sky, and estimate the sky parameters again.
; Repeat until converged or maximum number of interations...
for k=1,maxiter do begin

    if( k eq 1) then begin
        result = moment( array)
    endif else begin
        toplimit = mean + topsigma*sigma
        result = moment( array[ where( array le toplimit)])
    endelse

    ; check for illegal result values...
    if( result[1] le 0) then begin
        message, 'First moment of pixel distribution is <= 0', /continue
        return, 0L
    endif

    ; extract important field from result array...
    mean = result[0]
    sigma = sqrt(result[1])

    ; compute the fractional change in sigma...
    delta2 = ((sigma-oldsigma)/sigma)^2

    ; give feedback...
    if( debug) then $
      print, "mean=", mean, " sigma=", sigma, " toplimit=", toplimit, $
      " delta2=", delta2

    ; quit if the change was not significant...
    if( delta2 lt mindelta2 ) then return, 1L

    ; store value for comparison on next iteration...
    oldsigma = sigma
endfor

message, 'Maximum iterations exceeded', /continue
return, 0L

end
