;------------------------------------------------------------------------
function skystat, array, mean, stddev, $
                  MAXITER=maxiter, MINDELTA=mindelta, TOPSIGMA=topsigma, $
                  LOWER=lower, VERBOSE=verbose
;------------------------------------------------------------------------
;+
; NAME:
;	SKYSTAT
;
; PURPOSE:
;	This function computes and returns an interatively-derived
;	estimate of the mean value and standard deviation of the
;	passed array.  The return value indicates whether the solution
;	converged.
;
; CATEGORY:
;	Data analysis
;
; CALLING SEQUENCE:
;	Result = SKYSTAT( array, mean, stddev)
;
; INPUTS:
;	array:	array for which the mean and standard deviation are to
;	be computed.
;
; KEYWORD PARAMETERS:
;	MAXITER: Set this keyword to the maximum number of iterations
;	allowed [default=10] 
;
;	MINDELTA: Set this keyword to the minimum fractional change in
;	the derived standard deviation requiring an additional
;	iteration. [default=0.1]
;
;	TOPSIGMA: Set this keyword to the number of standard
;	deviations above the mean value at which to reject
;	values. [default=3] 
;
;	LOWER: Set this keyword to the lower threshold value above
;	which values are legal [default=minimum value in array - 1]
;
;	VERBOSE: Set this flag to cause updates to be printed during
;	execution. [default=n]
;
; OUTPUTS:
;	This function returns "true" (1B) if the computation completed
;	in fewer than the maximum number of iterations, else "false"
;	(0B).
;
; OPTIONAL OUTPUTS:
;	MEAN: will contain the iterative estimate of the mean value in
;	the array.
;
;	STDDEV: will contain the iterative estimate of the standard
;	deviation in the array.
;
; PROCEDURE:
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2004-Jun-04	GDW	Original version
;-
;------------------------------------------------------------------------

; defaults...
if ~ keyword_set(maxiter) then maxiter = 10
if ~ keyword_set(mindelta) then mindelta = 0.1
if ~ keyword_set(topsigma) then topsigma = 3.
if ~ keyword_set(lower) then lower = min(array)-1.

; Initialize parameters...
status = 0B
oldstd = 0.
mindelta2 = mindelta^2

; First guess at image parameters...
good = where( array gt lower)
result = moment(array[good], sdev=stddev) 
mean = result[0]

; Now truncate the pixels at the requested number of sigmas
; above the sky, and estimate the sky parameters again.
; Repeat until converged or maximum number of interations...
for k=1,maxiter do begin
    toplimit = mean + topsigma*stddev
    good = where( array gt lower and array lt toplimit)
    result = moment(array[good], sdev=stddev) 
    mean = result[0]

    if keyword_set(verbose) then $
      print, format='("iter=", i, " sky=", f15.5, " sigma=", f15.5)', $
             k, mean, stddev

    ; error if the STDDEV is not legal...
    if stddev le 0 then message, 'STDDEV is not a positive value'

    delta = ((stddev-oldstd)/stddev)^2
    if delta lt mindelta2 then begin
        status = 1B             ; converged!
        break
    endif 

    ; save value for next iteration...
    oldstd = stddev
endfor 

return, status

end
