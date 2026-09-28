;-----------------------------------------------------------------------
Function MAD, Array, SIGMA=sigma, MEDIAN=median
;-----------------------------------------------------------------------
;+
; NAME:
;	MAD
;
; PURPOSE:
;	Compute and return the Median Absolute Deviation of the passed
;	array.  Note that MAD/0.6745 is comparable to the standard
;	deviation for a random normal distribution.  This routine can
;	optionally return the MAD/0.6745 value.
;
; CATEGORY:
;	Statistics
;
; CALLING SEQUENCE:
;	Result = MAD( array)
;
; INPUTS:
;	Array:	A real-valued array for which the MAD is to be computed.
;
; KEYWORD PARAMETERS:
;	SIGMA:	Set this keyword to return MAD/0.6745, a robust
;	estimate of the standard deviation for the array.  The default 
;	is to return MAD.
;
;	MEDIAN: Set this keyword to a named variable to receive the
;	median value of the Array.  THis is provided as a convenience
;	so that the median does not need to be computed separately.
;
; OUTPUTS:
;	This function returns a scalar value which is either MAD (by
;	default) or MAD/0.6745 (if requested by use of the /SIGMA keyword)
;
; PROCEDURE:
;	Computes the MAD via
;		median( abs( Array - median( Array)))
;
; EXAMPLE:
;	1. To compute the median absolute deviation of array X:
;		value = mad(X)
;
;	2. To compute a mad-based estimate of the standard deviation
;	for the array X:
;		sigma = mad(X,/SIGMA)
;
; REFERENCE:
;	Beers, Flynn, & Gebhardt, 1990, AJ, 100, 32
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2002-Jan-18	GDW	Original version
;	2004-Jun-04	GDW	Added MEDIAN keyword
;-
;-----------------------------------------------------------------------

if keyword_set(SIGMA) then begin
    factor = 0.6745
endif else begin
    factor = 1.0
endelse

median = median( Array, /EVEN)
return, median( abs( Array - median))/factor

end
