;-----------------------------------------------------------------------
function wmean, data, sigma, sigma_mean
;-----------------------------------------------------------------------
;+
; NAME:
;	WMEAN
;
; PURPOSE:
;	This function returns the weighted mean value of the passed
;	values with uncertainties, and optionally the estimated
;	uncertainty of the mean.
;
; CATEGORY:
;	Statistics
;
; CALLING SEQUENCE:
;	Result = WMEAN( Data, Sigma, Sigma_Mean)
;
; INPUTS:
;	Data:	Array of values to average
;	Sigma: Array of cooresponding uncertainties
;
; OUTPUTS:
;	This function returns the weighted mean value of the passed array.
;
; OPTIONAL OUTPUTS:
;	Sigma_Mean:	The estimates uncertainty in the mean.
;
; EXAMPLE:
;	Compute the weighted mean value of the array X with
;	corresponding uncertainties Sigma:
;		Mean = WMEAN( X, Sigma, Sigma_mean)
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	1999-Nov-5	GDW	Original version
;-
;-----------------------------------------------------------------------

if( n_params() lt 2) then begin
    message, 'Missing parameter(s)'
endif

; count size of arrays...
n = n_elements(data)
n2 = n_elements(sigma)

if ( n2 ne n ) then begin
    message, 'data and sigma arrays are not the same size'
endif

; intialize...
sum1 = 0
sum2 = 0

; loop...
for i=0,n-1 do begin
    sum1 = sum1 + data[i]/sigma[i]^2
    sum2 = sum2 + 1./sigma[i]^2
end

result = sum1/sum2
sigma_mean = sqrt(1./sum2)
return, result
end


