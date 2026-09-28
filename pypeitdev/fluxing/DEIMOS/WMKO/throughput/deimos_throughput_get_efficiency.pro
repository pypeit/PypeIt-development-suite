;-----------------------------------------------------------------------
pro deimos_throughput_get_efficiency, wav, eff, params, VERBOSE=verbose
;-----------------------------------------------------------------------
;+
; NAME:
;	DEIMOS_THROUGHPUT_GET_EFFICIENCY
;
; PURPOSE:
;	This routine will compute the average efficiency over the 
;       specified wavelength range.
;
; CATEGORY:
;	Instrument performance monitoring
;
; CALLING SEQUENCE:
;	DEIMOS_THROUGHPUT_GET_EFFICIENCY, params
;
; INPUTS:
;       wav:    array of wavelengths
;       eff:    array of efficiencies
;	params: structure storing the parameter set, including:
;               lambda_eff  [I] - passband center
;               dlambda_eff [I] - passband width
;               efficiency  [O] - median efficiency within passband
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2000-Nov-16	GDW	Original version
;-
;-----------------------------------------------------------------------

;; loop over passbands...
n = n_elements( params)

for i=0,n-1 do begin

    ;; initial value is undefined...
    params[i].efficiency = !values.f_nan

    ;; generate endpoints for wavelength range...
    lambda1 = params[i].lambda_eff - 0.5 * params[i].dlambda_eff
    lambda2 = lambda1 + params[i].dlambda_eff

    ;; determine whether the passband is entirely contained within spectrum...
    good = where(wav le lambda1, count1)
    good = where(wav ge lambda2, count2)
    good = where(wav ge lambda1 and wav le lambda2, count3)
    if count1 lt 1 || count2 lt 1 || count3 lt 1 then begin
        if keyword_set(VERBOSE) then begin
            errmsg = 'passband at ' + stringify(params[i].lambda_eff) $
                  + ' not entirely contained within spectrum for dataset'
            message, errmsg, /info
        endif
        continue
    endif

    ;; get pixels in passband...
    good = where(wav ge lambda1 $
                 and wav le lambda2 $
                 and finite(eff), count)

    ;; skip if no good pixels...
    if count lt 1 then begin
        if keyword_set(VERBOSE) then begin
            errmsg = 'no good pixels in passband at ' $
                  + stringify(params[i].lambda_eff)
            message, errmsg , /info
        endif
        continue
    endif

    ;; take median efficiency within region...
    params[i].efficiency = median(eff[good])

endfor 

end
