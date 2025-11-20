;-----------------------------------------------------------------------
function fit_parabola, telfocus, fwhm, nimages, bestfocus, bestfwhm, yfit
;-----------------------------------------------------------------------
; fit_parabola / GDW / 14 Jan 1999
;
; Purpose:
;   Estimate the best focus value by analyzing the passed focus
;   images.
;
; Passed parameters:
;   real   telfocus[nimages]     I: array of focus settings for images
;   real   fwhm[nimages]         I: array of measured image sizes
;   int    nimages               I: number of images in image array
;   real   bestfocus             O: optimal focus setting
;   real   bestfwhm              O: optimal size value
;   real   yfit[nimages]         O: array of fitted FWHM at telfocus[i]
;-----------------------------------------------------------------------

; find the best of the focus values...
big = 1.e32
min_fwhm = big                  ; minimum FWHM
min_fwhm_focus = big            ; optimal focus
for i = 0,nimages-1 do begin
    if( fwhm[i] lt min_fwhm) then begin
        min_fwhm = fwhm[i]
        min_fwhm_focus = telfocus[i]
    endif
endfor

; define initial guesses...
a = [1., min_fwhm, min_fwhm_focus]

; set all weights to 1.
weights = replicate( 1.0, n_elements(telfocus))

; fit the parabola
yfit = gdwcurvefit( telfocus, fwhm, weights, a,  $
                    CHISQ=chisq, FUNCTION_NAME='parabola',  $
                    ITER=n_iter,STATUS=status)

if status then begin
    bestfwhm = a[1]
    bestfocus = a[2]
endif else begin
    bestfwhm = min_fwhm
    bestfocus = min_fwhm_focus
endelse

return, status

end
