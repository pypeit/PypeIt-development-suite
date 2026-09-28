;-----------------------------------------------------------------------
function itersmooth, y, WIDTH=width, MAXITER=maxiter, NSIGMA=nsigma, GROW=grow, REJECT=reject
;-----------------------------------------------------------------------
; Iteratively
;-----------------------------------------------------------------------

;; perform initial smoothing...
ysmooth = median( y, width)

;; iterate as requested...
for i=0,maxiter do begin

    ;; compute difference...
    diff = y - ysmooth
    sigma = mad( diff, /sigma)

    ;; locate points more than nsigma from the line...
    diff = abs(y - ysmooth)
    bad = where( diff gt nsigma*sigma, nbad)
    good = where( diff le nsigma*sigma, ngood)

    ;; we are done if no bad points exist...
    if nbad eq 0 then break



endfor

;; return smoothness...
return, ysmooth
end

