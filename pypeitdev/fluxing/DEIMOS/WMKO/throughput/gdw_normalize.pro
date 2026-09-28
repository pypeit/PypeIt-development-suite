;-----------------------------------------------------------------------
function gdw_normalize, x, y
;-----------------------------------------------------------------------

;; fit background with B-spline...
nord = 4
fullbkpt = bspline_bkpts( x, nord=nord, nbkpts=15)
sset = bspline_iterfit( x, y, nord=nord, maxiter=0, yfit=yfit, fullbkpt=fullbkpt)

;; normalize...
norm = y - yfit
return, norm

end
