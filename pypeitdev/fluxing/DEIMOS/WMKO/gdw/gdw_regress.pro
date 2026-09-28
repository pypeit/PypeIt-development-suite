;-----------------------------------------------------------------------
function gdw_regress, xbuf, ybuf, siglimit, STATUS=status, YFIT=yfit, $
  REJECTED=rejected
;-----------------------------------------------------------------------
; iterativelly run REGRESS with sigma rejection until it stabilizes
;-----------------------------------------------------------------------

x = xbuf
y = ybuf
n = n_elements(x)
keepgoing = 1B
rejected = bytarr(n) + 1B
good2 = where( x eq x)

;; loop...
while keepgoing do begin

    ;; fit line to data...
    m = REGRESS( X, Y, STATUS=status, YFIT=yfit2, CONST=b)

    ;; compute offsets from fit...
    ydiff = y - yfit2

    ;; determine variation...
    stmad = mad( ydiff, /sigma)

    ;; determine offset in sigma...
    ysigma = abs( ydiff/stmad)

    ;; determine number of bad points...
    bad = where( ysigma gt siglimit, nbad)
    if nbad eq 0 then keepgoing = 0B

    ;; retain good points...
    good = where( ysigma lt siglimit, ngood)
    if ngood lt 1 then begin
        message, 'no good data!', /info
        status = 4
        return, !values.f_nan
    endif 
    x = x[good]
    y = y[good]
    good2 = good2[good]

endwhile 

;; disable rejection for good data points...
rejected[good2] = 0B

yfit = m[0]*xbuf + b
return, m

end 
