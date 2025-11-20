;-----------------------------------------------------------------------
pro spectral_mask, x, y, xgood, ygood, threshold=threshold
;-----------------------------------------------------------------------
; Purpose:
;       Given a spectrum (x,y) and a threshold (thresh) defining the
;       minimum good pixel value (as a fraction of the peak value),
;       return the pixels which exceed this threshhold.
;-----------------------------------------------------------------------

;; convert fractional threshold to a value...
cut = threshold * max(y)

;; select pixels above threshold...
good = where( y ge cut, count)
if count le 0 then message, 'ERROR: no pixels above threshold -- abort!'

;; return good spectrum...
xgood = x[good]
ygood = y[good]

end
