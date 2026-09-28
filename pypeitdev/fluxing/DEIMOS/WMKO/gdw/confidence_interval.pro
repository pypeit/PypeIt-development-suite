;-----------------------------------------------------------------------
function confidence_interval, array, WIDTH=width
;-----------------------------------------------------------------------
; conf_intrvl.awk / GDW / 22 Jun 95
;
; Returns a the half-width of the symmetrical interval spanning
; the requested fraction of the data (defaults to 68%).  When I say
; symmetrical, I mean that the endpoints of the confidence interval are
; equally distant from the median value, rather than there necessarily
; being an equal number of points above and below the median within
; the stated range.
;-----------------------------------------------------------------------

if( not keyword_set( width)) then width = 0.68 ; default

; sort data...
x = array( sort( array))

; count data...
nr = n_elements(x)

; compute the equal-point interval...
npts = fix( width*nr+0.5)
i = fix((0.5 - width/2)*nr+0.5)
j = i + npts - 1
k = fix(nr/2)+1

; compute the median...
if nr mod 2 eq 1 then begin
    median = x[k-1] 
endif else begin
    median = 0.5*(x[k-1]+x[k])
    radius = 0.5*(x[j]-x[i])
endelse

; compute the symmetrical interval.  The number of points we must
; have in out confidence interval is "npts."  We begin at the first
; point, and we compute two things: d1 gives the distance from
; this point to the median value, and d2 gives the distance from
; the median to the other end of the interval containing "npts"
; points.  We are trying to balance d1 and d2, hence making the
; interval symmetrical.  The "delta" parameter measures how
; nearly balanced d1 and d2 are.  We choose the endpoints which
; give the minimum value of delta, "deltamin".
deltamin = 1.e32
for i=0,nr-npts do begin
    d1 = median - x[i]
    j = i + npts - 1
    d2 = x[j] - median
    delta = (d2-d1)^2
    if delta lt deltamin then begin 
        deltamin = delta
        top = j
        bot = i
    endif
endfor
radius = 0.5*(x[top]-x[bot])

; double-check: count the point within the range:
;;; ngood = 0
;;; for i=0,nr-1 do begin
;;;    if x[i] ge x[bot] and x[i] le x[top] then ngood = ngood + 1
;;; endfor

;;; print, median, radius, fix(100.*ngood/nr+0.5), "%"

return, radius
end

