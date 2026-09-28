;-----------------------------------------------------------------------
function pa, X, Y, T
;-----------------------------------------------------------------------
; Gregory D. Wirth / WMKO / 2005-Sep-13
;
; Computes the position angle [degress] of vector X w.r.t. Y, measured
; counter-clockwise from Y.  T is the vector from observer to target.
;
; X = vector defining direction (typically the slit direction)
; Y = vector defining PA=0 direction (typically toward north)
; T = vector from observer to target
; E = vector defining PA=90 direction 
;-----------------------------------------------------------------------

;; if X and Y are the same, then PA is zero by definition...
dp = dot(X,Y)
tiny = 1.d-5
if dp le tiny then begin
    pa_val = 0.d0
endif else if 180.d0-dp le tiny then begin
    pa_val = 180.d0
endif else begin

    ;; verify that X and Y are perp to T (these vectors should be in the
    ;; plane of the sky, which is tangent to T...
    A = cross(X,Y)
    B = cross(Y,X)
    tiny = 1.d-5
    if dot(A,T) gt tiny and dot(B,T) gt tiny then $
      message, 'vectors are not perpendiclar to T'

    ;; compute the dot product...
    pa_val = dot(X,Y)

    ;; correct dot product as needed...
    E = cross(Y,T)
    if dot(X,E) gt 90.d0 then pa_val = 360.d0 - pa_val
endelse
return, pa_val
end
