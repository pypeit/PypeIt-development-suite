;-----------------------------------------------------------------------
pro VNORM, Z
;-----------------------------------------------------------------------
;     normalizes vector to unit length.
;
;      Implicit None
;      Real Z(3)                 ! Three-vector
;      Real Zlen                 ! Length of Z vector
;      Integer i                 ! counter
;-----------------------------------------------------------------------
Zlen = Sqrt(Total(Z^2))
if Zlen le 0. then begin
    message, 'zero length vector'
endif else begin
    Z = Z/Zlen
endelse
End

