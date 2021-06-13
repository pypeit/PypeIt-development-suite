;-----------------------------------------------------------------------
function dot, X, Y
;-----------------------------------------------------------------------
;     Gregory D. Wirth / UC Santa Cruz / 24 May 1991
;
;     Calculates the angle between the UNIT 3-vectors X and Y.
;     Returns this angle Theta in degrees.
;
;      Real X(3)                 ! first unit vector
;      Real Y(3)                 ! first unit vector
;      Real Theta                ! angular separation
;      Real dp                   ! scalar dot product
;      Real Xlen                 ! length of the X vector
;      Real Ylen                 ! length of the Y vector
;      Real Pi                   ! pi
;      Real RaDeg                ! radians per degree
;-----------------------------------------------------------------------

dp = total(X*Y)

;; If abs(dp) exceeds 1, then we probably forgot to normalize the
;; unit vectors before calling DOT...

if abs(dp) gt 1 then begin
    
    xlen = Sqrt(Total(X^2))
    ylen = Sqrt(Total(Y^2))
    tiny = 1.e-12

    ;; First check whether the vectors are normalized...
    if xlen gt 1.d0+tiny or ylen gt 1.d0+tiny then begin
        message, "ERROR in DOT --- Vectors not normalized!"
        
        ;; Otherwise, chalk it up to roundoff error.  Reset to 1 and continue...
    endif else begin
        dp = dp/abs(dp)
    endelse
endIf

Theta = Acos(dp)*!radeg
return, theta

End

