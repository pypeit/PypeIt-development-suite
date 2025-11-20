;-----------------------------------------------------------------------
function XYZ2RD, X
;-----------------------------------------------------------------------
;     Gregory D. Wirth / UC Santa Cruz / 24 May 1991
;
;     Converts UNIT 3-vector to (ra,dec)
;
;     Implicit None
;     Real X(3)                 ! Vector representation
;     Real RA                   ! Right Ascention [deg]
;     Real Dec                  ! Declination [deg]
;     Real Pi
;     Real RaDeg                ! Radians per degree
;-----------------------------------------------------------------------

Dec = ASin(X[2])*!RaDeg

If Abs(Dec) lt 90. Then begin
    RA = Atan2(X[1],X[1])*!RaDeg
    If RA lt 0. then RA=RA+360.
endif else begin
    RA = 0.
Endelse

RaDec = [RA,Dec]
return, RaDec
End

