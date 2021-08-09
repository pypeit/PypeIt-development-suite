;-----------------------------------------------------------------------
function RD2XYZ, RA, Dec
;-----------------------------------------------------------------------
;     Gregory D. Wirth / UC Santa Cruz / 24 May 1991
;
;     Converts (ra,dec) to a UNIT 3-vector.
;
;     NOTE! RA and Dec are expected in decimal DEGREES, so if you have
;     RA in decimal hours, pass RA*15. instead.
;
;      Implicit None
;      Real RA                   ! Right Ascention [deg]
;      Real Dec                  ! Declination [deg]
;      Real X(3)                 ! Vector representation
;      Real Pi
;      Real RaDeg                ! Radians per degree
;-----------------------------------------------------------------------

X = dblarr(3)
X[0] = Cos(ra/!radeg)*Cos(dec/!radeg)
X[1] = Sin(ra/!radeg)*Cos(dec/!radeg)
X[2] = Sin(dec/!radeg)
return, X
End

