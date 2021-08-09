;-----------------------------------------------------------------------
function cross, X, Y
;-----------------------------------------------------------------------
;     Gregory D. Wirth / UC Santa Cruz / 24 May 1991
;
;     Computes the vector cross product of the unit vectors X and Y.
;
;      Implicit None
;      Real X(3)                 ! First vector in cross product
;      Real Y(3)                 ! Second vector in cross product
;      Real Z(3)                 ! Resultant vector in cross product
;-----------------------------------------------------------------------

VNORM,X
VNORM,Y

Z = Y
Z[0] = X[1]*Y[2] - Y[1]*X[2]
Z[1] = X[2]*Y[0] - Y[2]*X[0]
Z[2] = X[0]*Y[1] - Y[0]*X[1]

;; Vnormalize result...
VNORM, Z
return, Z
end

