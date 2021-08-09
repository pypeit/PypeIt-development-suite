;-----------------------------------------------------------------------
pro parabola, x, a, f, pder
;-----------------------------------------------------------------------
; parabola / GDW / 14 Jan 1999
;
; Purpose:
;   Compute the parabola of the form y=a[0]*(x-a[2])**2+a[1], plus
;   partial derivatives at all of the points in the x array.
;
; Passed parameters:
;   real   x                     I: independent variable
;   real   a[3]                  I: array of parameters
;   real   f                     O: function value at x
;   real   pder[n_elements(x),3] O: partial derivatives at x[1]
;-----------------------------------------------------------------------

; evaluate function...
f = a[0]*(x-a[2])^2+a[1]

; evaluate partials...
if n_params(0) le 3 then return
pder = fltarr(n_elements(x), 3)
pder[0,0] = (x-a[2])^2
pder[*,1] = 1.
pder[0,2] = -2*a[0]*(x-a[2])

end

