;-----------------------------------------------------------------------
pro ptype, n, istyle
;-----------------------------------------------------------------------
; Purpose:
;   Define point types as in Lick MONGO.
;
; Arguments:
;   N = number of vertices
;   Istyle = 0 for open symbols (vertices connected)
;            1 for skeletal symbols (vertices connected to center)
;            2 for stellated symbols
;            3 for filled symbols (with current color)
;
; Note:
;   - Define symbol size using SYMSIZE keyword to PLOT and OPLOT
;
; Usage:
;   ptype, 3, 3 ; create filled triangles
;   plot, x, y, symsize=2
;-----------------------------------------------------------------------

if keyword_set( istyle) then $
  style = istyle $
else $
  style = 0

scale = 0.2
!p.psym = 8
if n gt 20 then n=20
theta = !pi*2./n
theta0 = 1.5*!pi - 0.5*theta

; use dot if fewer than 2 vertices...
if n le 1 then begin
    usersym, [0,0], [0,0]
endif else if style eq 0 then begin ; open symbols...
    a = findgen(n+1) * theta + theta0
    usersym, cos(a), sin(a)
endif else if style eq 1 then begin  ; skeletal symbols...
    i = findgen(2*(n+1))
    a =  i * (!pi/n) + theta0
    x = cos(a)
    y = sin(a)
    odd = where( i mod 2 eq 1)
    x[odd] = 0.
    y[odd] = 0.
    usersym, x, y
endif else if style eq 2 then begin  ; stellated symbols...
    i = findgen(2*(n+1))
    a =  i * (!pi/n) + theta0
    x = cos(a)
    y = sin(a)
    odd = where( i mod 2 eq 1)
    x[odd] = scale*x[odd]
    y[odd] = scale*y[odd]
    usersym, x, y
endif else if style eq 3 then begin ; filled symbols...
    a = findgen(n+1) * theta + theta0
    usersym, cos(a), sin(a), /fill
endif else begin
    message, 'Illegal value for style'
endelse
end
