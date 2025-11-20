;-----------------------------------------------------------------------
pro tlabel, string, font=font, charsize=charsize
;-----------------------------------------------------------------------
; Purpose:
;   Plot only the main label, a la Lick MONGO
;-----------------------------------------------------------------------
save = !p
if keyword_set(font) then !p.font = font
if keyword_set(charsize) then !p.charsize = charsize
x1 = !x.crange[0]
x2 = !x.crange[1]
dx = x2-x1
y1 = !y.crange[0]
y2 = !y.crange[1]
dy = y2-y1
xyouts, x1+0.5*dx, y1+1.03*dy, string, alignment=0.5
!p = save
end
