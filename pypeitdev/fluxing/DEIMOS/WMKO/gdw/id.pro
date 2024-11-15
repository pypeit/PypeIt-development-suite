;-----------------------------------------------------------------------
pro id, string, font=font, charsize=charsize
;-----------------------------------------------------------------------
if not keyword_set(font) then font=!p.font
if not keyword_set(charsize) then charsize=!p.charsize
buf = systime(0)
if( n_params() gt 0) then buf = string + '  ' + buf
if font eq 1 then buf = '!11' + buf + '!3'
x1 = !x.crange[0]
x2 = !x.crange[1]
dx = x2-x1
y1 = !y.crange[0]
y2 = !y.crange[1]
dy = y2-y1
xyouts, x2, y1 + 1.005*dy, buf, alignment=1.0, charsize=charsize, font=font
end
