;-----------------------------------------------------------------------
pro putlabel, x, y, n, string, CHARSIZE=charsize, NORMAL=normal
;-----------------------------------------------------------------------
; Purpose: Annotate the current plot.
;-----------------------------------------------------------------------

; store the current graphics environment...
save = !p
if keyword_set(charsize) then !p.charsize = charsize

; plot at either the normalized or data coordinates...
if keyword_set( NORMAL) then begin
    x_plot = !x.crange[0] + x*(!x.crange[1]-!x.crange[0])
    y_plot = !y.crange[0] + y*(!y.crange[1]-!y.crange[0])
endif else begin
    x_plot = x
    y_plot = y
endelse

; set alignment based on n...
alignment = 0.5*((n-1) mod 3)
;;print, 'alignment=', alignment, ' xplot=', x_plot, ' yplot=', y_plot

; plot the text...
xyouts, x_plot, y_plot, string, ALIGNMENT=alignment

; restore graphics environment...
!p = save

end

