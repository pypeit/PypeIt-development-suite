;-----------------------------------------------------------------------
pro xlabel, string
;-----------------------------------------------------------------------
; Purpose:
;   Plot only the main label, a la Lick MONGO
;-----------------------------------------------------------------------
plot, [0,0], [1,1], /nodata, xstyle=4, ystyle=4, xtitle=string
end
