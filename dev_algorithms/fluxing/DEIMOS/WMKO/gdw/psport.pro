;-----------------------------------------------------------------------
pro psport, psfile
;-----------------------------------------------------------------------
; psland / GDW / 19 Jan 1999
;
; Purpose:
;   Set up printing for PostScript output to either:
;      1) the postscript printer
;      2) a named postscript file
;
; Note:
;   Output will NOT be sent to the printer unless/until you close the
;   device using "device, /close".
;-----------------------------------------------------------------------

if n_params(0) lt 1 then begin
    set_plot, 'printer'
    device, /portait, xsize=7, ysize=10, xoff=0.75, yoff=0.5, /inch
    ret = dialog_printersetup()
endif else begin
    set_plot, 'ps'
    device, filename=psfile, /portrait, xsize=7, ysize=10, xoff=0.75, yoff=0.5, /inch
endelse

end
