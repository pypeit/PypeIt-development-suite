;-----------------------------------------------------------------------
pro psland, psfile, ENCAPSULATED=encapsulated, COLOR=color
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
    device, /landscape
    ret = dialog_printersetup()
endif else begin
    set_plot, 'ps'
    device, filename=psfile, /landscape
endelse

if keyword_set(ENCAPSULATED) then begin
    device, /encapsul
endif

if keyword_set(COLOR) then begin
    device, bits_per_pixel=8, /color
endif

end
