pro psclose
;+
; NAME:
;PSCLOSE -- close ps, open X.
;     
; PURPOSE:
;       To close the Postscript device and set the graphics output
;       device back to X Windows.
;     
; CALLING SEQUENCE:
;       PSCLOSE
;     
; INPUTS:
;       None.
;     
; OUTPUTS:
;       None.
;
; KEYWORDS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       The device is changed.
;
; RESTRICTIONS:
;       A PostScript file must be open.
;
; RELATED PROCEDURES:
;       PSOPEN
;
; MODIFICATION HISTORY:
;       Written by Tim Robishaw in ancient times.
;-

; MAKE SURE WE HAVE POSTSCRIPT DEVICE OPEN...
if (!d.name ne 'PS') then begin
    message, 'DEVICE is not set to PS!', /INFO
    return
endif

; CLOSE THE POSTSCRIPT DEVICE...
device, /CLOSE_FILE

; SET THE GRAPHICS OUTPUT DEVICE TO X WINDOWS...
set_plot, 'X'

end; psclose
