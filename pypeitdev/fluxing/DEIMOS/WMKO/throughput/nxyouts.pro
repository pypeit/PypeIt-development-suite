;-----------------------------------------------------------------------
pro nxyouts, xnorm, ynorm, text, _EXTRA=ex
;-----------------------------------------------------------------------
;+
; NAME:
;	NXYOUTS
;
; PURPOSE:
;	This procedure will generate text output on the current plot
;	window in normalized coordinated.  Unlike xyouts, it will work
;	when !p.multi is in effect.
;
; CATEGORY:
;	Graphics.
;
; CALLING SEQUENCE:
;	nxyouts, xnorm, ynorm, text
;
; INPUTS:
;	xnorm:  normalized x position (0=left axis, 1=right axis)
;	ynorm:  normalized y position (0=bottom axis, 1=top axis)
;	text:   string to plot
;
; OUTPUTS:
;	This procedure will generate text output on the current plot
;	window in normalized coordinated.  Unlike xyouts, it will work
;	when !p.multi is in effect.
;
; EXAMPLE:
;       1) Put the string 'Have a nice day' in the center of the plot:
;       nxyouts, 0.5, 0.5, 'Have a nice day', align=0.5
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2011-Jan-20	GDW	Original version
;-
;-----------------------------------------------------------------------

x1 = !x.crange[1]
x0 = !x.crange[0]
dx = x1-x0

y1 = !y.crange[1]
y0 = !y.crange[0]
dy = y1-y0

x = x0 + xnorm*dx
y = y0 + ynorm*dy

xyouts, x, y, text, _EXTRA=ex

end
