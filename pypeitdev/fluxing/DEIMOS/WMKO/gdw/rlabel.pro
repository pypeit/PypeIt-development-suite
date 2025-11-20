;-----------------------------------------------------------------------
pro rlabel, string, CHARSIZE=charsize, FONT=font
;-----------------------------------------------------------------------
;+
; NAME:
;	RLABEL
;
; PURPOSE:
;	This procedure puts a label centered off the right-hand side of the
;	current plot, similar to the MONGO command 'rlabel'.
;
; CATEGORY:
;	Plotting.
;
; CALLING SEQUENCE:
;	RLABEL, string, CHARSIZE=charsize, FONT=font
;
; INPUTS:
;	string:	Text string to be plotted
;
; KEYWORD PARAMETERS:
;	CHARSIZE: character size to use for the label
;	FONT: font to use for the label
;
; OUTPUTS:
;	None.
;
; SIDE EFFECTS:
;	The specified string is added to the current plot.
;
; RESTRICTIONS:
;	None.
;
; EXAMPLE:
;	Add the label 'Some restrictions apply' to the right-hand side 
;	of the current plot.  Use TrueType fonts of twice normal size.
;		rlabel, 'Some restrictions apply', font=1, charsize=2
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2001-Mar-07	GDW	Original version
;-
;-----------------------------------------------------------------------

; store the current graphics environment...
save = !p
if keyword_set(charsize) then !p.charsize = charsize
if keyword_set(font) then !p.font = font

x_plot = !x.crange[0] + 1.025*(!x.crange[1]-!x.crange[0])
y_plot = !y.crange[0] + 0.5*(!y.crange[1]-!y.crange[0])

; plot the text...
xyouts, x_plot, y_plot, string, ALIGNMENT=0.5, ORIENT=-90

; restore graphics environment...
!p = save

end
