;-----------------------------------------------------------------------
function gdwcolor, name
;-----------------------------------------------------------------------
;+
; NAME:
;	GDWCOLOR
;
; PURPOSE:
;	This function adds a bunch of colors to the current color
;	table and returns the index correponding to the requested
;	color.
;
; CATEGORY:
;	Graphics
;
; CALLING SEQUENCE:
;       value = GDWCOLOR( name)
;
; INPUTS:
;	Name:	Name of the color for which to return the color
;               index.   Options are: black, white, red, lime, blue,
;               yellow, cyan, magenta, gray, olive, purple, teal, 
;               navy, green, maroon
;
; OUTPUTS:
;	Returns the color index corresponding to the named color.
;
; COMMON BLOCKS:
;	GDWCOLORLIST:	contains color table.  For internal use only.
;
; SIDE EFFECTS:
;	Modifies the color table.
;
; RESTRICTIONS:
;       The following commands should be executed to configure IDL
;       colors before using this routine:
;               device, true=24
;               device, decomp=0
;
; EXAMPLE:
;       1) plot something in navy blue
;       plot, x, y, color=gdwcolor('navy')
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	2010-Jun-26	GDW	Original version
;       2011-Nov-16     GDW     Add turquoise
;-
;-----------------------------------------------------------------------

common gdwcolorlist, color_table

;; initialize color table...
if n_elements(color_table) eq 0 then begin

    ;; define data...
    foo = {COLOR, name:'', red:0, green:0, blue:0, index:0}
    maxcolors = 100
    color_table = replicate({COLOR}, maxcolors)
    
    i = -1
    max_colors = 255
    color_table[++i] = { COLOR, name:'black',    red:000, green:000, blue:000, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'white',    red:255, green:255, blue:255, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'red',      red:255, green:000, blue:000, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'lime',     red:000, green:255, blue:000, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'blue',     red:000, green:000, blue:255, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'yellow',   red:255, green:211, blue:000, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'cyan',     red:000, green:255, blue:255, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'magenta',  red:255, green:000, blue:255, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'darkgray', red:0, green:0, blue:0, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'gray',     red:127, green:127, blue:127, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'lightgray',red:160, green:160, blue:160, index:max_colors-i-1}
;;    color_table[++i] = { COLOR, name:'lightgray',red:200, green:200, blue:200, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'olive',    red:127, green:127, blue:000, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'purple',   red:127, green:000, blue:127, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'turquoise',red:64,  green:224, blue:208, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'navy',     red:000, green:000, blue:127, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'green',    red:000, green:159, blue:107, index:max_colors-i-1}
    color_table[++i] = { COLOR, name:'maroon',   red:127, green:000, blue:000, index:max_colors-i-1}

    ncolors = i

    ;; loop over colors and load table...
    for i=0,ncolors-1 do begin
        c = color_table[i]
        tvlct, c.red, c.green, c.blue, c.index
    endfor

endif

;; locate color...
good = where( color_table.name eq name, n)
if n ne 1 then message, 'invalid color: '+name
return, color_table[good].index

end 
