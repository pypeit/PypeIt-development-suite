;-----------------------------------------------------------------------
pro display, image, frame, NOERASE=noerase, NOFILL=nofill, $
             XSIZE=xsize, YSIZE=ysize, $
             NOZSCALE=nozscale, CONTRAST=contrast, NOZRANGE=nozrange, $
             Z1=z1, Z2=z2, LOGARITHMIC=logarithmic, VERBOSE=verbose
;-----------------------------------------------------------------------
;+
; NAME:
;	DISPLAY
;
; PURPOSE:
;	This procedure renders the specified image to the specified
;	window, or the current one if not specified.
;
; CATEGORY:
;	Image processing
;
; CALLING SEQUENCE:
;	DISPLAY, Image, Frame
;
; INPUTS:
;	Image:	2-D array to display.
;
; OPTIONAL INPUTS:
;	Frame:	Index of the frame number on which to display.
;	
; KEYWORD PARAMETERS:
;	NOERASE:	erase frame before display?
;	NOFILL:		scale image to fit display window?
;	NOZSCALE:	display range of greylevels near median?
;	CONTRAST:	contrast adjustment for zscale algorithm [def=0.25]
;	NOZRANGE:	display full image intensity range?
;	Z1:		minimum greylevel to be displayed
;	Z2:		maximum greylevel to be displayed
;	LOGARITHMIC:	use logarithmic scale for color mapping?
;	VERBOSE:	print diagnostic messages?
;	XSIZE:		columns in window to create
;	YSIZE:		rows in window to create
;
; OUTPUTS:
;
; SIDE EFFECTS:
;	The image is rendered to the selected window.  If not
;	existing, the window is created.
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
;	1999-Nov-4	GDW	Original version
;-
;-----------------------------------------------------------------------

def_size = 512

; get keywords...
do_erase = not keyword_set(NOERASE)
fill = not keyword_set(NOFILL)
if( not keyword_set(XSIZE)) then xsize=def_size
if( not keyword_set(YSIZE)) then ysize=def_size
zscale = not keyword_set(NOZSCALE)
if( not keyword_set(CONTRAST)) then contrast=0.25
zrange = not keyword_set(NOZRANGE)
logarithmic = keyword_set(LOGARITHMIC)
verbose = keyword_set(VERBOSE)

; check args...
np = n_params()
if( np lt 1 or np gt 2) then begin
    message, 'Wrong number of parameters passed'
endif

; verify array...
buf = size( image)
n_dim = buf[0]
nx    = buf[1]
ny    = buf[2]
if( n_dim ne 2) then begin
    message, 'Image not 2-dimensional'
endif

if( contrast lt 0 or contrast gt 1) then begin
    message, 'CONTRAST parameter out of legal range 0-1.'
endif

; use the current device if none specified...
if( np lt 2) then begin
    frame = !d.window
endif

; check whether this window exists...
device, window_state=wexist
if( wexist[frame]) then begin
    wset, frame
    xsize = !d.x_size
    ysize = !d.y_size
endif else begin
    window, frame, xsize=xsize, ysize=ysize
endelse

; resample the image as needed...
if fill then begin
    display_image = congrid( image, xsize, ysize, /interp)
    xstart = 0
    ystart = 0
endif else begin

    if( nx le xsize) then begin
        ix1 = 0
        ix2 = nx - 1
        xstart = nint(0.5*(xsize-nx))
    endif else begin
        ixc = nint(0.5*nx)
        ix1 = ixc - nint(0.5*xsize)
        ix2 = ix1 + xsize - 1
        xstart = 0
    endelse

    if( ny le ysize) then begin
        iy1 = 0
        iy2 = ny - 1
        ystart = nint(0.5*(ysize-ny))
    endif else begin
        iyc = nint(0.5*ny)
        iy1 = iyc - nint(0.5*ysize)
        iy2 = iy1 + ysize - 1
        ystart = 0
    endelse

    display_image = image[ix1:ix2,iy1:iy2]

endelse

; determine lower and upper image scaling limits based on user settings...
if zscale then begin
    cdf_lo = 0.5*contrast
    lo = cdf_frac( image, cdf_lo, -0.05)

    cdf_hi = 1. - cdf_lo
    hi = cdf_frac( image, cdf_hi, -0.05)
    if verbose then print, 'zscale contrast=', contrast,' lo=', lo,' hi=', hi
endif else if zrange then begin
    lo = min( image, max=hi)
    if verbose then print, 'zrange lo=', lo,' hi=', hi
endif else begin
    lo = z1
    hi = z2
    if verbose then print, 'manual lo=', lo,' hi=', hi
endelse

; scale image appropriately...
if logarithmic then begin
    display_image = alog(temporary(display_image) - lo + 1. > 1.)
    lo = 0
    hi = alog(hi - lo + 1.)
    if verbose then print, 'ztrans=log lo=', lo,' hi=', hi,  $
      ' min=', min(display_image), $
      ' max=', max(display_image)
endif else begin
    if verbose then print, 'ztrans=linear'
endelse

if do_erase then erase
tv, bytscl( display_image, min=lo, max=hi), xstart, ystart

end
