;-----------------------------------------------------------------------
pro debias, array, xdata1, xdata2, xbias1, xbias2, y1, y2,  $
            ORDER=order, VERBOSE=verbose, DISPLAY=display, $
            INTERACTIVE=interactive
;-----------------------------------------------------------------------
; Gregory D. Wirth / W. M. Keck Observatory / 17 Mar 1999
;
; Purpose:
;   Fit the overscan region with a polynomial and subtract the
;   resulting fit from the data section.
;
; Passed Parameters:
;   real   array         I/O: 2-D image array
;   int    xdata1          I: first column of the data section
;   int    xdata2          I: last column of the data section
;   int    xbias1          I: first column of the bias section
;   int    xbias2          I: last column of the bias section
;   int    y1              I: first row of bias and data sections
;   int    y2              I: last row of bias and data sections
;   int    order           I: order of the polynomial to fit
;   bool   verbose         I: print diagnostics?
;   bool   display         I: graph results?
;   bool   interactive     I: prompt user for order change?
;
; Local variables:
;   int    width              number of columns in bias section
;   real   median_bias        smoothed version of 
;
; Internal Procedures Used:
;
; External Procedures Used:
;
; Modification history:
;   1999-Mar-17  gwirth   Original version
;-----------------------------------------------------------------------

if not keyword_set( ORDER) then $
  order = 0                     ; default is constant
verbose = keyword_set( VERBOSE) ; give feedback on terminal?
display = keyword_set( DISPLAY) ; make plot of fit?
interactive = keyword_set( INTERACTIVE) ; prompt user for change in order?

; verify the requested order...
min_order = 0
max_order = 10
if order lt min_order or order gt max_order then $
  message, 'received illegal value for order'

; mash the bias section...
dx = xbias2-xbias1+1
dy = y2-y1+1
median_bias = fltarr( dy)
for i = 0,dy-1 do begin
    median_bias[i] = median( array[xbias1:xbias2,y1+i])
endfor

; begin loop until user accepts result...
continue = 1L
while continue do begin

    ; fit the bias with a polynomial...
    x = findgen(dy)+y1
;;    coeff = svdfit( x, median_bias, order+1, /LEGENDRE, YFIT=yfit)
    coeff = poly_fit( x, median_bias, order, yfit)
    y = transpose(yfit)

    ; plot...
    if display then begin
        plot, x, median_bias, title="Bias Region and Fit",  $
          xtitle="Row Number", ytitle="DN", /ynozero, xstyle=1
        oplot, x, y, linestyle=1
    endif

    if interactive then begin
        bad = 1L
    endif else begin
        bad = 0L
        continue = 0L
    endelse

    ; prompt...
    while bad do begin
        buf = ''
        read, 'Enter new order or Q to quit: ', buf
        if strlowcase(strmid(buf,0,1)) eq 'q' then begin
            bad = 0L
            continue = 0L
        endif else begin
            reads, buf, new_order, format='(i)'
            if new_order ge min_order and new_order le max_order then begin
                order = new_order
                bad = 0L
            endif else begin
                message, 'legal values for order are between ' +  $
                  string(min_order) + ' and ' + string(max_order), /inf
            endelse
        endelse
    endwhile
endwhile

; subtract bias...
;;array[xdata1:xdata2,y1:y2] = array[xdata1:xdata2,y1:y2] - y
for i=xdata1,xdata2 do begin
    array[i,y1:y2] = array[i,y1:y2] - y
endfor

end
