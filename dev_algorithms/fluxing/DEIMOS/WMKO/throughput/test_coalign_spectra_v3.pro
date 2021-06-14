pro test_coalign_spectra_v3, infile1, infile2, PSFILE=psfile

;device, true=24
;device, decomp=0
;loadct, 13
red = gdwcolor('red')
white = gdwcolor('white')
green = gdwcolor('green')

;; verify files...
if ~ file_test(infile1) then begin
    if file_test(infile1+'.gz') then begin
        infile1 = infile1+'.gz'
    endif else begin
        message, 'ERROR: cannot locate file '+infile1
    endelse 
endif 

;; verify files...
if ~ file_test(infile2) then begin
    if file_test(infile2+'.gz') then begin
        infile2 = infile2+'.gz'
    endif else begin
        message, 'ERROR: cannot locate file '+infile2
    endelse 
endif 

;; read data...
meta1 = mrdfits( infile1, 1, /silent)
spec1 = mrdfits( infile1, 2, /silent)
x1 = spec1.wavelength
y1 = spec1.counts

meta2 = mrdfits( infile2, 1, /silent)
spec2 = mrdfits( infile2, 2, /silent)
x2 = spec2.wavelength
y2 = spec2.counts

;; start output file...
if keyword_set(PSFILE) then begin
    set_plot, 'ps'
    device, filename=psfile, /landscape, /color
    display = 0B
endif else begin
    set_plot, 'X'
    window, 0
    window, 1
    window, 2
    display = 1B
endelse 

;; cross-correlate...
x2_orig = x2
y2_orig = y2
delta = coalign_spectra( x1, y1, x2, y2, /shift, display=display)

;; plot...
if keyword_set(DISPLAY) then wset, 0
!p.font = 1
!p.charsize=1.5
plot, x1, y1, title=psfile, xtitle='Wavelength', ytitle='Flux'
oplot, x2_orig, y2_orig, color=red
oplot, x2, y2, color=green

xpos = !x.crange[0] + 0.05 * (!x.crange[1]-!x.crange[0])
ypos = !y.crange[0] + 0.95 * (!y.crange[1]-!y.crange[0])
xyouts, xpos, ypos, 'delta='+stringify(delta), align=0.

if keyword_set(DISPLAY) then cursor, xx, yy, 1

if keyword_set(PSFILE) then device, /close

end
