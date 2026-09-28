pro test_coalign_spectra

red = 255
green = 160
cyan = 100
blue = 50
yellow = 208
orange = 220

infile = '/s/sdata1001/dmoseng/kroot/throughput/extract/1200G/2004jun19_d0619_0109.fits'
psfile = 'test_coalign_spectra.ps'
  
;; read data...
meta = mrdfits( infile, 1)
spec = mrdfits( infile, 2)

x = spec.wavelength
y = spec.counts

set_plot, 'ps'
device, filename=psfile, /landscape, /color
plot, x, y, font=1, charsize=1.5

loadct, 13

;; shift...
nshift = 100
x2 = x
y2 = shift( y, nshift)
ny = n_elements(y)
if nshift gt 0 then begin
    y2[0:nshift-1] = y2[nshift]
endif else begin
    y2[ny-1+nshift:ny-1] = y2[ny-2+nshift]
endelse
oplot, x2, y2, color=red

;; cross-correlate...
coalign_spectra, x, y, x2, y2
plot, x, y, font=1, charsize=1.5
oplot, x2, y2, color=red


device, /close

end
