pro test_coalign_spectra_v2

red = 255
green = 160
cyan = 100
blue = 50
yellow = 208
orange = 220

;; infile = '2004jun19_d0619_0109.fits'
infile = '/s/sdata1001/dmoseng/kroot/throughput/extract/1200G/2004jun19_d0619_0109.fits'
psfile = 'test_coalign_spectra.ps'
  
;; read data...
meta = mrdfits( infile, 1)
spec = mrdfits( infile, 2)

x = spec.wavelength
y = spec.counts

set_plot, 'ps'
device, filename=psfile, /landscape, /color
;; !p.multi = [0, 1, 2]

plot, x, y, font=1, charsize=1.5
delta = 0.05*(!y.crange[1] - !y.crange[0])

loadct, 13

;; smooth with bspline...
nord = 4
fullbkpt = bspline_bkpts( x, nord=nord, nbkpts=15)
sset = bspline_iterfit( x, y, nord=nord, maxiter=0, yfit=yfit, fullbkpt=fullbkpt)
oplot, x, yfit, color=blue

;; smooth with bspline plus rejection...
nord = 4
fullbkpt = bspline_bkpts( x, nord=nord, nbkpts=15)
sset = bspline_iterfit( x, y, nord=nord, maxiter=10, upper=5., lower=1., $
                        yfit=yfit, fullbkpt=fullbkpt)
oplot, x, yfit, color=red

;; normalize...
norm = y / yfit
plot, x, norm, font=1, charsize=1.5

;; shift...
nshift = 100
n = n_elements(norm)
norm2 = shift( norm, nshift)
oplot, x, norm2, color=red

;; cross-correlate...
lag = findgen(n-6)+3
result = C_CORRELATE(norm, norm2, lag)
!p.multi = [0,1,1]
plot, lag, result, font=1, charsize=1.5
top = max( result, m)
mshift = lag[m]
oplot, [mshift, mshift], !y.crange, linestyle=1

print, nshift, mshift

device, /close

end
