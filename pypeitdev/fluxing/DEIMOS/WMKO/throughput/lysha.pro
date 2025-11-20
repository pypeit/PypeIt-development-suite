;------------------------------------------------------------------------
pro lysha
;------------------------------------------------------------------------

dir = '/s/sdata1001/dmoseng/kroot/throughput/raw/'

;; good image...
input = [dir + '2010nov05/d1105_0096.fits.gz']

do_deimos_throughput, input=input, /VERBOSE, /PSFILE, /CLOBBER


end
