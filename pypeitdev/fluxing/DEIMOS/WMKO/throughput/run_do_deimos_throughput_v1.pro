pro run_do_deimos_throughput_v1

input = ['2008oct02/d1002_0049.fits.gz', $
         '2008oct02/d1002_0050.fits.gz']

input = '/s/sdata1001/dmoseng/kroot/throughput/raw/' + input

;; do_deimos_throughput, input=input, /DISPLAY, /CLOBBER, /VERBOSE
do_deimos_throughput, input=input, /CLOBBER, /VERBOSE
end
