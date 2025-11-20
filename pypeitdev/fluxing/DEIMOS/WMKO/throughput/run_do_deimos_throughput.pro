;------------------------------------------------------------------------
pro run_do_deimos_throughput
;------------------------------------------------------------------------
; This is the general driver script.  Fill in names of input files and
; let 'er rip!
;
; REQUIREMENTS:
;       IDL must be invoked using the command
;               ~dmoseng/bin/do_deimos_throughput
;       in order to properly configure the IDL PATH
;------------------------------------------------------------------------

;input = ['2008oct02/d1002_0049.fits.gz', $
;         '2008oct02/d1002_0050.fits.gz']

;input = '/s/sdata1001/dmoseng/kroot/throughput/raw/' + input

;; do_deimos_throughput, input=input, /DISPLAY, /CLOBBER, /VERBOSE
;; do_deimos_throughput, input=input, /CLOBBER, /VERBOSE

input = [ $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2002aug27/d0827_0025.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2002dec30/d1230_0035.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2002jul26/d0726_0006.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2002jul26/d0726_0007.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2002jul26/d0726_0008.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2002jul26/d0726_0009.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003aug01/d0801_0066.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003aug01/d0801_0067.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003aug01/d0801_0068.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003aug21/d0821_0149.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003jan27/a0039.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003jan27/a0052.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003jan27/a0053.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003jul30/d0730_0123.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003jul31/d0731_0069.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003jul31/d0731_0070.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003mar26/d0326_0043.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct26/d1026_0113.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct26/d1026_0114.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct26/d1026_0115.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct26/d1026_0116.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct27/d1027_0124.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct28/d1028_0089.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct28/d1028_0090.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct28/d1028_0091.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct28/d1028_0093.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct29/d1029_0108.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct29/d1029_0109.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct30/d1030_0107.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2003oct30/d1030_0109.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2004apr21/d0421_0057.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2004apr21/d0421_0058.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2004apr21/d0421_0059.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2004aug16/d0816_0092.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2004aug20/d0820_0118.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2004aug20/d0820_0119.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2004aug20/d0820_0120.fits.gz', $
        '/s/sdata1001/dmoseng/kroot/throughput/raw/2004aug21/d0821_0108.fits.gz']

;;dir = '/s/sdata1001/dmoseng/kroot/throughput/raw/2011mar01/'
;;input = ['d0301_0065.fits.gz','d0301_0067.fits.gz','d0301_0068.fits.gz','d0301_0069.fits.gz','d0301_0070.fits.gz']

;; dir = '/s/sdata1001/dmoseng/kroot/throughput/raw/2012jul29/'
;; input = ['d0729_0158.fits.gz','d0729_0159.fits.gz','d0729_0161.fits.gz','d0729_0163.fits.gz','d0729_0165.fits.gz']
;;input = ['d0729_0158.fits.gz']

input = ['/s/sdata1001/dmoseng/kroot/throughput/raw/2011jan01/d0101_0045.fits.gz', $
         '/s/sdata1001/dmoseng/kroot/throughput/raw/2010oct06/d1006_0136.fits.gz' ]

;; good image...
;; input = [dir + '2010oct06/d1006_0140.fits.gz']

;; bad trace...
;; input = [dir + '2002aug11/d0811_0005.fits' ]

;; bad image...
;;;input = [dir + '2008jul07/d0707_0008.fits.gz']

;; use this to do selected files
;; do_deimos_throughput, input=input, /VERBOSE, /CLOBBER

;; use this to do EVERYTHING...
do_deimos_throughput, /VERBOSE

end
