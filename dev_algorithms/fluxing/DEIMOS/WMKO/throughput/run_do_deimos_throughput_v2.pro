;------------------------------------------------------------------------
pro run_do_deimos_throughput_v2
;------------------------------------------------------------------------
; This is the general driver script.  Fill in names of input files and
; let 'er rip!
;
; This version is for reducing a set of data with BD+284211 on 1200G
;
; REQUIREMENTS:
;       IDL must be invoked using the command
;               ~dmoseng/bin/do_deimos_throughput
;       in order to properly configure the IDL PATH
;------------------------------------------------------------------------

input = [ '2003may30/d0530_0126.fits.gz', $
;          '2003may30/d0530_0127.fits.gz', $
;          '2003may30/d0530_0128.fits.gz', $
;          '2003may30/d0530_0129.fits.gz', $
;          '2003may29/d0529_0112.fits.gz', $
;          '2003may28/d0528_0115.fits.gz', $
;          '2003may28/d0528_0116.fits.gz', $
;          '2003may28/d0528_0117.fits.gz', $
;          '2003aug27/d0827_0081.fits.gz', $
;          '2003aug27/d0827_0082.fits.gz', $
;          '2003aug27/d0827_0083.fits.gz', $
;          '2002aug27/d0827_0013.fits.gz', $
;          '2002aug27/d0827_0030.fits.gz', $
;          '2002aug27/d0827_0031.fits.gz', $
;          '2002aug27/d0827_0032.fits.gz', $
;          '2002aug27/d0827_0033.fits.gz', $
;          '2002jul26/d0726_0001.fits.gz', $
;          '2002jul26/d0726_0002.fits.gz', $
;          '2002jul26/d0726_0003.fits.gz', $
;          '2002jul26/d0726_0004.fits.gz', $
;          '2004may21/d0521_0088.fits.gz', $
;          '2004may21/d0521_0089.fits.gz', $
;          '2004may21/d0521_0090.fits.gz', $
;          '2004may23/d0523_0114.fits.gz', $
;          '2004may23/d0523_0115.fits.gz', $
;          '2004may23/d0523_0116.fits.gz', $
          '2004jun19/d0619_0109.fits.gz', $
          '2004jun19/d0619_0110.fits.gz', $
          '2004jun19/d0619_0111.fits.gz', $
          '2004jun20/d0620_0100.fits.gz', $
          '2004jun20/d0620_0101.fits.gz', $
          '2004jun20/d0620_0102.fits.gz', $
          '2004jul16/d0716_0107.fits.gz', $
          '2004jul16/d0716_0108.fits.gz', $
          '2004jul16/d0716_0109.fits.gz', $
          '2005apr03/d0403_0094.fits.gz', $
          '2005apr03/d0403_0095.fits.gz', $
          '2005apr03/d0403_0096.fits.gz', $
          '2005jun05/d0605_0111.fits.gz', $
          '2005jun05/d0605_0112.fits.gz', $
          '2005jun05/d0605_0113.fits.gz', $
          '2005may08/d0508_0115.fits.gz', $
          '2005may08/d0508_0116.fits.gz', $
          '2005may08/d0508_0117.fits.gz', $
          '2005may10/d0510_0108.fits.gz', $
          '2005may10/d0510_0109.fits.gz', $
          '2005may10/d0510_0110.fits.gz', $
          '2005aug27/d0827_0036.fits.gz', $
          '2005aug27/d0827_0037.fits.gz', $
          '2005aug27/d0827_0038.fits.gz', $
          '2005aug27/d0827_0039.fits.gz', $
          '2005sep05/d0905_0059.fits.gz', $
          '2005sep05/d0905_0060.fits.gz', $
          '2005sep05/d0905_0061.fits.gz', $
          '2005sep06/d0906_0070.fits.gz', $
          '2005sep06/d0906_0071.fits.gz', $
          '2005sep06/d0906_0072.fits.gz', $
          '2006may31/d0531_0083.fits.gz', $
          '2006may31/d0531_0084.fits.gz', $
          '2006may31/d0531_0085.fits.gz', $
          '2006sep18/d0918_0057.fits.gz', $
          '2006sep18/d0918_0058.fits.gz', $
          '2006sep18/d0918_0059.fits.gz', $
          '2007nov11/d1111_0054.fits.gz', $
          '2007nov11/d1111_0055.fits.gz', $
          '2007nov11/d1111_0056.fits.gz']

input = '/s/sdata1001/dmoseng/kroot/throughput/raw/' + input

;; do_deimos_throughput, input=input, /DISPLAY, /CLOBBER, /VERBOSE
do_deimos_throughput, input=input, /CLOBBER, /VERBOSE
end


