pro do_make_throughput_iraf

in = ['B0049.0001.fits','R0049.0001.fits']
in = '/h/scratch13/gwirth/deimos-throughput/' + in
outfile = 'spec0049.fits'
side=['B','R']
gratenam='830G'
cenlam=8000.
exptime=44.96680069  

make_throughput_iraf, IN=in, OUTFILE=outfile, $
  SIDE=side, GRATENAM=gratenam, CENLAM=cenlam, $
  EXPTIME=exptime, /DISPLAY

in = ['B0050.0001.fits','R0050.0001.fits']
in = '/h/scratch13/gwirth/deimos-throughput/' + in
outfile = 'spec0050.fits'
side=['B','R']
gratenam='830G'
cenlam=7000.
exptime = 45.06599808 

make_throughput_iraf, IN=in, OUTFILE=outfile, $
  SIDE=side, GRATENAM=gratenam, CENLAM=cenlam, $
  EXPTIME=exptime, /DISPLAY

end
