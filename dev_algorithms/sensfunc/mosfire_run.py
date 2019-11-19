import os
from flux_coadd_tell import flux_tell, stack_multinight

### sensfunction for MOSFIRE
instrument = 'MOSFIRE'

basedir = '/d2/Feige/'
#basedir = os.getenv('HOME')

std_path = os.path.join(basedir,'Dropbox/PypeIt_Redux/MOSFIRE/Nov19/nov_redux/Science')
stdfile = 'spec1d_m191118_0064-GD71_MOSFIRE_2019Nov18T104704.507.fits'

tell_method = 'qso'
#tell_method = 'poly'


z_qso = 7.6
sci_path = os.path.join(basedir,'Dropbox/PypeIt_Redux/MOSFIRE/Nov19/nov_redux/Science')
spec1dfiles = ['spec1d_d0527_0080-P261_OFF_DEIMOS_2017May27T102635.318.fits',
               'spec1d_d0527_0080-P261_OFF_DEIMOS_2017May27T102635.318.fits']
objids = ['SPAT0764-SLIT0000-DET03','SPAT0764-SLIT0000-DET07']
outroot = 'J1724+1901_MOSFIRE_2019nov18'
flux_tell(sci_path, stdfile, spec1dfiles=spec1dfiles, std_path=std_path, instrument=instrument,
          outroot=outroot, objids=objids,  z_qso=z_qso, tell_method=tell_method,
          fit_region_min=[9200.0], fit_region_max=[9900.0],
          do_sens=True, do_flux=False, do_stack=False, do_tell=False, disp=False, debug=False)