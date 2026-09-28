import os
from flux_coadd_tell import flux_tell, stack_multinight

### sensfunction for MOSFIRE
instrument = 'MOSFIRE'

#basedir = '/d2/Feige/'
basedir = os.getenv('HOME')

std_path = os.path.join(basedir,'Dropbox/PypeIt_Redux/MOSFIRE/Nov19/nov18_redux/Science/')
stdfile = 'spec1d_m191118_0065-GD71_MOSFIRE_2019Nov18T104910.987.fits'

tell_method = 'qso'
#tell_method = 'poly'


z_qso = 7.1
sci_path = os.path.join(basedir,'Dropbox/PypeIt_Redux/MOSFIRE/Nov19/nov18_redux/Science/')
spec1dfiles = ['spec1d_m191118_0059-J0038-0653_OFF_MOSFIRE_2019Nov18T092342.427.fits',
               'spec1d_m191118_0060-J0038-0653_OFF_MOSFIRE_2019Nov18T093452.567.fits']
objids = ['SPAT1092-SLIT0000-DET01','SPAT1121-SLIT0000-DET01']
outroot = 'J0038m0653_MOSFIRE_2019nov18'
flux_tell(sci_path, stdfile, spec1dfiles=spec1dfiles, std_path=std_path, instrument=instrument,
          outroot=outroot, objids=objids,  z_qso=z_qso, tell_method=tell_method,
          fit_region_min=[9200.0], fit_region_max=[9900.0],sens_polyorder=8,
          do_sens=False, do_flux=True, do_stack=True, do_tell=True, show=True, debug=True)