import os
from flux_coadd_tell import flux_tell, stack_multinight

### sensfunction for MOSFIRE
instrument = 'MOSFIRE'

basedir = '/d2/Feige/'
#basedir = os.getenv('HOME')

std_path = os.path.join(basedir,'Dropbox/PypeIt_Redux/MOSFIRE/Nov19/nov_redux/Science')
stdfile = 'spec1d_m191118_0064-GD71_MOSFIRE_2019Nov18T104704.507.fits'
tell_method = 'qso'

## J0313-1806
sci_path = os.path.join(basedir,'Dropbox/PypeIt_Redux/MOSFIRE/Nov19/nov_redux/Science')
z_qso = 7.6
fileroot = 'J0313-1806_FIRE'
flux_tell(sci_path, stdfile, std_path=std_path,instrument=instrument, fileroot=fileroot, z_qso=z_qso, tell_method=tell_method,
          do_sens=True, do_flux=False, do_stack=False, do_tell=False, disp=False, debug=False)
