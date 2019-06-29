import os
import numpy as np

import telluric
from flux1d import apply_sensfunc
from pypeit.core import coadd1d
from pypeit import msgs

dev_path = os.getenv('PYPEIT_DEV')
datapath = os.path.join(dev_path, 'VLT_XSHOOTER/NIR/Science')

spec1dfiles = ['spec1d_XSHOO.2016-08-02T09:57:17.147-STD,TELLURIC_XShooter_NIR_2016Aug02T095717.147.fits',
               'spec1d_XSHOO.2016-08-02T09:57:56.903-STD,TELLURIC_XShooter_NIR_2016Aug02T095756.903.fits']

fnames = os.path.join(datapath, spec1dfiles)