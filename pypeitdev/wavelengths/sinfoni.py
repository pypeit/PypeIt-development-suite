


import os
from pypeit.core.wavecal.templates import build_template

# Data Model
# FITS table
#  wave -- Wavelength values
#  flux -- Arc spectrum flux values
#
# Meta must include BINNING of the template with 1=native
if os.getenv('PYPEIT_DEV') is not None:
    template_path = os.path.join(os.getenv('PYPEIT_DEV'), 'dev_algorithms/wavelengths/template_files/')
else:
    # print("You may wish to set the PYPEIT_DEV environment variable")
    pass


# slits = [1-4]  # 5080 -- 7820
# slits = [1-7]  # 7820 -- 9170
binspec = 1
xidl_file = os.path.join(template_path, 'SINFONI', 'SINFONI_0.25_K.sav')
outroot = 'vlt_sinfoni_K.fits'
slits = [0]
lcut = [19000.0, 25000.0]
build_template(xidl_file, slits, lcut, binspec, outroot, lowredux=True, micron=True)

