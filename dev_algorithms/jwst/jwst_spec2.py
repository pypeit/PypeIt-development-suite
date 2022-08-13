

import os
from matplotlib import pyplot as plt
from astropy.io import fits

from pypeit.display import display
from jwst.pipeline import Spec2Pipeline

rawpath_level2 = '/Users/joe/jwst_redux/NIRSPEC/Raw/02736_ERO_SMACS0723_G395MG235M/level_2'

output_dir = '/Users/joe/jwst_redux/NIRSPEC/redux/calwebb'
asn_file = '/Users/joe/jwst_redux/NIRSPEC/redux/calwebb/asn_files/jwst_nrs_spec2_asn.json'


spec2 = Spec2Pipeline()
spec2.save_results = True
spec2.output_dir = output_dir
result = spec2(asn_file)
