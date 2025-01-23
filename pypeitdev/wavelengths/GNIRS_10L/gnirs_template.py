
import os
import numpy as np
from astropy import table
from pypeit.wavecalib import templates

redux_path = '/Users/joe/python/PypeIt/pypeit/data/arc_lines/reid_arxiv/'
archives =['gemini_gnirs_10L_K.fits', 'gemini_gnirs_10L_H.fits','gemini_gnirs_10L_J.fits',
           'gemini_gnirs_10L_Y.fits']
files = [os.path.join(redux_path, file) for file in archives]

wave = np.zeros((4,1022))
flux = np.zeros((4,1022))
order = np.arange(3,7,1,dtype=int)
for islit, file in enumerate(files):
    arch = table.Table.read(file)
    wave[islit,:] = arch['wave']
    flux[islit,:] = arch['flux']

outpath = '/Users/joe/python/PypeIt/pypeit/data/arc_lines/reid_arxiv/'
outroot = 'gemini_gnirs_10mm_LBSX.fits'
templates.write_template(wave,flux,1,outpath,outroot, order=order)