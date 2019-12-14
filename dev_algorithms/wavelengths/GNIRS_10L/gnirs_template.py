
import os
import numpy as np
from astropy import table
from pypeit.wavecalib import templates

redux_path = '/Users/joe/python/PypeIt/pypeit/data/arc_lines/reid_arxiv/'
archives =['gemini_gnirs_10L_K.fits', 'gemini_gnirs_10L_H.fits','gemini_gnirs_10L_J.fits',
           'gemini_gnirs_10L_Y.fits']
files = [os.path.join(redux_path, file) for file in archives]

arc_table = table.Table()
wave = np.zeros((6,1022))
flux = np.zeros((6,1022))
order = np.arange(6,2,-1,dtype=int)
arc_table['flux'] = np.zeros((6,1022))
arc_table['order'] = np.zeros(6,dtype=int)
for islit, file in enumerate(files):
    table.table.read(files)
    wave[islit,:] = table['wave']
    flux[islit,:] = table['flux']

outpath = '/Users/joe/python/PypeIt/pypeit/data/arc_lines/reid_arxiv/'
outroot = 'gemini_gnirs_10mm_LBSX.fits'
templates.write_template(wave,flux,1,outpath,outroot)