#!/usr/bin/env python3


import glob
import os
import time
import numpy

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    file_list = glob.glob('RAW_DATA/VLT_XSHOOTER/*/XSHO*.fits.gz')

    for p in file_list:
        d,f = os.path.split(p)
        os.rename(p, os.path.join(d,f.replace('_',':')))

    print('Elapsed time: {0} seconds'.format(time.clock() - t))



