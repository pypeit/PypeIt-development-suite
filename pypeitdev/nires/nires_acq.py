#!/usr/local/anaconda3-5.0.0.1/bin/1python

"""
Show difference MOSFIRE acquisition frames.
"""
import sys
import os
import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import argparse
#from IPython import embed

sys.path.append(r'/home/jwalawender/.local/lib/python3.6/site-packages/pyds9-1.9.dev0-py3.6-linux-x86_64.egg/')
from pyds9 import DS9

def parse_args(options=None):
    parser = argparse.ArgumentParser(description='Display NIRES SCAM cquisition images via DS9', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('rawpath', type = str, help = 'Raw data path, e.g. /s/sdata1500/nires3/2020dec26')
    parser.add_argument('prefix', type = str, help = 'Image prefix, e.g. v201226_')
    parser.add_argument('obj_frame', type = int, help = 'Object frame #, e.g. 146 for v201226_0146.fits')
    parser.add_argument('sky_frame', type = int, help = 'Sky frame #, e.g. 145 for v201226_0145.fits')
    parser.add_argument('--cut_min', type = float, default=None, help = 'Specify lower cut level in ADU. If not set uses sig_min')
    parser.add_argument('--cut_max', type = float, default=None, help = 'Specify upper cut level in ADU. If not set uses sig_max')
    parser.add_argument('--sig_min', type = float, default=3.0, help = 'Specify lower cut level in units of b/g std deviation')
    parser.add_argument('--sig_max', type = float, default=3.0, help = 'Specify upper cut level in units of b/g std deviation')
    return parser.parse_args() if options is None else parser.parse_args(options)



def main(args):

    obj_file = os.path.join(args.rawpath, '{:s}{:04d}.fits'.format(args.prefix, args.obj_frame))
    sky_file = os.path.join(args.rawpath, '{:s}{:04d}.fits'.format(args.prefix, args.sky_frame))

    obj = fits.getdata(obj_file)
    sky = fits.getdata(sky_file)
    diff = obj - sky

    mean_sky, med_sky, sigma_sky = sigma_clipped_stats(diff, sigma_lower=3.0, sigma_upper=3.0)
    cut_min = med_sky - args.sig_min*sigma_sky if args.cut_min is None else args.cut_min
    cut_max = med_sky + args.sig_max*sigma_sky if args.cut_max is None else args.cut_max
    print('cut_min={:5.3f}, cut_max={:5.3f}'.format(cut_min, cut_max))

    d = DS9()
    d.set_np2arr(diff)
    d.set('scale limits {:5.3} {:5.3}'.format(cut_min, cut_max))
    d.set('scale linear')
    d.set('zoom to fit')
    d.set('zoom 6')
    d.set('pan -300 -50')


if __name__ == '__main__':
    main(parse_args())
