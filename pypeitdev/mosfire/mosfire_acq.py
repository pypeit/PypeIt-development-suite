#!/usr/local/anaconda/bin/python3

"""
Show difference MOSFIRE acquisition frames.
"""
import sys
import os
import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import argparse

sys.path.append(r'/home/jwalawender/.local/lib/python3.6/site-packages/pyds9-1.9.dev0-py3.6-linux-x86_64.egg/')
from pyds9 import DS9

def parse_args(options=None):
    parser = argparse.ArgumentParser(description='Display MOSFIRE acquisition images via DS9', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('rawpath', type = str, help = 'Raw data path, e.g. /s/sdata1300/mosfire5/2020dec08')
    parser.add_argument('prefix', type = str, help = 'Image prefix, e.g. m201208_')
    parser.add_argument('sky_frame', type = int, help = 'Sky frame #, e.g. 145 for m201208_0145.fits')
    parser.add_argument('obj_frame', type = int, help = 'Object frame #, e.g. 146 for m201208_0146.fits')
    parser.add_argument('--cut_min', type = float, default=None, help = 'Specify lower cut level in ADU. If not set uses sig_min')
    parser.add_argument('--cut_max', type = float, default=None, help = 'Specify upper cut level in ADU. If not set uses sig_max')
    parser.add_argument('--sig_min', type = float, default=3.0, help = 'Specify lower cut level in units of b/g std deviation')
    parser.add_argument('--sig_max', type = float, default=3.0, help = 'Specify upper cut level in units of b/g std deviation')
    return parser.parse_args() if options is None else parser.parse_args(options)



def main(args):


    sky_file = os.path.join(args.rawpath, '{:s}{:04d}.fits'.format(args.prefix, args.sky_frame))
    obj_file = os.path.join(args.rawpath, '{:s}{:04d}.fits'.format(args.prefix, args.obj_frame))

    obj = fits.getdata(obj_file)
    sky = fits.getdata(sky_file)
    diff = obj - sky

    # Get the sky level in the window
    ny, nx = diff.shape
    # define crude y, x boundaries
    x_coord, y_coord = np.meshgrid(np.arange(nx), np.arange(ny))
    upper_left_y, upper_left_x = (1070, 1034)
    upper_right_y, upper_right_x = (1070, 1051)
    lower_left_y, lower_left_x = (1035, 1037)
    lower_right_y, lower_right_x = (1032, 1055)

    median_box = (x_coord > lower_left_x) & (x_coord < upper_right_x) & (y_coord > lower_left_y) & (y_coord < upper_right_y)
    mean_sky, med_sky, sigma_sky = sigma_clipped_stats(diff[median_box], sigma_lower=3.0, sigma_upper=3.0)
    cut_min = med_sky - args.sig_min*sigma_sky if args.cut_min is None else args.cut_min
    cut_max = med_sky + args.sig_max*sigma_sky if args.cut_max is None else args.cut_max

    d = DS9()
    d.set_np2arr(diff)
    d.set('scale limits {:5.3} {:5.3}'.format(cut_min, cut_max))
    d.set('scale linear')
    d.set('zoom to fit')
    d.set('zoom 10')


if __name__ == '__main__':
    main(parse_args())
