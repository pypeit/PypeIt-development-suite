#!/bin/env python2.4

'''
Name:
  mosfireGenerateSlitRegions

Purpose:
  Given a MOSFIRE FITS image, output a file with ds9 regions to indicate the 
  expected pixel positions of the various CSU bars on the direct image.

Syntax:
  mosfireGenerateSlitRegions <fits_file> <reg_file>

Parameters:
  fits_file = name of MOSFIRE FITS image to parse
  reg_file = name of ds9 regions file to generate

Author:
  Gregory D. Wirth, WMKO

Modification history:
  2013-May-30	GDW	Original version

'''

import logging as lg
import numpy as np
import os
import pyfits
import sys

tempscale = 0.99646 
mm = 1
demagnification = 7.24254
center_pix = (1042.986, 1035.879)

def mm_to_pix(mm):
	return mm/demagnification/0.018

# following code comes from old trunk of MOSFIRE code tree:
# http://code.google.com/p/mosfire/source/browse/MOSFIRE/CSU.py?r=823640967a0497bd7bafbb0f8147228c3d3993bd 
def csu_mm_to_pix(x_mm, slitno):
	'''Convert a slit's position into a pixel value. This is a linear approximation to a sixth order polynomial fit by ccs from cooldown 8 data.'''

	numslits = 46
	rotation = np.radians(0.2282)

	# _kfp is keck focal plane
	# The x_kfp has a "fudge factor"
	x_kfp = x_mm*tempscale  - 5.0*1.2 *mm
	y_kfp = 5.8*mm * tempscale * (numslits - slitno + 0.70)

	# _mfp is Mosfire focal plane
	# Convert the mms into pixels
	x_mfp = 2047 - mm_to_pix(x_kfp)
	y_mfp = mm_to_pix(y_kfp)

	# Rotate around the center
	x_mfp -= center_pix[0]
	y_mfp -= center_pix[1]

	x_mfp = np.cos(rotation)*x_mfp - np.sin(rotation)*y_mfp
	y_mfp = np.sin(rotation)*x_mfp + np.cos(rotation)*y_mfp
	x_mfp += center_pix[0]
	y_mfp += center_pix[1]

	return np.array([x_mfp, y_mfp])

# End of MOSFIRE.CSU code

lg.basicConfig(filename=None, level=lg.INFO, 
	format='%(asctime)s %(name)s %(levelname)s %(filename)s %(lineno)d: %(message)s')

# verify args...
args = sys.argv
args.pop(0)
if len(args) != 2 :
    print "Usage: %s fits_file reg_file" % args[0]
    sys.exit(1)

# parse input args...    
if len(args) >= 1 :
    fits_file = args.pop(0)
if len(args) >= 1 :
    reg_file = args.pop(0)

# confirm existence of input image...
if not os.path.exists(fits_file):
    lg.error( "Specified input file '%s' not found" % fits_file)
    sys.exit(1)

# open the FITS file for reading...
try: 
    f = pyfits.open(fits_file)
    header = f[0].header
except: 
    lg.error("Could not open '%s' for fits IO" % fits_file)
    sys.exit(1)

# confirm non-existence of output image...
if os.path.exists(reg_file):
    lg.error( "Operation would clobber existing file '%s'" % reg_file)
    sys.exit(1)

# open output file for writing...
try:
	o = open( reg_file, 'w')
except:
	lg.error( "Can't open '%s' for writing" % reg_file)

# write header info...
o.write( 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=1 edit=0 move=0 delete=1 include=1 source=1')
o.write( "\nphysical\n")

# loop over bars...
for i in xrange(1,93):

    # put marker on bar position...
    slit = int(i+1)/2
    pos = header["B%0.2iPOS" % i]	
    xpix, ypix = csu_mm_to_pix(pos, slit)
    o.write( "circle(%f, %f, 1)\n" % (xpix, ypix))

    # place label to left or right of marker, depending on whether this is an odd or even bar...
    if (i % 2 == 0):
        xtext = xpix - 20
    else:
        xtext = xpix + 20
    o.write( "text(%f, %f) # text = {%s}\n" % (xtext, ypix, i))
        
# close files...
f.close
o.close
