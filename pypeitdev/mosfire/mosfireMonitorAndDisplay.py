#!/bin/env python2.4
 
'''
Name: 
	mosfireMonitorAndDisplay

Purpose:
	Display new MOSFIRE images in the ds9 image display client.
	Display the most recent image in frame buffer 1, and display
	the difference between the previous and current images in
	frame buffer 2.  If the image was a direct image, overlay the
	predicted bar positions in frame buffer 1.

Syntax:
	mosfireMonitorAndDisplay [ds9-name]

Parameters:
	ds9-name: the name of the ds9 display to connect to 
		[default="AutoDisplay"]

Logging:
	Basic logging is written to STDERR, which may be redirected to
	a logfile by the startup script.  Complete logging is written to
	file mosfireMonitorAndDisplay.log in the user's home directory.

Modification history:
	Date unknown	NPK	Original version
	2014-Jan-22	GDW	- now set up to perform a "zoom-to-fit" 
				on both the primary and delta frame 
				buffers the first time a new image is 
				read out, but not on subsequent readouts
				(to preserve zoom factor)

				- no longer relies on the unreliable
				IMAGEDONE keyword to detect new
				images.  On each iteration, the script
				simply checks whether the value of the
				'lastimage' keyword changed and
				displays if so.  This appears to have
				increased its reliability; the script
				no longer skips images.

				- in addition to the basic logging on
				STDERR, the script now creates a more
				detailed log in the user's home
				directory.

'''

import logging as lg
import numpy as np
import inspect
import os
import pyfits
import shlex
import subprocess
import sys
import time
import ds9
import MOSFIRE
import mosfireMosaic

MOSFIRE.SIMULATE=False

# define a global param list...
params = {}

# The following code is from MOSFIRE.CSU
tempscale = 0.99646 
mm = 1
demagnification = 7.24254
center_pix = (1042.986, 1035.879)
barpitch_mm = (5.8 * mm * tempscale)
def mm_to_pix(mm):
	return mm/demagnification/0.018
barpitch_pix = mm_to_pix(5.8 * mm * tempscale)

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

# extract the program name (minus extension) from argv[0]
head, tail = os.path.split(sys.argv[0])
prog_name, ext = os.path.splitext(tail)

# configure logging...
logfile = os.environ['HOME'] + '/' + prog_name + '.log'
lg.basicConfig( filename=logfile,
                     level=lg.DEBUG,
                     format='%(asctime)s %(process)6d %(levelname)-8s %(filename)s %(lineno)d: %(message)s',
                     datefmt='%Y-%m-%dT%H:%M:%S')

# define a Handler which writes INFO messages or higher to the sys.stderr
console = lg.StreamHandler()
console.setLevel(lg.INFO)

# set a format which is simpler for console use
## formatter = lg.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
formatter = lg.Formatter('%(asctime)s %(levelname)s %(filename)s %(lineno)d: %(message)s')

# tell the handler to use this format
console.setFormatter(formatter)

# add the handler to the root logger
lg.getLogger('').addHandler(console)
lg.info( '---------- Starting up ----------')

mds = MOSFIRE.mds()
mgs = MOSFIRE.mosfire()
mosaic = mosfireMosaic.MosfireMosaic()

def path_to_lastfile():
	'''
	Purpose: Return the full disk name of the last file written
	'''

	try: 
		path = mds.getval("lastfile").rstrip()
		return path
	except: 
		lg.error("Could not poll mds for lastfile")
		print "MDS may be down?"
		return None

def path_to_previousfile():
	''' 
	Purpose: return the path of file to subtract from mgs.  if 
	value is "last", use delta file 
	'''

	global params

	try:
		path = mgs.deltafile()
		if (path == 'last'):
			return params["previous_file"]
		else:
			return path
	except:
		lg.error("Could not get deltafile from global server")
		print "Global server may be down?"
		return params["previous_file"]

def ds9_overplot(ds9, header):
	''' 
	Purpose: Draw ds9 regions on the edge of CSU bars 
	'''
	lg.info("generating ds9 overplot")
	
	reg = 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=1 edit=0 move=0 delete=1 include=1 source=1'
	reg += "\nphysical\n"
	y0 = 2004.9

	for i in xrange(1,93):
		if (i % 2 == 0):
			even = 1
		else:
			even = 0
		slit = int(i+1)/2
		pos = header["B%0.2iPOS" % i]	
		xpix, ypix = csu_mm_to_pix(pos, slit)

		reg += "circle(%f, %f, 1)\n" % (xpix, ypix)
		if (even == 1):
			xtext = xpix - 20
		else:
			xtext = xpix + 20
		reg += "text(%f, %f) # text = {%s}\n" % (xtext, ypix, i)
	ds9.xpapipe("regions", reg)

def display_difference_image():
	''' 
	Purpose: Display current_file - previous file in ds9 frame 2 
	'''

	global params
	lg.info("State %s" % (inspect.stack()[0][3]))
	lg.info(params)

	ds9 = params["ds9"]
        
	try:
		current_file = params["current_file"]
		zoomflag = params["zoom2"]
		previous_file = path_to_previousfile()
		lg.info("subtracting %s from %s" % (previous_file, current_file))
	except:
		return

	if (current_file is None) or (previous_file is None):
		return

	try:
		f1 = pyfits.open(current_file)
		f2 = pyfits.open(previous_file)
	except:
		lg.error("Could not open either '%s' or '%s'" % (current_file, 
			previous_file))
		return
	
	# send timestamp
	update_heartbeat()

	try:
		delt = f1[0].data - f2[0].data
		f1.close()
		f2.close()
	except:
		lg.error("Could not take difference image")
		return

	hdu = pyfits.PrimaryHDU(delt)
	hdu.header.update("object", "%s - %s" % (current_file.split("/")[-1], 
		previous_file.split("/")[-1]), "")

	# copy over WCS information
	to_copy = ["PSCALE", "CTYPE1", "CTYPE2", "WCSDIM", "CD1_1", "CD1_2", 
		"CD2_1", "CD2_2", "LTM1_1", "LTM2_2", "WAT0_001", "WAT1_001",
		"WAT2_001", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "RADECSYS"]

	for el in to_copy:
		hdu.header.update(el, f1[0].header[el], "")
		
	hdulist = pyfits.HDUList([hdu])
	deltname = "%s/delta.fits" % (os.getenv("HOME"))

	try: 
		os.remove(deltname)
	except: 
		pass

	hdulist.writeto(deltname)
	ds9.open(deltname, 2)
	lg.debug("zoom flag is %s" % zoomflag)
	if zoomflag == "to fit":
		lg.debug("Zoom frame 2 to fit")
		ds9.xpaset("zoom to fit")
		params["zoom2"] = None
	ds9.frameno(1)

def do_mosaic():
	''' 
	Purpose: Asks ds9 to display param "current_file" 
	'''
	global params
	lg.info("State %s" % (inspect.stack()[0][3]))
	lg.info(params)
	
	lg.debug("grating_mode is '%s'" % params["grating_mode"])
	if params["grating_mode"] != "mirror":
		lg.debug("NON-imaging mode")
		return

	lg.debug("imaging mode -- invoke mosaic")
	mosaic.queueFrame(params["current_file"])

def update_heartbeat():		
	'''
	Purpose: update heartbeat with current timestamp
	'''

	try:
		mgs.heartbeat("autodisplayhb")
	except:
		lg.error("Could not set autodisplay heartbeat")
		print "MGS may be down?"

def wait_for_next_image():
	'''
	Purpose: Iterate until the current frame changes, then return.
	The name of the new file will be stored in params["current_file"].
	'''

	lg.info("State %s" % (inspect.stack()[0][3]))

	# loop until change in new image name...
	lastImage = params["current_file"]
	while lastImage == params["current_file"]:	
		lg.debug("Waiting for new image")
		update_heartbeat()
		time.sleep(1)
		lastImage = path_to_lastfile()

	# get the lastfile..
	lg.info("New image %s has arrived" % lastImage )
	params["previous_file"] = params["current_file"]
	params["current_file"] = lastImage

def display_current_file( initialize=0):
	''' 
	Purpose: Asks ds9 to display param "current_file" 
	'''

	global params

	lg.info("State %s" % (inspect.stack()[0][3]))
	lg.info(params)

	ds9 = params["ds9"]
	current_file = params["current_file"]
        zoomflag = params["zoom1"]

	if current_file is None:
		return 0

	ds9.open(params["current_file"], 1)
	try: 
		f1 = pyfits.open(current_file)
	except: 
		lg.error("Could not open '%s' for fits IO" % current_file)
		return 0

	try:
		hdr = f1[0].header
		gtpos = hdr["MGTNAME"].rstrip()
		params["grating_mode"]=gtpos
	except:
		lg.error("Could not find MGTNAME in header")
		return 0
		
	if gtpos == "mirror":
		# imaging mode
		ds9_overplot(ds9, hdr)

	if zoomflag == "to fit":
		ds9.xpaset("zoom to fit")
		# change the value so that on subsequent displays we
		# do not override the zoom factor set by the
		# observer...
		if not initialize:
			params["zoom1"] = None

	# report success...
	return 1

def startup():
	'''
	Purpose: initialize the list of global variables (params)
	'''

	global params
	lg.info("State %s" % (inspect.stack()[0][3]))
	lg.info(params)

	D = 0
	if len(sys.argv) == 2:
		displayName = sys.argv[1]
	else:
		displayName = "Autodisplay"

	D = ds9.ds9(displayName)

	params["ds9"] = D
	path = path_to_lastfile()
	params["current_file"] = path
	params["previous_file"] = path
	params["grating_mode"] = None
	params["zoom1"] = "to fit" # frame 1
	params["zoom2"] = "to fit" # frame 2

	mosaic.setDs9(D);
	mosaic.start();

#------------------------------------------------------------------------
# Main Loop
#------------------------------------------------------------------------

# initialize...
startup()
initialize = 1
display_current_file( initialize)

while True:

	# loop for next image...
	wait_for_next_image()

	# display current file and difference image...
	if display_current_file():
		display_difference_image()
		do_mosaic()
