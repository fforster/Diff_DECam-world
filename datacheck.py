#!/usr/bin/python2.7

import os
import glob    # module to find pathnames
import sys
import numpy as np 
from astropy.io import fits


def make_a_table (instr, obsdate) :
	
	# make a path to scince images
	refdir = os.path.join('DATA', 'DATA_CMMPIPE')
	if instr == 'DuPont' :
		obsdir = os.path.join ('DATA', instr, obsdate)
	elif instr == 'SOI' :
	    obsdir = os.path.join ('DATA', 'SOAR', instr, obsdate)	
	else :
		sys.exit('instrument name is incorrect!')	
	
	# check if path exists
	if os.path.isdir(obsdir) :
	    print obsdir
	else :
	    sys.exit('observation date is incorrect!')	
	
	# make a table of objects observed at the given obsdate
	pathname = os.path.join(obsdir,'*.fits')
	files = glob.glob(pathname)
	print files
	
	
	
	



if __name__ == "__main__" :
    
    make_a_table (sys.argv[1], sys.argv[2])


#gethead image* OBJECT RA DEC FILTER1 FILTER2 AIRMASS EXPTIME
