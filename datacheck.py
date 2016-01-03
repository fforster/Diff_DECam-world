#!/usr/bin/python2.7

'''
This code is made of two functions:
table_of_fits 
    reads DuPont and SOI fits file and write tables of some header elements into files saved in obsdate folders.
'''


import os
import glob    # module to find pathnames
import sys
import numpy as np 
from astropy.io import fits


def table_of_fits (instr, obsdate) :
	
	### make paths to science images and output folder
	
	#refdir = os.path.join('DATA', 'DATA_CMMPIPE')
	if instr == 'DuPont' :
		instrdir = os.path.join ('DATA', instr)
	elif instr == 'SOI' :
	    instrdir = os.path.join ('DATA', 'SOAR', instr)	
	else :
		sys.exit('instrument name is incorrect!')	
	
	obsdir = os.path.join (instrdir, obsdate)
	
	# check if path exists
	if os.path.isdir(obsdir) :
	    print obsdir
	else :
	    sys.exit('observation date is incorrect!')

	outdir = os.path.join (obsdir, 'OUT')
	if not os.path.exists(outdir) :
		os.makedirs(outdir)

	### make a table of fits for a given obsdate
	
	# length of pathnames
	nchar1 = len(instrdir) + 1
	nchar2 = len(obsdir) + 1
	# list the fits files in obsdir
	pathname = os.path.join(obsdir,'*.fits')
	files = sorted(glob.glob(pathname))
	# open the obsdate file and write the headline
	outfile = os.path.join(outdir, '%s.dat' %obsdir[nchar1:])
	with open(outfile,'wb') as outf:
		outf.write(b'# FILENAME OBJECT RA DEC FILTER1 FILTER2 AIRMASS EXPTIME\n')
	
	for f in files :
		# open the header
		HDU = fits.open(f)
		priHDU = HDU[0].header
		# write header values into file
		with open (outfile,'a') as outf:
			outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(f[nchar2:], priHDU['OBJECT'], 
			priHDU['RA'], priHDU['DEC'], priHDU['FILTER1'], priHDU['FILTER2'], priHDU['AIRMASS'], priHDU['EXPTIME']))


def DECamSNlist (file2014) :
	
	# load list of DECam SNe
	SN, field, CCD, icoord, jcoord = np.loadtxt (file2014, dtype='string', delimiter='\t', skiprows=2, usecols=(1,4,5,6,7), unpack=True)
	# NOTE: epoch is missing in file2014!!! #
	
	# make a dictionary	
	d = {}
	for i in range(len(SN)) :
		d[SN[i]] = ['Blind14A-P_%s' %field[i], CCD[i], [int(icoord[i]), int(jcoord[i])]]
	print d
	
	# check if reference images are available
	
	# check if there are SOI or DuPont images of the same SNe ?
	
	return d
	

	
if __name__ == "__main__" :
    
    table_of_fits (sys.argv[1], sys.argv[2])
    SNlist = DECamSNlist ('SNHiTS2014.dat')


#gethead image* OBJECT RA DEC FILTER1 FILTER2 AIRMASS EXPTIME
