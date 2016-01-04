#!/usr/bin/python2.7

'''
This code is made of ??? functions:
InstrSNlist 
    reads DuPont and SOI fits file and write tables of some header elements into files saved in obsdate folders.
DECamSNlist
    reads the full list of HiTS SNe from 2014 and 2015 and make a dictionary where the key elements are SN nicknames   
'''


import os
import glob    # module to find pathnames
import sys
import numpy as np 
from astropy.io import fits


def DECamSNlist (file2014='SNHiTS2014.dat') :
	
	# load list of DECam SNe
	SN, field, CCD, icoord, jcoord = np.loadtxt (file2014, dtype='string', delimiter='\t', skiprows=2, usecols=(1,4,5,6,7), unpack=True)
	# NOTE: epoch is missing in file2014!!! #
	
	# make a dictionary	
	d = {}
	for i in range(len(SN)) :
		d[SN[i]] = ['Blind14A-P_%s' %field[i], CCD[i], [int(icoord[i]), int(jcoord[i])]]
	
	#print d
	
	return d


def InstrSNlist (instr, obsdate) :
	
	### make paths to science images and output folder
	
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
	
	SNlist = []
		
	for f in files :
		# open the header
		HDU = fits.open(f)
		priHDU = HDU[0].header
		# write header values into file
		with open (outfile,'a') as outf:
			outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(f[nchar2:], priHDU['OBJECT'], 
			priHDU['RA'], priHDU['DEC'], priHDU['FILTER1'], priHDU['FILTER2'], priHDU['AIRMASS'], priHDU['EXPTIME']))
		# make a list of SN available with given instr and at given obsdate
		if priHDU['OBJECT'] not in SNlist :
			SNlist.append (priHDU['OBJECT'])
	
	#print SNlist
	
	return SNlist		
			
	
	
def ifSNdataexist (SN, instr, obsdate, file2014='SNHiTS2014.dat') :
	
	refdir = os.path.join('DATA', 'DATA_CMMPIPE')
	# check if sn exist or data are available 
	DECamSN = DECamSNlist(file2014)
	if not DECamSN.has_key(SN) :
		sys.exit('%s does not correspond to any valid HiTS SN!' %SN)
	else :
		SNlist = InstrSNlist(instr, obsdate)
		if SN not in SNlist :
			# Note: resolve uppercase inconsistency
			sys.exit('%s was not observed in the observation date %s' %(SN,obsdate))
		else :
			# check if path exists
			sndir = os.path.join(refdir, DECamSN[SN][0], DECamSN[SN][1]) 
	        if os.path.isdir(sndir) :
				print 'Following files of SN with nickname %s are available in DECam folder\n%s\n' %(SN,sndir)
				for f in os.listdir(sndir) :
					print f
	        else :
	            sys.exit('DECam reference image not available')
	# check which of the SN observed at obsdate are available in DECam data
	for sn in SNlist :   
	    if DECamSN.has_key(sn) :        
	        sndir = os.path.join(refdir, DECamSN[sn][0], DECamSN[sn][1]) 
	        if os.path.isdir(sndir) :
				print 'Data of SN with nickname %s are available in DECam folder\n%s' %(sn,sndir)
            else :
				print 'Data of SN with nickname %s are NOT available' %sn
	
	
if __name__ == "__main__" :
    
    ifSNdataexist (sys.argv[1], sys.argv[2], sys.argv[3])
    

#gethead image* OBJECT RA DEC FILTER1 FILTER2 AIRMASS EXPTIME
