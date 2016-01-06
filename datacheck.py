#!/usr/bin/python2.7

'''
This code is made of one function that checks the availability of HiTS SN data on local disk.
A main function ifSNdataexist does it by calling other two functions:
- DECamSNlist
      reads the full list of HiTS SNe from 2014 and 2015 and make a dictionary where the key elements are SN nicknames
- InstrSNlist 
      reads DuPont and SOI fits file and write tables of some header elements into files saved in obsdate folders.

Note: this code can be implemented by downloading data of not available SN, i.e. passing (field, CCD, epoch) to filldata.py  
'''

### IMPORTANT: SOLVE THE ISSUE WITH EPOCH. SOLVE THE ISSUE WITH UPPERCASE (ELISE)

import os
import glob    # module to find pathnames
import sys
import numpy as np 
from astropy.io import fits


def DECamSNlist (file2014='SNHiTS2014.dat', file2015='SNHiTS2015.dat') :
	
    # load list of DECam SNe
    mat = np.loadtxt (file2014, dtype='string', delimiter='\t', skiprows=2, usecols=(1,4,5,6,7))
    SN = mat[:,0]
    field = mat[:,1]
    CCD = mat[:,2]
    icoord = mat[:,3]
    jcoord = mat[:,4]
    mat = np.loadtxt (file2015, dtype='string', delimiter='\t', skiprows=2, usecols=(1,4,5,6,7))
    SN = np.hstack((SN,mat[:,0]))    
    field = np.hstack((field,mat[:,1]))
    CCD = np.hstack((CCD,mat[:,2]))
    icoord = np.hstack((icoord,mat[:,3]))
    jcoord = np.hstack((jcoord,mat[:,4]))
	# NOTE: epoch is missing in file2014 and file2015!!! #
	
	# make a dictionary	
    d = {}
    for i in range(len(SN)) :
        SN[i] = SN[i].replace(" ","")
        d[SN[i].upper()] = ['Blind14A-P_%s' %field[i], CCD[i], [int(icoord[i]), int(jcoord[i])]]
    #print d
			
    return d


def InstrSNlist (instr, obsdate) :
	
    ### make paths to science images and output folder
    
    if instr == 'DuPont' :
        instrdir = os.path.join ('DATA', instr)
    elif instr == 'SOI' :
        instrdir = os.path.join ('DATA', 'SOAR', instr)
    else :
        sys.exit('\ninstrument name is incorrect!')	
    
    obsdir = os.path.join (instrdir, obsdate)
	
    # check if path exists
    if os.path.isdir(obsdir) :
        print '\n%s data for %s are in\n%s' %(instr,obsdate,obsdir)
    else :
        sys.exit('\nobservation date is incorrect!')
    
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
        if instr == 'DuPont' :
            outf.write(b'# FILENAME OBJECT RA DEC FILTER AIRMASS EXPTIME\n')
        else :
            outf.write(b'# FILENAME OBJECT RA DEC FILTER1 FILTER2 AIRMASS EXPTIME\n')
    
    SNlist = []
		
    for f in files :
        # open the fits
        try :
            HDU = fits.open(f)
        except IOError :    # error handling for empty/corrupted fits 
            print ('Error opening %s' %f)
            continue    
        priHDU = HDU[0].header
        # write header values into file
        with open (outfile,'a') as outf:
            if instr == 'DuPont' :
                outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(f[nchar2:], priHDU['OBJECT'], 
                priHDU['RA'], priHDU['DEC'], priHDU['FILTER'], priHDU['AIRMASS'], priHDU['EXPTIME']))
            else :
                outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(f[nchar2:], priHDU['OBJECT'], 
                priHDU['RA'], priHDU['DEC'], priHDU['FILTER1'], priHDU['FILTER2'], priHDU['AIRMASS'], priHDU['EXPTIME']))
        # make a list of SN available for given instr and at given obsdate
        objname = priHDU['OBJECT'].upper()
        objname = objname.replace(' ','')
        if (objname != '') and (objname[-2] == '_') :
            objname = objname[:-2]		
        if objname not in SNlist :
            SNlist.append(objname)
    #print SNlist

    return SNlist


def ifSNdataexist (SN, instr, obsdate, file2014='SNHiTS2014.dat', file2015='SNHiTS2015.dat') :
	
    SN = SN.upper()    # Note: upper/lowercase inconsistency solved everywhere in the code
	
    refdir = os.path.join('DATA', 'DATA_CMMPIPE')
    # check if sn exist or data are available 
    SNdict = DECamSNlist(file2014, file2015)
    if SNdict.has_key(SN) :
        SNlist = InstrSNlist(instr, obsdate)
        print '\nChecking whether DECam data for %s are available...\n' %SN
        if SN in SNlist :
            # check if path exists
            sndir = os.path.join(refdir, SNdict[SN][0], SNdict[SN][1])
            if os.path.isdir(sndir) :
                print '\nFollowing files of %s are available in the local folder %s' %(SN,sndir)
                for f in os.listdir(sndir) :
                    print f
            else :
                print 'DECam reference image not available'
        else :
            print '%s was not observed in the observation date %s' %(SN,obsdate)
    else :
        sys.exit('%s does not correspond to any valid HiTS SN!' %SN)
    # check which of the SN observed at obsdate are available in DECam data
    print '\nChecking DECam data available for the other SNe observed with %s in %s...\n' %(instr,obsdate)
    for sn in SNlist :   
        if SNdict.has_key(sn) :        
            sndir = os.path.join(refdir, SNdict[sn][0], SNdict[sn][1]) 
            if os.path.isdir(sndir) :
                print '%s SN data --> available in %s' %(sn,sndir)
            else :
                print '%s SN data --> NOT available' %sn
	
	
if __name__ == "__main__" :
    
    ifSNdataexist (sys.argv[1], sys.argv[2], sys.argv[3])
    

#gethead image* OBJECT RA DEC FILTER1 FILTER2 AIRMASS EXPTIME
