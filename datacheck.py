#!/usr/bin/python2.7
'''
A main function check_by_coordinates calls other two functions:
- DECamSNlist
      reads the full list of HiTS SNe from 2014 and 2015 and make a dictionary where the key elements are SN nicknames
- InstrSNlist 
      reads DuPont and SOI fits file and write tables of some header elements into files saved in obsdate folders.
The fuction check_by_coordinates search around the coordinates of HiTS SNe, in a circle of 1 arcmin radius,
to match with the coordinates in the DuPont and SOI fits file headers.
 
To run the code write on command line: python datacheck.py <instrument> <observation date>.
e.g. 
>>>python datacheck.py SOI 20140324
'''

import os
import glob    # module to find pathnames
import sys
import numpy as np 
import datetime as dt    # class for manipulating date and time
from shutil import copyfile
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky    # classes for matching catalogs by coordinates 


def DECamSNlist (file2014='SNHiTS2014.dat', file2015='SNHiTS2015.dat') :
	
    # load list of DECam SNe
    mat = np.loadtxt (file2014, dtype='string', delimiter='\t', skiprows=2, usecols=(0,2,3,4,5,6,7,8))
    SN = mat[:,0]
    ra = mat[:,1]
    dec = mat[:,2]
    field = mat[:,3]
    CCD = mat[:,4]
    epoch = mat[:,5]
    icoord = mat[:,6]
    jcoord = mat[:,7]
    mat = np.loadtxt (file2015, dtype='string', delimiter='\t', skiprows=2, usecols=(0,2,3,4,5,6,7,8))
    SN = np.hstack((SN,mat[:,0])) 
    ra = np.hstack((ra,mat[:,1]))
    dec = np.hstack((dec,mat[:,2]))   
    field = np.hstack((field,mat[:,3]))
    CCD = np.hstack((CCD,mat[:,4]))
    epoch = np.hstack((epoch,mat[:,5]))
    icoord = np.hstack((icoord,mat[:,6]))
    jcoord = np.hstack((jcoord,mat[:,7]))
		
	# make a dictionary	
    d = {}
    for i in range(len(SN)) :
        ### old issue with sn nickname syntax
        #SN[i] = SN[i].replace(" ","")
        #d[SN[i].upper()] = ['Blind14A_%s' %field[i], CCD[i], epoch[i], [int(icoord[i]), int(jcoord[i])]]
        ##################################
        d[SN[i]] = [ra[i], dec[i], 'Blind14A_%s' %field[i], CCD[i], epoch[i], [int(icoord[i]), int(jcoord[i])]]
    #print d
			
    return d


def InstrSNlist (instr, obsdate) :
	
    ### make paths to science images and output folder
    
    if instr == 'DuPont' :
        instrdir = os.path.join ('rawDATA', instr)
    elif instr == 'SOI' :
        instrdir = os.path.join ('rawDATA', 'SOAR', instr)
    else :
        sys.exit('\ninstrument name is incorrect!')	
    
    obsdir = os.path.join (instrdir, obsdate)
	
    # check if path exists
    if os.path.isdir(obsdir) :
        print '\n%s data for %s are in\n%s' %(instr,obsdate,obsdir)
    else :
        sys.exit('\nobservation date is incorrect!')
    
    outdir = os.path.join ('procDATA', 'SOAR', 'SOI', obsdate)
    if not os.path.exists(outdir) :
    	os.makedirs(outdir)
    
    ### make/read a table of fits for a given obsdate
	
    # list the fits files in obsdir
    pathname = os.path.join(obsdir,'*.fits')
    files = sorted(glob.glob(pathname))
    # open the obsdate file and write the headline
    outfile = os.path.join(outdir, '%s.dat' %obsdate)
    print ('\nWriting into %s' %outfile)
    with open(outfile,'wb') as outf:
        if instr == 'DuPont' :
            outf.write(b'# FILENAME OBJECT RA DEC FILTER AIRMASS EXPTIME\n')
        else :
            outf.write(b'# FILENAME OBJECT RA DEC FILTER1 FILTER2 AIRMASS EXPTIME\n')
    
    # writing the obsdate file and returning lists of fits, ra and dec
    list_of_fits = []
    ra = []
    dec = []
		
    for f in files :
        # open the fits
        try :
            HDU = fits.open(f)
        except IOError :    # error handling for empty/corrupted fits 
            print ('Error opening %s' %f)
            continue    
        priHDU = HDU[0].header
        # write header values into file
        #with open (outfile,'a') as outf:
        #    if instr == 'DuPont' :
        #        outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(os.path.basename(f), priHDU['OBJECT'], 
        #        priHDU['RA'], priHDU['DEC'], priHDU['FILTER'], priHDU['AIRMASS'], priHDU['EXPTIME']))
        #    else :
        #        outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(os.path.basename(f), priHDU['OBJECT'], 
        #        priHDU['RA'], priHDU['DEC'], priHDU['FILTER1'], priHDU['FILTER2'], priHDU['AIRMASS'], priHDU['EXPTIME']))
        # make a list of fits available for given instr and at given obsdate
        list_of_fits.append(f)
        ra.append(priHDU['RA'])
        dec.append(priHDU['DEC'])
        
    return list_of_fits, ra, dec, outfile


def check_by_coordinates (instr, obsdate, file2014='SNHiTS2014.dat', file2015='SNHiTS2015.dat') :
	
    # load the dictionary of DECam SNe 
    SNdict = DECamSNlist(file2014, file2015)
    SN = []
    raDECam = []
    decDECam = []
    for sn, val in SNdict.iteritems() :
        SN.append(sn)
        raDECam.append(val[0])
        decDECam.append(val[1]) 
    SN = np.array(SN)
    raDECam = np.array(raDECam)
    decDECam = np.array(decDECam)
    # load fits files and their coordinates at obsdate 
    fitsfile, raSOI, decSOI, outfile = InstrSNlist (instr, obsdate)
    fitsfile = np.array(fitsfile)
    raSOI = np.array(raSOI)
    decSOI = np.array(decSOI)
    
    # search matching coordinates around a separation of 1 arcsminute
    c = SkyCoord (ra=raSOI, dec=decSOI, unit=(u.hourangle, u.deg))
    catalog = SkyCoord (ra=raDECam, dec=decDECam, unit=(u.hourangle, u.deg))
    idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(c, 1*u.arcmin)
    
    # updating OBJECT card in fits file headers with SNHiTS name
    date = dt.date.today().strftime('%m/%d/%Y')
    print ('\nUpdating fits file headers on %s...' %date)
    i = 0
    for f in fitsfile[idxc] :    # loop only in the list of files that must be updated
        '''
        copying the SOAR raw file into doreduce input directory,
        where it can be processed.
        '''
        inf = f.replace ('rawDATA', 'procDATA')
        if not os.path.exists(inf) :
            print '\nCopying %s\nin %s' %(f,inf)
            copyfile (f, inf)
        else :
            print '\n%s already exist' %(inf)
        # open the fits (they are already selected without IOError)
        HDU = fits.open(inf, mode='update')    
        priHDU = HDU[0].header
        newname = SN[idxcatalog[i]]
        if priHDU['OBJECT'] != newname :    # update only if necessary
            print ('\nWriting OBJECT=%s into %s ...' %(newname, inf))
            priHDU.set ('OBJ_OLD', priHDU['OBJECT'], 'object name given during the observation', after='OBJECT')
            priHDU.set ('OBJECT', newname, 'object name updated on %s' %date)
            #priHDU.add_history('header updated on %s' %date)
        else :
            print ('\nSkipping %s ...' %inf)
        HDU.close()
        i += 1
    
    # saving header info of the fits file for the obsdata 
    for f in fitsfile :
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
                outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(os.path.basename(f), priHDU['OBJECT'], 
                priHDU['RA'], priHDU['DEC'], priHDU['FILTER'], priHDU['AIRMASS'], priHDU['EXPTIME']))
            else :
                outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(os.path.basename(f), priHDU['OBJECT'], 
                priHDU['RA'], priHDU['DEC'], priHDU['FILTER1'], priHDU['FILTER2'], priHDU['AIRMASS'], priHDU['EXPTIME']))   

	
if __name__ == "__main__" :
    
    check_by_coordinates (sys.argv[1], sys.argv[2])
    

#gethead image* OBJECT RA DEC FILTER1 FILTER2 AIRMASS EXPTIME
