#!/usr/local/bin/python
# coding: latin-1

import numpy as np
import math as m
import string as str
import astroTools 
import deg2sec 

def angsep (RAa, DECa, RAb, DECb) :
	
	# define functions
    hms = astroTools.convHMS    
    dms = astroTools.convDMS
    sec = deg2sec.conv
    
    dRA = hms(RAa) - hms(RAb)    # degrees
    avDEC = ((dms(DECa) + dms(DECb))/2.)*m.pi/180.    # radians
        
    deltaRA = dRA * m.cos (avDEC)
    deltaDEC = dms(DECa) - dms(DECb)
        
    return sec(deltaRA), sec(deltaDEC) 
    

if __name__ == "__main__" :
    
    import sys
    
    angsep (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]) 
       
