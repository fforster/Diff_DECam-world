#!/usr/bin/python2.7

'''
This module (may) contains the following functions:
- subsample
select manually a sub-sample from a given one
- jackknife 
accepts a sample and performs a left-one out resampling
- jackCI
accepts the flux computed for the original sample and all the fluxes computed for all the sub-samples 
to compute the confidence interval at 95%
- bootstrap
accepts a sample and performs a reampling with repetition
- bootCI 
accepts the flux computed for the original sample and all the fluxes computed for all the random samples 
to compute the confidence interval at 95%
'''


import os
import numpy as np
import numpy.random as npr
from astropy.stats import bootstrap

def jack_res (v) :
	
    print np.shape(v)
	
    if v.ndim == 1 :
        for i in range(len(v)) :
            if i==0 :
                matv = np.delete (v, i)
            else :
                vtemp = np.delete (v, i)
                matv = np.vstack ((matv, vtemp))
        print np.shape(matv)
        return matv
    
    elif v.ndim == 3 :
        for i in range(np.shape(v)[2]) :
            if i==0 :
                matv = np.delete (v, i, axis=2)
            else :
                vtemp = np.delete (v, i, axis=2)
                matv = np.vstack ((matv, vtemp))
        print np.shape(matv)
        print np.shape(matv[0,:,:])
        return matv

def jackknife (psf1s, psf2s, f1sel, e_f1sel, f2sel, e_f2sel, r1sel, r2sel) :
	
	mpsf1s = jack_res (np.array(psf1s))
	mpsf2s = jack_res (np.array(psf2s))
	mf1sel = jack_res (np.array(f1sel))
	me_f1sel = jack_res (np.array(e_f1sel))
	mf2sel = jack_res (np.array(f2sel))
	me_f2sel = jack_res (np.array(e_f2sel))
	mr1sel = jack_res (np.array(r1sel))
	mr2sel = jack_res (np.array(r2sel))
	
	return mpsf1s, mpsf2s, mf1sel, me_f1sel, mf2sel, me_f2sel, mr1sel, mr2sel


if __name__ == "__main__" :
	
	jackknife ()
